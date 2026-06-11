"""Packing-job construction, worker-side execution, and the
single-call :func:`create_packing` user entry point.

The search loop pre-samples diameters in the main process and ships them
to workers — function handles do not round-trip reliably under spawn /
unfamiliar pickle contexts. CONFIG drives every per-job knob, including
the optional `neighbor_max`, pre-screen, and `save_packings` controls.
"""
from __future__ import annotations

import multiprocessing as mp
import os
import pathlib
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

import rcpgenerator

from ._helpers import (
    auto_neighbor_max,
    get_final_phi,
    hypersphere_volume_from_diameter,
    safe_compute_phi,
)
from .model import sample_diameters
from .overlap import overlap_report, phi_corrected


# ============================================================================
# Single-packing primitives.
# ============================================================================
def run_one_custom_packing(
    diameters,
    *,
    seed=123,
    Ndim=2,
    phi_init=0.5,
    initializer_phi=0.11,
    box=None,
    walls=None,
    neighbor_max=0,
    fix_height=False,
    size_ratio_hint=None,
    verbose=False,
):
    """Run one packing with a user-supplied diameter list.

    Steps:
      1. Initialize a dilute monodisperse packing at `initializer_phi`
         (safe; high-phi mono init can hang).
      2. Overwrite diameters with `diameters` and rescale the box so the
         realized phi matches `phi_init`.
      3. Call `.pack()` and return the packed instance.

    `neighbor_max=0` (default) delegates to the engine's per-particle,
    diameter-ratio-scaled allocation. Pass a positive int to override the
    per-particle base.
    """
    diameters = np.asarray(diameters, dtype=float)
    diameters = diameters[np.isfinite(diameters)]
    diameters = np.maximum(diameters, 1e-12)
    N = len(diameters)

    # neighbor_max <= 0 delegates to the C++ per-particle layout (base 200
    # scaled by each particle's diameter ratio), which sizes rows by actual
    # need. The previous behavior resolved a GLOBAL value via
    # auto_neighbor_max — sized for the largest particle — which the C++
    # layer then used as the per-particle BASE, giving every small particle
    # the big particle's allocation (~15 GB at N=500K, SR>=28 in 3D vs
    # ~0.8 GB delegated). Explicit positive overrides pass through as the
    # C++ base, unchanged.
    neighbor_max = max(0, int(neighbor_max))

    if box is None:
        box = [1.0] * Ndim
    if walls is None:
        walls = [0] * Ndim

    packing = rcpgenerator.Packing(
        phi=float(initializer_phi),
        N=N,
        Ndim=Ndim,
        box=list(box),
        walls=list(walls),
        fix_height=fix_height,
        dist={"type": "mono", "d": 1.0},
        neighbor_max=neighbor_max,
        seed=int(seed),
        verbose=verbose,
    )

    box_vol = float(np.prod(np.asarray(packing.box, dtype=float)))
    particle_vol = hypersphere_volume_from_diameter(diameters, Ndim)
    current_phi = particle_vol / box_vol
    if current_phi <= 0:
        raise ValueError("Current custom-diameter phi is non-positive.")
    scale = (float(phi_init) / current_phi) ** (1.0 / Ndim)
    diameters_scaled = diameters * scale

    packing.diameters = diameters_scaled.tolist()
    packing.phi = float(phi_init)

    object.__setattr__(packing, "_initialized", True)
    object.__setattr__(packing, "_needs_initialize", False)
    object.__setattr__(packing, "_packed", False)
    object.__setattr__(packing, "_needs_pack", True)

    packing.pack()
    return packing


def _prescreen_reject(diameters, ndim, phi_target, threshold):
    """Geometric pre-screen: return (reject_bool, dmax, lmin)."""
    D = np.asarray(diameters, dtype=float)
    V = hypersphere_volume_from_diameter(D, ndim)
    L_min = (V / float(phi_target)) ** (1.0 / float(ndim))
    dmax = float(D.max())
    return (dmax > L_min / float(threshold)), dmax, L_min


def create_packing(theta, MODEL, *,
                   N=500, Ndim=None, seed=0,
                   phi_init=0.30, initializer_phi=0.11,
                   neighbor_max=0,
                   compute_voronoi=False,
                   compute_overlap=True,
                   prescreen=None,
                   speed=None,
                   verbose=False):
    """Single-call packing build.

    Returns ``(packing, stats)``. ``stats`` always has:
      - ``phi``           — recomputed from (positions, diameters, box)
      - ``phi_rcp``       — rcpgenerator's internal phi (diagnostic)
      - ``phi_corr``      — ``phi * (1 - max_overlap)**Ndim``
      - ``max_overlap``   — worst fractional overlap, overlap/(r_i+r_j)
      - ``mean_overlap``  — mean fractional overlap over overlapping pairs
      - ``n_overlapping`` — number of overlapping pairs

    With ``compute_voronoi=True``, also returns per-particle local-phi
    statistics and a sum-match sanity flag.

    ``prescreen`` (optional) is a dict ``{"enabled": bool, "phi_target":
    float, "threshold": float, "verbose": bool}``; if it rejects, ``phi=0.0``
    is returned without running the C++ packer.

    ``speed`` (optional) overrides the ``RCP_SPEED`` env var for this call.
    Accepts ``"immediate"`` / ``"quick"`` / ``"patient"`` / ``"forever"``.
    ``None`` leaves the existing env value untouched.
    """
    if speed is not None:
        os.environ["RCP_SPEED"] = str(speed)
    if Ndim is None:
        Ndim = int(MODEL.get("Ndim", 2))
    D = sample_diameters(theta, MODEL, N=N)
    rng = np.random.default_rng(int(seed))
    rng.shuffle(D)

    if prescreen and prescreen.get("enabled", False):
        reject, dmax, L_min = _prescreen_reject(
            D, Ndim,
            prescreen.get("phi_target", 0.85),
            prescreen.get("threshold", 2.5),
        )
        if reject:
            if prescreen.get("verbose", False):
                print(f"[prescreen] reject: D_max={dmax:.4f} > "
                      f"L_min/{prescreen['threshold']}={L_min/prescreen['threshold']:.4f}")
            return None, {
                "phi": 0.0, "phi_rcp": 0.0, "phi_corr": 0.0,
                "max_overlap": 0.0, "mean_overlap": 0.0, "n_overlapping": 0,
                "rejected": True, "reason": "prescreen",
            }

    packing = run_one_custom_packing(
        D, seed=int(seed), Ndim=int(Ndim),
        phi_init=float(phi_init),
        initializer_phi=float(initializer_phi),
        neighbor_max=int(neighbor_max),
        size_ratio_hint=MODEL.get("size_ratio"),
        verbose=verbose,
    )

    phi = safe_compute_phi(packing)
    stats = {
        "phi": phi,
        "phi_rcp": float(get_final_phi(packing)),
    }
    if compute_overlap:
        ov = overlap_report(packing, verbose=False)
        stats["max_overlap"] = ov["max_overlap"]
        stats["mean_overlap"] = ov["mean_overlap"]
        stats["n_overlapping"] = ov["n_overlapping"]
        stats["phi_corr"] = phi_corrected(phi, ov["max_overlap"], Ndim)
    else:
        stats["max_overlap"] = 0.0
        stats["mean_overlap"] = 0.0
        stats["n_overlapping"] = 0
        stats["phi_corr"] = phi

    if compute_voronoi:
        from .voronoi import voronoi_phi_local
        v = voronoi_phi_local(packing)
        lp = v["local_phi"]
        stats.update({
            "local_phi": lp,
            "local_phi_mean": float(lp.mean()),
            "local_phi_std": float(lp.std()),
            "local_phi_min": float(lp.min()),
            "local_phi_max": float(lp.max()),
            "voronoi_sum_match": bool(v["sum_match"]),
            "voronoi_sum_match_relative_error": float(v["sum_match_relative_error"]),
        })
    return packing, stats


# ============================================================================
# Job packaging + worker-side execution.
# ============================================================================
def make_packing_job(theta, MODEL, candidate, seed, config,
                     diameters=None, generation=None, stage="unknown", N=None):
    """Build one picklable job dict."""
    if N is None:
        N = int(config["N"])
    if diameters is None:
        diameters = sample_diameters(np.asarray(theta, dtype=float), MODEL, N=N)
    return {
        "diameters": np.asarray(diameters, dtype=float).tolist(),
        "theta": np.asarray(theta, dtype=float).tolist(),
        "candidate": int(candidate),
        "seed": int(seed),
        "generation": None if generation is None else int(generation),
        "stage": stage,
        "Ndim": int(config["Ndim"]),
        "phi_init": float(config["phi_init"]),
        "initializer_phi": float(config.get("initializer_phi", 0.11)),
        "box": config.get("box", None),
        "walls": config.get("walls", None),
        "neighbor_max": int(config.get("neighbor_max", 0)),
        "fix_height": bool(config.get("fix_height", False)),
        "size_ratio_hint": MODEL.get("size_ratio"),
        # Pre-screen, packing-save, and corrected-phi flags travel with each job.
        "prescreen": {
            "enabled": bool(config.get("prescreen_enabled", False)),
            "phi_target": float(config.get("prescreen_phi_target", 0.85)),
            "threshold": float(config.get("prescreen_threshold", 2.5)),
            "verbose": bool(config.get("prescreen_verbose", False)),
        },
        "save_packings": str(config.get("save_packings", "none")).lower(),
        "save_dir": str(config.get("save_dir", ".")),
    }


def _save_packing_npz(packing, *, save_dir, generation, candidate, seed, stage):
    """Write (positions, diameters, box) for one packing under save_dir/packings/."""
    out_dir = pathlib.Path(save_dir) / "packings"
    out_dir.mkdir(parents=True, exist_ok=True)
    gen_tag = "none" if generation is None else f"{int(generation):03d}"
    path = out_dir / (
        f"gen{gen_tag}_cand{int(candidate):03d}"
        f"_seed{int(seed):d}_{stage}.npz"
    )
    np.savez(
        path,
        positions=np.asarray(packing.positions, dtype=float),
        diameters=np.asarray(packing.diameters, dtype=float),
        box=np.asarray(packing.box, dtype=float),
    )
    return str(path)


def _estimate_pairs_gb(D, Ndim, base):
    """Predicted neighbor-matrix allocation (GB) for this packing —
    mirrors the engine's per-particle cap formula so pathological
    candidates can be refused BEFORE the OS out-of-memory killer
    silently destroys the worker process."""
    D = np.asarray(D, dtype=float)
    dmin = max(float(D.min()), 1e-300)
    scaling = Ndim / 3.0
    omega = 1.0 + Ndim / 6.0
    b = base if base > 0 else (300 if Ndim == 2 else 750)
    caps = np.ceil(scaling * (50.0 * (D / dmin) ** omega + b))
    caps = np.minimum(caps, len(D))
    return float((caps + 1).sum() * 4) / 1e9


def _check_worker_mem(D, Ndim, base):
    """Returns None if within budget, else an explanatory message."""
    cap_gb = float(os.environ.get("RCP_WORKER_MEM_GB", "8"))
    est = _estimate_pairs_gb(D, Ndim, base)
    if est <= cap_gb:
        return None
    return (f"rcp_search: refusing to pack — estimated neighbor-matrix "
            f"allocation {est:.1f} GB exceeds the per-worker budget "
            f"{cap_gb:.0f} GB (RCP_WORKER_MEM_GB). This candidate's size "
            f"distribution implies very large per-particle neighbor "
            f"capacities. Raise RCP_WORKER_MEM_GB if the host has room, "
            f"lower n_workers, or tighten the prescreen.")


def _pin_worker_threads():
    """Engine threads per worker process (default 1: the pool supplies the
    parallelism — 40+ workers each spawning machine-wide OpenMP teams
    oversubscribes the host ~40x). RCP_WORKER_THREADS overrides."""
    try:
        import rcpgenerator
        rcpgenerator.set_num_threads(
            max(1, int(os.environ.get("RCP_WORKER_THREADS", "1"))))
    except Exception:
        pass


def run_packing_job(job):
    """Worker-safe one-packing entry point."""
    try:
        _pin_worker_threads()
        D = np.asarray(job["diameters"], dtype=float)
        seed = int(job["seed"])
        rng = np.random.default_rng(seed)
        rng.shuffle(D)
        Ndim = int(job["Ndim"])

        # Pre-screen.
        ps = job.get("prescreen", {}) or {}
        if ps.get("enabled", False):
            reject, dmax, L_min = _prescreen_reject(
                D, Ndim, ps.get("phi_target", 0.85), ps.get("threshold", 2.5)
            )
            if reject:
                return {
                    "success": True, "rejected": True,
                    "candidate": job["candidate"], "seed": seed,
                    "generation": job["generation"], "stage": job["stage"],
                    "theta": job.get("theta"),
                    "phi": 0.0, "phi_rcp": 0.0, "phi_corr": 0.0,
                    "max_overlap": 0.0, "mean_overlap": 0.0, "n_overlapping": 0,
                    "steps": None, "force_magnitude": None, "max_min_dist": None,
                    "Dmin": float(D.min()) if len(D) else None,
                    "Dmax": float(D.max()) if len(D) else None,
                    "D_ratio": float(D.max() / D.min()) if len(D) and D.min() > 0 else None,
                    "packing_path": None,
                }

        mem_msg = _check_worker_mem(D, Ndim, int(job.get("neighbor_max", 0)))
        if mem_msg:
            return {
                "success": False, "error": mem_msg,
                "candidate": job["candidate"], "seed": seed,
                "generation": job["generation"], "stage": job["stage"],
                "theta": job.get("theta"),
            }
        packing = run_one_custom_packing(
            D, seed=seed, Ndim=Ndim,
            phi_init=float(job["phi_init"]),
            initializer_phi=float(job["initializer_phi"]),
            box=job.get("box"),
            walls=job.get("walls"),
            neighbor_max=int(job.get("neighbor_max", 0)),
            fix_height=bool(job.get("fix_height", False)),
            size_ratio_hint=job.get("size_ratio_hint"),
            verbose=False,
        )
        d_final = np.asarray(packing.diameters, dtype=float)

        phi = safe_compute_phi(packing)
        ov = overlap_report(packing, verbose=False)
        phi_corr = phi_corrected(phi, ov["max_overlap"], Ndim)

        packing_path = None
        save_mode = str(job.get("save_packings", "none")).lower()
        if save_mode == "all":
            packing_path = _save_packing_npz(
                packing, save_dir=job["save_dir"],
                generation=job["generation"], candidate=job["candidate"],
                seed=seed, stage=job["stage"],
            )

        return {
            "success": True,
            "rejected": False,
            "error": None,
            "candidate": job["candidate"],
            "seed": seed,
            "generation": job["generation"],
            "stage": job["stage"],
            "theta": job.get("theta"),
            "phi": float(phi),
            "phi_rcp": float(get_final_phi(packing)),
            "phi_corr": float(phi_corr),
            "max_overlap": float(ov["max_overlap"]),
            "mean_overlap": float(ov["mean_overlap"]),
            "n_overlapping": int(ov["n_overlapping"]),
            "steps": int(packing.steps) if packing.steps is not None else None,
            "force_magnitude": float(packing.force_magnitude) if packing.force_magnitude is not None else None,
            "max_min_dist": float(packing.max_min_dist) if packing.max_min_dist is not None else None,
            "Dmin": float(d_final.min()) if len(d_final) else None,
            "Dmax": float(d_final.max()) if len(d_final) else None,
            "D_ratio": float(d_final.max() / d_final.min()) if len(d_final) else None,
            "packing_path": packing_path,
        }
    except Exception as e:
        return {
            "success": False,
            "rejected": False,
            "error": repr(e),
            "traceback": traceback.format_exc(),
            "candidate": job.get("candidate"),
            "seed": job.get("seed"),
            "generation": job.get("generation"),
            "stage": job.get("stage"),
            "theta": job.get("theta"),
            "phi": None, "phi_rcp": None, "phi_corr": None,
            "max_overlap": None, "mean_overlap": None, "n_overlapping": None,
            "steps": None, "force_magnitude": None, "max_min_dist": None,
            "Dmin": None, "Dmax": None, "D_ratio": None,
            "packing_path": None,
        }


# ============================================================================
# Worker pool sizing and dispatch.
# ============================================================================
def resolve_n_workers(config):
    """Pick a pool size from CONFIG, honoring the RCP_MAX_WORKERS env cap."""
    if not bool(config.get("use_parallel", True)):
        return 1
    detected = os.cpu_count() or 1
    n_workers = config.get("n_workers", "auto")
    reserve = int(config.get("reserve_cpus", 0))
    if n_workers == "auto":
        n_workers = detected - reserve
    else:
        n_workers = int(n_workers)
    cap = os.environ.get("RCP_MAX_WORKERS")
    if cap is not None:
        n_workers = min(n_workers, int(cap))
    return max(1, min(int(n_workers), detected))


def run_jobs(jobs, config):
    """Run a list of jobs; serial if small / single-worker, otherwise parallel."""
    jobs = list(jobs)
    n_workers = resolve_n_workers(config)
    print(f"Detected CPUs: {os.cpu_count()}; using workers: {n_workers}")
    if n_workers <= 1 or len(jobs) <= 1:
        return [run_packing_job(job) for job in jobs]
    try:
        # forkserver: immune to parent OpenMP state (see scheduler.py).
        ctx = mp.get_context(os.environ.get("RCP_MP_CONTEXT", "forkserver"))
        results = []
        with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as ex:
            futures = [ex.submit(run_packing_job, job) for job in jobs]
            for fut in as_completed(futures):
                row = fut.result()
                if not row.get("success") and not row.get("rejected"):
                    err = row.get("error", "no error message recorded")
                    print(f"[heads-up] one packing run failed and was "
                          f"skipped (candidate={row.get('candidate')}, "
                          f"seed={row.get('seed')}, stage={row.get('stage')}). "
                          f"The engine said:\n  {err}", flush=True)
                results.append(row)
        return results
    except Exception as e:
        print("Parallel execution failed; falling back to serial.")
        print("Error:", repr(e))
        return [run_packing_job(job) for job in jobs]


def aggregate_candidate_results(job_results, theta_by_candidate, score_by="mean_phi"):
    """Group rows by candidate, compute scores, return list sorted desc."""
    grouped = {}
    for r in job_results:
        grouped.setdefault(int(r["candidate"]), []).append(r)
    out = []
    for candidate, rows in grouped.items():
        successes = [r for r in rows if r["success"] and r["phi"] is not None]
        failures = [r for r in rows if not r["success"]]
        rejects = [r for r in rows if r.get("rejected", False)]
        phis = [float(r["phi"]) for r in successes]
        phi_corrs = [float(r["phi_corr"]) for r in successes if r.get("phi_corr") is not None]
        max_ov = [float(r["max_overlap"]) for r in successes if r.get("max_overlap") is not None]
        steps = [r["steps"] for r in successes]
        forces = [r["force_magnitude"] for r in successes]
        mmds = [r["max_min_dist"] for r in successes]
        seeds = [int(r["seed"]) for r in rows]
        if phis:
            mean_phi = float(np.mean(phis)); median_phi = float(np.median(phis))
            best_phi = float(np.max(phis)); std_phi = float(np.std(phis))
        else:
            mean_phi = -np.inf; median_phi = -np.inf
            best_phi = -np.inf; std_phi = np.inf
        if phi_corrs:
            mean_phi_corr = float(np.mean(phi_corrs))
        else:
            mean_phi_corr = -np.inf
        theta = np.asarray(theta_by_candidate[candidate], dtype=float)
        finite = lambda xs: [x for x in xs if x is not None]
        out.append({
            "candidate": int(candidate),
            "theta": theta,
            "mean_phi": mean_phi, "median_phi": median_phi,
            "best_phi": best_phi, "std_phi": std_phi,
            "mean_phi_corr": mean_phi_corr,
            "mean_max_overlap": float(np.mean(max_ov)) if max_ov else float("inf"),
            "phis": phis, "phi_corrs": phi_corrs,
            "max_overlaps": max_ov,
            "steps": steps, "forces": forces, "max_min_dists": mmds, "seeds": seeds,
            "failures": len(failures),
            "rejects": len(rejects),
            "raw_rows": rows,
            "mean_force": float(np.nanmean(finite(forces))) if finite(forces) else float("inf"),
            "mean_steps": float(np.nanmean(finite(steps))) if finite(steps) else float("inf"),
            "mean_max_min_dist": float(np.nanmean(finite(mmds))) if finite(mmds) else float("inf"),
        })
    key = "mean_phi_corr" if score_by == "phi_corr" else "mean_phi"
    return sorted(out, key=lambda r: r[key], reverse=True)
