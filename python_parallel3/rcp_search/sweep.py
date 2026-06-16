"""Continuous-distribution baseline sweeps (power / lognormal / gaussian).

Generates the φ-vs-size-ratio (and φ-vs-CV for gaussian) baseline data for
the packing-distribution study. One packing per (family, shape parameter,
size, seed) "case"; results append to a census-format JSONL whose rows are
SELF-DESCRIBING (they embed the full make_universal_model spec incl grid_n,
so each row is reconstructable without a sidecar).

Design (agreed):
* The GRID is the sweep parameters (shape param × size). N is NOT swept —
  a finite-size tolerance on D_max/L imposes N per case; when the required
  N exceeds N_cap the case is INVALID (dropped), which is what caps the
  reachable size ratio.
* Resume/idempotent: a case whose key already exists in the output census
  is skipped, so you can add parameters or raise N_cap and only the new /
  newly-valid cases run.
* 2D/3D via Ndim (sets box volume L^Ndim in the autosizer, balance rule,
  and the φ used for the jamming-compression correction).

This module holds the LIBRARY pieces (picklable, importable) so the parallel
Colab harness can use them; the parallel/memory-monitored runner is added on
top of these. The local pilot uses :func:`run_sweep_sequential`.
"""
from __future__ import annotations

import json
import math
import os
import pathlib
import time

import numpy as np

from ._helpers import hypersphere_volume_from_diameter
from .model import (center_for_geometric_mean, initialize_theta,
                    make_universal_model, sample_diameters)
from .census import realized_distribution_stats

SWEEP_CENSUS_FILENAME = "baseline_census.jsonl"

# φ used to convert the dilute (φ_init) D_max/L to the JAMMED value when
# sizing N (the box is fixed at L=1, particles grow to φ_final, so
# D_max/L_final = D_max * (φ_final / V_sum)^(1/Ndim)). Per-dimension defaults.
_PHI_JAMMED = {2: 0.86, 3: 0.80}


# ============================================================================
# Model construction per family (returns MODEL, theta, and a recoverable spec)
# ============================================================================

def build_sweep_model(family, shape_param, size_ratio, Ndim,
                      grid_n=5000, alpha_tukey=0.0025):
    """Return (MODEL, theta, spec) for one continuous baseline.

    family: 'power' | 'lognormal' | 'gaussian'
    shape_param: exponent q (power), alpha (lognormal), CV (gaussian)
    size_ratio: D_max/D_min crop (power, lognormal). For gaussian it is
                DERIVED from CV (pass None) so the distribution is contained.
    `spec` is the exact make_universal_model kwargs (JSON-able) for recovery.
    """
    bal = "equal_volume" if int(Ndim) == 3 else "equal_area"
    family = str(family).lower()

    if family == "lognormal":
        SR = float(size_ratio)
        comps = [{"type": "lognormal",
                  "init": {"alpha": float(shape_param), "cutoff": 1.0}}]
        name = f"ln_a{shape_param:g}_S{SR:g}_d{Ndim}"
        spec = dict(components=comps, size_ratio=SR, Ndim=int(Ndim),
                    alpha_tukey=float(alpha_tukey), grid_n=int(grid_n),
                    balance_rule=bal, name=name)
        MODEL = make_universal_model(**spec)
        theta = initialize_theta(MODEL, rule="from_spec")
        # geometric mean (mu) at D_g = 1 (symmetric in log about the span)
        theta[MODEL["components"][0]["slice_start"] + 3] = \
            center_for_geometric_mean(MODEL, 1.0)

    elif family == "power":
        SR = float(size_ratio)
        comps = [{"type": "power",
                  "init": {"R_low_u": 0.0, "R_high_u": 1.0,
                           "exp": float(shape_param)}}]
        name = f"pw_q{shape_param:g}_S{SR:g}_d{Ndim}"
        spec = dict(components=comps, size_ratio=SR, Ndim=int(Ndim),
                    alpha_tukey=float(alpha_tukey), grid_n=int(grid_n),
                    balance_rule=bal, name=name)
        MODEL = make_universal_model(**spec)
        theta = initialize_theta(MODEL, rule="from_spec")

    elif family == "gaussian":
        # Gaussian in DIAMETER, mean D=1, sigma_D = CV. Sweep axis is CV;
        # SR is derived so the geometric span [1/sqrt(SR), sqrt(SR)] contains
        # the gaussian's ~±4σ (the low end is the binding constraint).
        CV = float(shape_param)
        if CV >= 0.24:
            # symmetric ±4σ would cross D=0; clean handling of broad gaussians
            # (asymmetric crop at a positive floor) is deferred — flag loudly.
            raise NotImplementedError(
                f"gaussian CV={CV} >= 0.24 needs asymmetric-crop handling "
                f"(not yet built); use CV < 0.24 for now.")
        SR = float(max((1.0 + 4.0 * CV) ** 2, 1.0 / (1.0 - 4.0 * CV) ** 2) * 1.05)
        comps = [{"type": "gaussian", "init": {}}]   # filled after we know span
        name = f"gs_cv{CV:g}_S{SR:g}_d{Ndim}"
        spec = dict(components=comps, size_ratio=SR, Ndim=int(Ndim),
                    alpha_tukey=float(alpha_tukey), grid_n=int(grid_n),
                    balance_rule=bal, name=name)
        MODEL = make_universal_model(**spec)
        R_min, R_max = MODEL["R_min"], MODEL["R_max"]
        span_R = R_max - R_min
        R_c = 0.5                              # D_c = 1
        R_center_u = (R_c - R_min) / span_R
        sigma_u = (CV / 2.0) / span_R          # sigma_R = sigma_D/2 = CV/2
        spec["components"] = [{"type": "gaussian",
                               "init": {"R_center_u": float(R_center_u),
                                        "sigma_u": float(sigma_u)}}]
        MODEL = make_universal_model(**spec)
        theta = initialize_theta(MODEL, rule="from_spec")
    else:
        raise ValueError(f"unknown family '{family}'")

    return MODEL, theta, spec


# ============================================================================
# Finite-size autosizer: D_max/L tolerance -> required N (or None if > cap)
# ============================================================================

def autosize_N(MODEL, theta, Ndim, dmaxL_tol, N_cap,
               phi_jammed=None, N_floor=2000):
    """Smallest N (in [N_floor, N_cap]) whose predicted JAMMED D_max/L <= tol.

    D_max/L_final = D_max * (phi_jammed / V_sum)^(1/Ndim), V_sum = sum of
    mean-normalized particle volumes in the unit box (box fixed at L=1; the
    pack grows particles to phi_jammed). Monotone decreasing in N, so bisect.
    Returns (N, predicted_dmaxL) or (None, predicted_dmaxL_at_cap) if even
    N_cap can't meet the tolerance.
    """
    Ndim = int(Ndim)
    phi_j = float(phi_jammed if phi_jammed is not None
                  else _PHI_JAMMED.get(Ndim, 0.80))

    def dmaxL(N):
        D = sample_diameters(theta, MODEL, N=int(N))
        Vsum = float(np.sum(hypersphere_volume_from_diameter(D, Ndim)))
        return float(D.max() * (phi_j / Vsum) ** (1.0 / Ndim))

    d_cap = dmaxL(N_cap)
    if d_cap > dmaxL_tol:
        return None, d_cap                    # infeasible within the cap
    lo, hi = int(N_floor), int(N_cap)
    if dmaxL(lo) <= dmaxL_tol:
        return lo, dmaxL(lo)
    # bisect for the smallest feasible N
    for _ in range(24):
        mid = (lo + hi) // 2
        if mid <= lo:
            break
        if dmaxL(mid) <= dmaxL_tol:
            hi = mid
        else:
            lo = mid
    return hi, dmaxL(hi)


# ============================================================================
# Case identity + resume
# ============================================================================

def case_key(family, shape_param, size_ratio, seed, Ndim, dmaxL_tol):
    """Canonical key for resume/skip. N is excluded (derived); size_ratio is
    rounded so a derived gaussian SR matches stably."""
    sz = "auto" if size_ratio is None else f"{float(size_ratio):.6g}"
    return (f"{family}|p={float(shape_param):.6g}|s={sz}|seed={int(seed)}"
            f"|d={int(Ndim)}|tol={float(dmaxL_tol):.4g}")


def census_path(save_dir):
    return pathlib.Path(save_dir) / SWEEP_CENSUS_FILENAME


def load_completed_keys(save_dir):
    """Set of case keys already present in the output census (for skip)."""
    p = census_path(save_dir)
    done = set()
    if not p.exists():
        return done
    with open(p) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                row = json.loads(line)
            except Exception:
                continue
            if "case_key" in row:
                done.add(row["case_key"])
    return done


# ============================================================================
# One case -> one census row (picklable; the parallel harness reuses this)
# ============================================================================

def run_one_case(family, shape_param, size_ratio, seed, Ndim, *,
                 dmaxL_tol=0.20, N_cap=250_000, grid_n=5000,
                 alpha_tukey=0.0025, phi_init=0.30, initializer_phi=0.11,
                 neighbor_max=0, speed="quick", N_override=None,
                 phi_jammed=None):
    """Build the model, size N (or use N_override), pack, return a
    self-describing census row dict. Returns a row with status 'invalid' if
    the autosizer can't meet the tolerance within N_cap (no packing run)."""
    from .jobs import create_packing
    t_build = time.perf_counter()
    MODEL, theta, spec = build_sweep_model(
        family, shape_param, size_ratio, Ndim,
        grid_n=grid_n, alpha_tukey=alpha_tukey)
    SR = float(spec["size_ratio"])

    if N_override is not None:
        N, pred_dmaxL = int(N_override), None
    else:
        N, pred_dmaxL = autosize_N(MODEL, theta, Ndim, dmaxL_tol, N_cap,
                                   phi_jammed=phi_jammed)

    key = case_key(family, shape_param, size_ratio, seed, Ndim, dmaxL_tol)
    base = {
        "case_key": key,
        "source": "continuous_baseline",        # provenance for the merged set
        "family": family,
        "model_name": spec["name"],             # mirrors packing_census schema
        "param_name": {"power": "exp", "lognormal": "alpha",
                       "gaussian": "CV"}[family],
        "param_value": float(shape_param),
        "size_ratio": SR, "Ndim": int(Ndim), "seed": int(seed),
        "dmaxL_tol": float(dmaxL_tol), "N_cap": int(N_cap),
        "model_spec": spec,
    }

    if N is None:
        base.update(status="invalid", reason="N>N_cap",
                    predicted_dmaxL_at_cap=pred_dmaxL,
                    wall_s=time.perf_counter() - t_build)
        return base

    # realized distribution stats (analytic moments of the sampled diameters)
    stats = realized_distribution_stats(theta, MODEL, N)
    sample_cv = float(math.sqrt(max(stats["m2"] - stats["m1"] ** 2, 0.0))
                      / stats["m1"]) if stats["m1"] else 0.0

    t0 = time.perf_counter()
    packing, pstats = create_packing(
        theta, MODEL, N=N, Ndim=int(Ndim), seed=int(seed),
        phi_init=phi_init, initializer_phi=initializer_phi,
        neighbor_max=neighbor_max, speed=speed, verbose=False)
    wall = time.perf_counter() - t0

    if not isinstance(pstats, dict) or pstats.get("phi") is None:
        base.update(status="failed", reason=str(pstats)[:200], N=int(N),
                    wall_s=time.perf_counter() - t_build)
        return base

    Dfin = np.asarray(packing.diameters, dtype=float)
    L = float(np.asarray(packing.box, dtype=float)[0])
    realized_dmaxL = float(Dfin.max() / L)
    steps = getattr(packing, "steps", None)

    base.update(
        status="ok", N=int(N), rejected=False, predicted_dmaxL=pred_dmaxL,
        realized_dmaxL=realized_dmaxL,
        dmaxL_overshoot=bool(realized_dmaxL > dmaxL_tol * 1.10),
        phi=float(pstats["phi"]), phi_corr=float(pstats["phi_corr"]),
        max_overlap=float(pstats["max_overlap"]),
        steps=(int(steps) if steps is not None else None),
        sample_CV=sample_cv,
        theta=[float(x) for x in np.asarray(theta, dtype=float)],
        wall_s=float(wall),
        **{k: stats[k] for k in stats},     # m1..m6, d_min, d_max
    )
    return base


# ============================================================================
# Sequential runner (LOCAL PILOT — conservative, no parallelism)
# ============================================================================

def run_sweep_sequential(cases, save_dir, *, dmaxL_tol=0.20, N_cap=250_000,
                         Ndim=3, grid_n=5000, alpha_tukey=0.0025,
                         speed="quick", local_N_hardcap=None, verbose=True):
    """Run a list of cases one at a time, appending self-describing rows.
    Idempotent: cases already in the census are skipped.

    cases: list of dicts {family, shape_param, size_ratio, seed}. (size_ratio
    is ignored/derived for gaussian -> pass None.)
    local_N_hardcap: if set, any autosized N above it is clamped down and the
    row flagged — a SAFETY rail for the local pilot so we never launch a
    heavy packing on the dev machine.
    """
    save_dir = str(save_dir)
    pathlib.Path(save_dir).mkdir(parents=True, exist_ok=True)
    done = load_completed_keys(save_dir)
    out = census_path(save_dir)
    rows = []
    for c in cases:
        fam = c["family"]
        p = c["shape_param"]
        sr = c.get("size_ratio")
        seed = int(c.get("seed", 0))
        key = case_key(fam, p, sr, seed, Ndim, dmaxL_tol)
        if key in done:
            if verbose:
                print(f"[skip] {key} (already in census)", flush=True)
            continue

        N_override, clamped = None, False
        if local_N_hardcap is not None:
            MODEL, theta, _ = build_sweep_model(fam, p, sr, Ndim,
                                                grid_n=grid_n,
                                                alpha_tukey=alpha_tukey)
            N_auto, _ = autosize_N(MODEL, theta, Ndim, dmaxL_tol, N_cap)
            if N_auto is None:
                N_override, clamped = int(local_N_hardcap), True
            else:
                N_override = int(min(N_auto, local_N_hardcap))
                clamped = N_auto > local_N_hardcap   # only flag a true clamp

        t0 = time.perf_counter()
        row = run_one_case(fam, p, sr, seed, Ndim, dmaxL_tol=dmaxL_tol,
                           N_cap=N_cap, grid_n=grid_n, alpha_tukey=alpha_tukey,
                           speed=speed, N_override=N_override)
        if clamped:
            row["local_N_hardcapped"] = True
        with open(out, "a") as f:
            f.write(json.dumps(row) + "\n")
        done.add(key)
        rows.append(row)
        if verbose:
            st = row.get("status")
            extra = (f"phi={row['phi']:.4f} N={row['N']} "
                     f"Dmax/L={row['realized_dmaxL']:.3f} "
                     f"CV={row['sample_CV']:.3f}" if st == "ok" else
                     row.get("reason", ""))
            print(f"[done] {fam} p={p} sr={sr} -> {st} {extra} "
                  f"({time.perf_counter()-t0:.1f}s)", flush=True)
    return rows


def make_grid(family, params, sizes=None, seeds=(0,)):
    """Build a `cases` list (cross product). For gaussian, `size_ratio` is
    derived from CV, so `sizes` is ignored."""
    cases = []
    for p in params:
        if str(family).lower() == "gaussian":
            for s in seeds:
                cases.append({"family": family, "shape_param": p,
                              "size_ratio": None, "seed": int(s)})
        else:
            for sr in (sizes or []):
                for s in seeds:
                    cases.append({"family": family, "shape_param": p,
                                  "size_ratio": sr, "seed": int(s)})
    return cases


# ============================================================================
# Parallel harness (Colab): picklable worker + dry-run planner + runner
# ============================================================================

def _sweep_worker(payload):
    """Pool worker: set engine threads to T, run one case, return its row.
    Module-level + picklable so the forkserver pool can use it."""
    try:
        import rcpgenerator
        rcpgenerator.set_num_threads(max(1, int(payload.get("T", 1))))
    except Exception:
        pass
    a = dict(payload)
    key = a.pop("key", None); a.pop("T", None)
    fam = a.pop("family"); p = a.pop("shape_param")
    sr = a.pop("size_ratio"); seed = a.pop("seed"); Ndim = a.pop("Ndim")
    row = run_one_case(fam, p, sr, seed, Ndim, **a)
    if key is not None:
        row["case_key"] = key
    return row


def _build_pending(cases, Ndim, dmaxL_tol, N_cap, grid_n, alpha_tukey, done):
    """Autosize per distinct (family, param, size); return (pending, invalid
    rows, skipped count). N is cached across seeds (seed-independent)."""
    from collections import deque
    pending, invalid_rows, skipped = deque(), [], 0
    autocache = {}
    for c in cases:
        fam = c["family"]; p = c["shape_param"]
        sr = c.get("size_ratio"); seed = int(c.get("seed", 0))
        key = case_key(fam, p, sr, seed, Ndim, dmaxL_tol)
        if key in done:
            skipped += 1
            continue
        akey = (fam, float(p), "auto" if sr is None else float(sr))
        if akey not in autocache:
            M, th, _ = build_sweep_model(fam, p, sr, Ndim, grid_n=grid_n,
                                         alpha_tukey=alpha_tukey)
            autocache[akey] = autosize_N(M, th, Ndim, dmaxL_tol, N_cap)
        N, pred = autocache[akey]
        if N is None:
            row = run_one_case(fam, p, sr, seed, Ndim, dmaxL_tol=dmaxL_tol,
                               N_cap=N_cap, grid_n=grid_n,
                               alpha_tukey=alpha_tukey)
            invalid_rows.append(row)
            continue
        pending.append({"family": fam, "shape_param": p, "size_ratio": sr,
                        "seed": seed, "Ndim": Ndim, "N_override": int(N),
                        "dmaxL_tol": dmaxL_tol, "N_cap": N_cap,
                        "grid_n": grid_n, "alpha_tukey": alpha_tukey,
                        "key": key})
    return pending, invalid_rows, skipped


def plan_sweep(cases, save_dir, *, dmaxL_tol=0.20, N_cap=250_000, Ndim=3,
               grid_n=5000, alpha_tukey=0.0025, T=4, total_cpus=None,
               s_per_particle=1.3e-3):
    """DRY RUN: autosize the whole grid (no packing), report valid/invalid/
    skipped, N distribution, dropped (capped) points, and a ROUGH cost
    estimate, so cost is approved before any compute. `s_per_particle` is the
    per-case wall ≈ s_per_particle * N (recalibrate after the first runs)."""
    done = load_completed_keys(save_dir)
    pending, invalid_rows, skipped = _build_pending(
        cases, Ndim, dmaxL_tol, N_cap, grid_n, alpha_tukey, done)
    Ns = sorted(c["N_override"] for c in pending)
    n_workers = max(1, int((total_cpus or os.cpu_count() or 4) // max(1, T)))
    print(f"=== sweep plan (Ndim={Ndim}, dmaxL_tol={dmaxL_tol}, "
          f"N_cap={N_cap}) ===")
    print(f"  pending (will run) : {len(pending)}")
    print(f"  already done (skip): {skipped}")
    print(f"  invalid (N>cap)    : {len(invalid_rows)}")
    for r in invalid_rows:
        print(f"      drop {r['family']} {r['param_name']}={r['param_value']} "
              f"SR={r['size_ratio']:g} (pred D_max/L@cap="
              f"{r.get('predicted_dmaxL_at_cap', float('nan')):.2f})")
    if Ns:
        tot = sum(Ns)
        core_h = tot * s_per_particle * T / 3600.0
        wall_h = (tot * s_per_particle) / max(1, n_workers) / 3600.0
        print(f"  N: min={Ns[0]} median={Ns[len(Ns)//2]} max={Ns[-1]} "
              f"sum={tot}")
        print(f"  ROUGH cost @ T={T}, {n_workers} concurrent: "
              f"~{core_h:.1f} core-h, ~{wall_h:.1f} wall-h "
              f"(s_per_particle={s_per_particle:g}; recalibrate)")
    return {"pending": list(pending), "invalid": invalid_rows,
            "skipped": skipped, "n_workers": n_workers}


def run_sweep_parallel(cases, save_dir, *, T=4, total_cpus=None,
                       mem_budget_gb=100.0, dmaxL_tol=0.20, N_cap=250_000,
                       Ndim=3, grid_n=5000, alpha_tukey=0.0025, speed="quick",
                       max_kills_per_case=2, mp_context="fork", verbose=True):
    """Run the grid in parallel: `total_cpus // T` cases concurrent, each at
    T engine threads. Idempotent (skips done), appends self-describing rows.

    Memory: a reactive guard polls system RAM each loop (~1 s); if usage
    exceeds `mem_budget_gb` it kills the largest-RSS worker, requeues that
    case (up to `max_kills_per_case`, then records it as killed_oom and
    advises), and reports — so a too-greedy parallel run self-throttles
    instead of crashing. Killed/missing cases refill on a later re-run.
    """
    from .queue_scheduler import PreemptivePool
    try:
        import psutil
    except Exception:
        psutil = None
        if verbose:
            print("[mem-monitor] psutil unavailable — RAM guard DISABLED.",
                  flush=True)

    save_dir = str(save_dir)
    pathlib.Path(save_dir).mkdir(parents=True, exist_ok=True)
    out = census_path(save_dir)
    done = load_completed_keys(save_dir)

    def _append(row):
        with open(out, "a") as f:
            f.write(json.dumps(row) + "\n")

    pending, invalid_rows, skipped = _build_pending(
        cases, Ndim, dmaxL_tol, N_cap, grid_n, alpha_tukey, done)
    for r in invalid_rows:
        _append(r); done.add(r["case_key"])
    counts = {"ok": 0, "failed": 0, "killed": 0,
              "invalid": len(invalid_rows), "skipped": skipped}
    for c in pending:
        c["speed"] = speed; c["T"] = T
    if verbose:
        print(f"[sweep] pending={len(pending)} skipped={skipped} "
              f"invalid={len(invalid_rows)}", flush=True)
    if not pending:
        return counts

    n_workers = max(1, int((total_cpus or os.cpu_count() or 4) // max(1, T)))
    n_workers = min(n_workers, len(pending))
    # 'fork' avoids the forkserver/spawn __main__ re-execution (children fork
    # from an OMP-clean parent here); works as a script and in a notebook.
    pool = PreemptivePool(n_workers, _sweep_worker, mp_context=mp_context)
    worker_aid = [None] * n_workers
    worker_case = [None] * n_workers
    aid_case = {}
    killed_aids = set()
    kills_per_key = {}
    next_aid = 0

    def dispatch_idle():
        nonlocal next_aid
        for w in range(n_workers):
            if worker_case[w] is not None or not pending:
                continue
            case = pending.popleft()
            aid = next_aid; next_aid += 1
            worker_aid[w] = aid; worker_case[w] = case
            aid_case[aid] = case
            try:
                pool.submit(w, aid, case)
            except Exception:
                worker_case[w] = None; worker_aid[w] = None
                pending.append(case)

    dispatch_idle()
    try:
        while True:
            if not pending and not any(worker_case[w] is not None
                                       for w in range(n_workers)):
                break
            # reactive RAM guard (single-threaded -> no race with kill)
            if psutil is not None:
                vm = psutil.virtual_memory()
                used = (vm.total - vm.available) / 1e9
                if used > mem_budget_gb:
                    best, best_rss = None, -1
                    for w in range(n_workers):
                        if worker_case[w] is None:
                            continue
                        try:
                            rss = psutil.Process(
                                pool._procs[w].pid).memory_info().rss
                        except Exception:
                            rss = 0
                        if rss > best_rss:
                            best_rss, best = rss, w
                    if best is not None:
                        case = worker_case[best]; aid = worker_aid[best]
                        pool.kill(best)
                        killed_aids.add(aid); aid_case.pop(aid, None)
                        worker_case[best] = None; worker_aid[best] = None
                        k = kills_per_key.get(case["key"], 0) + 1
                        kills_per_key[case["key"]] = k
                        counts["killed"] += 1
                        print(f"[mem-monitor] RAM {used:.0f}>{mem_budget_gb:.0f} "
                              f"GB: killed {case['family']} "
                              f"p={case['shape_param']} SR~{case['size_ratio']} "
                              f"N={case['N_override']} (RSS {best_rss/1e9:.0f} "
                              f"GB). Raise T (fewer concurrent) or lower N_cap.",
                              flush=True)
                        if k < max_kills_per_case:
                            pending.append(case)
                        else:
                            row = {"case_key": case["key"],
                                   "source": "continuous_baseline",
                                   "family": case["family"],
                                   "param_value": float(case["shape_param"]),
                                   "size_ratio": case["size_ratio"],
                                   "Ndim": Ndim, "seed": case["seed"],
                                   "N": case["N_override"],
                                   "status": "killed_oom",
                                   "reason": f"OOM-killed {k}x; raise T / lower "
                                             f"N_cap"}
                            _append(row); done.add(case["key"])
                        dispatch_idle()
                        continue
            got = pool.next_result(timeout=1.0)
            if got is None:
                continue
            w, aid, outrow = got
            if aid in killed_aids:
                killed_aids.discard(aid)
                continue                      # late result from a killed case
            case = aid_case.pop(aid, None)
            worker_case[w] = None; worker_aid[w] = None
            if case is None:
                continue
            if isinstance(outrow, dict):
                _append(outrow); done.add(case["key"])
                st = outrow.get("status", "ok")
                counts[st] = counts.get(st, 0) + 1
                if verbose and st == "ok":
                    print(f"[ok] {case['family']} p={case['shape_param']} "
                          f"SR~{case['size_ratio']} N={outrow.get('N')} "
                          f"phi={outrow.get('phi', float('nan')):.4f} "
                          f"Dmax/L={outrow.get('realized_dmaxL', float('nan')):.3f}",
                          flush=True)
            else:
                counts["failed"] += 1
            dispatch_idle()
    finally:
        pool.shutdown()
    if verbose:
        print(f"[sweep] done: {counts}", flush=True)
    return counts


def load_sweep_census(save_dir):
    """Load the baseline census into a pandas DataFrame."""
    import pandas as pd
    p = census_path(save_dir)
    if not p.exists():
        return pd.DataFrame()
    rows = []
    with open(p) as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return pd.DataFrame(rows)
