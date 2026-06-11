"""Evolutionary search loop over diameter-distribution parameters.

Single entry point: :func:`run_or_resume_search`. One CONFIG dict carries
every knob (packing, search, scheduler, pre-screen, time-aware, progress,
packing-save). MODEL is the parametric distribution from
:mod:`rcpgenerator.search.model`.

Two scheduler modes:

* ``"batch"`` (default) — cheap screen → barrier → elite confirm → barrier.
* ``"adaptive"`` — :class:`~rcpgenerator.search.scheduler.AdaptiveScheduler`
  with budget recycling on cheap drops / confirm drops. Activated by
  ``CONFIG['scheduler_mode'] == 'adaptive'``.
"""
from __future__ import annotations

import os
import traceback
from copy import deepcopy

import numpy as np

from ._helpers import get_final_phi, safe_compute_phi
from .jobs import (
    aggregate_candidate_results,
    make_packing_job,
    resolve_n_workers,
    run_jobs,
    run_one_custom_packing,
)
from .model import clip_theta, perturb_theta, sample_diameters
from .overlap import overlap_report, phi_corrected
from .persistence import (
    checkpoint_path,
    delete_checkpoint,
    load_checkpoint_dict,
    result_to_log_row,
    save_best_result,
    save_checkpoint,
    save_config_sidecar,
    save_model_sidecar,
    save_search_log,
    search_dir,
)


# ============================================================================
# Population sampling.
# ============================================================================
def random_theta(MODEL, rng=None):
    from .model import init_theta_random
    return init_theta_random(MODEL, rng=rng)


def propose_population(center, sigma, population_size, MODEL, rng=None):
    """Generate `population_size` candidates by Gaussian perturbation."""
    rng = np.random.default_rng(rng) if not isinstance(rng, np.random.Generator) else rng
    return [perturb_theta(center, sigma, MODEL, rng) for _ in range(population_size)]


# ============================================================================
# State construction.
# ============================================================================
def initialize_search_state(theta_start, config, MODEL):
    rng = np.random.default_rng(int(config.get("master_seed", 1234)))
    center = clip_theta(np.asarray(theta_start, dtype=float), MODEL)
    return {
        "config": deepcopy(config),
        "model_name": MODEL.get("name", "unnamed"),
        "center": center,
        "sigma": float(config["initial_sigma"]),
        "next_generation": 0,
        "best_result": None,
        "log_rows": [],
        "rng_state": rng.bit_generator.state,
    }


def print_config_summary(config, MODEL):
    print("\nSearch config:")
    keys = ["search_name", "Ndim", "N", "phi_init",
            "neighbor_max", "generations", "population_size",
            "cheap_seeds_per_candidate", "elite_count", "elite_extra_seeds",
            "initial_sigma", "sigma_decay", "scheduler_mode", "score_by",
            "save_packings", "use_parallel", "n_workers", "reserve_cpus",
            "save_dir"]
    for k in keys:
        if k in config:
            print(f"  {k}: {config[k]}")
    print(f"\nMODEL: {MODEL.get('name','unnamed')}  "
          f"(family={MODEL.get('family','?')}, theta_dim={MODEL['theta_dim']}, "
          f"size_ratio={MODEL.get('size_ratio')})")
    print("Detected CPUs:", os.cpu_count())
    print("Resolved workers:", resolve_n_workers(config))


# ============================================================================
# Batch-mode generation.
# ============================================================================
def _run_batch_generation(state, MODEL):
    config = state["config"]
    gen = int(state["next_generation"])
    center = np.asarray(state["center"], dtype=float)
    sigma = float(state["sigma"])
    score_by = str(config.get("score_by", "mean_phi"))
    save_mode = str(config.get("save_packings", "none")).lower()

    rng = np.random.default_rng()
    rng.bit_generator.state = state["rng_state"]

    print(f"\n=== Generation {gen + 1}/{config['generations']} "
          f"| sigma={sigma:.4f} | mode=batch | score_by={score_by} ===")

    candidates = [center] + propose_population(
        center, sigma, int(config["population_size"]) - 1, MODEL, rng)
    theta_by_candidate = {j: np.asarray(t, dtype=float) for j, t in enumerate(candidates)}

    N_packing = int(config["N"])
    diameters_by_candidate = {
        j: sample_diameters(t, MODEL, N=N_packing)
        for j, t in theta_by_candidate.items()
    }

    cheap_seeds = [int(100000 + 1000 * gen + i)
                   for i in range(int(config["cheap_seeds_per_candidate"]))]
    cheap_jobs = []
    for j, theta in theta_by_candidate.items():
        D = diameters_by_candidate[j]
        for seed in cheap_seeds:
            cheap_jobs.append(make_packing_job(
                theta, MODEL, j, seed, config,
                diameters=D, generation=gen, stage="cheap"))

    print(f"Cheap screen: {len(cheap_jobs)} packing jobs")
    cheap_results = run_jobs(cheap_jobs, config)
    cheap_agg = aggregate_candidate_results(
        cheap_results, theta_by_candidate, score_by=score_by)

    for r in cheap_agg:
        state["log_rows"].append(
            result_to_log_row(r, gen, "cheap", sigma, config, MODEL))

    print("Cheap screen top candidates:")
    for r in cheap_agg[:min(5, len(cheap_agg))]:
        print(f"  cand {r['candidate']:02d}: "
              f"mean_phi={r['mean_phi']:.6f}, "
              f"mean_phi_corr={r['mean_phi_corr']:.6f}, "
              f"n={len(r['phis'])}, failures={r['failures']}, "
              f"rejects={r['rejects']}")

    elite_count = int(config["elite_count"])
    elite_extra = int(config["elite_extra_seeds"])
    elites = cheap_agg[:elite_count]
    elite_theta_by_candidate = {int(r["candidate"]): np.asarray(r["theta"], dtype=float)
                                for r in elites}

    elite_jobs = []
    for rank, r in enumerate(elites):
        j = int(r["candidate"])
        theta = np.asarray(r["theta"], dtype=float)
        D = diameters_by_candidate[j]
        extra_seeds = [int(200000 + 1000 * gen + 100 * rank + i)
                       for i in range(elite_extra)]
        for seed in extra_seeds:
            elite_jobs.append(make_packing_job(
                theta, MODEL, j, seed, config,
                diameters=D, generation=gen, stage="elite_extra"))

    print(f"Elite confirmation: {len(elite_jobs)} extra packing jobs")
    elite_results = run_jobs(elite_jobs, config)
    cheap_for_elites = [row for row in cheap_results
                        if int(row["candidate"]) in elite_theta_by_candidate]
    confirmed_results = aggregate_candidate_results(
        cheap_for_elites + elite_results, elite_theta_by_candidate,
        score_by=score_by)

    for r in confirmed_results:
        state["log_rows"].append(
            result_to_log_row(r, gen, "confirmed", sigma, config, MODEL))

    print("Confirmed elites:")
    for r in confirmed_results:
        print(f"  cand {r['candidate']:02d}: "
              f"mean_phi={r['mean_phi']:.6f}, "
              f"mean_phi_corr={r['mean_phi_corr']:.6f}, "
              f"std={r['std_phi']:.6f}, "
              f"n={len(r['phis'])}")

    # Save the per-generation best packing if requested.
    if save_mode == "best" and confirmed_results:
        try:
            _save_best_per_gen(confirmed_results[0], MODEL, config, gen)
        except Exception as exc:
            print("Warning: best-per-gen save failed:", exc)

    if confirmed_results:
        best_this_gen = confirmed_results[0]
        cur_best = state["best_result"]
        best_key = "mean_phi_corr" if score_by == "phi_corr" else "mean_phi"
        if cur_best is None or best_this_gen[best_key] > cur_best.get(best_key, -np.inf):
            state["best_result"] = deepcopy(best_this_gen)
        elite_thetas = np.array(
            [np.asarray(r["theta"], dtype=float) for r in confirmed_results])
        ranks = np.arange(len(elite_thetas))
        weights = np.exp(-ranks); weights = weights / np.sum(weights)
        new_center = np.sum(elite_thetas * weights[:, None], axis=0)
        state["center"] = clip_theta(new_center, MODEL)

    state["sigma"] = sigma * float(config.get("sigma_decay", 0.80))
    state["next_generation"] = gen + 1
    state["rng_state"] = rng.bit_generator.state

    save_search_log(state["log_rows"], config)
    save_best_result(state["best_result"], config)
    ckpt = save_checkpoint(state, config)
    print("Checkpoint saved:", ckpt)
    if state["best_result"] is not None:
        print(f"Best overall mean_phi: {state['best_result']['mean_phi']:.6f}  "
              f"phi_corr: {state['best_result']['mean_phi_corr']:.6f}")
    return state


def _save_best_per_gen(best_candidate, MODEL, config, generation):
    """Re-run the candidate's seed-0 packing and dump it as a .npz."""
    theta = np.asarray(best_candidate["theta"], dtype=float)
    D = sample_diameters(theta, MODEL, N=int(config["N"]))
    rng = np.random.default_rng(0)
    rng.shuffle(D)
    packing = run_one_custom_packing(
        D, seed=0, Ndim=int(config["Ndim"]),
        phi_init=float(config["phi_init"]),
        initializer_phi=float(config.get("initializer_phi", 0.11)),
        neighbor_max=int(config.get("neighbor_max", 0)),
        size_ratio_hint=MODEL.get("size_ratio"),
        verbose=False,
    )
    import pathlib
    out_dir = pathlib.Path(config["save_dir"]) / "packings"
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"gen{generation:03d}_best_seed0.npz"
    np.savez(
        path,
        positions=np.asarray(packing.positions, dtype=float),
        diameters=np.asarray(packing.diameters, dtype=float),
        box=np.asarray(packing.box, dtype=float),
        theta=theta,
    )
    print(f"  saved best-per-gen packing: {path}")


# ============================================================================
# Adaptive-mode generation.
# ============================================================================
def _adaptive_packing_worker(payload):
    """Top-level (picklable) packing worker for the AdaptiveScheduler path."""
    try:
        from .jobs import _pin_worker_threads, _check_worker_mem
        _pin_worker_threads()
        _mem_msg = _check_worker_mem(
            payload["diameters"], int(payload["Ndim"]),
            int(payload.get("neighbor_max", 0)))
        if _mem_msg:
            return {"success": False, "error": _mem_msg,
                    "stage": payload.get("stage", "unknown"),
                    "seed_id": int(payload.get("seed_id", -1))}
        D = np.asarray(payload["diameters"], dtype=float)
        seed_id = int(payload["seed_id"])
        Ndim = int(payload["Ndim"])
        rng = np.random.default_rng(seed_id)
        rng.shuffle(D)
        packing = run_one_custom_packing(
            D, seed=seed_id, Ndim=Ndim,
            phi_init=float(payload["phi_init"]),
            initializer_phi=float(payload["initializer_phi"]),
            neighbor_max=int(payload.get("neighbor_max", 0)),
            size_ratio_hint=payload.get("size_ratio_hint"),
            verbose=False,
        )
        phi = safe_compute_phi(packing)
        from .jobs import engine_max_overlap
        mo = engine_max_overlap(packing)
        phi_corr = phi_corrected(phi, mo, Ndim)
        return {
            "success": True,
            "phi": float(phi),
            "phi_corr": float(phi_corr),
            "max_overlap": float(mo),
            "stage": payload.get("stage", "unknown"),
            "seed_id": seed_id,
        }
    except Exception as e:
        return {
            "success": False,
            "error": repr(e),
            "traceback": traceback.format_exc(),
            "stage": payload.get("stage", "unknown"),
            "seed_id": int(payload.get("seed_id", -1)),
        }


def _make_adaptive_propose_theta(state, MODEL, rng):
    center = np.asarray(state["center"], dtype=float)
    sigma = float(state["sigma"])
    def propose():
        return perturb_theta(center, sigma, MODEL, rng)
    return propose


def _candidate_to_result_dict(candidate, MODEL):
    phis = list(candidate.all_phi)
    return {
        "candidate": int(candidate.candidate_id),
        "theta": np.asarray(candidate.theta, dtype=float),
        "mean_phi": float(candidate.mean_phi) if candidate.mean_phi is not None else -np.inf,
        "median_phi": float(np.median(phis)) if phis else -np.inf,
        "best_phi": float(candidate.best_phi) if candidate.best_phi is not None else -np.inf,
        "std_phi": float(candidate.std_phi) if candidate.std_phi is not None else np.inf,
        "mean_phi_corr": float(candidate.mean_phi) if candidate.mean_phi is not None else -np.inf,
        "mean_max_overlap": 0.0,
        "phis": phis,
        "phi_corrs": [], "max_overlaps": [],
        "steps": [], "forces": [], "max_min_dists": [], "seeds": [],
        "failures": int(candidate.n_failed),
        "rejects": 0,
        "raw_rows": [],
        "mean_force": float("inf"),
        "mean_steps": float("inf"),
        "mean_max_min_dist": float("inf"),
    }


def _run_adaptive_generation(state, MODEL):
    from .scheduler import (
        AdaptiveScheduler,
        CandidateState,
        LogSink,
        RealRunner,
        SchedulerConfig,
    )

    if state["config"].get("time_aware_enabled", False) or \
            state["config"].get("progress_enabled", False):
        # Apply time-aware + progress monkey-patches if enabled.
        try:
            from . import time_aware  # noqa: F401  (import triggers patching)
            time_aware.apply(state["config"])
        except Exception as exc:
            print("Warning: time-aware patch not applied:", exc)

    config = state["config"]
    gen = int(state["next_generation"])
    sigma = float(state["sigma"])
    center = np.asarray(state["center"], dtype=float)
    score_by = str(config.get("score_by", "mean_phi"))
    save_mode = str(config.get("save_packings", "none")).lower()

    rng = np.random.default_rng()
    rng.bit_generator.state = state["rng_state"]

    print(f"\n=== Generation {gen + 1}/{config['generations']} "
          f"| sigma={sigma:.4f} | mode=adaptive | score_by={score_by} ===")

    initial_thetas = [center] + propose_population(
        center, sigma, int(config["population_size"]) - 1, MODEL, rng)

    sched_cfg = SchedulerConfig(
        population_size=int(config["population_size"]),
        cheap_seeds_per_candidate=int(config["cheap_seeds_per_candidate"]),
        elite_count=int(config["elite_count"]),
        elite_extra_seeds=int(config["elite_extra_seeds"]),
        cheap_delta=float(config.get("cheap_delta", 0.005)),
        phi_noise_floor=float(config.get("phi_noise_floor", 0.001)),
        n_workers=int(resolve_n_workers(config)),
        scheduler_seed=int(config.get("master_seed", 0)) + gen,
        throttle_local=bool(config.get("throttle_local", False)),
        worker_cooldown_seconds=float(config.get("worker_cooldown_seconds", 0.0)),
        real_smoke_max_wall_seconds=float(config.get("real_smoke_max_wall_seconds", 1e6)),
    )

    class _AdaptiveRealRunner(RealRunner):
        def __init__(self, *, n_workers, MODEL, packing_config,
                     throttle_local=False, worker_cooldown_seconds=0.0):
            super().__init__(
                n_workers=n_workers,
                packing_callable=_adaptive_packing_worker,
                throttle_local=throttle_local,
                worker_cooldown_seconds=worker_cooldown_seconds,
            )
            self._adaptive_MODEL = MODEL
            self._adaptive_N = int(packing_config["N"])
            self._adaptive_Ndim = int(packing_config["Ndim"])
            self._adaptive_phi_init = float(packing_config["phi_init"])
            self._adaptive_init_phi = float(packing_config.get("initializer_phi", 0.11))
            self._adaptive_neighbor_max = int(packing_config.get("neighbor_max", 0))
            self._adaptive_sr_hint = MODEL.get("size_ratio")

        def dispatch(self, attempt):
            w = self._pick_idle_worker()
            assert w is not None, "_AdaptiveRealRunner.dispatch: no idle worker"
            attempt.worker_id = w
            attempt.started_at = self.now()
            theta_arr = np.asarray(attempt.theta, dtype=float)
            D = sample_diameters(theta_arr, self._adaptive_MODEL,
                                 N=self._adaptive_N)
            stage_val = attempt.stage.value if hasattr(attempt.stage, "value") else str(attempt.stage)
            payload = {
                "diameters": D.tolist(),
                "seed_id": int(attempt.seed_id),
                "Ndim": self._adaptive_Ndim,
                "phi_init": self._adaptive_phi_init,
                "initializer_phi": self._adaptive_init_phi,
                "neighbor_max": self._adaptive_neighbor_max,
                "size_ratio_hint": self._adaptive_sr_hint,
                "stage": stage_val,
            }
            import time as _t
            start_wall = _t.monotonic()
            f = self.executor.submit(_adaptive_packing_worker, payload)
            self._future_to_data[f] = (attempt, start_wall, w)

    runner = _AdaptiveRealRunner(
        n_workers=sched_cfg.n_workers,
        MODEL=MODEL,
        packing_config=config,
        throttle_local=sched_cfg.throttle_local,
        worker_cooldown_seconds=sched_cfg.worker_cooldown_seconds,
    )
    sink_dir = search_dir(config) / f"adaptive_gen_{gen+1:03d}"
    sink_dir.mkdir(parents=True, exist_ok=True)
    sink = LogSink(sink_dir)
    propose_next = _make_adaptive_propose_theta(state, MODEL, rng)
    scheduler = AdaptiveScheduler(sched_cfg, runner, sink, propose_next)

    prior_best = (state["best_result"]["mean_phi"]
                  if state.get("best_result") is not None else float("-inf"))

    try:
        summary = scheduler.run_generation(initial_thetas, prior_best=prior_best)
    finally:
        runner.shutdown()

    sink.write_attempt_log()
    sink.write_dispatch_log()
    sink.write_drop_log()
    sink.write_budget_log()
    sink.write_candidate_log(scheduler.candidates, sched_cfg)
    sink.write_summary(summary)

    print(f"Adaptive scheduler stats:")
    print(f"  total packings dispatched: {summary['total_packings']}")
    print(f"  cheap dispatched: {summary['n_cheap_dispatched_total']}, "
          f"confirm dispatched: {summary['n_confirm_dispatched_total']}")
    print(f"  cheap drops: {summary['n_dropped_cheap']}, "
          f"confirm drops: {summary['n_dropped_confirm']}")
    print(f"  perturbations proposed: {summary['n_perturbations_proposed']}")
    print(f"  eligibles (full quota): {summary['n_eligibles']}")

    eligible_candidates = sorted(
        [c for c in scheduler.candidates.values()
         if c.state == CandidateState.CONFIRMED],
        key=lambda c: -c.mean_phi,
    )
    for c in eligible_candidates:
        r = _candidate_to_result_dict(c, MODEL)
        state["log_rows"].append(
            result_to_log_row(r, gen, "confirmed", sigma, config, MODEL))

    print("Confirmed elites:")
    for c in eligible_candidates:
        std_str = f"{c.std_phi:.6f}" if c.std_phi is not None else "n/a"
        print(f"  cand {c.candidate_id:02d}: mean={c.mean_phi:.6f}, "
              f"best={c.best_phi:.6f}, std={std_str}, n={c.n_valid_total}")

    if eligible_candidates:
        best_this_gen = _candidate_to_result_dict(eligible_candidates[0], MODEL)
        cur_best = state["best_result"]
        if cur_best is None or best_this_gen["mean_phi"] > cur_best.get("mean_phi", -np.inf):
            state["best_result"] = deepcopy(best_this_gen)
        elite_thetas = np.array(
            [np.asarray(c.theta, dtype=float) for c in eligible_candidates])
        ranks = np.arange(len(elite_thetas))
        weights = np.exp(-ranks); weights = weights / np.sum(weights)
        new_center = np.sum(elite_thetas * weights[:, None], axis=0)
        state["center"] = clip_theta(new_center, MODEL)

        if save_mode == "best":
            try:
                _save_best_per_gen(best_this_gen, MODEL, config, gen)
            except Exception as exc:
                print("Warning: best-per-gen save failed:", exc)

    state["sigma"] = sigma * float(config.get("sigma_decay", 0.80))
    state["next_generation"] = gen + 1
    state["rng_state"] = rng.bit_generator.state

    save_search_log(state["log_rows"], config)
    save_best_result(state["best_result"], config)
    ckpt = save_checkpoint(state, config)
    print("Checkpoint saved:", ckpt)
    if state["best_result"] is not None:
        print(f"Best overall mean_phi: {state['best_result']['mean_phi']:.6f}")
    return state


# ============================================================================
# Top-level entry point.
# ============================================================================
def run_one_generation(state, MODEL):
    """Dispatch one generation based on `config['scheduler_mode']`."""
    if str(state["config"].get("scheduler_mode", "batch")).lower() == "adaptive":
        return _run_adaptive_generation(state, MODEL)
    return _run_batch_generation(state, MODEL)


def run_or_resume_search(theta_start, config, MODEL, force_restart=False):
    """Main search entry point.

    On first call: writes MODEL + CONFIG sidecars next to the checkpoint, so
    the directory becomes self-describing. On resume: skips sidecar writes.

    If ``config['speed']`` is set, the corresponding ``RCP_SPEED`` env var
    is applied for this run (and inherited by forked workers).

    Returns ``(theta_best, state)`` where ``state`` is the live dict (also
    available as ``state.SearchState`` via :class:`SearchState`).
    """
    config = deepcopy(config)
    search_dir(config)
    speed = config.get("speed")
    if speed is not None:
        os.environ["RCP_SPEED"] = str(speed)
        print(f"RCP_SPEED set to '{speed}' for this run.")
    print_config_summary(config, MODEL)

    if force_restart:
        delete_checkpoint(config)

    state = None
    if bool(config.get("resume_if_checkpoint_exists", True)):
        state = load_checkpoint_dict(config)
        if state is not None:
            ckpt_model = state.get("model_name", "unnamed")
            if ckpt_model != MODEL.get("name", "unnamed"):
                print(f"\nWARNING: checkpoint was for model '{ckpt_model}', "
                      f"resuming with '{MODEL.get('name','unnamed')}'.")
            print("\nLoaded checkpoint:", checkpoint_path(config))
            state["config"] = deepcopy(config)

    if state is None:
        print("\nStarting new search.")
        state = initialize_search_state(theta_start, config, MODEL)

    # Sidecar persistence — always re-write to capture in-session edits.
    save_config_sidecar(config)
    save_model_sidecar(MODEL, config)

    while int(state["next_generation"]) < int(config["generations"]):
        state = run_one_generation(state, MODEL)

    print("\nSearch complete.")
    theta_best = None
    if state.get("best_result") is not None:
        theta_best = np.asarray(state["best_result"]["theta"], dtype=float)
    return theta_best, state


def run_validation(theta, config, MODEL, seeds=None, label="validation"):
    """Run an extra-seed validation batch for one theta. Returns aggregated row."""
    import pandas as pd
    if seeds is None:
        seeds = list(range(10000, 10020))
    D = sample_diameters(np.asarray(theta, dtype=float), MODEL, N=int(config["N"]))
    theta_by_candidate = {0: np.asarray(theta, dtype=float)}
    jobs = [make_packing_job(theta, MODEL, 0, int(seed), config,
                             diameters=D, generation=None, stage=label)
            for seed in seeds]
    results = run_jobs(jobs, config)
    agg = aggregate_candidate_results(
        results, theta_by_candidate,
        score_by=str(config.get("score_by", "mean_phi")))[0]
    print(f"\n{label}:")
    print(f"  n_success: {len(agg['phis'])}, failures: {agg['failures']}, "
          f"rejects: {agg['rejects']}")
    print(f"  mean_phi: {agg['mean_phi']:.6f}  "
          f"mean_phi_corr: {agg['mean_phi_corr']:.6f}")
    print(f"  median_phi: {agg['median_phi']:.6f}  std_phi: {agg['std_phi']:.6f}")
    return agg, pd.DataFrame(results)
