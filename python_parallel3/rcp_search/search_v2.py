"""run_or_resume_search_v2 — the v1 search loop on the v2 queue scheduler.

Drop-in alternative entry point: same signature, same CONFIG/MODEL, same
checkpoint format and path as :func:`rcp_search.search.run_or_resume_search`
(the state dict — center, sigma, generation counter, best_result — is
scheduler-agnostic), so a campaign can switch between v1 and v2 by
changing one function name, in either direction, mid-campaign.

What changes is only HOW a generation is evaluated: the no-phases
preemptive queue (rcp_search.queue_scheduler) replaces the cheap->elite
barrier — min_cheap floor, UCB-gated evictable elite bucket, racing
demotions with real in-flight kills, fill-and-stop, ~full worker
utilization.

Extra CONFIG keys (all optional; sensible defaults):
  min_cheap          candidates screened before elite seeds flow
                     (default max(6, population_size // 4))
  min_elite          racing-immune strongest elites (default 2)
  gate_z             one-sided confidence for the UCB gate (default 2.0)
  gate_epsilon       indifference zone in phi (default 5e-4)
  stop_eps / stop_G  search-level stopping rule: stop when the best
                     improves by < stop_eps for stop_G consecutive
                     generations (defaults 1e-3 / 2; set stop_G=0 to
                     disable and always run config["generations"])
  stop_cheap_when_filled / kill_cheap_when_filled   (defaults True)

Scientific semantics preserved: the center update uses the TOP-K
candidates by mean (the recombination set, decoupled from confirmation
depth as agreed), with the same exp-rank weights and sigma decay as v1;
best_result tracking uses confirmed elites only.
"""
from __future__ import annotations

import os
from copy import deepcopy
from functools import partial

import numpy as np

from .model import clip_theta, perturb_theta, sample_diameters
from .persistence import (
    checkpoint_path, delete_checkpoint, load_checkpoint_dict,
    save_best_result, save_checkpoint, save_config_sidecar,
    save_model_sidecar, save_search_log, search_dir)
from .jobs import (
    _check_worker_mem, _pin_worker_threads, _prescreen_reject,
    resolve_n_workers, run_one_custom_packing)
from .overlap import overlap_report, phi_corrected
from ._helpers import safe_compute_phi
from .search import initialize_search_state, print_config_summary
from .queue_scheduler import PreemptivePool, QSConfig, QueueScheduler


# ============================================================================
# Worker-side evaluation (module-level: picklable into the pool)
# ============================================================================

def _v2_packing_eval(payload, MODEL=None, pack_cfg=None):
    """One packing seed for the queue scheduler. Returns {"phi": <score>,
    ...record fields}. "phi" carries the SCORING value (mean_phi or
    phi_corr per config) so the scheduler's statistics drive on it."""
    _pin_worker_threads()
    theta = np.asarray(payload["theta"], dtype=float)
    cid, seed = int(payload["candidate"]), int(payload["seed"])
    N, Ndim = int(pack_cfg["N"]), int(pack_cfg["Ndim"])

    D = sample_diameters(theta, MODEL, N=N)
    rng = np.random.default_rng(cid * 100003 + seed)
    rng.shuffle(D)

    ps = pack_cfg.get("prescreen") or {}
    if ps.get("enabled", False):
        reject, dmax, L_min = _prescreen_reject(
            D, Ndim, ps.get("phi_target", 0.85), ps.get("threshold", 2.5))
        if reject:
            return {"success": True, "phi": 0.0, "raw_phi": 0.0,
                    "phi_corr": 0.0, "max_overlap": 0.0, "rejected": True}

    mem_msg = _check_worker_mem(D, Ndim, int(pack_cfg.get("neighbor_max", 0)))
    if mem_msg:
        return {"success": False, "error": mem_msg}

    packing = run_one_custom_packing(
        D, seed=seed, Ndim=Ndim,
        phi_init=float(pack_cfg["phi_init"]),
        initializer_phi=float(pack_cfg.get("initializer_phi", 0.11)),
        neighbor_max=int(pack_cfg.get("neighbor_max", 0)),
        size_ratio_hint=MODEL.get("size_ratio"),
        verbose=False)
    phi = float(safe_compute_phi(packing))
    from .jobs import engine_max_overlap
    mo = engine_max_overlap(packing)
    pc = float(phi_corrected(phi, mo, Ndim))
    score = pc if pack_cfg.get("score_by") == "phi_corr" else phi
    return {"success": True, "phi": score, "raw_phi": phi, "phi_corr": pc,
            "max_overlap": float(mo), "rejected": False}


# ============================================================================
# One generation on the queue scheduler
# ============================================================================

def _qs_config_from(config) -> QSConfig:
    pop = int(config["population_size"])
    cheap = int(config.get("cheap_seeds_per_candidate", 2))
    max_elite = int(config.get("elite_count", 2))
    return QSConfig(
        n_workers=resolve_n_workers(config),
        cheap_seeds=cheap,
        min_cheap=int(config.get("min_cheap", max(6, pop // 4))),
        max_elite=max_elite,
        min_elite=int(config.get("min_elite", min(2, max_elite))),
        elite_seeds=cheap + int(config.get("elite_extra_seeds", 6)),
        z=float(config.get("gate_z", 2.0)),
        epsilon=float(config.get("gate_epsilon", 5e-4)),
        stop_cheap_when_filled=bool(config.get("stop_cheap_when_filled", True)),
        kill_cheap_when_filled=bool(config.get("kill_cheap_when_filled", True)),
    )


def _cand_row(c, gen):
    """v1-compatible result row from a v2 candidate."""
    phis = list(c.results)
    raw_phis = [r.get("raw_phi", r.get("phi")) for r in c.raw
                if r.get("success", True)]
    pcs = [r.get("phi_corr") for r in c.raw
           if r.get("success", True) and r.get("phi_corr") is not None]
    n = len(phis)
    mean = float(np.mean(phis)) if phis else 0.0
    return {
        "candidate": int(c.cid),
        "generation": int(gen),
        "theta": np.asarray(c.theta, dtype=float),
        "phis": phis,
        "mean_phi": float(np.mean(raw_phis)) if raw_phis else mean,
        "mean_phi_corr": float(np.mean(pcs)) if pcs else mean,
        "median_phi": float(np.median(phis)) if phis else 0.0,
        "std_phi": float(np.std(phis, ddof=1)) if n >= 2 else 0.0,
        "score": mean,
        "n": n,
        "state": c.state.value,
    }


def _run_queue_generation(state, MODEL, pool, eval_fn):
    config = state["config"]
    gen = int(state["next_generation"])
    rng = np.random.default_rng()
    rng.bit_generator.state = state["rng_state"]
    center = np.asarray(state["center"], dtype=float)
    sigma = float(state["sigma"])
    pop = int(config["population_size"])

    print(f"\n=== Generation {gen + 1}/{config['generations']} "
          f"| sigma={sigma:.4f} | mode=queue-v2 ===")

    thetas = [perturb_theta(center, sigma, MODEL, rng) for _ in range(pop)]
    sched = QueueScheduler(_qs_config_from(config), thetas, eval_fn,
                           pool=pool, log=bool(config.get("verbose_scheduler",
                                                          False)))
    res = sched.run()

    rows = [_cand_row(c, gen) for c in res["candidates"] if c.n > 0]
    rows.sort(key=lambda r: r["score"], reverse=True)
    confirmed = [r for r in rows
                 if r["state"] in ("confirmed", "elite")]
    print(f"[gen {gen + 1}] seeds={res['total_seeds']} kills={res['kills']} "
          f"utilization={res['utilization']:.2f} "
          f"sigma_pooled={res['pooled_sigma']:.2e}")
    print("Elites:")
    for r in confirmed:
        print(f"  cand {r['candidate']:02d}: mean_phi={r['mean_phi']:.6f} "
              f"phi_corr={r['mean_phi_corr']:.6f} std={r['std_phi']:.6f} "
              f"n={r['n']} ({r['state']})")

    score_by = str(config.get("score_by", "mean_phi"))
    best_key = "mean_phi_corr" if score_by == "phi_corr" else "mean_phi"
    if confirmed:
        best_this_gen = confirmed[0]
        cur = state["best_result"]
        if cur is None or best_this_gen[best_key] > cur.get(best_key, -np.inf):
            state["best_result"] = deepcopy(best_this_gen)
        # center update: TOP-K BY MEAN (recombination set, decoupled from
        # confirmation depth) with v1's exp-rank weights
        recomb = rows[:int(config.get("elite_count", 2))]
        elite_thetas = np.array([np.asarray(r["theta"], dtype=float)
                                 for r in recomb])
        ranks = np.arange(len(elite_thetas))
        weights = np.exp(-ranks)
        weights = weights / np.sum(weights)
        state["center"] = clip_theta(
            np.sum(elite_thetas * weights[:, None], axis=0), MODEL)

    state["log_rows"].extend(
        {k: (v.tolist() if isinstance(v, np.ndarray) else v)
         for k, v in r.items() if k != "phis"} for r in rows)
    state["sigma"] = sigma * float(config.get("sigma_decay", 0.80))
    state["next_generation"] = gen + 1
    state["rng_state"] = rng.bit_generator.state

    save_search_log(state["log_rows"], config)
    save_best_result(state["best_result"], config)
    ckpt = save_checkpoint(state, config)
    print("Checkpoint saved:", ckpt)
    if state["best_result"] is not None:
        print(f"Best overall {best_key}: "
              f"{state['best_result'][best_key]:.6f}")
    return state


# ============================================================================
# Entry point
# ============================================================================

def run_or_resume_search_v2(theta_start, config, MODEL, force_restart=False):
    """v2 search entry point — same contract and checkpoint as
    run_or_resume_search; evaluation runs on the no-phases queue."""
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
        print("\nStarting new search (queue-v2).")
        state = initialize_search_state(theta_start, config, MODEL)

    save_config_sidecar(config)
    save_model_sidecar(MODEL, config)

    pack_cfg = {
        "N": int(config["N"]), "Ndim": int(config["Ndim"]),
        "phi_init": float(config["phi_init"]),
        "initializer_phi": float(config.get("initializer_phi", 0.11)),
        "neighbor_max": int(config.get("neighbor_max", 0)),
        "prescreen": config.get("prescreen"),
        "score_by": str(config.get("score_by", "mean_phi")),
    }
    eval_fn = partial(_v2_packing_eval, MODEL=MODEL, pack_cfg=pack_cfg)
    pool = PreemptivePool(resolve_n_workers(config), eval_fn)

    stop_eps = float(config.get("stop_eps", 1e-3))
    stop_G = int(config.get("stop_G", 2))
    best_key = ("mean_phi_corr"
                if str(config.get("score_by", "mean_phi")) == "phi_corr"
                else "mean_phi")
    flat = 0
    prev_best = (state["best_result"][best_key]
                 if state.get("best_result") else -np.inf)
    try:
        while int(state["next_generation"]) < int(config["generations"]):
            state = _run_queue_generation(state, MODEL, pool, eval_fn)
            new_best = (state["best_result"][best_key]
                        if state.get("best_result") else -np.inf)
            if stop_G > 0:
                if new_best - prev_best < stop_eps:
                    flat += 1
                    if flat >= stop_G:
                        print(f"\nStopping rule: best improved < {stop_eps} "
                              f"for {stop_G} consecutive generations.")
                        break
                else:
                    flat = 0
                prev_best = max(prev_best, new_best)
    finally:
        pool.shutdown()

    print("\nSearch complete.")
    theta_best = None
    if state.get("best_result") is not None:
        theta_best = np.asarray(state["best_result"]["theta"], dtype=float)
    return theta_best, state
