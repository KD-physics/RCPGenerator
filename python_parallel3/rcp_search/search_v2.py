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
    from .jobs import _MemTrace
    theta = np.asarray(payload["theta"], dtype=float)
    cid, seed = int(payload["candidate"]), int(payload["seed"])
    N, Ndim = int(pack_cfg["N"]), int(pack_cfg["Ndim"])
    mt = _MemTrace(f"c{cid}s{seed}")

    D = sample_diameters(theta, MODEL, N=N)
    rng = np.random.default_rng(cid * 100003 + seed)
    rng.shuffle(D)
    mt.mark("sample_D")

    ps = pack_cfg.get("prescreen") or {}
    if ps.get("enabled", False):
        reject, dmax, L_min = _prescreen_reject(
            D, Ndim, ps.get("phi_target", 0.85), ps.get("threshold", 2.5))
        if reject:
            return {"success": True, "phi": 0.0, "raw_phi": 0.0,
                    "phi_corr": 0.0, "max_overlap": 0.0, "rejected": True,
                    "seed": seed, "steps": 0, "wall_s": 0.0}

    mem_msg = _check_worker_mem(D, Ndim, int(pack_cfg.get("neighbor_max", 0)))
    if mem_msg:
        return {"success": False, "error": mem_msg, "seed": seed}

    import time as _time
    _t0 = _time.perf_counter()
    packing = run_one_custom_packing(
        D, seed=seed, Ndim=Ndim,
        phi_init=float(pack_cfg["phi_init"]),
        initializer_phi=float(pack_cfg.get("initializer_phi", 0.11)),
        neighbor_max=int(pack_cfg.get("neighbor_max", 0)),
        size_ratio_hint=MODEL.get("size_ratio"),
        verbose=False)
    _wall_s = _time.perf_counter() - _t0
    _steps = int(packing.steps) if getattr(packing, "steps", None) is not None else None
    mt.mark("pack")
    phi = float(safe_compute_phi(packing))
    mt.mark("compute_phi")
    from .jobs import engine_max_overlap
    mo = engine_max_overlap(packing)
    pc = float(phi_corrected(phi, mo, Ndim))
    mt.mark("overlap+corr")
    score = pc if pack_cfg.get("score_by") == "phi_corr" else phi
    out = {"success": True, "phi": score, "raw_phi": phi, "phi_corr": pc,
           "max_overlap": float(mo), "rejected": False,
           "seed": seed, "steps": _steps, "wall_s": float(_wall_s)}
    del packing
    mt.mark("teardown")
    mt.flush()
    return out


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
        progress_every=int(config.get("progress_every", 25)),
        progress_seconds=float(config.get("progress_seconds", 120.0)),
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


def _center_from_confirmed(confirmed_rows, config, MODEL):
    """Rank-weighted recombination center from the top-K confirmed rows
    (each a dict carrying 'theta', already ordered best-first). Returns None
    if empty.

    The exp-rank weights depend ONLY on rank, not on the mean values, so the
    center is bit-identical whenever the ordered top-K theta list is unchanged
    — even if more elite seeds refine the means. That is exactly the
    carry-over validity condition the speculator relies on.
    """
    recomb = confirmed_rows[:int(config.get("elite_count", 2))]
    if not recomb:
        return None
    elite_thetas = np.array([np.asarray(r["theta"], dtype=float)
                             for r in recomb])
    weights = np.exp(-np.arange(len(elite_thetas)))
    weights = weights / np.sum(weights)
    return clip_theta(np.sum(elite_thetas * weights[:, None], axis=0), MODEL)


class _Branch:
    """One hypothesised gen-(X+1) population: a candidate final recombination
    ordering, the center it implies, the gen+1 thetas drawn from it, and the
    cheap packings computed for it so far. `prob` is the Thompson-estimated
    probability this ordering is the one gen X actually lands on."""
    __slots__ = ("key", "prob", "center", "thetas", "pending",
                 "inflight", "results", "dispatched")

    def __init__(self, key, prob, center, thetas, cheap):
        self.key = key                 # tuple of recomb cids, ranked
        self.prob = float(prob)
        self.center = np.asarray(center, dtype=float)
        self.thetas = thetas
        pop = len(thetas)
        self.pending = [(cid, s) for cid in range(pop) for s in range(cheap)]
        self.inflight = set()
        self.results = {}              # (cid, seed) -> out
        self.dispatched = 0            # handed out (for proportional alloc)


class _Speculator:
    """Multi-branch compute-ahead: packs gen-(X+1) cheap seeds on gen-X's
    idle TAIL workers, HEDGED across the most likely final recombination
    orderings so the work survives the center moving.

    Why multiple branches: gen X's center is a rank-weighted sum of the
    top-K confirmed thetas, so it isn't pinned down until the elite means
    stop reordering. A single provisional center (phase 1) is discarded
    whenever that ordering shifts. Here we Thompson-sample the contenders'
    noisy means to estimate P(each ordering is final), build one gen+1
    population per high-probability ordering, and spend idle CPUs across
    them PROPORTIONAL to P(win). Whatever ordering gen X lands on, a branch
    almost surely already computed its cheaps.

    Correctness is identical to phase 1 per branch: each branch draws its
    thetas from a CLONE of the post-proposal RNG state, so the branch whose
    center equals gen X's FINAL center yields thetas bit-identical to what
    gen X+1 will draw, and its cheap packings are reusable. Non-matching
    branches are discarded. Strictly leftover capacity: real work always
    preempts (requeue() returns a killed seed to its branch).
    """

    def __init__(self, config, MODEL, rng_state_next, sigma_next, cheap_seeds,
                 branch_seed, n_branches=4, n_samples=400):
        self.config = config
        self.MODEL = MODEL
        self._rng_state_next = rng_state_next
        self._sigma_next = float(sigma_next)
        self._cheap = int(cheap_seeds)
        self._branch_seed = int(branch_seed)
        self._n_branches = int(config.get("compute_ahead_branches", n_branches))
        self._n_samples = int(config.get("compute_ahead_samples", n_samples))
        self.scheduler = None          # set by QueueScheduler.__init__
        self._branches = None          # list[_Branch] or None until init
        self.n_completed = 0

    # ---- branch construction ---------------------------------------------

    def _enumerate_branches(self, contenders, sigma):
        """Thompson-sample final recombination orderings from the contenders'
        means (each ~ N(mean, (sigma/sqrt(n))^2)); return the top-N most
        probable ordered top-K cid tuples with their probabilities."""
        K = int(self.config.get("elite_count", 2))
        rng = np.random.default_rng(self._branch_seed)
        cids = [c[0] for c in contenders]
        means = np.array([c[2] for c in contenders], dtype=float)
        ns = np.array([max(1, c[3]) for c in contenders], dtype=float)
        se = (sigma if sigma > 0 else 1e-9) / np.sqrt(ns)
        counts = {}
        for _ in range(self._n_samples):
            draw = means + se * rng.standard_normal(len(means))
            order = np.argsort(-draw)[:K]
            key = tuple(cids[i] for i in order)
            counts[key] = counts.get(key, 0) + 1
        ranked = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)
        ranked = ranked[:self._n_branches]
        tot = sum(c for _, c in ranked)
        return [(key, c / tot) for key, c in ranked]

    def _ensure_init(self) -> bool:
        if self._branches is not None:
            return len(self._branches) > 0
        # Lock branches only once gen X's screening floor is met (recomb set
        # settled enough to estimate). Before that _next_work() can idle in
        # the screening-drain gap with a partial, misleading bucket.
        sch = self.scheduler
        if sch._screened_count() < sch.cfg.min_cheap:
            return False
        contenders, sigma = sch.contenders_for_speculation()
        if not contenders:
            return False
        theta_by_cid = {c[0]: c[1] for c in contenders}
        branch_specs = self._enumerate_branches(contenders, sigma)
        pop = int(self.config["population_size"])
        branches = []
        for key, prob in branch_specs:
            rows = [{"theta": theta_by_cid[cid]} for cid in key]
            center = _center_from_confirmed(rows, self.config, self.MODEL)
            if center is None:
                continue
            rng = np.random.default_rng()
            rng.bit_generator.state = self._rng_state_next   # same start state
            thetas = [perturb_theta(center, self._sigma_next, self.MODEL, rng)
                      for _ in range(pop)]
            branches.append(_Branch(key, prob, center, thetas, self._cheap))
        self._branches = branches
        return len(branches) > 0

    # ---- work dispatch (proportional to branch probability) --------------

    def next_payload(self):
        if not self._ensure_init():
            return None
        # pick the branch with pending work that is most UNDER its proportional
        # share: minimise dispatched/prob (a branch with prob 0.5 should get
        # ~2x the seeds of a prob 0.25 branch).
        best, best_key = None, None
        for bi, b in enumerate(self._branches):
            if not b.pending:
                continue
            score = b.dispatched / b.prob if b.prob > 0 else float("inf")
            if best is None or score < best_key:
                best, best_key, best_i = b, score, bi
        if best is None:
            return None
        cid, seed = best.pending.pop(0)
        best.inflight.add((cid, seed))
        best.dispatched += 1
        return {"theta": best.thetas[cid], "candidate": cid, "seed": seed,
                "stage": "cheap", "branch": best_i}

    def on_result(self, payload, out):
        bi = int(payload.get("branch", 0))
        if bi >= len(self._branches):
            return
        b = self._branches[bi]
        key = (int(payload["candidate"]), int(payload["seed"]))
        b.inflight.discard(key)
        if out.get("success", True):
            b.results[key] = out
            self.n_completed += 1

    def requeue(self, payload):
        bi = int(payload.get("branch", 0))
        if bi >= len(self._branches):
            return
        b = self._branches[bi]
        key = (int(payload["candidate"]), int(payload["seed"]))
        if key in b.inflight:
            b.inflight.discard(key)
            b.pending.append(key)
            b.dispatched -= 1

    # ---- harvest ---------------------------------------------------------

    def harvest(self, real_center):
        """Return {'center','thetas','results'} for whichever branch's center
        matches gen X's FINAL center (at most one — the center is a
        deterministic function of the ordering); else None."""
        if not self._branches:
            return None
        rc = np.asarray(real_center, dtype=float)
        for b in self._branches:
            if b.results and np.allclose(b.center, rc, rtol=0.0, atol=1e-12):
                return {"center": b.center, "thetas": b.thetas,
                        "results": dict(b.results)}
        return None

    def branch_summary(self):
        if not self._branches:
            return "0 branches"
        return " ".join(f"b{i}(p={b.prob:.2f},n={len(b.results)})"
                        for i, b in enumerate(self._branches))


def _run_queue_generation(state, MODEL, pool, eval_fn, restore_snap=None):
    config = state["config"]
    gen = int(state["next_generation"])
    rng = np.random.default_rng()
    rng.bit_generator.state = state["rng_state"]
    center = np.asarray(state["center"], dtype=float)
    sigma = float(state["sigma"])
    pop = int(config["population_size"])

    restore_state = None
    if restore_snap is not None and int(restore_snap.get("generation", -1)) == gen:
        # Warm-start this generation from the realtime snapshot: reuse the
        # exact population thetas and re-seed completed results/states.
        sigma = float(restore_snap.get("sigma", sigma))
        thetas = [np.asarray(t, dtype=float) for t in restore_snap["thetas"]]
        restore_state = restore_snap.get("candidates")
        print(f"\n=== Generation {gen + 1}/{config['generations']} "
              f"| sigma={sigma:.4f} | mode=queue-v2 | RESTORED ===")
        done = sum(len(v.get("results", [])) for v in restore_state.values())
        bucket = sum(1 for v in restore_state.values()
                     if v.get("state") in ("elite", "confirmed"))
        leadvals = [sum(v["results"]) / len(v["results"])
                    for v in restore_state.values() if v.get("results")]
        best = max(leadvals) if leadvals else float("nan")
        print(f"[progress] RESTORED: completed_packings={done} "
              f"bucket={bucket} kills={int(restore_snap.get('kills', 0))} "
              f"best_so_far={best:.4f} | resuming unfinished + in-flight thetas",
              flush=True)
    else:
        print(f"\n=== Generation {gen + 1}/{config['generations']} "
              f"| sigma={sigma:.4f} | mode=queue-v2 ===")
        thetas = [perturb_theta(center, sigma, MODEL, rng) for _ in range(pop)]

    # ---- compute-ahead (opt-in: config['compute_ahead']) ----------------
    # Two halves: (a) CONSUME a cache the previous gen computed on its idle
    # tail (validated by re-deriving this gen's thetas and matching them);
    # (b) CREATE a speculator that draws gen+1's population from the
    # post-proposal RNG state and packs its cheaps on THIS gen's idle tail.
    precomputed = None
    speculator = None
    if config.get("compute_ahead", False) and restore_state is None:
        carry = state.pop("_carry", None)
        if carry is not None:
            cthetas = carry.get("thetas", [])
            pre = {}
            for (cid, seed), out in carry.get("results", {}).items():
                if cid < len(thetas) and cid < len(cthetas) and np.allclose(
                        np.asarray(thetas[cid], dtype=float),
                        np.asarray(cthetas[cid], dtype=float),
                        rtol=0.0, atol=1e-12):
                    pre[(cid, seed)] = out
            if pre:
                precomputed = pre
                print(f"[compute-ahead] carried {len(pre)} pre-screened "
                      f"packings into gen {gen + 1} (recompute skipped).",
                      flush=True)
            elif carry.get("results"):
                print("[compute-ahead] carry discarded (center moved); "
                      "gen recomputed from scratch.", flush=True)
        if gen + 1 < int(config["generations"]):
            speculator = _Speculator(
                config, MODEL,
                deepcopy(rng.bit_generator.state),
                sigma * float(config.get("sigma_decay", 0.80)),
                int(config.get("cheap_seeds_per_candidate", 2)),
                branch_seed=int(config.get("master_seed", 1234)) * 100003
                + gen * 7 + 99)
    else:
        state.pop("_carry", None)      # never let a stale carry linger

    # Realtime snapshot: written after each result+dispatch (worker-finish
    # cadence), atomic temp+rename. Only loss on a bump is in-flight compute.
    def _snap(sched_self):
        try:
            from .snapshot import write_snapshot
            write_snapshot(config["save_dir"], generation=gen, sigma=sigma,
                           thetas=thetas, candidates=sched_self.cands,
                           kills=sched_self.kills, config=config, model=MODEL)
        except Exception as _se:
            print(f"[snapshot] skipped (non-fatal): {_se!r}", flush=True)

    sched = QueueScheduler(_qs_config_from(config), thetas, eval_fn,
                           pool=pool,
                           log=bool(config.get("verbose_scheduler", False)),
                           on_tick=(_snap if config.get("snapshot", True)
                                    else None),
                           restore_state=restore_state,
                           speculator=speculator, precomputed=precomputed)
    res = sched.run()

    # Generation finished cleanly -> the checkpoint below supersedes the
    # snapshot; clear it so a later restore=True never grabs a done gen.
    try:
        from .snapshot import clear_snapshot
        clear_snapshot(config["save_dir"])
    except Exception:
        pass

    # Cycle 22 packing census: append one row per packing (parent-side,
    # single writer; reads results already collected — touches no worker).
    try:
        from .census import append_generation as _census_append
        _census_append(state, MODEL, res, config)
    except Exception as _census_exc:   # logging must never break a search
        print(f"[census] skipped (non-fatal): {_census_exc!r}", flush=True)

    rows = [_cand_row(c, gen) for c in res["candidates"] if c.n > 0]
    rows.sort(key=lambda r: r["score"], reverse=True)
    confirmed = [r for r in rows
                 if r["state"] in ("confirmed", "elite")]
    print(f"[gen {gen + 1}] seeds={res['total_seeds']} kills={res['kills']} "
          f"utilization={res['utilization']:.2f} "
          f"sigma_pooled={res['pooled_sigma']:.2e} "
          f"tail_idle={res.get('tail_idle_cpu_s', 0.0)/60.0:.1f} cpu-min")
    # Top-candidate distribution table (default on; set
    # config['print_top_distribution']=False to silence). Shows how the
    # winning distribution shape evolves generation to generation.
    if config.get("print_top_distribution", True) and rows:
        try:
            from .model import theta_to_dataframe
            # Show the top CONFIRMED candidate (the one that actually feeds
            # the decision), not the global max-mean — a 1-2 seed noise draw
            # can top the raw list without being trustworthy or selected.
            top = confirmed[0] if confirmed else rows[0]
            tag = "confirmed elite" if confirmed else "UNCONFIRMED (no elite this gen)"
            tdf = theta_to_dataframe(top["theta"], MODEL)
            print(f"Top candidate c{top['candidate']:02d} "
                  f"(mean_phi={top['mean_phi']:.6f}, n={top['n']}, "
                  f"state={top['state']}, {tag}) distribution:")
            print(tdf.to_string(index=False,
                                float_format=lambda x: f"{x:.6g}"))
        except Exception as _disp_exc:
            print(f"[top-dist] skipped (non-fatal): {_disp_exc!r}", flush=True)

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
        # center update: recombination set = top-K CONFIRMED candidates by
        # mean (NOT all rows — including dropped/demoted candidates would
        # let a single-seed noise outlier with a lucky high mean dominate
        # the center via the rank-0 weight, steering the search on noise).
        # Matches v1's batch path, which recombines from confirmed only.
        new_center = _center_from_confirmed(confirmed, config, MODEL)
        if new_center is not None:
            state["center"] = new_center

    # compute-ahead harvest: if the center held, carry the speculator's
    # idle-tail packings into gen+1; otherwise discard. Done AFTER the center
    # update so the validity check sees this gen's FINAL center.
    if speculator is not None:
        carry = speculator.harvest(state["center"])
        msg = (f"[compute-ahead] gen {gen + 1} idle tail: "
               f"computed={res.get('spec_completed', 0)} "
               f"preempted={res.get('spec_preempted', 0)} "
               f"| branches: {speculator.branch_summary()}")
        if carry is not None:
            state["_carry"] = carry
            msg += f" -> carrying {len(carry['results'])} into gen {gen + 2}"
        elif res.get("spec_completed", 0):
            msg += " -> discarded (no branch matched final center)"
        print(msg, flush=True)

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

def run_or_resume_search_v2(theta_start, config, MODEL, force_restart=False,
                            restore=False, eval_fn=None, return_restored=False):
    """v2 search entry point — same contract and checkpoint as
    run_or_resume_search; evaluation runs on the no-phases queue.

    restore=False (default): unchanged behavior — resume from the last
      COMPLETED generation in checkpoint.pkl, honoring the passed CONFIG/
      MODEL (so you can change generations/N/etc and continue).
    restore=True: true mid-generation restore from the realtime snapshot —
      reverts CONFIG/MODEL to the snapshot's, warm-starts the in-progress
      generation (only in-flight packings are lost), then continues.
    eval_fn: optional override of the per-packing evaluator (module-level,
      picklable). Default None -> the real rcpgenerator path. Used by tests.
    return_restored=True: return (theta_best, state, config, model) instead
      of (theta_best, state). The restored CONFIG/MODEL are also always on
      state['config'] / state['model'].
    """
    config = deepcopy(config)

    # restore=True: authoritatively revert CONFIG and MODEL to the snapshot
    # of the in-progress generation (the completed packings were computed
    # under them; a changed N etc. would be inconsistent).
    restore_snap = None
    if restore:
        from .snapshot import load_snapshot
        restore_snap = load_snapshot(config.get("save_dir"))
        if restore_snap is None:
            print("[restore] no snapshot found; falling back to normal "
                  "resume from the last completed generation.", flush=True)
        else:
            snap_cfg = deepcopy(restore_snap["config"])
            snap_cfg["save_dir"] = config.get("save_dir", snap_cfg.get("save_dir"))
            diffs = [k for k in snap_cfg
                     if k in config and config[k] != snap_cfg[k]]
            if diffs:
                print(f"[restore] using snapshot CONFIG/MODEL; ignoring "
                      f"passed-in differences in: {sorted(diffs)}", flush=True)
            config = snap_cfg
            MODEL = restore_snap["model"]

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
    if eval_fn is None:
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
        first = True
        while int(state["next_generation"]) < int(config["generations"]):
            snap_for_gen = restore_snap if first else None
            first = False
            state = _run_queue_generation(state, MODEL, pool, eval_fn,
                                          restore_snap=snap_for_gen)
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
    state["config"] = config
    state["model"] = MODEL
    if return_restored:
        return theta_best, state, config, MODEL
    return theta_best, state
