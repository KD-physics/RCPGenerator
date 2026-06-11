"""v2 "no-phases" preemptive scheduler for rcp_search.

One pool of workers, one priority queue of (candidate, seed) work items.
"Cheap screening" and "elite confirmation" are priorities, not phases, so
worker utilization stays ~100% whenever eligible work exists:

* min_cheap floor: elite seeds are withheld until at least `min_cheap`
  candidates have completed screening (an early, unrepresentative leader
  would otherwise set a miscalibrated qualification bar).
* UCB gate: a screened candidate qualifies for the elite bucket iff
  mean + z*sigma*sqrt(1/n) >= leader_mean - epsilon, with sigma pooled
  across all candidates' repeat seeds (prior until enough data).
* Evictable bucket: a later, stronger qualifier evicts the weakest
  provisional elite; the evictee's in-flight seeds are killed (real
  process termination, not censoring). Gate-passing candidates that
  find the bucket full wait on a ranked waitlist and are pulled in
  whenever a demotion/eviction opens a slot.
* Racing demotions: an elite whose UCB falls confidently below the
  bucket leader is demoted mid-confirmation and its in-flight seeds
  killed.
* Allocation: among bucket members still needing seeds, dispatch to the
  highest-UCB unresolved member (LUCB-flavored: effort concentrates on
  candidates whose fate is undecided).
* Fill-and-stop: once the bucket is full (and the floor met), no new
  cheap screening is dispatched and in-flight cheap seeds are killed
  (configurable).

Workers are dedicated processes (forkserver context) running a simple
task loop; kill == terminate + respawn, so preemption frees the CPU
immediately even mid-task. The evaluation callable is injected (any
picklable module-level function), which is how the mock validation in
tests/test_queue_v2.py drives the whole workflow without the engine.

Scientific semantics preserved by interface: the generation result
reports BOTH the confirmed winner (decision) and the top-K-by-mean
recombination set (search dynamics), per the agreed decoupling.
"""
from __future__ import annotations

import math
import os
import time
import multiprocessing as mp
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Callable, Optional


# =============================================================================
# Worker pool with real preemption
# =============================================================================

def _qs_worker_loop(conn, eval_fn):
    """Dedicated worker: one duplex pipe carries tasks in and results out.
    A pipe is private to its worker, so a kill (terminate) can corrupt at
    most this worker's stream — which the pool discards and replaces at
    respawn. (A SHARED queue would be corrupted for everyone by a kill
    landing mid-write; that is why this is a pipe-per-worker design.)
    Exceptions become clean failure results — the process survives."""
    while True:
        item = conn.recv()
        if item is None:
            return
        attempt_id, payload = item
        try:
            out = eval_fn(payload)
            if not isinstance(out, dict):
                out = {"phi": float(out)}
            out.setdefault("success", True)
        except Exception as e:  # noqa: BLE001 — worker must never die on a job
            out = {"success": False, "error": repr(e)}
        conn.send((attempt_id, out))


class PreemptivePool:
    """N dedicated worker processes, one duplex pipe each; targeted kill =
    terminate + respawn with a FRESH pipe (kill-safe by construction)."""

    def __init__(self, n_workers: int, eval_fn: Callable[[dict], dict],
                 mp_context: Optional[str] = None):
        ctx_name = mp_context or os.environ.get("RCP_MP_CONTEXT", "forkserver")
        self._ctx = mp.get_context(ctx_name)
        self._eval_fn = eval_fn
        self.n_workers = n_workers
        self._conns: list[Any] = [None] * n_workers   # parent ends
        self._procs: list[Any] = [None] * n_workers
        self.respawn_count = 0
        for w in range(n_workers):
            self._spawn(w)

    def _spawn(self, w: int):
        parent_conn, child_conn = self._ctx.Pipe(duplex=True)
        p = self._ctx.Process(
            target=_qs_worker_loop,
            args=(child_conn, self._eval_fn),
            daemon=True)
        p.start()
        child_conn.close()      # parent keeps only its end
        self._conns[w] = parent_conn
        self._procs[w] = p

    def submit(self, w: int, attempt_id: int, payload: dict):
        try:
            self._conns[w].send((attempt_id, payload))
        except (BrokenPipeError, OSError):
            # worker died since the last poll (e.g., OOM kill): respawn
            # with a fresh pipe and resubmit — the task is not lost.
            self.ensure_alive(w)
            self._conns[w].send((attempt_id, payload))

    def kill(self, w: int):
        """Terminate worker w mid-task and respawn with a fresh pipe."""
        p = self._procs[w]
        if p is not None and p.is_alive():
            p.terminate()
            p.join(timeout=10.0)
        try:
            self._conns[w].close()
        except Exception:
            pass
        self._spawn(w)
        self.respawn_count += 1

    def ensure_alive(self, w: int) -> bool:
        """True if the worker is alive; respawns (and returns False) if it
        died unexpectedly (e.g., OOM kill in real use)."""
        p = self._procs[w]
        if p is not None and p.is_alive():
            return True
        try:
            self._conns[w].close()
        except Exception:
            pass
        self._spawn(w)
        self.respawn_count += 1
        return False

    def next_result(self, timeout: float):
        """Wait for any worker's result: (worker_id, attempt_id, out), or
        None on timeout. A dead/corrupted pipe surfaces as (w, None, None)
        so the scheduler can fail the attempt and respawn."""
        from multiprocessing.connection import wait as _mpwait
        live = [c for c in self._conns if c is not None]
        ready = _mpwait(live, timeout=timeout)
        for conn in ready:
            w = self._conns.index(conn)
            try:
                aid, out = conn.recv()
                return (w, aid, out)
            except (EOFError, OSError):
                return (w, None, None)
        return None

    def shutdown(self):
        for w in range(self.n_workers):
            try:
                self._conns[w].send(None)
            except Exception:
                pass
        for p in self._procs:
            if p is not None:
                p.join(timeout=2.0)
                if p.is_alive():
                    p.terminate()


# =============================================================================
# Scheduler
# =============================================================================

class CandState(Enum):
    PENDING = "pending"        # not yet screened
    SCREENING = "screening"    # cheap seeds in flight
    WAITLIST = "waitlist"      # gate-passed, bucket full
    ELITE = "elite"            # in the bucket, confirmation seeds flow
    CONFIRMED = "confirmed"    # elite with full seed quota
    DROPPED = "dropped"        # gate-failed after screening
    DEMOTED = "demoted"        # raced out of the bucket
    EVICTED = "evicted"        # pushed out by a stronger qualifier
    ABANDONED = "abandoned"    # never screened (bucket filled first)


@dataclass
class QSConfig:
    n_workers: int = 5
    cheap_seeds: int = 2          # seeds per candidate to count as screened
    min_cheap: int = 10           # candidates screened before elite seeds flow
    max_elite: int = 4            # bucket capacity
    min_elite: int = 2            # always confirm at least this many
    elite_seeds: int = 6          # total seeds for a confirmed elite
    z: float = 2.0                # one-sided confidence for the gate
    epsilon: float = 5e-4         # indifference zone in phi
    sigma_prior: float = 2e-3     # pooled-sigma fallback until data exists
    stop_cheap_when_filled: bool = True
    kill_cheap_when_filled: bool = True
    result_timeout_s: float = 30.0  # watchdog for lost workers


@dataclass
class _Cand:
    cid: int
    theta: Any
    state: CandState = CandState.PENDING
    results: list = field(default_factory=list)
    inflight: set = field(default_factory=set)   # attempt ids

    @property
    def n(self) -> int:
        return len(self.results)

    @property
    def mean(self) -> float:
        return sum(self.results) / len(self.results) if self.results else float("-inf")


class QueueScheduler:
    """Runs ONE generation: a fixed population of thetas through the
    no-phases workflow. Returns the decision + diagnostics."""

    def __init__(self, cfg: QSConfig, thetas: list, eval_fn: Callable,
                 pool: Optional[PreemptivePool] = None, log: bool = False):
        self.cfg = cfg
        self.cands = [_Cand(cid=i, theta=t) for i, t in enumerate(thetas)]
        self.pool = pool or PreemptivePool(cfg.n_workers, eval_fn)
        self._own_pool = pool is None
        self.log = log
        # attempt bookkeeping
        self._next_attempt = 0
        self._attempt_info: dict[int, tuple[int, str, int]] = {}  # aid -> (cid, stage, worker)
        self._killed_attempts: set[int] = set()
        self._worker_busy: list[Optional[int]] = [None] * cfg.n_workers
        # diagnostics
        self.events: list[str] = []
        self._busy_time = [0.0] * cfg.n_workers
        self._attempt_start: dict[int, float] = {}
        self._t0: Optional[float] = None   # set at first dispatch
        self._idle_while_eligible = 0.0
        self.kills = 0

    # ---- statistics ------------------------------------------------------

    def pooled_sigma(self) -> float:
        ss, df = 0.0, 0
        for c in self.cands:
            if c.n >= 2:
                m = c.mean
                ss += sum((r - m) ** 2 for r in c.results)
                df += c.n - 1
        if df < 3:
            return self.cfg.sigma_prior
        return math.sqrt(ss / df)

    def _ucb(self, c: _Cand) -> float:
        if c.n == 0:
            return float("inf")
        return c.mean + self.cfg.z * self.pooled_sigma() / math.sqrt(c.n)

    def _leader_mean(self) -> float:
        means = [c.mean for c in self.cands if c.n >= self.cfg.cheap_seeds]
        return max(means) if means else float("-inf")

    def _gate(self, c: _Cand) -> bool:
        return self._ucb(c) >= self._leader_mean() - self.cfg.epsilon

    # ---- bucket ----------------------------------------------------------

    def _bucket(self) -> list:
        return [c for c in self.cands
                if c.state in (CandState.ELITE, CandState.CONFIRMED)]

    def _screened_count(self) -> int:
        return sum(1 for c in self.cands
                   if c.n >= self.cfg.cheap_seeds or c.state in
                   (CandState.DROPPED, CandState.DEMOTED, CandState.EVICTED))

    def _bucket_full(self) -> bool:
        return len(self._bucket()) >= self.cfg.max_elite

    def _promote(self, c: _Cand):
        c.state = CandState.ELITE
        self._log(f"promote c{c.cid} mean={c.mean:.5f}")

    def _kill_inflight(self, c: _Cand):
        for aid in list(c.inflight):
            info = self._attempt_info.get(aid)
            if info is None:
                continue
            _, _, w = info
            if self._worker_busy[w] == aid:
                self.pool.kill(w)
                self._account_busy(aid, w)
                self._worker_busy[w] = None
                self.kills += 1
            self._killed_attempts.add(aid)
            c.inflight.discard(aid)

    def _evict_weakest_for(self, newcomer: _Cand) -> bool:
        bucket = sorted(self._bucket(), key=lambda c: c.mean)
        weakest = bucket[0]
        if newcomer.mean <= weakest.mean:
            return False
        weakest.state = CandState.EVICTED
        self._kill_inflight(weakest)
        self._log(f"evict c{weakest.cid} (mean={weakest.mean:.5f}) "
                  f"for c{newcomer.cid} (mean={newcomer.mean:.5f})")
        self._promote(newcomer)
        return True

    def _pull_waitlist(self):
        while not self._bucket_full():
            wl = [c for c in self.cands if c.state == CandState.WAITLIST]
            if not wl:
                return
            best = max(wl, key=lambda c: c.mean)
            self._promote(best)

    def _racing_check(self):
        """Demote elites whose UCB falls confidently below the bucket
        leader. min_elite strongest members are immune."""
        bucket = sorted(self._bucket(), key=lambda c: c.mean, reverse=True)
        if len(bucket) <= self.cfg.min_elite:
            return
        lead_mean = bucket[0].mean
        for c in bucket[self.cfg.min_elite:]:
            if c.n >= self.cfg.cheap_seeds and \
               self._ucb(c) < lead_mean - self.cfg.epsilon:
                c.state = CandState.DEMOTED
                self._kill_inflight(c)
                self._log(f"demote c{c.cid} mean={c.mean:.5f} "
                          f"ucb={self._ucb(c):.5f} < lead-eps")
                self._pull_waitlist()

    # ---- screening completion -------------------------------------------

    def _on_screened(self, c: _Cand):
        if not self._gate(c):
            c.state = CandState.DROPPED
            self._log(f"drop c{c.cid} mean={c.mean:.5f} (gate)")
        elif not self._bucket_full():
            self._promote(c)
        elif not self._evict_weakest_for(c):
            c.state = CandState.WAITLIST
            self._log(f"waitlist c{c.cid} mean={c.mean:.5f}")
        # a new arrival can change the leader → re-gate the bucket
        self._racing_check()
        # bucket may have just filled (or the floor may have just been
        # met with a full bucket) → stop/kill cheap. Checked on EVERY
        # screening completion, including gate-failures.
        if self._bucket_full() and \
           self._screened_count() >= self.cfg.min_cheap and \
           self.cfg.stop_cheap_when_filled:
            self._halt_cheap()

    def _halt_cheap(self):
        for c in self.cands:
            if c.state == CandState.PENDING:
                c.state = CandState.ABANDONED
            elif c.state == CandState.SCREENING and \
                    self.cfg.kill_cheap_when_filled:
                self._kill_inflight(c)
                c.state = CandState.ABANDONED
                self._log(f"abandon in-flight cheap c{c.cid}")

    # ---- dispatch --------------------------------------------------------

    def _next_work(self) -> Optional[tuple[int, str]]:
        """Pick the next (cid, stage) per the priority rules."""
        floor_met = self._screened_count() >= self.cfg.min_cheap
        cheap_allowed = not (self._bucket_full() and floor_met
                             and self.cfg.stop_cheap_when_filled)
        # 1. cheap seeds (always until the floor; after that, while allowed)
        if cheap_allowed:
            for c in self.cands:
                if c.state in (CandState.PENDING, CandState.SCREENING) and \
                        c.n + len(c.inflight) < self.cfg.cheap_seeds:
                    return (c.cid, "cheap")
        # 2. elite seeds (withheld until the floor is met)
        if floor_met:
            needing = [c for c in self._bucket()
                       if c.n + len(c.inflight) < self.cfg.elite_seeds]
            if needing:
                pick = max(needing, key=self._ucb)   # LUCB-flavored
                return (pick.cid, "elite")
        return None

    def _dispatch(self, w: int, cid: int, stage: str):
        c = self.cands[cid]
        aid = self._next_attempt
        self._next_attempt += 1
        seed = c.n + len(c.inflight)
        payload = {"theta": c.theta, "candidate": cid, "seed": seed,
                   "stage": stage}
        self._attempt_info[aid] = (cid, stage, w)
        now = time.monotonic()
        if self._t0 is None:
            self._t0 = now
        self._attempt_start[aid] = now
        c.inflight.add(aid)
        if c.state == CandState.PENDING:
            c.state = CandState.SCREENING
        self._worker_busy[w] = aid
        self.pool.submit(w, aid, payload)

    def _account_busy(self, aid: int, w: int):
        t0 = self._attempt_start.pop(aid, None)
        if t0 is not None:
            self._busy_time[w] += time.monotonic() - t0

    # ---- main loop -------------------------------------------------------

    def run(self) -> dict:
        cfg = self.cfg
        while True:
            # fill idle workers
            dispatched = False
            for w in range(cfg.n_workers):
                if self._worker_busy[w] is not None:
                    continue
                work = self._next_work()
                if work is None:
                    break
                self._dispatch(w, *work)
                dispatched = True
            inflight = sum(1 for b in self._worker_busy if b is not None)
            if inflight == 0:
                if self._next_work() is None:
                    break       # generation complete
                if not dispatched:
                    break       # safety: nothing dispatchable
                continue
            # wait for a result on any worker pipe
            got = self.pool.next_result(timeout=cfg.result_timeout_s)
            if got is None:
                # watchdog: check for dead workers (external kill / OOM)
                for w in range(cfg.n_workers):
                    if self._worker_busy[w] is not None and \
                            not self.pool.ensure_alive(w):
                        dead_aid = self._worker_busy[w]
                        self._log(f"worker {w} died; attempt {dead_aid} lost")
                        self._fail_attempt(dead_aid, w)
                continue
            w, aid, out = got
            if aid is None:
                # pipe died (worker killed externally mid-send)
                self.pool.ensure_alive(w)
                dead_aid = self._worker_busy[w]
                if dead_aid is not None:
                    self._log(f"worker {w} pipe died; attempt {dead_aid} lost")
                    self._fail_attempt(dead_aid, w)
                continue
            self._on_result(aid, out)
        util = self.utilization()
        bucket = sorted(self._bucket(), key=lambda c: c.mean, reverse=True)
        confirmed = [c for c in bucket if c.n >= cfg.elite_seeds]
        winner = confirmed[0] if confirmed else (bucket[0] if bucket else None)
        by_mean = sorted([c for c in self.cands if c.n > 0],
                         key=lambda c: c.mean, reverse=True)
        out = {
            "winner": None if winner is None else
                {"cid": winner.cid, "theta": winner.theta,
                 "mean_phi": winner.mean, "n": winner.n},
            "recombination_set": [
                {"cid": c.cid, "theta": c.theta, "mean_phi": c.mean, "n": c.n}
                for c in by_mean[:cfg.max_elite]],
            "elites": [{"cid": c.cid, "mean_phi": c.mean, "n": c.n}
                       for c in bucket],
            "pooled_sigma": self.pooled_sigma(),
            "total_seeds": sum(c.n for c in self.cands),
            "kills": self.kills,
            "respawns": self.pool.respawn_count,
            "utilization": util,
            "states": {c.cid: c.state.value for c in self.cands},
            "events": list(self.events),
        }
        if self._own_pool:
            self.pool.shutdown()
        return out

    def _fail_attempt(self, aid: int, w: int):
        cid, _, _ = self._attempt_info[aid]
        c = self.cands[cid]
        c.inflight.discard(aid)
        self._killed_attempts.add(aid)
        self._account_busy(aid, w)
        self._worker_busy[w] = None

    def _on_result(self, aid: int, out: dict):
        if aid in self._killed_attempts:
            return                       # late result from a killed attempt
        cid, stage, w = self._attempt_info[aid]
        c = self.cands[cid]
        c.inflight.discard(aid)
        self._account_busy(aid, w)
        if self._worker_busy[w] == aid:
            self._worker_busy[w] = None
        if not out.get("success", True):
            self._log(f"[heads-up] run failed and skipped (c{cid} "
                      f"stage={stage}): {out.get('error', '?')}")
            return
        c.results.append(float(out["phi"]))
        if stage == "cheap" and c.n >= self.cfg.cheap_seeds and \
                c.state == CandState.SCREENING:
            self._on_screened(c)
        elif stage == "elite":
            self._racing_check()
            if c.state == CandState.ELITE and c.n >= self.cfg.elite_seeds:
                c.state = CandState.CONFIRMED
                self._log(f"confirm c{c.cid} mean={c.mean:.5f} n={c.n}")

    # ---- diagnostics -----------------------------------------------------

    def utilization(self) -> float:
        if self._t0 is None:
            return 1.0
        span = time.monotonic() - self._t0
        if span <= 0:
            return 1.0
        return sum(self._busy_time) / (self.cfg.n_workers * span)

    def _log(self, msg: str):
        self.events.append(msg)
        if self.log:
            print(f"[qs] {msg}", flush=True)


# =============================================================================
# Multi-generation driver with the stopping rule
# =============================================================================

def run_mock_search(propose_fn: Callable[[int], list], eval_fn: Callable,
                    cfg: QSConfig, max_generations: int = 12,
                    stop_eps: float = 1e-3, stop_G: int = 2,
                    log: bool = False) -> dict:
    """Generation loop: propose -> QueueScheduler -> track best; stop when
    the best improves by < stop_eps for stop_G consecutive generations."""
    best = float("-inf")
    flat = 0
    history = []
    pool = PreemptivePool(cfg.n_workers, eval_fn)
    gen_run = 0
    for gen in range(max_generations):
        thetas = propose_fn(gen)
        sched = QueueScheduler(cfg, thetas, eval_fn, pool=pool, log=log)
        res = sched.run()
        gen_run += 1
        w = res["winner"]
        gen_best = w["mean_phi"] if w else float("-inf")
        improved = gen_best - best
        history.append({"gen": gen, "best": gen_best, "improvement": improved,
                        "seeds": res["total_seeds"], "util": res["utilization"]})
        if improved < stop_eps:
            flat += 1
        else:
            flat = 0
        best = max(best, gen_best)
        if flat >= stop_G:
            break
    pool.shutdown()
    return {"best": best, "generations_run": gen_run, "history": history}
