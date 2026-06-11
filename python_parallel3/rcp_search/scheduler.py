"""Adaptive multi-fidelity scheduler.

Drives candidate dispatch under a fixed budget with three observable goals:

* **Cores fed** — no barrier between cheap screen and confirmation. Every
  completion event triggers re-dispatch; idle only when nothing decision-
  useful remains.
* **Cheap drops → new perturbations.** A candidate whose ``mean(samples)`` is
  worse than ``current_best - cheap_delta`` is dropped; its unused cheap
  seeds spawn a fresh ``propose_next_theta()`` call.
* **Confirm drops → survivors.** A confirming candidate whose
  ``mean + SE < max(prior_best, leader.mean - leader.SE)`` is dropped; its
  unused confirmation seeds are redistributed to surviving confirming
  elites.

The full audit set (per-attempt, dispatch, drop, budget, idle, decision,
worker, summary) is written by :class:`LogSink` under the configured
``out_dir``.

Optional. The search loop uses this module when
``CONFIG['scheduler_mode'] == 'adaptive'``.
"""
from __future__ import annotations


import csv
import enum
import hashlib
import heapq
import json
import math
import pathlib
import statistics
import time
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Callable, Optional


# =============================================================================
# Configuration
# =============================================================================


@dataclass
class SchedulerConfig:
    """All adaptive-scheduler tunable parameters with defaults sized for a
    small local smoke test. Real Colab/Cloud runs override these from the
    user's CONFIG dict."""

    # --- Population and budget ---
    population_size: int = 12
    cheap_seeds_per_candidate: int = 2
    elite_count: int = 2
    elite_extra_seeds: int = 6

    # --- Drop rules ---
    cheap_delta: float = 0.005        # mean_so_far < current_best - cheap_delta -> drop
    # Confirmation floor formula has no extra knob; uses SE directly.

    # --- Pathological detection ---
    phi_noise_floor: float = 0.001
    unstable_range_tol: float = 10.0
    unstable_std_tol: float = 5.0

    # --- Worker pool ---
    n_workers: int = 4
    scheduler_seed: int = 0

    # --- Local throttling (Colab leaves this False) ---
    throttle_local: bool = False
    worker_cooldown_seconds: float = 1.5

    # --- Safety stop ---
    real_smoke_max_wall_seconds: float = 600.0

    # --- Cap on samples per candidate (post-redistribution) ---
    max_extra_seeds_per_candidate: int = 32

    # --- Audit ---
    log_to_disk: bool = True


# =============================================================================
# Enumerations
# =============================================================================


class CandidateState(str, enum.Enum):
    PENDING = "pending"               # not yet dispatched
    CHEAP_RUNNING = "cheap_running"   # has at least one cheap attempt
    CHEAP_DECIDED = "cheap_decided"   # cheap quota met, awaiting promotion or post-cheap
    CONFIRMING = "confirming"         # confirmation samples in flight or completed
    CONFIRMED = "confirmed"           # reached elite_extra_seeds; eligible for "best"
    DROPPED = "dropped"               # dropped by cheap-delta or confirm-floor rule
    PATHOLOGICAL = "pathological"     # unstable / bimodal flagged


_FROZEN = {CandidateState.CONFIRMED, CandidateState.DROPPED,
           CandidateState.PATHOLOGICAL}


class Stage(str, enum.Enum):
    CHEAP = "cheap"
    CONFIRM = "confirm"


# =============================================================================
# Data classes
# =============================================================================


@dataclass
class Candidate:
    candidate_id: int
    theta: Any = None
    state: CandidateState = CandidateState.PENDING
    cheap_phi: list = field(default_factory=list)
    confirm_phi: list = field(default_factory=list)
    n_active: int = 0
    n_failed: int = 0
    n_censored: int = 0
    runtime_history: list = field(default_factory=list)
    state_history: list = field(default_factory=list)
    last_reason: str = ""
    used_seeds: set = field(default_factory=set)
    # Confirmation bookkeeping
    confirm_target: int = 0   # original elite_extra_seeds; eligibility check
    confirm_quota: int = 0    # current dispatch target; can grow with redistribution
    cheap_dispatched: int = 0
    confirm_dispatched: int = 0

    @property
    def all_phi(self):
        return self.cheap_phi + self.confirm_phi

    @property
    def n_cheap_valid(self): return len(self.cheap_phi)

    @property
    def n_confirm_valid(self): return len(self.confirm_phi)

    @property
    def n_valid_total(self): return self.n_cheap_valid + self.n_confirm_valid

    @property
    def mean_phi(self):
        ap = self.all_phi
        return statistics.fmean(ap) if ap else None

    @property
    def mean_cheap(self):
        return statistics.fmean(self.cheap_phi) if self.cheap_phi else None

    @property
    def std_phi(self):
        ap = self.all_phi
        return statistics.stdev(ap) if len(ap) >= 2 else None

    @property
    def se_phi(self):
        s = self.std_phi
        if s is None or self.n_valid_total < 2:
            return None
        return s / math.sqrt(self.n_valid_total)

    @property
    def best_phi(self):
        ap = self.all_phi
        return max(ap) if ap else None

    @property
    def phi_range(self):
        ap = self.all_phi
        return 0.0 if len(ap) < 2 else max(ap) - min(ap)

    def record_transition(self, new_state, reason, t):
        if new_state == self.state:
            return False
        self.state_history.append((t, self.state.value, new_state.value, reason))
        self.state = new_state
        self.last_reason = reason
        return True


@dataclass
class Attempt:
    attempt_id: int
    candidate_id: int
    theta: Any
    seed_id: int
    stage: Stage
    scheduled_at: float = 0.0
    started_at: Optional[float] = None
    ended_at: Optional[float] = None
    runtime: Optional[float] = None
    phi: Optional[float] = None
    failed: bool = False
    censored: bool = False
    completed: bool = False
    worker_id: int = -1
    decision_reason: str = ""


@dataclass
class Event:
    fake_time: float
    event_type: str  # "completed", "failed", "censored"
    attempt_id: int
    phi: Optional[float] = None
    payload: dict = field(default_factory=dict)


# =============================================================================
# Runner interface
# =============================================================================


class Runner:
    """Abstract worker pool."""

    n_workers: int

    def now(self) -> float: raise NotImplementedError
    def n_idle(self) -> int: raise NotImplementedError
    def dispatch(self, attempt: Attempt) -> None: raise NotImplementedError
    def next_event(self) -> Optional[Event]: raise NotImplementedError
    def in_flight_count(self) -> int: raise NotImplementedError
    def worker_timeline(self) -> list: return []
    def total_busy_time(self) -> float: return 0.0
    def shutdown(self): pass


# =============================================================================
# Log sink
# =============================================================================


class LogSink:
    """CSV + JSON log writers. Captures the three goal-signals explicitly:
    - dispatch_log: every dispatch with its reason
    - drop_log: every cheap/confirm drop with seeds_returned + redirected_to
    - worker_timeline: per-worker busy/idle segments with idle_reason
    """

    def __init__(self, out_dir: pathlib.Path):
        self.out_dir = pathlib.Path(out_dir)
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.events: list = []
        self.dispatch_log: list = []
        self.drop_log: list = []
        self.budget_log: list = []   # snapshots of budget state on key events
        self.idle_reasons: list = [] # (worker_id, start, end, reason)

    def log_event(self, t, etype, a: Attempt, candidate_state="", reason=""):
        self.events.append({
            "fake_time": t,
            "event_type": etype,
            "attempt_id": a.attempt_id,
            "candidate_id": a.candidate_id,
            "seed_id": a.seed_id,
            "stage": a.stage.value if isinstance(a.stage, Stage) else a.stage,
            "runtime": a.runtime if a.runtime is not None else "",
            "phi": a.phi if a.phi is not None else "",
            "failed": int(a.failed),
            "censored": int(a.censored),
            "worker_id": a.worker_id,
            "decision_reason": reason or a.decision_reason,
            "candidate_state_after": candidate_state,
        })

    def log_dispatch(self, t, a: Attempt, reason):
        self.dispatch_log.append({
            "fake_time": t,
            "attempt_id": a.attempt_id,
            "candidate_id": a.candidate_id,
            "stage": a.stage.value if isinstance(a.stage, Stage) else a.stage,
            "seed_id": a.seed_id,
            "decision_reason": reason,
        })

    def log_drop(self, t, candidate, kind, seeds_returned, redirected_to=""):
        self.drop_log.append({
            "fake_time": t,
            "candidate_id": candidate.candidate_id,
            "kind": kind,           # "cheap_drop" or "confirm_drop" or "pathological"
            "n_valid_at_drop": candidate.n_valid_total,
            "n_cheap": candidate.n_cheap_valid,
            "n_confirm": candidate.n_confirm_valid,
            "mean_phi": candidate.mean_phi if candidate.mean_phi is not None else "",
            "seeds_returned": seeds_returned,
            "redirected_to": redirected_to,
            "reason": candidate.last_reason,
        })

    def log_budget(self, t, label, cheap_remaining, confirm_remaining,
                    cheap_perturbations_total, confirm_dispatched_total):
        self.budget_log.append({
            "fake_time": t,
            "label": label,
            "cheap_budget_remaining": cheap_remaining,
            "confirm_budget_remaining": confirm_remaining,
            "cheap_perturbations_total": cheap_perturbations_total,
            "confirm_dispatched_total": confirm_dispatched_total,
        })

    def log_idle_segment(self, worker_id, start, end, reason):
        self.idle_reasons.append({
            "worker_id": worker_id,
            "segment_start": start,
            "segment_end": end,
            "idle_reason": reason,
        })

    def _write_rows(self, fname, rows):
        path = self.out_dir / fname
        if not rows:
            path.write_text("")
            return path
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        return path

    def write_attempt_log(self):
        return self._write_rows("per_job_event_log.csv", self.events)

    def write_dispatch_log(self):
        return self._write_rows("dispatch_log.csv", self.dispatch_log)

    def write_drop_log(self):
        return self._write_rows("drop_log.csv", self.drop_log)

    def write_budget_log(self):
        return self._write_rows("budget_log.csv", self.budget_log)

    def write_idle_reasons(self):
        return self._write_rows("idle_reasons.csv", self.idle_reasons)

    def write_candidate_log(self, candidates, cfg):
        path = self.out_dir / "per_candidate_decision_log.csv"
        rows = []
        for c in candidates.values():
            rows.append({
                "candidate_id": c.candidate_id,
                "n_cheap": c.n_cheap_valid,
                "n_confirm": c.n_confirm_valid,
                "n_total": c.n_valid_total,
                "n_failed": c.n_failed,
                "n_censored": c.n_censored,
                "mean_phi": c.mean_phi if c.mean_phi is not None else "",
                "std_phi": c.std_phi if c.std_phi is not None else "",
                "se_phi": c.se_phi if c.se_phi is not None else "",
                "best_phi": c.best_phi if c.best_phi is not None else "",
                "phi_range": c.phi_range,
                "confirm_target": c.confirm_target,
                "confirm_quota": c.confirm_quota,
                "final_state": c.state.value,
                "final_reason": c.last_reason,
                "state_history_json": json.dumps(c.state_history),
            })
        return self._write_rows("per_candidate_decision_log.csv", rows)

    def write_worker_timeline(self, runner: Runner):
        path = self.out_dir / "worker_timeline.csv"
        rows = runner.worker_timeline()
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["worker_id", "segment_start", "segment_end",
                        "attempt_id", "stage"])
            for r in rows:
                w.writerow(r)
        return path

    def write_decision_hash(self, candidates):
        items = []
        for cid in sorted(candidates):
            c = candidates[cid]
            items.append((cid, c.state.value,
                          round(c.mean_phi or -1.0, 6),
                          c.n_cheap_valid, c.n_confirm_valid,
                          c.n_failed, c.n_censored,
                          tuple(round(p, 6) for p in c.all_phi)))
        h = hashlib.sha256(repr(items).encode()).hexdigest()
        (self.out_dir / "decision_log_hash.txt").write_text(h + "\n")
        return h

    def write_summary(self, summary: dict):
        path = self.out_dir / "generation_summary.json"
        with open(path, "w") as f:
            json.dump(summary, f, indent=2, default=str)
        return path


# =============================================================================
# Adaptive scheduler
# =============================================================================


class AdaptiveScheduler:
    """Compute-allocation scheduler.

    Inputs:
      - cfg: SchedulerConfig
      - runner: Runner instance (Fake or Real)
      - log_sink: LogSink
      - propose_next_theta: callable () -> theta, called when cheap budget
                            has saved compute and a fresh perturbation can fit
      - initial_thetas: list of theta arrays for the bootstrap population
      - prior_best: best mean phi from previous generations (or -math.inf)
    """

    name = "adaptive_v2"

    def __init__(self, cfg: SchedulerConfig, runner: Runner,
                 log_sink: LogSink,
                 propose_next_theta: Callable[[], Any]):
        self.cfg = cfg
        self.runner = runner
        self.log = log_sink
        self.propose_next_theta = propose_next_theta
        self.candidates: dict = {}
        self.attempts: dict = {}
        self._next_aid = 0
        self._next_cid = 0
        # Compute budgets
        self.cheap_budget_remaining = (cfg.population_size *
                                        cfg.cheap_seeds_per_candidate)
        self.confirm_budget_remaining = cfg.elite_count * cfg.elite_extra_seeds
        # Reference
        self.prior_best = -math.inf
        self.current_best = -math.inf   # max(prior_best, best CONFIRMED this gen)
        # Phase
        self.phase = Stage.CHEAP

    # ---------- helpers ----------

    def _new_candidate(self, theta) -> Candidate:
        c = Candidate(candidate_id=self._next_cid, theta=theta)
        self._next_cid += 1
        self.candidates[c.candidate_id] = c
        return c

    def _new_attempt(self, c: Candidate, stage: Stage) -> Attempt:
        seed = 0
        while seed in c.used_seeds:
            seed += 1
        c.used_seeds.add(seed)
        a = Attempt(
            attempt_id=self._next_aid,
            candidate_id=c.candidate_id,
            theta=c.theta,
            seed_id=seed,
            stage=stage,
            scheduled_at=self.runner.now(),
        )
        self._next_aid += 1
        self.attempts[a.attempt_id] = a
        return a

    def _floor(self):
        """Confirmation floor = max(prior_best, leader.mean - leader.SE).

        Leader is the current confirming candidate with the highest mean among
        those that are not yet DROPPED/PATHOLOGICAL/CONFIRMED-without-leadership.
        We include CONFIRMED candidates as potential leaders too — their stats
        are real and contribute to the floor.
        """
        cfg = self.cfg
        eligible = [c for c in self.candidates.values()
                    if c.state in (CandidateState.CONFIRMING,
                                    CandidateState.CONFIRMED)
                    and c.mean_phi is not None]
        if not eligible:
            return self.prior_best
        leader = max(eligible, key=lambda c: c.mean_phi)
        se = leader.se_phi
        if se is None:
            se = cfg.phi_noise_floor
        return max(self.prior_best, leader.mean_phi - se)

    def _is_pathological(self, c: Candidate) -> bool:
        cfg = self.cfg
        if c.n_valid_total < 2:
            return False
        if c.phi_range > cfg.unstable_range_tol * cfg.phi_noise_floor:
            return True
        s = c.std_phi
        if s is not None and s > cfg.unstable_std_tol * cfg.phi_noise_floor:
            return True
        return False

    # ---------- bootstrap ----------

    def _bootstrap_population(self, initial_thetas):
        for theta in initial_thetas:
            self._new_candidate(theta)

    # ---------- cheap-phase logic ----------

    def _evaluate_cheap_after_completion(self, c: Candidate):
        """After each cheap completion, decide drop / continue / promote-ready."""
        cfg = self.cfg
        # Pathological first
        if self._is_pathological(c):
            saved = self._seeds_returned_on_cheap_drop(c)
            c.record_transition(CandidateState.PATHOLOGICAL,
                                 f"range={c.phi_range:.4f}>{cfg.unstable_range_tol}*"
                                 f"{cfg.phi_noise_floor}",
                                 self.runner.now())
            self.log.log_drop(self.runner.now(), c, "pathological_cheap",
                              saved, redirected_to="new_perturbation_queue")
            self.cheap_budget_remaining += saved
            return
        # Drop?
        if (self.current_best > -math.inf
                and c.mean_phi is not None
                and c.mean_phi < self.current_best - cfg.cheap_delta):
            saved = self._seeds_returned_on_cheap_drop(c)
            c.record_transition(
                CandidateState.DROPPED,
                f"cheap_drop: mean={c.mean_phi:.5f} < "
                f"current_best={self.current_best:.5f} - delta={cfg.cheap_delta}",
                self.runner.now(),
            )
            self.log.log_drop(self.runner.now(), c, "cheap_drop", saved,
                              redirected_to="new_perturbation_queue")
            self.cheap_budget_remaining += saved
            return
        # Cheap quota met?
        if c.n_cheap_valid >= cfg.cheap_seeds_per_candidate:
            c.record_transition(CandidateState.CHEAP_DECIDED,
                                "cheap quota met", self.runner.now())

    def _seeds_returned_on_cheap_drop(self, c: Candidate) -> int:
        """How many seeds were planned but not used."""
        planned = self.cfg.cheap_seeds_per_candidate
        used_or_active = c.cheap_dispatched
        return max(0, planned - used_or_active)

    # ---------- confirm-phase logic ----------

    def _evaluate_confirm_after_completion(self, c: Candidate):
        cfg = self.cfg
        # Pathological?
        if self._is_pathological(c):
            saved = self._seeds_returned_on_confirm_drop(c)
            c.record_transition(CandidateState.PATHOLOGICAL,
                                 f"range={c.phi_range:.4f} unstable in confirm",
                                 self.runner.now())
            self.log.log_drop(self.runner.now(), c, "pathological_confirm",
                              saved, redirected_to="confirm_survivors")
            self._redistribute_saved_confirm(saved)
            return
        # Eligibility check: did we reach elite_extra_seeds?
        if c.n_confirm_valid >= cfg.elite_extra_seeds and c.state != CandidateState.CONFIRMED:
            # Move to CONFIRMED only when its target quota is met
            if c.n_confirm_valid >= c.confirm_target:
                c.record_transition(
                    CandidateState.CONFIRMED,
                    f"n_extra={c.n_confirm_valid}>=target={c.confirm_target}",
                    self.runner.now(),
                )
                # Update current_best if this candidate's mean beats it
                if c.mean_phi is not None and c.mean_phi > self.current_best:
                    self.current_best = c.mean_phi
        # Floor rule (skip if just confirmed)
        if c.state in (CandidateState.CONFIRMING, CandidateState.CONFIRMED):
            floor = self._floor()
            se = c.se_phi or cfg.phi_noise_floor
            if c.mean_phi is not None and c.mean_phi + se < floor:
                # Don't drop a candidate that has already CONFIRMED — it earned its spot
                if c.state == CandidateState.CONFIRMED:
                    return
                saved = self._seeds_returned_on_confirm_drop(c)
                c.record_transition(
                    CandidateState.DROPPED,
                    f"confirm_drop: mean={c.mean_phi:.5f}+se={se:.5f} "
                    f"< floor={floor:.5f}",
                    self.runner.now(),
                )
                self.log.log_drop(self.runner.now(), c, "confirm_drop",
                                  saved, redirected_to="confirm_survivors")
                self._redistribute_saved_confirm(saved)

    def _seeds_returned_on_confirm_drop(self, c: Candidate) -> int:
        return max(0, c.confirm_target - c.confirm_dispatched)

    def _redistribute_saved_confirm(self, saved: int):
        if saved <= 0:
            return
        survivors = [c for c in self.candidates.values()
                     if c.state == CandidateState.CONFIRMING]
        if not survivors:
            # No survivors left; just return budget to pool (will go unspent)
            self.confirm_budget_remaining += saved
            return
        share = saved // len(survivors)
        leftover = saved - share * len(survivors)
        cap = self.cfg.max_extra_seeds_per_candidate
        for i, s in enumerate(survivors):
            extra = share + (1 if i < leftover else 0)
            new_quota = min(s.confirm_target + extra, cap)
            # Increase confirm_quota (not confirm_target — eligibility unchanged)
            s.confirm_quota = max(s.confirm_quota, new_quota)
        self.confirm_budget_remaining += saved

    # ---------- dispatch ----------

    def _pick_cheap_attempt(self) -> Optional[Attempt]:
        cfg = self.cfg
        # 1) Round-robin priority: dispatch first sample for every candidate
        #    before any second sample. This guarantees the drop rule can fire
        #    between sample #1 and sample #2, freeing budget for perturbations.
        pool = [c for c in self.candidates.values()
                if c.state not in _FROZEN
                and c.state not in (CandidateState.CHEAP_DECIDED,
                                     CandidateState.CONFIRMING,
                                     CandidateState.CONFIRMED)
                and c.cheap_dispatched < cfg.cheap_seeds_per_candidate]
        # Sort by (fewest dispatched, lowest candidate_id) — first samples
        # for every candidate first.
        pool.sort(key=lambda c: (c.cheap_dispatched, c.candidate_id))
        if pool:
            return self._dispatch_attempt(pool[0], Stage.CHEAP)
        # 2) Budget remaining AND no pending pool work -> ask for new theta
        if self.cheap_budget_remaining > 0:
            theta = None
            if self.propose_next_theta is not None:
                try:
                    theta = self.propose_next_theta()
                except Exception:
                    theta = None
            if theta is not None:
                c = self._new_candidate(theta)
                return self._dispatch_attempt(c, Stage.CHEAP)
        return None

    def _pick_confirm_attempt(self) -> Optional[Attempt]:
        """Choose the confirming candidate that benefits most from one more
        sample. Priority: those still under their (possibly enlarged) quota,
        ordered by `dispatched - n_valid` (least-dispatched-first) so all
        survivors stay roughly balanced; ties broken by higher current mean
        (favor the leader)."""
        cfg = self.cfg
        if self.confirm_budget_remaining <= 0:
            return None
        candidates = [c for c in self.candidates.values()
                      if c.state == CandidateState.CONFIRMING
                      and c.confirm_dispatched < c.confirm_quota]
        if not candidates:
            return None
        candidates.sort(key=lambda c: (
            c.confirm_dispatched - c.n_confirm_valid,   # least pending first
            -(c.mean_phi if c.mean_phi is not None else -math.inf),  # then leader bias
            c.candidate_id,                             # FIFO
        ))
        return self._dispatch_attempt(candidates[0], Stage.CONFIRM)

    def _dispatch_attempt(self, c: Candidate, stage: Stage) -> Attempt:
        a = self._new_attempt(c, stage)
        c.n_active += 1
        if stage == Stage.CHEAP:
            c.cheap_dispatched += 1
            self.cheap_budget_remaining -= 1
            if c.state == CandidateState.PENDING:
                c.record_transition(CandidateState.CHEAP_RUNNING,
                                    "first cheap dispatch", self.runner.now())
        else:
            c.confirm_dispatched += 1
            self.confirm_budget_remaining -= 1
        a.started_at = self.runner.now()
        reason = (f"{stage.value}_dispatch:cid={c.candidate_id}:"
                  f"state={c.state.value}")
        a.decision_reason = reason
        self.log.log_dispatch(self.runner.now(), a, reason)
        self.runner.dispatch(a)
        return a

    def _dispatch_eligible(self):
        """Fill all idle workers with work."""
        while self.runner.n_idle() > 0:
            if self.phase == Stage.CHEAP:
                a = self._pick_cheap_attempt()
                if a is None:
                    # No cheap work; if confirm phase is open, try that
                    if self._cheap_phase_complete():
                        self._enter_confirm_phase()
                        continue
                    break
            else:
                a = self._pick_confirm_attempt()
                if a is None:
                    break

    def _cheap_phase_complete(self):
        """Cheap phase is done when (a) cheap_budget_remaining <= 0 AND no
        candidate has cheap_dispatched < target, AND (b) no cheap attempts
        are in flight, AND (c) no propose_next_theta queue work remains.
        """
        cfg = self.cfg
        # Any candidate still needs cheap samples?
        any_pending = any(c.cheap_dispatched < cfg.cheap_seeds_per_candidate
                          and c.state not in _FROZEN
                          and c.state != CandidateState.CHEAP_DECIDED
                          and c.state not in (CandidateState.CONFIRMING,
                                               CandidateState.CONFIRMED)
                          for c in self.candidates.values())
        if any_pending:
            return False
        # Budget remaining and propose_next_theta could give new work?
        if self.cheap_budget_remaining > 0 and self.propose_next_theta is not None:
            # We need to check by trying. Safer: leave to the dispatcher to
            # discover.
            pass
        # In-flight cheap attempts?
        in_flight_cheap = any(
            (not a.completed) and (not a.failed) and (not a.censored)
            and a.stage == Stage.CHEAP
            for a in self.attempts.values()
        )
        if in_flight_cheap:
            return False
        return True

    def _enter_confirm_phase(self):
        """Select top elite_count by mean among non-frozen non-pathological
        candidates with cheap_decided state. Move them to CONFIRMING."""
        if self.phase == Stage.CONFIRM:
            return
        cfg = self.cfg
        survivors = [c for c in self.candidates.values()
                     if c.state == CandidateState.CHEAP_DECIDED
                     and c.mean_phi is not None]
        survivors.sort(key=lambda c: -c.mean_phi)
        for c in survivors[: cfg.elite_count]:
            c.confirm_target = cfg.elite_extra_seeds
            c.confirm_quota = cfg.elite_extra_seeds
            c.record_transition(CandidateState.CONFIRMING,
                                 "promoted to confirmation (top "
                                 f"{cfg.elite_count} by mean)",
                                 self.runner.now())
        # Remaining cheap-decided are simply 'cheap_decided' until generation end
        self.phase = Stage.CONFIRM
        self.log.log_budget(self.runner.now(), "enter_confirm_phase",
                            self.cheap_budget_remaining,
                            self.confirm_budget_remaining,
                            len(self.candidates),
                            sum(c.confirm_dispatched for c in self.candidates.values()))

    # ---------- event handling ----------

    def _handle_completed(self, ev: Event):
        a = self.attempts[ev.attempt_id]
        a.ended_at = ev.fake_time
        a.runtime = a.ended_at - (a.started_at or a.scheduled_at)
        a.phi = ev.phi
        a.completed = True
        c = self.candidates[a.candidate_id]
        c.n_active = max(0, c.n_active - 1)
        c.runtime_history.append(a.runtime)
        if a.stage == Stage.CHEAP:
            c.cheap_phi.append(a.phi)
            self._evaluate_cheap_after_completion(c)
        else:
            c.confirm_phi.append(a.phi)
            self._evaluate_confirm_after_completion(c)
        self.log.log_event(ev.fake_time, "completed", a,
                           candidate_state=c.state.value,
                           reason=c.last_reason)

    def _handle_failed(self, ev: Event):
        a = self.attempts[ev.attempt_id]
        a.ended_at = ev.fake_time
        a.runtime = a.ended_at - (a.started_at or a.scheduled_at)
        a.failed = True
        c = self.candidates[a.candidate_id]
        c.n_active = max(0, c.n_active - 1)
        c.n_failed += 1
        self.log.log_event(ev.fake_time, "failed", a,
                           candidate_state=c.state.value,
                           reason="worker failure")

    def _handle_censored(self, ev: Event):
        a = self.attempts[ev.attempt_id]
        if a.completed:
            return
        a.censored = True
        a.ended_at = ev.fake_time
        a.runtime = a.ended_at - (a.started_at or a.scheduled_at)
        c = self.candidates[a.candidate_id]
        c.n_active = max(0, c.n_active - 1)
        c.n_censored += 1
        self.log.log_event(ev.fake_time, "censored", a,
                           candidate_state=c.state.value,
                           reason="censored")

    # ---------- main loop ----------

    def run_generation(self, initial_thetas, prior_best=-math.inf) -> dict:
        self.prior_best = float(prior_best) if prior_best is not None else -math.inf
        self.current_best = self.prior_best
        self._bootstrap_population(initial_thetas)
        self.log.log_budget(0.0, "start",
                            self.cheap_budget_remaining,
                            self.confirm_budget_remaining,
                            len(self.candidates), 0)
        self._dispatch_eligible()
        start_wall = self.runner.now() if hasattr(self.runner, "_wall_start_real") else 0.0
        max_iters = 100000
        for _ in range(max_iters):
            if self.runner.in_flight_count() == 0:
                # No work in flight; check if we can dispatch more.
                self._dispatch_eligible()
                if self.runner.in_flight_count() == 0:
                    # If workers are all in throttle cooldown, wait for one
                    # to come back online and try again.
                    if (hasattr(self.runner, "wait_for_available_worker")
                            and self.runner.n_idle() == 0):
                        got = self.runner.wait_for_available_worker(max_wait=10.0)
                        if got:
                            self._dispatch_eligible()
                            if self.runner.in_flight_count() > 0:
                                pass  # work dispatched, continue
                            else:
                                # truly nothing eligible — proceed with phase logic
                                pass
                    # Try moving to confirm phase if cheap is done
                    if self.runner.in_flight_count() == 0:
                        if self.phase == Stage.CHEAP and self._cheap_phase_complete():
                            self._enter_confirm_phase()
                            self._dispatch_eligible()
                            if self.runner.in_flight_count() == 0:
                                break  # confirm phase has no work either
                        else:
                            break
            # Safety stop
            if (self.cfg.real_smoke_max_wall_seconds > 0
                    and self.runner.now() > self.cfg.real_smoke_max_wall_seconds):
                self.log.log_budget(self.runner.now(), "safety_stop_max_wall",
                                    self.cheap_budget_remaining,
                                    self.confirm_budget_remaining,
                                    len(self.candidates),
                                    sum(c.confirm_dispatched
                                        for c in self.candidates.values()))
                break
            ev = self.runner.next_event()
            if ev is None:
                break
            et = ev.event_type
            if et == "completed":
                self._handle_completed(ev)
            elif et == "failed":
                self._handle_failed(ev)
            elif et == "censored":
                self._handle_censored(ev)
            self._dispatch_eligible()

        # Final transition: any leftover CONFIRMING that hit its target?
        for c in self.candidates.values():
            if (c.state == CandidateState.CONFIRMING
                    and c.n_confirm_valid >= c.confirm_target):
                c.record_transition(CandidateState.CONFIRMED,
                                     "final: target met", self.runner.now())
                if c.mean_phi is not None and c.mean_phi > self.current_best:
                    self.current_best = c.mean_phi
        return self._build_summary()

    def _build_summary(self):
        cfg = self.cfg
        now = self.runner.now()
        total_busy = self.runner.total_busy_time()
        # Candidates eligible as "best": CONFIRMED with n_confirm_valid >= elite_extra_seeds
        eligibles = [c for c in self.candidates.values()
                     if c.state == CandidateState.CONFIRMED]
        eligibles.sort(key=lambda c: -c.mean_phi)
        # New best for this generation (only if eligible exists)
        new_best_phi = None
        new_best_cid = None
        if eligibles:
            top = eligibles[0]
            new_best_phi = top.mean_phi
            new_best_cid = top.candidate_id
        return {
            "policy": self.name,
            "n_candidates_total": len(self.candidates),
            "n_workers": cfg.n_workers,
            "wall_time": now,
            "core_utilization": (total_busy / (cfg.n_workers * now)) if now > 0 else 0.0,
            "idle_core_seconds": max(0.0, cfg.n_workers * now - total_busy),
            "prior_best": self.prior_best,
            "current_best": self.current_best,
            "new_best_phi": new_best_phi,
            "new_best_candidate_id": new_best_cid,
            "n_eligibles": len(eligibles),
            "eligible_candidates": [c.candidate_id for c in eligibles],
            "cheap_budget_remaining": self.cheap_budget_remaining,
            "confirm_budget_remaining": self.confirm_budget_remaining,
            "n_cheap_dispatched_total": sum(c.cheap_dispatched
                                            for c in self.candidates.values()),
            "n_confirm_dispatched_total": sum(c.confirm_dispatched
                                              for c in self.candidates.values()),
            "n_dropped_cheap": sum(1 for c in self.candidates.values()
                                    if c.state == CandidateState.DROPPED
                                    and c.n_confirm_valid == 0),
            "n_dropped_confirm": sum(1 for c in self.candidates.values()
                                      if c.state == CandidateState.DROPPED
                                      and c.n_confirm_valid > 0),
            "n_pathological": sum(1 for c in self.candidates.values()
                                   if c.state == CandidateState.PATHOLOGICAL),
            "n_perturbations_proposed": (len(self.candidates) -
                                          cfg.population_size),
            "total_packings": (sum(c.cheap_dispatched
                                    for c in self.candidates.values())
                               + sum(c.confirm_dispatched
                                      for c in self.candidates.values())),
            "scheduler_seed": cfg.scheduler_seed,
        }


# =============================================================================
# Real ProcessPoolExecutor adapter (with optional per-worker cooldown)
# =============================================================================


class RealRunner(Runner):
    """Wraps a ProcessPoolExecutor.

    If cfg.throttle_local is True, each worker waits `worker_cooldown_seconds`
    after every completion before being marked available again. The scheduler's
    n_idle() reads from that "available" count, so dispatch naturally throttles.

    On Colab/Cloud, throttle_local stays False and cooldown is skipped.
    """

    def wait_for_available_worker(self, max_wait: float = 5.0) -> bool:
        """Sleep until at least one worker is available (cooldown expired).
        Returns True if a worker becomes available, False on timeout.
        Used when in_flight==0 and all workers are simultaneously cooling
        down, so the scheduler doesn't think the phase is complete."""
        import time as _t
        deadline = _t.monotonic() + max_wait
        while True:
            if self.n_idle() > 0:
                return True
            if _t.monotonic() >= deadline:
                return False
            _t.sleep(0.05)

    def __init__(self, n_workers: int, packing_callable: Callable,
                 throttle_local: bool = False,
                 worker_cooldown_seconds: float = 0.0):
        from concurrent.futures import ProcessPoolExecutor
        import multiprocessing as mp
        self.n_workers = n_workers
        self.packing_callable = packing_callable
        self.throttle_local = throttle_local
        self.worker_cooldown_seconds = worker_cooldown_seconds
        ctx = mp.get_context("fork")
        self.executor = ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx)
        self._future_to_data = {}   # future -> (attempt, start_wall, worker_id)
        self._wall_start = time.monotonic()
        self._busy_total = 0.0
        self._cooldown_until = [0.0] * n_workers   # per-worker wall-time
        self._next_worker_round = 0
        self._segments = []   # (worker_id, start_rel, end_rel, attempt_id, stage)

    def now(self):
        return time.monotonic() - self._wall_start

    def n_idle(self):
        # Available = not in-flight AND past cooldown
        in_flight_workers = set(
            data[2] for data in self._future_to_data.values()
        )
        wall_now = time.monotonic()
        avail = 0
        for w in range(self.n_workers):
            if w in in_flight_workers:
                continue
            if self.throttle_local and wall_now < self._cooldown_until[w]:
                continue
            avail += 1
        return avail

    def in_flight_count(self):
        return len(self._future_to_data)

    def _pick_idle_worker(self):
        in_flight = set(data[2] for data in self._future_to_data.values())
        wall_now = time.monotonic()
        for w in range(self.n_workers):
            if w in in_flight:
                continue
            if self.throttle_local and wall_now < self._cooldown_until[w]:
                continue
            return w
        return None

    def dispatch(self, attempt: Attempt):
        w = self._pick_idle_worker()
        assert w is not None, "RealRunner.dispatch called with no idle worker"
        attempt.worker_id = w
        attempt.started_at = self.now()
        start_wall = time.monotonic()
        f = self.executor.submit(self.packing_callable, attempt.theta,
                                  attempt.seed_id, attempt.stage.value)
        self._future_to_data[f] = (attempt, start_wall, w)

    def next_event(self) -> Optional[Event]:
        if not self._future_to_data:
            return None
        from concurrent.futures import as_completed
        # Block on next completion. Use a long timeout from the iterator.
        it = as_completed(list(self._future_to_data.keys()), timeout=None)
        f = next(it, None)
        if f is None:
            return None
        attempt, start_wall, worker_id = self._future_to_data.pop(f)
        end_wall = time.monotonic()
        t_now = end_wall - self._wall_start
        elapsed = end_wall - start_wall
        self._busy_total += elapsed
        self._segments.append((worker_id, start_wall - self._wall_start,
                                t_now, attempt.attempt_id, attempt.stage.value))
        # Apply cooldown
        if self.throttle_local and self.worker_cooldown_seconds > 0:
            self._cooldown_until[worker_id] = end_wall + self.worker_cooldown_seconds
        try:
            result = f.result()
        except Exception as exc:
            print(f"[heads-up] one packing process failed and was skipped "
                  f"(stage={attempt.stage.value}, seed={attempt.seed_id}). "
                  f"The error was:\n  {exc!r}", flush=True)
            return Event(fake_time=t_now, event_type="failed",
                          attempt_id=attempt.attempt_id)
        if result.get("success"):
            return Event(fake_time=t_now, event_type="completed",
                          attempt_id=attempt.attempt_id,
                          phi=float(result["phi"]))
        else:
            # The worker survived and reported a clean failure (e.g. the
            # engine raised on neighbor-capacity overflow). Skip this run,
            # tell the user what the engine said, and keep going.
            err = result.get("error", "no error message recorded")
            print(f"[heads-up] one packing run failed and was skipped "
                  f"(stage={attempt.stage.value}, seed={attempt.seed_id}). "
                  f"The engine said:\n  {err}", flush=True)
            return Event(fake_time=t_now, event_type="failed",
                          attempt_id=attempt.attempt_id)

    def worker_timeline(self):
        return self._segments

    def total_busy_time(self):
        return self._busy_total

    def shutdown(self):
        self.executor.shutdown(wait=False, cancel_futures=True)
