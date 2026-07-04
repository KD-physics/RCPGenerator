"""Optional time-aware scheduler hardening + running progress reports.

Activates when ``CONFIG['time_aware_enabled']`` or ``CONFIG['progress_enabled']``
is true. Call :func:`apply` once per session (idempotent).

What it adds when active:

* **Per-packing dynamic timeout** with abandon-on-straggler semantics. The
  worker subprocess keeps running its slot until natural completion, but the
  scheduler stops waiting and treats the seed as failed. Compatible with the
  Python 3.12+ fork start method (no ``max_tasks_per_child`` required).
* **Decision-lock drain** — when the top confirmed elite's ``mean - SE``
  beats the runner-up's ``mean + SE``, the elite ranking is settled and every
  in-flight confirm seed is abandoned, ending the generation early.
* **Dropped-candidate drain** — as soon as a candidate hits DROPPED or
  PATHOLOGICAL, its in-flight sibling seeds are abandoned.
* **Running progress prints** every N scheduler events plus a periodic
  runner-side heartbeat during long quiet polls.
* **Phase-transition annotation** at the cheap→elite boundary.

All knobs live on the single CONFIG dict (``time_aware_*`` / ``progress_*``).
"""
from __future__ import annotations

import collections
import multiprocessing as _mp
import time as _time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from .scheduler import AdaptiveScheduler, CandidateState, Event, RealRunner, Stage

_PATCH_STATE = {
    "applied": False,
    "TIME_AWARE_CONFIG": None,
    "PROGRESS_CONFIG": None,
    "progress": {"last_print_wall": 0.0, "events_since_print": 0},
}


def _ta_cfg():
    return _PATCH_STATE["TIME_AWARE_CONFIG"] or {}


def _prog_cfg():
    return _PATCH_STATE["PROGRESS_CONFIG"] or {}


def _build_subconfigs(config):
    """Pull folded CONFIG keys into the legacy patch dicts."""
    return (
        {
            "enabled": bool(config.get("time_aware_enabled", False)),
            "initial_timeout_s": float(config.get("time_aware_initial_timeout_s", 3600.0)),
            "calibration_threshold": int(config.get("time_aware_calibration_threshold", 15)),
            "percentile": float(config.get("time_aware_percentile", 95.0)),
            "headroom": float(config.get("time_aware_headroom", 1.5)),
            "min_floor_fraction": float(config.get("time_aware_min_floor_fraction", 0.3)),
            "verbose": bool(config.get("time_aware_verbose", True)),
        },
        {
            "enabled": bool(config.get("progress_enabled", True)),
            "print_after_n_events": int(config.get("progress_print_after_n_events", 5)),
            "heartbeat_every_seconds": float(config.get("progress_heartbeat_every_seconds", 60.0)),
            "verbose_transitions": bool(config.get("progress_verbose_transitions", True)),
        },
    )


def apply(config):
    """Apply the patches. Re-callable to refresh config knobs."""
    ta, prog = _build_subconfigs(config)
    _PATCH_STATE["TIME_AWARE_CONFIG"] = ta
    _PATCH_STATE["PROGRESS_CONFIG"] = prog
    if _PATCH_STATE["applied"]:
        return

    # Preserve originals.
    if not hasattr(RealRunner, "_ta_original_init"):
        RealRunner._ta_original_init = RealRunner.__init__
        RealRunner._ta_original_next_event = RealRunner.next_event
    if not hasattr(AdaptiveScheduler, "_ta_original_handle_completed"):
        AdaptiveScheduler._ta_original_handle_completed = AdaptiveScheduler._handle_completed
        AdaptiveScheduler._ta_original_handle_failed = AdaptiveScheduler._handle_failed
        AdaptiveScheduler._ta_original_enter_confirm = AdaptiveScheduler._enter_confirm_phase

    RealRunner.__init__ = _ta_patched_init
    RealRunner.next_event = _ta_patched_next_event
    RealRunner._ta_update_timeout = _ta_update_timeout
    RealRunner._ta_heartbeat = _ta_heartbeat
    AdaptiveScheduler._handle_completed = _ta_patched_handle_completed
    AdaptiveScheduler._handle_failed = _ta_patched_handle_failed
    AdaptiveScheduler._enter_confirm_phase = _ta_patched_enter_confirm_phase

    _PATCH_STATE["applied"] = True
    print("[time-aware] patches applied.")


# ============================================================================
# RealRunner patches.
# ============================================================================
def _ta_patched_init(self, n_workers, packing_callable,
                     throttle_local=False, worker_cooldown_seconds=0.0):
    self.n_workers = n_workers
    self.packing_callable = packing_callable
    self.throttle_local = throttle_local
    self.worker_cooldown_seconds = worker_cooldown_seconds
    ctx = _mp.get_context("fork")
    self.executor = ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx)
    self._future_to_data = {}
    self._wall_start = _time.monotonic()
    self._busy_total = 0.0
    self._cooldown_until = [0.0] * n_workers
    self._next_worker_round = 0
    self._segments = []
    cfg = _ta_cfg()
    self._ta_packing_times = collections.deque(maxlen=50)
    self._ta_enabled = bool(cfg.get("enabled", False))
    self._ta_initial_timeout = float(cfg.get("initial_timeout_s", 3600.0))
    self._ta_calib_thresh = int(cfg.get("calibration_threshold", 15))
    self._ta_percentile = float(cfg.get("percentile", 95.0))
    self._ta_headroom = float(cfg.get("headroom", 1.5))
    self._ta_min_floor = float(cfg.get("min_floor_fraction", 0.3))
    self._ta_verbose = bool(cfg.get("verbose", True))
    self._ta_dynamic_timeout = self._ta_initial_timeout
    self._ta_n_timeouts = 0


def _ta_update_timeout(self):
    if len(self._ta_packing_times) >= self._ta_calib_thresh:
        t = np.array(self._ta_packing_times, dtype=float)
        p = float(np.percentile(t, self._ta_percentile))
        self._ta_dynamic_timeout = max(
            self._ta_initial_timeout * self._ta_min_floor,
            p * self._ta_headroom,
        )


def _ta_heartbeat(self):
    n = len(self._future_to_data)
    if n == 0:
        return
    now = _time.monotonic()
    cheap_n = sum(1 for (a, _, _) in self._future_to_data.values()
                  if getattr(a.stage, "value", str(a.stage)) == "cheap")
    confirm_n = n - cheap_n
    max_elapsed = max(now - sw for _, sw, _ in self._future_to_data.values())
    print(f"[heartbeat] in-flight={n} ({cheap_n} cheap, {confirm_n} elite), "
          f"max_elapsed={max_elapsed:.0f}s, "
          f"timeout_threshold={self._ta_dynamic_timeout:.0f}s, "
          f"timeouts_so_far={self._ta_n_timeouts}")
    _PATCH_STATE["progress"]["last_print_wall"] = now


def _ta_patched_next_event(self):
    if not self._future_to_data:
        return None
    if not self._ta_enabled:
        return RealRunner._ta_original_next_event(self)

    POLL_S = 1.0
    while self._future_to_data:
        now = _time.monotonic()
        timed_out = None
        for f, (attempt, start_wall, w) in self._future_to_data.items():
            if not f.done() and (now - start_wall) > self._ta_dynamic_timeout:
                timed_out = (f, attempt, start_wall, w)
                break
        if timed_out is not None:
            f, attempt, start_wall, worker_id = timed_out
            elapsed = now - start_wall
            del self._future_to_data[f]
            self._ta_n_timeouts += 1
            self._ta_packing_times.append(elapsed)
            self._ta_update_timeout()
            if self._ta_verbose:
                stage_val = getattr(attempt.stage, "value", str(attempt.stage))
                print(f"[time-aware] abandoned straggler attempt={attempt.attempt_id} "
                      f"seed={attempt.seed_id} stage={stage_val} elapsed={elapsed:.0f}s "
                      f"threshold={self._ta_dynamic_timeout:.0f}s")
            t_now = now - self._wall_start
            return Event(fake_time=t_now, event_type="failed",
                         attempt_id=attempt.attempt_id)

        done = [f for f in self._future_to_data if f.done()]
        if done:
            f = done[0]
            attempt, start_wall, worker_id = self._future_to_data.pop(f)
            end_wall = _time.monotonic()
            t_now = end_wall - self._wall_start
            elapsed = end_wall - start_wall
            self._busy_total += elapsed
            self._segments.append((worker_id, start_wall - self._wall_start,
                                   t_now, attempt.attempt_id, attempt.stage.value))
            if self.throttle_local and self.worker_cooldown_seconds > 0:
                self._cooldown_until[worker_id] = end_wall + self.worker_cooldown_seconds
            self._ta_packing_times.append(elapsed)
            self._ta_update_timeout()
            try:
                result = f.result()
            except Exception:
                return Event(fake_time=t_now, event_type="failed",
                             attempt_id=attempt.attempt_id)
            if result.get("success"):
                return Event(fake_time=t_now, event_type="completed",
                             attempt_id=attempt.attempt_id,
                             phi=float(result["phi"]))
            return Event(fake_time=t_now, event_type="failed",
                         attempt_id=attempt.attempt_id)

        prog = _prog_cfg()
        if prog.get("enabled", True):
            hb_interval = float(prog.get("heartbeat_every_seconds", 60.0))
            if (now - _PATCH_STATE["progress"]["last_print_wall"]) >= hb_interval:
                _ta_heartbeat(self)
        _time.sleep(POLL_S)
    return None


# ============================================================================
# AdaptiveScheduler patches.
# ============================================================================
def _ta_scheduler_progress_line(self):
    cands = self.candidates
    if not cands:
        return
    total = len(cands)
    n_cheap_done = sum(c.n_cheap_valid for c in cands.values())
    cheap_planned = total * int(self.cfg.cheap_seeds_per_candidate)
    n_dropped = sum(1 for c in cands.values() if c.state == CandidateState.DROPPED)
    n_path = sum(1 for c in cands.values() if c.state == CandidateState.PATHOLOGICAL)
    n_in_flight = self.runner.in_flight_count()

    if self.phase == Stage.CHEAP:
        means = [c.mean_phi for c in cands.values() if c.mean_phi is not None]
        best_str = f"{max(means):.4f}" if means else "n/a"
        thr = getattr(self.runner, "_ta_dynamic_timeout", None)
        thr_str = f" threshold={thr:.0f}s" if thr is not None else ""
        print(f"[progress] CHEAP {n_cheap_done}/{cheap_planned} done | "
              f"in-flight={n_in_flight} | dropped={n_dropped} pathological={n_path} | "
              f"best_so_far={best_str}{thr_str}")
    else:
        n_confirm_done = sum(c.n_confirm_valid for c in cands.values())
        confirm_planned = int(self.cfg.elite_count) * int(self.cfg.elite_extra_seeds)
        confirmed = sorted(
            [c for c in cands.values() if c.state == CandidateState.CONFIRMED],
            key=lambda c: -(c.mean_phi or -1e9),
        )
        top_str = run_str = "n/a"
        if confirmed:
            c = confirmed[0]
            top_str = f"{c.mean_phi:.4f}±{(c.se_phi or 0):.4f}(cand{c.candidate_id})"
        if len(confirmed) >= 2:
            c = confirmed[1]
            run_str = f"{c.mean_phi:.4f}±{(c.se_phi or 0):.4f}(cand{c.candidate_id})"
        print(f"[progress] ELITE {n_confirm_done}/{confirm_planned} done | "
              f"in-flight={n_in_flight} | top={top_str} | runner={run_str}")

    _PATCH_STATE["progress"]["last_print_wall"] = _time.monotonic()
    _PATCH_STATE["progress"]["events_since_print"] = 0


def _ta_after_event(self):
    prog = _prog_cfg()
    if not prog.get("enabled", True):
        return
    _PATCH_STATE["progress"]["events_since_print"] += 1
    if _PATCH_STATE["progress"]["events_since_print"] >= int(prog.get("print_after_n_events", 5)):
        _ta_scheduler_progress_line(self)


def _ta_abandon_for_candidate(scheduler, candidate_id):
    runner = scheduler.runner
    to_abandon = [f for f, (a, _, _) in runner._future_to_data.items()
                  if getattr(a, "candidate_id", None) == candidate_id]
    for f in to_abandon:
        del runner._future_to_data[f]
    return len(to_abandon)


def _ta_drain_dropped_candidates(scheduler):
    n_total = 0
    runner = scheduler.runner
    in_flight_cands = {getattr(a, "candidate_id", None)
                       for (a, _, _) in runner._future_to_data.values()}
    terminal = {CandidateState.DROPPED, CandidateState.PATHOLOGICAL}
    for cand_id in list(in_flight_cands):
        if cand_id is None:
            continue
        c = scheduler.candidates.get(cand_id)
        if c is None:
            continue
        if c.state in terminal:
            n = _ta_abandon_for_candidate(scheduler, cand_id)
            if n > 0:
                n_total += n
                if _prog_cfg().get("verbose_transitions", True):
                    print(f"[drain] candidate {cand_id} ({c.state.value}); "
                          f"abandoned {n} in-flight seeds")
    return n_total


def _ta_check_elite_decision_lock(scheduler):
    if scheduler.phase != Stage.CONFIRM:
        return 0
    confirmed = sorted(
        [c for c in scheduler.candidates.values()
         if c.state == CandidateState.CONFIRMED],
        key=lambda c: -(c.mean_phi if c.mean_phi is not None else -1e9),
    )
    if len(confirmed) < 2:
        return 0
    top, run = confirmed[0], confirmed[1]
    top_se = top.se_phi if top.se_phi is not None else 0.0
    run_se = run.se_phi if run.se_phi is not None else 0.0
    if top.mean_phi is None or run.mean_phi is None:
        return 0
    if (top.mean_phi - top_se) <= (run.mean_phi + run_se):
        return 0
    runner = scheduler.runner
    n_abandoned = len(runner._future_to_data)
    if n_abandoned == 0:
        return 0
    runner._future_to_data.clear()
    if _prog_cfg().get("verbose_transitions", True):
        print(f"[decision-lock] elite ranking settled; "
              f"abandoned {n_abandoned} in-flight elite seeds")
    return n_abandoned


def _ta_patched_handle_completed(self, ev):
    AdaptiveScheduler._ta_original_handle_completed(self, ev)
    _ta_after_event(self)
    _ta_drain_dropped_candidates(self)
    _ta_check_elite_decision_lock(self)


def _ta_patched_handle_failed(self, ev):
    AdaptiveScheduler._ta_original_handle_failed(self, ev)
    _ta_after_event(self)
    _ta_drain_dropped_candidates(self)
    _ta_check_elite_decision_lock(self)


def _ta_patched_enter_confirm_phase(self):
    prog = _prog_cfg()
    if prog.get("verbose_transitions", True):
        n_promoted = sum(1 for c in self.candidates.values()
                         if c.state == CandidateState.CHEAP_DECIDED)
        n_in_flight = self.runner.in_flight_count()
        n_dropped = sum(1 for c in self.candidates.values()
                        if c.state == CandidateState.DROPPED)
        print(f"[transition] CHEAP → ELITE: promoted={n_promoted} | "
              f"dropped={n_dropped} | cheap_in_flight_to_abandon={n_in_flight}")
    AdaptiveScheduler._ta_original_enter_confirm(self)
    if prog.get("verbose_transitions", True):
        confirm_planned = int(self.cfg.elite_count) * int(self.cfg.elite_extra_seeds)
        print(f"[transition] elite phase: {self.cfg.elite_count} elites × "
              f"{self.cfg.elite_extra_seeds} seeds = {confirm_planned} planned")
