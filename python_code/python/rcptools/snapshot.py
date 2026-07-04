"""Realtime mid-generation snapshot — crash recovery for long searches.

Separate from ``checkpoint.pkl`` (which is written only at generation
boundaries). The snapshot captures the FULL state of the in-progress
generation and is overwritten after each worker finishes (and the next
theta is launched), so a Colab/VM bump costs only the compute of
whatever packings were mid-flight at the instant of the kill — every
completed (theta, phi) and every in-flight theta survives.

Used only when the search is resumed with ``restore=True``; otherwise it
is written but ignored (a normal resume picks up from the last completed
generation via the checkpoint).

Crash-safety: written with a temp file + atomic rename, so a kill during
a write can never corrupt the live snapshot — the previous complete one
survives until the new one is fully on disk. Single writer (the main
process), so no locking. The whole worker / pool / kill machinery is
untouched by this module.
"""
from __future__ import annotations

import os
import pickle
import pathlib


SNAPSHOT_FILENAME = "snapshot.pkl"


def snapshot_path(save_dir) -> pathlib.Path:
    return pathlib.Path(save_dir) / SNAPSHOT_FILENAME


def write_snapshot(save_dir, *, generation, sigma, thetas, candidates,
                   kills, config, model):
    """Atomically overwrite the snapshot with the current generation state.

    `candidates` is the scheduler's live candidate list. We store each
    candidate's completed results + raw dicts + state; in-flight seeds are
    implicitly recovered because the candidate's theta is in `thetas` and
    its unfinished seeds simply get re-dispatched on restore.
    """
    import numpy as np
    cand_state = {}
    for c in candidates:
        st = c.state.value if hasattr(c.state, "value") else str(c.state)
        cand_state[int(c.cid)] = {
            "results": list(c.results),
            "raw": list(c.raw),
            "state": st,
        }
    payload = {
        "generation": int(generation),
        "sigma": float(sigma),
        "thetas": [list(map(float, np.asarray(t, dtype=float))) for t in thetas],
        "candidates": cand_state,
        "kills": int(kills),
        "config": config,
        "model": model,
    }
    path = snapshot_path(save_dir)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "wb") as f:
        pickle.dump(payload, f)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)        # atomic on POSIX; survives a mid-write kill


def load_snapshot(save_dir):
    """Return the snapshot dict, or None if absent/unreadable (a torn
    final write degrades to None -> caller falls back to checkpoint)."""
    path = snapshot_path(save_dir)
    if not path.exists():
        return None
    try:
        with open(path, "rb") as f:
            return pickle.load(f)
    except Exception:
        return None


def clear_snapshot(save_dir):
    """Remove the snapshot at a clean generation boundary (the checkpoint
    now covers completed generations). Deletes exactly one file by literal
    path — never a directory, never a wildcard."""
    path = snapshot_path(save_dir)
    try:
        if path.exists():
            path.unlink()
    except Exception:
        pass


def restore_from_snapshot(save_dir):
    """Helper: return ``(CONFIG, MODEL)`` as they were when the in-progress
    generation started. Raises FileNotFoundError if there is no snapshot.
    """
    snap = load_snapshot(save_dir)
    if snap is None:
        raise FileNotFoundError(
            f"no readable {SNAPSHOT_FILENAME} in {save_dir}")
    return snap["config"], snap["model"]
