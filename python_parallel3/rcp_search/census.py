"""Packing census — a lightweight, append-only log of EVERY packing the
search generates (not just the winners).

Motivation
----------
An rcp_search run evaluates thousands–tens of thousands of packings while
optimizing distribution shape, then discards everything but the best. But
each packing is a free (distribution -> phi) data point spanning a far
wider range of polydispersity / skewness than any deliberate sweep — and
the search's *losers* (extreme distributions that pack poorly, including
phi=0 rejects) are exactly the signal for theory-comparison work like the
SkewnessStudyWithEric study. This module captures all of it.

Design
------
* Parent-side, single writer: rows are written by the main process at the
  end of each generation from results the scheduler already collected. The
  workers, the pool, and the C++ engine are never touched, so this cannot
  affect the (carefully tuned) worker call / kill / respawn machinery.
* Append-per-generation JSONL: each completed generation is durably on
  disk, so an interrupted run loses nothing prior. One JSON object per
  packing per line; trivially loaded into pandas later.
* Sufficient-statistics, not derived quantities: we log the raw moments
  m1..m6 of the realized diameter distribution plus its min/max, and the
  full theta. Polydispersity delta, skewness S (in any convention —
  number- vs volume-weighted, radius vs diameter), and distribution
  family are CLOSED-FORM functions of these, computed on demand at
  analysis time. The system/distribution definitions live in the
  model.pkl / config.json sidecars already saved next to this file.

The realized moments are computed here from `sample_diameters`, which is
deterministic (quantile sampler), so they reproduce exactly the diameters
that were packed — and are recorded now so the record is independent of
any future change to the sampler code.
"""
from __future__ import annotations

import json
import pathlib

import numpy as np

from .model import sample_diameters


CENSUS_FILENAME = "packing_census.jsonl"
_N_MOMENTS = 6


def census_path(save_dir) -> pathlib.Path:
    return pathlib.Path(save_dir) / CENSUS_FILENAME


def realized_distribution_stats(theta, MODEL, N):
    """Raw moments m1..m6 of the realized (mean-normalized) diameter
    distribution, plus min/max. Recomputed from the deterministic quantile
    sampler, so identical to the diameters actually packed (the per-seed
    shuffle is moment-invariant)."""
    D = sample_diameters(np.asarray(theta, dtype=float), MODEL, N=int(N))
    D = np.asarray(D, dtype=float)
    out = {f"m{k}": float(np.mean(D ** k)) for k in range(1, _N_MOMENTS + 1)}
    out["d_min"] = float(D.min())
    out["d_max"] = float(D.max())
    return out


def append_generation(state, MODEL, res, config):
    """Append one row per packing produced this generation. Driven entirely
    by `res` (the QueueScheduler result) + the candidates' theta; reads
    nothing from the workers directly. No-op if config['census'] is False.

    A "packing" = one completed seed of one candidate. Every candidate that
    produced at least one result is logged, regardless of its final state
    (confirmed / dropped / demoted / evicted), and prescreen-rejected
    packings are logged too (phi == 0, rejected flag set). Only hard
    failures (worker crash / OOM, which produce no result) are absent.
    """
    if not config.get("census", True):
        return 0

    path = census_path(config["save_dir"])
    gen = int(state.get("next_generation", 0))
    N = int(config["N"])
    Ndim = int(config["Ndim"])
    search_name = config.get("search_name", "unnamed")
    model_name = MODEL.get("name", "unnamed") if MODEL else "unnamed"
    size_ratio = (MODEL.get("size_ratio") if MODEL else None)

    n_written = 0
    with open(path, "a") as f:
        for c in res.get("candidates", []):
            raw = getattr(c, "raw", None)
            if not raw:
                continue
            # distribution stats are per-candidate (same diameters for all
            # of its seeds) — compute once.
            stats = realized_distribution_stats(c.theta, MODEL, N)
            theta_list = [float(x) for x in np.asarray(c.theta, dtype=float)]
            state_str = c.state.value if hasattr(c.state, "value") else str(c.state)
            for ordinal, r in enumerate(raw):
                row = {
                    "search_name": search_name,
                    "model_name": model_name,
                    "generation": gen,
                    "candidate": int(c.cid),
                    "seed": int(r.get("seed", ordinal)),
                    "seed_ordinal": ordinal,
                    "state": state_str,
                    "rejected": bool(r.get("rejected", False)),
                    "N": N,
                    "Ndim": Ndim,
                    "size_ratio": size_ratio,
                    # outputs (as returned by the worker)
                    "phi": float(r.get("raw_phi", r.get("phi", 0.0))),
                    "phi_corr": float(r.get("phi_corr", 0.0)),
                    "max_overlap": float(r.get("max_overlap", 0.0)),
                    "steps": (int(r["steps"]) if r.get("steps") is not None else None),
                    "wall_s": (float(r["wall_s"]) if r.get("wall_s") is not None else None),
                    # realized distribution: raw moments + extremes
                    **stats,
                    # full reconstruction handle
                    "theta": theta_list,
                }
                f.write(json.dumps(row) + "\n")
                n_written += 1
    return n_written


def load_census(save_dir):
    """Load the census as a pandas DataFrame (one row per packing).

    Derived quantities (delta, skewness in your chosen convention, family
    code) are computed from the logged moments at analysis time — see
    `moments_to_delta_S` for the moment-based forms, or use `theta` +
    the saved MODEL for full-PDF-based definitions.
    """
    import pandas as pd
    rows = []
    p = census_path(save_dir)
    if not p.exists():
        return pd.DataFrame()
    with open(p) as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return pd.DataFrame(rows)


def moments_to_delta_S(df, weighting="number"):
    """Convenience: number- or volume-weighted polydispersity (delta = std/mean)
    and skewness (third standardized moment) from the logged raw moments.

    `weighting`:
      - "number": moments of D as logged (m1..m6)
      - "volume": moments under weight proportional to D**3, i.e.
        E_v[D^k] = m_{k+3} / m_3

    Returns a DataFrame with added `delta` and `skewness` columns. This is
    ONE convention; adjust to match the target study's exact definition
    (radius vs diameter, which weighting) — all are closed-form in m1..m6.
    """
    import pandas as pd
    out = df.copy()
    if weighting == "volume":
        def mk(k):  # volume-weighted raw moment E_v[D^k]
            return out[f"m{k + 3}"] / out["m3"]
        e1, e2, e3 = mk(1), mk(2), mk(3)
    else:
        e1, e2, e3 = out["m1"], out["m2"], out["m3"]
    var = e2 - e1 ** 2
    std = np.sqrt(var.clip(lower=0))
    out["delta"] = std / e1
    mu3 = e3 - 3 * e1 * e2 + 2 * e1 ** 3
    out["skewness"] = mu3 / std.pow(3).replace(0, np.nan)
    return out
