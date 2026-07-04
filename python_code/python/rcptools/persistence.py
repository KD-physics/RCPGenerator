"""Checkpoint and sidecar persistence + the :class:`SearchState` accessor wrapper.

Each search directory holds:

* ``checkpoint.pkl``  — the live search state (pickled)
* ``model.pkl``       — the MODEL used to produce the state
* ``config.json``     — the CONFIG used to produce the state
* ``search_log.csv``  — flat per-candidate-per-stage log rows
* ``best_result.json``— summary of the global best candidate
* ``packings/``       — optional per-packing ``.npz`` files (only if
  ``CONFIG['save_packings']`` is "best" or "all")

The sidecars make a saved directory self-describing: opening one in §4
analysis reconstructs the search without the originating notebook.
"""
from __future__ import annotations

import datetime
import json
import pathlib
import pickle
from copy import deepcopy

import numpy as np
import pandas as pd


# ============================================================================
# Path helpers.
# ============================================================================
def search_dir(config) -> pathlib.Path:
    p = pathlib.Path(config["save_dir"])
    p.mkdir(parents=True, exist_ok=True)
    return p


def checkpoint_path(config) -> pathlib.Path:
    return search_dir(config) / "checkpoint.pkl"


def model_path(config) -> pathlib.Path:
    return search_dir(config) / "model.pkl"


def config_path(config) -> pathlib.Path:
    return search_dir(config) / "config.json"


def log_csv_path(config) -> pathlib.Path:
    return search_dir(config) / "search_log.csv"


def best_json_path(config) -> pathlib.Path:
    return search_dir(config) / "best_result.json"


# ============================================================================
# Sidecar I/O.
# ============================================================================
def _json_safe(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj


def save_config_sidecar(config):
    """Write CONFIG to ``config.json`` next to the checkpoint."""
    path = config_path(config)
    payload = _json_safe(config)
    with open(path, "w") as f:
        json.dump(payload, f, indent=2, default=str)
    return path


def save_model_sidecar(MODEL, config):
    """Pickle MODEL (function refs survive pickle round-trip in the same kernel)."""
    path = model_path(config)
    with open(path, "wb") as f:
        pickle.dump(MODEL, f)
    return path


def load_model_sidecar(save_dir):
    path = pathlib.Path(save_dir) / "model.pkl"
    if not path.exists():
        return None
    with open(path, "rb") as f:
        return pickle.load(f)


def load_config_sidecar(save_dir):
    path = pathlib.Path(save_dir) / "config.json"
    if not path.exists():
        return None
    with open(path) as f:
        return json.load(f)


# ============================================================================
# State checkpoint + log + best.
# ============================================================================
def save_checkpoint(state, config):
    path = checkpoint_path(config)
    with open(path, "wb") as f:
        pickle.dump(state, f)
    return path


def load_checkpoint_dict(config):
    path = checkpoint_path(config)
    if not path.exists():
        return None
    with open(path, "rb") as f:
        return pickle.load(f)


def delete_checkpoint(config):
    path = checkpoint_path(config)
    if path.exists():
        path.unlink()
        print("Deleted checkpoint:", path)


def save_search_log(log_rows, config):
    df = pd.DataFrame(log_rows)
    df.to_csv(log_csv_path(config), index=False)
    return df


def save_best_result(best_result, config):
    if best_result is None:
        return
    path = best_json_path(config)
    payload = {k: v for k, v in best_result.items() if k != "raw_rows"}
    payload = _json_safe(payload)
    with open(path, "w") as f:
        json.dump(payload, f, indent=2, default=str)


def result_to_log_row(result, generation, stage, sigma, config, MODEL):
    """Flatten one aggregated candidate result into a CSV-friendly row."""
    from .model import sample_diameters

    theta = np.asarray(result["theta"], dtype=float)
    D = sample_diameters(theta, MODEL, N=int(config["N"]))
    return {
        "timestamp": datetime.datetime.now().isoformat(),
        "search_name": config.get("search_name", "unnamed"),
        "model_name": MODEL.get("name", "unnamed"),
        "generation": int(generation),
        "stage": stage,
        "candidate": int(result["candidate"]),
        "sigma": float(sigma),
        "Ndim": int(config["Ndim"]),
        "N": int(config["N"]),
        "n_success": len(result["phis"]),
        "failures": int(result["failures"]),
        "rejects": int(result.get("rejects", 0)),
        "mean_phi": float(result["mean_phi"]),
        "median_phi": float(result["median_phi"]),
        "best_phi": float(result["best_phi"]),
        "std_phi": float(result["std_phi"]),
        "mean_phi_corr": float(result.get("mean_phi_corr", result["mean_phi"])),
        "mean_max_overlap": float(result.get("mean_max_overlap", 0.0)),
        "mean_force": float(result["mean_force"]),
        "mean_steps": float(result["mean_steps"]),
        "mean_max_min_dist": float(result["mean_max_min_dist"]),
        "seeds_json": json.dumps(result["seeds"]),
        "phis_json": json.dumps(result["phis"]),
        "phi_corrs_json": json.dumps(result.get("phi_corrs", [])),
        "max_overlaps_json": json.dumps(result.get("max_overlaps", [])),
        "steps_json": json.dumps(result["steps"]),
        "forces_json": json.dumps(result["forces"]),
        "max_min_dists_json": json.dumps(result["max_min_dists"]),
        "theta_json": json.dumps(theta.tolist()),
        "Dmin": float(D.min()),
        "Dmax": float(D.max()),
        "D_ratio": float(D.max() / D.min()),
    }


# ============================================================================
# SearchState — accessor wrapper around the dict-form state.
# ============================================================================
class SearchState:
    """Accessor wrapper around the dict produced by the search loop.

    Pass either a live dict or a `save_dir` to load from sidecars. All
    methods take a MODEL argument when the answer requires the PDF or the
    diameter sampler (e.g. regenerate_packing, distribution_at).

    The accessors all read the per-row CSV log (the durable record);
    `state["best_result"]` is used for `best_*` shortcuts.
    """

    def __init__(self, state=None, MODEL=None, config=None, save_dir=None):
        if state is not None:
            self._state = state
            self._model = MODEL
            self._config = config or state.get("config", {})
        elif save_dir is not None:
            self._config = load_config_sidecar(save_dir)
            self._model = load_model_sidecar(save_dir)
            self._state = load_checkpoint_dict({"save_dir": str(save_dir)})
            if self._state is None:
                raise FileNotFoundError(
                    f"no checkpoint.pkl found in {save_dir}"
                )
            if self._config is None:
                self._config = self._state.get("config", {})
        else:
            raise ValueError("Need state= or save_dir=")
        self._df_cache = None

    # -- properties -------------------------------------------------------
    @property
    def config(self):
        return self._config

    @property
    def model(self):
        return self._model

    @property
    def state(self):
        return self._state

    @property
    def save_dir(self):
        return pathlib.Path(self._config["save_dir"])

    @property
    def df(self):
        """The flattened log DataFrame (cached). Falls back to CSV on disk
        if the state dict's `log_rows` is empty."""
        if self._df_cache is not None:
            return self._df_cache
        rows = self._state.get("log_rows") or []
        if rows:
            self._df_cache = pd.DataFrame(rows)
        else:
            path = self.save_dir / "search_log.csv"
            if path.exists():
                self._df_cache = pd.read_csv(path)
            else:
                self._df_cache = pd.DataFrame()
        return self._df_cache

    @property
    def n_generations(self):
        df = self.df
        return 0 if df.empty else int(df["generation"].max()) + 1

    # -- accessors -------------------------------------------------------
    def best_theta(self):
        """The global best theta (across all generations)."""
        br = self._state.get("best_result")
        if br is None:
            return None
        return np.asarray(br["theta"], dtype=float)

    def best_theta_at(self, generation, by="mean_phi", stage="confirmed"):
        """Best theta within a single generation. `by` may be 'mean_phi' or
        'phi_corr'."""
        df = self.df
        if df.empty:
            return None
        sub = df[(df["generation"] == int(generation)) & (df["stage"] == stage)]
        if sub.empty:
            return None
        col = "mean_phi_corr" if by == "phi_corr" else "mean_phi"
        if col not in sub.columns:
            col = "mean_phi"
        row = sub.loc[sub[col].idxmax()]
        return np.asarray(json.loads(row["theta_json"]), dtype=float)

    def top_n_at(self, generation, n=5, by="mean_phi", stage="confirmed"):
        """Top-n candidates within a single generation."""
        df = self.df
        if df.empty:
            return pd.DataFrame()
        sub = df[(df["generation"] == int(generation)) & (df["stage"] == stage)]
        col = "mean_phi_corr" if by == "phi_corr" else "mean_phi"
        if col not in sub.columns:
            col = "mean_phi"
        return sub.sort_values(col, ascending=False).head(n)

    def all_thetas_at(self, generation, stage=None):
        """Every candidate theta (and its row metrics) in one generation."""
        df = self.df
        if df.empty:
            return pd.DataFrame()
        sub = df[df["generation"] == int(generation)]
        if stage is not None:
            sub = sub[sub["stage"] == stage]
        return sub

    def phi_stats_at(self, generation, stage="confirmed"):
        """Aggregate phi statistics for a single generation."""
        sub = self.all_thetas_at(generation, stage=stage)
        if sub.empty:
            return {}
        phis = np.concatenate([np.asarray(json.loads(s), dtype=float)
                               for s in sub["phis_json"] if isinstance(s, str)])
        if phis.size == 0:
            return {}
        return {
            "n_candidates": int(len(sub)),
            "n_seeds": int(phis.size),
            "mean": float(np.mean(phis)),
            "median": float(np.median(phis)),
            "std": float(np.std(phis)),
            "min": float(np.min(phis)),
            "max": float(np.max(phis)),
            "p10": float(np.quantile(phis, 0.10)),
            "p90": float(np.quantile(phis, 0.90)),
        }

    def phi_trajectory(self, stage="confirmed", by="mean_phi"):
        """Per-generation best score + per-generation min/max/percentiles."""
        df = self.df
        if df.empty:
            return pd.DataFrame()
        sub = df[df["stage"] == stage].copy()
        if sub.empty:
            return pd.DataFrame()
        col = "mean_phi_corr" if by == "phi_corr" else "mean_phi"
        if col not in sub.columns:
            col = "mean_phi"
        agg = sub.groupby("generation")[col].agg(
            best="max", mean="mean", std="std",
            p10=lambda v: float(np.quantile(v, 0.10)),
            p90=lambda v: float(np.quantile(v, 0.90)),
        ).reset_index()
        return agg

    def overlap_trajectory(self, stage="confirmed"):
        """Per-generation mean / max of `mean_max_overlap`."""
        df = self.df
        if df.empty or "mean_max_overlap" not in df.columns:
            return pd.DataFrame()
        sub = df[df["stage"] == stage]
        agg = sub.groupby("generation")["mean_max_overlap"].agg(
            mean="mean", max="max"
        ).reset_index()
        return agg

    # -- convenience -----------------------------------------------------
    def summary(self):
        cfg = self._config
        mdl = self._model
        br = self._state.get("best_result") or {}
        print(f"search:  {cfg.get('search_name', '?')}")
        print(f"  Ndim={cfg.get('Ndim')}  N={cfg.get('N')}  "
              f"generations={cfg.get('generations')}")
        if mdl is not None:
            print(f"model:   {mdl.get('name', '?')}  family={mdl.get('family')}  "
                  f"size_ratio={mdl.get('size_ratio')}  Ng/Nc={len(mdl.get('components', []))}")
        print(f"state:   next_generation={self._state.get('next_generation')}  "
              f"sigma={self._state.get('sigma'):.4f}")
        print(f"best:    mean_phi={br.get('mean_phi', float('nan')):.6f}  "
              f"best_phi={br.get('best_phi', float('nan')):.6f}")
        print(f"log:     {len(self.df)} rows")
        print(f"save_dir: {self.save_dir}")
