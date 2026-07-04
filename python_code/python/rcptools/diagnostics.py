"""Post-search visualizations operating on a :class:`SearchState`.

Plots are designed for the universal-mixture family produced by
:func:`rcpgenerator.search.model.make_universal_model`.
"""
from __future__ import annotations

import json

import numpy as np

from .jobs import create_packing
from .model import evaluate_model_pdf, sample_diameters
from .persistence import SearchState


def _to_state(state_or_dir):
    if isinstance(state_or_dir, SearchState):
        return state_or_dir
    if hasattr(state_or_dir, "get") and "config" in state_or_dir:
        return SearchState(state=state_or_dir,
                           MODEL=state_or_dir.get("MODEL_for_plotting"),
                           config=state_or_dir["config"])
    return SearchState(save_dir=state_or_dir)


def plot_distribution_from_theta(theta, MODEL, *, density_space="R", N=None,
                                 grid_n=None, bins=30, title=None,
                                 xlog=False, ylog=False, xlim=None, ylim=None):
    """Plot the PDF implied by `theta` plus a sampled-diameter histogram."""
    import matplotlib.pyplot as plt

    grid_n = int(grid_n) if grid_n is not None else int(MODEL.get("grid_n", 5000))
    N = 250 if N is None else int(N)
    R_min, R_max = MODEL["R_min"], MODEL["R_max"]
    R_grid = np.linspace(R_min, R_max, grid_n)
    pR_raw = evaluate_model_pdf(R_grid, theta, MODEL)
    Z = float(np.trapezoid(pR_raw, R_grid))
    pR = pR_raw / Z if Z > 0 else pR_raw
    D_grid = 2.0 * R_grid
    pD = pR / 2.0

    D_sample = sample_diameters(theta, MODEL, N=N)
    R_sample = D_sample / 2.0
    title = title or MODEL.get("name", "universal mixture")

    ds = str(density_space).lower()
    if ds == "r":
        x, y, samp = R_grid, pR, R_sample
        xlabel, ylabel, hist_x = "radius R", "P(R)", "sampled R"
    elif ds == "d":
        x, y, samp = D_grid, pD, D_sample
        xlabel, ylabel, hist_x = "diameter D", "P(D)", "sampled D"
    elif ds == "logd":
        y = R_grid * pR
        x, samp = R_grid, R_sample
        xlabel, ylabel, hist_x = "radius R", "P(logD) = R * P(R)", "sampled R"
    else:
        raise ValueError("density_space must be 'R', 'D', or 'logD'")

    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
    axes[0].plot(x, y)
    axes[0].set_xlabel(xlabel); axes[0].set_ylabel(ylabel); axes[0].set_title(title)
    axes[1].hist(samp, bins=bins, density=True)
    axes[1].set_xlabel(hist_x); axes[1].set_ylabel("density")
    axes[1].set_title("sampled diameters")
    for ax in axes:
        if xlog: ax.set_xscale("log")
        if ylog: ax.set_yscale("log")
        if xlim: ax.set_xlim(xlim)
        if ylim: ax.set_ylim(ylim)
    plt.tight_layout(); plt.show()


def plot_search_history(state, *, by="mean_phi"):
    """Per-generation candidate scores + best-so-far trajectory."""
    import matplotlib.pyplot as plt

    s = _to_state(state)
    df = s.df
    if df.empty:
        print("No log rows yet.")
        return
    col = "mean_phi_corr" if by == "phi_corr" else "mean_phi"
    if col not in df.columns:
        col = "mean_phi"

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    cheap = df[df["stage"] == "cheap"]
    confirmed = df[df["stage"] == "confirmed"]
    if not cheap.empty:
        axes[0].scatter(cheap["generation"], cheap[col], alpha=0.35, label="cheap")
    if not confirmed.empty:
        axes[0].scatter(confirmed["generation"], confirmed[col], alpha=0.8,
                        label="confirmed")
    axes[0].set_xlabel("generation"); axes[0].set_ylabel(col)
    axes[0].set_title(f"candidates by generation ({col})")
    axes[0].legend()

    if not confirmed.empty:
        best_by_gen = confirmed.groupby("generation")[col].max()
        axes[1].plot(best_by_gen.index, best_by_gen.values, marker="o")
        axes[1].set_xlabel("generation"); axes[1].set_ylabel(f"best {col}")
        axes[1].set_title("best confirmed per generation")
    plt.tight_layout(); plt.show()


def plot_distribution_evolution(state, *, stage="confirmed", by="mean_phi",
                                xlog=False, ylog=False, xlim=None, ylim=None):
    """Overlay the best-θ P(R) across generations."""
    import matplotlib.pyplot as plt

    s = _to_state(state)
    if s.model is None:
        raise ValueError("SearchState has no MODEL — pass MODEL when loading.")
    df = s.df
    sub = df[df["stage"] == stage]
    if sub.empty:
        print(f"No '{stage}' rows yet."); return
    col = "mean_phi_corr" if by == "phi_corr" else "mean_phi"
    if col not in sub.columns:
        col = "mean_phi"
    idx = sub.groupby("generation")[col].idxmax()
    best_gen = sub.loc[idx].sort_values("generation")
    MODEL = s.model
    R_grid = np.linspace(MODEL["R_min"], MODEL["R_max"],
                         MODEL.get("grid_n", 5000))

    plt.figure(figsize=(8, 4))
    for _, row in best_gen.iterrows():
        theta = np.asarray(json.loads(row["theta_json"]), dtype=float)
        pdf = evaluate_model_pdf(R_grid, theta, MODEL)
        Z = float(np.trapezoid(pdf, R_grid))
        if Z > 0:
            pdf = pdf / Z
        plt.plot(R_grid, pdf, alpha=0.7, label=f"gen {int(row['generation'])}")
    plt.xlabel("R"); plt.ylabel("P(R)")
    plt.title("Best distribution per generation")
    plt.legend(fontsize=8, ncol=2)
    if xlog: plt.xscale("log")
    if ylog: plt.yscale("log")
    if xlim: plt.xlim(xlim)
    if ylim: plt.ylim(ylim)
    plt.tight_layout(); plt.show()


def plot_population_at_generation(state, generation, *, stage="confirmed",
                                  xlog=False, ylog=False, xlim=None, ylim=None):
    """Overlay every candidate's P(R) within a single generation."""
    import matplotlib.pyplot as plt

    s = _to_state(state)
    if s.model is None:
        raise ValueError("SearchState has no MODEL.")
    sub = s.all_thetas_at(generation, stage=stage)
    if sub.empty:
        print(f"No '{stage}' rows at generation {generation}."); return
    MODEL = s.model
    R_grid = np.linspace(MODEL["R_min"], MODEL["R_max"],
                         MODEL.get("grid_n", 5000))
    plt.figure(figsize=(8, 4))
    for _, row in sub.iterrows():
        theta = np.asarray(json.loads(row["theta_json"]), dtype=float)
        pdf = evaluate_model_pdf(R_grid, theta, MODEL)
        Z = float(np.trapezoid(pdf, R_grid))
        if Z > 0:
            pdf = pdf / Z
        plt.plot(R_grid, pdf, alpha=0.4)
    plt.xlabel("R"); plt.ylabel("P(R)")
    plt.title(f"Population at generation {generation} (stage={stage})")
    if xlog: plt.xscale("log")
    if ylog: plt.yscale("log")
    if xlim: plt.xlim(xlim)
    if ylim: plt.ylim(ylim)
    plt.tight_layout(); plt.show()


def regenerate_packing(state, *, generation=None, seed=0,
                       theta=None, compute_voronoi=False,
                       compute_overlap=True):
    """Recreate a packing from a saved state.

    With ``theta=None`` uses the best-theta-at-generation (default: global
    best). Returns ``(packing, stats)`` from :func:`create_packing`.
    """
    s = _to_state(state)
    if s.model is None:
        raise ValueError("SearchState has no MODEL.")
    if theta is None:
        if generation is None:
            theta = s.best_theta()
        else:
            theta = s.best_theta_at(generation)
        if theta is None:
            raise ValueError("No theta available.")
    cfg = s.config
    return create_packing(
        theta, s.model,
        N=int(cfg["N"]),
        Ndim=int(cfg["Ndim"]),
        seed=int(seed),
        phi_init=float(cfg["phi_init"]),
        initializer_phi=float(cfg.get("initializer_phi", 0.11)),
        neighbor_max=int(cfg.get("neighbor_max", 0)),
        compute_voronoi=compute_voronoi,
        compute_overlap=compute_overlap,
        verbose=False,
    )


def compare_runs(states, *, by="mean_phi", labels=None):
    """Side-by-side phi trajectory for multiple runs."""
    import matplotlib.pyplot as plt

    if labels is None:
        labels = [s.config.get("search_name", str(i))
                  if hasattr(s, "config") else str(i)
                  for i, s in enumerate(states)]
    plt.figure(figsize=(8, 4))
    for s_in, lbl in zip(states, labels):
        s = _to_state(s_in)
        traj = s.phi_trajectory(by=by)
        if traj.empty:
            continue
        plt.plot(traj["generation"], traj["best"], marker="o", label=lbl)
    col = "mean_phi_corr" if by == "phi_corr" else "mean_phi"
    plt.xlabel("generation"); plt.ylabel(f"best {col}")
    plt.title("Best-per-generation across runs")
    plt.legend(fontsize=8)
    plt.tight_layout(); plt.show()


def run_all_diagnostics(state, *, plot_kwargs=None, by="mean_phi"):
    """Convenience: history + per-gen distribution evolution + best-θ plot."""
    if plot_kwargs is None:
        plot_kwargs = {}
    s = _to_state(state)
    plot_search_history(s, by=by)
    plot_distribution_evolution(s, by=by, **plot_kwargs)
    theta = s.best_theta()
    if theta is not None and s.model is not None:
        plot_distribution_from_theta(
            theta, s.model, N=int(s.config["N"]),
            title=f"best (mean_phi={s.state['best_result']['mean_phi']:.6f})",
            **plot_kwargs,
        )
    return s.df, theta
