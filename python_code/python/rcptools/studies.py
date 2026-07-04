"""Pre-populated study menu.

Each study returns ``(model, theta_start, config_overrides)``. The notebook
combines the overrides with :func:`default_config` to produce a fully
populated CONFIG.

Conventions
-----------

* 2D studies use ``N=1225`` by default and `phi_init=0.6`.
* 3D studies use a smaller per-study ``N`` baseline because the C++ packer
  scales worse in memory/time at high SR and high N. Override as needed.
* The diameter distribution is built with :func:`make_universal_model`.
  Gaussian-only studies use ``Ng`` Gaussian components with random centers;
  hybrid studies use a power-law tail plus a single Gaussian; tri-mix uses
  a power-law plus a low-R Gaussian plus a high-R Gaussian.

Use :func:`list_studies` to discover available keys.
"""
from __future__ import annotations

import pathlib

import numpy as np

from .model import make_universal_model, set_perturb_weights


def default_config(*, search_name, save_dir, Ndim, N,
                   generations=16, population_size=84,
                   cheap_seeds_per_candidate=2,
                   elite_count=10, elite_extra_seeds=16,
                   initial_sigma=0.35, sigma_decay=0.88,
                   master_seed=1234,
                   neighbor_max=0,
                   phi_init=0.60, initializer_phi=0.11):
    """A fully populated CONFIG with reasonable defaults.

    Override any field at the call site or by ``CONFIG[k] = v`` after the
    fact. Every per-section knob (scheduler, pre-screen, time-aware,
    progress, packing-save) is included so the user has one dict to scan.
    """
    return {
        # --- Identity ---
        "search_name": str(search_name),
        "save_dir": str(save_dir),
        "resume_if_checkpoint_exists": True,

        # --- Packing ---
        "N": int(N),
        "Ndim": int(Ndim),
        "phi_init": float(phi_init),
        "initializer_phi": float(initializer_phi),
        "neighbor_max": int(neighbor_max),   # 0 = Ndim/SR-scaled default
        "fix_height": False,
        "box": None,
        "walls": None,

        # --- Search loop ---
        "generations": int(generations),
        "population_size": int(population_size),
        "cheap_seeds_per_candidate": int(cheap_seeds_per_candidate),
        "elite_count": int(elite_count),
        "elite_extra_seeds": int(elite_extra_seeds),
        "initial_sigma": float(initial_sigma),
        "sigma_decay": float(sigma_decay),
        "master_seed": int(master_seed),

        # --- Scoring ---
        "score_by": "mean_phi",   # or "phi_corr"

        # --- Parallel ---
        "use_parallel": True,
        "n_workers": "auto",
        "reserve_cpus": 1,

        # --- Scheduler ---
        "scheduler_mode": "adaptive",   # or "batch"
        "cheap_delta": 0.012,
        "phi_noise_floor": 0.004,
        "throttle_local": False,
        "worker_cooldown_seconds": 0.0,

        # --- Pre-screen (D_max geometric reject before C++ call) ---
        "prescreen_enabled": True,
        "prescreen_phi_target": 0.85,
        "prescreen_threshold": 2.5,    # D_max < L_min / threshold
        "prescreen_verbose": False,

        # --- Time-aware scheduler (straggler abandon, dynamic timeout) ---
        "time_aware_enabled": False,
        "time_aware_initial_timeout_s": 3600.0,
        "time_aware_calibration_threshold": 15,
        "time_aware_percentile": 95.0,
        "time_aware_headroom": 1.5,
        "time_aware_min_floor_fraction": 0.3,
        "time_aware_verbose": True,

        # --- Progress reports ---
        "progress_enabled": True,
        "progress_print_after_n_events": 25,
        "progress_heartbeat_every_seconds": 900.0,
        "progress_verbose_transitions": True,

        # --- Packing save ---
        "save_packings": "none",   # "none" / "best" / "all"

        # --- Packer convergence budget (sets RCP_SPEED env var per run) ---
        # "immediate" / "quick" / "patient" / "forever".
        # None = leave whatever the session has already set.
        "speed": "quick",
    }


# ============================================================================
# Study builders. Each returns (model, theta_start, config_overrides).
# ============================================================================
def _gaussian_only_model(Ng, size_ratio, Ndim, *, name):
    """Equally-spaced Gaussians in R, equal-volume balance."""
    centers = np.linspace(0.1, 0.9, Ng).tolist() if Ng > 1 else [0.5]
    components = [
        {"type": "gaussian", "init": {"R_center_u": float(c), "sigma_u": 0.05}}
        for c in centers
    ]
    rule = "equal_volume" if Ndim == 3 else "equal_area"
    MODEL = make_universal_model(
        components=components,
        size_ratio=float(size_ratio),
        Ndim=int(Ndim),
        alpha_tukey=0.05,
        balance_rule=rule,
        name=name,
    )
    theta = MODEL["init_rules"]["balanced"](MODEL)
    return MODEL, theta


def _powlaw_gauss_model(size_ratio, Ndim, *, name, exp=-2.3,
                        R_low_u=0.01, R_high_u=0.25,
                        gauss_R_u=0.90, gauss_sigma_u=0.05):
    """Power-law tail (the fines) plus one Gaussian (the coarse peak)."""
    components = [
        {"type": "power", "init": {
            "R_low_u": R_low_u, "R_high_u": R_high_u, "exp": exp}},
        {"type": "gaussian", "init": {
            "R_center_u": gauss_R_u, "sigma_u": gauss_sigma_u}},
    ]
    rule = "equal_volume" if Ndim == 3 else "equal_area"
    MODEL = make_universal_model(
        components=components,
        size_ratio=float(size_ratio),
        Ndim=int(Ndim),
        alpha_tukey=0.01,
        balance_rule=rule,
        name=name,
    )
    # Slow the power-law exponent perturbations — search shouldn't yank it
    # around step-by-step.
    set_perturb_weights(MODEL, {"power_0.exp": 0.05})
    theta = MODEL["init_rules"]["balanced"](MODEL)
    return MODEL, theta


def _trimix_model(size_ratio, Ndim, *, name, exp=-3.47):
    """Power-law fines + low-R Gaussian (mid) + high-R Gaussian (coarse)."""
    components = [
        {"type": "power", "init": {
            "R_low_u": 0.01, "R_high_u": 0.10, "exp": exp}},
        {"type": "gaussian", "init": {
            "R_center_u": 0.09, "sigma_u": 0.05}},
        {"type": "gaussian", "init": {
            "R_center_u": 0.90, "sigma_u": 0.05}},
    ]
    rule = "equal_volume" if Ndim == 3 else "equal_area"
    MODEL = make_universal_model(
        components=components,
        size_ratio=float(size_ratio),
        Ndim=int(Ndim),
        alpha_tukey=0.0025,
        balance_rule=rule,
        grid_n=10000,
        name=name,
    )
    set_perturb_weights(MODEL, {"power_0.exp": 0.05})
    theta = MODEL["init_rules"]["balanced"](MODEL)
    return MODEL, theta


# ============================================================================
# Study registry — pass each entry to build_study(...) to materialize it.
# ============================================================================
STUDIES_2D = {
    # --- Gaussian-only sweeps (Ng × SR) ---
    "gauss_Ng2_SR20":  {"family": "gaussian", "Ng": 2,  "size_ratio": 20},
    "gauss_Ng4_SR20":  {"family": "gaussian", "Ng": 4,  "size_ratio": 20},
    "gauss_Ng4_SR40":  {"family": "gaussian", "Ng": 4,  "size_ratio": 40},
    "gauss_Ng4_SR60":  {"family": "gaussian", "Ng": 4,  "size_ratio": 60},
    "gauss_Ng12_SR40": {"family": "gaussian", "Ng": 12, "size_ratio": 40},
    "gauss_Ng12_SR60": {"family": "gaussian", "Ng": 12, "size_ratio": 60},

    # --- power-law + Gaussian hybrid sweep ---
    "powlaw_gauss_SR40":  {"family": "powlaw_gauss", "size_ratio": 40},
    "powlaw_gauss_SR60":  {"family": "powlaw_gauss", "size_ratio": 60},
    "powlaw_gauss_SR80":  {"family": "powlaw_gauss", "size_ratio": 80},
    "powlaw_gauss_SR100": {"family": "powlaw_gauss", "size_ratio": 100},

    # --- tri-mix at extreme SR ---
    "trimix_SR150": {"family": "trimix", "size_ratio": 150},
}

STUDIES_3D = {
    "gauss_Ng2_SR20": {"family": "gaussian", "Ng": 2, "size_ratio": 20, "N": 30000},
    "gauss_Ng4_SR20": {"family": "gaussian", "Ng": 4, "size_ratio": 20, "N": 30000},
    "gauss_Ng4_SR30": {"family": "gaussian", "Ng": 4, "size_ratio": 30, "N": 30000},

    "powlaw_gauss_SR25": {"family": "powlaw_gauss", "size_ratio": 25, "N": 30000},
    "powlaw_gauss_SR30": {"family": "powlaw_gauss", "size_ratio": 30, "N": 30000},

    "trimix_SR30": {"family": "trimix", "size_ratio": 30, "N": 30000},
    "trimix_SR50": {"family": "trimix", "size_ratio": 50, "N": 30000},
}


def list_studies(Ndim=None):
    """Return the registry keyed by Ndim. Pass `Ndim=2` or `3` to filter."""
    if Ndim == 2:
        return dict(STUDIES_2D)
    if Ndim == 3:
        return dict(STUDIES_3D)
    return {"2D": dict(STUDIES_2D), "3D": dict(STUDIES_3D)}


def build_study(key, *, Ndim, save_root, N=None, generations=None,
                config_overrides=None):
    """Materialize a study by key.

    Returns ``(model, theta_start, config)``. ``config['save_dir']`` is set
    to ``save_root / <study_key>``.
    """
    registry = STUDIES_2D if Ndim == 2 else STUDIES_3D
    if key not in registry:
        raise KeyError(f"unknown study '{key}'; available {Ndim}D: {list(registry)}")
    spec = dict(registry[key])
    family = spec.pop("family")
    size_ratio = float(spec.pop("size_ratio"))
    N_baseline = int(spec.pop("N", 1225 if Ndim == 2 else 30000))
    if N is not None:
        N_baseline = int(N)

    name = f"{Ndim}d_{key}_sr{size_ratio:g}"

    if family == "gaussian":
        Ng = int(spec.pop("Ng"))
        model, theta = _gaussian_only_model(Ng, size_ratio, Ndim, name=name)
    elif family == "powlaw_gauss":
        model, theta = _powlaw_gauss_model(size_ratio, Ndim, name=name)
    elif family == "trimix":
        model, theta = _trimix_model(size_ratio, Ndim, name=name)
    else:
        raise ValueError(f"unknown family '{family}'")

    save_dir = pathlib.Path(save_root) / f"{Ndim}d_{key}"
    cfg = default_config(
        search_name=name,
        save_dir=str(save_dir),
        Ndim=Ndim,
        N=N_baseline,
        generations=int(generations) if generations is not None else 16,
    )
    if config_overrides:
        cfg.update(config_overrides)
    return model, theta, cfg
