"""Universal mixture model for the diameter distribution.

A MODEL is built from a roster of components — `power`, `gaussian`, or
`delta` — via :func:`make_universal_model`. Each component contributes a
normalized PDF on physical radius R; component weights come from a softmax
over per-component amplitude slots in `theta`.

Public surface:

- :func:`make_universal_model`
- :func:`balance_amplitudes`
- :func:`set_perturb_weights`
- :func:`evaluate_model_pdf`
- :func:`sample_diameters`
- :func:`initialize_theta`
- :func:`perturb_theta`
- :func:`clip_theta`
- :data:`COMPONENT_REGISTRY`
"""
from __future__ import annotations

import numpy as np


# ============================================================================
# Tukey window applied at the framework level (PDF * window).
# ============================================================================
def tukey_window_on_interval(x, xmin, xmax, alpha=0.25):
    x = np.asarray(x)
    w = np.zeros_like(x, dtype=float)
    inside = (x >= xmin) & (x <= xmax)
    if alpha <= 0:
        w[inside] = 1.0
        return w
    L = xmax - xmin
    taper = alpha * L / 2.0
    flat_left = xmin + taper
    flat_right = xmax - taper
    flat = inside & (x >= flat_left) & (x <= flat_right)
    w[flat] = 1.0
    left = inside & (x < flat_left)
    if np.any(left):
        w[left] = 0.5 * (1.0 - np.cos(np.pi * (x[left] - xmin) / taper))
    right = inside & (x > flat_right)
    if np.any(right):
        w[right] = 0.5 * (1.0 - np.cos(np.pi * (xmax - x[right]) / taper))
    return w


# ============================================================================
# Component registry. Each entry registers the param list + bounds + a
# raw (un-normalized) PDF function on physical R.
# ============================================================================
def _power_pdf_R(R, params, MODEL):
    """Power-law on a normalized [R_low_u, R_high_u] sub-interval."""
    R_min, R_max = MODEL["R_min"], MODEL["R_max"]
    span = R_max - R_min
    R_low = R_min + float(params[1]) * span
    R_high = R_min + float(params[2]) * span
    if R_high < R_low:
        R_low, R_high = R_high, R_low
    R = np.asarray(R, dtype=float)
    out = np.zeros_like(R)
    inside = (R >= R_low) & (R <= R_high) & (R > 0)
    if np.any(inside):
        out[inside] = R[inside] ** float(params[3])
    return out


def _gaussian_pdf_R(R, params, MODEL):
    """Gaussian centered at R_center_u with width sigma_u (in [0,1] of span)."""
    R_min, R_max = MODEL["R_min"], MODEL["R_max"]
    span = R_max - R_min
    R_c = R_min + float(params[1]) * span
    sigma = max(float(params[2]) * span, 1e-9)
    R = np.asarray(R, dtype=float)
    return np.exp(-0.5 * ((R - R_c) / sigma) ** 2)


def _delta_pdf_R(R, params, MODEL):
    """Narrow spike at R_center_u (fixed sub-grid width = 0.005 * span)."""
    R_min, R_max = MODEL["R_min"], MODEL["R_max"]
    span = R_max - R_min
    R_c = R_min + float(params[1]) * span
    sigma = 0.005 * span
    R = np.asarray(R, dtype=float)
    return np.exp(-0.5 * ((R - R_c) / sigma) ** 2)


def _lognormal_pdf_R(R, params, MODEL):
    """Truncated lognormal in DIAMETER, expressed on R (D = 2R).

    A lognormal of log-width ``alpha`` (std of ln D) whose GEOMETRIC MEAN
    (mu = exp(mu_log) = median) sits at a center placed LINEARLY across the
    span — identical convention to the gaussian's ``R_center_u``:

        R_c = R_min + center * (R_max - R_min)        # D_c = 2 R_c = mu

    Truncated symmetrically in log about D_c to an effective size ratio
    ``S_eff = cutoff * S`` (S = MODEL['size_ratio']), i.e. the local support
    is [D_c * S_eff^-1/2, D_c * S_eff^+1/2]. Anything past the global span
    [R_min, R_max] is clipped by the framework grid (same as the gaussian
    when its center/width overrun the span). The returned density is
    un-normalized (the framework normalizes over the grid); the 1/D factor
    becomes 1/R up to a constant, which normalization absorbs.

    params = [amp, alpha, cutoff, center]  (amp handled by the mixer).
    """
    R_min, R_max = MODEL["R_min"], MODEL["R_max"]
    span = R_max - R_min
    alpha = max(float(params[1]), 1e-6)
    cutoff = min(max(float(params[2]), 1e-6), 1.0)
    center = min(max(float(params[3]), 0.0), 1.0)
    S = float(MODEL.get("size_ratio", 1.0))
    R_c = R_min + center * span                     # geometric mean (in R)
    S_eff = max(cutoff * S, 1.0 + 1e-12)            # local size ratio >= 1
    h = 0.5 * np.log(S_eff)                          # half-width in log
    R_lo, R_hi = R_c * np.exp(-h), R_c * np.exp(h)
    R = np.asarray(R, dtype=float)
    out = np.zeros_like(R)
    inside = (R >= R_lo) & (R <= R_hi) & (R > 0)
    if np.any(inside):
        x = np.log(R[inside] / R_c)                  # = ln(D/D_c)
        out[inside] = np.exp(-0.5 * (x / alpha) ** 2) / R[inside]
    return out


COMPONENT_REGISTRY = {
    "power": {
        "n_params": 4,
        "param_names": ["amp", "R_low_u", "R_high_u", "exp"],
        "bounds_lower": [-30.0, 0.0, 0.0, -6.0],
        "bounds_upper": [30.0, 1.0, 1.0, 12.0],
        "perturb_scale": [1.0, 0.10, 0.10, 0.05],
        "defaults": {"amp": 0.0, "R_low_u": 0.0, "R_high_u": 1.0, "exp": 0.0},
        "pdf_fn": _power_pdf_R,
    },
    "gaussian": {
        "n_params": 3,
        "param_names": ["amp", "R_center_u", "sigma_u"],
        "bounds_lower": [-30.0, 0.0, 0.001],
        "bounds_upper": [30.0, 1.0, 0.5],
        "perturb_scale": [1.0, 0.10, 0.05],
        "defaults": {"amp": 0.0, "R_center_u": 0.5, "sigma_u": 0.05},
        "pdf_fn": _gaussian_pdf_R,
    },
    "delta": {
        "n_params": 2,
        "param_names": ["amp", "R_center_u"],
        "bounds_lower": [-30.0, 0.0],
        "bounds_upper": [30.0, 1.0],
        "perturb_scale": [1.0, 0.05],
        "defaults": {"amp": 0.0, "R_center_u": 0.5},
        "pdf_fn": _delta_pdf_R,
    },
    "lognormal": {
        # Truncated lognormal in diameter. shape params:
        #   alpha  — log-width (std of ln D)
        #   cutoff — effective size ratio as a LINEAR fraction of the global
        #            S: S_eff = cutoff * size_ratio  (cutoff=1 -> full span)
        #   center — geometric mean (mu) position, LINEAR across the span
        #            (same convention as gaussian R_center_u; default 0.5 is
        #            the linear midpoint, NOT D_g=1 — see lognormal_diagnostics)
        "n_params": 4,
        "param_names": ["amp", "alpha", "cutoff", "center"],
        "bounds_lower": [-30.0, 0.01, 0.05, 0.0],
        "bounds_upper": [30.0, 3.0, 1.0, 1.0],
        "perturb_scale": [1.0, 0.05, 0.05, 0.05],
        "defaults": {"amp": 0.0, "alpha": 0.9, "cutoff": 1.0, "center": 0.5},
        "pdf_fn": _lognormal_pdf_R,
    },
}


# ============================================================================
# Internal helpers used by the wrapper PDF and balance routines.
# ============================================================================
def _u_component_params(theta, MODEL, k):
    c = MODEL["components"][k]
    return np.asarray(theta[c["slice_start"]:c["slice_end"]], dtype=float)


def _u_norm_grid(MODEL):
    return np.linspace(MODEL["R_min"], MODEL["R_max"], int(MODEL.get("grid_n", 5000)))


def _u_component_moment(MODEL, k, params, q):
    """Per-component normalized moment ⟨R^q⟩ used by balance_amplitudes."""
    R = _u_norm_grid(MODEL)
    pdf = COMPONENT_REGISTRY[MODEL["components"][k]["type"]]["pdf_fn"](R, params, MODEL)
    Z = float(np.trapezoid(pdf, R))
    if Z < 1e-300:
        return 0.0
    return float(np.trapezoid((R ** q) * pdf, R) / Z)


# ============================================================================
# Universal-mixture PDF on physical R (no Tukey here; framework applies it).
# ============================================================================
def universal_pdf_R(R, theta, MODEL):
    R = np.asarray(R, dtype=float)
    theta = np.asarray(theta, dtype=float)
    components = MODEL["components"]
    amps = np.array([theta[c["slice_start"]] for c in components], dtype=float)
    a_shift = amps - np.max(amps)
    w = np.exp(a_shift)
    w_sum = float(np.sum(w))
    if w_sum <= 0.0 or not np.isfinite(w_sum):
        return np.zeros_like(R)
    w = w / w_sum

    R_norm = _u_norm_grid(MODEL)
    total = np.zeros_like(R)
    for k, c in enumerate(components):
        params = _u_component_params(theta, MODEL, k)
        fn = COMPONENT_REGISTRY[c["type"]]["pdf_fn"]
        norm_pdf = fn(R_norm, params, MODEL)
        Z = float(np.trapezoid(norm_pdf, R_norm))
        if Z < 1e-300:
            continue
        total = total + (w[k] / Z) * fn(R, params, MODEL)
    return total


def _universal_perturb_scale_func(theta, MODEL):
    """Theta-step scale = MODEL['theta_perturb_scale'] * MODEL['theta_perturb_weights']."""
    base = np.asarray(MODEL["theta_perturb_scale"], dtype=float)
    weights = MODEL.get("theta_perturb_weights")
    if weights is None:
        return base
    w = np.asarray(weights, dtype=float)
    if w.shape != base.shape:
        raise ValueError(
            f"theta_perturb_weights shape {w.shape} != theta_dim {base.shape}"
        )
    return base * w


# ============================================================================
# Framework primitives — model-agnostic; the search loop uses these.
# ============================================================================
def evaluate_model_pdf(R, theta, MODEL):
    """User PDF * Tukey window. This is what samplers and plots use."""
    R = np.asarray(R, dtype=float)
    raw = MODEL["pdf"](R, theta, MODEL)
    win = tukey_window_on_interval(R, MODEL["R_min"], MODEL["R_max"],
                                   alpha=MODEL.get("alpha_tukey", 0.25))
    return np.maximum(raw * win, 0.0)


def sample_diameters(theta, MODEL, N, rng=None):
    """Inverse-CDF sampling of N diameters. Deterministic quantile by default;
    pass `rng` for random uniform draws. Mean-normalized for numerical
    consistency with the search loop's downstream phi rescaling."""
    R_min = MODEL["R_min"]
    R_max = MODEL["R_max"]
    grid_n = MODEL.get("grid_n", 5000)
    R_grid = np.linspace(R_min, R_max, grid_n)
    pdf = evaluate_model_pdf(R_grid, theta, MODEL) + 1e-300
    pdf = pdf / np.trapezoid(pdf, R_grid)
    dx = R_grid[1] - R_grid[0]
    cdf = np.cumsum(pdf) * dx
    cdf = cdf / cdf[-1]
    if rng is None:
        u = (np.arange(N) + 0.5) / N
    else:
        rng_obj = np.random.default_rng(rng) if not isinstance(rng, np.random.Generator) else rng
        u = rng_obj.uniform(0, 1, N)
    R_sample = np.interp(u, cdf, R_grid)
    D_sample = 2.0 * R_sample
    D_sample = D_sample / np.mean(D_sample)
    return D_sample


def clip_theta(theta, MODEL):
    lo = MODEL["theta_bounds_lower"]
    hi = MODEL["theta_bounds_upper"]
    return np.clip(np.asarray(theta, dtype=float).copy(), lo, hi)


def perturb_theta(theta, sigma, MODEL, rng=None):
    rng = np.random.default_rng(rng) if not isinstance(rng, np.random.Generator) else rng
    if MODEL.get("perturb_scale_func") is not None:
        scale = MODEL["perturb_scale_func"](theta, MODEL)
    else:
        scale = MODEL.get("theta_perturb_scale", np.ones(MODEL["theta_dim"]))
    step = rng.normal(0.0, sigma, size=MODEL["theta_dim"]) * np.asarray(scale, dtype=float)
    return clip_theta(np.asarray(theta, dtype=float) + step, MODEL)


def init_theta_random(MODEL, rng=None, **_):
    rng = np.random.default_rng(rng) if not isinstance(rng, np.random.Generator) else rng
    lo = MODEL["theta_bounds_lower"]
    hi = MODEL["theta_bounds_upper"]
    return clip_theta(rng.uniform(lo, hi), MODEL)


# ============================================================================
# Module-level init rules (picklable). The factory binds these into
# MODEL['init_rules'] so they survive a pickle round-trip.
# ============================================================================
def init_theta_from_spec(MODEL, rng=None, **_):
    """Initialize theta from each component's `init` dict, defaulting any
    missing parameter to the registry's `defaults`."""
    theta = np.zeros(MODEL["theta_dim"], dtype=float)
    for c in MODEL["components"]:
        reg = COMPONENT_REGISTRY[c["type"]]
        for j, pname in enumerate(reg["param_names"]):
            if pname in c["init"]:
                theta[c["slice_start"] + j] = float(c["init"][pname])
            else:
                theta[c["slice_start"] + j] = float(reg["defaults"][pname])
    return clip_theta(theta, MODEL)


def init_theta_balanced(MODEL, rng=None, **_):
    """`init_theta_from_spec` then apply :func:`balance_amplitudes`."""
    theta = init_theta_from_spec(MODEL, rng=rng)
    return balance_amplitudes(theta, MODEL, rule=MODEL.get("balance_rule"))


def initialize_theta(MODEL, rule="balanced", rng=None, **rule_kwargs):
    """Dispatch to a registered init rule. 'random' always works; other rules
    are pulled from `MODEL['init_rules']` (factory-provided)."""
    if rule == "random":
        return init_theta_random(MODEL, rng=rng)
    rules = MODEL.get("init_rules", {})
    if rule in rules:
        return rules[rule](MODEL, rng=rng, **rule_kwargs)
    raise ValueError(
        f"unknown rule '{rule}'; MODEL supports ['random'] + {list(rules)}"
    )


# ============================================================================
# Balance + perturb-weight helpers (component-aware).
# ============================================================================
def balance_amplitudes(theta, MODEL, rule=None, fractions=None):
    """Set per-component amp slots so the mixture matches a balance rule.

    rules:
      'equal_count'  — each component contributes 1/K of the mass
      'equal_area'   — weights ∝ 1/⟨R^2⟩_k   (use for Ndim=2)
      'equal_volume' — weights ∝ 1/⟨R^3⟩_k   (use for Ndim=3)
      'none'         — leave amp slots unchanged

    `fractions` overrides `rule` with explicit target mass fractions.
    """
    if rule is None:
        rule = MODEL.get("balance_rule", "equal_count")
    theta = np.asarray(theta, dtype=float).copy()
    K = len(MODEL["components"])

    if fractions is None:
        if rule == "none":
            return clip_theta(theta, MODEL)
        if rule == "equal_count":
            fractions = np.full(K, 1.0 / K)
        elif rule in ("equal_area", "equal_volume"):
            q = 2 if rule == "equal_area" else 3
            inv = np.zeros(K, dtype=float)
            for k in range(K):
                params = _u_component_params(theta, MODEL, k)
                mq = _u_component_moment(MODEL, k, params, q)
                inv[k] = 1.0 / mq if mq > 1e-300 else 0.0
            tot = float(np.sum(inv))
            fractions = inv / tot if tot > 0 else np.full(K, 1.0 / K)
        else:
            raise ValueError(f"unknown balance rule '{rule}'")

    fractions = np.asarray(fractions, dtype=float)
    if fractions.shape[0] != K:
        raise ValueError(
            f"fractions length {fractions.shape[0]} != n_components {K}"
        )
    fractions = np.clip(fractions, 1e-300, None)
    fractions = fractions / float(np.sum(fractions))

    log_w = np.log(fractions)
    log_w = log_w - float(np.max(log_w))
    for k, c in enumerate(MODEL["components"]):
        theta[c["slice_start"]] = log_w[k]
    return clip_theta(theta, MODEL)


def set_perturb_weights(MODEL, spec):
    """Lock or rescale per-parameter perturbation steps.

    `spec` forms (anything unset stays at 1.0):
      - flat array of length theta_dim
      - list of per-component dicts: [{'exp': 0.0}, {'sigma_u': 0.1}, {}]
      - dict keyed by '<label>.<param>': {'power_0.exp': 0.0, 'gaussian_0.sigma_u': 0.1}

    Weight 0 = locked; 0.1 = slow drift; >1 = accelerated.
    """
    w = np.ones(MODEL["theta_dim"], dtype=float)

    if isinstance(spec, np.ndarray):
        arr = np.asarray(spec, dtype=float)
        if arr.shape != (MODEL["theta_dim"],):
            raise ValueError(
                f"flat weights length {arr.shape[0]} != theta_dim {MODEL['theta_dim']}"
            )
        w = arr.copy()
    elif isinstance(spec, (list, tuple)) and all(np.isscalar(x) for x in spec):
        arr = np.asarray(spec, dtype=float)
        if arr.shape != (MODEL["theta_dim"],):
            raise ValueError(
                f"flat weights length {arr.shape[0]} != theta_dim {MODEL['theta_dim']}"
            )
        w = arr.copy()
    elif isinstance(spec, (list, tuple)):
        if len(spec) != len(MODEL["components"]):
            raise ValueError(
                f"per-component list length {len(spec)} != "
                f"n_components {len(MODEL['components'])}"
            )
        for k, sub in enumerate(spec):
            c = MODEL["components"][k]
            reg = COMPONENT_REGISTRY[c["type"]]
            sub = dict(sub) if sub else {}
            for j, pname in enumerate(reg["param_names"]):
                if pname in sub:
                    w[c["slice_start"] + j] = float(sub[pname])
    elif isinstance(spec, dict):
        pos = {}
        for c in MODEL["components"]:
            reg = COMPONENT_REGISTRY[c["type"]]
            for j, pname in enumerate(reg["param_names"]):
                pos[(c["label"], pname)] = c["slice_start"] + j
                pos[f"{c['label']}.{pname}"] = c["slice_start"] + j
        for key, val in spec.items():
            if key not in pos:
                avail = [c["label"] for c in MODEL["components"]]
                raise KeyError(
                    f"unknown key '{key}'; component labels: {avail}"
                )
            w[pos[key]] = float(val)
    else:
        raise TypeError(f"unsupported spec type {type(spec)}")

    MODEL["theta_perturb_weights"] = w
    return w


# ============================================================================
# Factory.
# ============================================================================
def make_universal_model(components, *, size_ratio=5.0, Ndim=2,
                         alpha_tukey=0.25, grid_n=5000,
                         balance_rule="equal_count", name=None):
    """Build a MODEL dict from a roster spec.

    Parameters
    ----------
    components : list of {'type': str, 'init': dict}
        Each entry's 'type' must be a key of :data:`COMPONENT_REGISTRY`. The
        'init' dict overrides per-parameter defaults (see the registry).
    size_ratio : float
        Dmax / Dmin. The model's R range spans Dmin/2 to Dmax/2 with
        Dmin = 1/sqrt(size_ratio) and Dmax = sqrt(size_ratio).
    Ndim : int
        Dimension passed to balance_amplitudes when the search calls
        balance_rule='equal_area' or 'equal_volume'.
    alpha_tukey : float
        Tukey-window edge fraction.
    grid_n : int
        Grid resolution for numerical normalization and the CDF sampler.
    balance_rule : str
        Default rule used by the 'balanced' init rule and bare
        :func:`balance_amplitudes` calls.
    name : str | None
        Optional human label used in checkpoints / logs.
    """
    Dmin = 1.0 / np.sqrt(size_ratio)
    Dmax = np.sqrt(size_ratio)
    R_min = Dmin / 2.0
    R_max = Dmax / 2.0

    component_records = []
    bounds_lower, bounds_upper, perturb_scale = [], [], []
    theta_names, component_labels = [], []
    type_counts = {}
    offset = 0
    for spec in components:
        t = spec["type"]
        if t not in COMPONENT_REGISTRY:
            raise ValueError(
                f"unknown component type '{t}'; available: "
                f"{list(COMPONENT_REGISTRY)}"
            )
        reg = COMPONENT_REGISTRY[t]
        n = int(reg["n_params"])
        idx_of_type = type_counts.get(t, 0)
        label = f"{t}_{idx_of_type}"
        type_counts[t] = idx_of_type + 1
        component_records.append({
            "type": t,
            "label": label,
            "slice_start": offset,
            "slice_end": offset + n,
            "param_names": list(reg["param_names"]),
            "init": dict(spec.get("init", {})),
        })
        component_labels.append(label)
        bounds_lower.extend(reg["bounds_lower"])
        bounds_upper.extend(reg["bounds_upper"])
        perturb_scale.extend(reg["perturb_scale"])
        theta_names.extend(f"{label}.{p}" for p in reg["param_names"])
        offset += n

    theta_dim = offset
    bounds_lower = np.asarray(bounds_lower, dtype=float)
    bounds_upper = np.asarray(bounds_upper, dtype=float)
    perturb_scale = np.asarray(perturb_scale, dtype=float)
    theta_perturb_weights = np.ones(theta_dim, dtype=float)

    return {
        "name": name or "universal_" + "_".join(component_labels),
        "family": "universal_mixture",
        "pdf": universal_pdf_R,
        "components": component_records,
        "component_labels": component_labels,
        "size_ratio": float(size_ratio),
        "distribution_space": "R",
        "alpha_tukey": float(alpha_tukey),
        "R_min": float(R_min),
        "R_max": float(R_max),
        "grid_n": int(grid_n),
        "Ndim": int(Ndim),
        "balance_rule": str(balance_rule),
        "theta_dim": int(theta_dim),
        "theta_names": theta_names,
        "theta_bounds_lower": bounds_lower,
        "theta_bounds_upper": bounds_upper,
        "theta_perturb_scale": perturb_scale,
        "theta_perturb_weights": theta_perturb_weights,
        "perturb_scale_func": _universal_perturb_scale_func,
        "init_rules": {
            "from_spec": init_theta_from_spec,
            "balanced": init_theta_balanced,
            "random": init_theta_random,
        },
    }


def theta_to_dataframe(theta, MODEL):
    """Decode a universal-mixture theta into a physical-parameter table.

    Returns a long-format DataFrame, one row per (component, parameter),
    with the raw theta value, its physical value, and units. The physical
    transforms mirror the model's own PDF code exactly:
      - amp slots -> softmax mixture weights across components
      - R_low_u / R_high_u / R_center_u -> R = R_min + raw * (R_max-R_min)
      - sigma_u -> raw * (R_max - R_min)
      - exp -> raw (1-to-1, no transform)
    """
    import pandas as pd
    theta = np.asarray(theta, dtype=float)
    R_min, R_max = MODEL["R_min"], MODEL["R_max"]
    span = R_max - R_min

    amps = np.array([theta[c["slice_start"]] for c in MODEL["components"]])
    w = np.exp(amps - amps.max())
    w = w / w.sum()

    rows = []
    for k, c in enumerate(MODEL["components"]):
        params = theta[c["slice_start"]:c["slice_end"]]
        label, ctype = c["label"], c["type"]
        rows.append({"component": label, "type": ctype, "param": "mixture_weight",
                     "raw_theta": float(amps[k]), "physical": float(w[k]),
                     "units": "fraction of mass"})
        for j, pname in enumerate(c["param_names"]):
            raw = float(params[j])
            if pname == "amp":
                continue
            elif pname in ("R_low_u", "R_high_u", "R_center_u"):
                R = R_min + raw * span
                rows.append({"component": label, "type": ctype, "param": pname,
                             "raw_theta": raw, "physical": R,
                             "units": f"R (D = {2 * R:.4f})"})
            elif pname == "sigma_u":
                rows.append({"component": label, "type": ctype, "param": pname,
                             "raw_theta": raw, "physical": raw * span,
                             "units": "R (width)"})
            elif pname == "exp":
                rows.append({"component": label, "type": ctype, "param": pname,
                             "raw_theta": raw, "physical": raw,
                             "units": "P(R) ∝ R^exp"})
            elif pname == "alpha":
                rows.append({"component": label, "type": ctype, "param": pname,
                             "raw_theta": raw, "physical": raw,
                             "units": "log-width (std of ln D)"})
            elif pname == "cutoff":
                S = float(MODEL.get("size_ratio", 1.0))
                rows.append({"component": label, "type": ctype, "param": pname,
                             "raw_theta": raw, "physical": raw * S,
                             "units": f"S_eff = cutoff*S (S={S:g})"})
            elif pname == "center":   # geometric mean (mu), linear across span
                R_c = R_min + raw * span
                rows.append({"component": label, "type": ctype, "param": pname,
                             "raw_theta": raw, "physical": 2.0 * R_c,
                             "units": "D_c = mu (geometric mean)"})
            else:
                rows.append({"component": label, "type": ctype, "param": pname,
                             "raw_theta": raw, "physical": raw, "units": ""})
        if ctype == "delta":
            rows.append({"component": label, "type": ctype, "param": "sigma_fixed",
                         "raw_theta": 0.005, "physical": 0.005 * span,
                         "units": "R (width, internal)"})
    return pd.DataFrame(rows)


def center_for_geometric_mean(MODEL, D_target=1.0):
    """The linear `center` value that places the lognormal geometric mean
    (mu) at diameter `D_target`. center=0.5 is the linear midpoint of the
    span, which is NOT D_g=1; this returns the center giving D_c=D_target.
    """
    R_min, R_max = MODEL["R_min"], MODEL["R_max"]
    span = R_max - R_min
    return float((0.5 * float(D_target) - R_min) / span)


def lognormal_diagnostics(theta, MODEL, component=None):
    """Analytic diagnostics for a 'lognormal' component (exact, independent
    of the grid/Tukey sampling path). `component` is an index into
    MODEL['components']; default picks the first lognormal.

    Returns a dict: alpha, cutoff, center, S, S_eff, D_c (=mu, geometric
    mean), D_min, D_max, a (=ln(S_eff)/(2 alpha)), parent_sD (untruncated
    coefficient of variation sqrt(exp(alpha^2)-1)). For the realized sample
    CV, std/mean of `sample_diameters(theta, MODEL, N)` (scale-invariant).
    """
    theta = np.asarray(theta, dtype=float)
    comps = MODEL["components"]
    idx = None
    if isinstance(component, int) and 0 <= component < len(comps) \
            and comps[component]["type"] == "lognormal":
        idx = component
    else:
        for i, c in enumerate(comps):
            if c["type"] == "lognormal":
                idx = i
                break
    if idx is None:
        raise ValueError("no 'lognormal' component found in MODEL")
    c = comps[idx]
    p = theta[c["slice_start"]:c["slice_end"]]
    alpha = float(p[1])
    cutoff = float(p[2])
    center = float(p[3])
    S = float(MODEL.get("size_ratio", 1.0))
    R_min, R_max = MODEL["R_min"], MODEL["R_max"]
    span = R_max - R_min
    D_c = 2.0 * (R_min + center * span)          # geometric mean (mu)
    S_eff = max(cutoff * S, 1.0)
    h = 0.5 * np.log(S_eff)
    D_min = float(D_c * np.exp(-h))
    D_max = float(D_c * np.exp(h))
    a = float(h / max(alpha, 1e-12))             # = ln(S_eff)/(2 alpha)
    parent_sD = float(np.sqrt(np.exp(alpha ** 2) - 1.0))
    return {"component": c["label"], "alpha": alpha, "cutoff": cutoff,
            "center": center, "S": S, "S_eff": float(S_eff),
            "D_c": float(D_c), "D_min": D_min, "D_max": D_max,
            "a": a, "parent_sD": parent_sD}
