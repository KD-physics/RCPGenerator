"""Internal helpers shared by the search package."""
from __future__ import annotations

import math

import numpy as np


def hypersphere_volume_from_diameter(diameters, ndim):
    """Total d-dimensional measure for hyperspheres with given diameters."""
    diameters = np.asarray(diameters, dtype=float)
    radii = diameters / 2.0
    if ndim == 2:
        return float(np.sum(np.pi * radii ** 2))
    if ndim == 3:
        return float(np.sum((4.0 / 3.0) * np.pi * radii ** 3))
    unit_ball = np.pi ** (ndim / 2.0) / math.gamma(ndim / 2.0 + 1.0)
    return float(np.sum(unit_ball * radii ** ndim))


def get_final_phi(packing):
    """Prefer the packed phi_final; fall back to the requested phi."""
    if hasattr(packing, "phi_final") and packing.phi_final is not None:
        return float(packing.phi_final)
    if hasattr(packing, "phi") and packing.phi is not None:
        return float(packing.phi)
    raise AttributeError("packing has no phi_final or phi")


def safe_compute_phi(packing):
    """Recompute phi from (positions, diameters, box). Returns 0.0 on any
    sanity failure (non-finite positions/diameters, particles outside box,
    box <= 0). Independent of rcpgenerator's internal phi_final."""
    pos = None
    for attr in ("positions", "points", "x", "coords"):
        if hasattr(packing, attr):
            pos = np.asarray(getattr(packing, attr), dtype=float)
            break
    if pos is None or pos.size == 0:
        return 0.0
    D = np.asarray(packing.diameters, dtype=float)
    if D.size != pos.shape[0]:
        return 0.0
    N, ndim = pos.shape

    box = None
    for attr in ("box", "L", "box_size"):
        if hasattr(packing, attr):
            b = np.asarray(getattr(packing, attr), dtype=float)
            box = (b[1::2] - b[0::2]) if b.size == 2 * ndim else b
            break
    if box is None:
        return 0.0
    box = np.atleast_1d(box).astype(float)
    if box.size == 1:
        box = np.full(ndim, float(box.item()))
    if box.size != ndim:
        return 0.0

    if not np.all(np.isfinite(pos)):
        return 0.0
    if not np.all(np.isfinite(D)) or np.any(D <= 0):
        return 0.0
    tol = 1e-6 * np.max(box)
    if np.any(pos < -tol) or np.any(pos > box + tol):
        return 0.0
    if np.any(box <= 0):
        return 0.0

    V_particles = hypersphere_volume_from_diameter(D, ndim)
    V_box = float(np.prod(box))
    return V_particles / V_box if V_box > 0 else 0.0


def auto_neighbor_max(N, Ndim, size_ratio=1.0):
    """Default per-particle neighbor allocation when CONFIG['neighbor_max']==0.

    Scales as the C++ default does: with Ndim and the diameter ratio.
    Capped at 7500 to keep memory bounded on large N.
    """
    base = 200.0
    scaling = Ndim / 3.0
    omega = 1.0 + Ndim / 6.0
    cap = scaling * (50.0 * (size_ratio ** omega) + base)
    cap = max(50.0, min(cap, 7500.0))
    return int(round(cap))


def summarize_packing(label, packing):
    """Compact text summary for ad-hoc inspection."""
    D = np.asarray(packing.diameters, dtype=float)
    print(f"\ncase: {label}")
    print("  N:", packing.N)
    print("  Ndim:", packing.Ndim)
    print("  box:", packing.box)
    print("  initial phi:", packing.phi)
    print("  final phi:", packing.phi_final)
    print("  steps:", packing.steps)
    print("  force_magnitude:", packing.force_magnitude)
    print("  max_min_dist:", packing.max_min_dist)
    if len(D):
        print("  D_min:", float(D.min()))
        print("  D_max:", float(D.max()))
        print("  D_ratio:", float(D.max() / D.min()))
