"""Per-particle local packing fraction via periodic radical Voronoi.

Requires the `pyvoro-mmalahe` package (drop-in replacement for `pyvoro`
that handles radii correctly). Importing this module without the package
installed will raise ImportError.
"""
from __future__ import annotations

import numpy as np


def voronoi_phi_local(packing):
    """Per-particle local packing fraction with periodic BC.

    Returns dict:
      local_phi               : per-particle local phi
      cell_volumes            : per-particle Voronoi volume (3D) or area (2D)
      box_volume              : nominal box volume / area
      sum_volumes             : sum of cell volumes
      sum_match               : True iff sum matches box_volume within 1e-9 rel
      sum_match_relative_error: actual relative mismatch
    """
    import pyvoro

    positions = np.asarray(packing.positions, dtype=float)
    radii = np.asarray(packing.diameters, dtype=float) / 2.0
    Ndim = int(packing.Ndim)
    box = list(packing.box)
    limits = [[0.0, float(L)] for L in box]
    dispersion = float(2.0 * radii.max())

    if Ndim == 2:
        cells = pyvoro.compute_2d_voronoi(
            positions.tolist(), limits, dispersion,
            radii=radii.tolist(), periodic=[True, True],
        )
        particle_vol = np.pi * radii ** 2
        box_vol = box[0] * box[1]
    elif Ndim == 3:
        cells = pyvoro.compute_voronoi(
            positions.tolist(), limits, dispersion,
            radii=radii.tolist(), periodic=[True, True, True],
        )
        particle_vol = (4.0 / 3.0) * np.pi * radii ** 3
        box_vol = box[0] * box[1] * box[2]
    else:
        raise ValueError(f"Voronoi only for Ndim=2 or 3, got Ndim={Ndim}")

    cell_vols = np.array([c["volume"] for c in cells], dtype=float)
    local_phi = particle_vol / cell_vols
    sum_vols = float(cell_vols.sum())
    rel_err = abs(sum_vols - box_vol) / box_vol
    return {
        "local_phi": local_phi,
        "cell_volumes": cell_vols,
        "box_volume": box_vol,
        "sum_volumes": sum_vols,
        "sum_match": rel_err < 1e-9,
        "sum_match_relative_error": rel_err,
    }


def plot_local_phi_histogram(packing, *, bins=50, title=None, figsize=(8, 4)):
    """Histogram of per-particle local phi."""
    import matplotlib.pyplot as plt

    v = voronoi_phi_local(packing)
    lp = v["local_phi"]
    plt.figure(figsize=figsize)
    plt.hist(lp, bins=bins, density=True, alpha=0.7)
    plt.axvline(lp.mean(), color="k", linestyle="--",
                label=f"mean = {lp.mean():.4f}")
    plt.xlabel("local phi")
    plt.ylabel("density")
    plt.title(title or "per-particle local packing fraction")
    plt.legend()
    plt.show()
    return v


def plot_local_phi_vs_diameter(packing, *, title=None, figsize=(8, 4),
                               xlog=False, alpha=0.4):
    """Scatter of per-particle local phi against particle diameter."""
    import matplotlib.pyplot as plt

    v = voronoi_phi_local(packing)
    D = np.asarray(packing.diameters, dtype=float)
    plt.figure(figsize=figsize)
    plt.scatter(D, v["local_phi"], alpha=alpha, s=8)
    plt.xlabel("particle diameter D")
    plt.ylabel("local phi")
    if xlog:
        plt.xscale("log")
    plt.title(title or "local phi vs. particle diameter")
    plt.show()
    return v
