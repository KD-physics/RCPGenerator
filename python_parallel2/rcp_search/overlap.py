"""Pair-overlap diagnostics for a packed configuration."""
from __future__ import annotations

import numpy as np


def overlap_report(packing, frac_def="min", verbose=False):
    """Distribution of fractional overlaps over all overlapping pairs.

    Returns a dict with:
      - n_overlapping       : number of overlapping pairs
      - n_candidate_pairs   : number of close pairs queried from the tree
      - fractions_min       : overlap / min(r_i, r_j)  (DW09 convention)
      - fractions_sum       : overlap / (r_i + r_j)    (contraction calc)
      - max_overlap         : worst pair's `fractions_min`
      - mean_overlap        : mean of `fractions_min` over overlapping pairs
      - X_shrink_to_clear   : uniform contraction needed to clear all overlaps
      - ndim
      - pairs               : (m, 2) index array of overlapping pairs

    `verbose=True` prints a per-quantile summary.
    """
    from scipy.spatial import cKDTree

    pos = None
    for attr in ("positions", "points", "x", "coords"):
        if hasattr(packing, attr):
            pos = np.asarray(getattr(packing, attr), dtype=float)
            break
    if pos is None:
        raise AttributeError("packing has no positions/points/x/coords attribute")
    D = np.asarray(packing.diameters, dtype=float)
    R = D / 2.0
    N, ndim = pos.shape

    boxsize = None
    for attr in ("box", "L", "box_size"):
        if hasattr(packing, attr):
            b = np.asarray(getattr(packing, attr), dtype=float)
            boxsize = (b[1::2] - b[0::2]) if b.size == 2 * ndim else b
            break

    if N == 0 or R.size == 0:
        return {
            "n_overlapping": 0, "n_candidate_pairs": 0,
            "fractions_min": np.array([]), "fractions_sum": np.array([]),
            "max_overlap": 0.0, "mean_overlap": 0.0,
            "X_shrink_to_clear": 0.0, "ndim": ndim, "pairs": np.zeros((0, 2), int),
        }

    tree = cKDTree(pos, boxsize=boxsize) if boxsize is not None else cKDTree(pos)
    pairs = tree.query_pairs(r=2.0 * R.max(), output_type="ndarray")
    if pairs.shape[0] == 0:
        if verbose:
            print("No candidate pairs found.")
        return {
            "n_overlapping": 0, "n_candidate_pairs": 0,
            "fractions_min": np.array([]), "fractions_sum": np.array([]),
            "max_overlap": 0.0, "mean_overlap": 0.0,
            "X_shrink_to_clear": 0.0, "ndim": ndim, "pairs": np.zeros((0, 2), int),
        }

    i_idx, j_idx = pairs[:, 0], pairs[:, 1]
    delta = pos[i_idx] - pos[j_idx]
    if boxsize is not None:
        delta -= np.round(delta / boxsize) * boxsize
    dist = np.linalg.norm(delta, axis=1)
    rsum = R[i_idx] + R[j_idx]
    overlap = rsum - dist
    mask = overlap > 0
    overlap = overlap[mask]
    Ri, Rj = R[i_idx][mask], R[j_idx][mask]
    rsum = rsum[mask]
    rmin = np.minimum(Ri, Rj)

    f_min = overlap / rmin
    f_sum = overlap / rsum
    X_clear = float(f_sum.max()) if f_sum.size else 0.0

    out = {
        "n_overlapping": int(mask.sum()),
        "n_candidate_pairs": int(pairs.shape[0]),
        "fractions_min": f_min,
        "fractions_sum": f_sum,
        "max_overlap": float(f_min.max()) if f_min.size else 0.0,
        "mean_overlap": float(f_min.mean()) if f_min.size else 0.0,
        "X_shrink_to_clear": X_clear,
        "ndim": ndim,
        "pairs": pairs[mask],
    }

    if verbose:
        rep = f_min if frac_def == "min" else f_sum
        label = "overlap / min(r_i,r_j)" if frac_def == "min" else "overlap / (r_i+r_j)"
        print(f"N = {N}, candidate close-pairs = {pairs.shape[0]}, "
              f"overlapping = {out['n_overlapping']} "
              f"({100 * out['n_overlapping'] / (N * (N - 1) / 2):.4f}% of all pairs)")
        print(f"fractional overlap ({label}):")
        if rep.size:
            q = np.quantile(rep, [0.5, 0.9, 0.99, 1.0])
            print(f"  median = {q[0]:.2e}")
            print(f"  p90    = {q[1]:.2e}")
            print(f"  p99    = {q[2]:.2e}")
            print(f"  max    = {q[3]:.2e}")
            print(f"shrink-to-clear: X = {X_clear:.2e}; "
                  f"phi_new ≈ phi * (1 - X)**{ndim} = phi * {(1 - X_clear) ** ndim:.6f}")
        else:
            print("  no overlapping pairs.")
    return out


def phi_corrected(phi, max_overlap, Ndim):
    """Quality-corrected phi: phi - (Ndim/2) * max_overlap.

    max_overlap is the worst pair's fractional overlap (already
    dimensionless). Subtracts the leading-order phi loss that a uniform
    contraction would incur to clear the worst contact.
    """
    return float(phi) - 0.5 * float(Ndim) * float(max_overlap)
