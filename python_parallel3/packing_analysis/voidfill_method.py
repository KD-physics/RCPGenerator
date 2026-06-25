#!/usr/bin/env python3
"""Void-fill method (Gate 1): quantify how completely the voids of a coarse
scaffold are filled by finer particles, with a scale-normalized per-pore metric.

Pipeline (see Claude/VOIDFILL_METHOD.md):
  - scaffold: coarse = particles with D > D_c.
  - voids:    radical (Laguerre) Voronoi of the COARSE set only (periodic); each
              coarse particle -> one cell that tiles the box. void W_c = V_c - v_big.
  - fill:     eta_c = (fine solid volume assigned to cell c) / W_c, fines assigned
              by radical (power-distance) membership to coarse cells.
  - normalize:rho_c = eta_c / phi_fill*, phi_fill* = FG RCP of the fine sub-set.
  - summaries:<rho> (vol-weighted), f_under (void-vol fraction with rho<rho_crit),
              rho_min, and integrated underfill Delta_phi_recon.

All quantities dimensionless and invariant under uniform rescaling of (D, box).
Importable; CLI: `python voidfill_method.py <packing.npz> [Dc ...]`.
"""
import sys, pathlib, numpy as np
import pyvoro
from scipy.spatial import cKDTree

# ----------------------------------------------------------------------- FG ref
def _greedy_phi(chords, f=0.7654):
    """1D greedy rod packing -> calibrated 3D phi (f fixes mono -> 0.6435)."""
    import heapq
    r = sorted(chords, reverse=True)
    used = r[0]; gap = f * r[0]; h = [-(f * r[0])]
    for L in r[1:]:
        g = -heapq.heappop(h); gap -= g; c = (1 + f) * L
        a, b = (f * L, max(g - c, f * L)) if g >= c else (f * L, f * L)
        heapq.heappush(h, -a); heapq.heappush(h, -b); gap += a + b; used += L
    return used / (used + gap)

def fg_rcp(diameters, n=20000, runs=6, seed=12345, f=0.7654):
    """Farr-Groot RCP of an explicit diameter set: D^2-weighted chords L=D*sqrt(U),
    greedy 1D, calibrated. Scale-invariant (depends only on diameter ratios)."""
    D = np.asarray(diameters, float)
    if len(D) == 0:
        return np.nan
    if np.allclose(D, D[0]):
        return 0.6435  # monodisperse anchor
    w = D**2; w = w / w.sum()
    rng = np.random.default_rng(seed); vals = []
    for _ in range(runs):
        idx = rng.choice(len(D), size=min(n, 50_000_000), p=w)
        L = D[idx] * np.sqrt(rng.random(len(idx)))
        vals.append(_greedy_phi(L, f))
    return float(np.mean(vals))

# ------------------------------------------------------------------- geometry
def sphere_vol(d):
    return (np.pi / 6.0) * np.asarray(d, float)**3

def radical_voronoi_volumes(pos, radii, box):
    """Per-particle radical-Voronoi cell volume (periodic). Borrowed call pattern
    from the rcp_search notebook `voronoi_phi_local`."""
    pos = np.asarray(pos, float); radii = np.asarray(radii, float)
    L = [float(b) for b in box]
    limits = [[0.0, L[i]] for i in range(3)]
    disp = float(2.0 * radii.max())
    cells = pyvoro.compute_voronoi((pos % np.array(L)).tolist(), limits, disp,
                                   radii=radii.tolist(), periodic=[True, True, True])
    return np.array([c["volume"] for c in cells], float)

def assign_fines_to_coarse(fine_pos, coarse_pos, coarse_r, box, kq=12):
    """Radical (power-distance) membership: each fine -> coarse j minimizing
    |x-x_j|^2 - r_j^2. Periodic. Query k nearest coarse centers, pick min power
    distance among them (exact when the true owner is within the k nearest)."""
    L = np.asarray(box, float)
    tree = cKDTree(coarse_pos % L, boxsize=L)
    k = min(kq, len(coarse_pos))
    _, nn = tree.query(fine_pos % L, k=k)
    nn = np.atleast_2d(nn.T).T if k > 1 else nn.reshape(-1, 1)
    best = np.full(len(fine_pos), -1, int); bestpd = np.full(len(fine_pos), np.inf)
    for col in range(k):
        j = nn[:, col]
        dr = fine_pos - coarse_pos[j]; dr -= L * np.round(dr / L)
        pd = np.einsum("ij,ij->i", dr, dr) - coarse_r[j]**2
        upd = pd < bestpd
        best[upd] = j[upd]; bestpd[upd] = pd[upd]
    return best  # index into the coarse arrays

# ------------------------------------------------------------------- the metric
def void_fill(pos, dia, box, Dc, phi_fill=None, fg_kw=None):
    """Core metric at a single cutoff Dc. Returns a dict of per-cell arrays and
    scalar summaries."""
    pos = np.asarray(pos, float); dia = np.asarray(dia, float); box = np.asarray(box, float)
    Vbox = float(np.prod(box))
    coarse = dia > Dc; fine = ~coarse
    nc = int(coarse.sum())
    if nc < 4 or fine.sum() == 0:
        return None
    cpos, cdia = pos[coarse], dia[coarse]
    cvol = radical_voronoi_volumes(cpos, cdia / 2.0, box)          # cells tile box
    vbig = sphere_vol(cdia)
    W = cvol - vbig                                                # per-pore void vol
    # assign fines -> coarse cell, accumulate fine solid per cell
    owner = assign_fines_to_coarse(pos[fine], cpos, cdia / 2.0, box)
    fvol = sphere_vol(dia[fine])
    S = np.zeros(nc); np.add.at(S, owner, fvol)
    eta = np.divide(S, W, out=np.zeros_like(S), where=W > 0)
    if phi_fill is None:
        phi_fill = fg_rcp(dia[fine], **(fg_kw or {}))
    rho = eta / phi_fill
    # summaries
    rho_mean = float(S.sum() / (phi_fill * W.sum()))              # vol-weighted <rho>
    dphi_recon = float((phi_fill * W.sum() - S.sum()) / Vbox)     # integrated underfill
    return dict(Dc=float(Dc), n_coarse=nc, n_fine=int(fine.sum()),
                cell_vol=cvol, void=W, fine_solid=S, eta=eta, rho=rho,
                phi_fill=float(phi_fill),
                phi_scaffold=float(vbig.sum() / Vbox),            # coarse solid / box
                void_frac=float(W.sum() / Vbox),
                rho_mean=rho_mean, rho_min=float(rho.min()),
                dphi_recon=dphi_recon)

def f_under(res, rho_crit):
    """Void-volume fraction in underfilled pores (rho < rho_crit)."""
    m = res["rho"] < rho_crit
    return float(res["void"][m].sum() / res["void"].sum())

def sweep_cutoff(pos, dia, box, Dc_list, **kw):
    return [r for Dc in Dc_list for r in [void_fill(pos, dia, box, Dc, **kw)] if r]

# ------------------------------------------------------------------- ablation
def ablate_uniform(dia, fine_mask, q, seed=0):
    """Remove a fraction q of fine *volume* at random. Returns (keep_mask,
    removed_volume)."""
    rng = np.random.default_rng(seed)
    fidx = np.where(fine_mask)[0]
    fvol = sphere_vol(dia[fidx])
    order = rng.permutation(len(fidx))
    target = q * fvol.sum(); cum = np.cumsum(fvol[order])
    ncut = int(np.searchsorted(cum, target)) + 1
    drop = fidx[order[:ncut]]
    keep = np.ones(len(dia), bool); keep[drop] = False
    return keep, float(fvol[order[:ncut]].sum())

def ablate_graded(pos, dia, box, Dc, cell_indices, fracs, seed=0):
    """Remove a known fraction `frac` of the fine *volume* in each listed cell.
    Returns (keep_mask, total_removed_volume, {cell: removed_volume})."""
    coarse = dia > Dc; fine = ~coarse
    cpos, cdia = pos[coarse], dia[coarse]
    owner = assign_fines_to_coarse(pos[fine], cpos, cdia / 2.0, box)
    fidx = np.where(fine)[0]; rng = np.random.default_rng(seed)
    keep = np.ones(len(dia), bool); removed = 0.0; imposed = {}
    for jc, g in zip(cell_indices, fracs):
        members = fidx[owner == jc]
        if len(members) == 0:
            imposed[int(jc)] = 0.0; continue
        fv = sphere_vol(dia[members]); order = rng.permutation(len(members))
        tgt = g * fv.sum(); cum = np.cumsum(fv[order])
        ncut = int(np.searchsorted(cum, tgt)) + 1
        drop = members[order[:ncut]]; keep[drop] = False
        rv = float(fv[order[:ncut]].sum()); removed += rv; imposed[int(jc)] = rv
    return keep, float(removed), imposed

def ablate_pore(pos, dia, box, Dc, cell_index, seed=0):
    """Empty the fines of ONE coarse cell. Returns (keep_mask, removed_volume,
    cell_index)."""
    coarse = dia > Dc; fine = ~coarse
    cpos, cdia = pos[coarse], dia[coarse]
    owner = assign_fines_to_coarse(pos[fine], cpos, cdia / 2.0, box)
    fidx = np.where(fine)[0]
    drop = fidx[owner == cell_index]
    keep = np.ones(len(dia), bool); keep[drop] = False
    return keep, float(sphere_vol(dia[drop]).sum()), int(cell_index)

# ------------------------------------------------------------------- CLI
def _load(npz):
    d = np.load(npz, allow_pickle=True)
    return d["positions"], d["diameters"], d["box"]

def main():
    npz = sys.argv[1]
    pos, dia, box = _load(npz)
    Dmin, Dmax = dia.min(), dia.max()
    if len(sys.argv) > 2:
        Dcs = [float(x) for x in sys.argv[2:]]
    else:
        Dcs = list(np.geomspace(Dmin * 2, Dmax * 0.6, 8))
    print(f"# {pathlib.Path(npz).name}  N={len(dia)} Dmin={Dmin:.4f} Dmax={Dmax:.4f} "
          f"phi={(sphere_vol(dia).sum()/np.prod(box)):.4f}")
    print(f"# {'Dc':>8} {'Dc/Dmax':>7} {'n_co':>6} {'phi_sc':>7} {'voidfr':>7} "
          f"{'phi_fill':>8} {'<rho>':>6} {'rho_min':>7} {'dphi_recon':>10}")
    for Dc in Dcs:
        r = void_fill(pos, dia, box, Dc)
        if r is None:
            continue
        print(f"  {r['Dc']:8.4f} {r['Dc']/Dmax:7.3f} {r['n_coarse']:6d} "
              f"{r['phi_scaffold']:7.4f} {r['void_frac']:7.4f} {r['phi_fill']:8.4f} "
              f"{r['rho_mean']:6.3f} {r['rho_min']:7.3f} {r['dphi_recon']:10.4f}")

if __name__ == "__main__":
    main()
