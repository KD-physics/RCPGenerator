#!/usr/bin/env python3
"""Chord-gap decomposition of phi (Farr-Groot's language).

A random line alternates solid chord / void gap; exactly phi = solid_len / line_len
(mean-chord theorem). We cut AXIS-ALIGNED rays (exact periodic wrap; x/y/z also gives
the anisotropy check; statistically identical to isotropic lines for an isotropic
packing). Records: per-sphere chords L_s (each line-sphere intersection -> FG-marginal
check) and void gaps L_v (between merged solid runs -> the structural object).

phi-recovery (mean solid fraction over rays vs true phi) is the built-in correctness
gate. Usage: python chord_cut.py <packing.npz> [n_rays_per_axis]
"""
import sys, numpy as np
import voidfill_method as vf

def _runs_and_gaps(centers, halfs, Lax):
    """1D periodic interval union. Returns (void_gap_lengths, solid_total)."""
    ivs = []
    for c, hh in zip(centers % Lax, halfs):
        lo, hi = c - hh, c + hh
        if lo < 0:        ivs += [(0.0, hi), (lo + Lax, Lax)]
        elif hi > Lax:    ivs += [(lo, Lax), (0.0, hi - Lax)]
        else:             ivs.append((lo, hi))
    if not ivs:
        return np.array([Lax]), 0.0
    ivs.sort()
    merged = [list(ivs[0])]
    for lo, hi in ivs[1:]:
        if lo <= merged[-1][1] + 1e-12:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    solid = sum(h - l for l, h in merged)
    gaps = [merged[k + 1][0] - merged[k][1] for k in range(len(merged) - 1)]
    gaps.append((merged[0][0] + Lax) - merged[-1][1])          # wrap gap
    return np.array([g for g in gaps if g > 1e-9]), solid

def cut_axis(pos, dia, box, axis, n_rays, rng):
    r = dia / 2.0; t = [a for a in range(len(box)) if a != axis]
    Lax = box[axis]; per_sphere, voids = [], []; solids = []
    qy = rng.random(n_rays) * box[t[0]]; qz = rng.random(n_rays) * box[t[1]]
    for j in range(n_rays):
        dy = pos[:, t[0]] - qy[j]; dy -= box[t[0]] * np.round(dy / box[t[0]])
        dz = pos[:, t[1]] - qz[j]; dz -= box[t[1]] * np.round(dz / box[t[1]])
        d2 = dy * dy + dz * dz; hit = d2 < r * r
        if not hit.any():
            voids.append(Lax); solids.append(0.0); continue
        half = np.sqrt(r[hit] ** 2 - d2[hit])
        per_sphere.append(2.0 * half)
        g, sol = _runs_and_gaps(pos[hit, axis], half, Lax)
        voids.append(g); solids.append(sol)
    return (np.concatenate([np.atleast_1d(x) for x in per_sphere]) if per_sphere else np.array([]),
            np.concatenate([np.atleast_1d(x) for x in voids]),
            np.array(solids) / Lax)

def main():
    npz = sys.argv[1]; nray = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    pos, dia, box = vf._load(npz); d = len(box)
    phi_true = float(vf.sphere_vol(dia).sum() / np.prod(box))
    rng = np.random.default_rng(7)
    print(f"# {npz}  N={len(dia)}  S={dia.max()/dia.min():.0f}  phi_true={phi_true:.4f}")
    allLs, allLv = [], []
    for ax in range(d):
        Ls, Lv, sf = cut_axis(pos, dia, box, ax, nray, rng)
        allLs.append(Ls); allLv.append(Lv)
        print(f"  axis {ax}: phi_recov={sf.mean():.4f}  (resid {abs(sf.mean()-phi_true):.4f})  "
              f"<Ls>={Ls.mean():.4f} <Lv>={Lv.mean():.4f}  n_chord={len(Ls)} n_gap={len(Lv)}")
    Ls = np.concatenate(allLs); Lv = np.concatenate(allLv)
    print(f"  ALL: <Ls>={Ls.mean():.4f}  <Lv>={Lv.mean():.4f}  "
          f"phi_from_means~ n/a (per-ray mean is the phi)")

if __name__ == "__main__":
    main()
