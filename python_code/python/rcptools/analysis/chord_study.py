#!/usr/bin/env python3
"""Chord-gap decomposition across the S-ladder: the deficit = excess void-gap, and
the gap-distribution SHAPE = mechanism (fat tail -> voids ; uniform shift -> uniform
inefficiency). Real (line-cut) vs FG (greedy). Cutoff-free.

Usage: python chord_study.py [n_rays]
"""
import glob, pathlib, heapq, numpy as np
import voidfill_method as vf, chord_cut as cc

def greedy_gaps(chords, f=0.7654):
    """FG greedy 1D pack -> (phi, final void-gap lengths). Mirrors vf._greedy_phi
    but also returns the leftover gap segments (FG's predicted void distribution)."""
    r = sorted(chords, reverse=True)
    used = r[0]; h = [-(f * r[0])]
    for L in r[1:]:
        g = -heapq.heappop(h); c = (1 + f) * L
        a, b = (f * L, max(g - c, f * L)) if g >= c else (f * L, f * L)
        heapq.heappush(h, -a); heapq.heappush(h, -b); used += L
    gaps = np.array([-x for x in h])
    return used / (used + gaps.sum()), gaps

def fg_chords_and_gaps(dia, n=20000, seed=12345, f=0.7654):
    D = np.asarray(dia, float); w = D ** 2; w = w / w.sum()
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(D), size=min(n, len(D) * 50), p=w)
    L = D[idx] * np.sqrt(rng.random(len(idx)))
    phi, gaps = greedy_gaps(L, f)
    return L, phi, gaps

def shape(gaps):
    """scale-free gap-shape stats: CV and fraction of total void in the top-10% gaps."""
    g = np.sort(gaps); n = len(g)
    cv = g.std() / g.mean()
    top = g[int(0.9 * n):].sum() / g.sum()
    return cv, top

def analyze(npz, nray):
    pos, dia, box = vf._load(npz); d = len(box)
    phi_true = float(vf.sphere_vol(dia).sum() / np.prod(box))
    rng = np.random.default_rng(7)
    Ls, Lv, sf = [], [], []
    for ax in range(d):
        a, b, s = cc.cut_axis(pos, dia, box, ax, nray, rng)
        Ls.append(a); Lv.append(b); sf.append(s.mean())
    Ls = np.concatenate(Ls); Lv = np.concatenate(Lv)
    fgL, phi_fg, fgGap = fg_chords_and_gaps(dia)
    # normalize gaps by mean per-sphere chord (geometric scale) -> dimensionless
    lv_n = Lv / Ls.mean(); fg_n = fgGap / fgL.mean()
    cv_r, top_r = shape(lv_n); cv_f, top_f = shape(fg_n)
    phi_gm = greedy_gaps(Ls)[0]                       # greedy on MEASURED rods
    return dict(S=dia.max() / dia.min(), phi_true=phi_true, phi_recov=np.mean(sf),
                phi_fg=phi_fg, phi_gm=phi_gm, deficit=phi_fg - phi_true,
                lvn_real=lv_n.mean(), lvn_fg=fg_n.mean(), cv_r=cv_r, cv_f=cv_f,
                top_r=top_r, top_f=top_f, Ls=Ls, fgL=fgL, Lv=lv_n, fg_n=fg_n)

def make_figs(res, FIG, tag="pl"):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rs = [r for _, r in res]; Ss = sorted({round(r["S"]) for r in rs})
    cmap = plt.cm.viridis; bins = np.linspace(0, 1.5, 60)
    fig, ax = plt.subplots(2, 2, figsize=(13, 10))
    # (A) gap distribution real (solid) vs FG (dashed), per S
    for k, S in enumerate(Ss):
        col = cmap(k / max(1, len(Ss) - 1))
        lv = np.concatenate([r["Lv"] for r in rs if round(r["S"]) == S])
        fg = np.concatenate([r["fg_n"] for r in rs if round(r["S"]) == S])
        ax[0, 0].hist(lv, bins=bins, density=True, histtype="step", color=col, lw=2, label=f"real S={S}")
        ax[0, 0].hist(fg, bins=bins, density=True, histtype="step", color=col, lw=1, ls="--")
    ax[0, 0].set_xlabel("void gap / <L_s>"); ax[0, 0].set_yscale("log")
    ax[0, 0].set_title("gap distribution: real (solid) vs FG-greedy (dashed)"); ax[0, 0].legend(fontsize=7)
    # (B) deficit and excess mean gap vs S
    Sx = [r["S"] for r in rs]
    ax[0, 1].scatter(Sx, [r["deficit"] for r in rs], label="deficit phiFG-phiT")
    ax[0, 1].scatter(Sx, [r["lvn_real"] - r["lvn_fg"] for r in rs], marker="^", label="excess mean gap (real-FG)")
    ax[0, 1].set_xlabel("S"); ax[0, 1].set_title("deficit & excess gap vs S"); ax[0, 1].legend(fontsize=8)
    # (C) gap-shape vs S
    ax[1, 0].scatter(Sx, [r["cv_r"] for r in rs], label="CV real")
    ax[1, 0].scatter(Sx, [r["cv_f"] for r in rs], marker="x", label="CV FG")
    ax[1, 0].scatter(Sx, [r["top_r"] for r in rs], marker="^", label="top10% void real")
    ax[1, 0].scatter(Sx, [r["top_f"] for r in rs], marker="v", label="top10% void FG")
    ax[1, 0].set_xlabel("S"); ax[1, 0].set_title("gap-shape: fat tail (real) vs uniform (FG)"); ax[1, 0].legend(fontsize=8)
    # (D) solid-chord sanity: measured P(Ls) vs FG fgL, one packing
    r0 = rs[-1]
    b2 = np.linspace(0, np.percentile(r0["Ls"], 99), 50)
    ax[1, 1].hist(r0["Ls"], bins=b2, density=True, histtype="step", lw=2, label="measured L_s")
    ax[1, 1].hist(r0["fgL"], bins=b2, density=True, histtype="step", lw=2, ls="--", label="FG D^2 sqrt(U)")
    ax[1, 1].set_xlabel("solid chord L_s"); ax[1, 1].set_title("solid-chord marginal: measured vs FG (geometry)")
    ax[1, 1].legend(fontsize=8)
    fig.suptitle(f"Chord-gap decomposition ({tag}): deficit = excess gap; mechanism = fat gap tail", fontsize=13)
    fig.tight_layout(); out = FIG / f"chord_decomp_{tag}.png"
    fig.savefig(out, dpi=130); plt.close(); print(f"\nfigure -> {out}")

def main():
    import sys
    nray = int(sys.argv[1]) if len(sys.argv) > 1 else 250
    DATA = pathlib.Path(__file__).resolve().parent.parent / "Data" / "voidfill"
    print(f"# chord-gap decomposition (n_rays/axis={nray}); gaps normalized by <Ls>")
    print(f"# {'pack':13s} {'S':>4} {'phiT':>6} {'phiRec':>6} {'phiFG':>6} {'defic':>6} "
          f"{'<Lv>_r':>6} {'<Lv>_FG':>7} {'CV_r':>5} {'CV_FG':>5} {'top_r':>5} {'top_FG':>6}")
    res = []
    for fam in ("pl_p355", "pl_p363"):
        for f in sorted(glob.glob(str(DATA / f"{fam}_S*.npz")),
                        key=lambda x: int(x.split('_S')[-1].split('.')[0])):
            r = analyze(f, nray); res.append((pathlib.Path(f).stem, r))
            print(f"  {pathlib.Path(f).stem:13s} {r['S']:4.0f} {r['phi_true']:6.3f} "
                  f"{r['phi_recov']:6.3f} {r['phi_fg']:6.3f} {r['deficit']:6.3f} "
                  f"{r['lvn_real']:6.3f} {r['lvn_fg']:7.3f} {r['cv_r']:5.2f} {r['cv_f']:5.2f} "
                  f"{r['top_r']:5.2f} {r['top_f']:6.2f}")
    print("\n# greedy(measured rods) should ~= phiFG (>phiT) -> rods fine, gaps are the issue:")
    for name, r in res:
        print(f"  {name:13s} phiT={r['phi_true']:.3f} greedy(meas L_s)={r['phi_gm']:.3f} "
              f"phiFG={r['phi_fg']:.3f}")
    make_figs(res, pathlib.Path(__file__).resolve().parent.parent / "Figures", "pl")
    return res

if __name__ == "__main__":
    main()
