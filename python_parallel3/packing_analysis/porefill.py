#!/usr/bin/env python3
"""Pore-fill study: per-CAVITY fill distribution vs size ratio S (and the
accessibility split), using the verified pore-network tool (porenet, frac=0.35).

Per scaffold cavity (merged Delaunay-cell pore): void W, contained fine solid Sf,
fill eta=Sf/W and rho=eta/phi_fill (phi_fill = FG-of-fines). Plus accessibility:
the cavity's widest throat (max Delaunay-pair gap) vs its own size -> constriction.

The point (vs Gate-2's mean): look at the LOW-fill TAIL and whether it grows with S.
Usage: python porefill.py   (loops the existing p=-3.55/-3.63 S-ladders)
"""
import glob, pathlib, numpy as np
import porenet as pn, voidfill_method as vf

FRAC = 0.35; DCFRAC = 0.35           # gap-merge frac ; scaffold cutoff Dc/Dmax
ETA_STARVE = 0.5                     # cavity "starved" if eta < this * phi_fill (i.e. rho<0.5)

def per_cavity(npz):
    pos, dia, box = vf._load(npz); d = len(box); Vbox = float(np.prod(box))
    Dc = DCFRAC * dia.max()
    info = pn.assign_fines(pn.pore_network(pos, dia, box, Dc))
    th = pn.throats(info)
    info = pn.merge_pores_gap(info, th, frac=FRAC)
    bodies = info["merge"]["bodies"]
    void = info["void"]; fine = info["fine_solid"]
    phi_fill = vf.fg_rcp(dia[dia <= Dc])
    # per-body aggregates
    Wb, Sb, sizeb = [], [], []
    body_of = {}                                      # simplex -> body id
    for bid, ix in bodies.items():
        for i in ix: body_of[i] = bid
        Wb.append(void[ix].sum()); Sb.append(fine[ix].sum())
        sizeb.append(max(void[ix].sum(), 0.0)**(1.0 / d))   # cavity linear size
    Wb = np.array(Wb); Sb = np.array(Sb); sizeb = np.array(sizeb)
    eta = np.divide(Sb, Wb, out=np.zeros_like(Sb), where=Wb > 0)
    rho = eta / phi_fill
    # accessibility: widest throat (gate) per body, from inter-body throats
    bid_list = list(bodies.keys()); idx = {b: k for k, b in enumerate(bid_list)}
    gate = np.zeros(len(bid_list))
    for a, b, area, gap, ref in th:
        ba, bb = body_of.get(int(a)), body_of.get(int(b))
        if ba is None or bb is None or ba == bb:
            continue
        for B in (ba, bb):
            gate[idx[B]] = max(gate[idx[B]], gap)     # widest open gate into the cavity
    constr = np.divide(gate, sizeb, out=np.zeros_like(gate), where=sizeb > 0)
    phi_sim = float(vf.sphere_vol(dia).sum() / Vbox); phi_fg = vf.fg_rcp(dia)
    return dict(S=float(dia.max() / dia.min()), n_cav=len(rho), rho=rho, eta=eta,
                W=Wb, size=sizeb, gate=gate, constr=constr, phi_fill=phi_fill,
                deficit=phi_fg - phi_sim, Vbox=Vbox,
                starved_underfill=float(((phi_fill * Wb - Sb)[rho < ETA_STARVE]).sum() / Vbox))

def make_figs(res, FIG):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    Ss = sorted({r["S"] for r in res})
    cmap = plt.cm.viridis
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    # (A) P(rho) per S, pooled over both exponents
    for k, S in enumerate(Ss):
        rho = np.concatenate([r["rho"] for r in res if r["S"] == S])
        ax[0].hist(rho, bins=np.linspace(0, 1.3, 40), histtype="step", density=True,
                   color=cmap(k / max(1, len(Ss) - 1)), lw=2, label=f"S={S:.0f} (n={len(rho)})")
    ax[0].axvline(0.5, color="r", ls=":", lw=1); ax[0].set_xlabel("per-cavity fill rho = eta / phi_fill")
    ax[0].set_ylabel("density"); ax[0].set_title("P(rho) vs S  (tail = starved cavities)")
    ax[0].legend(fontsize=8)
    # (B) accessibility: rho vs constriction ratio (pooled)
    rho = np.concatenate([r["rho"] for r in res]); con = np.concatenate([r["constr"] for r in res])
    m = np.isfinite(rho) & np.isfinite(con) & (con > 0)
    cc = np.corrcoef(rho[m], con[m])[0, 1]
    ax[1].scatter(con[m], rho[m], s=5, alpha=0.25)
    ax[1].axhline(0.5, color="r", ls=":", lw=1)
    ax[1].set_xlabel("constriction ratio  (widest gate gap / cavity size)")
    ax[1].set_ylabel("rho"); ax[1].set_title(f"accessibility: rho vs gate  (corr={cc:+.2f})")
    # (C) trends vs S
    fS = [float((np.concatenate([r["rho"] for r in res if r["S"] == S]) < 0.5).mean()) for S in Ss]
    def_S = [np.mean([r["deficit"] for r in res if r["S"] == S]) for S in Ss]
    uf_S = [np.mean([r["starved_underfill"] for r in res if r["S"] == S]) for S in Ss]
    ax[2].plot(Ss, fS, "o-", label="f_starved (rho<0.5)")
    ax[2].plot(Ss, def_S, "s-", label="FG deficit (phi_FG-phi_sim)")
    ax[2].plot(Ss, uf_S, "^-", label="starved under-fill (localized)")
    ax[2].set_xlabel("size ratio S"); ax[2].set_title("trends vs S"); ax[2].legend(fontsize=8)
    fig.tight_layout(); out = FIG / "porefill_distributions.png"
    fig.savefig(out, dpi=130); plt.close()
    print(f"\nfigure -> {out}\n  accessibility corr(rho, constriction) = {cc:+.3f}")
    return cc

def main():
    HERE = pathlib.Path(__file__).resolve().parent
    DATA = HERE.parent / "Data" / "voidfill"; FIG = HERE.parent / "Figures"
    print(f"# per-cavity fill vs S (frac={FRAC}, Dc/Dmax={DCFRAC})")
    print(f"# {'pack':14s} {'S':>4} {'ncav':>5} {'<rho>':>6} {'med':>5} {'f_rho<.5':>8} "
          f"{'f_empty':>7} {'deficit':>7} {'starve_uf':>9}  {'corr_rc':>7}")
    res = []
    for fam in ("pl_p355", "pl_p363"):
        for f in sorted(glob.glob(str(DATA / f"{fam}_S*.npz")),
                        key=lambda x: int(x.split('_S')[-1].split('.')[0])):
            r = per_cavity(f); res.append(r); rho = r["rho"]
            fstarve = float((rho < 0.5).mean()); fempty = float((rho < 0.05).mean())
            mc = (r["constr"] > 0) & np.isfinite(rho)
            crc = np.corrcoef(rho[mc], r["constr"][mc])[0, 1] if mc.sum() > 3 else float("nan")
            print(f"  {pathlib.Path(f).stem:14s} {r['S']:4.0f} {r['n_cav']:5d} "
                  f"{rho.mean():6.3f} {np.median(rho):5.3f} {fstarve:8.3f} {fempty:7.3f} "
                  f"{r['deficit']:7.4f} {r['starved_underfill']:9.4f}  {crc:+7.2f}")
    make_figs(res, FIG)

if __name__ == "__main__":
    main()
