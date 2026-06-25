#!/usr/bin/env python3
"""Secondary, scale-free void check: density-fluctuation spectrum sigma^2_phi(ell).
Tile the box at window size ell, local phi per cell (particle volume -> center cell;
valid for ell >> particle size, i.e. the LARGE-scale void-organization question),
variance across cells vs ell, swept. Excess / growth with S at some ell* = void
organization; smooth power-law decay = homogeneous. Caveat: S<=100 (below saturation).

Usage: python chord_boxbin.py
"""
import glob, pathlib, numpy as np
import voidfill_method as vf

def sigma2_phi(pos, dia, box, ms):
    vol = vf.sphere_vol(dia); Vbox = float(np.prod(box)); d = len(box)
    phi = vol.sum() / Vbox; out = []
    for m in ms:
        idx = np.floor((pos % box) / (box / m)).astype(int)
        idx = np.clip(idx, 0, m - 1)
        flat = np.ravel_multi_index([idx[:, k] for k in range(d)], (m,) * d)
        cellvol = np.bincount(flat, weights=vol, minlength=m ** d)
        cphi = cellvol / (Vbox / m ** d)
        out.append((box[0] / m, cphi.std() / phi))     # ell, normalized rms fluctuation
    return phi, np.array(out)

def main():
    DATA = pathlib.Path(__file__).resolve().parent.parent / "Data" / "voidfill"
    FIG = pathlib.Path(__file__).resolve().parent.parent / "Figures"
    ms = np.array([2, 3, 4, 5, 6, 8, 10, 13, 16, 20, 25, 32])
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(13, 5)); cmap = plt.cm.viridis
    print(f"# sigma_phi/phi at fixed large window (ell~L/4) and slope, vs S")
    print(f"# {'pack':13s} {'S':>4} {'sig@L/4':>8} {'sig@L/8':>8} {'slope':>6}")
    fams = ["pl_p355", "pl_p363"]
    for fam in fams:
        files = sorted(glob.glob(str(DATA / f"{fam}_S*.npz")),
                       key=lambda x: int(x.split('_S')[-1].split('.')[0]))
        for f in files:
            pos, dia, box = vf._load(f); S = dia.max() / dia.min()
            phi, sf = sigma2_phi(pos, dia, box, ms)
            ell, sig = sf[:, 0], sf[:, 1]
            slope = np.polyfit(np.log(ell), np.log(sig), 1)[0]
            s4 = np.interp(box[0] / 4, ell[::-1], sig[::-1]); s8 = np.interp(box[0] / 8, ell[::-1], sig[::-1])
            print(f"  {pathlib.Path(f).stem:13s} {S:4.0f} {s4:8.4f} {s8:8.4f} {slope:6.2f}")
            col = cmap(files.index(f) / max(1, len(files) - 1))
            axi = ax[fams.index(fam)]
            axi.loglog(ell / box[0], sig, "o-", color=col, label=f"S={S:.0f}")
        ax[fams.index(fam)].set_xlabel("window ell / L"); ax[fams.index(fam)].set_ylabel("sigma_phi/phi")
        ax[fams.index(fam)].set_title(f"{fam}: density fluctuations"); ax[fams.index(fam)].legend(fontsize=8)
    fig.suptitle("Box-binning sigma^2_phi(ell): voids would show as growing large-scale excess vs S", fontsize=12)
    fig.tight_layout(); out = FIG / "chord_boxbin.png"; fig.savefig(out, dpi=130); plt.close()
    print(f"figure -> {out}")

if __name__ == "__main__":
    main()
