#!/usr/bin/env python3
"""Generate a 250k 2D lognormal RCP packing with the broadest width that still keeps
Dmax/L ~ 0.40 (largest particle ~40% of the box), then save an .npz. A fixed normal
sample is rescaled by sigma so Dmax/L is monotonic in sigma -> binary-search sigma."""
import sys, time, pathlib, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "Scripts"))
import rcpgenerator as r

N = 250_000; PHI_EST = 0.84; TARGET = 0.40
OUT = pathlib.Path(__file__).resolve().parent.parent / "Data" / "voidfill" / "logn2d_250k.npz"

def main():
    # Broad lognormal, truncated to a finite size ratio S. NOTE: a 250k lognormal
    # gives Dmax/L ~ 0.003 regardless of width (volume is shared across all N) -- the
    # requested Dmax/L=0.4 needs a power-law/bimodal (volume in a few bigs), not a
    # lognormal. This bundle is the broad-lognormal zoom demo (size range = S).
    S = float(sys.argv[1]) if len(sys.argv) > 1 else 60.0
    sig = float(sys.argv[2]) if len(sys.argv) > 2 else 0.85
    rng = np.random.default_rng(7)
    D = np.exp(sig * rng.standard_normal(N)); D = np.clip(D / D.min(), 1.0, S)
    L = np.sqrt((np.pi / 4 * D ** 2).sum() / PHI_EST)
    print(f"[gen] sigma={sig}  S(trunc)={S:.0f}  realized S={D.max()/D.min():.1f}  "
          f"Dmax/L_est={D.max()/L:.4f}", flush=True)
    t = time.time()
    pk = r.Packing(N=N, Ndim=2, phi=0.05, seed=0, walls=[0, 0],
                   dist={"type": "custom", "custom": list(map(float, D))})
    pk.initialize(); print(f"[gen] initialized in {time.time()-t:.1f}s", flush=True)
    t = time.time(); pk.pack(verbose=False, progress_interval=5000)
    pos = np.asarray(pk.positions, float); dia = np.asarray(pk.diameters, float); box = np.asarray(pk.box, float)
    np.savez(OUT, positions=pos, diameters=dia, box=box, Ndim=2,
             phi_final=float(pk.phi_final), phi_corr=float(pk.phi_final), seed=0,
             note="2D lognormal 250k, broad, Dmax/L~0.4")
    print(f"[gen] packed in {time.time()-t:.1f}s  phi={pk.phi_final:.4f}  "
          f"S={dia.max()/dia.min():.1f}  Dmax/L={dia.max()/box[0]:.3f}  N={len(dia)}", flush=True)
    print(f"[gen] saved -> {OUT}", flush=True)

if __name__ == "__main__":
    main()
