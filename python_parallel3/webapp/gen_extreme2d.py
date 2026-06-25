#!/usr/bin/env python3
"""Extreme multiscale CONTINUOUS lognormal 2D packing using the NATIVE lognormal dist
(dist={"type":"lognormal","mu","sigma"}, neighbor_max=0) -- the efficient path the 3D
validation used (N~1e5, S~1e4 at <1GB). The earlier `custom` explicit-diameter route
was what blew memory at large S, not the generator.
  python gen_extreme2d.py [N] [sigma]   (defaults N=250000, sigma=1.0 -> S~1e4)
"""
import sys, time, pathlib, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "Scripts"))
import rcpgenerator as r

# working config: N=100k, sigma=1.3 -> S~8400, phi~0.91 (2D memory cliff is N*S~1.5e9,
# so S~1e4 needs N<=~120k; native lognormal dist, NOT the `custom` route which blew up).
N = int(sys.argv[1]) if len(sys.argv) > 1 else 100_000
SIGMA = float(sys.argv[2]) if len(sys.argv) > 2 else 1.3
OUT = pathlib.Path(__file__).resolve().parent.parent / "Data" / "voidfill" / "logn2d_extreme.npz"

def main():
    print(f"[gen] native lognormal  N={N}  Ndim=2  sigma={SIGMA}", flush=True)
    t = time.time()
    pk = r.Packing(N=N, Ndim=2, phi=0.05, seed=0, walls=[0, 0], fix_height=False,
                   dist={"type": "lognormal", "mu": 0.0, "sigma": SIGMA}, neighbor_max=0)
    pk.initialize(); print(f"[gen] initialized in {time.time()-t:.1f}s", flush=True)
    t = time.time(); pk.pack(verbose=False, progress_interval=5000)
    pos = np.asarray(pk.positions, float); dia = np.asarray(pk.diameters, float); box = np.asarray(pk.box, float)
    np.savez(OUT, positions=pos, diameters=dia, box=box, Ndim=2,
             phi_final=float(pk.phi_final), phi_corr=float(pk.phi_final), seed=0,
             note=f"2D native lognormal sigma={SIGMA}")
    print(f"[gen] packed in {time.time()-t:.1f}s  phi={pk.phi_final:.4f}  "
          f"S={dia.max()/dia.min():.0f}  Dmax/L={dia.max()/box[0]:.3f}  N={len(dia)}", flush=True)
    print(f"[gen] saved -> {OUT}", flush=True)

if __name__ == "__main__":
    main()
