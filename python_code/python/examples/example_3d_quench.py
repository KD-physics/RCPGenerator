"""Quench to force balance on a 3-D power-law packing.

1. Pack (with diameter permutation on) and quench the result at the default settings.
2. Quench the same packed state in two stages: a loose tolerance first, then
   resume=True with a tight one.
3. One quasi-static shear step: displace a slab of particles, hold it and a second slab
   half a box away, quench, and resolve the force the free particles exert on each slab
   (equal and opposite: the shear stress is carried through the packing).

Run:  python example_3d_quench.py
Takes about a minute at one thread.
"""
import warnings

import numpy as np

import rcpgenerator

from _common import build_summary
from _common import print_case_summary


CASE_NAME = "3d_quench"
N, S, P = 2000, 8.0, -3.0      # particles, D_max/D_min, number-density exponent


def powerlaw_diameters(n: int, size_ratio: float, exponent: float, seed: int = 0) -> list[float]:
    """n diameters with number density ~ D**exponent between 1 and size_ratio (deterministic quantiles)."""
    a = exponent + 1.0
    u = (np.arange(n) + 0.5) / n
    d = (1.0 + u * (size_ratio ** a - 1.0)) ** (1.0 / a)
    rng = np.random.default_rng(seed)
    return rng.permutation(d).tolist()


def report(label: str, packing: rcpgenerator.Packing) -> None:
    print(f"{label:<28} status = {packing.quench_status:<10} residual = {packing.quench_residual:.2e}"
          f"   iterations = {packing.quench_iters}")


def force_from_free(packing: rcpgenerator.Packing, members: np.ndarray, held: np.ndarray) -> np.ndarray:
    """Total contact force the free (non-held) particles exert on `members`, periodic box.

    A contact's force is its overlap, (D_i + D_j)/2 - r_ij, along the line of centres
    (the engine's force law). Packing.quench_held_force gives only the total over ALL held
    particles; resolving one slab needs this sum over its contacts."""
    x = np.asarray(packing.positions, float)
    d = np.asarray(packing.diameters, float)
    L = np.asarray(packing.box, float)
    free = np.ones(len(d), bool)
    free[held] = False
    others = np.flatnonzero(free)
    total = np.zeros(x.shape[1])
    for i in members:
        dx = x[others] - x[i]
        dx -= np.round(dx / L) * L                        # minimum image
        dist = np.linalg.norm(dx, axis=1)
        overlap = (d[others] + d[i]) / 2.0 - dist
        touching = overlap > 0.0
        total -= ((overlap[touching] / dist[touching])[:, None] * dx[touching]).sum(axis=0)
    return total


def shear_step(balanced: rcpgenerator.Packing, strain: float = 1e-5) -> tuple[np.ndarray, np.ndarray]:
    """Displace the slab x < 0.1 L by strain * L along y, hold it and the slab 0.5 L <= x < 0.6 L,
    quench the rest, and return the force from the free particles on each slab."""
    sheared = balanced.copy()
    x = np.asarray(sheared.positions, float)
    L = np.asarray(sheared.box, float)
    moved = np.flatnonzero(x[:, 0] < 0.1 * L[0])
    fixed = np.flatnonzero((x[:, 0] >= 0.5 * L[0]) & (x[:, 0] < 0.6 * L[0]))
    held = np.concatenate([moved, fixed])
    x[moved, 1] = (x[moved, 1] + strain * L[1]) % L[1]
    sheared.positions = x.tolist()
    sheared.quench(resume=True, held=held.tolist())   # particle i keeps index i through the quench
    report("shear step", sheared)
    return force_from_free(sheared, moved, held), force_from_free(sheared, fixed, held)


def run_case(verbose: bool = False) -> dict[str, object]:
    packing = rcpgenerator.Packing(
        phi=0.05,
        N=N,
        Ndim=3,
        box=[1.0, 1.0, 1.0],
        walls=[0, 0, 0],
        dist={"type": "custom", "custom": powerlaw_diameters(N, S, P)},
        seed=1,
        verbose=False,
        diameter_permute=True,
    )
    packing.pack()
    packed = packing.copy()                       # keep the packed state for the other demonstrations

    # 1. default quench: diameters x 1.003 once, then FIRE to a residual of 1e-10
    packing.quench()
    report("default quench", packing)
    if verbose:
        print("   ", packing.quench_message)

    # 2. staged: the second call continues from the first without rescaling the diameters
    staged = packed.copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)   # stage 1 may stop early; we continue below
        staged.quench(tol=1e-4)
    report("stage 1 (tol 1e-4)", staged)
    staged.quench(tol=1e-10, resume=True)
    report("stage 2 (tol 1e-10, resume)", staged)

    # 3. a shear step on the balanced state: the force on the displaced slab is the transmitted shear load
    f_moved, f_fixed = shear_step(packing)
    stress = f_moved[1] / (2.0 * packing.box[0] * packing.box[2])   # the slab has two faces
    if verbose:
        print(f"    y-force on displaced slab {f_moved[1]:+.4e}, on fixed slab {f_fixed[1]:+.4e}; shear stress {stress:+.3e}")

    print_case_summary(CASE_NAME, packing)
    summary = build_summary(CASE_NAME, packing)
    summary["quench_status"] = packing.quench_status
    summary["quench_residual"] = packing.quench_residual
    return summary


def main() -> None:
    run_case(verbose=True)


if __name__ == "__main__":
    main()
