"""Diameter permutation on a 3-D power-law distribution: the same particles packed with the
permutation off (pristine algorithm) and on (validated default recipe), and once more with a
reduced budget, to show the main dial (permute_budget) at work.

Run:  python example_3d_diameter_permute.py
Takes a few minutes at one thread.
"""
import math

import numpy as np

import rcpgenerator

from _common import build_summary
from _common import print_case_summary


CASE_NAME = "3d_diameter_permute"
N, S, P = 4000, 20.0, -3.45      # particles, D_max/D_min, number-density exponent


def powerlaw_diameters(n: int, size_ratio: float, exponent: float, seed: int = 0) -> list[float]:
    """n diameters with number density ~ D**exponent between 1 and size_ratio (deterministic quantiles)."""
    a = exponent + 1.0
    u = (np.arange(n) + 0.5) / n
    d = (1.0 + u * (size_ratio ** a - 1.0)) ** (1.0 / a)
    rng = np.random.default_rng(seed)
    return rng.permutation(d).tolist()


def build(**permute_kwargs) -> rcpgenerator.Packing:
    return rcpgenerator.Packing(
        phi=0.05,
        N=N,
        Ndim=3,
        box=[1.0, 1.0, 1.0],
        walls=[0, 0, 0],
        dist={"type": "custom", "custom": powerlaw_diameters(N, S, P)},
        neighbor_max=0,
        seed=1,
        verbose=False,
        **permute_kwargs,
    )


def phi_estimate(packing: rcpgenerator.Packing) -> float:
    """Packing fraction after shrinking every diameter by the residual maximum overlap."""
    return packing.phi_history[-1] * (1.0 - packing.max_overlap_history[-1]) ** 3


def run_case(verbose: bool = False) -> dict[str, object]:
    arms = {
        "off": dict(),
        "on (default: budget 0.25, start 4000)": dict(diameter_permute=True),
        "on, budget 0.03": dict(diameter_permute=True, permute_budget=0.03),
    }
    results = {}
    for name, kw in arms.items():
        packing = build(**kw)
        packing.pack()
        results[name] = phi_estimate(packing)
        if verbose:
            print(f"{name:<40} phi = {results[name]:.4f}   steps = {packing.steps}")
        last = packing
    print_case_summary(CASE_NAME, last)
    summary = build_summary(CASE_NAME, last)
    summary["phi_by_arm"] = results
    return summary


def main() -> None:
    run_case(verbose=True)


if __name__ == "__main__":
    main()
