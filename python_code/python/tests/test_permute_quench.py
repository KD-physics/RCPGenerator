"""Functional tests for diameter permutation and the quench to force balance.

These check stability and intended behaviour across dimensions and boundary conditions;
they are not packing-fraction studies. Every case uses a few hundred particles at one
thread.

    python tests/test_permute_quench.py           # full battery (a few minutes)
    python tests/test_permute_quench.py --quick   # representative subset (well under a minute)
    pytest tests/test_permute_quench.py           # the same tests under pytest; RCP_QUICK=1 for the subset
"""

from __future__ import annotations

import math
import os
import sys
import time
import warnings

import numpy as np

import rcpgenerator

rcpgenerator.set_num_threads(1)   # one thread: bit-reproducible runs (the determinism test relies on it)

N = 500
QUICK = os.environ.get("RCP_QUICK", "") not in ("", "0")

# Quench fixtures use an explicit compression of x1.01 rather than the 2D/3D default x1.003. A quench gives up when its
# worst residual has not halved for stall_window iterations; at x1.003 that happens for a few percent of realisations
# whose minimisation is still converging (3D periodic N=500 S=5: 2 of 40 seeds; resumed, they converge). A different
# compiler or CPU (FMA contraction, sort tie order) is a different realisation of the same fixture, so these
# tests must not depend on one realisation. At x1.01 the same fixtures converged in every one of 380 seeded runs
# (single quench, resume and shear paths). The default value itself is checked where convergence does not matter.
FIXTURE_KAPPA = 1.01

# name -> (Ndim, walls). walls[d] = 1: hard walls at 0 and L along d; walls[0] = -t: the first t
# dimensions form a disk (t=2) or sphere (t=3) of diameter box[0]; 0: periodic.
GEOMETRIES = {
    "2d_periodic": (2, [0, 0]),
    "2d_walls_1dir": (2, [0, 1]),
    "2d_walls_2dir": (2, [1, 1]),
    "2d_disk": (2, [-2, 0]),
    "3d_periodic": (3, [0, 0, 0]),
    "3d_walls_1dir": (3, [0, 0, 1]),
    "3d_walls_2dir": (3, [0, 1, 1]),
    "3d_walls_3dir": (3, [1, 1, 1]),
    "3d_cylinder": (3, [-2, 0, 0]),
    "3d_cylinder_capped": (3, [-2, 0, 1]),
    "3d_sphere": (3, [-3, 0, 0]),
}
QUICK_GEOMETRIES = ("2d_walls_2dir", "3d_periodic", "3d_sphere")
ARMS = {
    "off": {},
    "permute_default": {"diameter_permute": True},
    "permute_low_budget": {"diameter_permute": True, "permute_budget": 0.01},
}


def powerlaw(n: int, size_ratio: float, exponent: float = -2.5, seed: int = 0) -> list[float]:
    """n deterministic power-law diameters on [1, size_ratio], shuffled."""
    a = exponent + 1.0
    u = (np.arange(n) + 0.5) / n
    d = (1.0 + u * (size_ratio ** a - 1.0)) ** (1.0 / a)
    return np.random.default_rng(seed).permutation(d).tolist()


def make(ndim: int = 3, walls=None, n: int = N, size_ratio: float = 5.0, seed: int = 3, **kw) -> rcpgenerator.Packing:
    walls = [0] * ndim if walls is None else list(walls)
    return rcpgenerator.Packing(phi=0.05 if ndim == 3 else 0.1, N=n, Ndim=ndim, box=[1.0] * ndim, walls=walls,
                                dist={"type": "custom", "custom": powerlaw(n, size_ratio)}, seed=seed, **kw)


def same_multiset_up_to_scale(a, b, rel: float = 1e-12) -> float:
    """Return the common scale factor between two diameter multisets; fail if there is none."""
    ra = np.sort(np.asarray(a, float)) / np.sort(np.asarray(b, float))
    assert np.ptp(ra) <= rel * ra.mean(), f"diameter multiset changed beyond a common scale (spread {np.ptp(ra):.2e})"
    return float(ra.mean())


def quench_quiet(pk: rcpgenerator.Packing, **kw) -> list:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        pk.quench(**kw)
    return [w for w in caught if issubclass(w.category, RuntimeWarning)]


def straight_line_balance(pk: rcpgenerator.Packing) -> float:
    """The quench's residual recomputed independently: max |F_i| / S_i over particles with >= Ndim + 1 contacts and a
    contact-force sum >= 4e-2 x the mean contact force. Pair distances wrap only along periodic dimensions, never
    through a wall; wall contacts (planar or curved container) push with 2 x overlap."""
    x, d, L = np.asarray(pk.positions, float), np.asarray(pk.diameters, float), np.asarray(pk.box, float)
    walls = list(pk.walls)
    n, nd = x.shape
    t = -walls[0] if walls[0] < 0 else 0
    periodic = np.array([w == 0 and k >= t for k, w in enumerate(walls)])
    F = np.zeros_like(x)
    S = np.zeros(n)
    z = np.zeros(n, int)
    for i in range(n):
        dx = x[i] - x
        dx -= np.round(dx / L) * L * periodic
        r = np.linalg.norm(dx, axis=1)
        r[i] = np.inf
        ov = (d[i] + d) / 2.0 - r
        m = ov > 0.0
        F[i] = ((ov[m] / r[m])[:, None] * dx[m]).sum(axis=0)
        S[i], z[i] = ov[m].sum(), m.sum()
        for k in range(nd):
            if walls[k] == 1:
                for side, sign in ((0.0, 1.0), (L[k], -1.0)):
                    o = d[i] / 2.0 - abs(x[i, k] - side)
                    if o > 0.0:
                        F[i, k] += sign * 2.0 * o
                        S[i] += 2.0 * o
                        z[i] += 1
        if t:
            u = x[i, :t] - L[0] / 2.0
            dc = np.linalg.norm(u)
            o = d[i] / 2.0 - (L[0] / 2.0 - dc)
            if o > 0.0:
                F[i, :t] -= 2.0 * o * u / dc
                S[i] += 2.0 * o
                z[i] += 1
    fmean = S[S > 0].sum() / z[S > 0].sum()
    backbone = (z >= nd + 1) & (S >= 4e-2 * fmean)
    assert backbone.any()
    return float((np.linalg.norm(F, axis=1)[backbone] / S[backbone]).max())


def check_state(pk: rcpgenerator.Packing) -> None:
    x = np.asarray(pk.positions, float)
    d = np.asarray(pk.diameters, float)
    assert x.shape == (len(d), pk.Ndim)
    assert np.isfinite(x).all() and np.isfinite(d).all() and (d > 0).all()
    box = np.asarray(pk.box, float)
    assert ((x >= 0.0) & (x < box + 1e-12)).all(), "a particle left the box"


# --------------------------------------------------------------------------- the matrix

def run_case(geometry: str, arm: str) -> dict:
    ndim, walls = GEOMETRIES[geometry]
    pk = make(ndim, walls, **ARMS[arm])
    d_init = list(pk.diameters)
    t0 = time.time()
    pk.pack()
    t_pack = time.time() - t0
    check_state(pk)
    same_multiset_up_to_scale(pk.diameters, d_init)

    passes, exchanges, counts = pk.permute_passes, pk.permute_exchanges, list(pk.permute_counts)
    if arm == "off":
        assert passes == 0 and exchanges == 0, (passes, exchanges)
    else:
        assert passes > 0 and exchanges > 0, f"{geometry}/{arm}: permutation did not run"
        assert len(counts) == len(d_init) and sum(counts) > 0
        cap = max(1, int(pk.permute_budget * len(d_init)))
        assert exchanges <= 2 * cap * passes, "more exchanges than the budget allows (pairwise + 3-cycle caps)"

    d_packed = list(pk.diameters)
    x_packed = np.asarray(pk.positions, float)
    t0 = time.time()
    caught = quench_quiet(pk, kappa_mult=FIXTURE_KAPPA)
    t_quench = time.time() - t0
    check_state(pk)
    assert pk.quench_status == "converged", f"{geometry}/{arm}: {pk.quench_message}"
    assert pk.quench_converged == 1 and pk.quench_residual <= 1e-10 and not caught
    if any(w != 0 for w in walls):
        # no contact through a wall: the balance holds with straight-line distances across walled dimensions
        bal = straight_line_balance(pk)
        assert bal <= 1e-8, f"{geometry}/{arm}: straight-line force balance {bal:.2e} (quench reports {pk.quench_residual:.2e})"
    scale = same_multiset_up_to_scale(pk.diameters, d_packed)
    assert abs(scale - FIXTURE_KAPPA) < 1e-12, f"diameters rescaled by {scale}, expected {FIXTURE_KAPPA} once"
    # particle i keeps index i through the quench, and the pack's permutation diagnostics survive
    assert np.allclose(np.asarray(pk.diameters) / np.asarray(d_packed), FIXTURE_KAPPA, rtol=0, atol=1e-12)
    box = np.asarray(pk.box, float)
    dx = np.asarray(pk.positions) - x_packed
    dx -= np.round(dx / box) * box                    # minimum image across periodic seams
    assert np.max(np.abs(dx)) < 0.1 * box.min(), "a particle's index changed across the quench"
    assert (pk.permute_passes, pk.permute_exchanges, list(pk.permute_counts)) == (passes, exchanges, counts)
    return {"pack_s": t_pack, "quench_s": t_quench, "quench_iters": pk.quench_iters, "exchanges": exchanges}


def test_matrix():
    geos = QUICK_GEOMETRIES if QUICK else tuple(GEOMETRIES)
    arms = ("off", "permute_default") if QUICK else tuple(ARMS)
    for g in geos:
        for a in arms:
            run_case(g, a)


# --------------------------------------------------------------------------- quench API

def test_quench_resume_and_give_up():
    base = make()
    base.pack()
    saved = base.to_dict()
    d_packed = list(base.diameters)

    pk = rcpgenerator.Packing.from_dict(saved)
    caught = quench_quiet(pk, max_iter=200, kappa_mult=FIXTURE_KAPPA)
    assert pk.quench_status == "max_iter" and pk.quench_iters == 200 and caught
    assert "resume=True" in pk.quench_message
    assert abs(same_multiset_up_to_scale(pk.diameters, d_packed) - FIXTURE_KAPPA) < 1e-12
    d_first = list(pk.diameters)
    caught = quench_quiet(pk, resume=True)
    assert pk.quench_status == "converged" and not caught, pk.quench_message
    assert list(pk.diameters) == d_first, "resume=True must not rescale"

    try:
        pk.quench(resume=True, kappa_mult=1.01)
        raise AssertionError("resume=True with kappa_mult must raise")
    except ValueError:
        pass

    pk = rcpgenerator.Packing.from_dict(saved)       # shrunk: nothing left to balance
    x0 = np.asarray(pk.positions)
    caught = quench_quiet(pk, kappa_mult=0.99)
    assert pk.quench_status == "unjammed" and caught and pk.quench_converged == 0
    assert np.isfinite(np.asarray(pk.positions)).all() and np.asarray(pk.positions).shape == x0.shape

    pk = rcpgenerator.Packing.from_dict(saved)       # a stall window of one iteration gives up at once
    x0 = np.asarray(pk.positions)
    caught = quench_quiet(pk, stall_window=1)
    assert pk.quench_status == "stalled" and caught and pk.quench_iters < 10, (pk.quench_status, pk.quench_iters)
    assert not np.array_equal(np.asarray(pk.positions), x0), "the stopped state must be kept"
    assert pk.quench_message.startswith("gave up")
    # the default compression in 2D/3D is x1.003, applied once (checked here: it does not depend on convergence)
    assert abs(same_multiset_up_to_scale(pk.diameters, d_packed) - 1.003) < 1e-12


def test_quench_held_particles():
    for ndim, walls in ((3, [0, 0, 0]), (2, [1, 1])):
        pk = make(ndim, walls)
        pk.pack()
        x0 = np.asarray(pk.positions, float)
        held = [i for i in range(len(x0)) if x0[i, 0] < 0.15 * pk.box[0]]
        assert held
        quench_quiet(pk, held=held)
        x1 = np.asarray(pk.positions, float)
        assert np.array_equal(x1[held], x0[held]), "held particles moved"
        assert len(pk.quench_held_force) == ndim and np.isfinite(pk.quench_held_force).all()
        free = np.setdiff1d(np.arange(len(x0)), held)
        assert not np.array_equal(x1[free], x0[free])


def slab_force_from_free(pk: rcpgenerator.Packing, members, held) -> np.ndarray:
    """Contact force (overlap along the line of centres) from the free particles on `members`, periodic box."""
    x, d, L = np.asarray(pk.positions, float), np.asarray(pk.diameters, float), np.asarray(pk.box, float)
    free = np.ones(len(d), bool)
    free[held] = False
    others = np.flatnonzero(free)
    total = np.zeros(x.shape[1])
    for i in members:
        dx = x[others] - x[i]
        dx -= np.round(dx / L) * L
        dist = np.linalg.norm(dx, axis=1)
        ov = (d[others] + d[i]) / 2.0 - dist
        m = ov > 0.0
        total -= ((ov[m] / dist[m])[:, None] * dx[m]).sum(axis=0)
    return total


def test_quench_shear_response():
    """Held-slab shear on a balanced packing: the load is transmitted through the packing (equal and opposite
    on the two slabs) and is linear in the displacement."""
    pk = make(n=N)
    pk.pack()
    quench_quiet(pk, kappa_mult=FIXTURE_KAPPA)
    assert pk.quench_status == "converged", pk.quench_message
    base = pk.to_dict()
    x0 = np.asarray(pk.positions, float)
    moved = np.flatnonzero(x0[:, 0] < 0.15)
    fixed = np.flatnonzero((x0[:, 0] >= 0.5) & (x0[:, 0] < 0.65))
    held = np.concatenate([moved, fixed])
    fy = []
    for strain in (1e-5, 2e-5):
        s = rcpgenerator.Packing.from_dict(base)
        x = x0.copy()
        x[moved, 1] = (x[moved, 1] + strain) % 1.0
        s.positions = x.tolist()
        caught = quench_quiet(s, resume=True, held=held.tolist())
        assert s.quench_status == "converged" and not caught, s.quench_message
        f_m, f_f = slab_force_from_free(s, moved, held)[1], slab_force_from_free(s, fixed, held)[1]
        assert f_m < 0.0, "the displaced slab must be pushed back"
        assert abs(f_m + f_f) <= 1e-3 * abs(f_m), f"slab forces not equal and opposite: {f_m:.4e} vs {f_f:.4e}"
        fy.append(f_m)
    assert abs(fy[1] / fy[0] - 2.0) < 0.1, f"response not linear: ratio {fy[1] / fy[0]:.3f}"


def test_no_contact_through_walls():
    """Hard walls are not periodic. In 1.0.0 every pair distance was wrapped through walls too, so particles touching
    opposite walls could push on each other through them, and a 'converged' quench was out of balance by up to ~0.6
    on those particles once distances were measured in a straight line. Here the balance must hold that way."""
    for ndim, walls in ((3, [0, 0, 1]), (3, [1, 1, 1])):
        pk = make(ndim, walls)
        pk.pack()
        caught = quench_quiet(pk, kappa_mult=FIXTURE_KAPPA)
        assert pk.quench_status == "converged" and not caught, f"walls={walls}: {pk.quench_message}"
        bal = straight_line_balance(pk)
        assert bal <= 1e-8, f"walls={walls}: straight-line force balance {bal:.2e} (quench reports {pk.quench_residual:.2e})"


def test_quench_needs_a_state_and_validates():
    pk = make(n=200)
    # the constructor already places particles: quenching that dilute state is allowed and reports unjammed
    caught = quench_quiet(pk)
    assert pk.quench_status == "unjammed" and caught
    pk.reset()                                        # no realized state at all
    try:
        pk.quench()
        raise AssertionError("quench() without a realized state must raise")
    except RuntimeError:
        pass
    pk.pack()
    for kwargs in ({"mix_mode": 2}, {"held": [10 ** 6]}, {"stall_window": -1}, {"tol": 0.0}, {"kappa_mult": -1.0}):
        try:
            pk.quench(**kwargs)
            raise AssertionError(f"quench({kwargs}) must raise ValueError")
        except ValueError:
            pass


# --------------------------------------------------------------------------- permutation API

def test_permute_validation():
    for kwargs in ({"permute_budget": 0.0}, {"permute_budget": 0.6}, {"permute_interval": 0}, {"permute_start": -1}):
        pk = make(n=100, diameter_permute=True, **kwargs)
        try:
            pk.pack()
            raise AssertionError(f"{kwargs} must raise ValueError")
        except ValueError:
            pass
    try:
        make(n=100, diameter_permutation=True)        # misspelt keyword
        raise AssertionError("unknown keyword must raise TypeError")
    except TypeError:
        pass


def test_off_ignores_permute_knobs_and_is_deterministic():
    runs = []
    for kw in ({}, {"diameter_permute": False, "permute_budget": 0.1, "permute_interval": 2, "permute_start": 10}, {}):
        pk = make(n=300, **kw)
        pk.pack()
        runs.append((pk.steps, np.asarray(pk.positions).tobytes(), np.asarray(pk.diameters).tobytes()))
    assert runs[0] == runs[1] == runs[2], "permutation off must not depend on the permutation knobs"
    on = []
    for _ in range(2):
        pk = make(n=300, diameter_permute=True)
        pk.pack()
        on.append((pk.steps, pk.permute_exchanges, np.asarray(pk.positions).tobytes()))
    assert on[0] == on[1], "permutation must be deterministic for a fixed seed at one thread"


def test_copy_and_serialisation():
    pk = make(n=300, diameter_permute=True, permute_budget=0.1)
    pk.pack()
    for clone in (pk.copy(), rcpgenerator.Packing.from_dict(pk.to_dict())):
        assert clone.diameter_permute and clone.permute_budget == 0.1
        assert clone.to_dict()["pack_config"] == pk.to_dict()["pack_config"]
        quench_quiet(clone, kappa_mult=FIXTURE_KAPPA)
        assert clone.quench_status == "converged", clone.quench_message
    off = pk.copy()
    off.diameter_permute = False
    assert off._collect_pack_config()["swap"] == 0.0, "turning permutation off on a copy must switch it off"
    # a dictionary written by version 1.0 (no permutation or quench keys) still loads, packs and quenches
    legacy = pk.to_dict()
    legacy["pack_config"] = {k: v for k, v in legacy["pack_config"].items()
                             if k in ("box", "walls", "neighbor_max", "seed", "fix_height")}
    old = rcpgenerator.Packing.from_dict(legacy)
    assert old.diameter_permute is False
    old.pack()
    quench_quiet(old, kappa_mult=FIXTURE_KAPPA)
    assert old.quench_status == "converged", old.quench_message


# --------------------------------------------------------------------------- runner

def test_quench_high_dimensions():
    """4D and 5D at the default compression (x1.01 above 3D): a hard-walled 4D quench converges and balances with
    straight-line distances, and a periodic 5D quench keeps every particle in the box (its translation-mean removal
    once overflowed a 3-slot array above 3D). The walled case uses a size ratio of 2: at a size ratio of 5, 4D hard
    walls give up early (the stall rule) in about 1 single quench call in 6 (33/40 seeds converged; 30/30 at ratio 2),
    which a cross-platform test cannot absorb."""
    cases = [(4, [1] * 4, 1000, 2.0), (5, [0] * 5, 1000, 1.0)]
    for ndim, walls, n, ratio in cases:
        pk = make(ndim, walls, n=n, size_ratio=ratio)
        pk.pack()
        d_packed = list(pk.diameters)
        caught = quench_quiet(pk)
        check_state(pk)
        tag = f"{ndim}D walls={walls} ratio={ratio}"
        assert pk.quench_status == "converged" and pk.quench_residual <= 1e-10 and not caught, f"{tag}: {pk.quench_message}"
        assert np.allclose(np.asarray(pk.diameters) / np.asarray(d_packed), 1.01, rtol=0, atol=1e-12), tag
        if walls[0] != 0:
            bal = straight_line_balance(pk)
            assert bal <= 1e-8, f"{tag}: straight-line force balance {bal:.2e}"


TESTS = [test_matrix, test_quench_resume_and_give_up, test_quench_held_particles, test_quench_shear_response,
         test_no_contact_through_walls, test_quench_needs_a_state_and_validates, test_permute_validation,
         test_off_ignores_permute_knobs_and_is_deterministic, test_copy_and_serialisation, test_quench_high_dimensions]


def main(argv: list[str]) -> int:
    global QUICK
    QUICK = QUICK or "--quick" in argv
    rcpgenerator.set_num_threads(1)
    failures = 0
    for test in TESTS:
        t0 = time.time()
        try:
            test()
            print(f"PASS  {test.__name__}  ({time.time() - t0:.1f} s)", flush=True)
        except Exception as exc:  # report every test, then fail
            failures += 1
            print(f"FAIL  {test.__name__}: {type(exc).__name__}: {exc}", flush=True)
    print(f"{len(TESTS) - failures}/{len(TESTS)} passed" + (" (quick)" if QUICK else ""))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
