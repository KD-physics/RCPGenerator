"""Fast installed-wheel smoke test for ellipsoids and rods.

Imports only the installed wheel and NumPy. One thread, fixed seeds. Each shape is packed once in a periodic box
and then checked with geometry written here, independently of the package:

* sphere                  pair distances
* prolate / oblate        the Perram-Wertheim contact function of two ellipsoids
* rod with flat ends, 2D  the separating-axis test of two rectangles
* rod with round ends, 3D the distance between two segments

For every shape: at the reported contact state (diameters x shrink_factor) no pair overlaps, a slightly larger
state does overlap (the shrink is tight; asked only when the packed state overlapped, i.e. shrink_factor < 1: a
run that ends with no force left has nothing to shrink and reports the state as packed), phi_contact equals the particle volume computed here, arrays have the
documented shapes and are finite, orientations are rotation matrices, and the state survives a JSON round trip.
One case is packed twice to show that a single-thread run is reproducible bit for bit.

    python wheel_smoke_shapes.py
"""

import json
import math

import numpy as np

import rcpgenerator

TIGHT = 1e-4          # a state this much larger than the reported contact state must overlap


# --------------------------------------------------------------------------- independent geometry ----
def unit_ball_volume(n: int) -> float:
    return math.pi ** (n / 2) / math.gamma(n / 2 + 1)


def close_pairs(x, reach, box):
    """All pairs (i < j) whose centres, by the minimum image, are closer than reach_i + reach_j."""
    n_part = len(x)
    i, j = np.triu_indices(n_part, 1)
    r = x[j] - x[i]
    r -= box * np.round(r / box)
    near = np.linalg.norm(r, axis=1) < reach[i] + reach[j]
    return i[near], j[near], r[near]


def sphere_gap(x, radius, box):
    """Smallest (distance - sum of radii) over all pairs."""
    i, j, r = close_pairs(x, 2.0 * radius, box)
    return float(np.min(np.linalg.norm(r, axis=1) - radius[i] - radius[j]))


def ellipsoid_contact(x, semi, rot, box):
    """Smallest Perram-Wertheim contact value over all pairs: max over lam in (0, 1) of
    lam (1 - lam) r . [(1 - lam) A + lam B]^-1 r, with A = R diag(semi^2) R^T. Two ellipsoids are apart iff it
    is above 1, touch at 1 and overlap below. The function of lam is concave: golden-section search."""
    shape = np.einsum("nij,nj,nkj->nik", rot, semi ** 2, rot)
    i, j, r = close_pairs(x, semi.max(axis=1), box)
    a, b = shape[i], shape[j]

    def value(lam):
        m = (1.0 - lam)[:, None, None] * a + lam[:, None, None] * b
        return lam * (1.0 - lam) * np.einsum("pi,pi->p", r, np.linalg.solve(m, r[:, :, None])[:, :, 0])

    lo, hi = np.zeros(len(r)), np.ones(len(r))
    g = (math.sqrt(5.0) - 1.0) / 2.0
    for _ in range(60):
        c, d = hi - g * (hi - lo), lo + g * (hi - lo)
        left = value(c) > value(d)
        hi = np.where(left, d, hi)
        lo = np.where(left, lo, c)
    return float(np.min(value(0.5 * (lo + hi))))


def segment_distance(r, u, hu, v, hv):
    """Distance between the segments {s u, |s| <= hu} and {r + t v, |t| <= hv} (unit u, v), for arrays of pairs.
    The squared distance is a convex quadratic in (s, t) on a rectangle, so its minimum is at the stationary point
    when that lies inside, and otherwise on one of the four edges, where the other variable has a closed form."""
    uv = np.einsum("pi,pi->p", u, v)
    ru, rv = np.einsum("pi,pi->p", r, u), np.einsum("pi,pi->p", r, v)

    def dist(s, t):
        return np.linalg.norm(r + t[:, None] * v - s[:, None] * u, axis=1)

    den = 1.0 - uv ** 2
    safe = np.where(den > 1e-12, den, 1.0)
    s0, t0 = (ru - uv * rv) / safe, (uv * ru - rv) / safe
    inside = (den > 1e-12) & (np.abs(s0) <= hu) & (np.abs(t0) <= hv)
    best = np.where(inside, dist(s0, t0), np.inf)
    for sign in (-1.0, 1.0):
        t = sign * hv                                   # an end of the second segment against the first
        best = np.minimum(best, dist(np.clip(ru + t * uv, -hu, hu), t))
        s = sign * hu                                   # an end of the first segment against the second
        best = np.minimum(best, dist(s, np.clip(s * uv - rv, -hv, hv)))
    return best


def capped_rod_gap(x, axis, half, radius, box):
    """Smallest (segment distance - sum of radii) over all pairs of rods with round ends."""
    i, j, r = close_pairs(x, half + radius, box)
    return float(np.min(segment_distance(r, axis[i], half[i], axis[j], half[j]) - radius[i] - radius[j]))


def rectangle_gap(x, axis, half, radius, box):
    """2D rods with flat ends are rectangles: the largest separation over the four edge normals (separating-axis
    test) is positive iff two rectangles are apart. Returns the smallest such value over all pairs."""
    i, j, r = close_pairs(x, np.hypot(half, radius), box)
    perp = np.stack([-axis[:, 1], axis[:, 0]], axis=1)
    best = np.full(len(r), -np.inf)
    for normal in (axis[i], perp[i], axis[j], perp[j]):
        def extent(k):
            return half[k] * np.abs(np.einsum("pi,pi->p", axis[k], normal)) + radius[k] * np.abs(np.einsum("pi,pi->p", perp[k], normal))
        best = np.maximum(best, np.abs(np.einsum("pi,pi->p", r, normal)) - extent(i) - extent(j))
    return float(np.min(best))


# -------------------------------------------------------------------------------------- the cases ----
def check_arrays(packing, n_part: int, ndim: int, shaped: bool) -> None:
    x, d = np.asarray(packing.positions), np.asarray(packing.diameters)
    assert x.shape == (n_part, ndim) and x.dtype == np.float64, (x.shape, x.dtype)
    assert d.shape == (n_part,) and d.dtype == np.float64, (d.shape, d.dtype)
    assert np.isfinite(x).all() and np.isfinite(d).all() and (d > 0).all()
    assert packing.steps > 0 and not packing.capped
    assert np.isfinite(np.asarray(packing.phi_history)).all()
    assert np.isfinite(np.asarray(packing.max_overlap_history)).all()
    if shaped:
        rot = np.asarray(packing.orientations)
        assert rot.shape == (n_part, ndim, ndim) and rot.dtype == np.float64, (rot.shape, rot.dtype)
        eye = np.einsum("nij,nkj->nik", rot, rot)
        assert np.abs(eye - np.eye(ndim)).max() < 1e-9, "orientations are not orthonormal"
        assert np.abs(np.linalg.det(rot) - 1.0).max() < 1e-9, "orientations are not proper rotations"
        assert 0.0 < packing.shrink_factor <= 1.0


def check_round_trip(packing, shaped: bool) -> None:
    restored = rcpgenerator.Packing.from_dict(json.loads(json.dumps(packing.to_dict())))
    assert np.array_equal(np.asarray(restored.positions), np.asarray(packing.positions))
    assert np.array_equal(np.asarray(restored.diameters), np.asarray(packing.diameters))
    if shaped:
        assert np.array_equal(np.asarray(restored.orientations), np.asarray(packing.orientations))
        assert restored.phi_contact == packing.phi_contact and restored.shrink_factor == packing.shrink_factor
    twin = packing.copy()
    assert np.array_equal(np.asarray(twin.positions), np.asarray(packing.positions))


def main() -> None:
    rcpgenerator.set_num_threads(1)
    box3, box2 = np.ones(3), np.ones(2)

    # ---- sphere: the unchanged path ------------------------------------------------------------------------
    spheres = rcpgenerator.Packing(N=200, Ndim=3, box=[1.0] * 3, walls=[0] * 3, seed=1)
    spheres.pack()
    check_arrays(spheres, 200, 3, shaped=False)
    state = spheres.shrink_to_contact()
    radius = 0.5 * np.asarray(spheres.diameters) * state["scale"]
    gap = sphere_gap(np.asarray(spheres.positions), radius, box3)
    assert gap > -1e-9, "spheres overlap at the contact state: %g" % gap
    assert state["scale"] == 1.0 or sphere_gap(np.asarray(spheres.positions), radius * (1 + TIGHT), box3) < 0.0, "sphere shrink is not tight"
    volume = float(np.sum(unit_ball_volume(3) * radius ** 3))
    assert abs(state["phi"] / volume - 1.0) < 1e-12
    assert 0.60 < state["phi"] < 0.66, state["phi"]
    check_round_trip(spheres, shaped=False)

    # ---- ellipsoids: prolate (one long axis) and oblate (two long axes) -------------------------------------
    first = None
    for name, ratios, band in (("prolate", [0.5, 0.5], (0.64, 0.72)), ("oblate", [1.0, 0.5], (0.64, 0.72))):
        ell = rcpgenerator.Packing(N=150, Ndim=3, box=[1.0] * 3, walls=[0] * 3, seed=2, aspect_ratios=ratios)
        ell.pack()
        check_arrays(ell, 150, 3, shaped=True)
        semi = 0.5 * np.asarray(ell.diameters)[:, None] * ell.shrink_factor * np.array([1.0] + ratios)[None, :]
        x, rot = np.asarray(ell.positions), np.asarray(ell.orientations)
        assert 2.0 * semi.max() < 0.5, "an ellipsoid is longer than half the box"
        contact = ellipsoid_contact(x, semi, rot, box3)
        assert contact > 1.0 - 1e-7, "%s ellipsoids overlap at the contact state: %g" % (name, contact)
        assert ell.shrink_factor == 1.0 or ellipsoid_contact(x, semi * (1 + TIGHT), rot, box3) < 1.0, "%s shrink is not tight" % name
        volume = float(np.sum(unit_ball_volume(3) * np.prod(semi, axis=1)))
        assert abs(ell.phi_contact / volume - 1.0) < 1e-12, (name, ell.phi_contact, volume)
        assert band[0] < ell.phi_contact < band[1], (name, ell.phi_contact)
        check_round_trip(ell, shaped=True)
        if first is None:
            first = ell

    # the same inputs again: a single-thread run is reproducible bit for bit
    again = rcpgenerator.Packing(N=150, Ndim=3, box=[1.0] * 3, walls=[0] * 3, seed=2, aspect_ratios=[0.5, 0.5])
    again.pack()
    assert again.steps == first.steps
    assert np.array_equal(np.asarray(again.positions), np.asarray(first.positions))
    assert np.array_equal(np.asarray(again.orientations), np.asarray(first.orientations))
    assert np.array_equal(np.asarray(again.diameters), np.asarray(first.diameters))

    # ---- rods with flat ends, 2D (rectangles) ---------------------------------------------------------------
    flat = rcpgenerator.Packing(N=100, Ndim=2, box=[1.0] * 2, walls=[0] * 2, seed=3,
                                rod_length=2.0, rod_radius=0.5, rod_capped=False)
    flat.pack()
    check_arrays(flat, 100, 2, shaped=True)
    d = np.asarray(flat.diameters) * flat.shrink_factor
    half, radius = 0.5 * d * 2.0, d * 0.5            # half the straight length, and half the width
    x, axis = np.asarray(flat.positions), np.asarray(flat.orientations)[:, :, 0]
    assert 2.0 * np.hypot(half, radius).max() < 0.5, "a rod is longer than half the box"
    gap = rectangle_gap(x, axis, half, radius, box2)
    assert gap > -1e-9, "flat-ended rods overlap at the contact state: %g" % gap
    assert flat.shrink_factor == 1.0 or rectangle_gap(x, axis, half * (1 + TIGHT), radius * (1 + TIGHT), box2) < 0.0, "flat-ended shrink is not tight"
    volume = float(np.sum(4.0 * half * radius))
    assert abs(flat.phi_contact / volume - 1.0) < 1e-12, (flat.phi_contact, volume)
    assert 0.80 < flat.phi_contact < 0.93, flat.phi_contact
    check_round_trip(flat, shaped=True)

    # ---- rods with round ends, 3D (spherocylinders) ----------------------------------------------------------
    capped = rcpgenerator.Packing(N=150, Ndim=3, box=[1.0] * 3, walls=[0] * 3, seed=4,
                                  rod_length=2.0, rod_radius=0.5, rod_capped=True)
    capped.pack()
    check_arrays(capped, 150, 3, shaped=True)
    d = np.asarray(capped.diameters) * capped.shrink_factor
    half, radius = 0.5 * d * 2.0, d * 0.5            # half the straight length; the caps add one radius at each end
    x, axis = np.asarray(capped.positions), np.asarray(capped.orientations)[:, :, 0]
    assert 2.0 * (half + radius).max() < 0.5, "a rod is longer than half the box"
    gap = capped_rod_gap(x, axis, half, radius, box3)
    assert gap > -1e-9, "capped rods overlap at the contact state: %g" % gap
    assert capped.shrink_factor == 1.0 or capped_rod_gap(x, axis, half * (1 + TIGHT), radius * (1 + TIGHT), box3) < 0.0, "capped shrink is not tight"
    volume = float(np.sum(math.pi * radius ** 2 * 2.0 * half + unit_ball_volume(3) * radius ** 3))
    assert abs(capped.phi_contact / volume - 1.0) < 1e-12, (capped.phi_contact, volume)
    assert 0.55 < capped.phi_contact < 0.72, capped.phi_contact
    check_round_trip(capped, shaped=True)

    # ---- input validation: a wrong shape input raises, it is not silently accepted ---------------------------
    for bad in (dict(aspect_ratios=[0.5]),                                  # 3D needs two ratios
                dict(aspect_ratios=[0.5, 1.5]),                             # ratios are at most 1 and decreasing
                dict(rod_length=2.0),                                       # rod_length without rod_radius
                dict(rod_length=2.0, rod_radius=0.5, aspect_ratios=[0.5, 0.5])):   # rods and ellipsoids together
        try:
            rcpgenerator.Packing(N=20, Ndim=3, box=[1.0] * 3, walls=[0] * 3, seed=1, **bad)
        except (ValueError, TypeError, NotImplementedError):
            continue
        raise AssertionError("accepted without an error: %r" % (bad,))

    print("wheel_smoke_shapes: OK   phi at contact: sphere %.4f, last ellipsoid %.4f, flat-ended rod (2D) %.4f, capped rod %.4f"
          % (state["phi"], ell.phi_contact, flat.phi_contact, capped.phi_contact))


if __name__ == "__main__":
    main()
