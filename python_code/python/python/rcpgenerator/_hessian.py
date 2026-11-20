"""The Hessian of the contact energy of a packing, and its eigenvalues: a check that a small system is jammed.

The matrix is built from force evaluations of the packing's own force law at the state as it is (sizes fixed; ``mu``
plays no part): every coordinate of a particle is displaced by a small amount in both directions and the change of
the forces on that particle and on its neighbours is recorded. Particles that are neither neighbours nor share a
neighbour are displaced together, so the number of force evaluations does not grow with the number of particles.
It works for every shape without second derivatives of the contact function.

Coordinates, per particle: its Ndim position components; for ellipsoids and rods also Ndim (Ndim - 1) / 2 rotation
components, each a rotation in the plane of two axes, measured as the arc length on the particle's bounding sphere
(angle x bounding radius), so that every entry of the matrix has the units of a stiffness (pair stiffness 1)."""
from __future__ import annotations

import numpy as np

from ._rcpgenerator import contact_stress as _contact_stress
from ._rcpgenerator import evaluate_once as _evaluate_once

HESSIAN_DOF_CAP = 15000


class HessianTooLarge(RuntimeError):
    pass


def check_size(n_dof: int, force: bool, what: str) -> None:
    if n_dof > HESSIAN_DOF_CAP and not force:
        gb = 8.0 * n_dof * n_dof / 1e9
        raise HessianTooLarge(
            f"{what}: this packing has {n_dof:,} degrees of freedom, above the limit of {HESSIAN_DOF_CAP:,}. All "
            f"eigenvalues of a system this size need a dense matrix of about {gb:.1f} GB (8 x n^2 bytes), and the "
            "time grows as the cube of the number of degrees of freedom. Pass force=True to go ahead.")


class Hessian:
    """Sparse Hessian in coordinate form: entry k is ``values[k]`` at (``rows[k]``, ``cols[k]``); ``n`` rows and
    columns; ``dof_per_particle`` coordinates per particle (coordinate c of particle i is index i x dof_per_particle
    + c). ``contacts``: per particle, the number of particles it exchanges force with (from the matrix itself).
    ``to_dense()`` gives the symmetric dense matrix."""

    def __init__(self, n, dpp, rows, cols, values, n_evaluations, step):
        self.n = int(n); self.dof_per_particle = int(dpp)
        self.rows = rows; self.cols = cols; self.values = values
        self.n_evaluations = int(n_evaluations); self.step = float(step)

    @property
    def shape(self):
        return (self.n, self.n)

    @property
    def nnz(self):
        return int(len(self.values))

    def to_dense(self, particles=None) -> np.ndarray:
        """The dense matrix, made exactly symmetric (the mean of the matrix and its transpose). ``particles``: only
        the coordinates of these particles (a sub-block)."""
        H = np.zeros((self.n, self.n))
        np.add.at(H, (self.rows, self.cols), self.values)
        H = 0.5 * (H + H.T)
        if particles is not None:
            d = self.dof_per_particle
            idx = (np.asarray(particles, dtype=np.int64)[:, None] * d + np.arange(d)[None, :]).ravel()
            H = H[np.ix_(idx, idx)]
        return H


def _rotate(R, n, a, b, angle):
    """exp(angle G) R with G[a, b] = +1, G[b, a] = -1: the rotation whose conjugate force is torque component (a, b)."""
    c, s = np.cos(angle), np.sin(angle)
    ra, rb = R[a].copy(), R[b].copy()
    out = R.copy()
    out[a] = c * ra + s * rb
    out[b] = -s * ra + c * rb
    return out


def build(inp: dict, cfg: dict, *, shape: bool, mean_overlap: float | None, step: float | None = None) -> Hessian:
    x0 = np.asarray(inp["positions"], dtype=float)
    N, n = x0.shape
    comps = [(a, b) for a in range(n) for b in range(a + 1, n)] if shape else []
    nrot = len(comps); dpp = n + nrot
    cs = _contact_stress(inp, cfg, True)
    pi = np.asarray(cs["pair_i"], dtype=np.int64); pj = np.asarray(cs["pair_j"], dtype=np.int64)
    A = 0.5 * np.asarray(cs["bounding"], dtype=float)          # bounding radius: the length scale of a rotation
    adj = [[] for _ in range(N)]
    for i, j in zip(pi.tolist(), pj.tolist()):
        adj[i].append(j); adj[j].append(i)
    # colours: two particles of a colour are neither neighbours nor share a neighbour
    colour = np.full(N, -1, dtype=np.int64)
    order = np.argsort([-len(a) for a in adj], kind="stable")
    for i in order.tolist():
        taken = set()
        for j in adj[i]:
            if colour[j] >= 0: taken.add(int(colour[j]))
            for k in adj[j]:
                if colour[k] >= 0: taken.add(int(colour[k]))
        c = 0
        while c in taken: c += 1
        colour[i] = c
    ncol = int(colour.max()) + 1 if N else 0
    Dm = float(np.mean(2.0 * A))
    if step is None:
        # 1e-8 mean diameters: far below the overlaps of a packing (a contact that opens or closes inside the step
        # gives a wrong entry: measured, a step of 1e-6 diameters produced a spurious negative eigenvalue at a mean
        # overlap of 1e-3), far above rounding (the forces are piecewise linear in a displacement for spheres, so a
        # central difference is exact up to rounding).
        step = 1e-8 * Dm
    O0 = np.asarray(inp["orientations"], dtype=float).reshape(N, n, n) if shape else None

    def gradient(x, O):
        d = dict(inp); d["positions"] = x.tolist()
        if shape: d["orientations"] = O.reshape(N, n * n).tolist() if isinstance(inp["orientations"][0], (list, tuple)) else O.ravel().tolist()
        ev = _evaluate_once(d, cfg, 0.0)
        g = np.empty((N, dpp))
        g[:, :n] = -np.asarray(ev["forces"], dtype=float).reshape(N, n)
        if shape: g[:, n:] = -np.asarray(ev["torques"], dtype=float).reshape(N, nrot) / A[:, None]
        return g

    rows, cols, vals = [], [], []
    for c in range(ncol):
        members = np.where(colour == c)[0]
        tgt = np.concatenate([np.array([i] + adj[i], dtype=np.int64) for i in members.tolist()])   # rows that member i touches
        src = np.concatenate([np.full(1 + len(adj[i]), i, dtype=np.int64) for i in members.tolist()])
        for a in range(dpp):
            g = []
            for sgn in (+1.0, -1.0):
                x = x0.copy(); O = O0.copy() if shape else None
                if a < n:
                    x[members, a] += sgn * step
                else:
                    ca, cb = comps[a - n]
                    for i in members.tolist():
                        O[i] = _rotate(O0[i], n, ca, cb, sgn * step / A[i])
                g.append(gradient(x, O))
            dg = (g[0] - g[1]) / (2.0 * step)                    # (N, dpp): d gradient / d coordinate a of the members
            blk = dg[tgt]                                         # (m, dpp)
            rows.append((tgt[:, None] * dpp + np.arange(dpp)[None, :]).ravel())
            cols.append(np.repeat(src * dpp + a, dpp))
            vals.append(blk.ravel())
    rows = np.concatenate(rows) if rows else np.zeros(0, dtype=np.int64)
    cols = np.concatenate(cols) if cols else np.zeros(0, dtype=np.int64)
    vals = np.concatenate(vals) if vals else np.zeros(0)
    keep = vals != 0.0
    H = Hessian(N * dpp, dpp, rows[keep], cols[keep], vals[keep], 2 * ncol * dpp + 1, step)
    H.n_colours = ncol
    H.z = np.asarray(cs["z"], dtype=np.int64)   # contacts of the force law at this state: pairs, plus one per boundary pressed on
    # contacts per particle from the matrix: a neighbour j of i with a non-zero block
    contact = np.zeros(pi.shape[0], dtype=bool)
    if pi.shape[0]:
        key = {}
        blockmax = np.zeros(pi.shape[0])
        for k, (i, j) in enumerate(zip(pi.tolist(), pj.tolist())): key[(i, j)] = k; key[(j, i)] = k
        ri = H.rows // dpp; ci = H.cols // dpp
        off = ri != ci
        for r, cc, v in zip(ri[off].tolist(), ci[off].tolist(), np.abs(H.values[off]).tolist()):
            k = key.get((r, cc))
            if k is not None and v > blockmax[k]: blockmax[k] = v
        contact = blockmax > 1e-9
    H.contact_i = pi[contact]; H.contact_j = pj[contact]
    return H


def spin_modes(n: int, axes) -> int:
    """Rotations that leave a particle unchanged: for semi-axes with equal values, the rotations among them."""
    vals = {}
    for v in axes:
        k = round(float(v), 12); vals[k] = vals.get(k, 0) + 1
    return sum(m * (m - 1) // 2 for m in vals.values())


def verdict(H: Hessian, *, n_dim: int, walls, spin, tol: float = 1e-8) -> dict:
    """All eigenvalues of the Hessian (every particle), and what they say about the particles that are held.

    Held: what is left after taking out the rattlers one after another (a particle with fewer than Ndim + 1
    contacts, walls counted, cannot be held). A rattler may rest against held particles, so it stays in the matrix;
    its own free motions are zero modes that do not count against the packing, and so is a negative eigenvalue
    whose mode lies mostly on rattlers (the quench does not balance a rattler: one squeezed between two neighbours
    sits on a saddle). The zero modes are therefore
    separated: those that move held particles (the rank of the zero eigenvectors on the held coordinates) must be
    exactly the allowed ones: a uniform translation along each periodic direction, a rotation of the whole packing
    in a curved container, and, per held particle, the rotations that leave its shape unchanged (``spin``)."""
    d = H.dof_per_particle; N = H.n // d
    zc = np.asarray(H.z, dtype=np.int64).copy()
    nb = [[] for _ in range(N)]
    for i, j in zip(H.contact_i.tolist(), H.contact_j.tolist()):
        nb[i].append(j); nb[j].append(i)
    held = np.ones(N, dtype=bool)
    stack = [i for i in range(N) if zc[i] < n_dim + 1]
    while stack:
        i = stack.pop()
        if not held[i]: continue
        held[i] = False
        for j in nb[i]:
            if held[j]:
                zc[j] -= 1
                if zc[j] < n_dim + 1: stack.append(j)
    idx = np.where(held)[0]
    ev, V = np.linalg.eigh(H.to_dense()) if H.n else (np.zeros(0), np.zeros((0, 0)))
    zero = np.abs(ev) <= tol
    n_neg = int((ev < -tol).sum()); n_zero = int(zero.sum())
    n_zero_held = 0; n_neg_held = 0
    rows = (idx[:, None] * d + np.arange(d)[None, :]).ravel() if len(idx) else np.zeros(0, dtype=np.int64)
    if len(idx) and n_zero:
        sv = np.linalg.svd(V[np.ix_(rows, np.where(zero)[0])], compute_uv=False)
        n_zero_held = int((sv > 1e-3).sum())
    if len(idx) and n_neg:   # a negative mode counts against the packing when most of it lies on held particles
        wt = (V[np.ix_(rows, np.where(ev < -tol)[0])] ** 2).sum(0)
        n_neg_held = int((wt >= 0.5).sum())
    tc = -int(walls[0]) if int(walls[0]) < 0 else 0
    rigid = sum(1 for k in range(tc, n_dim) if int(walls[k]) == 0) + tc * (tc - 1) // 2
    allowed = (rigid + int(sum(int(spin[i]) for i in idx.tolist()))) if len(idx) else 0
    return dict(eigenvalues=ev, n_negative=n_neg, n_negative_held=n_neg_held, n_negative_rattlers=n_neg - n_neg_held, n_zero=n_zero, n_zero_held=n_zero_held, n_zero_allowed=int(allowed),
                n_zero_rattlers=n_zero - n_zero_held, n_rattlers=int(N - len(idx)), n_held=int(len(idx)), rattlers=np.where(~held)[0].tolist(),
                tol=float(tol), lowest_nonzero=float(ev[ev > tol].min()) if (ev > tol).any() else None,
                jammed=bool(len(idx) > 0 and n_neg_held == 0 and n_zero_held == allowed))
