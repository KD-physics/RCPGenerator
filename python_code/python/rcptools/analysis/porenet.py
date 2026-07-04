#!/usr/bin/env python3
"""Per-pore network tool (2D & 3D) for a coarse scaffold of a jammed packing.

Backbone: periodic Delaunay of the scaffold (D>Dc) centers. Each simplex is a pore
ELEMENT; its void = simplex measure - the scaffold-particle caps inside it, where
each particle's measure is split among its incident simplices by the subtended
angle (2D) / solid angle (3D). That accounting guarantees
    sum(simplex void) + sum(scaffold particle measure) = box   (EXACT closure),
because the angle fractions around any interior vertex sum to 1. Fines (D<=Dc) are
assigned to the simplex containing them -> per-element fill eta, rho.

(This is the simplex-element layer with exact closure + throats-as-faces. Merging
elements into physical pore bodies and assigning throat volumes is layered on top.)

Dimension-agnostic (Ndim in {2,3}); grid-free. CLI renders 2D.
"""
import sys, itertools, math, pathlib, numpy as np
from scipy.spatial import Delaunay, cKDTree
import voidfill_method as vf

# ----------------------------------------------------------------- geometry
def simplex_measure(V):
    """Area (2D) / volume (3D) of simplices. V: (m, d+1, d)."""
    d = V.shape[2]
    M = V[:, 1:, :] - V[:, :1, :]            # (m, d, d)
    return np.abs(np.linalg.det(M)) / math.factorial(d)

def vertex_fraction(V):
    """Per-vertex angle fraction (of full 2pi / 4pi) at each simplex vertex.
    Returns (m, d+1). Sum over simplices incident to an interior vertex -> 1."""
    m, k, d = V.shape
    frac = np.zeros((m, k))
    for a in range(k):
        others = [i for i in range(k) if i != a]
        vecs = V[:, others, :] - V[:, a:a+1, :]      # (m, d, d) edge vectors from a
        if d == 2:
            u = vecs[:, 0, :]; w = vecs[:, 1, :]
            cu = u / np.linalg.norm(u, axis=1, keepdims=True)
            cw = w / np.linalg.norm(w, axis=1, keepdims=True)
            ang = np.arccos(np.clip((cu * cw).sum(1), -1, 1))
            frac[:, a] = ang / (2 * np.pi)
        else:
            A = vecs[:, 0, :]; B = vecs[:, 1, :]; C = vecs[:, 2, :]
            na = np.linalg.norm(A, axis=1); nb = np.linalg.norm(B, axis=1); nc = np.linalg.norm(C, axis=1)
            triple = np.abs(np.einsum("ij,ij->i", A, np.cross(B, C)))
            den = (na * nb * nc + np.einsum("ij,ij->i", A, B) * nc
                   + np.einsum("ij,ij->i", A, C) * nb + np.einsum("ij,ij->i", B, C) * na)
            Omega = 2 * np.arctan2(triple, den)          # solid angle in [0, 2pi)
            frac[:, a] = Omega / (4 * np.pi)
    return frac

def ball_measure(r, d):
    return (np.pi * r**2) if d == 2 else (4.0 / 3.0 * np.pi * r**3)

# ----------------------------------------------------------------- periodic Delaunay
def periodic_delaunay(pts, box):
    """Return (P, orig, simplices) where simplices index into the replicated P and
    each distinct periodic simplex is kept exactly once (primary-centroid rule)."""
    d = pts.shape[1]; L = np.asarray(box, float)
    pts = pts % L
    reps, origs = [], []
    for o in itertools.product((-1, 0, 1), repeat=d):
        reps.append(pts + np.array(o) * L); origs.append(np.arange(len(pts)))
    P = np.vstack(reps); orig = np.concatenate(origs)
    tri = Delaunay(P)
    cent = P[tri.simplices].mean(axis=1)
    inbox = np.all((cent >= 0) & (cent < L), axis=1)
    return P, orig, tri, np.where(inbox)[0]      # tri + indices of kept (primary) simplices

# ----------------------------------------------------------------- pore network
def pore_network(pos, dia, box, Dc, phi_fill=None):
    pos = np.asarray(pos, float); dia = np.asarray(dia, float); box = np.asarray(box, float)
    d = len(box); Vbox = float(np.prod(box))
    coarse = dia > Dc; fine = ~coarse
    cpos = pos[coarse] % box; crad = dia[coarse] / 2.0
    P, orig, tri, keep = periodic_delaunay(cpos, box)
    simp = tri.simplices[keep]
    Prad = crad[orig]
    V = P[simp]                                          # (m, d+1, d)
    meas = simplex_measure(V)
    frac = vertex_fraction(V)
    caps = (frac * ball_measure(Prad[simp], d)).sum(axis=1)   # scaffold measure inside each simplex
    void = meas - caps
    # assign fines to the simplex containing them (Delaunay.find_simplex on replicated tri is
    # awkward; use centroid-nearest as a robust proxy, then exact point-in-simplex refine)
    cent = V.mean(axis=1)
    info = dict(P=P, orig=orig, simp=simp, V=V, meas=meas, void=void, cent=cent,
                tri=tri, keep=keep, pore_bigs=orig[simp], crad=crad,   # bounding bigs + scaffold radii
                Dc=Dc, d=d, Vbox=Vbox, coarse_meas=float(ball_measure(crad, d).sum()),
                box=box, fine_pos=pos[fine] % box, fine_meas=ball_measure(dia[fine] / 2.0, d))
    # closure check (scaffold-only): sum(void) + sum(coarse particle measure) == box
    info["closure_resid"] = abs(void.sum() + info["coarse_meas"] - Vbox) / Vbox
    info["n_simplices"] = len(simp); info["n_coarse"] = int(coarse.sum()); info["n_fine"] = int(fine.sum())
    return info

def assign_fines(info):
    """Assign each fine to the kept simplex containing it, via periodic min-image
    barycentric test against the nearest simplex centroids. Robust at the periodic
    boundary (where the kept primary-centroid simplices poke outside [0,L))."""
    V = info["V"]; cent = info["cent"]; box = info["box"]; d = info["d"]
    fp = info["fine_pos"]; m = len(cent)
    tree = cKDTree(cent % box, boxsize=box)
    Kq = min(16, m)
    _, nn = tree.query(fp % box, k=Kq)
    nn = nn.reshape(len(fp), -1)
    V0 = V[:, 0, :]
    T = np.transpose(V[:, 1:, :] - V[:, :1, :], (0, 2, 1))    # (m,d,d)
    Tinv = np.linalg.inv(T)
    Sfine = np.zeros(m); assigned = np.zeros(len(fp), bool)
    for col in range(Kq):
        j = nn[:, col]
        p = fp - box * np.round((fp - cent[j]) / box)        # min-image toward simplex j
        rel = p - V0[j]
        lam1 = np.einsum("nij,nj->ni", Tinv[j], rel)
        lam0 = 1.0 - lam1.sum(1)
        inside = (lam0 >= -1e-9) & np.all(lam1 >= -1e-9, axis=1)
        take = inside & (~assigned)
        np.add.at(Sfine, j[take], info["fine_meas"][take])
        assigned[take] = True
    info["fine_solid"] = Sfine
    info["assigned_fine_meas"] = float(Sfine.sum())
    info["total_fine_meas"] = float(info["fine_meas"].sum())
    info["frac_assigned"] = float(assigned.mean())
    info["eta"] = np.divide(Sfine, info["void"], out=np.zeros_like(Sfine), where=info["void"] > 0)
    return info

def throats(info):
    """Throats = faces shared by two adjacent pore-elements. Returns an array of
    (pore_i, pore_j, face_area, aperture). face_area is edge length (2D) / triangle
    area (3D); aperture = min pairwise surface gap among the face's scaffold bigs
    (the constriction). Periodic adjacency via the replicated triangulation."""
    tri = info["tri"]; keep = info["keep"]; orig = info["orig"]; P = info["P"]
    crad = info["crad"]; d = info["d"]; S = tri.simplices
    kept_by_key = {tuple(sorted(orig[S[ti]])): row for row, ti in enumerate(keep)}
    seen = set(); out = []
    for row, ti in enumerate(keep):
        vi = S[ti]; oi = set(orig[vi])
        for nb in tri.neighbors[ti]:
            if nb < 0:
                continue
            r2 = kept_by_key.get(tuple(sorted(orig[S[nb]])), -1)
            if r2 < 0 or r2 == row:
                continue
            a, b = (row, r2) if row < r2 else (r2, row)
            if (a, b) in seen:
                continue
            seen.add((a, b))
            common = oi & set(orig[S[nb]])
            shared = [(orig[v], P[v]) for v in vi if orig[v] in common]
            fc = np.array([s[1] for s in shared]); fr = crad[np.array([s[0] for s in shared])]
            area = (np.linalg.norm(fc[1] - fc[0]) if d == 2
                    else 0.5 * np.linalg.norm(np.cross(fc[1] - fc[0], fc[2] - fc[0])))
            ap = np.inf; ref = 0.0                          # min pair gap + that pair's (r_i+r_j)
            for p in range(len(fc)):
                for q in range(p + 1, len(fc)):
                    g = np.linalg.norm(fc[p] - fc[q]) - (fr[p] + fr[q])
                    if g < ap:
                        ap = g; ref = fr[p] + fr[q]
            out.append((a, b, float(area), float(ap), float(ref)))
    return np.array(out) if out else np.empty((0, 5))

def merge_pores_gap(info, th, frac=0.5):
    """Merge two pore-elements into one pore body iff the throat between them is a
    WIDE passage: min Delaunay-pair gap > frac * (r_i + r_j) of the defining pair.
    Narrow throats (gap <= frac*(r_i+r_j)) stay as pore boundaries, so disconnected
    pores remain disconnected. Pure geometric criterion (Ken's): the throat gap
    comes straight from the Delaunay pair. Union-find over simplices; every simplex
    lands in exactly one body -> body_void + particles = box (exact)."""
    m = len(info["void"]); parent = list(range(m))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for a, b, area, gap, ref in th:
        if gap > frac * ref:                                 # wide -> same pore body
            parent[find(int(a))] = find(int(b))
    bodies = {}
    for i in range(m):
        bodies.setdefault(find(i), []).append(i)
    throats = [(int(a), int(b), area, gap, ref) for a, b, area, gap, ref in th
               if not (gap > frac * ref)]                    # narrow faces = real throats
    info["merge"] = dict(bodies=bodies, n_bodies=len(bodies),
                         is_throat=np.zeros(m, bool), n_throat_simplices=0,
                         n_throat_faces=len(throats), throats=throats, frac=frac,
                         body_void=float(info["void"].sum()), throat_void=0.0)
    return info

def circumcenters(V):
    """Circumcenter of each simplex. V: (m, d+1, d) -> (m, d)."""
    A = 2.0 * (V[:, 1:, :] - V[:, :1, :])
    b = (V[:, 1:, :]**2).sum(2) - (V[:, :1, :]**2).sum(2)
    return np.linalg.solve(A, b[..., None])[..., 0]

def merge_pores(info, th, h=0.10, apf=0.4):
    """Watershed on simplex void-depth: deep simplices grow into pore BODIES;
    saddle simplices (adjacent to >=2 distinct basins) become THROATS with volume.
    Each simplex -> exactly one body or one throat, so the void split is exact:
        body_void + throat_void = total void  =>  +particles = box  (3-way closure).
    depth_i = circumradius - largest scaffold radius (the inscribed empty-sphere
    radius at the simplex's Voronoi vertex)."""
    V = info["V"]; void = info["void"]; m = len(void)
    cc = circumcenters(V); R = np.linalg.norm(cc - V[:, 0, :], axis=1)
    bigr = info["crad"][info["pore_bigs"]]                    # (m, d+1)
    depth = R - bigr.max(axis=1)
    adj = [[] for _ in range(m)]
    for row in th:
        a, b = int(row[0]), int(row[1]); adj[a].append(b); adj[b].append(a)
    label = np.full(m, -1, int); is_th = np.zeros(m, bool); nxt = 0
    for i in np.argsort(-depth):                              # flood from deepest
        labs = {label[j] for j in adj[i] if label[j] >= 0 and not is_th[j]}
        if not labs:
            label[i] = nxt; nxt += 1
        elif len(labs) == 1:
            label[i] = labs.pop()
        else:
            is_th[i] = True                                  # provisional saddle -> throat
    # --- HYBRID merge: absorb a saddle (merge its basins) only if BOTH
    #   (i) shallow barrier  : depth[t] >= (1-h)*min(peak)   [persistence; fixes
    #       one physical pore split by a minor internal max], AND
    #   (ii) wide enough     : depth[t] >= apf*rbig[t]        [aperture floor; never
    #       merge across a NARROW constriction -> disconnected pores stay separate].
    rbig = bigr.max(axis=1)
    parent = list(range(nxt)); peak = np.zeros(nxt)
    for i in range(m):
        if not is_th[i]:
            peak[label[i]] = max(peak[label[i]], depth[i])
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for t in sorted((i for i in range(m) if is_th[i]), key=lambda i: -depth[i]):
        if depth[t] < apf * rbig[t]:                         # narrow throat -> keep boundary
            continue
        bas = {find(label[j]) for j in adj[t] if label[j] >= 0 and not is_th[j]}
        for A in list(bas):
            for B in list(bas):
                A, B = find(A), find(B)
                if A != B and depth[t] >= (1 - h) * min(peak[A], peak[B]):
                    parent[A] = B; peak[B] = max(peak[A], peak[B])
        nb = [find(label[j]) for j in adj[t] if label[j] >= 0 and not is_th[j]]
        if len(set(nb)) <= 1:
            is_th[t] = False; label[t] = nb[0] if nb else label[t]
    bodies = {}
    for i in range(m):
        if not is_th[i]:
            bodies.setdefault(find(label[i]), []).append(i)
    body_void = float(sum(void[ix].sum() for ix in bodies.values()))
    throat_void = float(void[is_th].sum())
    info["merge"] = dict(label=label, is_throat=is_th, depth=depth, bodies=bodies,
                         n_bodies=len(bodies), n_throat_simplices=int(is_th.sum()),
                         body_void=body_void, throat_void=throat_void)
    return info

# ----------------------------------------------------------------- 2D render
def render_2d(info, pos, dia, Dc, out, nlabel=8):
    """Render the MERGED pore network: pore bodies colored by volume, throat
    simplices in grey, scaffold in red, fines as dots, a few pore volumes labeled."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.patches import Circle
    box = info["box"]; V = info["V"]; mg = info.get("merge")
    fig, ax = plt.subplots(figsize=(9, 9))
    if mg is not None:
        isth = mg["is_throat"]; bodyvol = np.zeros(len(V)); bodycent = {}
        # categorical color per pore body (random shuffle so neighbors differ)
        labs = list(mg["bodies"].keys())
        rng = np.random.default_rng(1)
        colmap = {lab: np.r_[0.25 + 0.7 * rng.random(3), 1.0] for lab in labs}  # distinct RGB
        face = np.zeros((len(V), 4))
        for lab, ix in mg["bodies"].items():
            bodyvol[ix] = info["void"][ix].sum()
            bodycent[lab] = info["cent"][ix].mean(axis=0)
            for i in ix:
                face[i] = colmap[lab]
        ax.add_collection(PolyCollection(V[~isth], facecolors=face[~isth],
                                         edgecolors="0.6", linewidths=0.2))
        ax.add_collection(PolyCollection(V[isth], facecolors="0.82", edgecolors="0.55",
                                         linewidths=0.2, alpha=0.85))
        big = sorted(mg["bodies"], key=lambda L: -bodyvol[mg["bodies"][L][0]])[:nlabel]
        for lab in big:
            c = bodycent[lab]
            ax.text(c[0], c[1], f"{bodyvol[mg['bodies'][lab][0]]:.3f}", fontsize=6,
                    ha="center", color="black")
        pc = None; cbl = None
        ttl = (f"Merged pores (distinct color per body) + throats (grey); "
               f"{mg['n_bodies']} bodies, {mg['n_throat_simplices']} throats; 3-way closure exact")
    else:
        pc = PolyCollection(V, array=info["eta"], cmap="viridis", edgecolors="0.5", linewidths=0.3)
        ax.add_collection(pc); cbl = "eta"; ttl = "Pore network (simplex layer)"
    fine = dia <= Dc; coarse = dia > Dc
    ax.scatter((pos[fine] % box)[:, 0], (pos[fine] % box)[:, 1], s=1.2, c="white", alpha=0.45)
    for x, r in zip(pos[coarse] % box, dia[coarse] / 2):
        ax.add_patch(Circle(x, r, fill=False, ec="red", lw=0.8))
    ax.set_xlim(0, box[0]); ax.set_ylim(0, box[1]); ax.set_aspect("equal")
    ax.set_title(ttl, fontsize=9)
    if pc is not None:
        plt.colorbar(pc, ax=ax, shrink=0.7, label=cbl)
    plt.tight_layout(); plt.savefig(out, dpi=140); plt.close()

def main():
    npz = sys.argv[1]; frac = float(sys.argv[2]) if len(sys.argv) > 2 else 0.30
    gfrac = float(sys.argv[3]) if len(sys.argv) > 3 else 0.10   # gap-merge fraction
    pos, dia, box = vf._load(npz); Dc = frac * dia.max()
    info = assign_fines(pore_network(pos, dia, box, Dc))
    print(f"# {pathlib.Path(npz).name}  Ndim={len(box)} N={len(dia)} Dc/Dmax={frac}")
    print(f"  n_coarse={info['n_coarse']} n_fine={info['n_fine']} n_simplices={info['n_simplices']}")
    print(f"  CLOSURE residual = {info['closure_resid']:.2e}   (pores+scaffold == box)")
    print(f"  fines assigned = {info['frac_assigned']:.4f} (count), "
          f"{info['assigned_fine_meas']/info['total_fine_meas']:.4f} (measure)")
    print(f"  void fraction = {info['void'].sum()/info['Vbox']:.4f}  "
          f"mean eta = {np.average(info['eta'], weights=info['void']):.4f}")
    th = throats(info)
    print(f"  throats = {len(th)} (faces between pores); "
          f"total face {'length' if len(box)==2 else 'area'} = {th[:,2].sum():.4f}; "
          f"aperture median = {np.median(th[:,3]):.4f}" if len(th) else "  throats = 0")
    print(f"  each pore-element bounded by {info['d']+1} scaffold bigs")
    if len(th):
        info = merge_pores_gap(info, th, frac=gfrac); mg = info["merge"]
        clo = (mg["body_void"] + info["coarse_meas"]) / info["Vbox"]
        print(f"  GAP-MERGE (frac={gfrac}): {mg['n_bodies']} pore bodies; "
              f"{mg['n_throat_faces']} throat faces (narrow Delaunay gaps)")
        print(f"  CLOSURE (pore bodies + particles)/box = {clo:.12f}  resid={abs(clo-1):.2e}")
    if len(box) == 2:
        FIG = pathlib.Path(__file__).resolve().parent.parent / "Figures"
        out = FIG / f"porenet_{pathlib.Path(npz).stem}.png"
        render_2d(info, pos, dia, Dc, out)
        print(f"  figure -> {out.name}")

if __name__ == "__main__":
    main()
