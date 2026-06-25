"""Analysis registry: name -> function(pos, dia, box) -> a bundle "layer" dict.

A layer dict is:
  {"meta": {"name","type","file"}, "data": ndarray}      # scalar_per_particle
  {"meta": {"name","type","file"}, "json": {...}}         # polygons
and may carry an optional "summary": {...} of scalar metrics (saved into the
bundle manifest for the science).

Each analysis lazy-imports its own dependencies so importing this module is cheap
and dependency-free. EXPANDABLE: register a new analysis by adding it to ANALYSES;
nothing downstream changes.
"""
import numpy as np


def local_phi(pos, dia, box):
    """Per-particle local packing fraction from a radical (power) Voronoi
    tessellation. Needs pyvoro."""
    import voidfill_method as vf
    box = np.asarray(box, float)
    cell = vf.radical_voronoi_volumes(pos % box, dia / 2.0, box)
    d = len(box)
    pv = (np.pi / 6 if d == 3 else np.pi / 4) * dia ** d
    lp = np.divide(pv, cell, out=np.zeros_like(pv), where=cell > 0)
    mean = float(np.average(lp, weights=cell)) if cell.sum() > 0 else None
    return {"meta": {"name": "local_phi", "type": "scalar_per_particle", "file": "local_phi.f32"},
            "data": lp, "summary": {"mean_local_phi": mean}}


def pores(pos, dia, box, dcfrac=0.35, frac=0.35):
    """Segmented pore map (2D): Delaunay pore network, gap-merged into bodies,
    per-cell fill fraction eta. Needs scipy. 2D-only."""
    box = np.asarray(box, float)
    if len(box) != 2:
        raise ValueError("pores layer is 2D-only")
    import porenet as pn
    dia = np.asarray(dia, float)
    Dc = dcfrac * dia.max()
    info = pn.assign_fines(pn.pore_network(pos, dia, box, Dc))
    th = pn.throats(info)
    info = pn.merge_pores_gap(info, th, frac=frac)
    V = info["V"]; void = info["void"]; fine = info["fine_solid"]
    bodies = info["merge"]["bodies"]
    body_of, eta_of = {}, {}
    for bid, ix in bodies.items():
        w = float(void[ix].sum()); e = float(fine[ix].sum() / w) if w > 0 else 0.0
        for i in ix:
            body_of[i] = int(bid); eta_of[i] = e
    tris, bid_l, eta_l = [], [], []
    for i in range(len(V)):
        t = V[i][:, :2]
        tris.append([round(float(c), 5) for p in t for c in p])
        bid_l.append(body_of.get(i, -1)); eta_l.append(round(eta_of.get(i, 0.0), 4))
    return {"meta": {"name": "pores", "type": "polygons", "file": "pores.json"},
            "json": {"tris": tris, "body": bid_l, "eta": eta_l, "n_bodies": len(bodies)},
            "summary": {"n_bodies": len(bodies), "n_cells": len(V)}}


# name -> analysis function. Add new analyses here (e.g. "chord", "depletion_shell").
ANALYSES = {
    "local_phi": local_phi,
    "pores": pores,
}


def available():
    """Names of registered analyses."""
    return sorted(ANALYSES)


def run(names, pos, dia, box):
    """Run the requested analyses; skip-with-warning on unknown name / missing
    dependency / failure (never raises). Returns a list of layer dicts."""
    pos = np.asarray(pos, float); dia = np.asarray(dia, float)
    out = []
    for nm in (names or []):
        fn = ANALYSES.get(nm)
        if fn is None:
            print(f"  [analysis '{nm}' not in registry {available()}; skipped]")
            continue
        try:
            out.append(fn(pos, dia, box))
        except Exception as e:
            print(f"  [analysis '{nm}' skipped: {type(e).__name__}: {str(e)[:90]}]")
    return out
