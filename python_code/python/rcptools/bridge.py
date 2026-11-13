#!/usr/bin/env python3
"""npz -> viz bundle bridge. A bundle = a folder with manifest.json + raw float32
arrays (positions, diameters) + optional metadata layers (computed via the analysis
toolbox). The webapp reads any bundle (drag-drop folder/.zip or CLI). Modular: a new
layer = a new manifest entry + a renderer module; no core change.

Usage:
  python bridge.py <packing.npz> <out_bundle_dir> [--layers]
  python bridge.py --ground-truth <out_bundle_dir>
"""
import sys, os, json, pathlib
import numpy as np

# make the sibling `packing_analysis` package importable when bridge is run
# standalone; in Colab (cwd = python_parallel3) it is already on sys.path.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

def _load_npz(npz):
    """Load pos/dia/box from a packing .npz (keys positions/diameters/box, with
    pos/dia fallbacks). No heavy deps -> geometry-only bundles need nothing."""
    z = np.load(npz, allow_pickle=True)
    pos = np.asarray(z["positions"] if "positions" in z else z["pos"], float)
    dia = np.asarray(z["diameters"] if "diameters" in z else z["dia"], float)
    box = np.asarray(z["box"], float)
    return pos, dia, box

def _f32(arr, path):
    np.asarray(arr, dtype="<f4").tofile(path)

def _container_frac(walls, d):
    """Box fraction filled by a curved (hyperspherical) container (walls[0] = -t); 1.0 if flat."""
    import math
    t = -int(walls[0]) if (walls is not None and len(walls) and int(walls[0]) < 0) else 0
    if t < 2:
        return 1.0
    return math.pi ** (t / 2.0) / math.gamma(t / 2.0 + 1.0) * (0.5 ** t)


def write_bundle(out_dir, pos, dia, box, walls=None, layers=None, aspect_ratios=None, orientations=None,
                 rod_length=None, rod_radius=None, rod_capped=True):
    """Write a viewer bundle. Ellipsoids: also pass `aspect_ratios` (Ndim-1 numbers) and `orientations` (one
    Ndim x Ndim rotation matrix per particle, as `Packing.orientations`); both or neither. Without them the
    bundle is exactly the sphere bundle it always was.

    Rods (2D and 3D): pass `rod_length`, `rod_radius` (each one number or a list of N, in units of the
    particle's D as in `Packing`), `rod_capped` and `orientations`; rod inputs and orientations both or neither,
    and never together with `aspect_ratios`. What a rod bundle adds to the sphere bundle:
      manifest["rods"] = {"capped": true | false,
                          "length": <number>  or  {"file": "rodl.f32", "dtype": "float32", "shape": [N]},
                          "radius": <number>  or  {"file": "rodr.f32", "dtype": "float32", "shape": [N]}}
          (a number when one value holds for every particle, a file of N float32 when a list was given;
           straight section = D x length, radius = D x radius, capped = half-sphere ends, else flat ends)
      manifest["orientations"] = the same entry as for ellipsoids (ori.f32, row-major rotation per particle);
          the rod's axis is the first body axis = column 0 of each matrix
      manifest["phi"] = the true volume (area) fraction of the rods.
    Call: write_bundle(dir, p.positions, p.diameters, p.box, walls=p.walls, rod_length=p.rod_length,
                       rod_radius=p.rod_radius, rod_capped=p.rod_capped, orientations=p.orientations)."""
    out = pathlib.Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    pos = np.asarray(pos, float); dia = np.asarray(dia, float); box = list(map(float, box))
    d = len(box); N = len(dia)
    rod = rod_length is not None or rod_radius is not None
    if rod:
        if aspect_ratios:
            raise ValueError("rod bundle: rod_length / rod_radius were given together with aspect_ratios; a "
                             "packing is of rods or of ellipsoids, not both")
        if rod_length is None or rod_radius is None or orientations is None or len(orientations) == 0:
            raise ValueError("rod bundle: pass rod_length, rod_radius and orientations together (all three), or none")
        if d not in (2, 3):
            raise NotImplementedError(f"viewer bundle: rods are drawn in 2D and 3D only, not Ndim = {d}")
        rod_vals = {}
        for key, val in (("length", rod_length), ("radius", rod_radius)):
            arr = np.asarray(val, float)
            if arr.ndim == 0:
                rod_vals[key] = float(arr)
            elif arr.shape == (N,):
                rod_vals[key] = arr
            else:
                raise ValueError(f"rod_{key} must be one number or a list of N = {N} numbers, got shape {arr.shape}")
        rod_ori = np.asarray(orientations, float).reshape(N, d * d)
    ell = (not rod) and (bool(aspect_ratios) or (orientations is not None and len(orientations) > 0))
    if ell:
        if not aspect_ratios or orientations is None or len(orientations) == 0:
            raise ValueError("ellipsoid bundle: pass both aspect_ratios and orientations, or neither")
        if isinstance(aspect_ratios[0], (list, tuple, np.ndarray)):
            raise NotImplementedError("viewer bundle: a per-particle aspect_ratios table is not supported yet; the "
                                      "viewer draws one shape for all particles (Packing.render draws tables)")
        ratios = [float(a) for a in aspect_ratios]
        ori = np.asarray(orientations, float).reshape(N, d * d)
        if len(ratios) != d - 1:
            raise ValueError(f"aspect_ratios needs Ndim-1 = {d - 1} entries, got {len(ratios)}")
    _f32(pos.ravel(), out / "pos.f32"); _f32(dia, out / "dia.f32")
    man = {"schema": 1, "ndim": d, "box": box, "N": N,
           # boundary conditions per dimension: 0 = periodic, 1 = hard wall; a leading
           # negative -t = one hyperspherical (disk/cylinder/sphere) hard boundary over
           # the first t dims. None -> viewer falls back to a dashed periodic box
           # (backward compatible with older bundles that carry no walls).
           "walls": ([int(w) for w in walls] if walls is not None else None),
           "positions": {"file": "pos.f32", "dtype": "float32", "shape": [N, d]},
           "diameters": {"file": "dia.f32", "dtype": "float32", "shape": [N]},
           "phi": float((np.pi/6 if d == 3 else np.pi/4) * (dia**d).sum() / np.prod(box) / _container_frac(walls, d)),
           "layers": []}
    if ell:   # shape and orientation; phi is the true particle volume fraction
        _f32(ori.ravel(), out / "ori.f32")
        man["aspect_ratios"] = ratios
        man["orientations"] = {"file": "ori.f32", "dtype": "float32", "shape": [N, d * d],
                               "layout": "row-major rotation matrix per particle, columns = body axes"}
        man["phi"] = float(man["phi"] * np.prod(ratios))
    if rod:   # shape (per packing or per particle) and orientation; phi is the true particle volume fraction
        _f32(rod_ori.ravel(), out / "ori.f32")
        man["rods"] = {"capped": bool(rod_capped)}
        for key, fname in (("length", "rodl.f32"), ("radius", "rodr.f32")):
            if isinstance(rod_vals[key], float):
                man["rods"][key] = rod_vals[key]
            else:
                _f32(rod_vals[key], out / fname)
                man["rods"][key] = {"file": fname, "dtype": "float32", "shape": [N]}
        man["orientations"] = {"file": "ori.f32", "dtype": "float32", "shape": [N, d * d],
                               "layout": "row-major rotation matrix per particle, columns = body axes; the rod's axis is column 0"}
        L_, R_ = dia * rod_vals["length"], dia * rod_vals["radius"]      # actual straight length and radius
        if d == 2:
            vol = 2.0 * R_ * L_ + (np.pi * R_ ** 2 if rod_capped else 0.0)
        else:
            vol = np.pi * R_ ** 2 * L_ + (4.0 / 3.0 * np.pi * R_ ** 3 if rod_capped else 0.0)
        man["phi"] = float(np.sum(vol) / np.prod(box) / _container_frac(walls, d))
    man["analysis_summary"] = {}
    for lyr in (layers or []):
        man["layers"].append(lyr["meta"])
        if lyr.get("summary"):
            man["analysis_summary"][lyr["meta"]["name"]] = lyr["summary"]
        if "data" in lyr:
            _f32(lyr["data"], out / lyr["meta"]["file"])
        elif "json" in lyr:
            (out / lyr["meta"]["file"]).write_text(json.dumps(lyr["json"]))
    (out / "manifest.json").write_text(json.dumps(man, indent=2))
    print(f"bundle -> {out}  (N={N}, ndim={d}, phi={man['phi']:.4f}, layers={[l['name'] for l in man['layers']]})")
    return man

def from_npz(npz, out_dir, with_layers=False, analyses=None):
    """Build a bundle from a packing .npz. `analyses` = list of registry names
    (e.g. ["local_phi","pores"]); if None and with_layers, use all registered."""
    pos, dia, box = _load_npz(npz)
    layers = build_layers(pos, dia, box, analyses) if (with_layers or analyses) else []
    return write_bundle(out_dir, pos, dia, box, layers=layers)

def build_layers(pos, dia, box, analyses=None):
    """Compute analysis layers via the packing_analysis registry. `analyses` =
    list of names; None -> all registered analyses. Each missing dep / failure
    skips-with-warning (the registry handles it), so a bundle never breaks."""
    from rcptools import analysis as pa
    names = pa.available_analyses() if analyses is None else analyses
    return pa.run_analyses(names, pos, dia, box)

def build_pore_layer(pos, dia, box, dcfrac=0.35, frac=0.35):
    """Back-compat shim — the pore-map analysis now lives in packing_analysis.registry."""
    from rcptools import analysis as pa
    return pa.registry.pores(pos, dia, box, dcfrac=dcfrac, frac=frac)

def ground_truth(out_dir):
    """Tiny known packing for pixel-exact rendering tests. Box 10x10, 4 circles at
    known centers with known radii, non-overlapping."""
    pos = np.array([[2.5, 2.5], [7.5, 2.5], [2.5, 7.5], [7.5, 7.5]], float)
    dia = np.array([3.0, 2.0, 2.0, 1.0], float)
    box = [10.0, 10.0]
    return write_bundle(out_dir, pos, dia, box)

if __name__ == "__main__":
    a = sys.argv[1:]
    if a and a[0] == "--ground-truth":
        ground_truth(a[1])
    else:
        from_npz(a[0], a[1], with_layers=("--layers" in a))
