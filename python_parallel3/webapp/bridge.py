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

def write_bundle(out_dir, pos, dia, box, layers=None):
    out = pathlib.Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    pos = np.asarray(pos, float); dia = np.asarray(dia, float); box = list(map(float, box))
    d = len(box); N = len(dia)
    _f32(pos.ravel(), out / "pos.f32"); _f32(dia, out / "dia.f32")
    man = {"schema": 1, "ndim": d, "box": box, "N": N,
           "positions": {"file": "pos.f32", "dtype": "float32", "shape": [N, d]},
           "diameters": {"file": "dia.f32", "dtype": "float32", "shape": [N]},
           "phi": float((np.pi/6 if d == 3 else np.pi/4) * (dia**d).sum() / np.prod(box)),
           "layers": []}
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
    return write_bundle(out_dir, pos, dia, box, layers)

def build_layers(pos, dia, box, analyses=None):
    """Compute analysis layers via the packing_analysis registry. `analyses` =
    list of names; None -> all registered analyses. Each missing dep / failure
    skips-with-warning (the registry handles it), so a bundle never breaks."""
    import packing_analysis as pa
    names = pa.available_analyses() if analyses is None else analyses
    return pa.run_analyses(names, pos, dia, box)

def build_pore_layer(pos, dia, box, dcfrac=0.35, frac=0.35):
    """Back-compat shim — the pore-map analysis now lives in packing_analysis.registry."""
    import packing_analysis as pa
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
