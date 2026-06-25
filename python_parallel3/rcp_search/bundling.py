"""Optional save/bundle layer for generated packings.

A packing sweep by default only tabulates a census row (saving every packing
would be ~100s of GB). This module adds *targeted* saving, behind flags, with no
effect on the default path:

  * save_packing -> {save_dir}/packings/{name}.npz   (pos/dia/box + provenance)
  * save_bundle  -> {save_dir}/bundles/{name}/       (webapp bundle, optional
                    analysis layers via the packing_analysis registry)

Everything is LAZY-imported (webapp, packing_analysis) so importing rcp_search
never hard-depends on them, and a bundle failure is recorded, never raised, so a
sweep is never lost. Environment-neutral: all paths are caller-supplied.
"""
import json
import pathlib
import re
import numpy as np


def _safe_name(name):
    """Filesystem-safe name (the case_key uses '|', invalid on Windows where the
    files often land). Replace reserved chars with '_'."""
    return re.sub(r'[<>:"/\\|?*]', '_', str(name))


def _subdir(save_dir, sub):
    p = pathlib.Path(save_dir) / sub
    p.mkdir(parents=True, exist_ok=True)
    return p


def save_npz(out_path, pos, dia, box, provenance=None):
    """Write a packing .npz (positions/diameters/box) with a JSON `provenance`
    string (family, params, N, seed, phi, ...). The toolbox + webapp read it."""
    out_path = pathlib.Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    prov = json.dumps(provenance or {}, default=str)
    np.savez(out_path,
             positions=np.asarray(pos, float),
             diameters=np.asarray(dia, float),
             box=np.asarray(box, float),
             provenance=np.array(prov))
    return str(out_path)


def bundle_from_npz(npz_path, out_dir, analyses=None):
    """Standalone: (re)build a webapp bundle from a saved .npz, attaching any
    requested analyses (registry names, e.g. ["local_phi","pores"]). Callable
    anywhere, anytime — decouples analysis from generation."""
    from webapp import bridge   # lazy
    return bridge.from_npz(str(npz_path), str(out_dir),
                           with_layers=bool(analyses), analyses=analyses)


def save_outputs(packing, provenance, bundle_opts):
    """Write the artifacts requested by `bundle_opts` for a packed object.

    bundle_opts (dict, all optional):
      save_packing : bool   -> write the .npz
      save_bundle  : bool   -> write the webapp bundle
      analyses     : list   -> analysis layer names (registry); [] = geometry only
      save_dir     : path   -> root for packings/ and bundles/ (default ".")
      name         : str    -> base filename (default provenance["case_key"])

    Returns a dict of written paths (and any bundle_error). Never raises."""
    if not bundle_opts:
        return {}
    save_dir = bundle_opts.get("save_dir") or "."
    name = _safe_name(bundle_opts.get("name") or provenance.get("case_key") or "packing")
    analyses = bundle_opts.get("analyses") or []
    pos = np.asarray(packing.positions, float)
    dia = np.asarray(packing.diameters, float)
    box = np.asarray(packing.box, float)
    out = {}
    npz_path = None

    if bundle_opts.get("save_packing"):
        try:
            npz_path = _subdir(save_dir, "packings") / f"{name}.npz"
            save_npz(npz_path, pos, dia, box, provenance=provenance)
            out["npz"] = str(npz_path)
        except Exception as e:
            out["npz_error"] = f"{type(e).__name__}: {str(e)[:120]}"
            npz_path = None

    if bundle_opts.get("save_bundle"):
        try:
            from webapp import bridge   # lazy
            bdir = _subdir(save_dir, "bundles") / name
            if npz_path is not None:
                bridge.from_npz(str(npz_path), str(bdir),
                                with_layers=bool(analyses), analyses=analyses)
            else:
                layers = bridge.build_layers(pos, dia, box, analyses) if analyses else []
                bridge.write_bundle(str(bdir), pos, dia, box, layers)
            out["bundle"] = str(bdir)
        except Exception as e:
            out["bundle_error"] = f"{type(e).__name__}: {str(e)[:120]}"

    return out
