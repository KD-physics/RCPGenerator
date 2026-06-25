#!/usr/bin/env python3
"""CLI launcher. Point it at a packing (.npz) or a bundle dir; it produces a
self-contained HTML with the data EMBEDDED (no server, standalone) and opens it.
  python launch.py <packing.npz | bundle_dir> [--no-open] [--out FILE]
"""
import sys, os, json, base64, pathlib, tempfile, webbrowser
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build as B
import bridge as BR

def make_bundle(arg):
    p = pathlib.Path(arg)
    if p.is_dir():
        return p
    out = pathlib.Path(tempfile.mkdtemp(prefix="vizbundle_"))
    BR.from_npz(str(p), str(out), with_layers=True)   # .npz -> bundle (+local-phi)
    return out

def embed(bundle_dir, out_file=None):
    man = json.loads((bundle_dir / "manifest.json").read_text())
    files = {"manifest.json": None}
    need = [man["positions"]["file"], man["diameters"]["file"]] + [L["file"] for L in man.get("layers", [])]
    enc = {}
    for f in need:
        enc[f] = base64.b64encode((bundle_dir / f).read_bytes()).decode()
    payload = {"manifest": man, "files": enc}
    standalone = B.build().read_text()
    inject = "<script>window.__EMBEDDED_BUNDLE=" + json.dumps(payload) + ";</script>\n"
    html = standalone.replace("<body>", "<body>\n" + inject, 1)
    out = pathlib.Path(out_file) if out_file else (B.DIST / ("packing_" + bundle_dir.name + ".html"))
    out.write_text(html)
    print(f"standalone -> {out}  ({len(html)//1024} KB, data embedded)")
    return out

if __name__ == "__main__":
    a = [x for x in sys.argv[1:] if not x.startswith("--")]
    opts = [x for x in sys.argv[1:] if x.startswith("--")]
    out_file = None
    if "--out" in sys.argv: out_file = sys.argv[sys.argv.index("--out") + 1]
    bundle = make_bundle(a[0])
    html = embed(bundle, out_file)
    if "--no-open" not in opts:
        webbrowser.open("file://" + str(html.resolve()))
