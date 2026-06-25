"""Standalone deliverable (B4) + drag-drop load (C1): single self-contained HTML from
file:// with embedded data (no server), and the drag-drop folder-load path."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))   # viz/ for build,launch
import harness as H
import build, launch

VIZ = pathlib.Path(__file__).resolve().parent.parent

def run():
    res = []
    def check(n, c, d=""): res.append((n, c)); print(f"  [{'PASS' if c else 'FAIL'}] {n}  {d}")

    # --- standalone embedded HTML from file:// (no server) ---
    html = launch.embed(VIZ / "bundles" / "ground_truth", str(VIZ / "dist" / "packing_gt_embed.html"))
    from playwright.sync_api import sync_playwright
    with sync_playwright() as p:
        b = p.chromium.launch(args=H.GL_ARGS + ["--allow-file-access-from-files"])
        pg = b.new_page(viewport={"width": 900, "height": 700}); errs = []
        pg.on("console", lambda m: errs.append(m.text) if m.type == "error" else None)
        pg.goto("file://" + str(pathlib.Path(html).resolve()))
        pg.wait_for_function("window.__diag && window.__diag.rendered===true", timeout=15000)
        d = pg.evaluate("window.__diag")
        check("B4 standalone file:// renders (no server)", d["rendered"] and d["N"] == 4, f"N={d['N']}")
        check("B4 standalone no console errors", len(errs) == 0, str(errs[:2]))
        b.close()

    # --- C1 bundle + layers load (folder/URL path; same normalize() drag-drop uses) ---
    with H.serve() as base, H.page() as pg:
        diag = H.load(pg, base, "pl_S40_3d")     # bundle with a local_phi layer
        ok = pg.evaluate("!!(window.__viz.currentData().layers && window.__viz.currentData().layers.local_phi)")
        check("C1 bundle manifest + layers parsed (real URL load)", diag["N"] > 0 and ok, f"N={diag['N']} layer={ok}")
        # (real picker/drop-extraction load paths are verified end-to-end in test_load.py)

    npass = sum(1 for _, c in res if c)
    print(f"\n{npass}/{len(res)} checks passed")
    return all(c for _, c in res)

if __name__ == "__main__":
    sys.exit(0 if run() else 1)
