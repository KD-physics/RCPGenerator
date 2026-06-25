#!/usr/bin/env python3
"""Capture a curated gallery of views for human (Ken) acceptance + real-RTX check."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent / "tests"))
import harness as H
from PIL import Image, ImageDraw

OUT = pathlib.Path(__file__).resolve().parent / "gallery"; OUT.mkdir(exist_ok=True)

def cap(pg, name, label):
    p = OUT / f"{name}.png"; pg.screenshot(path=str(p))
    return (p, label)

def run():
    shots = []
    with H.serve() as base, H.page(viewport=(900, 800)) as pg:
        # 2D bidisperse, full
        H.load(pg, base, "bi2d_r10"); shots.append(cap(pg, "g_bi2d", "2D bidisperse packing (5k)"))
        # deep zoom into a pore
        H.set_view(pg, {"target": [0.5, 0.5, 0], "zoom": 9}); pg.wait_for_timeout(200)
        shots.append(cap(pg, "g_zoom", "zoom into a pore (crisp)"))
        # pore-map overlay
        H.set_view(pg, {"target": [0.5, 0.5, 0], "zoom": 7}); pg.wait_for_timeout(100)
        pg.evaluate("window.__viz.state.activeOverlays=['pores']; window.__viz.refresh()"); pg.wait_for_timeout(300)
        shots.append(cap(pg, "g_pores", "segmented-pore map (eta)"))
        # 3D slice colored by local-phi
        H.load(pg, base, "pl_S40_3d"); pg.evaluate("VIZ.applyColorBy(window.__viz,'local_phi','viridis')"); pg.wait_for_timeout(300)
        shots.append(cap(pg, "g_localphi", "3D slice: Voronoi local-phi"))
        # 3D slice at another depth
        info = pg.evaluate("window.__viz.sliceInfo()"); pg.evaluate("z=>window.__viz.setSlice(z)", info["box"][2] * 0.25)
        pg.evaluate("VIZ.applyColorBy(window.__viz,'local_phi','turbo')"); pg.wait_for_timeout(300)
        shots.append(cap(pg, "g_slice2", "3D slice (z=0.25, turbo)"))
        # full UI
        H.load(pg, base, "pl_S40_3d", ui=True); pg.wait_for_timeout(300)
        pg.evaluate("VIZ.applyColorBy(window.__viz,'local_phi','viridis'); VIZ.stats.draw()"); pg.wait_for_timeout(300)
        shots.append(cap(pg, "g_ui", "full app: render + controls + FOV stats"))
    # montage 2x3
    ims = [Image.open(p).convert("RGB").resize((440, 391)) for p, _ in shots]
    labels = [l for _, l in shots]
    W, Hh, cols = 440, 391, 3; rows = (len(ims) + cols - 1) // cols
    canvas = Image.new("RGB", (W * cols, (Hh + 22) * rows), (20, 20, 24))
    dr = ImageDraw.Draw(canvas)
    for i, (im, lab) in enumerate(zip(ims, labels)):
        x, y = (i % cols) * W, (i // cols) * (Hh + 22)
        canvas.paste(im, (x, y + 22)); dr.text((x + 6, y + 6), lab, fill=(230, 230, 235))
    out = OUT / "gallery.png"; canvas.save(out); print(f"gallery -> {out}")

if __name__ == "__main__":
    run()
