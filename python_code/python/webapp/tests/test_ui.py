"""UI audit — every control driven by the REAL DOM action; PNG export verified by a
REAL download. Stats over the real viewport. No proxies."""
import sys, pathlib, io
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import harness as H
from PIL import Image
import numpy as np

def run():
    res = []
    def check(n, c, d=""): res.append((n, c)); print(f"  [{'PASS' if c else 'FAIL'}] {n}  {d}")
    with H.serve() as base:
        # ---- 3D bundle: stats, color-by, palette, slice slider, export ----
        with H.page() as pg:
            diag = H.load(pg, base, "pl_S40_3d", ui=True); pg.wait_for_timeout(300)
            check("UI builds, no console errors", len(pg._errs) == 0, str(pg._errs[:2]))
            st = pg.evaluate("window.__stats")
            check("D1 stats over real FOV (n==visible)", st and st["n"] == diag["visible"], f"{st and st['n']} vs {diag['visible']}")
            check("D1 phi-hist sums to n", sum(st["phi_hist"]["counts"]) == st["n"])
            # D2 real zoom -> stats shrink
            H.set_view(pg, {"target": diag["target"], "zoom": diag["zoom"] + 2}); pg.wait_for_timeout(300)
            st2 = pg.evaluate("window.__stats")
            check("D2 stats update on zoom", st2["n"] < st["n"], f"n {st['n']}->{st2['n']}")
            # B1 color-by via REAL select
            pg.select_option("#colorBy", "local_phi"); pg.wait_for_timeout(250)
            check("B1 color-by select (real)", pg.evaluate("window.__diag.colorBy") == "local_phi")
            # B1 palette via REAL select -> pixels change
            a1, _ = H.shot(pg, "ui_viridis"); pg.select_option("#palette", "turbo"); pg.wait_for_timeout(250)
            a2, _ = H.shot(pg, "ui_turbo")
            check("B1 palette select changes colors (real)", float(np.abs(a1.astype(int) - a2.astype(int)).mean()) > 3)
            # B1 slice slider via REAL input event -> slice z changes
            z0 = pg.evaluate("window.__diag.slice.z")
            pg.evaluate("const s=document.getElementById('sliceSlider'); s.value=200; s.dispatchEvent(new Event('input'))")
            pg.wait_for_timeout(300)
            check("B1 slice slider (real input)", pg.evaluate("window.__diag.slice.z") != z0, f"z {z0:.2f}->{pg.evaluate('window.__diag.slice.z'):.2f}")
            # E1 PNG export — REAL download triggered by the button
            with pg.expect_download() as di:
                pg.click("#exportBtn")
            path = di.value.path(); im = Image.open(path)
            exp_w = round((800 - 300) * 150 / 96)
            check("E1 export downloads a real PNG at DPI", im.format == "PNG" and abs(im.size[0] - exp_w) < 80, f"{im.format} {im.size} expect_w~{exp_w}")
        # ---- 2D bundle: overlay toggle via REAL checkbox ----
        with H.page() as pg:
            H.load(pg, base, "bi2d_r10", ui=True); pg.wait_for_timeout(300)
            a0, _ = H.shot(pg, "ui_pores_off")
            pg.check("#ov_pores"); pg.wait_for_timeout(350)
            a1, _ = H.shot(pg, "ui_pores_on")
            on = "pores" in pg.evaluate("window.__diag.activeOverlays")
            changed = float(np.abs(a0.astype(int) - a1.astype(int)).mean()) > 5
            check("B1 pore-overlay checkbox (real click)", on and changed, f"active={on} dimg>{changed}")
    npass = sum(1 for _, c in res if c)
    print(f"\n{npass}/{len(res)} checks passed")
    return all(c for _, c in res)

if __name__ == "__main__":
    sys.exit(0 if run() else 1)
