"""Render-class palettes in the color-by dropdown (real DOM selection -> the palette's
colors appear on screen), matching rcpgenerator/render.py _PALETTES."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import harness as H
import numpy as np

# the 12 render palettes (must match app/js/palettes.js)
PAL5 = [[255,99,71],[255,165,0],[255,223,0],[192,241,192],[255,250,205],[50,205,50]]

def present(arr, rgb, tol=40):
    d = np.abs(arr.astype(int) - np.array(rgb)).sum(axis=2)
    return float((d < tol).mean())

def run():
    res = []
    def check(n, c, d=""): res.append((n, c)); print(f"  [{'PASS' if c else 'FAIL'}] {n}  {d}")
    with H.serve() as base, H.page() as pg:
        H.load(pg, base, "bi2d_r10", ui=True); pg.wait_for_timeout(300)
        # the dropdown should contain palette options
        opts = pg.eval_on_selector_all("#colorBy option", "els => els.map(e=>e.value)")
        check("color-by dropdown lists render palettes", "p1" in opts and "p12" in opts, f"{[o for o in opts if o.startswith('p')][:13]}")
        # select palette 5 via the REAL dropdown
        pg.select_option("#colorBy", "p5"); pg.wait_for_timeout(300)
        check("palette select active", pg.evaluate("window.__diag.colorBy") == "p5")
        a5, _ = H.shot(pg, "pal5")
        # palette-5 colors should appear on screen (random assignment -> several present)
        hit = sum(present(a5, c) > 0.005 for c in PAL5)
        check("palette-5 colors render", hit >= 3, f"{hit}/6 palette colors present")
        # a different palette -> different pixels
        pg.select_option("#colorBy", "p9"); pg.wait_for_timeout(300)
        a9, _ = H.shot(pg, "pal9")
        check("different palette -> different colors", float(np.abs(a5.astype(int) - a9.astype(int)).mean()) > 5,
              f"mean|p5-p9|={float(np.abs(a5.astype(int)-a9.astype(int)).mean()):.1f}")
        check("no console errors", len(pg._errs) == 0, str(pg._errs[:2]))
    npass = sum(1 for _, c in res if c)
    print(f"\n{npass}/{len(res)} checks passed")
    return all(c for _, c in res)

if __name__ == "__main__":
    sys.exit(0 if run() else 1)
