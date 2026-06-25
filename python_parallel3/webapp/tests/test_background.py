"""Background-colour control: real dropdown selection changes the render-area bg and
auto-flips the outline/boundary ink for contrast."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import harness as H
import numpy as np

def run():
    res = []
    def check(n, c, d=""): res.append((n, c)); print(f"  [{'PASS' if c else 'FAIL'}] {n}  {d}")
    with H.serve() as base, H.page() as pg:
        H.load(pg, base, "bi2d_r10", ui=True); pg.wait_for_timeout(300)
        opts = pg.eval_on_selector_all("#bgSelect option", "els=>els.map(e=>e.textContent)")
        check("background control present with presets", "white" in opts and "black" in opts, f"{opts}")
        # default dark -> outline light
        ink0 = pg.evaluate("window.__viz.state.outlineColor")
        check("dark bg -> light outline ink", ink0[0] > 150, f"ink={ink0}")
        # select white via the REAL dropdown
        pg.select_option("#bgSelect", label="white"); pg.wait_for_timeout(300)
        bg = pg.evaluate("getComputedStyle(document.getElementById('deck')).backgroundColor")
        check("render-area bg changed to white", "255, 255, 255" in bg, f"bg={bg}")
        ink1 = pg.evaluate("window.__viz.state.outlineColor")
        check("light bg -> dark outline ink (auto-contrast)", ink1[0] < 80, f"ink={ink1}")
        a, _ = H.shot(pg, "bg_white")
        # the voids (between particles) are now white -> a real fraction of near-white pixels
        white = float((a.min(axis=2) > 240).mean())
        check("render shows white background in voids", white > 0.02, f"white_frac={white:.3f}")
        check("no console errors", len(pg._errs) == 0, str(pg._errs[:2]))
    npass = sum(1 for _, c in res if c)
    print(f"\n{npass}/{len(res)} checks passed")
    return all(c for _, c in res)

if __name__ == "__main__":
    sys.exit(0 if run() else 1)
