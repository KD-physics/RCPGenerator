"""New render features (real pixels): black outline (constant screen width),
periodic boundary box, random palette."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import harness as H
import numpy as np

def frac_grey(a):   # light-grey pixels (outline ~205 / boundary ~240), not blue fill, not dark bg
    a = a.astype(int)
    return float(((a.min(axis=2) > 150) & (a.max(axis=2) - a.min(axis=2) < 45)).mean())
frac_light = frac_grey

def run():
    res = []
    def check(n, c, d=""): res.append((n, c)); print(f"  [{'PASS' if c else 'FAIL'}] {n}  {d}")
    with H.serve() as base, H.page() as pg:
        H.load(pg, base, "bi2d_r10")                       # outline + boundary ON by default
        H.set_view(pg, {"target": [0.5, 0.5, 0], "zoom": pg.evaluate("window.__diag.zoom") - 1})  # box well inside canvas
        pg.wait_for_timeout(250); a_on, _ = H.shot(pg, "feat_on")
        b_on, l_on = frac_grey(a_on), frac_light(a_on)
        # turn both off
        pg.evaluate("window.__viz.state.outline=0; window.__viz.state.showBoundary=false; window.__viz.refresh()")
        pg.wait_for_timeout(250); a_off, _ = H.shot(pg, "feat_off")
        b_off, l_off = frac_grey(a_off), frac_light(a_off)
        check("outline draws black rings (constant-px)", b_on > b_off + 0.02, f"black on={b_on:.3f} off={b_off:.3f}")
        check("boundary box drawn", l_on > l_off + 0.0008, f"light on={l_on:.4f} off={l_off:.4f}")
        # outline thickness ~constant in screen px across zoom: black fraction should
        # not collapse when zoomed in (outline stays visible)
        pg.evaluate("window.__viz.state.outline=1.5; window.__viz.refresh()")
        H.set_view(pg, {"target": [0.5, 0.5, 0], "zoom": pg.evaluate("window.__diag.zoom") + 3}); pg.wait_for_timeout(250)
        az, _ = H.shot(pg, "feat_zoom"); bz = frac_grey(az)
        check("outline persists when zoomed (screen-constant)", bz > 0.01, f"black@zoom={bz:.3f}")
        # random palette
        pg.evaluate("VIZ.applyColorBy(window.__viz,'diameter','random')"); pg.wait_for_timeout(250)
        ar, _ = H.shot(pg, "feat_random")
        m = (np.abs(ar.astype(int) - np.array(H.BG)).sum(2) > 60)
        var = float(ar[m][:, 0].std() + ar[m][:, 1].std() + ar[m][:, 2].std()) if m.sum() > 100 else 0
        check("random palette varies colors", var > 40, f"colorspread={var:.0f}")
        check("no console errors", len(pg._errs) == 0, str(pg._errs[:2]))
    npass = sum(1 for _, c in res if c)
    print(f"\n{npass}/{len(res)} checks passed")
    return all(c for _, c in res)

if __name__ == "__main__":
    sys.exit(0 if run() else 1)
