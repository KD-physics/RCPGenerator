"""Segmented-pore-map overlay (C3): porenet polygons render + per-pore stats."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import harness as H
import numpy as np

def run():
    res = []
    def check(n, c, d=""): res.append((n, c)); print(f"  [{'PASS' if c else 'FAIL'}] {n}  {d}")
    with H.serve() as base, H.page() as pg:
        H.load(pg, base, "bi2d_r10")
        has = pg.evaluate("!!(window.__viz.currentData().layers.pores)")
        check("C3 pores layer present", has)
        arr0, _ = H.shot(pg, "pores_off")
        pg.evaluate("window.__viz.state.activeOverlays=['pores']; window.__viz.refresh()")
        pg.wait_for_timeout(400); arr1, _ = H.shot(pg, "pores_on")
        diff = float(np.abs(arr0.astype(int) - arr1.astype(int)).mean())
        check("C3 pore overlay renders (changes view)", diff > 8, f"mean|on-off|={diff:.1f}")
        ps = pg.evaluate("VIZ.poreStats(window.__viz.currentData())")
        # C3 cross-check vs toolbox: porenet on same packing/frac gives same body count
        sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent.parent / "Scripts"))
        import voidfill_method as vf, porenet as pn
        pos, dia, box = vf._load(str(pathlib.Path(__file__).resolve().parent.parent.parent / "Data/voidfill/bi2d_r10.npz"))
        info = pn.assign_fines(pn.pore_network(pos, dia, box, 0.35 * dia.max()))
        info = pn.merge_pores_gap(info, pn.throats(info), frac=0.35)
        nb = info["merge"]["n_bodies"]
        check("C3 pore count matches toolbox", ps["n_bodies"] == nb, f"app={ps['n_bodies']} toolbox={nb}")
        check("C3 no console errors", len(pg._errs) == 0, str(pg._errs[:2]))
    npass = sum(1 for _, c in res if c)
    print(f"\n{npass}/{len(res)} checks passed")
    return all(c for _, c in res)

if __name__ == "__main__":
    sys.exit(0 if run() else 1)
