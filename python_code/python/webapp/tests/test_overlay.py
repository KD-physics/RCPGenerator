"""Overlays: color-by Voronoi local-phi (C2), palette switch (C4)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import harness as H
import numpy as np

def color_variance(arr, bg=H.BG, tol=60):
    diff = np.abs(arr.astype(int) - np.array(bg)).sum(axis=2)
    m = diff > tol
    if m.sum() < 100: return 0.0
    px = arr[m].astype(float)
    return float(px[:, 0].std() + px[:, 2].std())   # spread across R and B channels

def run():
    res = []
    def check(n, c, d=""): res.append((n, c)); print(f"  [{'PASS' if c else 'FAIL'}] {n}  {d}")
    with H.serve() as base, H.page() as pg:
        H.load(pg, base, "pl_S40_3d")          # 3D slice, carries local_phi layer
        arr0, _ = H.shot(pg, "ov_solid"); v_solid = color_variance(arr0)
        # color by local_phi
        dom = pg.evaluate("VIZ.applyColorBy(window.__viz,'local_phi','viridis')")
        pg.wait_for_timeout(250); arr1, _ = H.shot(pg, "ov_localphi"); v_vir = color_variance(arr1)
        d1 = pg.evaluate("window.__diag")
        check("C2 colorBy local_phi active", d1["colorBy"] == "local_phi", f"domain={dom['domain']}")
        check("C2 coloring varies (vs solid)", v_vir > v_solid + 15, f"solid={v_solid:.1f} viridis={v_vir:.1f}")
        check("C2 no console errors", len(pg._errs) == 0, str(pg._errs[:2]))
        # switch palette -> different pixels
        pg.evaluate("VIZ.applyColorBy(window.__viz,'local_phi','turbo')")
        pg.wait_for_timeout(250); arr2, _ = H.shot(pg, "ov_turbo")
        diff = float(np.abs(arr1.astype(int) - arr2.astype(int)).mean())
        check("C4 palette switch changes colors", diff > 3, f"mean|viridis-turbo|={diff:.1f}")
        # color by diameter works
        pg.evaluate("VIZ.applyColorBy(window.__viz,'diameter','viridis')")
        pg.wait_for_timeout(200); d3 = pg.evaluate("window.__diag")
        check("C4 color-by diameter", d3["colorBy"] == "diameter")
        # revert to solid
        pg.evaluate("VIZ.applyColorBy(window.__viz,null)")
        pg.wait_for_timeout(200); arrR, _ = H.shot(pg, "ov_revert"); vR = color_variance(arrR)
        check("revert to solid fill", vR < v_solid + 10 and pg.evaluate("window.__diag").get("colorBy") is None, f"var={vR:.1f}")
    npass = sum(1 for _, c in res if c)
    print(f"\n{npass}/{len(res)} checks passed")
    return all(c for _, c in res)

if __name__ == "__main__":
    sys.exit(0 if run() else 1)
