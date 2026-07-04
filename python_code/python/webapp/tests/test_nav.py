"""Navigation correctness: zoom scaling (A4), pan fidelity (A5), deep-zoom crispness."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import harness as H

def run():
    res = []
    def check(n, c, d=""): res.append((n, c)); print(f"  [{'PASS' if c else 'FAIL'}] {n}  {d}")
    with H.serve() as base, H.page() as pg:
        H.load(pg, base, "ground_truth")   # 10x10 box, BL circle at (2.5,2.5) r=1.5
        cx, cy = 400, 400                  # canvas center
        # A4 zoom scaling: center BL, measure colored radius at zoom Z and Z+1 -> ~2x
        H.set_view(pg, {"target": [2.5, 2.5, 0], "zoom": 6})
        arr, _ = H.shot(pg, "nav_z6"); r6 = H.colored_radius(arr, cx, cy)
        H.set_view(pg, {"target": [2.5, 2.5, 0], "zoom": 7})
        arr, _ = H.shot(pg, "nav_z7"); r7 = H.colored_radius(arr, cx, cy)
        ratio = r7 / r6 if r6 > 0 else 0
        check("A4 zoom scaling x2", abs(ratio - 2.0) < 0.15, f"r6={r6} r7={r7} ratio={ratio:.2f}")
        # expected absolute: r ~ 1.5 * 2^Z (px). At Z=6: ~96
        check("A4 absolute radius", abs(r6 - 1.5 * 64) < 12, f"r6={r6} (expect ~96)")
        # A5 pan: at fixed zoom, center BL -> center colored; pan +5 world in x -> BR at center
        H.set_view(pg, {"target": [2.5, 2.5, 0], "zoom": 6})
        arr, _ = H.shot(pg, "nav_panA"); a = H.is_colored(arr, cx, cy)
        H.set_view(pg, {"target": [7.5, 2.5, 0], "zoom": 6})   # pan to BR circle
        arr, _ = H.shot(pg, "nav_panB"); b = H.is_colored(arr, cx, cy)
        check("A5 pan fidelity", a and b, f"BL@center={a} BR@center={b}")
        # deep zoom on bi2d: zoom way in, still renders (not blank), crisp
        d2 = H.load(pg, base, "bi2d_r10")
        H.set_view(pg, {"target": [d2["box"][0] / 2, d2["box"][1] / 2, 0], "zoom": 11})
        arr, _ = H.shot(pg, "nav_deep"); frac = H.colored_fraction(arr)
        check("deep-zoom renders (not blank)", 0.05 < frac < 0.999, f"colored_frac={frac:.3f}")
    npass = sum(1 for _, c in res if c)
    print(f"\n{npass}/{len(res)} checks passed")
    return all(c for _, c in res)

if __name__ == "__main__":
    sys.exit(0 if run() else 1)
