"""Core rendering correctness: render, area-fraction = phi, ground-truth pixels, count."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import harness as H

def run():
    results = []
    def check(name, cond, detail=""):
        results.append((name, cond, detail)); print(f"  [{'PASS' if cond else 'FAIL'}] {name}  {detail}")

    with H.serve() as base, H.page() as pg:
        # --- ground truth (4 known circles in a 10x10 box, phi=0.1414) ---
        diag = H.load(pg, base, "ground_truth")
        H.no_outline(pg)
        arr, _ = H.shot(pg, "gt")
        check("GT load+render", diag["ready"] and diag["rendered"], f"N={diag['N']}")
        check("GT no console errors", len(pg._errs) == 0, str(pg._errs[:2]))
        af = H.coverage_fraction(arr)
        check("GT area-fraction = phi", abs(af - diag["phi"]) < 0.025, f"area={af:.3f} phi={diag['phi']:.3f}")
        # four circles at the quadrant centers; box center (5,5) is empty
        # canvas 800x800, box 10x10 fills it -> world (x,y)->px (x*80, (10-y)*80) [y down]
        def w2p(x, y): return int(x * 80), int((10 - y) * 80)
        for (wx, wy, lbl) in [(2.5, 2.5, "BL"), (7.5, 2.5, "BR"), (2.5, 7.5, "TL"), (7.5, 7.5, "TR")]:
            px, py = w2p(wx, wy)
            check(f"GT circle {lbl} colored", H.is_colored(arr, px, py), f"@({px},{py})")
        cx, cy = w2p(5, 5)
        check("GT center empty", not H.is_colored(arr, cx, cy), f"@({cx},{cy})")
        check("GT visible count = N", diag["visible"] == diag["N"], f"vis={diag['visible']}")

        # --- bi2d_r10 (5000 circles, phi=0.9165) ---
        diag2 = H.load(pg, base, "bi2d_r10")
        H.no_outline(pg)
        arr2, _ = H.shot(pg, "bi2d")
        check("bi2d load+render", diag2["rendered"], f"N={diag2['N']}")
        af2 = H.coverage_fraction(arr2)
        check("bi2d area-fraction = phi", abs(af2 - diag2["phi"]) < 0.025, f"area={af2:.3f} phi={diag2['phi']:.3f}")
        check("bi2d visible count = N", diag2["visible"] == diag2["N"], f"vis={diag2['visible']}")

    npass = sum(1 for _, c, _ in results if c)
    print(f"\n{npass}/{len(results)} checks passed")
    return all(c for _, c, _ in results)

if __name__ == "__main__":
    sys.exit(0 if run() else 1)
