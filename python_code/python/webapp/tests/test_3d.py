"""3D cross-section (A6) + large-N load/perf (B3, software-GL floor)."""
import sys, pathlib, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import harness as H

def run():
    res = []
    def check(n, c, d=""): res.append((n, c)); print(f"  [{'PASS' if c else 'FAIL'}] {n}  {d}")
    with H.serve() as base, H.page() as pg:
        # --- 3D slice correctness on pl_S40_3d (53k, 3D) ---
        diag = H.load(pg, base, "pl_S40_3d")
        check("3D detected, renders a 2D slice", diag["ndim"] == 2 and diag["sourceNdim"] == 3,
              f"ndim={diag['ndim']} src={diag['sourceNdim']} slice={diag['slice']}")
        check("slice has particles", diag["visible"] > 100, f"vis={diag['visible']}")
        H.no_outline(pg)
        arr, _ = H.shot(pg, "slice_mid")
        # area-fraction = phi holds in EXPECTATION over planes -> average several interior slices
        Lz = diag["box"]  # slice box is 2D; get source Lz from sliceInfo
        info = pg.evaluate("window.__viz.sliceInfo()"); Lz = info["box"][info["axis"]]
        afs = []
        for frac in (0.25, 0.4, 0.55, 0.7, 0.85):
            pg.evaluate("z => window.__viz.setSlice(z)", frac * Lz); pg.wait_for_timeout(200)
            a, _ = H.shot(pg, f"slice_{int(frac*100)}"); afs.append(H.coverage_fraction(a))
        mean_af = sum(afs) / len(afs)
        check("A6 mean slice area-fraction = phi", abs(mean_af - diag["phi"]) < 0.03,
              f"mean_area={mean_af:.3f} phi={diag['phi']:.3f} slices={[round(x,3) for x in afs]}")
        check("3D no console errors", len(pg._errs) == 0, str(pg._errs[:2]))
        # move the slice -> content changes
        z0 = diag["slice"]["z"]; box = diag["box"]
        v0 = diag["visible"]
        pg.evaluate("z => window.__viz.setSlice(z)", z0 * 0.4)
        pg.wait_for_timeout(300); diag2 = pg.evaluate("window.__diag")
        check("slice slider changes content", diag2["slice"]["z"] != z0 and diag2["visible"] != v0,
              f"z {z0:.2f}->{diag2['slice']['z']:.2f}, vis {v0}->{diag2['visible']}")

        # --- large-N load/perf on pl_S100_3d (~461k, 3D) ---
        t0 = time.time()
        d3 = H.load(pg, base, "pl_S100_3d", timeout=60000)
        load_t = time.time() - t0
        check("large-N loads + renders", d3["rendered"], f"sourceN(slice)={d3['N']} load={load_t:.1f}s (sw-GL floor)")
        check("large-N load under 30s (sw-GL floor)", load_t < 30, f"{load_t:.1f}s")
    npass = sum(1 for _, c in res if c)
    print(f"\n{npass}/{len(res)} checks passed")
    return all(c for _, c in res)

if __name__ == "__main__":
    sys.exit(0 if run() else 1)
