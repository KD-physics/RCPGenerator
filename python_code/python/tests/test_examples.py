"""Every snippet in rcpgenerator.example() runs as printed; boundary-input checks (walls, fix_height, small box).

    python tests/test_examples.py            # all snippets (a few minutes)
    python tests/test_examples.py --quick    # catalogue checks and one snippet per group
    pytest tests/test_examples.py            # same; RCP_QUICK=1 for the subset

Each snippet runs in its own subprocess, in a temporary directory, with a non-interactive
Matplotlib backend and one thread. The 'voronoi' snippet is skipped when pyvoro is not installed.
"""

from __future__ import annotations

import io
import os
import subprocess
import sys
import tempfile
import time
from contextlib import redirect_stdout

import rcpgenerator
from rcpgenerator.examples import _CATALOGUE, example_names, example_source

QUICK = os.environ.get("RCP_QUICK", "") not in ("", "0")
MAX_LINES = 12


def _have(module: str) -> bool:
    try:
        __import__(module)
        return True
    except ImportError:
        return False


def run_snippet(name: str) -> tuple[bool, str, float]:
    code = example_source(name)
    env = dict(os.environ, MPLBACKEND="Agg", OMP_NUM_THREADS="1")
    with tempfile.TemporaryDirectory() as tmp:
        if name == "from_file":
            with open(os.path.join(tmp, "diameters.txt"), "w") as f:
                f.write("\n".join(str(1.0 + 0.002 * k) for k in range(500)))
        t0 = time.time()
        proc = subprocess.run([sys.executable, "-c", code], cwd=tmp, env=env, capture_output=True, text=True,
                              timeout=900)
        return proc.returncode == 0, (proc.stdout + proc.stderr)[-2000:], time.time() - t0


def test_catalogue_format():
    buf = io.StringIO()
    with redirect_stdout(buf):
        rcpgenerator.example()
    text = buf.getvalue()
    for title, entries in _CATALOGUE:
        assert title in text
        for name, desc, code in entries:
            assert f"  - '{name}': {desc}" in text
            assert 1 <= len(desc.split()) <= 6, f"{name}: descriptor should be a few words"
            n = len(code.strip().splitlines())
            assert 3 <= n <= MAX_LINES, f"{name}: {n} lines"
            assert "import rcpgenerator" in code or "rcpgenerator" in code.splitlines()[0], f"{name}: not self-contained"
    assert len(set(example_names())) == len(example_names()), "duplicate example names"


def test_lookup():
    buf = io.StringIO()
    with redirect_stdout(buf):
        rcpgenerator.example("disk")
    assert buf.getvalue() == example_source("disk")
    try:
        rcpgenerator.example("cylindr")
        raise AssertionError("an unknown name must raise")
    except ValueError as exc:
        assert "cylinder" in str(exc), str(exc)


def test_walls_validation():
    bad = {(2, (0, -2)): "FIRST", (3, (0, 0, -2)): "FIRST", (2, (2, 0)): "not 0 or 1", (2, (-2,)): "entries but",
           (3, (-4, 0, 0)): "container in 3", (2, (0.5, 0)): "integers"}
    for (nd, walls), fragment in bad.items():
        try:
            rcpgenerator.Packing(N=50, Ndim=nd, phi=0.05, box=[1.0] * nd, walls=list(walls), seed=1)
            raise AssertionError(f"walls={walls} must raise")
        except ValueError as exc:
            assert fragment in str(exc), (walls, str(exc))
    for nd, walls in ((2, [-1, 0]), (2, [-2, 1]), (3, [-2, -1, 1]), (3, [-3, 0, 0]), (3, [0, 1, 1])):
        rcpgenerator.Packing(N=50, Ndim=nd, phi=0.05, box=[1.0] * nd, walls=walls, seed=1)   # accepted as before


def test_fix_height():
    """fix_height keeps the last box length a fixed multiple of the largest diameter through a run (broken in
    1.0.0: the height compounded to inf)."""
    import math
    for nd, box, walls in ((3, [1.0, 1.0, 4.0], [0, 0, 1]), (2, [1.0, 6.0], [0, 1])):
        p = rcpgenerator.Packing(N=300, Ndim=nd, phi=0.05 if nd == 3 else 0.1, box=box, walls=walls,
                                 fix_height=True, dist={"type": "mono", "d": 1.0}, seed=1)
        p.pack()
        ratio = p.box[-1] / max(p.diameters)
        assert abs(ratio - box[-1]) < 1e-9 * box[-1], f"height / largest diameter = {ratio}, expected {box[-1]}"
        assert math.isfinite(p.phi_final) and (0.4 if nd == 3 else 0.7) < p.phi_final < (0.7 if nd == 3 else 0.9)
    p = rcpgenerator.Packing(N=300, Ndim=3, phi=0.05, box=[1.0, 1.0, 4.0], walls=[0, 0, 1], fix_height=True,
                             dist={"type": "mono", "d": 1.0}, seed=1)
    p.relax(n_steps=8000, target_phi=0.5)
    assert abs(p.phi_final - 0.5) < 1e-9 and abs(p.box[-1] / max(p.diameters) - 4.0) < 1e-9


def test_small_periodic_box_warning():
    """A periodic axis too short for the largest particle warns after the run (6D equal sizes need about 3 largest
    diameters; this box has 2.1); the same box with hard walls, and a 2D box, stay silent."""
    import warnings
    def warned(**kw):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            rcpgenerator.Packing(N=200, phi=0.01, seed=3, **kw).pack()
        return [str(w.message) for w in caught if "periodic box is small" in str(w.message)]
    msgs = warned(Ndim=6, box=[1.0] * 6, walls=[0] * 6)
    assert len(msgs) == 1 and "axis 5" in msgs[0], msgs
    assert not warned(Ndim=6, box=[1.0] * 6, walls=[1] * 6)
    assert not warned(Ndim=2, box=[1.0, 1.0], walls=[0, 0])


def test_snippets_run():
    names = example_names()
    if QUICK:
        names = [entries[0][0] for _, entries in _CATALOGUE]
    failed = []
    for name in names:
        if name == "voronoi" and not _have("pyvoro"):
            print(f"  skip  {name} (pyvoro not installed)")
            continue
        ok, out, dt = run_snippet(name)
        if ok and "periodic box is small" in out:
            ok, out = False, "a snippet should not trigger the small periodic box warning\n" + out
        print(f"  {'ok  ' if ok else 'FAIL'}  {name} ({dt:.1f} s)", flush=True)
        if not ok:
            failed.append(f"{name}:\n{out}")
    assert not failed, "\n\n".join(failed)


TESTS = [test_catalogue_format, test_lookup, test_walls_validation, test_fix_height, test_small_periodic_box_warning,
         test_snippets_run]


def main(argv: list[str]) -> int:
    global QUICK
    QUICK = QUICK or "--quick" in argv
    rcpgenerator.set_num_threads(1)
    failures = 0
    for test in TESTS:
        t0 = time.time()
        try:
            test()
            print(f"PASS  {test.__name__}  ({time.time() - t0:.1f} s)", flush=True)
        except Exception as exc:
            failures += 1
            print(f"FAIL  {test.__name__}: {type(exc).__name__}: {exc}", flush=True)
    print(f"{len(TESTS) - failures}/{len(TESTS)} passed" + (" (quick)" if QUICK else ""))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
