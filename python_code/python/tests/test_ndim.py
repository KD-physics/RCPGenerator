"""4D+ neighbour search: the k-d tree backend reproduces the generic sort+walk search exactly.

For Ndim >= 4 the default pair search is a dimension-generic k-d tree; RCP_PAIRS_BACKEND=sortwalk selects the
previous generic search. Both must produce identical packings (same steps, bit-identical positions and diameters),
and the automatic neighbour capacity must hold every neighbour list. The backend is read once per process, so
each run is a subprocess.

    python tests/test_ndim.py      # ~1 minute
    pytest tests/test_ndim.py
"""

from __future__ import annotations

import json
import os
import subprocess
import sys

RUNNER = r"""
import hashlib, json, sys
import numpy as np, rcpgenerator
rcpgenerator.set_num_threads(1)
c = json.loads(sys.argv[1])
p = rcpgenerator.Packing(N=c["N"], Ndim=c["Ndim"], phi=0.02, box=[1.0] * c["Ndim"], walls=c["walls"],
                         dist=c["dist"], seed=1)
p.pack()
print(json.dumps(dict(steps=p.steps, phi=p.phi_final,
      sha=hashlib.sha1(np.asarray(p.positions).tobytes() + np.asarray(p.diameters).tobytes()).hexdigest())))
"""

CASES = [
    dict(Ndim=4, N=400, walls=[0, 0, 0, 0], dist={"type": "bidisperse", "d1": 1.0, "d2": 1.4, "p": 0.5}),
    dict(Ndim=4, N=400, walls=[0, 0, 0, 1], dist={"type": "mono", "d": 1.0}),
    dict(Ndim=5, N=400, walls=[0] * 5, dist={"type": "mono", "d": 1.0}),
]


def run(case: dict, backend: str | None) -> dict:
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("RCP_PAIRS_BACKEND", None)
    if backend:
        env["RCP_PAIRS_BACKEND"] = backend
    proc = subprocess.run([sys.executable, "-c", RUNNER, json.dumps(case)], env=env, capture_output=True,
                          text=True, timeout=900)
    assert proc.returncode == 0, proc.stderr[-1500:]
    return json.loads(proc.stdout.strip().splitlines()[-1])


def test_kdtree_matches_sortwalk():
    for case in CASES:
        kd, sw = run(case, None), run(case, "sortwalk")
        assert kd == sw, f"{case['Ndim']}D {case['walls']}: k-d tree {kd} vs sort+walk {sw}"
        lo, hi = (0.35, 0.6) if case["Ndim"] == 4 else (0.2, 0.4)
        assert lo < kd["phi"] < hi, f"implausible packing fraction {kd['phi']}"


if __name__ == "__main__":
    test_kdtree_matches_sortwalk()
    print("PASS  test_kdtree_matches_sortwalk")
