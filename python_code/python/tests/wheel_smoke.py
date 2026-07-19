"""Fast installed-wheel smoke test used by cibuildwheel."""

import math

import numpy as np

import rcpgenerator


def main() -> None:
    rcpgenerator.set_num_threads(1)
    packing = rcpgenerator.Packing(
        phi=0.08,
        N=16,
        Ndim=2,
        box=[1.0, 1.0],
        walls=[0, 0],
        dist={"type": "mono", "d": 1.0},
        seed=20240719,
    )
    result = packing.pack()

    positions = np.asarray(result["positions"], dtype=float)
    diameters = np.asarray(result["diameters"], dtype=float)
    assert positions.shape == (16, 2), positions.shape
    assert diameters.shape == (16,), diameters.shape
    assert np.isfinite(positions).all()
    assert np.isfinite(diameters).all()
    assert (diameters > 0.0).all()
    assert math.isfinite(float(result["phi"]))
    assert math.isfinite(float(result["force_magnitude"]))
    assert int(result["steps"]) > 0
    for history_name in ("phi_history", "force_history", "energy_history"):
        history = np.asarray(result[history_name], dtype=float)
        assert history.shape == (int(result["steps"]),), (
            history_name,
            history.shape,
            result["steps"],
        )
        if history_name == "force_history":
            # Contact-free phases use NaN for the undefined mean force.
            assert not np.isinf(history).any()
            assert math.isfinite(float(history[-1]))
        else:
            assert np.isfinite(history).all(), history_name


if __name__ == "__main__":
    main()
