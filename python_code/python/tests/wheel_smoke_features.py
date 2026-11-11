"""Fast installed-wheel smoke test for the public rcpgenerator 1.1 features."""

import io
import json
import warnings
from contextlib import redirect_stdout

import numpy as np

import rcpgenerator


def main() -> None:
    rcpgenerator.set_num_threads(1)

    # The snippet catalogue must be packaged and available from the public API.
    output = io.StringIO()
    with redirect_stdout(output):
        rcpgenerator.example()
        rcpgenerator.example("permute_4d")
    catalogue = output.getvalue()
    assert "Higher Dimensions" in catalogue
    assert "diameter_permute=True" in catalogue

    # Permutation must run while preserving the diameter multiset up to its
    # common packing scale.
    initial_diameters = np.linspace(1.0, 4.0, 300)
    packing = rcpgenerator.Packing(
        N=300,
        Ndim=3,
        phi=0.05,
        box=[1.0] * 3,
        walls=[0] * 3,
        seed=1,
        dist={"type": "custom", "custom": initial_diameters.tolist()},
        diameter_permute=True,
    )
    packing.pack()
    scales = np.sort(packing.diameters) / np.sort(initial_diameters)
    assert packing.permute_exchanges > 0
    assert np.ptp(scales) <= 1e-12 * scales.mean()

    # Quench, serialize the realized state, then resume to the tighter target.
    packing.quench(tol=1e-4)
    restored = rcpgenerator.Packing.from_dict(
        json.loads(json.dumps(packing.to_dict()))
    )
    restored.quench(tol=1e-10, resume=True)
    assert restored.quench_status == "converged", restored.quench_message
    assert restored.quench_residual <= 1e-10

    # Exercise the 4D neighbor search and quench path used by the public wheel.
    higher_dimensional = rcpgenerator.Packing(
        N=300,
        Ndim=4,
        phi=0.05,
        box=[1.0] * 4,
        walls=[1] * 4,
        seed=1,
        dist={
            "type": "powerlaw",
            "d_min": 1.0,
            "d_max": 2.0,
            "exponent": -3.0,
        },
    )
    higher_dimensional.pack()
    higher_dimensional.quench()
    assert higher_dimensional.quench_status == "converged", (
        higher_dimensional.quench_message
    )

    # Malformed walls must fail clearly.
    try:
        rcpgenerator.Packing(
            N=50,
            Ndim=3,
            phi=0.05,
            box=[1.0] * 3,
            walls=[0, -2, 0],
            seed=1,
        )
    except ValueError:
        pass
    else:
        raise AssertionError("malformed walls must raise ValueError")

    # A very small periodic high-dimensional cell must report its limitation.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        rcpgenerator.Packing(
            N=200,
            Ndim=6,
            phi=0.01,
            box=[1.0] * 6,
            walls=[0] * 6,
            seed=3,
        ).pack()
    assert any("periodic box is small" in str(item.message) for item in caught)

    print("rcpgenerator 1.1 feature smoke OK")


if __name__ == "__main__":
    main()
