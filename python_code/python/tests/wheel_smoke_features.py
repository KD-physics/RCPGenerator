"""Fast installed-wheel smoke test for established public rcpgenerator features."""

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
    packing.quench(tol=1e-4, kappa_mult=1.003)
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
    higher_dimensional.quench(kappa_mult=1.01)
    assert higher_dimensional.quench_status == "converged", (
        higher_dimensional.quench_message
    )

    # A longer first stage (extend), then the adaptive quench: the diameters follow
    # the pressure mu until the forces balance. It is run at ten times the pressure
    # the packing run ended on, so that the overlaps are large enough for the plain
    # quench below to reach its tolerance well above floating-point rounding.
    adaptive = rcpgenerator.Packing(
        N=300,
        Ndim=3,
        phi=0.05,
        box=[1.0] * 3,
        walls=[0] * 3,
        seed=2,
        dist={
            "type": "bidisperse",
            "d_small": 1.0,
            "d_large": 1.4,
            "fraction_small": 0.5,
        },
    )
    adaptive.pack(extend=1.5)
    assert adaptive.mu_final is not None and adaptive.mu_final > 0.0
    pressure = 10.0 * adaptive.mu_final
    packed_diameters = np.array(adaptive.diameters)
    adaptive.adaptive_quench(mu=pressure)
    assert adaptive.quench_status == "converged", adaptive.quench_message
    assert adaptive.mu_final == pressure
    scale = np.array(adaptive.diameters) / packed_diameters
    assert np.ptp(scale) <= 1e-12, "adaptive quench must scale all diameters alike"

    # A plain quench of that balanced state: it must converge without a warning
    # and must leave the diameters exactly as they are (no default expansion).
    balanced_diameters = np.array(adaptive.diameters)
    adaptive.quench(tol=1e-8)
    assert adaptive.quench_status == "converged", adaptive.quench_message
    assert adaptive.quench_residual <= 1e-8
    assert np.array_equal(np.array(adaptive.diameters), balanced_diameters)

    # Contact numbers: one integer per particle (a wall would count as one), and
    # the three summaries that follow from them.
    z = adaptive.z
    assert isinstance(z, list) and len(z) == adaptive.N
    assert all(isinstance(v, int) for v in z)
    assert min(z) >= 0 and max(z) <= 24, (min(z), max(z))
    need = adaptive.Ndim + 1
    held = [v for v in z if v >= need]
    assert abs(adaptive.z_mean - sum(z) / len(z)) <= 1e-12
    assert adaptive.n_rattlers == len(z) - len(held)
    assert held and abs(adaptive.z_mean_no_rattlers - sum(held) / len(held)) <= 1e-12
    # a force-balanced packing: the particles that are not rattlers have at least
    # 2 * Ndim contacts on average
    assert adaptive.z_mean_no_rattlers >= 2 * adaptive.Ndim, adaptive.z_mean_no_rattlers

    # Mean and largest relative overlap of the pair contacts.
    assert 0.0 < adaptive.mean_overlap <= adaptive.max_overlap < 0.01
    assert adaptive.max_overlap <= 1.0 - adaptive.shrink_factor + 1e-9

    # They are not saved; loading a packed state recomputes them from the state.
    loaded = rcpgenerator.Packing.from_dict(json.loads(json.dumps(adaptive.to_dict())))
    assert loaded.mu_final == adaptive.mu_final
    assert loaded.phi_contact == adaptive.phi_contact
    assert loaded.z == z
    assert loaded.n_rattlers == adaptive.n_rattlers
    assert loaded.z_mean == adaptive.z_mean
    assert loaded.z_mean_no_rattlers == adaptive.z_mean_no_rattlers
    assert abs(loaded.mean_overlap - adaptive.mean_overlap) <= 1e-12
    assert loaded.max_overlap == adaptive.max_overlap

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

    # After 1.2.0, second part: lr_method, pack(jam=...), pressure and stress, the Hessian check, and the build itself.
    # A wheel must not contain the developer switches (CMake option RCP_DEV_SWITCHES, off by default).
    import os
    if not os.environ.get("RCP_ALLOW_DEV_BUILD"):
        assert rcpgenerator._rcpgenerator.DEV_BUILD is False, "this build has the developer switches compiled in"
        try:
            rcpgenerator.Packing(N=10, Ndim=2, adam_alpha_mode=5)
        except TypeError:
            pass
        else:
            raise AssertionError("a developer keyword must be refused by a normal build")
    spec = dict(N=300, Ndim=3, box=[1.0] * 3, walls=[0, 0, 1], seed=3,
                dist={"type": "bidisperse", "d_small": 1.0, "d_large": 1.4, "fraction_small": 0.5})
    default = rcpgenerator.Packing(**spec)
    default.pack()
    named = rcpgenerator.Packing(lr_method="bounded_feedback", **spec)
    named.pack()
    assert named.positions == default.positions and named.diameters == default.diameters   # the default by its name
    active = rcpgenerator.Packing(lr_method="highly_active", **spec)
    active.pack(jam=5e-4)                                   # no warning suppressed: the quench must converge
    assert active.lr_method == "highly_active" and active.quench_status == "converged"
    assert active.positions != default.positions and 0.5 < active.phi_contact < 0.7
    n = 3
    assert active.pressure > 0.0 and len(active.particle_virial) == spec["N"] and min(active.particle_virial) >= 0.0
    assert abs(sum(active.particle_virial) / (n * active.stress_volume) / active.pressure - 1.0) < 1e-12
    assert abs(sum(active.stress[d][d] for d in range(n)) / n / active.pressure - 1.0) < 1e-12
    assert active.boundaries == [(2, 0), (2, 1)]            # the two walls of the third axis
    for wall_pressure in active.boundary_pressure:           # momentum balance: wall pressure = normal stress
        assert abs(wall_pressure / active.stress[2][2] - 1.0) < 1e-6
    values = active.eigenvalues()                            # built at the end of the quench; all eigenvalues
    check = active.jam_check
    assert len(values) == n * spec["N"] and check["jammed"] and check["n_negative_held"] == 0
    assert check["n_zero_held"] == check["n_zero_allowed"] == 2      # translations along the two periodic axes
    assert check["n_held"] + check["n_rattlers"] == spec["N"]
    matrix = active.hessian()
    assert matrix.shape == (n * spec["N"], n * spec["N"]) and matrix.nnz > 0
    reloaded = rcpgenerator.Packing.from_dict(active.to_dict())
    assert reloaded.lr_method == "highly_active" and abs(reloaded.pressure / active.pressure - 1.0) < 1e-9
    permuted = rcpgenerator.Packing(N=400, Ndim=3, box=[1.0] * 3, walls=[0] * 3, seed=5, lr_method="highly_active", diameter_permute=True,
                                    dist={"type": "powerlaw", "d_min": 1.0, "d_max": 3.0, "exponent": -3.0})
    permuted.pack()
    assert permuted.permute_exchanges > 0

    print("rcpgenerator feature smoke OK")


if __name__ == "__main__":
    main()
