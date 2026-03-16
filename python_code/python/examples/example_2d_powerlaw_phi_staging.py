import rcpgenerator

from _common import build_summary
from _common import print_case_summary
from _common import render_case


CASE_NAME = "2d_powerlaw_phi_staging"


def run_case(verbose: bool = False) -> dict[str, object]:
    packing = rcpgenerator.Packing(
        phi=0.11,
        N=250,
        Ndim=2,
        box=[1.0, 1.0],
        walls=[0, 0],
        fix_height=False,
        dist={"type": "powerlaw", "d_min": 0.3, "d_max": 1.8, "exponent": -2.5},
        neighbor_max=0,
        seed=131,
        verbose=verbose,
    )

    packing.relax(n_steps=5000, target_phi=0.82, verbose=verbose)
    for target_phi in [0.84, 0.86, 0.88, 0.90, 0.92, 0.94]:
        packing.update_phi(target_phi)
        packing.relax(n_steps=1000, fix_diameter=True, verbose=verbose)

    print_case_summary(CASE_NAME, packing)
    image_path = render_case(CASE_NAME, packing)
    print(f"image: {image_path}")
    summary = build_summary(CASE_NAME, packing)
    summary["image_path"] = image_path
    return summary


def main() -> None:
    run_case(verbose=True)


if __name__ == "__main__":
    main()
