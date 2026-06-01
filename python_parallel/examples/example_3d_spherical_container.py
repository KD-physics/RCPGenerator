import rcpgenerator

from _common import build_summary
from _common import print_case_summary
from _common import render_case


CASE_NAME = "3d_spherical_container"


def run_case(verbose: bool = False) -> dict[str, object]:
    packing = rcpgenerator.Packing(
        phi=0.25,
        N=500,
        Ndim=3,
        box=[1.0, 1.0, 1.0],
        walls=[-3, 0, 0],
        fix_height=False,
        dist={"type": "mono", "d": 1.0},
        neighbor_max=0,
        seed=130,
        verbose=verbose,
    )

    packing.pack()
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
