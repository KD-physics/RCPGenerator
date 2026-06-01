from __future__ import annotations

from example_2d_circular_container import run_case as run_2d_circle
from example_2d_bidisperse_box import run_case as run_2d_bidisperse
from example_2d_polydisperse_box import run_case as run_2d_polydisperse
from example_2d_powerlaw_phi_staging import run_case as run_2d_powerlaw_staging
from example_3d_bidisperse_box import run_case as run_3d_bidisperse
from example_3d_cylindrical_container import run_case as run_3d_cylinder
from example_3d_monodisperse_box import run_case as run_3d_mono
from example_3d_spherical_container import run_case as run_3d_sphere
from minimal import run_case as run_2d_mono


def main() -> None:
    runs = [
        run_2d_mono,
        run_2d_bidisperse,
        run_2d_polydisperse,
        run_2d_circle,
        run_2d_powerlaw_staging,
        run_3d_mono,
        run_3d_bidisperse,
        run_3d_cylinder,
        run_3d_sphere,
    ]

    summaries = []
    for run in runs:
        summaries.append(run(verbose=False))
        print("-" * 60)

    print("summary:")
    for summary in summaries:
        print(
            f"{summary['case']}: phi={summary['phi']}, steps={summary['steps']}, "
            f"force={summary['force_magnitude']}, max_min_dist={summary['max_min_dist']}, "
            f"image={summary['image_path']}"
        )


if __name__ == "__main__":
    main()
