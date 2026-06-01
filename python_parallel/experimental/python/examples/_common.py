from __future__ import annotations

from typing import Any
from pathlib import Path

import rcpgenerator


def distribution_summary(dist: dict[str, Any]) -> str:
    dist_type = dist.get("type", "mono")
    if dist_type == "mono":
        return f"mono(d={dist.get('d')})"
    if dist_type == "bidisperse":
        return (
            f"bidisperse(d1={dist.get('d1')}, d2={dist.get('d2')}, "
            f"p={dist.get('p')})"
        )
    if dist_type == "lognormal":
        return f"lognormal(mu={dist.get('mu')}, sigma={dist.get('sigma')})"
    if dist_type == "flat":
        return f"flat(d_min={dist.get('d_min')}, d_max={dist.get('d_max')})"
    if dist_type == "custom":
        values = dist.get("custom", [])
        return f"custom(len={len(values)})"
    return dist_type


def walls_summary(walls: list[int]) -> str:
    return "periodic" if all(int(w) == 0 for w in walls) else str(walls)


def print_case_summary(case_name: str, packing: rcpgenerator.Packing) -> None:
    print(f"case: {case_name}")
    print(f"N: {packing.N}")
    print(f"Ndim: {packing.Ndim}")
    print(f"box: {list(packing.box)}")
    print(f"walls: {walls_summary(list(packing.walls))}")
    print(f"distribution: {distribution_summary(dict(packing.dist))}")
    print(f"final phi: {packing.phi_final}")
    print(f"steps: {packing.steps}")
    print(f"force magnitude: {packing.force_magnitude}")
    print(f"max_min_dist: {packing.max_min_dist}")


def render_case(case_name: str, packing: rcpgenerator.Packing) -> str:
    output_dir = Path(__file__).resolve().parent / "output"
    output_path = output_dir / f"{case_name}.png"
    palette_map = {
        "2d_monodisperse_box": 1,
        "2d_bidisperse_box": 5,
        "2d_polydisperse_box": 10,
        "2d_circular_container": 12,
        "2d_powerlaw_phi_staging": 8,
        "3d_monodisperse_box": 2,
        "3d_bidisperse_box": 6,
        "3d_cylindrical_container": 9,
        "3d_spherical_container": 4,
        "starter": 7,
    }
    return packing.savefig(
        path=output_path,
        palette_choice=palette_map.get(case_name, 1),
    )


def build_summary(case_name: str, packing: rcpgenerator.Packing) -> dict[str, Any]:
    return {
        "case": case_name,
        "N": packing.N,
        "Ndim": packing.Ndim,
        "box": list(packing.box),
        "walls": walls_summary(list(packing.walls)),
        "distribution": distribution_summary(dict(packing.dist)),
        "phi": packing.phi_final,
        "steps": packing.steps,
        "force_magnitude": packing.force_magnitude,
        "max_min_dist": packing.max_min_dist,
    }
