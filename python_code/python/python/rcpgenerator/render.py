"""Python-side rendering helpers for rcpgenerator packings."""

from __future__ import annotations

from pathlib import Path
import os
from typing import Iterable

import matplotlib


def _should_force_agg() -> bool:
    try:
        from IPython import get_ipython
    except Exception:
        get_ipython = None

    if get_ipython is not None and get_ipython() is not None:
        return False
    if os.name == "nt":
        return False
    return not bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))


if _should_force_agg():
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

_PALETTES = {
    1: np.array([[204, 85, 85], [255, 179, 71], [255, 255, 128], [119, 221, 119], [127, 106, 217], [0, 225, 225]], dtype=float) / 255.0,
    2: np.array([[0, 128, 128], [236, 224, 200], [112, 128, 144], [240, 128, 128]], dtype=float) / 255.0,
    3: np.array([[170, 189, 145], [222, 210, 186], [198, 134, 103], [176, 182, 186]], dtype=float) / 255.0,
    4: np.array([[143, 170, 180], [244, 221, 220], [196, 177, 162], [60, 76, 89]], dtype=float) / 255.0,
    5: np.array([[255, 99, 71], [255, 165, 0], [255, 223, 0], [192, 241, 192], [255, 250, 205], [50, 205, 50]], dtype=float) / 255.0,
    6: np.array([[0, 191, 255], [30, 124, 245], [255, 105, 180], [235, 215, 230], [220, 235, 235], [255, 69, 0]], dtype=float) / 255.0,
    7: np.array([[255, 110, 255], [0, 255, 127], [173, 255, 47], [255, 235, 250], [155, 80, 210]], dtype=float) / 255.0,
    8: np.array([[255, 20, 147], [0, 255, 255], [238, 130, 238], [124, 252, 0]], dtype=float) / 255.0,
    9: np.array([[176, 224, 230], [173, 216, 230], [255, 250, 205], [234, 154, 86], [240, 255, 255], [70, 130, 180]], dtype=float) / 255.0,
    10: np.array([[255, 162, 173], [142, 251, 142], [212, 251, 212], [173, 216, 230], [255, 235, 250], [90, 150, 200]], dtype=float) / 255.0,
    11: np.array([[255, 165, 0], [255, 99, 71], [186, 184, 108], [255, 250, 205], [135, 206, 235]], dtype=float) / 255.0,
    12: np.array([[205, 92, 92], [213, 181, 110], [156, 111, 68], [210, 105, 30], [234, 154, 86], [245, 218, 171]], dtype=float) / 255.0,
}


def _particle_colors(count: int, palette_choice: int) -> np.ndarray:
    rng = np.random.default_rng(42)
    palette = _PALETTES.get(palette_choice, _PALETTES[1])
    idx = rng.integers(0, len(palette), size=count)
    return palette[idx]


def _sphere_facecolors(base_color: np.ndarray, sx: np.ndarray, sy: np.ndarray, sz: np.ndarray) -> np.ndarray:
    normals = np.stack((sx, sy, sz), axis=-1)
    light_dir = np.array([-0.3, -0.18, 0.94], dtype=float)
    light_dir /= np.linalg.norm(light_dir)
    view_dir = np.array([0.18, -0.12, 0.98], dtype=float)
    view_dir /= np.linalg.norm(view_dir)
    half_vec = light_dir + view_dir
    half_vec /= np.linalg.norm(half_vec)

    diffuse = np.clip(np.tensordot(normals, light_dir, axes=([-1], [0])), 0.0, 1.0)
    specular = np.clip(np.tensordot(normals, half_vec, axes=([-1], [0])), 0.0, 1.0) ** 28
    ambient = 0.5
    intensity = ambient + 0.48 * diffuse
    lifted = 0.7 * base_color + 0.3
    rgb = np.clip(lifted[None, None, :] * intensity[..., None] + 0.3 * specular[..., None], 0.0, 1.0)
    alpha = np.ones((*rgb.shape[:2], 1), dtype=float)
    return np.concatenate((rgb, alpha), axis=-1)


def _normalized_shape(walls: Iterable[int], ndim: int) -> tuple[int, ...]:
    values = list(int(v) for v in walls)
    if not values:
        return tuple(0 for _ in range(ndim))
    if values[0] < 0:
        k = min(-values[0], ndim)
        normalized = values[:]
        if len(normalized) < ndim:
            normalized.extend([0] * (ndim - len(normalized)))
        for i in range(1, k):
            normalized[i] = -1
        return tuple(normalized[:ndim])
    if len(values) < ndim:
        values.extend([0] * (ndim - len(values)))
    return tuple(values[:ndim])


def _container_kind(box: list[float], walls: list[int]) -> str:
    ndim = len(box)
    normalized = _normalized_shape(walls, ndim)
    if ndim == 2 and normalized[0] < 0:
        return "circle"
    if ndim == 3 and normalized[0] == -2:
        return "cylinder"
    if ndim == 3 and normalized[0] == -3:
        return "sphere"
    return "box"


def _draw_2d_boundary(ax: plt.Axes, box: list[float], walls: list[int]) -> None:
    kind = _container_kind(box, walls)
    if kind == "circle":
        ax.add_patch(Circle((box[0] / 2.0, box[0] / 2.0), box[0] / 2.0, fill=False, ec="black", lw=2.0))
    else:
        x = [0, box[0], box[0], 0, 0]
        y = [0, 0, box[1], box[1], 0]
        ax.plot(x, y, color="black", linewidth=2.0)


def _periodic_shifts_2d(center: np.ndarray, radius: float, box: list[float], walls: list[int]) -> list[tuple[float, float]]:
    shifts_x = [0.0]
    shifts_y = [0.0]

    if int(walls[0]) == 0:
        if center[0] - radius < 0.0:
            shifts_x.append(box[0])
        if center[0] + radius > box[0]:
            shifts_x.append(-box[0])

    if int(walls[1]) == 0:
        if center[1] - radius < 0.0:
            shifts_y.append(box[1])
        if center[1] + radius > box[1]:
            shifts_y.append(-box[1])

    shifts = []
    for dx in shifts_x:
        for dy in shifts_y:
            shifts.append((dx, dy))
    return shifts


def _draw_box_wireframe(ax: plt.Axes, box: list[float]) -> None:
    lx, ly, lz = box
    corners = np.array(
        [
            [0, 0, 0],
            [lx, 0, 0],
            [lx, ly, 0],
            [0, ly, 0],
            [0, 0, lz],
            [lx, 0, lz],
            [lx, ly, lz],
            [0, ly, lz],
        ]
    )
    edges = [
        (0, 1), (1, 2), (2, 3), (3, 0),
        (4, 5), (5, 6), (6, 7), (7, 4),
        (0, 4), (1, 5), (2, 6), (3, 7),
    ]
    for a, b in edges:
        pts = corners[[a, b]]
        ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], color="black", linewidth=1.0, alpha=0.8)


def _draw_cylinder_wireframe(ax: plt.Axes, box: list[float]) -> None:
    theta = np.linspace(0.0, 2.0 * np.pi, 100)
    radius = box[0] / 2.0
    center = radius
    z0, z1 = 0.0, box[2]
    x = center + radius * np.cos(theta)
    y = center + radius * np.sin(theta)
    ax.plot(x, y, np.full_like(x, z0), color="black", linewidth=1.0)
    ax.plot(x, y, np.full_like(x, z1), color="black", linewidth=1.0)
    for angle in np.linspace(0.0, 2.0 * np.pi, 8, endpoint=False):
        xc = center + radius * np.cos(angle)
        yc = center + radius * np.sin(angle)
        ax.plot([xc, xc], [yc, yc], [z0, z1], color="black", linewidth=0.8, alpha=0.6)


def _draw_sphere_wireframe(ax: plt.Axes, box: list[float]) -> None:
    radius = box[0] / 2.0
    center = np.array([radius, radius, radius])
    u = np.linspace(0.0, 2.0 * np.pi, 40)
    v = np.linspace(0.0, np.pi, 20)
    x = center[0] + radius * np.outer(np.cos(u), np.sin(v))
    y = center[1] + radius * np.outer(np.sin(u), np.sin(v))
    z = center[2] + radius * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_wireframe(x, y, z, color="black", linewidth=0.4, alpha=0.35, rstride=2, cstride=2)


def _draw_3d_boundary(ax: plt.Axes, box: list[float], walls: list[int]) -> None:
    kind = _container_kind(box, walls)
    if kind == "cylinder":
        _draw_cylinder_wireframe(ax, box)
    elif kind == "sphere":
        _draw_sphere_wireframe(ax, box)
    else:
        _draw_box_wireframe(ax, box)


def _plot_2d(packing, palette_choice: int, figsize=(5, 5), dpi=120) -> plt.Figure:
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    positions = np.asarray(packing.positions, dtype=float)
    diameters = np.asarray(packing.diameters, dtype=float)
    colors = _particle_colors(len(diameters), palette_choice)

    for i in range(len(diameters)):
        radius = diameters[i] / 2.0
        for dx, dy in _periodic_shifts_2d(positions[i], radius, list(packing.box), list(packing.walls)):
            circle = Circle(
                (positions[i, 0] + dx, positions[i, 1] + dy),
                radius,
                facecolor=colors[i],
                edgecolor="black",
                linewidth=0.6,
            )
            ax.add_patch(circle)

    _draw_2d_boundary(ax, list(packing.box), list(packing.walls))
    ax.set_xlim(0.0, packing.box[0])
    ax.set_ylim(0.0, packing.box[1])
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    fig.tight_layout(pad=0.05)
    return fig


def _plot_3d(packing, palette_choice: int, figsize=(5, 5), dpi=120) -> plt.Figure:
    fig = plt.figure(figsize=figsize, dpi=dpi)    
    ax = fig.add_subplot(111, projection="3d")
    positions = np.asarray(packing.positions, dtype=float)
    diameters = np.asarray(packing.diameters, dtype=float)
    colors = _particle_colors(len(diameters), palette_choice)
    u = np.linspace(0.0, 2.0 * np.pi, 36)
    v = np.linspace(0.0, np.pi, 24)
    sx = np.outer(np.cos(u), np.sin(v))
    sy = np.outer(np.sin(u), np.sin(v))
    sz = np.outer(np.ones_like(u), np.cos(v))

    order = np.argsort(positions[:, 2])
    for idx in order:
        radius = diameters[idx] / 2.0
        x = radius * sx + positions[idx, 0]
        y = radius * sy + positions[idx, 1]
        z = radius * sz + positions[idx, 2]
        facecolors = _sphere_facecolors(colors[idx], sx, sy, sz)
        ax.plot_surface(
            x,
            y,
            z,
            rstride=1,
            cstride=1,
            facecolors=facecolors,
            edgecolor="none",
            linewidth=0.0,
            shade=False,
            antialiased=True,
        )

    _draw_3d_boundary(ax, list(packing.box), list(packing.walls))
    ax.set_xlim(0.0, packing.box[0])
    ax.set_ylim(0.0, packing.box[1])
    ax.set_zlim(0.0, packing.box[2])
    ax.set_box_aspect(tuple(packing.box))
    ax.view_init(elev=22, azim=35)
    ax.set_axis_off()
    ax.set_facecolor("white")
    fig.tight_layout(pad=0.05)
    return fig


def render_packing(
    packing,
    path: str | Path | None = None,
    show: bool = False,
    palette_choice: int = 1,
    figsize: tuple[float, float] = (4, 4),
    dpi: int = 110,
) -> str | None:
    output_path = Path(path) if path is not None else None
    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)

    if packing.Ndim == 2:
        fig = _plot_2d(packing, palette_choice, figsize=figsize, dpi=dpi)
    elif packing.Ndim == 3:
        fig = _plot_3d(packing, palette_choice, figsize=figsize, dpi=dpi)
    else:
        raise ValueError(f"Unsupported rendering dimension: {packing.Ndim}")

    if output_path is not None:
        fig.savefig(output_path, bbox_inches="tight", facecolor="white")
    if show:
        plt.show()
    plt.close(fig)
    return str(output_path) if output_path is not None else None

def animate_packing_2d(
    packing,
    path: str | Path | None = None,
    show: bool = True,
    palette_choice: int = 1,
    interval_ms: int = 80,
    repeat: bool = False,
):
    if packing.Ndim != 2:
        raise ValueError("animate_packing_2d only supports 2D packings")
    if not getattr(packing, "trajectory_positions", None):
        raise ValueError("No trajectory_positions recorded; call pack(capture_trajectory=True, ...) first")
    if not getattr(packing, "trajectory_diameters", None):
        raise ValueError("No trajectory_diameters recorded; call pack(capture_trajectory=True, ...) first")

    positions_frames = [np.asarray(frame, dtype=float) for frame in packing.trajectory_positions]
    diameters_frames = [np.asarray(frame, dtype=float) for frame in packing.trajectory_diameters]
    colors = _particle_colors(len(diameters_frames[0]), palette_choice)

    fig, ax = plt.subplots(figsize=(8, 8), dpi=180)
    ax.set_xlim(0.0, packing.box[0])
    ax.set_ylim(0.0, packing.box[1])
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    _draw_2d_boundary(ax, list(packing.box), list(packing.walls))

    circles = []
    for i in range(len(diameters_frames[0])):
        circle = Circle(
            (positions_frames[0][i, 0], positions_frames[0][i, 1]),
            diameters_frames[0][i] / 2.0,
            facecolor=colors[i],
            edgecolor="black",
            linewidth=0.5,
        )
        ax.add_patch(circle)
        circles.append(circle)

    def update(frame_index: int):
        positions = positions_frames[frame_index]
        diameters = diameters_frames[frame_index]
        for idx, circle in enumerate(circles):
            circle.center = (positions[idx, 0], positions[idx, 1])
            circle.radius = diameters[idx] / 2.0
        title = ""
        if getattr(packing, "trajectory_steps", None):
            step = packing.trajectory_steps[frame_index]
            phi = packing.trajectory_phi[frame_index] if frame_index < len(packing.trajectory_phi) else None
            title = f"step={step}"
            if phi is not None:
                title += f", phi={phi:.6f}"
        ax.set_title(title)
        return circles

    anim = animation.FuncAnimation(
        fig,
        update,
        frames=len(positions_frames),
        interval=interval_ms,
        blit=False,
        repeat=repeat,
    )

    output_path = Path(path) if path is not None else None
    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        if output_path.suffix.lower() != ".gif":
            raise ValueError("Only GIF saving is currently supported for animate_packing_2d")
        anim.save(output_path, writer=animation.PillowWriter(fps=max(1, int(round(1000 / interval_ms)))))

    if show:
        try:
            from IPython.display import HTML, display

            display(HTML(anim.to_jshtml()))
        except Exception:
            plt.show()
    plt.close(fig)
    return anim
