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
from matplotlib.patches import Circle, Ellipse
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


# ---- ellipsoids (additive). A packing with aspect ratios that are not all 1 is drawn by the functions below;
# everything else (spheres, and ellipsoids whose ratios are all exactly 1) goes through the sphere functions
# above, untouched, so those pictures are pixel-identical to what they always were.
def _ellipsoid_ratios(packing) -> list[float] | None:
    """The packing's aspect ratios if it is to be drawn as ellipsoids, else None."""
    ratios = getattr(packing, "aspect_ratios", None)
    if not ratios:
        return None
    if isinstance(ratios[0], (list, tuple)):   # a shape per particle: N rows
        table = np.asarray(ratios, dtype=float)
        return None if np.all(table == 1.0) else table
    ratios = [float(a) for a in ratios]
    if all(a == 1.0 for a in ratios):
        return None
    return ratios


def _axes_rel(ratios, n: int) -> np.ndarray:
    """Axis lengths in units of D for every particle, N x d: (1, ratios...) from one list or from a table of N rows."""
    r = np.asarray(ratios, dtype=float)
    r = np.broadcast_to(r, (n, r.shape[-1])) if r.ndim == 1 else r
    return np.concatenate((np.ones((n, 1)), r), axis=1)


def _shape_matrices(diameters: np.ndarray, orientations: np.ndarray, ratios: list[float]) -> np.ndarray:
    """A_i = R_i diag(semi-axes^2) R_i^T (N x d x d); semi-axis k = D/2 x (1, ratios...)[k]."""
    semi = 0.5 * diameters[:, None] * _axes_rel(ratios, len(diameters))
    return np.einsum("iak,ik,ibk->iab", orientations, semi ** 2, orientations)


def _ellipse_params_2d(diameters: np.ndarray, orientations: np.ndarray, ratios: list[float]):
    """For matplotlib's Ellipse: width (full length along body axis 1 = D), height (D x ratio), angle in degrees
    of body axis 1 from the x axis, and the half-extents along x and y (sqrt of the diagonal of A)."""
    width = diameters
    height = diameters * _axes_rel(ratios, len(diameters))[:, 1]
    angle = np.degrees(np.arctan2(orientations[:, 1, 0], orientations[:, 0, 0]))
    A = _shape_matrices(diameters, orientations, ratios)
    return width, height, angle, np.sqrt(A[:, 0, 0]), np.sqrt(A[:, 1, 1])


def _periodic_shifts_2d_extents(center, hx: float, hy: float, box: list[float], walls: list[int]) -> list[tuple[float, float]]:
    """Periodic images needed for a shape with half-extents hx, hy (same rule as _periodic_shifts_2d)."""
    shifts_x = [0.0]
    shifts_y = [0.0]
    if int(walls[0]) == 0:
        if center[0] - hx < 0.0:
            shifts_x.append(box[0])
        if center[0] + hx > box[0]:
            shifts_x.append(-box[0])
    if int(walls[1]) == 0:
        if center[1] - hy < 0.0:
            shifts_y.append(box[1])
        if center[1] + hy > box[1]:
            shifts_y.append(-box[1])
    return [(dx, dy) for dx in shifts_x for dy in shifts_y]


def _plot_2d_ellipses(packing, ratios: list[float], palette_choice: int, figsize=(5, 5), dpi=120) -> plt.Figure:
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    positions = np.asarray(packing.positions, dtype=float)
    diameters = np.asarray(packing.diameters, dtype=float)
    orientations = np.asarray(packing.orientations, dtype=float)
    colors = _particle_colors(len(diameters), palette_choice)
    width, height, angle, hx, hy = _ellipse_params_2d(diameters, orientations, ratios)

    for i in range(len(diameters)):
        for dx, dy in _periodic_shifts_2d_extents(positions[i], hx[i], hy[i], list(packing.box), list(packing.walls)):
            ax.add_patch(Ellipse(
                (positions[i, 0] + dx, positions[i, 1] + dy),
                width[i],
                height[i],
                angle=angle[i],
                facecolor=colors[i],
                edgecolor="black",
                linewidth=0.6,
            ))

    _draw_2d_boundary(ax, list(packing.box), list(packing.walls))
    ax.set_xlim(0.0, packing.box[0])
    ax.set_ylim(0.0, packing.box[1])
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    fig.tight_layout(pad=0.05)
    return fig


def _ellipsoid_surface(center: np.ndarray, rotation: np.ndarray, semi: np.ndarray,
                       sx: np.ndarray, sy: np.ndarray, sz: np.ndarray):
    """Surface points of one ellipsoid and its outward unit normals, from the unit-sphere mesh (sx, sy, sz):
    point = centre + R (semi * u), normal ~ R (u / semi)."""
    unit = np.stack((sx, sy, sz), axis=-1)
    pts = (unit * semi) @ rotation.T + center
    nrm = (unit / semi) @ rotation.T
    nrm /= np.linalg.norm(nrm, axis=-1, keepdims=True)
    return pts[..., 0], pts[..., 1], pts[..., 2], nrm[..., 0], nrm[..., 1], nrm[..., 2]


def _plot_3d_ellipsoids(packing, ratios: list[float], palette_choice: int, figsize=(5, 5), dpi=120) -> plt.Figure:
    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111, projection="3d")
    positions = np.asarray(packing.positions, dtype=float)
    diameters = np.asarray(packing.diameters, dtype=float)
    orientations = np.asarray(packing.orientations, dtype=float)
    colors = _particle_colors(len(diameters), palette_choice)
    u = np.linspace(0.0, 2.0 * np.pi, 36)
    v = np.linspace(0.0, np.pi, 24)
    sx = np.outer(np.cos(u), np.sin(v))
    sy = np.outer(np.sin(u), np.sin(v))
    sz = np.outer(np.ones_like(u), np.cos(v))
    axes_rel = _axes_rel(ratios, len(diameters))

    order = np.argsort(positions[:, 2])
    for idx in order:
        semi = 0.5 * diameters[idx] * axes_rel[idx]
        x, y, z, nx, ny, nz = _ellipsoid_surface(positions[idx], orientations[idx], semi, sx, sy, sz)
        facecolors = _sphere_facecolors(colors[idx], nx, ny, nz)
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


# ---- rods (additive). A packing with rod_length / rod_radius is drawn by the functions below. Conventions
# (BUILD.md section 9): both numbers are in units of the particle's D, so the straight section is D x rod_length
# long and the radius is D x rod_radius; the rod's axis is the first body axis = column 0 of the orientation
# matrix; rod_capped = half-sphere ends (2D: a stadium), otherwise flat ends (2D: a rectangle, 3D: a cylinder).
_ROD_ARC = 64      # polygon segments per half-circle end of a 2D stadium
_ROD_AROUND = 24   # 3D: facets around the rod's axis
_ROD_CAP = 6       # 3D: rings per hemispherical end


def _rod_spec(packing):
    """(half_length, radius, capped) of every particle in real units (two arrays of N and a flag) if the packing
    is of rods, else None. rod_length and rod_radius may each be one number or a list of N."""
    length = getattr(packing, "rod_length", None)
    radius = getattr(packing, "rod_radius", None)
    if length is None and radius is None:
        return None
    if length is None or radius is None:
        raise ValueError("rods need both rod_length and rod_radius to be drawn")
    diameters = np.asarray(packing.diameters, dtype=float)
    n = len(diameters)
    length = np.broadcast_to(np.asarray(length, dtype=float), (n,))
    radius = np.broadcast_to(np.asarray(radius, dtype=float), (n,))
    orientations = getattr(packing, "orientations", None)
    if orientations is None or len(orientations) != n:
        raise ValueError("rods need one orientation matrix per particle to be drawn (the rod's axis is its first column)")
    return 0.5 * diameters * length, diameters * radius, bool(getattr(packing, "rod_capped", True))


def _rod_outline_2d(center, axis, half_length: float, radius: float, capped: bool) -> np.ndarray:
    """Vertices (counter-clockwise) of one 2D rod: a stadium (two half-circles of `radius` centred at
    centre +- half_length x axis, joined by straight sides) or the rectangle 2 half_length by 2 radius."""
    c = np.asarray(center, dtype=float)
    u = np.asarray(axis, dtype=float)
    v = np.array([-u[1], u[0]])
    if not capped:
        return np.array([c + half_length * u - radius * v, c + half_length * u + radius * v,
                         c - half_length * u + radius * v, c - half_length * u - radius * v])
    a = np.linspace(-0.5 * np.pi, 0.5 * np.pi, _ROD_ARC + 1)
    arc = radius * (np.cos(a)[:, None] * u + np.sin(a)[:, None] * v)
    return np.concatenate((c + half_length * u + arc, c - half_length * u - arc))


def _rod_half_extents_2d(axis, half_length, radius, capped: bool):
    """Half-widths of a 2D rod's bounding box along x and y (for the periodic images)."""
    ux, uy = abs(axis[0]), abs(axis[1])
    if capped:
        return half_length * ux + radius, half_length * uy + radius
    return half_length * ux + radius * uy, half_length * uy + radius * ux


def _plot_2d_rods(packing, rods, palette_choice: int, figsize=(5, 5), dpi=120) -> plt.Figure:
    from matplotlib.patches import Polygon
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    positions = np.asarray(packing.positions, dtype=float)
    orientations = np.asarray(packing.orientations, dtype=float)
    half_length, radius, capped = rods
    colors = _particle_colors(len(half_length), palette_choice)

    for i in range(len(half_length)):
        axis = orientations[i][:, 0]
        hx, hy = _rod_half_extents_2d(axis, half_length[i], radius[i], capped)
        for dx, dy in _periodic_shifts_2d_extents(positions[i], hx, hy, list(packing.box), list(packing.walls)):
            ax.add_patch(Polygon(
                _rod_outline_2d(positions[i] + np.array([dx, dy]), axis, half_length[i], radius[i], capped),
                closed=True,
                facecolor=colors[i],
                edgecolor="black",
                linewidth=0.6,
            ))

    _draw_2d_boundary(ax, list(packing.box), list(packing.walls))
    ax.set_xlim(0.0, packing.box[0])
    ax.set_ylim(0.0, packing.box[1])
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    fig.tight_layout(pad=0.05)
    return fig


def _rod_profile(half_length: float, radius: float, capped: bool):
    """A rod as a surface of revolution about its axis: points (t, r) (t along the axis, r from it) from one end
    to the other, with the outward normal at each as (axial, radial) components. Flat ends meet the side at a
    sharp edge, which is two profile points at the same place with different normals."""
    n_len = int(min(12, max(1, np.ceil(2.0 * half_length / max(radius, 1e-300)))))
    side_t = np.linspace(-half_length, half_length, n_len + 1)
    if capped:
        a = np.linspace(0.0, 0.5 * np.pi, _ROD_CAP + 1)
        t = np.concatenate((-half_length - radius * np.cos(a), side_t[1:], half_length + radius * np.cos(a[::-1][1:])))
        r = np.concatenate((radius * np.sin(a), np.full(n_len, radius), radius * np.sin(a[::-1][1:])))
        nt = np.concatenate((-np.cos(a), np.zeros(n_len), np.cos(a[::-1][1:])))
        nr = np.concatenate((np.sin(a), np.ones(n_len), np.sin(a[::-1][1:])))
    else:
        # the end faces take the axial normal; the side takes the radial one (the profile point that starts each
        # segment gives that segment its normal, see _rod_quads)
        t = np.concatenate(([-half_length, -half_length], side_t, [half_length, half_length]))
        r = np.concatenate(([0.0, radius], np.full(n_len + 1, radius), [radius, 0.0]))
        nt = np.concatenate(([-1.0, -1.0], np.zeros(n_len), [1.0, 1.0, 1.0]))
        nr = np.concatenate(([0.0, 0.0], np.ones(n_len), [0.0, 0.0, 0.0]))
    return t, r, nt, nr


def _rod_quads(center: np.ndarray, rotation: np.ndarray, half_length: float, radius: float, capped: bool):
    """Surface of one 3D rod as quadrilaterals (M x 4 x 3) and one outward unit normal per quadrilateral (M x 3).
    Every vertex lies on the rod: on the side at distance `radius` from the axis, on a half-sphere end at
    distance `radius` from the end point of the axis segment, or on a flat end at |t| = half_length."""
    t, r, nt, nr = _rod_profile(half_length, radius, capped)
    phi = np.linspace(0.0, 2.0 * np.pi, _ROD_AROUND + 1)
    u, e2, e3 = rotation[:, 0], rotation[:, 1], rotation[:, 2]
    ring = np.cos(phi)[:, None] * e2 + np.sin(phi)[:, None] * e3                    # unit radial directions
    pts = center + t[:, None, None] * u + r[:, None, None] * ring[None, :, :]       # (profile, around, 3)
    keep = (np.abs(np.diff(t)) + np.abs(np.diff(r))) > 0.0                          # drop zero-length profile steps
    quads = np.stack((pts[:-1, :-1], pts[1:, :-1], pts[1:, 1:], pts[:-1, 1:]), axis=2)[keep]
    if capped:   # smooth surface: the normal at the middle of each profile step
        mt, mr = 0.5 * (nt[:-1] + nt[1:]), 0.5 * (nr[:-1] + nr[1:])
    else:        # faceted along the profile: each step has one normal (end face or side)
        mt, mr = nt[:-1], nr[:-1]
    mid = 0.5 * (ring[:-1] + ring[1:])
    mid /= np.linalg.norm(mid, axis=1, keepdims=True)
    normals = mt[:, None, None] * u + mr[:, None, None] * mid[None, :, :]
    normals = normals[keep]
    normals /= np.linalg.norm(normals, axis=-1, keepdims=True)
    return quads.reshape(-1, 4, 3), normals.reshape(-1, 3)


def _plot_3d_rods(packing, rods, palette_choice: int, figsize=(5, 5), dpi=120) -> plt.Figure:
    """3D rods: cylinders with half-sphere or flat ends. All rods go into ONE polygon collection, so that the
    depth sort is per facet and not per particle (long rods pass in front of and behind each other). No
    periodic images, as for spheres and ellipsoids."""
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111, projection="3d")
    positions = np.asarray(packing.positions, dtype=float)
    orientations = np.asarray(packing.orientations, dtype=float)
    half_length, radius, capped = rods
    colors = _particle_colors(len(half_length), palette_choice)

    all_quads, all_colors = [], []
    for idx in range(len(half_length)):
        quads, normals = _rod_quads(positions[idx], orientations[idx], half_length[idx], radius[idx], capped)
        shaded = _sphere_facecolors(colors[idx], normals[:, None, 0], normals[:, None, 1], normals[:, None, 2])
        all_quads.append(quads)
        all_colors.append(shaded[:, 0, :])
    if all_quads:
        facecolors = np.concatenate(all_colors)
        ax.add_collection3d(Poly3DCollection(np.concatenate(all_quads), facecolors=facecolors, edgecolors=facecolors,
                                             linewidths=0.2))

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

    ratios = _ellipsoid_ratios(packing)
    rods = _rod_spec(packing)
    if rods is not None and packing.Ndim == 2:
        fig = _plot_2d_rods(packing, rods, palette_choice, figsize=figsize, dpi=dpi)
    elif rods is not None and packing.Ndim == 3:
        fig = _plot_3d_rods(packing, rods, palette_choice, figsize=figsize, dpi=dpi)
    elif packing.Ndim == 2:
        fig = (_plot_2d_ellipses(packing, ratios, palette_choice, figsize=figsize, dpi=dpi) if ratios is not None
               else _plot_2d(packing, palette_choice, figsize=figsize, dpi=dpi))
    elif packing.Ndim == 3:
        fig = (_plot_3d_ellipsoids(packing, ratios, palette_choice, figsize=figsize, dpi=dpi) if ratios is not None
               else _plot_3d(packing, palette_choice, figsize=figsize, dpi=dpi))
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
    is_rod = getattr(packing, "rod_length", None) is not None or getattr(packing, "rod_radius", None) is not None
    if not getattr(packing, "trajectory_positions", None):
        raise ValueError("No trajectory_positions recorded; call pack(capture_trajectory=True, ...) first")
    if not getattr(packing, "trajectory_diameters", None):
        raise ValueError("No trajectory_diameters recorded; call pack(capture_trajectory=True, ...) first")

    positions_frames = [np.asarray(frame, dtype=float) for frame in packing.trajectory_positions]
    diameters_frames = [np.asarray(frame, dtype=float) for frame in packing.trajectory_diameters]
    # The core renumbers particles while it runs, so slot k of one frame is not slot k of the next. With the
    # recorded identities (trajectory_origin) every frame is put into ONE order, that of the last frame: patch k
    # is then the same particle throughout (its colour stays, and its shape is the shape slot k has at the end).
    origin_frames = getattr(packing, "trajectory_origin", None) or []
    slot_frames = None
    if len(origin_frames) == len(positions_frames) and len(origin_frames) > 0:
        final = np.asarray(origin_frames[-1], dtype=int)
        slot_frames = []
        for o in origin_frames:
            where = np.empty(len(final), dtype=int)
            where[np.asarray(o, dtype=int)] = np.arange(len(final))      # identity -> slot in this frame
            slot_frames.append(where[final])                             # slot in this frame of the particle in final slot k
        positions_frames = [p[s] for p, s in zip(positions_frames, slot_frames)]
        diameters_frames = [d[s] for d, s in zip(diameters_frames, slot_frames)]
    colors = _particle_colors(len(diameters_frames[0]), palette_choice)

    fig, ax = plt.subplots(figsize=(8, 8), dpi=180)
    ax.set_xlim(0.0, packing.box[0])
    ax.set_ylim(0.0, packing.box[1])
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    _draw_2d_boundary(ax, list(packing.box), list(packing.walls))

    n_particles = len(diameters_frames[0])
    ratios = None if is_rod else _ellipsoid_ratios(packing)
    orientation_frames = None
    # A frame records position, D and orientation of every slot, not its shape. With one shape for all particles
    # that is enough; with a shape per particle it is enough as long as no exchange moved shapes between slots.
    per_particle = (ratios is not None and np.ndim(ratios) == 2) or (is_rod and (
        isinstance(getattr(packing, "rod_length", None), (list, tuple)) or isinstance(getattr(packing, "rod_radius", None), (list, tuple))))
    if per_particle and getattr(packing, "diameter_permute", False):
        raise NotImplementedError("animate_2d with a shape per particle AND diameter_permute is not supported: "
                                  "exchanges move shapes between particles and the trajectory does not record "
                                  "which shape sat where at each frame")
    height_rel = None
    if ratios is not None:   # ellipses: height / D of every particle (one shape for all, or a table of N rows)
        height_rel = _axes_rel(ratios, n_particles)[:, 1]
    rod_rel = None
    if is_rod:               # rods: half-length / D and radius / D of every particle, and the kind
        if getattr(packing, "rod_length", None) is None or getattr(packing, "rod_radius", None) is None:
            raise ValueError("rods need both rod_length and rod_radius to be drawn")
        rod_rel = (0.5 * np.broadcast_to(np.asarray(packing.rod_length, dtype=float), (n_particles,)),
                   np.broadcast_to(np.asarray(packing.rod_radius, dtype=float), (n_particles,)),
                   bool(getattr(packing, "rod_capped", True)))
    if height_rel is not None or rod_rel is not None:   # one patch per particle, turned frame by frame (no periodic images, as for discs)
        if not getattr(packing, "trajectory_orientations", None):
            raise ValueError("No trajectory_orientations recorded; call pack(capture_trajectory=True, ...) first")
        orientation_frames = [np.asarray(frame, dtype=float) for frame in packing.trajectory_orientations]
        if slot_frames is not None:
            orientation_frames = [o[s] for o, s in zip(orientation_frames, slot_frames)]
        elif per_particle:
            raise ValueError("this trajectory does not record which particle is in which slot (it was captured by an "
                             "older version); pack again with capture_trajectory=True to animate a shape per particle")

    def rod_outline(frame_index: int, i: int) -> np.ndarray:
        d = diameters_frames[frame_index][i]
        return _rod_outline_2d(positions_frames[frame_index][i], orientation_frames[frame_index][i][:, 0],
                               d * rod_rel[0][i], d * rod_rel[1][i], rod_rel[2])

    circles = []
    for i in range(n_particles):
        if rod_rel is not None:
            from matplotlib.patches import Polygon
            circle = Polygon(rod_outline(0, i), closed=True, facecolor=colors[i], edgecolor="black", linewidth=0.5)
        elif height_rel is not None:
            circle = Ellipse(
                (positions_frames[0][i, 0], positions_frames[0][i, 1]),
                diameters_frames[0][i],
                diameters_frames[0][i] * height_rel[i],
                angle=float(np.degrees(np.arctan2(orientation_frames[0][i, 1, 0], orientation_frames[0][i, 0, 0]))),
                facecolor=colors[i],
                edgecolor="black",
                linewidth=0.5,
            )
        else:
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
            if rod_rel is not None:
                circle.set_xy(rod_outline(frame_index, idx))
                continue
            if height_rel is not None:
                rot = orientation_frames[frame_index][idx]
                circle.set_center((positions[idx, 0], positions[idx, 1]))
                circle.width = diameters[idx]
                circle.height = diameters[idx] * height_rel[idx]
                circle.angle = float(np.degrees(np.arctan2(rot[1, 0], rot[0, 0])))
                continue
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
