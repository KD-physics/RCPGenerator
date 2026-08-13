# rcpgenerator

`rcpgenerator` generates dense packings of polydisperse particles with a
Python interface backed by a compiled C++/OpenMP engine. It supports disks,
spheres, and higher-dimensional hyperspheres in periodic boxes, hard-wall
boxes, and curved containers. The stateful `Packing` class covers the usual
workflow from initialization through relaxation, inspection, persistence, and
visualization. Public dictionary-based functions are also available for callers
that need direct access to the compiled layer.

The distribution also contains `rcptools`, a research toolkit for packing
analysis, searches, job scheduling, persistence, and viewer-bundle export.

Project resources: [repository](https://github.com/KD-physics/RCPGenerator) ·
[Getting Started notebook](https://github.com/KD-physics/RCPGenerator/blob/main/getting_started.ipynb) ·
[examples](https://github.com/KD-physics/RCPGenerator/tree/main/python_code/python/examples) ·
[hosted packing viewer](https://kd-physics.github.io/RCPGenerator/webapp/index.html) ·
[issue tracker](https://github.com/KD-physics/RCPGenerator/issues)

## Installation

```bash
pip install rcpgenerator
```

Version 1.0.0 wheels are built and tested for CPython 3.10–3.14 on:

- Windows x86-64
- Linux x86-64
- macOS 15 or later on Apple Silicon
- macOS 14 or later on Intel

GitHub Actions has verified wheel building, dependency repair, clean
installation, and installed-package smoke tests for every listed Python version
and platform. Wheels have also been installed and run manually on Windows
x86-64, macOS Apple Silicon, and Linux x86-64. The macOS Intel wheels are
verified through GitHub Actions; no manual test on a physical Intel Mac is
currently claimed.

Installing a compatible wheel does not require Git, CMake, Visual Studio,
Xcode, a local C++ compiler, or a separately configured OpenMP installation.
The compiled extension and its required OpenMP runtime are handled by the wheel.

## Quickstart

```python
import numpy as np
import rcpgenerator

# The installed extension defaults to one thread unless OMP_NUM_THREADS is set.
rcpgenerator.set_num_threads(1)

packing = rcpgenerator.Packing(
    phi=0.08,
    N=64,
    Ndim=2,
    box=[1.0, 1.0],
    walls=[0, 0],              # 0 = periodic, 1 = hard wall
    dist={"type": "mono", "d": 1.0},
    neighbor_max=0,            # choose neighbor capacity automatically
    seed=20240719,
)
result = packing.pack()

positions = np.asarray(packing.positions)   # shape (N, Ndim)
diameters = np.asarray(packing.diameters)   # shape (N,)

print(packing.summary())
print("final packing fraction:", packing.phi_final)
print("completed steps:", packing.steps)

# In an interactive session:
# packing.show_packing()
```

Construction initializes a non-overlapping dilute configuration immediately;
`pack()` runs the expansion–relaxation algorithm on that realized state. The
[`examples` directory](https://github.com/KD-physics/RCPGenerator/tree/main/python_code/python/examples)
contains complete box, bidisperse, polydisperse, circular, cylindrical,
spherical, and staged-packing-fraction cases. The
[Getting Started notebook](https://github.com/KD-physics/RCPGenerator/blob/main/getting_started.ipynb)
is the easiest entry point for adapting a worked case.

## Public API at a glance

The `rcpgenerator` package exports exactly these public names:

| Name | Purpose |
| --- | --- |
| `Packing` | Recommended stateful initialization, packing, editing, and visualization API. |
| `initialize_particles(config)` | Low-level dictionary API for initialization only. |
| `run_packing(input, config)` | Low-level dictionary API for a standard packing run. |
| `set_num_threads(n)` | Set the OpenMP thread count; `n` must be at least 1. |
| `get_num_threads()` | Return the current OpenMP thread count. |
| `render_packing(packing, ...)` | Render a 2D or 3D packing. |
| `animate_packing_2d(packing, ...)` | Animate a previously captured 2D trajectory. |

Use `Packing` unless you specifically need the dictionary interface.

## `Packing` configuration

```python
packing = rcpgenerator.Packing(
    phi=0.05,
    N=1000,
    Ndim=3,
    box=[1.0, 1.0, 1.0],
    walls=[0, 0, 0],
    fix_height=False,
    dist={"type": "mono", "d": 1.0},
    neighbor_max=0,
    seed=0,
    verbose=False,
)
```

| Argument | Type and default | Meaning |
| --- | --- | --- |
| `phi` | `float`, `0.05` | Initial particle-volume fraction. It must be between 0 and 1 and low enough for non-overlapping random placement. |
| `N` | `int`, no usable default | Number of particles; must be positive. |
| `Ndim` | `int`, no usable default | Number of spatial dimensions; must be at least 2. |
| `box` | sequence of `float` | Box length in each dimension. Supply exactly `Ndim` positive entries. |
| `walls` | sequence of `int` | Boundary mode in each dimension. Defaults to periodic entries when `box` is supplied. See [Boundaries and containers](#boundaries-and-containers). |
| `fix_height` | `bool`, `False` | Treat the final box dimension as a fixed multiple of the largest initial particle diameter and scale that dimension with particle growth. |
| `dist` | `dict` | Diameter-distribution specification. See [Diameter distributions](#diameter-distributions). |
| `neighbor_max` | nonnegative `int`, `0` | Base used to preallocate each particle's neighbor row; `0` selects the dimension-specific automatic value. See [Neighbor capacity and memory](#neighbor-capacity-and-memory). |
| `seed` | nonnegative `int`, `0` | Deterministic seed used for diameter generation, initial placement, and packing internals. Zero is a valid deterministic seed. |
| `verbose` | `bool`, `False` | Default progress-printing choice used by `pack()` when its `verbose` argument is omitted. |

Unknown constructor keywords raise `TypeError`. Because construction calls
`initialize()` immediately, callers must provide valid `N`, `Ndim`, `box`, and
related configuration rather than relying on the placeholder defaults for `N`
and `Ndim`.

The initializer generates diameters, rescales them together to the requested
initial `phi`, sorts particles by decreasing diameter, and places them without
overlap. Distribution parameters therefore specify relative sizes or the shape
of a distribution; the realized absolute diameters are set by `phi` and `box`.

### Diameter distributions

Distribution names are lowercase and case-sensitive.

| `dist["type"]` | Required or relevant keys | Meaning |
| --- | --- | --- |
| `"mono"` | `d` (default `1.0`) | All particles begin with one relative diameter. |
| `"gaussian"` | `mu`, `sigma` (defaults `1.0`, `0.2`) | Absolute value of normally distributed samples. |
| `"bigaussian"` | `mu1`, `sigma1`, `mu2`, `sigma2`, `p` | Two Gaussian populations; `p` is the fraction drawn from the first population. |
| `"bidisperse"` | `d1`, `d2`, `p` (defaults `1.0`, `2.0`, `0.5`) | Two discrete relative diameters; `p` is the fraction assigned `d1`. |
| `"lognormal"` | `mu`, `sigma` | Lognormal samples using the supplied underlying-normal parameters. |
| `"flat"` | `d_min`, `d_max` (defaults `0.5`, `1.5`) | Uniform samples, with the realized sample extrema mapped to the requested bounds. |
| `"powerlaw"` | `d_min`, `d_max`, `exponent` (default `-2.5`) | Truncated power-law samples, including the logarithmic special case at exponent `-1`. |
| `"exponential"` | `d_min`, `d_max` | Bounded exponential samples, with realized extrema mapped to the requested bounds. |
| `"weibull"` | `shape`, `scale` (defaults `2.0`, `1.0`) | Weibull-distributed relative diameters. |
| `"custom"` | `custom` | Explicit list of exactly `N` relative diameters. |

Examples:

```python
bidisperse = {"type": "bidisperse", "d1": 1.0, "d2": 1.4, "p": 0.5}
powerlaw = {
    "type": "powerlaw",
    "d_min": 1.0,
    "d_max": 10.0,
    "exponent": -2.5,
}
custom = {"type": "custom", "custom": [1.0, 1.0, 1.4, 1.4]}
```

For `custom`, the association between list position and final particle index is
not preserved because particles are sorted by diameter during initialization.

### Boundaries and containers

For ordinary rectangular boxes, each `walls` entry controls one axis:

| Value | Boundary |
| --- | --- |
| `0` | Periodic along that axis. |
| `1` | Hard wall at both ends of that axis. |

Mixed boundaries are supported, for example `walls=[0, 0, 1]` for periodic
`x` and `y` with hard walls in `z`.

A negative first entry selects one curved hard boundary over the first `t`
dimensions, where `t = -walls[0]`. The first box length is the curved
container's diameter:

| Configuration | Container |
| --- | --- |
| `Ndim=2`, `walls=[-2, -1]` | Circular 2D container. |
| `Ndim=3`, `walls=[-2, -1, 1]` | Cylindrical cross-section in `x-y`, with hard end walls in `z`. |
| `Ndim=3`, `walls=[-3, -1, -1]` | Spherical 3D container. |

The core normalizes the subordinate curved entries, so the examples commonly
use forms such as `[-2, 0]`, `[-2, 0, 1]`, or `[-3, 0, 0]`. Explicit forms are
shown above to make the geometry clear. Curved-container `Packing.phi_final`
and `Packing.phi_history` are reported relative to the actual circle, cylinder,
or sphere volume rather than its rectangular bounding box.

With `fix_height=True`, the last initial `box` entry is interpreted as a
multiple of the largest initial particle diameter. This mode is useful for
confined geometries whose height should grow consistently with the particle
scale.

## Running and controlling a packing

### `Packing.initialize()`

```python
initial = packing.initialize()
```

Regenerates diameters and non-overlapping positions from the current
configuration. It returns a dictionary with:

- `positions`: nested list with shape `(N, Ndim)`;
- `diameters`: list of length `N`;
- `box` and normalized `walls`;
- `diameter_scale_factor`: common scale applied to generated diameters;
- `phi_modifier`: curved-container volume fraction, or `1.0` for an ordinary box.

Assigning `phi`, `N`, `Ndim`, `box`, `walls`, `fix_height`, or `dist` marks the
object for reinitialization. The next `pack()` calls `initialize()`
automatically. Assigning `seed`, `neighbor_max`, `positions`, or `diameters`
marks only the packing result stale.

### `Packing.pack()`

```python
result = packing.pack(
    verbose=False,
    progress_interval=1000,
    capture_trajectory=False,
    trajectory_interval=1000,
)
```

| Argument | Meaning |
| --- | --- |
| `verbose` | Print progress. `None` uses the constructor's `verbose` value. |
| `progress_interval` | Step interval between printed updates when verbose output is enabled. |
| `capture_trajectory` | Record sampled positions, diameters, and scalar diagnostics for animation or analysis. |
| `trajectory_interval` | Step interval between captured trajectory samples. |

`pack()` returns the raw result dictionary and also updates the corresponding
`Packing` attributes. Trajectory capture consumes additional memory roughly in
proportion to the number of samples times `N * Ndim`.

### `Packing.relax()`

`relax()` runs the current realized state with explicit step and growth
controls:

```python
result = packing.relax(
    n_steps=5000,
    mu=None,
    fix_diameter=False,
    target_phi=None,
    verbose=False,
    progress_interval=1000,
    capture_trajectory=False,
    trajectory_interval=1000,
)
```

| Argument | Meaning |
| --- | --- |
| `n_steps` | Positive maximum number of steps for this run. It is a cap, not a promise that all steps will execute; convergence may terminate earlier. |
| `mu` | Optional override for the packing algorithm's chemical-potential parameter. `None` uses the internal schedule. |
| `fix_diameter` | Hold all diameters fixed and relax positions only. |
| `target_phi` | Grow until this packing fraction is reached, then hold diameters fixed while relaxing. |
| Remaining arguments | Same progress and trajectory controls as `pack()`. |

`target_phi` and `fix_diameter=True` cannot be combined. A common staged
workflow is:

```python
packing.relax(n_steps=5000, target_phi=0.82)
packing.update_phi(0.84)
packing.relax(n_steps=1000, fix_diameter=True)
```

## Results and diagnostics

After `pack()` or `relax()`, the most commonly used attributes are:

| Attribute | Meaning |
| --- | --- |
| `positions` | Final particle centers, shape `(N, Ndim)`. Particle order can change during spatial reordering. |
| `diameters` | Final diameters, length `N`. |
| `box`, `walls` | Realized container configuration. |
| `steps` | Number of completed steps. |
| `phi_final` | Final packing fraction; corrected to actual container volume for curved boundaries. |
| `max_min_dist` | Core convergence/separation diagnostic. |
| `force_magnitude` | Final reported mean-force diagnostic. |
| `phi_history` | Packing fraction after every completed step. |
| `force_history` | Mean-force history. Contact-free phases may contain `NaN` because mean contact force is undefined there. |
| `energy_history` | Energy history. |
| `max_overlap_history` | Maximum pair-overlap history. |
| `mu_history` | Chemical-potential value by step. |
| `alpha_history` | Learning-rate value by step. |
| `mu_flag_history` | Chemical-potential rung indicator (`-1`, `0`, or `1`) by step. |

The step histories have length `steps`; they are not padded to a requested
maximum. Check outputs using the validation appropriate to your scientific
application. For example:

```python
import math
import numpy as np

assert np.asarray(packing.positions).shape == (packing.N, packing.Ndim)
assert np.asarray(packing.diameters).shape == (packing.N,)
assert np.isfinite(packing.positions).all()
assert np.isfinite(packing.diameters).all()
assert len(packing.phi_history) == packing.steps
assert len(packing.force_history) == packing.steps
assert len(packing.energy_history) == packing.steps
assert math.isfinite(packing.phi_final)
assert not np.isinf(packing.force_history).any()
```

The raw result returned by `pack()` or `relax()` contains all fields above plus
advanced engine diagnostics:

- `particle_origin`, mapping current particle order to the order at run start;
- `pp_F_sq_sum`, `pp_F_sign_flips`, `pp_pos_range_max`, `pp_contact_sum`,
  `pp_signflip_resets`, and `pp_window_steps`, which summarize the final
  diagnostic window per particle;
- `spike_log`, a flat array with eight values per event:
  `[step, force_ratio, force, max_overlap, mu_flag, alpha,
  step_mod_reorder, step_mod_500]`;
- `reorder_log`, a flat array with five values per event:
  `[step, mu_flag, pairs_before, pairs_after, force]`.

These advanced fields are available in the returned dictionary but are not all
copied to named `Packing` attributes or included by `to_dict()`. Also note that
raw `result["phi"]` is the compiled core's bounding-box fraction. For curved
containers, use `packing.phi_final` and `packing.phi_history`, which apply the
container-volume correction.

### Captured trajectory

When `capture_trajectory=True`, these attributes contain sampled values:

- `trajectory_steps`;
- `trajectory_positions` and `trajectory_diameters`;
- `trajectory_phi`, `trajectory_force`, and `trajectory_energy`;
- `trajectory_max_min_dist`.

They are cleared at the start of every new packing execution.

## Editing, copying, and serialization

| Method | Behavior |
| --- | --- |
| `update_phi(target_phi)` | Isotropically rescale diameters to the requested current packing fraction. In fixed-height mode, also scales the last box coordinate and particle coordinates. Returns the realized fraction. |
| `update_box(box)` | Rescale positions by axis and diameters isotropically so the current packing fraction is preserved. Returns the new box list. |
| `update_height(height)` | Convenience form of `update_box()` for the final dimension. |
| `copy()` / `clone()` | Deep-copy the serializable state into an independent `Packing`. |
| `to_dict()` | Return initialization config, packing config, realized state, principal final results, trajectory, and status flags. |
| `Packing.from_dict(data)` | Reconstruct a `Packing` without running initialization. |
| `reset()` | Clear realized state, results, and trajectory; the next run reinitializes. |
| `summary()` | Return a concise human-readable state summary. |

Example JSON persistence:

```python
import json

with open("packing.json", "w", encoding="utf-8") as stream:
    json.dump(packing.to_dict(), stream)

with open("packing.json", encoding="utf-8") as stream:
    restored = rcpgenerator.Packing.from_dict(json.load(stream))
```

`to_dict()` preserves the principal histories and captured trajectory, but not
every advanced raw diagnostic listed above. Save the raw result separately if
those fields are required.

State flags are available as read-only properties:

- `initialized`: a realized initial state exists;
- `packed`: the current realized state has a completed packing result;
- `needs_initialize`: configuration changes require regeneration;
- `needs_pack`: the realized state has no current packing result.

## Neighbor capacity and memory

The pair list uses a preallocated row for each particle. `neighbor_max` is the
base term in that capacity estimate; it is not a strict global maximum number
of physical neighbors. With particle diameter `D_i`, minimum diameter `D_min`,
and dimension `d`, the implementation allocates approximately

```text
ceil((d / 3) * (50 * (D_i / D_min) ** (1 + d / 6) + base))
```

slots for particle `i`, capped at `N`. With `neighbor_max=0`, the automatic
base is currently dimension-specific (300 in 2D, 750 in 3D, and 5500 in 4D and
higher). A larger particle-size ratio can therefore increase memory use
substantially even when the automatic setting succeeds.

If capacity is exhausted, the core raises a `RuntimeError` similar to:

```text
rcpgenerator: particle 17 exceeded its neighbor capacity (750 slots)
during pair-list construction (...).
```

Increase the base and rerun, for example:

```python
packing.neighbor_max = 1500
result = packing.pack()
```

The assignment marks the current result stale but does not regenerate the
particle positions. If memory use is the problem, reducing an extreme
`D_max / D_min` ratio can be more effective than increasing `neighbor_max`.

## Threads and reproducibility

The extension defaults to one OpenMP thread unless the `OMP_NUM_THREADS`
environment variable is already set when it is imported.

```python
import rcpgenerator

rcpgenerator.set_num_threads(4)
print(rcpgenerator.get_num_threads())
```

`set_num_threads(n)` requires `n >= 1`. Set the thread count before a run. A
fixed seed and one thread provide the strongest run-to-run reproducibility.
Multi-threaded reductions need not be bit-identical because floating-point
summation order can vary. For reproducible research, record the package
version, seed, thread count, inputs, compiler/wheel provenance, and platform.

## Rendering and animation

`Packing.render()` and the equivalent module function `render_packing()` support
2D and 3D packings:

```python
packing.render(
    path="packing.png",       # omit to avoid saving
    show=False,
    palette_choice=1,         # palettes 1 through 12; other values fall back to 1
    figsize=(6, 6),
    dpi=150,
)

packing.savefig("packing.svg", palette_choice=3)
packing.show_packing(palette_choice=3)
```

`render()` returns the output path as a string when saving, otherwise `None`.
It raises `ValueError` outside 2D and 3D. `savefig()` is a convenience wrapper
with the renderer's default figure size and DPI.

To animate 2D relaxation, capture a trajectory during the run:

```python
packing.pack(capture_trajectory=True, trajectory_interval=250)
packing.animate_2d(
    path="packing.gif",       # GIF is the only supported saved format
    show=False,
    palette_choice=1,
    interval_ms=80,
    repeat=False,
)
```

`animate_2d()` returns a Matplotlib `FuncAnimation`. It requires a 2D packing
and a captured trajectory. In a notebook, leaving `show=True` displays an HTML
animation when possible.

## Low-level dictionary API

The public low-level functions are useful for stateless integrations. The
high-level class performs state tracking and curved-container reporting for you.

### `initialize_particles(config)`

The configuration keys are `phi`, `N`, `Ndim`, `box`, `walls`, `fix_height`,
`dist`, and `seed`, with the same meanings as the `Packing` constructor.

```python
initial = rcpgenerator.initialize_particles({
    "phi": 0.08,
    "N": 64,
    "Ndim": 2,
    "box": [1.0, 1.0],
    "walls": [0, 0],
    "fix_height": False,
    "dist": {"type": "bidisperse", "d1": 1.0, "d2": 1.4, "p": 0.5},
    "seed": 20240719,
})
```

The return fields are documented under [`Packing.initialize()`](#packinginitialize).

### `run_packing(input, config)`

```python
result = rcpgenerator.run_packing(
    {
        "positions": initial["positions"],
        "diameters": initial["diameters"],
    },
    {
        "box": initial["box"],
        "walls": initial["walls"],
        "neighbor_max": 0,
        "seed": 20240719,
        "fix_height": False,
    },
)
```

The `input` dictionary requires `positions` and `diameters`. It also accepts an
advanced `fixed` mask of length `N`, where `1` exempts that particle's diameter
from global diameter scaling and `0` leaves it active. The `config` dictionary
accepts `box`, `walls`, `neighbor_max`, `seed`, and `fix_height`.

This function performs a standard run with the engine defaults. Explicit
maximum-step, target-fraction, progress-observer, and trajectory controls are
provided by `Packing.relax()` and `Packing.pack()`, not by public
`run_packing()` arguments.

## Common errors and remedies

| Symptom | Likely cause and remedy |
| --- | --- |
| `Missing or invalid phi/N` | Use `0 < phi < 1` and a positive `N`. |
| `--Ndim must be >=2` | Set `Ndim` to 2 or higher. |
| `--box length mismatch` | Supply exactly `Ndim` box lengths. |
| `--walls expects Ndim entries` | Supply exactly `Ndim` wall entries to the packing run. |
| `custom distribution needs N entries` | Provide exactly one custom diameter per particle. |
| Initial placement exceeds its attempt limit | The requested initial `phi`, size ratio, or confined geometry is too difficult for non-overlapping random placement. Lower initial `phi` or revise the geometry/distribution. |
| `exceeded its neighbor capacity` | Increase `neighbor_max` or reduce the particle-size ratio; see [Neighbor capacity and memory](#neighbor-capacity-and-memory). |
| No trajectory recorded | Call `pack(capture_trajectory=True, ...)` before animation. |
| Unsupported rendering dimension | Built-in rendering is limited to 2D and 3D, although the packing core supports higher dimensions. |

Use the exact lowercase distribution names documented above. An unsupported
distribution name is not a recoverable Python validation path in the current
initializer and should be avoided.

## Interactive packing viewer

`rcptools.bridge.write_bundle` creates a viewer bundle containing
`manifest.json`, `pos.f32`, and `dia.f32`:

```python
from rcptools.bridge import write_bundle

write_bundle(
    "my_packing",
    packing.positions,
    packing.diameters,
    packing.box,
    walls=packing.walls,
)
```

Open the [hosted packing viewer](https://kd-physics.github.io/RCPGenerator/webapp/index.html)
and select the generated bundle folder. Solid boundaries represent hard walls
and dashed boundaries represent periodic boundaries; the viewer can inspect 2D
packings directly and 3D packings by cross-section.

The viewer application itself is not included in the Python wheel. Its
development files and `launch.py` helper are repository-only and live in the
[`python_code/python/webapp` directory](https://github.com/KD-physics/RCPGenerator/tree/main/python_code/python/webapp).

## `rcptools`

The installed distribution also provides the `rcptools` namespace. It is a
research-workflow toolkit rather than part of the small core generator API. Its
modules include:

| Area | Modules and purpose |
| --- | --- |
| Data model and persistence | `model`, `persistence`, `snapshot`, and `bundling` for representing, saving, restoring, and grouping runs. |
| Execution and searches | `jobs`, `scheduler`, `queue_scheduler`, `search`, `search_v2`, `sweep`, and `time_aware`. |
| Diagnostics and quality | `diagnostics`, `overlap`, and `census`. |
| Viewer export | `bridge`, including `write_bundle`. |
| Scientific analyses | `analysis.chord_*`, `analysis.porefill`, `analysis.porenet`, and `analysis.voidfill_method`. |
| Study helpers | `studies` and optional Voronoi-related functionality in `voronoi`. |

Some `rcptools` workflows require method-specific inputs or optional external
libraries. Treat their module docstrings and the
[repository examples](https://github.com/KD-physics/RCPGenerator/tree/main/python_code/python/examples)
as the authoritative interface for those specialized research workflows.

## Interpretation and limitations

Packing is an iterative numerical procedure. Convergence behavior, memory use,
and runtime depend on particle count and size distribution, dimensionality,
boundaries, initial state, thread count, and packing controls. Reaching a step
cap means the run stopped at that cap; it does not by itself establish the same
condition as an earlier convergence termination.

Inspect the returned histories and diagnostics and apply validation appropriate
to the intended scientific use. The package does not label every returned state
as perfectly force-balanced or mechanically stable. Built-in visual rendering
is limited to 2D and 3D, while the numerical core supports higher dimensions.

## Building from source

Source builds require Git, a C++17 compiler, CMake-compatible build tooling, and
a usable OpenMP development/runtime installation. Platform-specific compiler
setup is the builder's responsibility.

```bash
git clone https://github.com/KD-physics/RCPGenerator.git
cd RCPGenerator/python_code/python
pip install -v .
```

Normal source builds use portable CPU instructions. A local build may explicitly
tune the compiled core for the current machine:

```bash
CMAKE_ARGS="-DRCP_NATIVE_OPTIMIZATION=ON" pip install -v .
```

In PowerShell:

```powershell
$env:CMAKE_ARGS = "-DRCP_NATIVE_OPTIMIZATION=ON"
pip install -v .
```

GCC native builds enable `-march=native` and the configured 512-bit vector-width
preference; Clang native builds enable `-march=native`. Other compilers emit a
CMake warning and retain their default architecture. A native-optimized build
is intended only for the machine on which it is compiled. Do not redistribute
it as a portable wheel.

## Method and citation

The package uses an iterative expansion–relaxation method based on the approach
described by Desmond and Weeks, with a compiled neighbor-search and relaxation
implementation.

- Desmond, K. “Random close packing at extreme size ratios with an Adam-based
  inflation protocol.” arXiv:2608.12235 (2026).
  [arXiv:2608.12235](https://arxiv.org/abs/2608.12235)
- Desmond, K. W. and Weeks, E. R. “Random close packing of disks and spheres in
  confined geometries.” *Physical Review E* **80**, 051305 (2009).
  [arXiv:0903.0864](https://arxiv.org/abs/0903.0864)
- Desmond, K. W. and Weeks, E. R. “Influence of particle size distribution on
  random close packing of spheres.” *Physical Review E* **90**, 022204 (2014).
  [arXiv:1303.4627](https://arxiv.org/abs/1303.4627)

Use the RCPGenerator arXiv manuscript as the preferred scientific citation and
record the `rcpgenerator` version used. The two earlier papers provide method
background.

`rcpgenerator` is distributed under the
[MIT License](https://github.com/KD-physics/RCPGenerator/blob/main/python_code/python/LICENSE).
