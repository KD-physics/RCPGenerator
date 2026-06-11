# rcpgenerator

Fast random close packing (RCP) generator for polydisperse spheres in N
dimensions. A parallel C++/OpenMP core with SIMD-accelerated inner loops,
wrapped in a stateful Python API. The user supplies particle positions and
diameters (or a built-in distribution) and receives a jammed configuration.
Includes rendering helpers and optional trajectory capture.

This release (`python_parallel3`) is a performance release: the Python API
is unchanged from `python_parallel2`, while the core is substantially
faster — large high-polydispersity packings (hundreds of thousands of
particles, size ratios of 30–60×) complete several times faster than the
previous release through improved neighbor management, vectorized overlap
screening, fused integrator passes, and a shortened (statistically
validated) compression schedule.

## Install

```bash
git clone https://github.com/KD-physics/RCPGenerator.git
cd RCPGenerator/python_parallel3
pip install -v .
```

Build requirements pulled automatically from `pyproject.toml`:
`scikit-build-core`, `pybind11`. System: C++17 compiler with OpenMP and
CMake 3.18+. The build targets the host CPU (`-march=native`); AVX-512
fast paths are used when available, with scalar fallbacks otherwise.

```python
import rcpgenerator
```

## Quick start

```python
import rcpgenerator

p = rcpgenerator.Packing(
    phi=0.11,
    N=250,
    Ndim=2,
    box=[1.0, 1.0],
    walls=[0, 0],
    fix_height=False,
    dist={"type": "mono", "d": 1.0},
    neighbor_max=0,
    seed=123,
)

p.pack()

print(p.phi_final, p.steps, p.force_magnitude, p.max_min_dist)
```

## `rcpgenerator.Packing`

Stateful pure-Python wrapper. Construction initializes particle positions
immediately.

### Methods

| Method | Signature |
| --- | --- |
| Constructor | `Packing(phi, N, Ndim, box, walls, fix_height, dist, neighbor_max=0, seed=None, verbose=False)` |
| `initialize()` | Re-run initialization with current config |
| `pack(...)` | `pack(verbose=None, progress_interval=1000, capture_trajectory=False, trajectory_interval=1000)` |
| `relax(...)` | `relax(n_steps, mu=None, fix_diameter=False, target_phi=None, verbose=False, progress_interval=1000, capture_trajectory=False, trajectory_interval=1000)` |
| `update_phi(target_phi)` | Rescale diameters to hit a given phi |
| `update_box(box)` | Update box; rescale positions and diameters consistently |
| `update_height(height)` | Update last box dimension |
| `render(...)` | `render(path=None, show=False, palette_choice=1, figsize=(4,4), dpi=110)` |
| `savefig(path, ...)` | Render to file |
| `show_packing(palette_choice=1, ...)` | Inline display |
| `animate_2d(...)` | `animate_2d(path=None, show=True, palette_choice=1, interval_ms=80, repeat=False)` — requires recorded 2D trajectory |
| `summary()` | Short status string |
| `reset()` | Keep config, clear realized state |
| `copy()`, `clone()` | Independent copies |
| `to_dict()`, `Packing.from_dict(d)` | Serialize / reconstruct |

### Attributes

- **Init config:** `phi`, `N`, `Ndim`, `box`, `walls`, `fix_height`, `dist`
- **Pack config:** `neighbor_max`, `seed`
- **Realized state:** `positions`, `diameters`
- **Results:** `steps`, `phi_final`, `max_min_dist`, `force_magnitude`,
  `phi_history`, `force_history`, `energy_history`
- **Flags:** `initialized`, `packed`, `needs_initialize`, `needs_pack`
- **Trajectory (if captured):** `trajectory_positions`, `trajectory_diameters`,
  `trajectory_steps`, `trajectory_phi`, `trajectory_force`,
  `trajectory_energy`, `trajectory_max_min_dist`

### Stateful behavior

- Construction runs initialization.
- Mutating an init-config attribute marks `needs_initialize = True`.
- Mutating a pack-config or realized-state attribute marks `needs_pack = True`.
- `pack()` re-initializes only when needed.
- `reset()` keeps config, clears realized state and results.

## Constructor parameters

| Param | Type | Meaning |
| --- | --- | --- |
| `phi` | float | Initial volume fraction used for placement |
| `N` | int | Number of particles |
| `Ndim` | int | Dimension |
| `box` | list[float] | Box / container size parameters (length `Ndim`) |
| `walls` | list[int] | Boundary convention per axis (see Geometries) |
| `fix_height` | bool | Fixed-height semantics for the last dimension |
| `dist` | dict | Diameter distribution (see Distributions) |
| `neighbor_max` | int | Per-particle neighbor allocation hint; `0` uses dimension- and size-ratio-scaled defaults |
| `seed` | int \| None | Seed carried through the API |
| `verbose` | bool | Default verbosity for `pack()` / `relax()` |

### Geometries (`walls`)

| Geometry | `Ndim` | `walls` |
| --- | --- | --- |
| Periodic box | 2 | `[0, 0]` |
| Periodic box | 3 | `[0, 0, 0]` |
| Circle | 2 | `[-2, 0]` |
| Cylinder | 3 | `[-2, 0, 0]` |
| Sphere | 3 | `[-3, 0, 0]` |

### Distributions (`dist`)

| `dist["type"]` | Required fields |
| --- | --- |
| `mono` | `d` |
| `bidisperse` | `d1`, `d2`, `p` |
| `lognormal` | `mu`, `sigma` |
| `flat` | `d_min`, `d_max` |
| `powerlaw` | `d_min`, `d_max`, `exponent` |
| `custom` | `custom` (list of `N` diameters) |

## `pack()` parameters

| Param | Default | Meaning |
| --- | --- | --- |
| `verbose` | `None` | Overrides instance `verbose` |
| `progress_interval` | `1000` | Printed-progress step interval |
| `capture_trajectory` | `False` | Record sampled snapshots |
| `trajectory_interval` | `1000` | Snapshot step interval |

## `relax()` parameters

| Param | Default | Meaning |
| --- | --- | --- |
| `n_steps` | required | Hard upper bound; convergence may exit earlier |
| `mu` | `None` | Override the diameter-growth coupling for this call |
| `fix_diameter` | `False` | Freeze diameter evolution (true fixed-diameter mode, not `mu=0`) |
| `target_phi` | `None` | Evolve toward target phi; on reach, lock diameters and continue |
| `verbose` | `False` | Per-call verbosity |
| `progress_interval` | `1000` | Printed-progress step interval |
| `capture_trajectory` | `False` | Record sampled snapshots |
| `trajectory_interval` | `1000` | Snapshot step interval |

`target_phi` and `fix_diameter=True` are mutually exclusive.

## Staged target-phi workflow

```python
p.relax(n_steps=5000, target_phi=0.82)
for target_phi in [0.84, 0.86, 0.88, 0.90, 0.92, 0.94]:
    p.update_phi(target_phi)
    p.relax(n_steps=1000, fix_diameter=True)
```

`update_phi` rescales diameters deterministically; the subsequent
`fix_diameter=True` relax settles positions at that phi.

## Runtime controls

### Threads

```python
import rcpgenerator
rcpgenerator.set_num_threads(4)
rcpgenerator.get_num_threads()
```

Set before `pack()` / `relax()`.

### `RCP_SPEED` (environment variable)

Controls the inner-loop budget for the mu-oscillation convergence
schedule.

| Value | Budget |
| --- | --- |
| `immediate` | Smallest; fastest exit |
| `quick` *(default)* | Balanced |
| `patient` | Larger; tighter convergence |
| `forever` | Very large |

```bash
RCP_SPEED=patient python my_script.py
```

The preset sets the rung-1 deadline as
`deadline = Y * 4000 + X * sqrt(N * Ndim) - 2000` and the number of
oscillation cycles. `N` enters automatically. The early compression
phase runs a front-loaded diameter-growth schedule (5x, relaxing over a
few hundred steps), which is what makes the shortened deadline safe:
final packing fractions are statistically indistinguishable from the
long-schedule result while total step counts drop ~13%.

Exact overrides:

- `RCP_RUNG1_DEADLINE` — integer step count; overrides the preset's
  deadline directly (no implicit `-2000`)
- `RCP_OSC_CYCLES` — integer cycle count; overrides the preset's cycle count
- `RCP_MU1_BOOST` — initial compression boost factor (default `5`; `1`
  disables and restores the flat schedule)
- `RCP_MU1_BOOST_TAU` — boost relaxation timescale in steps (default `350`)

### Diagnostics (optional)

```bash
RCP_PROFILE=1 python my_script.py
```

Prints a per-kernel timing breakdown and periodic progress snapshots to
stderr. Off by default.

### Self-check (optional)

```bash
RCP_SENTINEL=1000 python my_script.py
```

Every `RCP_SENTINEL` steps the engine rebuilds the neighbor structure
from scratch with an independent method and verifies that no
force-generating pair is missing from the live lists, printing a
diagnostic on any discrepancy. Useful when validating new parameter
regimes; adds overhead proportional to the check frequency.

## Trajectory and animation

```python
p.pack(capture_trajectory=True, trajectory_interval=500)
p.animate_2d(path="packing.gif", show=False, palette_choice=3)
```

`animate_2d` requires a previously recorded 2D trajectory.

## Thin binding API

Dict-based functions exposed for callers that don't want the wrapper:

```python
import rcpgenerator

state = rcpgenerator.initialize_particles({
    "phi": 0.11, "N": 4, "Ndim": 2,
    "box": [1.0, 1.0], "walls": [0, 0], "fix_height": False,
    "dist": {"type": "mono", "d": 1.0},
})

result = rcpgenerator.run_packing(
    {"positions": state["positions"], "diameters": state["diameters"]},
    {"box": state["box"], "walls": state["walls"],
     "neighbor_max": 0, "seed": 123, "fix_height": False},
)
```

Exports from `rcpgenerator`:

- `Packing`
- `initialize_particles(config)`
- `run_packing(input, config)`
- `set_num_threads(n)`, `get_num_threads()`
- `render_packing(...)`, `animate_packing_2d(...)`

## Performance notes

- The core parallelizes with OpenMP; thread count is set via
  `rcpgenerator.set_num_threads(n)`. Throughput is strong at modest
  thread counts (2–4); the force loop and neighbor management both
  scale.
- Neighbor lists use a kd-tree backend in 2D/3D and a dimension-generic
  sort-and-sweep backend for `Ndim > 3`, selected automatically.
- Overlap screening in 2D/3D uses 8-wide AVX-512 batching when the host
  supports it; results are bit-identical to the scalar path.
- Indicative scale (4 threads, default `quick` preset): 500,000
  particles at ~38x size ratio in 3D packs in roughly 20–25 minutes on a
  modern desktop CPU; 100,000 particles in ~3 minutes.

## `rcp_search/`

The `rcp_search` package colocated in this folder is a separate research
study (search protocols over diameter distributions) that uses
`rcpgenerator` as its engine. It is not part of the `rcpgenerator`
package and will move to its own repository in the future; it lives here
for import convenience (`sys.path` addition or running notebooks from
this folder).

## `tests/`

C++ parity tests for the initializer and packer, plus a hard-wall
regression test. Build with `-DBUILD_PARITY_TESTS=ON` and run via
`ctest`.

## License

MIT — see `LICENSE`.
