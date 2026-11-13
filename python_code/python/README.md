# rcpgenerator

`rcpgenerator` generates dense packings of polydisperse particles with a
Python interface backed by a compiled C++/OpenMP engine. It supports disks,
spheres, and higher-dimensional hyperspheres in periodic boxes, hard-wall
boxes, and curved containers, and it packs ellipsoids and rods in the same boxes
and containers (see [Ellipsoids and rods](#ellipsoids-and-rods)). The stateful `Packing` class covers the usual
workflow from initialization through relaxation, inspection, persistence, and
visualization. Public dictionary-based functions are also available for callers
that need direct access to the compiled layer.

The distribution also contains `rcptools`, a research toolkit for packing
analysis, searches, job scheduling, persistence, and viewer-bundle export.

Project resources: [repository](https://github.com/KD-physics/RCPGenerator) ·
[Getting Started notebook](https://github.com/KD-physics/RCPGenerator/blob/main/getting_started.ipynb) ·
[examples](https://github.com/KD-physics/RCPGenerator/tree/main/python_code/python/examples) ·
[interactive packing app](https://kd-physics.github.io/RCPGenerator/webapp/index.html) ·
[large-packing visualizer](https://kd-physics.github.io/RCPGenerator/python_code/python/webapp/app/index.html) ·
[issue tracker](https://github.com/KD-physics/RCPGenerator/issues)

## Installation

```bash
pip install rcpgenerator
```

The version 1.2.0 wheel workflow builds and tests CPython 3.10–3.14 on:

- Windows x86-64
- Linux x86-64
- macOS 15 or later on Apple Silicon
- macOS 14 or later on Intel

GitHub Actions verifies wheel building, dependency repair, clean installation,
and installed-package smoke tests for every listed Python version and platform.
Previous wheel builds have also been installed and run manually on Windows
x86-64, macOS Apple Silicon, and Linux x86-64. macOS Intel is verified through
GitHub Actions; no manual test on a physical Intel Mac is claimed.

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

## Built-in examples

`rcpgenerator.example()` lists short, complete snippets by topic: 2D, 3D and higher-dimensional packings,
ellipsoids (one shape or a shape per particle, orientations, permutation, quench), rods (round and flat ends,
mixed lengths, quench, a movie, the viewer), every boundary and container type, arbitrary size distributions, diameter permutation, the quench,
contact forces, plotting, saving, the web viewer, and Voronoi analysis. `example("name")` prints
one, ready to copy, paste and modify:

```python
import rcpgenerator

rcpgenerator.example()                    # the catalogue
rcpgenerator.example("capped_cylinder")   # one snippet
```

```text
Boundary Conditions Examples
  - 'walls_one_axis': hard walls along z
  - 'walls_all_axes': closed hard-wall box
  - 'channel_2d': 2D channel, walls in y
  - 'disk': 2D circular container
  - 'cylinder': periodic-axis cylinder
  - 'capped_cylinder': cylinder with end walls
  - 'tall_box': height in particle diameters
  - 'sphere': spherical 3D container
...
```

Every snippet runs as printed at a small particle count; the test suite executes each one.

## Public API at a glance

The `rcpgenerator` package exports exactly these public names:

| Name | Purpose |
| --- | --- |
| `Packing` | Recommended stateful initialization, packing, editing, and visualization API. |
| `initialize_particles(config)` | Low-level dictionary API for initialization only. |
| `run_packing(input, config)` | Low-level dictionary API for a standard packing run. |
| `set_num_threads(n)` | Set the OpenMP thread count; `n` must be at least 1. |
| `get_num_threads()` | Return the current OpenMP thread count. |
| `render_packing(packing, ...)` | Render a 2D or 3D packing of spheres, ellipsoids, or rods. |
| `animate_packing_2d(packing, ...)` | Animate a previously captured 2D trajectory (discs, ellipses, or rods). |
| `example(name=None)` | List the built-in snippets, or print one to copy. See [Built-in examples](#built-in-examples). |

Use `Packing` unless you specifically need the dictionary interface. Ellipsoids and rods add no
new names: they are selected by `Packing` keywords (`aspect_ratios`, or `rod_length` and `rod_radius`),
and their results are `Packing` attributes (`orientations`, `phi_contact`, `shrink_factor`) and the method
`shrink_to_contact()`. See [Ellipsoids and rods](#ellipsoids-and-rods).

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
| `diameter_permute` | `bool`, `False` | Turn on diameter permutation (which particle carries which diameter is re-assigned during compression). `False` reproduces the original algorithm bit for bit. See [Diameter permutation](#diameter-permutation). |
| `permute_budget` | `float` in (0, 0.5], `0.25` | Exchanges per permutation pass as a fraction of N. Smaller = less uplift. |
| `permute_interval` | `int` ≥ 1, `5` | Steps between permutation passes. Longer = fewer passes and less uplift. |
| `permute_start` | `int`, `4000` | Step at which permutation begins. Later = less uplift. |
| `aspect_ratios` | list of `Ndim - 1` `float`, or `N` such lists; default none | Pack ellipsoids. `D` is the longest axis and axis `k` has length `D × ratio_k`; the ratios are written in decreasing order, each at most 1. One list for all particles, or one per particle. See [Ellipsoids](#ellipsoids). |
| `orientations` | `N` rotation matrices (`Ndim × Ndim`); default generated | Starting orientations for ellipsoids or rods, one per particle; columns are the body axes. Omit to have random ones generated. Not accepted without a shape. |
| `rod_length` | `float` or list of `N`; default none | Pack rods. Length of the straight section in units of `D` (capped rods), or the full length (flat-ended rods). Requires `rod_radius`. See [Rods](#rods). |
| `rod_radius` | `float` or list of `N`; default none | Rod radius in units of `D`; `0.5` makes `D` the rod's thickness. Requires `rod_length`. |
| `rod_capped` | `bool`, `True` | `True`: a half-sphere on each end. `False`: flat ends (2D and 3D only). One value for the packing. |
| `ell_permute_orient` | `int`, `0` | With `diameter_permute=True` and ellipsoids or rods: `0` exchanges diameter and shape and leaves orientations in place; `1` makes the orientation travel with the particle. See [Permutation with ellipsoids and rods](#permutation-with-ellipsoids-and-rods). |
| `ell_iters`, `ell_track_tol` | `int`, `5`; `float`, `1e-6` | Advanced controls of the ellipsoid contact solve during a run: at most `ell_iters + 1` evaluations per pair per step, and the tolerance at which a pair's solution carried over from the previous step is accepted after one evaluation (also used for 3D flat-ended rods). `phi_contact` is always computed from a converged solve. |

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

Only periodic axes wrap. Particles never interact through a hard wall or a curved
container boundary. In 1.0.0 pair distances were also wrapped across walls, so two
particles touching opposite walls at nearly the same position could push on each other
through them: a few such contacts per walled packing, and after a quench those particles
were out of balance by tens of percent of their contact load when distances were measured
in a straight line. Walled packings made with 1.0.0 can therefore differ slightly from the
same seed now. Periodic packings are unchanged.

A negative first entry selects one curved hard boundary over the first `t`
dimensions, where `t = -walls[0]`. The first box length is the curved
container's diameter:

| Configuration | Container |
| --- | --- |
| `Ndim=2`, `walls=[-2, -1]` | Circular 2D container. |
| `Ndim=3`, `walls=[-2, -1, 1]` | Cylindrical cross-section in `x-y`, with hard end walls in `z`. |
| `Ndim=3`, `walls=[-3, -1, -1]` | Spherical 3D container. |

The core normalizes the subordinate curved entries, so the examples commonly
use forms such as `[-2, 0]`, `[-2, 0, 1]`, or `[-3, 0, 0]`. A `walls` list with the
wrong length, a value other than `0` or `1`, or a negative entry anywhere but first
raises `ValueError` naming the valid forms (a curved container must be the first
entry; `rcpgenerator.example()` shows each container). Explicit forms are
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

Assigning `phi`, `N`, `Ndim`, `box`, `walls`, `fix_height`, `dist`, `aspect_ratios`,
`rod_length`, `rod_radius`, or `rod_capped` marks the
object for reinitialization. The next `pack()` calls `initialize()`
automatically. Assigning `seed`, `neighbor_max`, `positions`, `diameters`, or `orientations`
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

A run stops when both the largest pair overlap and the net force on the particles are small (or
at the step cap). Results can differ slightly from versions 1.0.0 and 1.1.0: in rare cases
(a box too small for its largest particles, very elongated shapes) a run takes more steps.

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

To bring a finished packing to mechanical force balance, see
[Quench to force balance](#quench-to-force-balance) (`packing.quench()`).

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
| `orientations` | Ellipsoids and rods: final rotation matrices, shape `(N, Ndim, Ndim)`, in the particle order of `positions`; columns are the body axes, column 0 the longest axis or the rod's axis. Empty for spheres. |
| `aspect_ratios`, `rod_length`, `rod_radius` | The shape. When given per particle, the table or list in the particle order of `positions`. |
| `phi_contact` | Ellipsoids and rods: packing fraction after the exact shrink to zero overlap; the value to report. `None` for spheres. See [Reading the result](#reading-the-result-phi_contact-and-phi_final). |
| `shrink_factor` | Ellipsoids and rods: the common factor of that shrink; multiply every diameter by it. `None` for spheres. |
| `shrink_to_contact()` | Method: the same shrink computed on the current state, for spheres too (periodic boxes, flat walls and curved containers). Returns `scale`, `phi`, `n_overlapping`, `pair`; does not modify the packing. |

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
- `trajectory_max_min_dist`;
- `trajectory_origin`: for each frame, the index at the start of the run of the particle in each slot (the
  packer renumbers particles while it runs; see [Rendering and animation](#rendering-and-animation));
- `trajectory_orientations`: ellipsoids and rods only, one `Ndim × Ndim` rotation matrix per particle per frame.

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

Sphere states saved by versions 1.0 and 1.1 load in version 1.2. Ellipsoid and
rod states require RCPGenerator 1.2 or newer; older releases do not understand
their shape and orientation fields.

`to_dict()` preserves the principal histories and captured trajectory, but not
every advanced raw diagnostic listed above. Save the raw result separately if
those fields are required.

State flags are available as read-only properties:

- `initialized`: a realized initial state exists;
- `packed`: the current realized state has a completed packing result;
- `needs_initialize`: configuration changes require regeneration;
- `needs_pack`: the realized state has no current packing result.

## Diameter permutation

`diameter_permute=True` enables an optional move that permutes *which particle carries
which diameter* while the packing is being compressed. Positions are never moved by the
permutation and the multiset of relative diameters is preserved exactly (only the
assignment changes), so the requested distribution is what the final packing carries.
All continuous motion still comes from the packing algorithm, which is otherwise
unchanged. `diameter_permute=False` (the default) runs the original algorithm bit for bit.

```python
packing = rcpgenerator.Packing(..., diameter_permute=True)      # validated default recipe
packing = rcpgenerator.Packing(..., diameter_permute=True,
                               permute_budget=0.01)              # a gentler setting
result = packing.pack()
packing.permute_counts     # accepted exchanges per particle
packing.permute_passes     # passes run
packing.permute_exchanges  # accepted exchanges, pairwise + 3-cycle
```

**What it does.** From step `permute_start` on, every `permute_interval` steps, a pass
ranks the particles by size, scores each particle's *capacity* (the size its local cage
could hold) against the diameter it carries, and exchanges diameters between particles
within a 5 % rank window when the exchange lowers the overlap energy; 3-cycles
(A→B→C→A) are then tried to escape the pairwise local optimum. Each pass visits every
particle once, greedily, and stops after `permute_budget · N` exchanges (the 3-cycle
search has the same cap). Under-loaded particles grow and over-loaded ones shed size, so
the small particles fill the voids the large ones leave and the packing jams denser. It
works in two and three dimensions, in periodic boxes, with hard walls in any subset of
directions, and in curved containers.

**Knobs.**

| Argument | Default | Meaning |
| --- | --- | --- |
| `diameter_permute` | `False` | Turn the permutation on. Off is the original algorithm. |
| `permute_budget` | `0.25` | Cap on exchanges per pass as a fraction of N, in (0, 0.5]. The main dial: smaller = less uplift, on a smooth logarithmic scale (table below). |
| `permute_interval` | `5` | Steps between passes. Longer = fewer passes, less work, less uplift. |
| `permute_start` | `4000` | Step at which permutation begins. Later = less uplift; rarely worth changing. |

The defaults are the setting that reached the highest packing fractions in validation
(3-D power-law distributions with number density ∝ D^−3.45; φ is the packing fraction at
the end of the run with the residual overlap removed):

| size ratio S (N) | off | on (default) | gain |
| --- | --- | --- | --- |
| 40 (25 000) | 0.819 | 0.866 | +0.047 |
| 80 (242 481) | 0.849 | 0.903 | +0.054 |
| 100 (413 696) | 0.855 | 0.911 | +0.056 |

The budget dials smoothly from the original algorithm to the full gain (S = 40,
N = 25 000, one seed):

| `permute_budget` | off | 4×10⁻⁵ (1/N) | 10⁻⁴ | 10⁻³ | 10⁻² | 3×10⁻² | 0.1 | 0.25 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| φ | 0.819 | 0.821 | 0.822 | 0.840 | 0.859 | 0.863 | 0.866 | 0.866 |

The cap binds on every pass up to about 3 % of N and stops binding near 10 %, where the
gain saturates. The pass interval is the other rate knob, and wall time follows the
number of passes (same case):

| `permute_interval` | 1 | 2 | 5 (default) | 10 | 20 | 50 |
| --- | --- | --- | --- | --- | --- | --- |
| φ | 0.869 | 0.869 | 0.866 | 0.865 | 0.862 | 0.859 |
| wall time vs default | 3.3× | 1.8× | 1× | 0.7× | 0.57× | 0.48× |

`permute_start` moves the other way: 0.866 at 4000, 0.865 at 8000, 0.851 at 12000. The
permutation is deterministic for a given `seed` at one thread.

**Cost.** A pass costs a sort and O(N) capacity scores. At S = 80 to 100 a run with
permutation on takes roughly 3–4× the wall time of the run with it off, most of it in the
3-cycle search; `permute_budget=0.1` keeps nearly all of the gain at lower cost.

**What it cannot do.** Permutation reassigns diameters; it does not move particles. The
arrangement of the largest particles is set while the packing is still mobile and is not
improved by later exchanges, which is why the gain saturates about 0.04 below the density
an ideal largest-first (Apollonian-type) construction of the same size distribution
reaches at S = 80 to 100.

**Diagnostics.** After `pack()`, `packing.permute_counts` holds the accepted exchanges each
particle took part in (in the current particle order), `permute_passes` the number of
passes and `permute_exchanges` the total accepted exchanges. With permutation on, the raw
result's `spike_log` carries these per-particle counts. The recipe is assembled from
internal fields with a legacy `swap_` prefix; they are not part of the public API, and the
knobs above are the supported controls.

## Quench to force balance

`packing.quench()` takes the packing's current state, with or without diameter
permutation, and drives it to mechanical force balance at fixed diameters: every
particle's contact forces cancel. It is a FIRE minimisation (Bitzek et al., Phys. Rev.
Lett. 97, 170201, 2006) of the overlap energy; the diameters are frozen, the growth
schedule is suspended and no permutation runs.

```python
packing.pack()
start = packing.copy()                  # optional: keep the packed state
packing.quench()                        # compress 0.3 %, then balance to a residual of 1e-10
print(packing.quench_status)            # "converged", or why it stopped
print(packing.quench_message)

# tighten in stages: a quick first pass, then continue from where it stopped
packing = start.copy()
packing.quench(tol=1e-4)
packing.quench(tol=1e-10, resume=True)
```

**The small compression.** A packing run stops at the edge of jamming. At exactly its own
diameters the particles barely touch, a quench relaxes the remaining overlaps away, and
nothing is left to balance. `kappa_mult` multiplies every diameter once before the quench
(default 1.003, i.e. 0.3 % larger; 1.01 in 4D and above, where broad size distributions need
the larger compression to converge reliably), which gives a load-bearing contact network with small
overlaps. The quenched packing is therefore slightly compressed: in 3D its packing
fraction is about 0.9 % above the run's. `resume=True` continues from the previous
quench's state without rescaling again; combining it with `kappa_mult` is an error.

**Convergence measure.** Each particle's residual is |Σ f| / Σ |f| over its contacts (wall
contacts included): dimensionless and scale-free. `quench_residual` is the maximum over
the load-bearing particles (Ndim + 1 or more contacts and at least 1 % of the typical
contact load; below that, the ratio is limited by floating-point rounding, not by
balance). `quench_global_residual` is the RMS net force over the mean contact force.
`quench_log` holds, every `report_every` iterations, `[iteration, max residual, global
residual, energy, load-bearing count, time step]`.

**How it ends.** The quench stops as soon as one of these holds and reports which:

| `quench_status` | Meaning | What to try |
| --- | --- | --- |
| `"converged"` | Max residual ≤ `tol`. | — |
| `"max_iter"` | Iteration cap reached while still improving. | `quench(resume=True)` |
| `"stalled"` | The best residual has not halved within `stall_window` iterations. | Quench the packed state again with a lower `dt_max` or a slightly larger `kappa_mult` (see below); `quench(resume=True)` only if the residual was still falling. Below about 1e-8 the stop is the floating-point floor (roughly 6e-14 divided by the mean relative overlap), and the state is balanced as far as it can be. |
| `"unjammed"` | No load-bearing particles remain: the diameters are too small for the structure to carry load. | A larger `kappa_mult`. |
| `"step_collapse"` | FIRE's time step collapsed: it found no step it could follow downhill (typical of loose, under-compressed states). | Quench the packed state again with a lower `dt_max`, or a larger `kappa_mult`. |
| `"non_finite"` | A non-finite force or energy appeared. | Report it; do not use the state. |

Every status other than `"converged"` raises a `RuntimeWarning`, and `quench_message`
explains the outcome and the likely cause. Whatever the outcome, the packing keeps the
state the quench stopped in, and particle `i` keeps index `i` (positions, diameters,
`particle_origin` and `permute_counts` stay aligned).

**If a quench stalls or fails, retry it from the packed state** with either a lower `dt_max`
or slightly more compression (`kappa_mult`). `dt_max` is the largest time step the minimiser
may take (default 0.5 for spheres, 0.2 for ellipsoids and rods); `kappa_mult` is the factor
applied once to every diameter before the minimisation (default 1.003, and 1.01 for
ellipsoids, rods and 4D and above). A quench changes the packing in place, so keep a copy of
the packed state:

```python
import warnings
import rcpgenerator
packing = rcpgenerator.Packing(N=300, Ndim=2, phi=0.1, box=[1.0, 1.0], walls=[0, 0],
                               aspect_ratios=[0.5], seed=1)
packing.pack()
packed = packing.copy()                # the packed state, before the first quench
with warnings.catch_warnings():
    warnings.simplefilter("ignore", RuntimeWarning)   # a quench that does not converge warns
    packing.quench()
    if packing.quench_status != "converged":
        packing = packed.copy()
        packing.quench(dt_max=0.1)     # or: packing.quench(kappa_mult=1.02)
print(packing.quench_status, packing.quench_iters)
```

An example from the development runs: 400 ellipses of elongation 4 (`aspect_ratios=[0.25]`),
whose quench at the defaults (`kappa_mult=1.01`, `dt_max=0.2`) was still at a residual of 0.4
after 200,000 iterations.

| Retry of the packed state | Outcome | Iterations |
| --- | --- | --- |
| `quench(dt_max=0.1)` | converged | 35,248 |
| `quench(dt_max=0.05)` | converged | 157,347 |
| `quench(kappa_mult=1.02)` | converged | 15,830 |
| default settings, then `quench(resume=True)` | not converged | 200,000 more |

What to take from it:

- **Either remedy works.** Half the time step or 1 % more compression each brought the packing
  to a residual of 1e-10. More compression changes the state that is balanced (the particles
  are larger); a lower `dt_max` does not.
- **Do not go much lower on `dt_max` than needed.** Going from 0.1 to 0.05 cost 4.5 times the
  iterations. Halve it once, and again only if that fails.
- **Resuming with unchanged settings does not help a quench that is stuck.** `resume=True` is
  for a run that was still improving when it stopped (status `"max_iter"`, or `"stalled"`
  while the residual was still falling slowly).
- **When a lower `dt_max` does not help, more compression is the next thing to try** (second
  example below).

A second example: **flat-ended rods in 3D inside a curved container often need a retry.** A
flat end standing on a curved wall is held by many points of its rim at once, a much stiffer
contact than any other in the packing. In a survey of 150 flat-ended rods (`rod_length=1.5`,
three seeds per container):

| Container | Converged at the defaults | After a retry with `dt_max=0.1` | The remaining case |
| --- | --- | --- | --- |
| periodic box | 3 of 3 | | |
| sphere | 1 of 3 | 3 of 3 | |
| cylinder | 2 of 3 | 2 of 3 | `dt_max=0.05`: not converged after 200,000 iterations; `kappa_mult=1.02`: converged in 10,275 |

Ellipsoids, capped rods and 2D rectangles in curved containers converged at the defaults in
the cases tested.

The defaults for ellipsoids and rods are already gentler than for spheres for this reason.

**Where it works.** On 3-D power-law packings with the default settings the quench converges
to 1e-10 for size ratios up to about 10 (N = 4000: 40 000–56 000 iterations, about 1 ms per
iteration per 4000 particles at one thread). It also converges in 2D, with hard walls in
one, two or three directions, and in disk, cylinder and sphere containers (N = 500, size
ratio 5). The iteration count grows with N, and steeply for compressions close to the
jamming point. A small fraction of runs stop early as `"stalled"` while still converging,
and `quench(resume=True)` usually finishes them. Each new seed or machine is a different
realisation and can land either way. With N = 500 and size ratio 5 this happened for 2 of
40 seeds in 3D periodic boxes at the default compression, and for 2 of 660 runs across all
2D/3D geometries at `kappa_mult=1.01`. In 4D with hard walls and size ratio 5 it is about 1
call in 6 (33 of 40 seeds converged in one call; 30 of 30 at size ratio 2).
Above a size ratio of about 20 the method usually cannot
balance the packing as it is: most particles are rattlers, the load is carried through
small particles with only two or three contacts, and the quench converges only after a
substantial compression (S = 25 converges with `kappa_mult=1.01`; at S = 40 even the
smallest compression that can be balanced, about 2.5 % in diameter, raises the packing
fraction from 0.79 to 0.85). That is a compressed packing, not force balance of the
original; at the default setting the quench gives up with a message rather than run on. Gradient,
quasi-Newton, Newton and active-set minimisers were also tried on these states and none
finished faster than FIRE: the last stage of convergence is a sequence of contact changes,
which only the dynamics resolve efficiently.

**Holding particles.** `held` lists indices into `positions` of particles that must not
move. Their contacts still act on the others, and `quench_held_force` is the total force
on all held particles. A shear step is: displace a slab of particles, hold it and a second
slab, quench, and sum the contact forces the free particles exert on the displaced slab
(a contact's force is its overlap along the line of centres; `examples/example_3d_quench.py`
does this); the stress is that force over the slab's area. Keep the two slabs half a box
apart so that they do not touch through a periodic boundary. On the validation cases, three steps of 10⁻⁵ L gave
forces in the ratio 1 : 2.00 : 2.99, equal and opposite on the two slabs to five digits.

**Settings.**

| Argument | Default | Meaning |
| --- | --- | --- |
| `tol` | `1e-10` | Stop when the max residual is at or below this. |
| `max_iter` | `200000` | Iteration cap for this call. |
| `kappa_mult` | `1.003` (2D/3D), `1.01` (4D+) | One-time rescaling of all diameters before the quench. Not applied with `resume=True`. |
| `resume` | `False` | Continue from the previous quench's state (no rescaling). |
| `held` | `None` | Indices of particles that must not move. |
| `stall_window` | `max(20000, 20·N)` | Give up when the best residual has not halved in this many iterations; `0` disables. |
| `verbose` | packing's `verbose` | Print the outcome. |
| `mass_scaled`, `dt0`, `dt_max`, `alpha0`, `mix_mode`, `report_every` | `True`, `0.05`, `0.5`, `0.03`, `1`, `500` | Advanced FIRE settings; the defaults are the measured best. `mix_mode=0` uses global velocity mixing as in the original FIRE, `1` per-particle mixing. |

With ellipsoids and rods the defaults of `kappa_mult` and `dt_max` are `1.01` and `0.2`, and rods take
one more argument, `rod_round`; see [Quench](#quench-with-ellipsoids-and-rods) under Ellipsoids and rods.

## Ellipsoids and rods

The same algorithm packs particles that are not spheres: ellipses and ellipsoids, and rods with round
or flat ends. Every particle then also has an orientation, which is relaxed together with its position
while the particles grow. The shape is chosen by constructor keywords: `aspect_ratios` for ellipsoids,
or `rod_length` and `rod_radius` for rods (one or the other, never both). Without these keywords
the sphere code path is the one of 1.1.0; results differ from 1.1.0 only where the two repairs listed in
`CHANGELOG.md` act (when a run stops, and the step size in curved containers).

Everything else on this page still applies, with the exceptions listed under
[What works with them and what refuses them](#what-works-with-them-and-what-refuses-them).
`rcpgenerator.example()` has a group of snippets for each shape.

### Ellipsoids

```python
import numpy as np
import rcpgenerator

packing = rcpgenerator.Packing(
    N=300,
    Ndim=2,
    box=[1.0, 1.0],
    walls=[0, 0],
    aspect_ratios=[0.5],       # second axis = 0.5 D: 2:1 ellipses
    seed=1,
)
packing.pack()

print("phi at contact:", packing.phi_contact)
orientations = np.asarray(packing.orientations)   # shape (N, Ndim, Ndim)
print("long axis of particle 0:", orientations[0][:, 0])
packing.render(path="ellipses.png")
```

**`aspect_ratios` and what `D` means.** `D`, the diameter drawn from `dist`, is the full length of a
particle along its first body axis, and that axis is its **longest**. `aspect_ratios` lists the other
`Ndim - 1` axes in units of `D`: axis `k` has length `D × ratio_k`.

| `Ndim` | `aspect_ratios` | Axes | Shape |
| --- | --- | --- | --- |
| 2 | `[0.5]` | `D`, `0.5 D` | 2:1 ellipse |
| 3 | `[0.5, 0.5]` | `D`, `0.5 D`, `0.5 D` | Prolate (one long axis), 2:1 |
| 3 | `[1, 0.5]` | `D`, `D`, `0.5 D` | Oblate (one short axis), 2:1 |
| 3 | `[0.8, 0.6]` | `D`, `0.8 D`, `0.6 D` | Three different axes |

**The decreasing-order rule.** The ratios are written longest first and none exceeds 1:
`1 ≥ ratio_2 ≥ ratio_3 ≥ … > 0`. Every shape then has exactly one way of being written. Input that
breaks the rule raises `ValueError`; it is not reordered for you, because that would rescale `D` and
change what the size distribution means. To convert a shape written another way, divide `D` by the
largest ratio and list the ratios longest first: axes `D × (1, 2)` become `aspect_ratios=[0.5]` with
`D` twice as large. The size distribution is therefore the distribution of longest-axis lengths, and a
particle's volume is the volume of a sphere of diameter `D` times the product of its ratios.

**One list, or one per particle.** A single list applies to every particle. A list of `N` lists gives
each particle its own shape; every row obeys the order rule. Which particle gets which row depends on where the diameters
come from:

- **Diameters generated by the package** (any `dist` but `"custom"`): row `i` goes to the `i`-th
  **largest** particle, because the generated particles are sorted by decreasing diameter before the
  shapes are attached. A table drawn independently of size is thus a random assignment.
- **Your own diameter list** (`dist={"type": "custom", "custom": [...]}`): row `i` belongs to entry
  `i` of your diameter list and stays with it through the sort. This is the way to tie a shape to a
  size.

After `pack()`, `packing.aspect_ratios` is the table in
the packing's own particle order: row `k` is the shape of the particle at `positions[k]` with
`diameters[k]`.

```python
import rcpgenerator

N = 250
shapes = [[0.8, 0.6]] * 125 + [[0.5, 0.5]] * 125    # one row per particle
packing = rcpgenerator.Packing(
    N=N, Ndim=3, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
    dist={"type": "lognormal", "mu": 0.0, "sigma": 0.2},
    aspect_ratios=shapes,      # the 125 largest particles get [0.8, 0.6]
    seed=1,
)
packing.pack()
print("phi at contact:", packing.phi_contact)
print("particle 0:", packing.diameters[0], packing.aspect_ratios[0])
```

**Orientations in and out.** `packing.orientations` holds one `Ndim × Ndim` rotation matrix per
particle, in the same particle order as `positions`. The columns of a matrix are the particle's body
axes in the box frame, so column 0 is the direction of its longest axis. The initializer generates
random orientations, as it generates positions, and `pack()` returns the final ones. To start from
your own, pass `orientations=` to the constructor: one matrix per particle (nested, or flat as
`Ndim × Ndim` numbers row by row), in the order of `packing.diameters` after construction (largest
first). Orientations without `aspect_ratios` (or the rod keywords) raise `ValueError`: supply a shape
and orientations, or a shape alone, or neither.

```python
import numpy as np
import rcpgenerator

N, angle = 200, np.radians(30.0)
R = [[np.cos(angle), -np.sin(angle)],
     [np.sin(angle), np.cos(angle)]]          # columns = body axes; column 0 = long axis
packing = rcpgenerator.Packing(
    N=N, Ndim=2, box=[1.0, 1.0], walls=[0, 0],
    aspect_ratios=[0.5],
    orientations=[R] * N,                     # every ellipse starts at 30 degrees
    seed=1,
)
packing.pack()                                # orientations evolve with the positions
long_axis = np.asarray(packing.orientations)[:, :, 0]
degrees = np.degrees(np.arctan2(long_axis[:, 1], long_axis[:, 0])) % 180.0
print("final angles of the first three particles:", degrees[:3].round(1))
```

### Rods

Two kinds of rod, chosen with `rod_capped`, one value for the whole packing:

| `rod_capped` | Shape | What `rod_length` measures |
| --- | --- | --- |
| `True` (default) | A straight section with a half-sphere on each end (a spherocylinder; in 2D a stadium). | The **straight section only**. The rod is `rod_length + 2 × rod_radius` from tip to tip, and a capped rod of length 0 is exactly a sphere. |
| `False` | Flat ends: a cylinder (in 2D a rectangle). 2D and 3D only. | The **full length** of the rod. It must be positive. |

**How the two numbers combine with `D`.** Both `rod_length` and `rod_radius` are in units of the
particle's `D` from the size distribution: the straight length is `D × rod_length` and the radius is
`D × rod_radius`. With `rod_radius=0.5`, `D` is the rod's thickness and `rod_length` is its length in
thicknesses; a capped rod with `rod_length=2` is then 3 thicknesses from tip to tip. The particles
grow by scaling `D`, so a rod keeps its shape while it grows. Both keywords are required together, and
there is no ordering rule.

Flat-ended rods stop with larger leftover forces than other shapes; if force balance matters, use
[`quench()`](#quench-to-force-balance) on the packed state.

**Per-particle lists.** Each of `rod_length` and `rod_radius` is one number for every particle or a
list of `N`. As for ellipsoids, entry `i` goes to the `i`-th largest particle when the package generates the
diameters, and with entry `i` of your own diameter list when you supply one (`dist` of type
`"custom"`); after `pack()`
`packing.rod_length` and `packing.rod_radius` are in the packing's own particle order.

**Orientations** work exactly as for ellipsoids. Only the direction of the rod's axis matters; it is
the first body axis, column 0 of the particle's matrix.

Capped rods:

```python
import rcpgenerator

packing = rcpgenerator.Packing(
    N=250, Ndim=3, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
    rod_length=1.5,            # straight section: 1.5 thicknesses, 2.5 from tip to tip
    rod_radius=0.5,            # radius 0.5 D, so D is the thickness
    rod_capped=True,           # half-sphere ends
    seed=1,
)
packing.pack()
print("phi at contact:", packing.phi_contact)
```

Flat-ended rods, here 2 × 1 rectangles between hard walls:

```python
import rcpgenerator

packing = rcpgenerator.Packing(
    N=300, Ndim=2, box=[1.0, 1.0],
    walls=[0, 1],              # periodic in x, hard walls in y
    rod_length=2.0,            # full length: 2 thicknesses
    rod_radius=0.5,
    rod_capped=False,          # flat ends
    seed=2,
)
packing.pack()
print("phi at contact:", packing.phi_contact)
packing.render(path="rectangles.png")
```

Two lengths in one packing:

```python
import rcpgenerator

N = 300
packing = rcpgenerator.Packing(
    N=N, Ndim=2, box=[1.0, 1.0], walls=[0, 0],
    rod_length=[0.5] * (N // 2) + [3.0] * (N - N // 2),   # one entry per particle
    rod_radius=0.5,
    rod_capped=True,
    seed=1,
)
packing.pack()
print("phi at contact:", packing.phi_contact)
print("particle 0 has straight length", packing.rod_length[0], "x D")
```

### Reading the result: `phi_contact` and `phi_final`

An ellipsoid or rod packing carries two packing fractions:

| Attribute | Meaning |
| --- | --- |
| `phi_final` | The raw value at the end of the run, as for spheres. The particles still overlap slightly. |
| `phi_contact` | The packing fraction after every particle has been shrunk by one common factor to exact contact, so that no pair and no wall overlaps. **This is the value to report.** |
| `shrink_factor` | That common factor (at most 1). Multiply every diameter by it to obtain the non-overlapping configuration. |

Both fractions use the true particle volume (ellipsoid or rod, not the bounding sphere).
`pack()` fills `phi_contact` and `shrink_factor` for ellipsoids and rods; for spheres they stay
`None`.

`packing.shrink_to_contact()` computes the same shrink on the current state at any time, for spheres
as well, in periodic boxes, with flat walls and in curved containers. It returns a dictionary with `scale` (the factor, with a
safety margin of 1e-12), `phi` (the packing fraction at that scale), `n_overlapping` (contacts that
overlapped before the shrink) and `pair` (the tightest pair). It does not modify the packing.

```python
import numpy as np
import rcpgenerator

packing = rcpgenerator.Packing(N=300, Ndim=2, box=[1.0, 1.0], walls=[0, 0],
                               aspect_ratios=[0.5], seed=1)
packing.pack()
print("raw:", packing.phi_final, " at contact:", packing.phi_contact)
contact_diameters = np.asarray(packing.diameters) * packing.shrink_factor   # no overlap left

spheres = rcpgenerator.Packing(N=300, Ndim=2, box=[1.0, 1.0], walls=[0, 0], seed=1)
spheres.pack()
shrunk = spheres.shrink_to_contact()          # works for spheres too
print("spheres raw:", spheres.phi_final, " at contact:", shrunk["phi"], " factor:", shrunk["scale"])
```

The shrink is exact for ellipsoids, capped rods and 2D flat-ended rods. For 3D flat-ended rods the
factor is a lower bound, so `phi_contact` can be slightly too low but the shrunk configuration never
overlaps.

### What works with them and what refuses them

Ellipsoids and rods work wherever spheres do: every boundary and container, `fix_height`, permutation,
the quench, rendering, the viewer and animation. The few things that are not extended raise an error;
nothing runs on the assumption that the particles are spheres. They are listed under
[What still refuses them](#what-still-refuses-them).

| Feature | Ellipsoids | Rods |
| --- | --- | --- |
| Periodic boxes; hard walls on any subset of axes (`walls` entries `0`, `1`) | Yes | Yes |
| Diameter distributions (`dist`) | Yes: the distribution of `D` | Yes: the distribution of `D` |
| Dimensions | 2 to 8 | Capped: 2 to 8. Flat-ended: 2 and 3; above, `NotImplementedError` |
| Curved containers (negative first `walls` entry): disk, cylinder, sphere, and their higher-dimensional forms | Yes, in every dimension | Capped: every dimension. Flat-ended: 2D and 3D |
| `fix_height=True` | Yes: the height is a multiple of the largest `D` | Yes: a multiple of the largest rod's full extent (tip to tip, or the diagonal of a flat-ended rod) |
| `diameter_permute=True` | Yes | Yes |
| `quench()` | Yes | 2D and 3D; above, `NotImplementedError` |
| `copy()`, `to_dict()`, `from_dict()` | Yes: shapes and orientations are carried | Yes |
| `render()`, `savefig()`, `show_packing()` | 2D and 3D, true axes and orientations | 2D and 3D |
| `animate_2d()` | Yes | Yes |
| Viewer bundle (`write_bundle`) | 2D and 3D, one shape for all particles; a per-particle table raises `NotImplementedError` | 2D and 3D, per-particle lists included |

**Curved containers** use the same `walls` values as for spheres (`[-2, 0]`: a disk; `[-2, 0, 0]`: a
cylinder, with `1` instead of the last `0` for end walls; `[-3, 0, 0]`: a sphere), and `phi_contact` is
the particle volume over the container's own volume. A particle touches the wall where its farthest
point from the container's axis reaches it.

```python
import rcpgenerator

# oblate ellipsoids in a sphere
packing = rcpgenerator.Packing(N=200, Ndim=3, box=[1.0, 1.0, 1.0], walls=[-3, 0, 0],
                               aspect_ratios=[1.0, 0.5], seed=1)
packing.pack()
print("ellipsoids in a sphere:   phi =", round(packing.phi_contact, 4))

# flat-ended rods (rectangles) in a disk
packing = rcpgenerator.Packing(N=200, Ndim=2, box=[1.0, 1.0], walls=[-2, 0],
                               rod_length=2.0, rod_radius=0.5, rod_capped=False, seed=1)
packing.pack()
print("rectangles in a disk:     phi =", round(packing.phi_contact, 4))

# capped rods in a slab whose height stays 3 rods (tip to tip) tall
packing = rcpgenerator.Packing(N=200, Ndim=3, box=[1.0, 1.0, 3.0], walls=[0, 0, 1], fix_height=True,
                               rod_length=1.0, rod_radius=0.5, rod_capped=True, seed=1)
packing.pack()
print("capped rods, fix_height:  phi =", round(packing.phi_contact, 4))
```

Examples at N = 300 (one seed): 2:1 ellipses in a disk 0.871, 2 x 1 rectangles in a disk 0.870, capped
rods in a disk 0.879; prolate ellipsoids in a cylinder 0.66, oblate ellipsoids in a sphere 0.65. As for
spheres, the fraction in a container is lower than in a periodic box of the same N, because of the layer
at the wall.

#### What still refuses them

| What | Ellipsoids | Rods | Why |
| --- | --- | --- | --- |
| Flat-ended rods above 3D | – | `NotImplementedError` | The contact of two flat-ended cylinders is worked out for 2D and 3D only. Capped rods work in every dimension. |
| `quench()` of rods above 3D | works | `NotImplementedError` | Rods are balanced as slightly rounded rods, and that contact is written for 2D and 3D. |
| More than 8 dimensions | error | error | The fixed-size arrays of the shape code (the sphere kernel has the same limit). |
| Viewer bundle with a per-particle `aspect_ratios` table | `NotImplementedError` | lists work | The viewer draws one ellipsoid shape for all particles; `render()` draws tables. |
| Viewer bundle above 3D | as spheres | `NotImplementedError` | The viewer shows 2D packings and slices of 3D ones. |
| `animate_2d()` with a shape per particle and `diameter_permute=True` | `NotImplementedError` | `NotImplementedError` | Exchanges move shapes between particles and the frames do not record which shape sat where. |
| Research inputs of the core: the distribution morph (`morph_target`), fixed particles, the internal `swap_*` settings and swap injection, the quench's two-level, activity-damping, warm-up and per-particle-reset variants | error | error | Undocumented research features that work from sphere quantities; only the public `diameter_permute` and `quench()` are extended. |

Not a refusal, but a limit to know: **the quench of 3D flat-ended rods in a curved container often needs
a retry** at the default settings (lower `dt_max`, then slightly more compression); see
[Quench to force balance](#quench-to-force-balance).

#### Permutation with ellipsoids and rods

`diameter_permute=True` works with both shapes, with the same controls as for spheres
(`permute_budget`, `permute_interval`, `permute_start`; see
[Diameter permutation](#diameter-permutation)). An exchange trades the **diameter and the shape** of
two particles; positions are not moved. By default the **orientations stay in place**
(`ell_permute_orient=0`): each site keeps its orientation and receives the other particle's size and
shape. `ell_permute_orient=1` makes the orientation travel with the particle, so that the two whole
particles trade places. The keyword has this name for rods too.

Two particles may be exchanged only when they are alike in every direction: every axis length within
5 % for ellipsoids, straight length and radius each within 5 % for rods. Shapes further apart than
that are never exchanged with each other, so in a mixture of two well-separated shapes each shape is
permuted on its own. There is no 3-cycle pass for these shapes.

```python
import rcpgenerator

dist = {"type": "powerlaw", "d_min": 1.0, "d_max": 5.0, "exponent": -3.0}
for permute in (False, True):
    packing = rcpgenerator.Packing(
        N=300, Ndim=2, box=[1.0, 1.0], walls=[0, 0], dist=dist, seed=1,
        aspect_ratios=[0.8],
        diameter_permute=permute,      # add ell_permute_orient=1 to make orientations travel
    )
    packing.pack()
    print(f"permute {permute}: phi at contact = {packing.phi_contact:.4f}, "
          f"exchanges = {packing.permute_exchanges}")
```

Measured gain in `phi_contact` on power-law size distributions, permutation on against off, same
seeds:

| Shape | Orientations in place (default) | Orientations travel (`ell_permute_orient=1`) |
| --- | --- | --- |
| Ellipses and ellipsoids | +0.007 to +0.025 (2D), +0.007 to +0.018 (3D) | 0.000 to +0.002 for `aspect_ratios` of 0.8 and below |
| Capped rods | +0.012 to +0.029 | No gain |
| Flat-ended rods | +0.023 to +0.071 (2D rectangles most) | No gain, except +0.006 to +0.023 for 2D rectangles |

For ellipsoids the default gains a quarter to a half of what spheres gain on the same sizes, most for
shapes nearest a sphere. A mixture of two well-separated ellipsoid shapes gained no more than its
less favourable shape alone, and a mixture of two rod lengths no more than its better length alone. The cost is higher than for spheres: at N = 1000 on one thread a run with permutation took 3.6
to 4.6 times as long as the run without it for ellipsoids with orientations in place (4.5 to 6.3 times
with travelling orientations), against 1.7 to 2.0 times for spheres. These timings are from two seeds
and the pass has not been optimised.

#### Quench with ellipsoids and rods

`packing.quench()` moves positions and orientations together and brings the packing to force **and
torque** balance. The residual of a particle combines its net force and its net torque (divided by
the particle's bounding radius) relative to the sum of its contact forces. Statuses, `resume`, `held` and the
other arguments are as described under [Quench to force balance](#quench-to-force-balance). Two
defaults differ from spheres:

| Argument | Spheres | Ellipsoids and rods | Why |
| --- | --- | --- | --- |
| `kappa_mult` | `1.003` (2D/3D) | `1.01` | At 0.3 % compression the quench did not converge for the more elongated or flattened ellipsoids (2D at 4:1, 3D prolate and oblate at 3:1, and the three-axis shape `[0.8, 0.64]`); at 1 % all converged. |
| `dt_max` | `0.5` | `0.2` | At 0.5 the time step keeps overshooting and being cut back on elongated or flattened shapes: in an earlier survey 50 of 56 cases converged and 2D ellipses at 4:1 never did, where 0.2 brought all 56 to balance. With an exact contact solution in both cases, 0.5 still failed on 4 of 6 of the hardest cases and took 7 to 29 times the iterations on the other 2. |

**If a quench of ellipsoids or rods stalls or fails, quench the packed state again with a lower
`dt_max` (0.1, say) or slightly more compression (`kappa_mult=1.02`)**, as shown with measured
numbers under [Quench to force balance](#quench-to-force-balance). Resuming with unchanged settings
does not help, and a `dt_max` much lower than needed costs many more iterations. Keep a copy of the packed state for this,
because a quench changes the packing in place.

```python
import warnings
import rcpgenerator

packing = rcpgenerator.Packing(N=200, Ndim=2, box=[1.0, 1.0], walls=[0, 0],
                               aspect_ratios=[0.5], seed=1)
packing.pack()
packed = packing.copy()                # keep the packed state
with warnings.catch_warnings():
    warnings.simplefilter("ignore", RuntimeWarning)
    packing.quench()                   # kappa_mult=1.01, dt_max=0.2 for this shape
    if packing.quench_status != "converged":
        packing = packed.copy()        # start again from the packed state
        packing.quench(dt_max=0.1)     # retry the packed state: lower dt_max, or kappa_mult=1.02
print(packing.quench_status, packing.quench_residual, packing.quench_iters)
```

**Rods are quenched as slightly rounded rods.** Where two flat sides or ends rest on each other, the
contact energy of a sharp rod has a kink and a minimiser finds no balance point: trial quenches of
sharp rods stalled in 4 of 4 cases. For the quench each rod is therefore replaced by a rod whose edges
are rounded by 1 % (`rod_round=0.01`, an argument of `quench()` for rods). The balanced state is that
of the rounded rods at 1 % compression. The packing itself is always made with the sharp shape.
`rod_round=0` quenches the sharp rods and usually stalls. The rounded contact exists in 2D and 3D, so
rods cannot be quenched in 4D and above.

**Where it has been run.** With the defaults, 52 of 53 ellipsoid cases converged to a residual of
1e-10 (the one that did not, 2D ellipses at 4:1, is the retry example under
[Quench to force balance](#quench-to-force-balance)): 2D ellipses from nearly round to 4:1, 3D prolate and oblate shapes to 3:1 and a three-axis
shape, hard walls, two size distributions, 4D, and N up to 2000, with two or three seeds each. The
median was 14,100 iterations, against roughly 2,000 to 5,000 for spheres. For rods, 78 of 78 cases
converged to 1e-10: 2D and 3D, capped and flat-ended, `rod_length` from 0.25 to 4, lognormal sizes,
hard walls and mixed lengths. More elongated shapes than these have not been tried. The quenched state
is not a small perturbation of the packed one: at 1 % compression the median particle moved 0.1 to
0.2 `D` in the ellipsoid cases (N = 400).

#### Rendering, the viewer and animation

`render()`, `savefig()` and `show_packing()` draw ellipses, ellipsoids and rods with their true axes
and orientations, in 2D and 3D, including a shape per particle. As for spheres, 2D pictures include
the periodic images and 3D pictures do not. Drawing a 3D packing takes much longer than packing it.

`animate_2d()` draws ellipses and rods at their recorded positions, sizes and orientations; capture
the trajectory as for discs (see [Rendering and animation](#rendering-and-animation)).
`trajectory_orientations[f][k]` is the rotation matrix of the particle in slot `k` of frame `f`. One
combination is refused: a shape per particle (an `aspect_ratios` table, or a list for `rod_length` or
`rod_radius`) together with `diameter_permute=True` raises `NotImplementedError`, because exchanges
move shapes between particles and the frames do not record which shape sat where.

`rcptools.bridge.write_bundle` writes the shape and the orientations into the viewer bundle. The
viewer then shows a 2D packing directly and a 3D packing through its slice, where each particle
appears as the exact outline that the slicing plane cuts from it. Pass the shape and the orientations
together:

```python
import rcpgenerator
from rcptools.bridge import write_bundle

ellipses = rcpgenerator.Packing(N=80, Ndim=2, box=[1.0, 1.0], walls=[0, 0],
                                aspect_ratios=[0.5], seed=1)
ellipses.pack(capture_trajectory=True, trajectory_spacing="log", trajectory_per_decade=8)
ellipses.render(path="ellipses.png")
ellipses.animate_2d(path="ellipses.gif", show=False)
write_bundle("my_ellipses", ellipses.positions, ellipses.diameters, ellipses.box,
             walls=ellipses.walls,
             aspect_ratios=ellipses.aspect_ratios, orientations=ellipses.orientations)

rods = rcpgenerator.Packing(N=80, Ndim=2, box=[1.0, 1.0], walls=[0, 0],
                            rod_length=2.0, rod_radius=0.5, rod_capped=True, seed=1)
rods.pack()
write_bundle("my_rods", rods.positions, rods.diameters, rods.box, walls=rods.walls,
             rod_length=rods.rod_length, rod_radius=rods.rod_radius,
             rod_capped=rods.rod_capped, orientations=rods.orientations)
```

Leaving out the shape arguments and the orientations writes an ordinary sphere bundle of diameter `D`
without an error, so pass them whenever the packing has a shape. The viewer draws each ellipse as a
polygon, one object per particle; this has been used up to N = 2500 and is not meant for the
hundreds of thousands of particles that the disc display handles.

### Cost and limits

All numbers in this section are wall times of a complete `pack()` on one thread against a sphere
packing with the same N, dimension and size distribution. *Elongation* is the longest axis divided by
the shortest.

**Ellipsoids** (N = 1000, equal sizes, one seed):

| Shape | `aspect_ratios` | Time against spheres | `phi_contact` |
| --- | --- | --- | --- |
| 2D, elongation 1.25 to 2 | `[0.8]` to `[0.5]` | 2.1 to 2.3 times | 0.881 to 0.892 |
| 2D, elongation 3 | `[0.333]` | 2.7 times | 0.858 |
| 2D, elongation 4 | `[0.25]` | 3.2 times | 0.841 |
| 2D, elongation 6 | `[0.167]` | 3.8 times | 0.802 |
| 3D prolate, elongation 1.25 to 2 | `[0.8, 0.8]` to `[0.5, 0.5]` | 2.4 to 3.7 times | 0.695 to 0.709 |
| 3D prolate, elongation 3 | `[0.333, 0.333]` | 6 times | 0.661 |
| 3D prolate, elongation 4 | `[0.25, 0.25]` | 9 times | 0.621 |
| 3D prolate, elongation 8 | `[0.125, 0.125]` | 23 times | 0.467 |
| 3D oblate, elongation 1.25 to 2 | `[1, 0.8]` to `[1, 0.5]` | 2.2 to 2.7 times | 0.692 to 0.705 |
| 3D oblate, elongation 4 | `[1, 0.25]` | 3.5 times | 0.623 |
| 3D, three axes | `[0.8, 0.64]` | 2.7 times | 0.725 |

Up to an elongation of 2 an ellipsoid packing takes 2.1 to 3.7 times as long as the sphere packing
(3.0 times at N = 4000 on four threads, for prolate ellipsoids of elongation 1.5). Beyond that the
cost is set by the neighbour search, which works on the particles' bounding spheres of diameter `D`:
the number of pairs whose bounding spheres overlap grows roughly as the square of the elongation for
prolate shapes, and the time follows it. Oblate shapes are much cheaper than prolate ones at the same
elongation. Peak memory was 100 to 160 MB at N = 1000 whatever the shape. The packing fraction falls
smoothly as the elongation grows.

Where ellipsoids stop working (one seed each):

- **3D prolate:** elongation 8 packs; elongation 12 fails with the `exceeded its neighbor capacity`
  error, because a particle's bounding sphere then has more neighbours than its row holds. A larger
  `neighbor_max` was not tried.
- **2D:** elongation 12 packs; elongation 20 did not reach a jammed packing within the step limit.
- **Oblate** shapes pack down to at least `[1, 0.08]`. Size ratios of 50 in 2D and 10 in 3D pack
  (elongation 1.5), as do 4D and 5D ellipsoids. The shape code has a hard limit of 8 dimensions.

**Rods** (N = 1000, `rod_radius=0.5`, equal sizes, mean of two seeds):

| `rod_length` | 2D capped | 2D flat-ended | 3D capped | 3D flat-ended |
| --- | --- | --- | --- | --- |
| 0.5 | 2.1 times | — | 2.3 times | — |
| 1 | 2.3 times | 2.6 times | 3.2 times | 7.0 times |
| 2 | 3.0 times | 2.7 times | 6.2 times | 6.7 times |
| 4 | 4.0 times | 3.4 times | 12 times | 9.5 times |

Short rods, up to 2 thicknesses from tip to tip, take 2.1 to 2.7 times as long as spheres in 2D and
2.3 to 3.2 times for capped rods in 3D. **3D flat-ended rods are the slowest shape: 6 to 7 times the
sphere time even when short.** Two flat faces resting on each other have no closed-form contact, so
each such contact is solved iteratively inside every step. As for ellipsoids, the time for longer rods
follows the number of neighbours inside the bounding spheres: for 3D capped rods at N = 500, 12
neighbours per rod and 5 s at `rod_length=0.5`, 384 neighbours and 79 s at `rod_length=16`, against
2 s for spheres. Memory stayed between 78 and 125 MB in those runs.

Packing fractions at contact from single runs at N = 1000 (`rod_radius=0.5`): 3D capped 0.708, 0.696,
0.666 and 0.594 at `rod_length` 0.5, 1, 2 and 4; 3D flat-ended 0.684, 0.659 and 0.583 at 1, 2 and 4;
2D capped 0.894 to 0.900, with stacks of parallel rods; 2D flat-ended 0.936, 0.887 and 0.855 at 1, 2
and 4.

Where rods stop working (one seed each, N = 500):

- **Rods longer than half the box.** A rod that long can meet a periodic copy of its neighbour or of
  itself, which the forces do not account for (`rod_length=16` at N = 500 is such a case). The
  `periodic box is small` warning takes the rod's full extent into account; use more particles.
- **3D flat-ended rods at `rod_length=16`** ended loose, with 0.6 contacts per rod.
- **Capped rods in 5D** failed with the `exceeded its neighbor capacity` error at the automatic
  setting; 4D packs.

## Neighbor capacity and memory

The pair list uses a preallocated row for each particle. `neighbor_max` is the
base term in that capacity estimate; it is not a strict global maximum number
of physical neighbors. With particle diameter `D_i`, minimum diameter `D_min`,
and dimension `d`, the implementation allocates approximately

```text
ceil((d / 3) * (50 * (D_i / D_min) ** (1 + d / 6) + base))
```

slots for particle `i`, capped at `N`. With `neighbor_max=0`, the automatic
base is dimension-specific (300 in 2D, 750 in 3D). A larger particle-size ratio
can therefore increase memory use substantially even when the automatic setting
succeeds.

In 4D and higher with `neighbor_max=0`, each particle instead gets four times the
number of neighbours expected inside its search shell in a jammed packing, plus
100 (about 250 slots per particle in 4D, 400 in 5D, 700 in 6D for equal sizes;
larger particles get proportionally more). The number of neighbours roughly
doubles with each added dimension, so this keeps memory proportional to what the
packing needs: about 3 GB for two million particles in 5D. Passing `neighbor_max`
explicitly restores the formula above in any dimension.

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
with the renderer's default figure size and DPI. Ellipsoids and rods are drawn with their
true shapes and orientations, by the same calls; see
[Rendering, the viewer and animation](#rendering-the-viewer-and-animation).

To animate 2D relaxation, capture a trajectory during the run:

```python
packing.pack(capture_trajectory=True, trajectory_spacing="log")   # frames evenly in log(step)
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

Almost all of the motion happens in the first few hundred of some 16,000 steps.
`trajectory_spacing="log"` records frames evenly in the logarithm of the step
(`trajectory_per_decade` frames per factor of 10, default 20, about 85 frames in all), so
the movie shows the particles spreading and growing throughout. The default,
`trajectory_spacing="linear"`, records one frame every `trajectory_interval` steps, which
is nearly still after the first frames.

The packer renumbers particles while it runs, so `trajectory_positions[f][k]` is not the
same particle in every frame `f`. `trajectory_origin[f][k]` gives the identity (the index
at the start of the run) of the particle in slot `k` of frame `f`; `animate_2d()` uses it
so that a particle keeps its colour for the whole movie.

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

For ellipsoids, add `orientations` to `input` (one flat list of `Ndim × Ndim` numbers per particle, a
rotation matrix row by row) and `aspect_ratios` to `config`. For rods, add `orientations`, `rod_length`
and `rod_radius` (lists of one number or of `N`) to `input` and `rod_capped` to `config`. At this level
nothing generates orientations for you: a shape without orientations, or orientations without a shape,
raises `RuntimeError`, and with neither the sphere code runs unchanged. The result then also carries
`orientations`.

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
| `ValueError: walls=...` | The boundary specification is malformed; the message lists the valid forms, and `rcpgenerator.example()` has a snippet for each boundary type. |
| `UserWarning: periodic box is small for the largest particle` | A periodic box length is under about 2 largest diameters (2.1 in 2D/3D; in 4D and above more, up to about 3 for equal sizes, as the message states), so particles could interact through more than one periodic copy, which the forces do not account for. Common with few particles in 4D and above, or with broad size distributions. Use more particles, or hard walls (`walls` entry `1`) on that axis. |
| `custom distribution needs N entries` | Provide exactly one custom diameter per particle. |
| Initial placement exceeds its attempt limit | The requested initial `phi`, size ratio, or confined geometry is too difficult for non-overlapping random placement. Lower initial `phi` or revise the geometry/distribution. |
| `exceeded its neighbor capacity` | Increase `neighbor_max` or reduce the particle-size ratio; see [Neighbor capacity and memory](#neighbor-capacity-and-memory). |
| No trajectory recorded | Call `pack(capture_trajectory=True, ...)` before animation. |
| Unsupported rendering dimension | Built-in rendering is limited to 2D and 3D, although the packing core supports higher dimensions. |
| `ValueError: aspect_ratios = [...] breaks the order rule` | The ratios are not in decreasing order or one exceeds 1. `D` is the longest axis: divide `D` by the largest ratio and list the ratios longest first, for example `[0.5]` and not `[2.0]`. The input is never reordered for you. |
| `aspect_ratios needs Ndim-1 = ... entries` | Give one ratio for each axis after the first: one number in 2D, two in 3D. |
| `a per-particle aspect_ratios table needs one row per particle`, `rod_length as a list needs one entry per particle` | A per-particle table or list must have exactly `N` entries. |
| `orientations were given without aspect_ratios (or rod_length and rod_radius)` | Orientations need a shape. Supply the shape as well, or leave both out for spheres. |
| `orientations needs one rotation per particle` | Pass exactly `N` rotation matrices, each `Ndim × Ndim`. |
| `rods need both rod_length and rod_radius` | Give both keywords; `rod_radius=0.5` makes `D` the rod's thickness. |
| `rod_length / rod_radius were given together with aspect_ratios` | A packing is of rods or of ellipsoids, not both. |
| `rod_length must be positive: an uncapped rod of length 0 has no volume`, `rod_radius must be positive` | A flat-ended rod needs a positive length (a capped rod may have length 0, which is a sphere); every rod needs a positive radius. |
| `NotImplementedError: uncapped (flat-ended) rods are implemented in 2D and 3D only` | Use capped rods above 3D. |
| `NotImplementedError: ... is not supported with aspect_ratios (ellipsoids) yet` (or `with rods (rod_length / rod_radius) yet`) | The named feature (a research input such as the distribution morph or fixed particles) has not been extended to these shapes; see [What works with them and what refuses them](#what-works-with-them-and-what-refuses-them). |
| `NotImplementedError: quench() of rods uses the rounded-rod contact` | Rods can be quenched in 2D and 3D only. |
| `NotImplementedError: animate_2d with a shape per particle AND diameter_permute is not supported` | Animate a run without permutation, or with one shape for all particles. |
| `NotImplementedError: viewer bundle: a per-particle aspect_ratios table is not supported yet` | The viewer draws one ellipsoid shape for all particles; `Packing.render()` draws tables. |
| `exceeded its neighbor capacity` with long ellipsoids or rods | The neighbour search works on bounding spheres, so very elongated particles have many neighbours. See [Cost and limits](#cost-and-limits). |
| `RuntimeWarning: quench() did not converge` | Read `packing.quench_message`; the [quench status table](#quench-to-force-balance) lists the causes and remedies. With ellipsoids or rods, first quench the packed state again with a lower `dt_max`. The packing keeps the state the quench stopped in. |

Use the exact lowercase distribution names documented above. An unsupported
distribution name is not a recoverable Python validation path in the current
initializer and should be avoided.

## Interactive packing app

🧪 [Open the interactive packing app](https://kd-physics.github.io/RCPGenerator/webapp/index.html)

The hosted app is a self-contained 2D packing experiment. Choose the particle
count, initial packing fraction, and diameter distribution, generate a periodic
initial state, and run or stop its browser-based ADAM relaxation. The μ and
learning-rate sliders can be adjusted while it runs to explore how those
controls affect overlaps, energy, packing fraction, and the evolving
configuration. Positions and diameters can be downloaded as CSV.

This educational app is separate from the maintained `rcpgenerator` Python
package and compiled engine.

## Large-packing visualizer

🔍 [Open the large-packing visualizer](https://kd-physics.github.io/RCPGenerator/python_code/python/webapp/app/index.html)

The separate visualizer loads packings produced by `rcpgenerator`. It provides
map-like pan and zoom from the full system down to individual particles while
keeping the rendering crisp. It displays 2D packings directly and 3D packings
through an adjustable cross-sectional slice, with periodic wrapping, hard and
curved boundaries, color and analysis layers, field-of-view statistics, and
PNG export.

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

For ellipsoids and rods, pass the shape and the orientations as well; the call is shown under
[Rendering, the viewer and animation](#rendering-the-viewer-and-animation).

Open the visualizer and select the generated bundle folder or ZIP file. Solid
boundaries represent hard walls and dashed boundaries represent periodic
boundaries. Repository tools can also create a self-contained HTML file with
the packing data embedded for offline inspection.

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
as perfectly force-balanced or mechanically stable; `Packing.quench()` brings a
state to force balance where that is reachable and reports when it is not. Built-in visual rendering
is limited to 2D and 3D, while the numerical core supports higher dimensions. Ellipsoids and rods
cost more than spheres and stop working at large elongations; the measured numbers are under
[Cost and limits](#cost-and-limits).

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
