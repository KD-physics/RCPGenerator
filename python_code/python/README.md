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
[interactive packing app](https://kd-physics.github.io/RCPGenerator/webapp/index.html) ·
[large-packing visualizer](https://kd-physics.github.io/RCPGenerator/python_code/python/webapp/app/index.html) ·
[issue tracker](https://github.com/KD-physics/RCPGenerator/issues)

## Installation

```bash
pip install rcpgenerator
```

Version 1.1.0 wheels are built and tested for CPython 3.10–3.14 on:

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

## Built-in examples

`rcpgenerator.example()` lists short, complete snippets by topic: 2D, 3D and higher-dimensional packings, every
boundary and container type, arbitrary size distributions, diameter permutation, the quench,
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
| `render_packing(packing, ...)` | Render a 2D or 3D packing. |
| `animate_packing_2d(packing, ...)` | Animate a previously captured 2D trajectory. |
| `example(name=None)` | List the built-in snippets, or print one to copy. See [Built-in examples](#built-in-examples). |

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
| `diameter_permute` | `bool`, `False` | Turn on diameter permutation (which particle carries which diameter is re-assigned during compression). `False` reproduces the original algorithm bit for bit. See [Diameter permutation](#diameter-permutation). |
| `permute_budget` | `float` in (0, 0.5], `0.25` | Exchanges per permutation pass as a fraction of N. Smaller = less uplift. |
| `permute_interval` | `int` ≥ 1, `5` | Steps between permutation passes. Longer = fewer passes and less uplift. |
| `permute_start` | `int`, `4000` | Step at which permutation begins. Later = less uplift. |

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
| `"stalled"` | The best residual has not halved within `stall_window` iterations. | `quench(resume=True)` to keep trying, or quench the packed state with a larger `kappa_mult`. Below about 1e-8 the stop is the floating-point floor (roughly 6e-14 divided by the mean relative overlap), and the state is balanced as far as it can be. |
| `"unjammed"` | No load-bearing particles remain: the diameters are too small for the structure to carry load. | A larger `kappa_mult`. |
| `"step_collapse"` | FIRE's time step collapsed: it found no step it could follow downhill (typical of loose, under-compressed states). | A larger `kappa_mult`. |
| `"non_finite"` | A non-finite force or energy appeared. | Report it; do not use the state. |

Every status other than `"converged"` raises a `RuntimeWarning`, and `quench_message`
explains the outcome and the likely cause. Whatever the outcome, the packing keeps the
state the quench stopped in, and particle `i` keeps index `i` (positions, diameters,
`particle_origin` and `permute_counts` stay aligned).

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
| `ValueError: walls=...` | The boundary specification is malformed; the message lists the valid forms, and `rcpgenerator.example()` has a snippet for each boundary type. |
| `UserWarning: periodic box is small for the largest particle` | A periodic box length is under about 2 largest diameters (2.1 in 2D/3D; in 4D and above more, up to about 3 for equal sizes, as the message states), so particles could interact through more than one periodic copy, which the forces do not account for. Common with few particles in 4D and above, or with broad size distributions. Use more particles, or hard walls (`walls` entry `1`) on that axis. |
| `custom distribution needs N entries` | Provide exactly one custom diameter per particle. |
| Initial placement exceeds its attempt limit | The requested initial `phi`, size ratio, or confined geometry is too difficult for non-overlapping random placement. Lower initial `phi` or revise the geometry/distribution. |
| `exceeded its neighbor capacity` | Increase `neighbor_max` or reduce the particle-size ratio; see [Neighbor capacity and memory](#neighbor-capacity-and-memory). |
| No trajectory recorded | Call `pack(capture_trajectory=True, ...)` before animation. |
| Unsupported rendering dimension | Built-in rendering is limited to 2D and 3D, although the packing core supports higher dimensions. |
| `RuntimeWarning: quench() did not converge` | Read `packing.quench_message`; the [quench status table](#quench-to-force-balance) lists the causes and remedies. The packing keeps the state the quench stopped in. |

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
