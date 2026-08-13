# rcpgenerator

`rcpgenerator` is a Python package for generating dense packings of polydisperse
particles with a compiled C++/OpenMP engine. It supports 2D disks, 3D spheres, and
higher-dimensional hyperspheres in periodic boxes, hard-wall boxes, and curved
containers. The high-level `rcpgenerator.Packing` API manages initialization,
relaxation, results, and rendering. The distribution also includes `rcptools` for
packing analysis, parameter searches, persistence, and viewer-bundle export.

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

GitHub Actions has verified wheel building, dependency repair, clean installation, and
the installed-package smoke test for every listed Python version and platform. Wheels
have also been installed and run manually on Windows x86-64, macOS Apple Silicon, and
Linux x86-64. The macOS Intel wheels are verified through GitHub Actions; no manual test
on a physical Intel Mac is currently claimed.

Installing a compatible wheel does not require Git, CMake, Visual Studio, Xcode, a
local C++ compiler, or a separately configured OpenMP installation. The compiled
extension and its required OpenMP runtime are handled by the platform wheel.

## Quickstart

```python
import rcpgenerator

rcpgenerator.set_num_threads(1)  # optional; the default uses all available threads

packing = rcpgenerator.Packing(
    phi=0.08,
    N=64,
    Ndim=2,
    box=[1.0, 1.0],
    walls=[0, 0],  # 0 = periodic, 1 = hard wall
    dist={"type": "mono", "d": 1.0},
    neighbor_max=0,  # automatic capacity
    seed=20240719,
)
result = packing.pack()

print("phi =", packing.phi_final)
print("steps =", packing.steps)
print("positions =", len(result["positions"]))

# Display the packing in an interactive Python session:
# packing.show_packing()
```

The [`examples` directory](https://github.com/KD-physics/RCPGenerator/tree/main/python_code/python/examples)
contains complete box, bidisperse, polydisperse, circular, cylindrical, spherical,
and target-packing-fraction cases.

## Containers and boundaries

For a rectangular `box`, each entry in `walls` controls the corresponding dimension:

- `0` selects a periodic boundary.
- `1` selects a hard wall.
- A negative first entry, `-t`, selects one hyperspherical hard boundary over the
  first `t` dimensions. Its diameter is the first component of `box`: `-2` gives a
  disk in 2D or a cylindrical cross-section in 3D, while `-3` gives a sphere in 3D.
- `fix_height=True` makes the final box dimension a fixed multiple of the first
  particle diameter.

See the [Getting Started notebook](https://github.com/KD-physics/RCPGenerator/blob/main/getting_started.ipynb)
for configuration details and worked examples.

## Interactive packing viewer

`rcptools.bridge.write_bundle` creates a viewer bundle containing `manifest.json`,
`pos.f32`, and `dia.f32`:

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
and select the generated bundle folder. Solid boundaries represent hard walls and
dashed boundaries represent periodic boundaries; the viewer can inspect 2D packings
directly and 3D packings by cross-section.

The viewer application itself is not included in the Python wheel. Its development
files and `launch.py` helper are available only in the repository's
[`python_code/python/webapp` directory](https://github.com/KD-physics/RCPGenerator/tree/main/python_code/python/webapp).

## `rcptools`

The installed distribution also provides the `rcptools` package. It contains packing
search and scheduling tools, persistence helpers, bundle export, and analyses including
chord, pore, and void-fill workflows. Some analysis workflows have additional
data- or method-specific requirements; consult their module documentation and the
[repository examples](https://github.com/KD-physics/RCPGenerator/tree/main/python_code/python/examples)
before using them in an automated study.

## Interpretation and limitations

Packing is an iterative numerical procedure. Convergence behavior and runtime depend on
the particle distribution, dimensionality, boundaries, initial state, and packing
settings. Inspect the returned diagnostics—such as `steps`, `phi`, `max_min_dist`,
`force_magnitude`, and the recorded histories—and apply validation appropriate to the
intended scientific use. For reproducibility studies, record the package version, seed,
thread count, inputs, and relevant platform details.

## Building from source

Source builds require Git, a C++17 compiler, CMake-compatible build tooling, and a usable
OpenMP development/runtime installation. Platform-specific compiler setup is the
builder's responsibility.

```bash
git clone https://github.com/KD-physics/RCPGenerator.git
cd RCPGenerator/python_code/python
pip install -v .
```

Normal source builds use portable CPU instructions. A local build may explicitly tune
the compiled core for the current machine:

```bash
CMAKE_ARGS="-DRCP_NATIVE_OPTIMIZATION=ON" pip install -v .
```

In PowerShell, set the option before installing:

```powershell
$env:CMAKE_ARGS = "-DRCP_NATIVE_OPTIMIZATION=ON"
pip install -v .
```

GCC native builds enable `-march=native` and the configured 512-bit vector-width
preference; Clang native builds enable `-march=native`. Other compilers emit a CMake
warning and retain their default architecture. A native-optimized build is intended only
for the machine on which it is compiled. Do not redistribute it as a portable wheel.

## Method and citation

The package uses an iterative expansion–relaxation method based on the approach described
by Desmond and Weeks, with a compiled neighbor-search and relaxation implementation.

- Desmond, K. “Random close packing at extreme size ratios with an Adam-based
  inflation protocol.” arXiv:2608.12235 (2026).
  [arXiv:2608.12235](https://arxiv.org/abs/2608.12235)
- Desmond, K. W. and Weeks, E. R. “Random close packing of disks and spheres in
  confined geometries.” *Physical Review E* **80**, 051305 (2009).
  [arXiv:0903.0864](https://arxiv.org/abs/0903.0864)
- Desmond, K. W. and Weeks, E. R. “Influence of particle size distribution on random
  close packing of spheres.” *Physical Review E* **90**, 022204 (2014).
  [arXiv:1303.4627](https://arxiv.org/abs/1303.4627)

Use the RCPGenerator arXiv manuscript as the preferred scientific citation and record
the `rcpgenerator` version used. The two earlier papers provide method background.

`rcpgenerator` is distributed under the
[MIT License](https://github.com/KD-physics/RCPGenerator/blob/main/python_code/python/LICENSE).
