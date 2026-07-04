# rcpgenerator

Fast random close packing (RCP) generator for polydisperse spheres in N dimensions — a
parallel C++/OpenMP core with SIMD-accelerated inner loops, wrapped in a stateful Python
API. Supply a size distribution (or explicit diameters) and a container, and get back a
jammed configuration. 2D disks, 3D spheres, hyperspheres; periodic, hard-wall, and curved
(disk / cylinder / sphere) boundaries; size ratios exceeding 100 in 2D.

This is a **performance release** — several times faster than the original with improved
convergence — while the Python API and the packings it produces are unchanged. One install
provides two importable packages: **`rcpgenerator`** (the engine) and **`rcptools`** (search,
per-packing analysis, and bundle / visualisation helpers).

## Install

```bash
git clone https://github.com/KD-physics/RCPGenerator.git
cd RCPGenerator/python_code/python
pip install -v .
```

Requires a C++17 compiler with OpenMP. A guided tour is in `getting_started.ipynb` (repo root).

## Quickstart

```python
import rcpgenerator
rcpgenerator.set_num_threads(4)                 # optional; default = all cores

p = rcpgenerator.Packing(phi=0.11, N=1000, Ndim=3, box=[1., 1., 1.],
                         walls=[0, 0, 0],       # 0 = periodic, 1 = hard wall
                         dist={"type": "lognormal", "mu": 0.0, "sigma": 0.3},
                         neighbor_max=0, seed=0)
p.pack()
print("phi =", p.phi_final, " steps =", p.steps)
p.show_packing()                                # inline render
```

## Containers & boundaries (`walls`)

- `0` = periodic, `1` = hard wall, per dimension.
- A leading negative `-t` = one hyperspherical hard boundary over the first `t` dims
  (its diameter is the first `box` component): `-2` → disk (2D) / cylinder (3D),
  `-3` → sphere (3D).
- `fix_height=True` makes the last box dimension a fixed multiple of the first particle
  diameter.

`examples/` covers box (mono / bidisperse / lognormal / power-law), circular, cylindrical,
spherical, and target-φ staging cases; run them with `python examples/run_all_examples.py`.

## Interactive web viewer

Save any packing as a **bundle** (`manifest.json` + raw `pos.f32` / `dia.f32`) and open it in
the viewer to rotate, slice, and colour it:

```python
from rcptools.bridge import write_bundle
write_bundle("bundles/my_packing", p.positions, p.diameters, list(p.box), walls=list(p.walls))
```

The viewer ships in the package (`webapp/`). Build a self-contained HTML for your bundle and
open it in any browser:

```bash
python webapp/launch.py bundles/my_packing
```

(Or serve `webapp/` with a local web server — e.g. `python -m http.server` in `webapp/` — and
load the `bundles/my_packing` folder via the **Choose bundle folder** button.)

In the viewer a **solid** boundary is a hard wall and a **dashed** boundary is periodic; a
sphere's cross-section circle grows and shrinks as you scrub through z. (Dragging a whole
*folder* onto a `file://` page is unreliable across browsers — use `launch.py`, the button, or
multiselect the three files.)

## `rcptools`

A companion toolbox: search over size-distribution families, per-packing analyses
(`rcptools.analysis`: pore network, chords, void fill), and bundle export (`rcptools.bridge`).
It grew out of a density-search effort; explore it with `import rcptools; dir(rcptools)`.

## Method & citation

Iterative expansion–relaxation (Desmond & Weeks 2009) driven by an ADAM-based inflation
schedule, with a kd-tree neighbour manager for high polydispersity.

- Desmond, K. W. & Weeks, E. R. *Random close packing of disks and spheres in confined
  geometries.* Phys. Rev. E **80**, 051305 (2009). [arXiv:0903.0864](https://arxiv.org/abs/0903.0864)
- Desmond, K. W. & Weeks, E. R. *Influence of particle size distribution on random close
  packing of spheres.* Phys. Rev. E **90**, 022204 (2014). [arXiv:1303.4627](https://arxiv.org/abs/1303.4627)

Released under the MIT License.
