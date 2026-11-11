"""Built-in catalogue of short, copy-paste snippets.

    rcpgenerator.example()            # list every example, grouped
    rcpgenerator.example("cylinder")  # print one snippet

Every snippet is complete (imports included) and runs as printed at a small particle
count; change the numbers and go. The package's test suite executes each one.
"""

from __future__ import annotations

import difflib

# (group title, [(name, descriptor, snippet), ...]) in display order
_CATALOGUE: list[tuple[str, list[tuple[str, str, str]]]] = [
    ("2D Examples", [
        ("mono_2d", "monodisperse disks, periodic", '''\
import rcpgenerator
packing = rcpgenerator.Packing(N=500, Ndim=2, phi=0.1, box=[1.0, 1.0], walls=[0, 0],
                               dist={"type": "mono", "d": 1.0}, seed=1)
packing.pack()
print("phi =", packing.phi_final)
packing.show_packing()
'''),
        ("bidisperse_2d", "two disk sizes", '''\
import rcpgenerator
dist = {"type": "bidisperse", "d1": 1.0, "d2": 1.4, "p": 0.5}   # p = fraction with diameter d1
packing = rcpgenerator.Packing(N=500, Ndim=2, phi=0.1, box=[1.0, 1.0], walls=[0, 0],
                               dist=dist, seed=1)
packing.pack()
print("phi =", packing.phi_final)
packing.show_packing()
'''),
        ("polydisperse_2d", "lognormal disk sizes", '''\
import rcpgenerator
dist = {"type": "lognormal", "mu": 0.0, "sigma": 0.3}   # parameters of the underlying normal
packing = rcpgenerator.Packing(N=500, Ndim=2, phi=0.1, box=[1.0, 1.0], walls=[0, 0],
                               dist=dist, seed=1)
packing.pack()
print("phi =", packing.phi_final)
packing.show_packing()
'''),
        ("target_phi_2d", "stop at chosen fraction", '''\
import rcpgenerator
packing = rcpgenerator.Packing(N=500, Ndim=2, phi=0.1, box=[1.0, 1.0], walls=[0, 0],
                               dist={"type": "bidisperse", "d1": 1.0, "d2": 1.4, "p": 0.5}, seed=1)
packing.relax(n_steps=5000, target_phi=0.80)    # grow to phi = 0.80, then relax at fixed size
packing.update_phi(0.82)                        # rescale the diameters to a new fraction
packing.relax(n_steps=1000, fix_diameter=True)  # relax positions only
print("phi =", packing.phi_final)
'''),
        ("animate_2d", "movie of the growth", '''\
import rcpgenerator
packing = rcpgenerator.Packing(N=300, Ndim=2, phi=0.1, box=[1.0, 1.0], walls=[0, 0],
                               dist={"type": "bidisperse", "d1": 1.0, "d2": 1.4, "p": 0.5}, seed=1)
packing.pack(capture_trajectory=True, trajectory_interval=500)   # record every 500 steps
packing.animate_2d(path="packing.gif", show=False)               # GIF is the saved format
'''),
    ]),
    ("3D Examples", [
        ("mono_3d", "monodisperse spheres, periodic", '''\
import rcpgenerator
packing = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
                               dist={"type": "mono", "d": 1.0}, seed=1)
packing.pack()
print("phi =", packing.phi_final)   # packing.show_packing() draws it (3D rendering is slow for large N)
'''),
        ("bidisperse_3d", "two sphere sizes", '''\
import rcpgenerator
dist = {"type": "bidisperse", "d1": 1.0, "d2": 1.4, "p": 0.5}   # p = fraction with diameter d1
packing = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
                               dist=dist, seed=1)
packing.pack()
print("phi =", packing.phi_final)
'''),
        ("powerlaw_3d", "power-law sphere sizes", '''\
import rcpgenerator
dist = {"type": "powerlaw", "d_min": 1.0, "d_max": 10.0, "exponent": -2.5}   # number density ~ D**exponent
packing = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
                               dist=dist, seed=1)
packing.pack()
print("phi =", packing.phi_final, " size ratio =", max(packing.diameters) / min(packing.diameters))
'''),
        ("arrays_3d", "positions and diameters as arrays", '''\
import numpy as np, rcpgenerator
packing = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
                               dist={"type": "lognormal", "mu": 0.0, "sigma": 0.2}, seed=1)
packing.pack()
x = np.asarray(packing.positions)     # shape (N, 3), each coordinate in [0, box)
d = np.asarray(packing.diameters)     # shape (N,), same particle order as x
np.savetxt("packing.csv", np.column_stack([x, d]), delimiter=",", header="x,y,z,diameter")
'''),
    ]),
    ("Higher Dimensions", [
        ("permute_4d", "4D packing with permutation", '''\
import rcpgenerator
packing = rcpgenerator.Packing(N=2000, Ndim=4, phi=0.1, box=[1.0] * 4, walls=[0] * 4,
                               dist={"type": "powerlaw", "d_min": 1.0, "d_max": 3.0, "exponent": -3.0},
                               seed=1, diameter_permute=True)
packing.pack()                # any Ndim >= 2 packs; rendering is 2D and 3D only
print(packing.phi_final, packing.permute_exchanges)
'''),
        ("quench_5d_walls", "5D quench, hard walls", '''\
import rcpgenerator
# few particles in 5D: hard walls, since a periodic box would be too narrow for the largest particles
packing = rcpgenerator.Packing(N=1000, Ndim=5, phi=0.05, box=[1.0] * 5, walls=[1] * 5,
                               dist={"type": "powerlaw", "d_min": 1.0, "d_max": 5.0, "exponent": -3.0}, seed=1)
packing.pack()
packing.quench()              # above 3D the default compression is kappa_mult=1.01
print(packing.quench_status, packing.quench_iters)
'''),
    ]),
    ("Boundary Conditions Examples", [
        ("walls_one_axis", "hard walls along z", '''\
import rcpgenerator
# walls: one entry per axis; 0 = periodic, 1 = hard walls at both ends of that axis
packing = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 1],
                               dist={"type": "mono", "d": 1.0}, seed=1)
packing.pack()
print("phi =", packing.phi_final)
'''),
        ("walls_all_axes", "closed hard-wall box", '''\
import rcpgenerator
packing = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[1, 1, 1],
                               dist={"type": "mono", "d": 1.0}, seed=1)
packing.pack()
print("phi =", packing.phi_final)
'''),
        ("channel_2d", "2D channel, walls in y", '''\
import rcpgenerator
packing = rcpgenerator.Packing(N=500, Ndim=2, phi=0.1, box=[2.0, 1.0], walls=[0, 1],   # periodic in x
                               dist={"type": "bidisperse", "d1": 1.0, "d2": 1.4, "p": 0.5}, seed=1)
packing.pack()
packing.show_packing()
'''),
        ("disk", "2D circular container", '''\
import rcpgenerator
# a curved container is a NEGATIVE FIRST entry: walls[0] = -2 makes axes 1-2 a disk of diameter box[0]
packing = rcpgenerator.Packing(N=500, Ndim=2, phi=0.1, box=[1.0, 1.0], walls=[-2, 0],
                               dist={"type": "mono", "d": 1.0}, seed=1)
packing.pack()
print("phi (relative to the disk area) =", packing.phi_final)
packing.show_packing()
'''),
        ("cylinder", "periodic-axis cylinder", '''\
import rcpgenerator
# walls[0] = -2: x-y form a circle of diameter box[0]; walls[2] = 0: the cylinder axis z is periodic
packing = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 2.0], walls=[-2, 0, 0],
                               dist={"type": "mono", "d": 1.0}, seed=1)
packing.pack()
print("phi (relative to the cylinder volume) =", packing.phi_final)
'''),
        ("capped_cylinder", "cylinder with end walls", '''\
import rcpgenerator
# walls[0] = -2: circular cross-section in x-y; walls[2] = 1: hard end caps along z
packing = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 2.0], walls=[-2, 0, 1],
                               dist={"type": "mono", "d": 1.0}, seed=1)
packing.pack()
print("phi (relative to the cylinder volume) =", packing.phi_final)
'''),
        ("tall_box", "height in particle diameters", '''\
import rcpgenerator
# fix_height=True: the LAST box entry is a multiple of the largest particle diameter and stays so while
# the particles grow; here a slab 4 diameters thick between hard walls, periodic in x and y
packing = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 4.0], walls=[0, 0, 1],
                               fix_height=True, dist={"type": "mono", "d": 1.0}, seed=1)
packing.pack()
print("height / largest diameter =", packing.box[-1] / max(packing.diameters), " phi =", packing.phi_final)
'''),
        ("sphere", "spherical 3D container", '''\
import rcpgenerator
# walls[0] = -3: all three axes form a sphere of diameter box[0] (the other entries are ignored)
packing = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[-3, 0, 0],
                               dist={"type": "mono", "d": 1.0}, seed=1)
packing.pack()
print("phi (relative to the sphere volume) =", packing.phi_final)
'''),
    ]),
    ("Arbitrary Distribution", [
        ("custom_list", "your own diameter list", '''\
import rcpgenerator
diameters = [1.0] * 300 + [1.5] * 150 + [2.5] * 50   # any N values; only their RATIOS matter
packing = rcpgenerator.Packing(N=len(diameters), Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
                               dist={"type": "custom", "custom": diameters}, seed=1)
packing.pack()        # the absolute size is set by the packing; the multiset of ratios is kept
print("phi =", packing.phi_final)
'''),
        ("from_inverse_cdf", "any shape via quantiles", '''\
import numpy as np, rcpgenerator
N, S, p = 1000, 10.0, -3.0                    # count, size ratio D_max/D_min, number density ~ D**p
u = (np.arange(N) + 0.5) / N                  # evenly spaced quantiles: no sampling noise
D = (1.0 + u * (S ** (p + 1) - 1.0)) ** (1.0 / (p + 1))   # inverse CDF; swap in any distribution's
packing = rcpgenerator.Packing(N=N, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
                               dist={"type": "custom", "custom": D.tolist()}, seed=1)
packing.pack()
print("phi =", packing.phi_final)
'''),
        ("from_file", "measured sizes from file", '''\
import numpy as np, rcpgenerator
D = np.loadtxt("diameters.txt")               # one diameter per line, any units
packing = rcpgenerator.Packing(N=len(D), Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
                               dist={"type": "custom", "custom": D.tolist()}, seed=1)
packing.pack()
print("phi =", packing.phi_final)
'''),
        ("mixture", "mix of two distributions", '''\
import numpy as np, rcpgenerator
rng = np.random.default_rng(0)
D = np.concatenate([rng.lognormal(0.0, 0.1, 700),            # 70 % small population
                    rng.lognormal(np.log(3.0), 0.1, 300)])   # 30 % large population
packing = rcpgenerator.Packing(N=len(D), Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
                               dist={"type": "custom", "custom": D.tolist()}, seed=1)
packing.pack()
print("phi =", packing.phi_final)
'''),
        ("builtin_types", "every built-in distribution", '''\
import rcpgenerator
dists = [{"type": "mono", "d": 1.0}, {"type": "bidisperse", "d1": 1.0, "d2": 1.4, "p": 0.5},
         {"type": "gaussian", "mu": 1.0, "sigma": 0.2}, {"type": "lognormal", "mu": 0.0, "sigma": 0.3},
         {"type": "flat", "d_min": 0.5, "d_max": 1.5}, {"type": "exponential", "d_min": 0.5, "d_max": 2.0},
         {"type": "powerlaw", "d_min": 1.0, "d_max": 10.0, "exponent": -2.5}, {"type": "weibull", "shape": 2.0, "scale": 1.0},
         {"type": "bigaussian", "mu1": 1.0, "sigma1": 0.1, "mu2": 2.0, "sigma2": 0.2, "p": 0.5}]
for dist in dists:   # construction draws the diameters; pack() would then pack them
    p = rcpgenerator.Packing(N=300, Ndim=2, phi=0.1, box=[1.0, 1.0], walls=[0, 0], dist=dist, seed=1)
    print(f"{dist['type']:<12} size ratio {max(p.diameters) / min(p.diameters):.2f}")
'''),
    ]),
    ("Implementing Diameter Permute", [
        ("permute", "switch on, default recipe", '''\
import rcpgenerator
dist = {"type": "powerlaw", "d_min": 1.0, "d_max": 10.0, "exponent": -3.0}
packing = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
                               dist=dist, seed=1, diameter_permute=True)   # off by default
packing.pack()
phi = packing.phi_history[-1] * (1 - packing.max_overlap_history[-1]) ** 3   # overlaps removed
print(f"phi = {phi:.4f} after {packing.permute_exchanges} diameter exchanges")
'''),
        ("permute_dial", "dial from little to lot", '''\
import rcpgenerator
dist = {"type": "powerlaw", "d_min": 1.0, "d_max": 10.0, "exponent": -3.0}
for budget in (1e-3, 1e-2, 0.25):             # exchanges per pass as a fraction of N; 0.25 = default
    packing = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
                                   dist=dist, seed=1, diameter_permute=True, permute_budget=budget)
    packing.pack()
    phi = packing.phi_history[-1] * (1 - packing.max_overlap_history[-1]) ** 3
    print(f"budget {budget:<6} phi = {phi:.4f}")
'''),
        ("permute_vs_off", "same particles, off and on", '''\
import rcpgenerator
dist = {"type": "powerlaw", "d_min": 1.0, "d_max": 10.0, "exponent": -3.0}
for on in (False, True):
    packing = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
                                   dist=dist, seed=1, diameter_permute=on)
    packing.pack()
    phi = packing.phi_history[-1] * (1 - packing.max_overlap_history[-1]) ** 3
    print(f"permute {'on ' if on else 'off'}  phi = {phi:.4f}")
'''),
        ("permute_diagnostics", "who exchanged, how often", '''\
import numpy as np, matplotlib.pyplot as plt, rcpgenerator
packing = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
                               dist={"type": "powerlaw", "d_min": 1.0, "d_max": 10.0, "exponent": -3.0},
                               seed=1, diameter_permute=True)
packing.pack()
print(packing.permute_passes, "passes,", packing.permute_exchanges, "exchanges")
plt.scatter(packing.diameters, packing.permute_counts, s=4)   # per particle, same order as positions
plt.xscale("log"); plt.xlabel("diameter"); plt.ylabel("exchanges taken part in"); plt.show()
'''),
    ]),
    ("Quenching to Force Balance", [
        ("quench", "balance a finished packing", '''\
import rcpgenerator
packing = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
                               dist={"type": "lognormal", "mu": 0.0, "sigma": 0.3}, seed=1)
packing.pack()
packing.quench()          # diameters x 1.003 once, then minimise until every net force < 1e-10 of its load
print(packing.quench_status, packing.quench_residual, packing.quench_iters)
'''),
        ("quench_staged", "loose pass, then tighten", '''\
import rcpgenerator
packing = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
                               dist={"type": "lognormal", "mu": 0.0, "sigma": 0.3}, seed=1)
packing.pack()
packing.quench(tol=1e-4)                  # quick first pass
packing.quench(tol=1e-10, resume=True)    # continue from there; resume does not rescale again
print(packing.quench_status, packing.quench_residual)
'''),
        ("quench_give_up", "handle a quench failure", '''\
import warnings, rcpgenerator
packing = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
                               dist={"type": "powerlaw", "d_min": 1.0, "d_max": 5.0, "exponent": -2.5}, seed=3)
packing.pack()
start = packing.copy()                    # keep the packed state: a failed quench keeps where it stopped
with warnings.catch_warnings():
    warnings.simplefilter("ignore", RuntimeWarning)
    packing.quench(max_iter=2000)         # deliberately too few iterations
if packing.quench_status != "converged":
    print(packing.quench_status, "-", packing.quench_message)   # the reason, and what to try
    packing.quench(resume=True)           # here: continue where it stopped
print(packing.quench_status)
'''),
        ("quench_hold", "hold particles in place", '''\
import numpy as np, rcpgenerator
packing = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
                               dist={"type": "lognormal", "mu": 0.0, "sigma": 0.3}, seed=1)
packing.pack()
x = np.asarray(packing.positions)
bottom = np.flatnonzero(x[:, 2] < 0.1 * packing.box[2])     # indices into positions, as they are now
packing.quench(held=bottom.tolist())       # held particles do not move; particle i keeps index i
print(packing.quench_status, "net force on the held set:", packing.quench_held_force)
'''),
        ("quench_shear", "one shear step, stress", '''\
import numpy as np, rcpgenerator
p = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
                         dist={"type": "lognormal", "mu": 0.0, "sigma": 0.3}, seed=1)
p.pack(); p.quench()
x, d = np.asarray(p.positions), np.asarray(p.diameters)
moved, fixed = np.flatnonzero(x[:, 0] < 0.1), np.flatnonzero((x[:, 0] >= 0.5) & (x[:, 0] < 0.6))  # half a box apart
x[moved, 1] = (x[moved, 1] + 1e-5) % 1.0; p.positions = x.tolist()      # displace one slab along y
p.quench(resume=True, held=np.concatenate([moved, fixed]).tolist())
print(p.quench_status, "- slab forces: see example('contact_forces') to resolve each slab")
'''),
    ]),
    ("Other", [
        ("contact_forces", "per-particle contact forces", '''\
import numpy as np, rcpgenerator
from scipy.spatial import cKDTree
p = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
                         dist={"type": "lognormal", "mu": 0.0, "sigma": 0.3}, seed=1)
p.pack(); p.quench()
x, d, L = np.asarray(p.positions), np.asarray(p.diameters), np.asarray(p.box)
i, j = cKDTree(x, boxsize=L).query_pairs(d.max(), output_type="ndarray").T   # periodic box only
dx = x[j] - x[i]; dx -= np.round(dx / L) * L; r = np.linalg.norm(dx, axis=1)
ov = (d[i] + d[j]) / 2 - r; c = ov > 0                    # contacts; force = overlap along the line of centres
F = np.zeros_like(x); f = (ov[c] / r[c])[:, None] * dx[c]; np.add.at(F, j[c], f); np.add.at(F, i[c], -f)
print(c.sum(), "contacts; largest net force", np.abs(F).max(), "vs mean contact force", ov[c].mean())
'''),
        ("phi_no_overlap", "phi with overlaps removed", '''\
import rcpgenerator
from rcptools import overlap_report, phi_corrected
packing = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
                               dist={"type": "mono", "d": 1.0}, seed=1)
packing.pack()
report = overlap_report(packing)          # pairwise overlaps of the final state
print("raw phi", packing.phi_final, " after removing overlaps", phi_corrected(packing.phi_final, report["max_overlap"], 3))
'''),
        ("plot_history", "phi and energy versus step", '''\
import matplotlib.pyplot as plt, rcpgenerator
packing = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
                               dist={"type": "mono", "d": 1.0}, seed=1)
packing.pack()
fig, (a, b) = plt.subplots(2, 1, sharex=True)
a.plot(packing.phi_history); a.set_ylabel("phi")
b.semilogy(packing.energy_history); b.set_ylabel("energy"); b.set_xlabel("step")
plt.show()
'''),
        ("save_image", "PNG or SVG figure", '''\
import rcpgenerator
packing = rcpgenerator.Packing(N=500, Ndim=2, phi=0.1, box=[1.0, 1.0], walls=[0, 0],
                               dist={"type": "lognormal", "mu": 0.0, "sigma": 0.3}, seed=1)
packing.pack()
packing.render(path="packing.png", palette_choice=3, figsize=(6, 6), dpi=150)   # palettes 1-12
packing.savefig("packing.svg")
'''),
        ("save_and_load", "save, reload as JSON", '''\
import json, rcpgenerator
packing = rcpgenerator.Packing(N=500, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
                               dist={"type": "mono", "d": 1.0}, seed=1)
packing.pack()
with open("packing.json", "w") as f:
    json.dump(packing.to_dict(), f)
with open("packing.json") as f:
    restored = rcpgenerator.Packing.from_dict(json.load(f))   # no re-initialisation; continue or quench it
print(restored.phi_final)
'''),
        ("threads", "use more CPU threads", '''\
import rcpgenerator
rcpgenerator.set_num_threads(4)     # default is one thread; set before the run
print("threads:", rcpgenerator.get_num_threads())
# one thread is bit-for-bit reproducible; several are faster at large N but may differ in the last digits
'''),
        ("viewer_bundle", "open in web viewer", '''\
import rcpgenerator
from rcptools.bridge import write_bundle
packing = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 1],
                               dist={"type": "lognormal", "mu": 0.0, "sigma": 0.3}, seed=1)
packing.pack()
write_bundle("my_packing", packing.positions, packing.diameters, packing.box, walls=packing.walls)
# open https://kd-physics.github.io/RCPGenerator/webapp/index.html and select the my_packing folder
'''),
        ("voronoi", "local phi (needs pyvoro)", '''\
import rcpgenerator                  # first: pip install pyvoro-mmalahe
from rcptools.voronoi import voronoi_phi_local
packing = rcpgenerator.Packing(N=1000, Ndim=3, phi=0.05, box=[1.0, 1.0, 1.0], walls=[0, 0, 0],
                               dist={"type": "lognormal", "mu": 0.0, "sigma": 0.3}, seed=1)
packing.pack()
v = voronoi_phi_local(packing)       # periodic box; per-particle Voronoi cell volumes
print("mean local phi", v["local_phi"].mean(), " cells tile the box:", v["sum_match"])
'''),
    ]),
]

_SNIPPETS = {name: code for _, entries in _CATALOGUE for name, _, code in entries}


def example(name: str | None = None) -> None:
    """Print the catalogue of examples, or the snippet called `name`.

    ``rcpgenerator.example()`` lists every example by group; ``rcpgenerator.example("disk")``
    prints a short, complete script to copy, paste and modify.
    """
    if name is None:
        lines = []
        for title, entries in _CATALOGUE:
            lines.append(title)
            lines.extend(f"  - '{n}': {desc}" for n, desc, _ in entries)
        lines.append("")
        lines.append("Print one with rcpgenerator.example('name').")
        print("\n".join(lines))
        return
    key = str(name).strip().strip("'\"")
    if key not in _SNIPPETS:
        close = difflib.get_close_matches(key, list(_SNIPPETS), n=3, cutoff=0.5)
        hint = f" Did you mean {', '.join(repr(c) for c in close)}?" if close else ""
        raise ValueError(f"no example named {name!r}.{hint} Call rcpgenerator.example() for the list.")
    print(_SNIPPETS[key], end="")


def example_names() -> list[str]:
    """All example names, in catalogue order."""
    return list(_SNIPPETS)


def example_source(name: str) -> str:
    """The snippet called `name`, as a string (for tooling and tests)."""
    return _SNIPPETS[name]
