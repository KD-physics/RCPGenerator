# Changelog

## 1.2.2

### New
- **`extend`**: `Packing(extend=x)` or `pack(extend=x)` lengthens the first stage of a run (the stage at the starting
  pressure `mu`) by the factor `x` and leaves the rest of the run as it is. A longer first stage gives a denser
  packing: 5,000 spheres with a size ratio of 10 reach `phi_contact` 0.7372 in 17,821 steps without it, 0.7405 in
  32,145 steps with `extend=3`, 0.7426 in 67,843 steps with `extend=8.04`. Not set: the run is unchanged.
- **`adaptive_quench(mu=..., tol=..., max_iter=...)`**: force balance with the diameters free. One common scale
  factor on all diameters follows the pressure term `mu` until the contact forces carry it. Default `mu`: the new
  attribute `mu_final`, the value the packing run ended on. Spheres, ellipsoids and rods (rods in 2D and 3D).
  README, "Adaptive quench".
- `mu_final`: the `mu` a run ended on. `to_dict()` has one more `final_result` key for it.
- **Contact numbers**: `z` (contacts per particle, in the order of `positions`; a wall the particle touches counts
  as one), `z_mean`, `n_rattlers` (fewer than `Ndim + 1` contacts) and `z_mean_no_rattlers`. Filled by `pack()`,
  `quench()` and `adaptive_quench()` and on loading a packed state, for every shape (`None` for spheres above 8
  dimensions). They are counted in the pass that computes `phi_contact`, at no extra cost. Not saved by `to_dict()`.
- **Overlaps**: `mean_overlap` and `max_overlap`, the mean and the largest relative overlap over the pairs that
  overlap (walls not included), filled with the contact numbers. Use `mean_overlap` to compare packing fractions
  at the same degree of overlap.
- **Spheres get `phi_contact` and `shrink_factor` too**: `pack()` fills them for every shape (`None` for spheres
  above 8 dimensions); `phi_final` is unchanged. A sphere `to_dict()` therefore has two more `final_result` keys.
- The exact shrink behind them examines only pairs that can overlap (a cell grid) instead of every pair, with the
  same result bit for bit: 2.5 s for 200,000 equal spheres on one thread.
- `quench(coast=...)` and `adaptive_quench(coast=...)`: whether a particle that loses its last contact keeps its
  velocity (see "Changed: quench").
- `rcpgenerator.example()`: snippets `adaptive_quench` and `extend`.

### Changed: quench
- **`quench()` no longer expands the diameters by default.** `kappa_mult` now defaults to `1.0`; 1.1.0 and 1.2.0
  used 1.003 (spheres in 2D and 3D) or 1.01 (4D and above, ellipsoids, rods). Pass the value to get the former
  behaviour: `quench(kappa_mult=1.003)`. Reason: `quench()` is also the step that restores balance after each change
  to a packing (particles moved or held, a shear step), and repeated calls must not grow the diameters. A packing
  straight from `pack()` quenched without an expansion usually ends as `"unjammed"`; the message says what to pass.
- **Coasting, on by default for spheres.** During a quench a particle that loses its last contact used to stop where
  the contact ended; it now keeps its velocity until the minimiser next resets. Sphere quenches therefore take a
  different path than before and usually fewer iterations (11 of 12 packings measured; two that gave up without it
  converged). `quench(coast=False)` is the former behaviour, bit for bit. For ellipsoids and rods the default is
  off (no gain measured), so their quench results are unchanged.

### Changed: speed presets
The environment variable `RCP_SPEED` now changes only the length of the first stage: `immediate` is `extend=0.5`,
`quick` (the default) is unchanged, `patient` is `extend=4`, `forever` is `extend=10`. Before, the three non-default
presets also changed the number and length of the later cycles; runs with those presets give different results now.
Runs with the default preset are not affected.

### Changed: loading saved states
`Packing.from_dict()` gives a loaded state every attribute a freshly built packing has. A state saved by an earlier
version gets `phi_contact`, `shrink_factor` and the contact numbers computed on loading; `mu_final` is `None` for
it. `copy()` goes through the same path, so a copy now also has the history and diagnostic attributes (empty) that
it lacked before.

## 1.2.0

### New
- **Ellipsoids** (`aspect_ratios`) and **rods** (`rod_length`, `rod_radius`, `rod_capped`; rounded or flat ends)
  pack wherever spheres do: periodic boxes, flat walls, the curved containers (disk, cylinder, sphere, and their
  higher-dimensional forms) and `fix_height`; with diameter permutation, `quench()`, rendering, the viewer and
  `animate_2d()`. One shape for all particles or one per particle. See README, "Ellipsoids and rods".
- `phi_contact`, `shrink_factor` and `shrink_to_contact()`: the packing fraction after an exact uniform shrink to
  zero overlap. `shrink_to_contact()` now also works for spheres in curved containers.
- `pack(capture_trajectory=True, trajectory_spacing="log")`: frames evenly in the logarithm of the step.
  `trajectory_origin`: which particle is in which slot of each frame; `animate_2d()` uses it, so a particle keeps
  its colour through a movie (the packer renumbers particles while it runs).
- `rcpgenerator.example()`: groups of ellipsoid and rod snippets.

### Changed: when a run stops
A run stops when the largest pair overlap is under 5e-4 **and** the normalised force is under its threshold, three
checks in a row (or at the step cap).

- **What was wrong in 1.0.0 and 1.1.0.** The stop test read the force of a step that had not been computed yet,
  which was always zero, so only the overlap was tested. The test now reads the force of the previous step, the
  one the overlap value belongs to.
- **No force left.** The normalised force is the mean net force over the particles that carry one, divided by the
  mean contact force; it is undefined when no particle carries a net force, as when the particles have shrunk out
  of contact. An undefined force counts as the force condition met. (Without this, a run that could not meet the
  force test idled to the step cap.)
- **Who is affected.** With no ellipsoid or rod input, 9 of 11 reference cases are bit-identical to 1.1.0. The
  stopping rule changes one of them: a box too small for its largest particles (500 spheres, size ratio 20),
  whose force at the old stop was 33 times the threshold: 16,098 -> 26,446 steps, phi 0.7641 -> 0.7635. (The other
  is a cylinder case, changed by the step-size correction below.) In the result tables of this project, 7 of 869
  runs moved: 2D ellipses at 12:1 now run to the step cap (39,666 -> 60,000 steps); 3D prolate ellipsoids at 8:1
  take 16,712 -> 26,734 steps; five others take 500 to 12,000 more steps; phi changes between -0.0006 and +0.0007.
  Sphere runs with diameter permutation in 2D (power law 1:5, N = 1000): 11 of 12 unchanged, one takes 17,304 ->
  20,554 steps with phi +0.0002 (seed scatter 0.002). Small systems whose force sits near the threshold can take
  about twice the steps with the same packing fraction: 2D capped rods of length 1 at N = 100 (four of six runs,
  15,700 -> about 29,000 steps, phi within 0.0005); at N = 300 and above none of them moves.
- **Flat-ended rods have their own force thresholds: 2e-2 in 2D and 1e-1 in 3D** (5e-3 for every other shape).
  Their normalised force does not fall below about 6e-3 (2D) and 1e-2 to 5e-2 (3D) with this method, however
  tightly the contact is solved, because of their flat contacts. The thresholds lie above that level, so the force
  test rarely decides when these shapes stop and their packings stop with larger leftover forces than other
  shapes. Use `quench()` when force balance matters.
- "Force" at the stop is the net translational force. Torques are not tested, for any shape; `quench()` balances
  forces and torques.

### Changed: step size in curved containers
In a curved container the step-size feedback compared the packing fraction with a history stored on a different
normalisation, so it always saw the fraction "falling" and cut the step size by 1.5 every 500 steps for the whole
run. It now compares like with like. Effect on spheres in a disk, a cylinder and a sphere (N = 500, 6 seeds, three
size distributions, and a 4D container): the mean packing fraction moves by -0.0002 to +0.0009, less than the seed
scatter in every case (largest 0.78 of it); equal-size disks in a disk finish in 17,110 steps instead of 20,526.
Boxes are not affected.

### Changed: quench
- Rods are balanced as slightly rounded rods (`rod_round=0.01`). The solver of that rounded contact could stop
  before it had converged for some crossed rods, which made two of 78 survey quenches fail; repaired (78 of 78).
- A quench that stalls or fails: quench the packed state again with a lower `dt_max`, or with slightly more
  compression (`kappa_mult`); resuming with unchanged settings does not help a quench that is stuck. Flat-ended
  rods in 3D inside a curved container often need this at the defaults. README, "Quench to force balance".

### Changed: per-particle lists
With your own diameter list (`dist={"type": "custom", ...}`), entry `i` of a per-particle shape table or rod list
goes with entry `i` of that diameter list. With generated diameters, row `i` goes to the `i`-th largest particle.
