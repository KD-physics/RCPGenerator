# Changelog

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
