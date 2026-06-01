# Skew Boundary Plan

## Purpose

This document is the working plan for adding skewed-cell boundary support in a new experimental Python refactor tree while preserving a direct baseline against the current `python/` implementation.

It is intended to serve as the restart point if the session is lost.

## Scope

- No changes to the top-level legacy reference files.
- Baseline implementation remains the current `python/` tree.
- Experimental work will be done in `python/experimental/python/`.
- The goal is to add a skew-boundary capability with as few behavioral changes as possible outside the boundary math.

## High-Level Strategy

1. Preserve the current `python/` tree as the local baseline.
2. Use `python/experimental/python/` as the local experimental copy.
3. Add skew-boundary support in `python/experimental/python/` only.
4. Keep a zero-skew mode in `python/experimental/python/` that should reproduce `python/` behavior as closely as possible.
5. Compare correctness and runtime between `python/` and `python/experimental/python/` at zero skew.
6. Add finite-skew tests and examples only after zero-skew parity is established.

## Design Intent

The first implementation should prefer the smallest coherent geometry change rather than a broad redesign.

Expected design direction:

- Add an explicit skew/basis representation rather than encoding skew indirectly through ad hoc trigonometric corrections.
- Preserve the existing axis-aligned fast path when skew is disabled or zero.
- Route to generalized cell math only when skew is active.
- Keep Python-facing APIs as close as possible to the current interface.

## Working Assumptions

- Skew support is expected to matter first for periodic boundaries.
- Hard-wall skew containers may be added, but they should be treated as a separate geometric problem from periodic skew wrapping.
- The current code should remain the performance baseline for axis-aligned cases.
- Zero-skew comparison against the baseline does not need to be bitwise identical.
- The practical baseline standard is:
  - same qualitative behavior
  - reasonably close final `phi`
  - reasonably close runtime
  - no obvious convergence or stability regression
- User-facing seed reproducibility should be preserved as much as possible.

## Proposed Execution Phases

### Phase 0: Baseline Freeze

Goal:
- Establish the current `python/` behavior and runtime as the comparison target.

Tasks:
- Record the exact baseline examples and tests that currently pass.
- Record the current build and test commands.
- Record representative runtimes for selected examples.
- Record representative output metrics such as final `phi`, `steps`, and `force_magnitude`.

Required baseline cases:
- 2D periodic box
- 2D one-hard-wall case
- 2D two-hard-wall case
- 2D circular container
- 3D periodic box
- Any existing parity tests

Artifacts to preserve:
- test results
- runtime measurements
- selected example outputs
- environment notes if relevant

### Phase 1: Experimental Tree Setup

Goal:
- Use `python/experimental/python/` as a local experimental copy of `python/`.

Tasks:
- Confirm `python/experimental/python/` builds and runs identically before any edits.
- Run the baseline test suite in both trees and compare outputs.

Success criterion:
- `python/experimental/python/` at zero changes matches `python/` behavior and is within expected runtime variation.

### Phase 2: Zero-Skew Infrastructure

Goal:
- Introduce the new boundary/skew API and internal representation in `python/experimental/python/` without changing behavior when skew is zero.

Tasks:
- Decide on the user-facing skew representation.
- Decide whether skew is represented by:
  - cell basis vectors
  - a reduced skew parameter
  - tilt/shear factors
- Add parsing, state storage, and validation for skew configuration.
- Keep the existing axis-aligned code path active when skew is disabled or zero.

Success criterion:
- All zero-skew tests still match baseline behavior.

### Phase 3: Periodic Skew Cell

Goal:
- Implement skewed periodic boundary handling in `python/experimental/python/`.

Tasks:
- Add cell-matrix and inverse-cell handling.
- Convert wrapping and minimum-image calculations to fractional-coordinate logic for skew-enabled runs.
- Update neighbor search logic if current assumptions depend on axis alignment.
- Add focused low-particle tests where expected wrapping behavior is easy to inspect.

Success criterion:
- Finite-skew periodic cases run stably.
- Zero-skew cases remain baseline-compatible.

### Phase 4: Hard-Wall Skew Container

Goal:
- Decide whether skewed hard-wall boundaries are in scope for the first version.

Tasks:
- If included, define wall equations and inside/outside tests.
- Add particle-wall distance logic for skewed walls.
- Validate against simple hand-checkable cases.

Success criterion:
- Hard-wall skew logic is stable and does not perturb zero-skew behavior.

### Phase 5: Performance Evaluation

Goal:
- Measure the runtime overhead of the generalized geometry path.

Tasks:
- Benchmark `python/` vs `python/experimental/python/` on zero-skew cases.
- Benchmark `python/experimental/python/` with zero skew using:
  - axis-aligned fast path
  - generalized skew-capable path if both exist
- Benchmark `python/experimental/python/` with representative finite skew.

Performance questions to answer:
- Is zero-skew overhead negligible?
- Is a dedicated axis-aligned fast path still justified?
- Does skew materially affect neighbor search cost?

### Phase 6: Decision Point

Goal:
- Decide whether the experimental implementation is ready to replace `python/`.

Criteria:
- zero-skew baseline behavior is preserved
- finite-skew cases are stable
- performance cost is acceptable
- tests and examples are sufficient for maintenance

## Logging Requirements

Future sessions should maintain a simple append-only log under:

- `python/planning/logs/`

Suggested files:

- `python/planning/logs/session_log.md`
- `python/planning/logs/baseline_results.md`
- `python/planning/logs/performance_results.md`
- `python/planning/logs/skew_cases.md`

### Session Log Format

Each work session should append:

- date/time
- branch or folder context used
- goals for the session
- commands run
- files changed
- tests run
- results
- blockers or open questions
- next recommended step

Suggested template:

```md
## Session YYYY-MM-DD HH:MM

Context:
- baseline tree:
- experimental tree:

Goals:
- ...

Commands run:
- ...

Files changed:
- ...

Results:
- ...

Blockers:
- ...

Next step:
- ...
```

### Baseline Results Format

For each baseline case, record:

- case name
- tree used
- command run
- runtime
- final phi
- steps
- force magnitude
- pass/fail
- notes

### Performance Results Format

For each benchmark, record:

- case name
- tree used
- skew setting
- runtime
- repeated-run statistics if measured
- notes on variance

## Minimum Test Matrix

### Baseline Compatibility

- `python/` current parity tests
- `python2/` parity tests at zero skew
- representative Python examples in both trees at zero skew

### Zero-Skew Regression Cases

- 2D `walls=[0,0]`
- 2D `walls=[0,1]`
- 2D `walls=[1,1]`
- 2D curved container case
- 3D periodic case

### Finite-Skew Cases

Start with small deterministic tests:

- 2D periodic box with mild skew
- 2D periodic box with stronger skew
- simple particle pairs near wrap boundaries
- optional 3D skew case if design generalizes beyond 2D immediately

### Hard-Wall Finite-Skew Cases

Only if included in first implementation:

- 2D parallelogram hard container with a small number of particles
- particle near each wall
- particle near corner intersections

## Open Design Questions

- What exact Python API should enable skew?
- Is skew periodic-only in the first version, or should hard-wall skew also be included immediately?
- Should the skew representation be angle-based or basis-vector-based?
- Should zero skew route to existing code or always use generalized code?
- Will neighbor-list logic need a full rewrite for skewed periodic geometry, or only a local adaptation?

## Initial API Proposal

This section captures the current recommended direction before implementation.

### First-Version Scope

Recommended first version:

- periodic skew only
- 2D first
- preserve the existing non-skew code path
- do not include skewed hard-wall containers in the first pass

Reasoning:

- periodic skew is the cleaner and more standard geometry problem
- hard-wall skew adds a separate wall-distance and corner-handling problem
- keeping the axis-aligned fast path reduces regression risk and simplifies performance comparison

### Public Configuration Shape

Recommended first-pass user-facing API:

```python
packing = rcpgenerator.Packing(
    ...,
    box_skew=[a_deg, b_deg, c_deg, ...],
)
```

Interpretation:

- `box_skew`
  - list of per-basis-direction angles in degrees
  - length must match `box`
  - default is `[0, 0, ..., 0]`
- each entry is interpreted as a counter-clockwise rotation relative to the corresponding axis direction
- the design is not restricted to 2D
- `box[0]` / x is the anchor basis direction used for guardrail checks and auto-correction
- for 2D:
  - `box[0]` defines the magnitude of the first basis vector, whose reference direction is the x-axis
  - `box[1]` defines the magnitude of the second basis vector, whose reference direction is the y-axis

Rationale:

- both basis directions are intentionally user-adjustable
- this is not limited to simple shear
- the purpose includes frustrating crystalline order along both directions
- the same per-direction idea is intended to carry across dimensions, with x acting as the anchor

### Zero-Skew Behavior

When `box_skew` is omitted or all entries are zero:

the code should route to the current axis-aligned logic and preserve present behavior as closely as practical.

### Internal Representation

Recommended internal representation:

- convert the user-facing `box` and `box_skew` values into a cell basis
- precompute the cell matrix and its inverse
- do periodic wrapping and minimum-image logic in fractional coordinates only when skew is enabled

Suggested 2D basis form:

- `a = Lx * [cos(ax), sin(ax)]`
- `b = Ly * [cos(90 deg + ay), sin(90 deg + ay)]`

This keeps the public API simple while preserving a clean path to future generalization.

### Guardrail on Basis Angle

To avoid extremely sharp cells that may be numerically awkward or difficult to converge:

- in any dimension, the angle between the x anchor basis vector and each other basis vector should be constrained to be at least `30 degrees`
- if the requested geometry violates that guardrail:
  - keep the x-direction basis vector fixed
  - adjust the violating basis direction just enough to make the basis angle `30 degrees`
  - emit a warning to the user

Requested behavior after correction:

- the corrected value should be pushed back into the stored `box_skew`
- once packing is done, the user should see the latest corrected `box_skew` values rather than the invalid original request
- rendering methods, including `show_packing`, should display the corrected skewed cell geometry consistently

### Future Extension Path

If the first version succeeds, a later generalized API could become:

```python
box_basis=[[ax, ay], [bx, by]]
```

or a 3D analogue.

The first implementation should not jump to that broader form unless it is required immediately.

### Multidimensional Clarification

`box_skew` is intended to be available in any supported packing dimension, not only in 2D.

Important limitation:

- one scalar angle per basis direction does not describe the most general skewed cell in dimensions above 2
- instead, it defines a constrained family of direction adjustments relative to the original axis-aligned basis

That limitation is acceptable for the current objective as long as it is documented clearly.

## Recommended Next Step

Before any implementation:

1. Review and approve this plan.
2. Decide the initial scope:
   - periodic skew only
   - or periodic plus hard-wall skew
3. Decide the preferred public representation of skew:
   - angle/tilt parameter
   - or explicit basis vectors
4. After approval, validate `python/experimental/python/` as the experimental baseline and continue with zero-skew comparison work.
