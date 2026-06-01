# Skew Boundary Handoff

Use this file as the restart point for a new Codex session.

## Repo Context

- Repo root: `C:\dev\RCPGenerator_build\c++`
- Baseline implementation: `python/`
- Experimental implementation: `python/experimental/python/`
- All skew-boundary work is in `python/experimental/python/`

## Goal

Add experimental skewed periodic boundary support with minimal disruption to the current `python/` baseline.

Current public API direction:

- `box_skew=[a, b, ...]`
- length matches `box`
- default all zeros
- x is the anchor direction
- guardrail: each non-x basis direction must remain at least `30 degrees` from x
- if violated:
  - keep x fixed
  - adjust the violating direction
  - warn the user
  - store the corrected `box_skew` back into the result/state

## Files Already Changed

Experimental code changes already exist in:

- `python/experimental/python/include/rcpgenerator/initialize_particles.hpp`
- `python/experimental/python/include/rcpgenerator/rcp_generator.hpp`
- `python/experimental/python/src/core/initialize_particles.cpp`
- `python/experimental/python/src/core/rcp_generator.cpp`
- `python/experimental/python/src/bindings/module.cpp`
- `python/experimental/python/python/rcpgenerator/_wrapper.py`
- `python/experimental/python/python/rcpgenerator/render.py`
- `python/experimental/python/CMakeLists.txt`
- `python/experimental/python/tests/skew_periodic_regression.cpp`

Planning/log files already updated:

- `python/planning/skew_boundary_plan.md`
- `python/planning/logs/session_log.md`
- `python/planning/logs/performance_results.md`
- `python/planning/logs/skew_cases.md`

## What Works

The experimental runtime path and the fresh experimental build/test path have both now been validated.

Validated from Python after installing the experimental package:

- zero skew:
  - `walls=[0,0]`
  - `box_skew=[0.0, 0.0]`
  - finite and converged

- finite skew:
  - `walls=[0,0]`
  - `box_skew=[0.0, -20.0]`
  - finite and converged

- guardrail:
  - requested `box_skew=[0.0, -75.0]`
  - corrected to about `box_skew=[0.0, -59.9]`
  - finite and converged

Observed Python results from the working experimental install:

- zero skew:
  - `phi_final=0.8346338128981675`
  - `force_magnitude=0.0006616441923446396`
  - `steps=5378`

- mild skew:
  - `phi_final=0.8692422996623266`
  - `force_magnitude=0.0006202817109434772`
  - `steps=10744`

- guardrail case:
  - `phi_final=0.861896653478047`
  - `force_magnitude=0.004548571856850665`
  - `steps=21518`

Freshly reverified from Python after reinstalling the experimental package:

- zero skew:
  - `phi_final=0.8289107602609177`
  - `force_magnitude=0.0007183038147209226`
  - `steps=8230`

- mild skew:
  - `phi_final=0.8692422996468798`
  - `force_magnitude=0.0005308890768692847`
  - `steps=5764`

- guardrail case:
  - corrected to `box_skew=[0.0, -59.90000000000057]`
  - `phi_final=0.8618598323195438`
  - `force_magnitude=0.0003941107040262731`
  - `steps=9956`

## Important Debugging Detail

At one point Python was loading a stale installed package, which caused confusing results.

After reinstalling the experimental package with:

```powershell
python -m pip install -e python/experimental/python --no-build-isolation
```

the finite-skew runtime path behaved correctly.

## Build Status

The previously reported configure blocker turned out to be a sandbox restriction, not a demonstrated skew-code defect.

Fresh configure/build/test succeeded when run unsandboxed in:

- `python/experimental/python/build-exp-live`

Verified commands:

```powershell
cmake -S python/experimental/python -B python/experimental/python/build-exp-live -G Ninja -DBUILD_PARITY_TESTS=ON -Dpybind11_DIR=C:/Users/kend1/anaconda3/envs/rcp/Lib/site-packages/pybind11/share/cmake/pybind11
cmake --build python/experimental/python/build-exp-live --target _rcpgenerator rcpgenerator_skew_periodic_regression rcpgenerator_initializer_parity rcpgenerator_packer_parity rcpgenerator_hard_wall_regression
ctest --test-dir python/experimental/python/build-exp-live --output-on-failure
```

Verified `ctest` result:

- `4/4` tests passed
- `rcpgenerator_skew_periodic_regression` is now confirmed in the generated test manifest and passes

The stale `python/experimental/python/build-py/` directory should not be used as evidence for the experimental skew build state because its generated files still point at the baseline `python/` source tree.

## Recommended Next Steps For New Session

1. Read:
   - `python/planning/skew_boundary_plan.md`
   - `python/planning/logs/session_log.md`
   - `python/planning/skew_boundary_handoff.md`

2. Inspect the experimental source changes in:
   - `python/experimental/python/src/core/initialize_particles.cpp`
   - `python/experimental/python/src/core/rcp_generator.cpp`
   - `python/experimental/python/python/rcpgenerator/_wrapper.py`
   - `python/experimental/python/python/rcpgenerator/render.py`

3. Use the verified experimental build directory as the current reference point:

```powershell
cmake -S python/experimental/python -B python/experimental/python/build-exp-live -G Ninja -DBUILD_PARITY_TESTS=ON -Dpybind11_DIR=C:/Users/kend1/anaconda3/envs/rcp/Lib/site-packages/pybind11/share/cmake/pybind11
cmake --build python/experimental/python/build-exp-live --target _rcpgenerator rcpgenerator_skew_periodic_regression
ctest --test-dir python/experimental/python/build-exp-live --output-on-failure
```

4. If future configure/build attempts fail only inside the Codex sandbox with temp-file deletion errors, treat that as a sandbox/environment issue first, not a skew-code issue.

5. Continue from the verified state:
   - decide whether to add more skew-specific regression cases beyond the current mild-skew and guardrail coverage
   - review whether the current generalized neighbor search path is sufficient for larger skew sweeps
   - decide whether any documentation/examples should be added for the experimental API

## Suggested Prompt For New Session

Use something like:

`Look at python/planning/skew_boundary_handoff.md and familiarize yourself with the code base, where we are in the build, and what to do next.`
