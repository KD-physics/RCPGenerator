# Baseline Results

This file records baseline correctness and representative output metrics for the current `python/` tree before skew-boundary work.

## Commands

- Build:
  - `cmake -S python -B python/build-py -G Ninja -DBUILD_PARITY_TESTS=ON -Dpybind11_DIR=...`
- Test:
  - `ctest --test-dir python/build-py --output-on-failure`
- Python package:
  - `python -m pip install -e python --no-build-isolation`

## Current Recorded Results

### C++ Tests

- tree: `python/`
- command: `ctest --test-dir python/build-py --output-on-failure`
- result: `3/3 tests passed`
- tests:
  - `rcpgenerator_initializer_parity`
  - `rcpgenerator_packer_parity`
  - `rcpgenerator_hard_wall_regression`

### Python 2D Boundary Probe

- tree: `python/`
- configuration: `phi=0.11, N=64, Ndim=2, box=[1.0,1.0], dist=mono(d=1.0), neighbor_max=0, seed=123`

Recorded outputs:

- `walls=[0,0]`
  - initialization finite: yes
  - pack finite: yes
  - `phi=0.8397127786333479`
  - `force=0.0022241976896253905`
  - `steps=36008`

- `walls=[-1,0]`
  - initialization finite: yes
  - pack finite: yes
  - `phi=0.6171180902126232`
  - `force=nan`
  - `steps=60000`
  - note: this probe uses `[-1,0]`, while the shipped circular example uses `[-2,0]`

- `walls=[0,1]`
  - initialization finite: yes
  - pack finite: yes
  - `phi=0.8226183402306758`
  - `force=0.00036328990442205904`
  - `steps=8737`

- `walls=[1,1]`
  - initialization finite: yes
  - pack finite: yes
  - `phi=0.7734388611881311`
  - `force=0.004753444541347523`
  - `steps=16202`

### Python Example-Scale Hard-Wall Probe

- tree: `python/`
- configuration: `phi=0.11, N=250, Ndim=2, box=[1.0,1.0], walls=[1,1], dist=mono(d=1.0), neighbor_max=0, seed=123`
- initialization finite: yes
- final finite: yes
- `phi=0.8294609996921389`
- `force=0.00028931377560504354`
- `steps=19470`

## Baseline Timing Snapshot

### 2026-03-18

- tree: `python/`
- measurement method: end-to-end `Packing(...)` construction plus `pack()`

- case: `2d_periodic_64`
  - command shape: Python one-off probe
  - skew: none
  - runtime_s: `0.238824`
  - `phi=0.835359496846431`
  - `force=0.0008072806565487041`
  - `steps=16978`

- case: `2d_two_hard_64`
  - command shape: Python one-off probe
  - skew: none
  - runtime_s: `0.407780`
  - `phi=0.7950366058386629`
  - `force=0.0014628739013737887`
  - `steps=26774`

- case: `2d_two_hard_250`
  - command shape: Python one-off probe
  - skew: none
  - runtime_s: `1.639498`
  - `phi=0.8248721719803589`
  - `force=0.0004036556248876518`
  - `steps=25609`

- case: `3d_periodic_64`
  - command shape: Python one-off probe
  - skew: none
  - runtime_s: `0.396154`
  - `phi=0.625645841533886`
  - `force=0.0014977944725883155`
  - `steps=10438`

## Experimental Tree Validation Snapshot

### 2026-03-18

- tree: `python/experimental/python/`
- test command: `ctest --test-dir python/experimental/python/build-exp --output-on-failure`
- result: `3/3 tests passed`

Direct Python probes from the experimental tree:

- case: `exp_2d_periodic_64`
  - runtime_s: `0.097186`
  - `phi=0.8371066713193166`
  - `force=0.0005357810263300632`
  - `steps=7851`
  - finite: yes

- case: `exp_2d_two_hard_64`
  - runtime_s: `0.186161`
  - `phi=0.789772797120563`
  - `force=0.00400507785100868`
  - `steps=15456`
  - finite: yes

Notes:

- These direct Python probe values differ from the baseline snapshot because initialization remains non-deterministic through `std::random_device`.
- The C++ parity tests remain the reliable zero-change regression gate for deterministic behavior.
