# Skew Cases

This file records finite-skew validation cases for `python/experimental/python/`.

## 2026-03-18

- tree: `python/experimental/python/`
- scope covered so far:
  - 2D periodic mild skew
  - 2D periodic guardrail correction
  - 2D rendering path for skewed periodic cells

### Case: 2D periodic mild skew

- config:
  - `N=32`
  - `Ndim=2`
  - `phi=0.11`
  - `box=[1.0, 1.0]`
  - `walls=[0, 0]`
  - `box_skew=[0.0, -20.0]`
  - `dist={'type': 'mono', 'd': 1.0}`
  - `seed=123`
- result:
  - converged without NaNs
  - corrected `box_skew` remained `[0.0, -20.0]`
  - `phi_final=0.8692422996623266`
  - `force_magnitude=0.0006202817109434772`
  - `steps=10744`

### Case: 2D periodic guardrail

- config:
  - same as mild skew except requested `box_skew=[0.0, -75.0]`
- result:
  - initializer warning emitted
  - corrected `box_skew=[0.0, -59.90000000000057]`
  - packing remained finite
  - `phi_final=0.861896653478047`
  - `force_magnitude=0.004548571856850665`
  - `steps=21518`

### Notes

- zero-skew periodic remains finite in the same experimental tree
- the first skew probe failed before reinstall because Python was still loading a stale installed extension module
- an experimental C++ regression test now exists at `python/experimental/python/tests/skew_periodic_regression.cpp`
