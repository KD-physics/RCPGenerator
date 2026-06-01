# Performance Results

This file records runtime comparisons between the baseline tree and the future experimental skew-capable tree.

## Status

- No baseline-vs-experimental comparison has been recorded yet.
- Baseline runtime capture should be added before modifying the experimental tree.

## Intended Table

For each case, record:

- date/time
- tree
- case name
- command
- skew setting
- runtime
- repeated-run notes
- key output metrics

## Initial Baseline Snapshot

### 2026-03-18

- tree: `python/`
- skew setting: zero / disabled

- case: `2d_periodic_64`
  - runtime_s: `0.238824`
  - `phi=0.835359496846431`
  - `force=0.0008072806565487041`
  - `steps=16978`

- case: `2d_two_hard_64`
  - runtime_s: `0.407780`
  - `phi=0.7950366058386629`
  - `force=0.0014628739013737887`
  - `steps=26774`

- case: `2d_two_hard_250`
  - runtime_s: `1.639498`
  - `phi=0.8248721719803589`
  - `force=0.0004036556248876518`
  - `steps=25609`

- case: `3d_periodic_64`
  - runtime_s: `0.396154`
  - `phi=0.625645841533886`
  - `force=0.0014977944725883155`
  - `steps=10438`

## Experimental Zero-Change Snapshot

### 2026-03-18

- tree: `python/experimental/python/`
- skew setting: zero / disabled

- case: `exp_2d_periodic_64`
  - runtime_s: `0.097186`
  - `phi=0.8371066713193166`
  - `force=0.0005357810263300632`
  - `steps=7851`

- case: `exp_2d_two_hard_64`
  - runtime_s: `0.186161`
  - `phi=0.789772797120563`
  - `force=0.00400507785100868`
  - `steps=15456`

Notes:

- Do not interpret the direct Python runtime and result differences as a geometry regression by themselves.
- Current end-to-end Python initialization is non-deterministic, so these measurements are suitable only as rough runtime snapshots.
- Deterministic regression should continue to rely on the C++ parity tests unless initialization seeding is made explicit in the experimental API.

## Experimental Skew Snapshot

### 2026-03-18

- tree: `python/experimental/python/`
- command: targeted Python probe via installed experimental package

- case: `exp_2d_periodic_32_zero_skew`
  - skew setting: `box_skew=[0.0, 0.0]`
  - runtime_s:
    - init: `0.000410`
    - pack: `0.036751`
  - `phi=0.8346338128981675`
  - `force=0.0006616441923446396`
  - `steps=5378`

- case: `exp_2d_periodic_32_mild_skew`
  - skew setting: `box_skew=[0.0, -20.0]`
  - runtime_s:
    - init: `0.000162`
    - pack: `0.217313`
  - `phi=0.8692422996623266`
  - `force=0.0006202817109434772`
  - `steps=10744`

- case: `exp_2d_periodic_32_guardrail`
  - skew setting:
    - requested: `box_skew=[0.0, -75.0]`
    - corrected: `box_skew=[0.0, -59.90000000000057]`
  - runtime_s:
    - init: `0.000529`
    - pack: `0.424291`
  - `phi=0.861896653478047`
  - `force=0.004548571856850665`
  - `steps=21518`

Notes:

- the first finite-skew Python probe initially failed because the installed package was stale relative to the rebuilt experimental tree
- reinstalling `python/experimental/python/` resolved that mismatch

## Experimental Skew Reverification

### 2026-03-18

- tree: `python/experimental/python/`
- command: targeted Python probe after fresh editable reinstall of the experimental package

- case: `exp_2d_periodic_32_zero_skew_reverified`
  - skew setting: `box_skew=[0.0, 0.0]`
  - `phi=0.8289107602609177`
  - `force=0.0007183038147209226`
  - `steps=8230`

- case: `exp_2d_periodic_32_mild_skew_reverified`
  - skew setting: `box_skew=[0.0, -20.0]`
  - `phi=0.8692422996468798`
  - `force=0.0005308890768692847`
  - `steps=5764`

- case: `exp_2d_periodic_32_guardrail_reverified`
  - skew setting:
    - requested: `box_skew=[0.0, -75.0]`
    - corrected: `box_skew=[0.0, -59.90000000000057]`
  - `phi=0.8618598323195438`
  - `force=0.0003941107040262731`
  - `steps=9956`

Notes:

- these values differ from earlier snapshots because end-to-end initialization is still non-deterministic
- the important verified behavior is that all cases remain finite and converged, and the guardrail correction is preserved after a fresh build/install cycle

## Baseline Vs Experimental Timing Comparison

### 2026-03-18

- objective:
  - compare original baseline `python/` runtime against the experimental tree at zero skew and finite skew
- case:
  - `N=400`
  - `Ndim=2`
  - periodic box
  - `phi=0.11`
  - `dist={"type": "mono", "d": 1.0}`
  - `walls=[0,0]`
  - repeated `5` times per case

- case: `baseline_zero_skew_400`
  - tree: `python/`
  - skew setting: `box_skew=[0.0, 0.0]`
  - `mean_init_s=0.000865`
  - `mean_pack_s=1.149757`
  - `mean_total_s=1.150622`
  - `total_stdev_s=0.338242`
  - `mean_steps=12817`

- case: `experimental_zero_skew_400`
  - tree: `python/experimental/python/`
  - skew setting: `box_skew=[0.0, 0.0]`
  - `mean_init_s=0.000906`
  - `mean_pack_s=1.813930`
  - `mean_total_s=1.814836`
  - `total_stdev_s=1.033650`
  - `mean_steps=20644.6`
  - versus baseline:
    - `total_ratio=1.5773x`
    - `pack_ratio=1.5777x`

- case: `experimental_mild_skew_400`
  - tree: `python/experimental/python/`
  - skew setting: `box_skew=[0.0, -20.0]`
  - `mean_init_s=0.012072`
  - `mean_pack_s=13.793923`
  - `mean_total_s=13.805995`
  - `total_stdev_s=2.364408`
  - `mean_steps=17105.2`
  - versus baseline:
    - `total_ratio=11.9987x`
    - `pack_ratio=11.9972x`
  - versus experimental zero skew:
    - `total_ratio=7.6074x`
    - `pack_ratio=7.6044x`

Notes:

- this is the first direct baseline-vs-experimental timing comparison recorded in the log
- the numbers are useful for coarse runtime comparison, but they are not a strict algorithmic microbenchmark
- initialization is still non-deterministic, so the realized particle state and resulting step count vary noticeably between runs
- the current data nevertheless suggests:
  - the experimental zero-skew path is materially slower than the original baseline on this `N=400` case
  - enabling mild skew is much more expensive again than experimental zero skew
