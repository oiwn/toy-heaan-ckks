# Current Task Context: Sensor std-dev end-to-end (CKKS green + native demo)

State: in progress. Design steps below come first, then task files 0001 -> 0002 -> 0003. Roadmap M1 + M2.

## Plan

- [x] Range contract: fix the fields, the demo values, and how out-of-range readings are clamped.
- [x] Degree vs accuracy: choose the Chebyshev degree and set the std error budget for the demo.
- [x] Level/prime budget: recount levels for the chosen degree and set the prime count for N = 8192.
- [x] Record the decisions in `specs/overview.md` (gotchas) and `specs/demo.md`; move tasks 0002-0003 from draft to ready.
- [x] Task 0001-green-baseline (fmt, clippy, scratch cleanup).
- [x] Task 0004-encoder-wide-staging (no silent i64 cap in encode/decode).
- [x] Task 0002-range-contract-std-dev (fixes the red test).
- [ ] Task 0003-sensor-demo-example (native N = 8192 demo).

## Findings

- **Why `encrypted_std_dev_end_to_end` was red** (wrong sign: -2.2256 vs 2.2214): the direct fit on variance [1, 100] has c0 ~ 6.5, and `6.5 * 2^61 > 2^63` saturated in the encoder's old i64 staging. Confirmed: task 0004 alone turns the test green. The "wraps past q0" diagnosis was wrong.
- **Encoder fixed (task 0004):** coefficients are staged in f64 and reduced straight into RNS; the limit is now |coeff| < 2^(total_bits - 1) (~Q/2), with `try_encode` returning an error instead of saturating. Precision is ~2^-52 relative to the largest coefficient.
- **Small values cost approximation accuracy, not CKKS precision.** sqrt near 0 is hard for a polynomial; CKKS noise stays ~1e-12 absolute.
- **Decided 2026-10-12 — sensor range contract** (replaces the 1/B normalization):
  - `SensorContract { lo: -40.0, hi: 60.0, std_min: 1.0, std_max: 8.0 }`, `W = hi - lo = 100`.
  - Client clamps readings to `[lo, hi]`, sends `x_n = (x - lo) / W` (in [0, 1], under the encoder's i64 cap), and asserts its plaintext std is inside `[std_min, std_max]`.
  - Fit sqrt on `[(std_min/W)^2, (std_max/W)^2]` = [1e-4, 6.4e-3]; max coefficient ~0.05, so nothing wraps.
  - Client restores `std = std_n * W` in f64 after decrypt: no 1/B level, no encrypted x10.
- **Why a std envelope, not the 0.25 Popoviciu bound:** sqrt is homogeneous, so the error in degC depends only on `std_max/std_min` and the degree, not on `W`. Sweep (max error over the envelope): [1, 8] degree 7 -> 0.055, degree 9 -> 0.025; [1, 10] degree 7 -> 0.107; [0.5, 8] degree 7 -> 0.170.
- **Realistic data:** outdoor daily swing 3-20 degC, daily std ~0.35 x swing = 1-7 degC; indoor 0.2-1.5 degC (outside the envelope).
- **Population std (divide by n)** is the computed statistic, kept deliberately (operator, 2026-10-12); sample std would only change the `1/n` constant to `1/(n-1)`.
- **Degree 7, budget:** std error < 0.1 degC (expect <= 0.055), intermediate variance error < 1e-9.
- **Levels:** 3 + (degree + 1) = 11; the final Clenshaw mul needs >= 3 channels, so 13 x 61-bit primes at N = 8192.
- **For 0002:** the affine slope `4 / (v_hi - v_lo)` ~ 630 no longer hits an encoder cap; drop `clamped_scalar_bits` and the i64 asserts in `pipeline.rs`, and keep the range contract for sqrt accuracy (the std envelope), not for encoder limits.

## Context

Next task: `specs/tasks/0003-sensor-demo-example.md`. Demo code: `examples/sensor_std_dev/`; generic fitter `fit_chebyshev` in `src/crypto/pipeline.rs`.

## Next

Start task 0003 (`specdev task advance 3`).
