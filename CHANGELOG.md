# Changelog

Completed tasks, newest entry at the top, dated. `specdev task done` adds an entry for each finished task file; work done directly in `specs/ctx.md` is moved here once every checkbox in its plan is done. Coexists with any conventional release notes already here.

<!-- Entry template:
## YYYY-MM-DD — <task title>

- what shipped
- files touched / decisions locked
-->

## 2026-10-12 — range contract std dev

- Sensor std-dev demo moved to `examples/sensor_std_dev/` (`SensorContract` + `encrypted_std_dev`); library keeps generic ops with `fit_chebyshev`; i64-era scalar clamp removed.
- Task `0002-range-contract-std-dev`; source: user request 2026-10-12: revive project, CKKS green + native sensor demo with range normalization

## 2026-10-12 — encoder wide staging

- Encoder stages coefficients in f64 and reduces them straight into RNS (no i64 cap, no silent saturation); `try_encode` errors past Q/2; decode reconstructs in f64.
- Task `0004-encoder-wide-staging`; source: user request 2026-10-12: encoder i64 staging caps |x| < 4 and saturates silently; widen like SEAL, keep f64 FFT

## 2026-10-12 — green baseline

- Green baseline on rustc 1.99: removed scratch example, fixed 2 clippy errors, applied cargo fmt, fixed stale AGENTS.md commands.
- Task `0001-green-baseline`; source: user request 2026-10-12: revive project, CKKS green + native sensor demo with range normalization

## 2026-08-17 — Demo pre-requisites, Phase A steps 1–4

- `DIGIT_BITS` 15 -> 20 (`src/crypto/engine.rs`); `std_dev_8` variance error 2.5e-13 -> 7.1e-12, inside the 1e-9 budget.
- `mul_plaintext` engine op (ct x pt, no relin) and test `mul_plaintext_matches_ct_mul`.
- Pipeline ops moved into `src/crypto/pipeline.rs` (`sum_slots`, `mul_plain_scalar`, `align_level`, `reduce_secret_key`, `decode_wide`, `fit_sqrt_chebyshev`, `eval_chebyshev_clenshaw`); `std_dev_8` uses them; `RnsBasis::reconstruct_centered_coeff_f64`; `Ciphertext: Clone`.
- Encrypted Chebyshev eval (Clenshaw) and `LeveledEvaluator` (per-level gadget-key cache); matches plaintext Clenshaw within 1e-3 at toy params.
- Left red: `encrypted_std_dev_end_to_end` (direct fit on [1, 100] wraps constants past q0); fixed by task 0002.
