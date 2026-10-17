# Changelog

Completed tasks, newest entry at the top, dated. `specdev task done` adds an entry for each finished task file; work done directly in `specs/ctx.md` is moved here once every checkbox in its plan is done. Coexists with any conventional release notes already here.

<!-- Entry template:
## YYYY-MM-DD — <task title>

- what shipped
- files touched / decisions locked
-->

## 2026-10-16 — generic key switching

- One generic `KeySwitchingKey<P>` (gadget encryption of an old key, `switch(d)`) now backs `RelinearizationKey` and `RotationKey`; the legacy single-pair keys, `RnsGadget*` types and `mul_ciphertexts` are gone, and ops, evaluator, algorithms and engine are generic over `PolyRing`.
- Task `0012-generic-key-switching`; source: user request 2026-10-13: pluggable ring backends, scheme over a small trait (overview.md Target architecture)

## 2026-10-16 — backend trait

- One backend trait `PolyRing` (arithmetic, sampling, automorphism, modulus chain, encoding boundary, gadget) replaces the three ring traits; `RnsPoly` rests in NTT form behind a private flag, and the encoder, `Ciphertext` level methods and gadget key switching use only the trait. Output bit-identical, sensor pipeline 24 s -> 9 s.
- Task `0011-backend-trait`; source: user request 2026-10-13: pluggable ring backends, scheme over a small trait (overview.md Target architecture)

## 2026-10-15 — remove commented code

- Removed dead commented-out code from the old `RnsPolyRing` API: `keys/common.rs` and `benches/basic.rs` (both moved to `tmp/`), stale tests in `math/utils.rs` and two `CkksError` variants.
- Task `0014-remove-commented-code`; source: operator 2026-10-14: huge commented-out sections in keys/common.rs

## 2026-10-14 — chebyshev module

- `pipeline.rs` retired: `sum_slots` and the Chebyshev code (encrypted Clenshaw, fit, plaintext reference) now live in `ckks/algorithms/{sum_slots, chebyshev}.rs`, built only on ops, level methods and the evaluator; new `Ciphertext::zero_like` op.
- Task `0010-chebyshev-module`; source: user request 2026-10-13: pipeline.rs is a std-dev grab bag; decompose into composable ops like fhe.rs / tfhe-rs

## 2026-10-14 — ops composition

- Homomorphic ops moved to `ckks/ops.rs`: `std::ops` Add/Sub/Neg plus `Ciphertext::{mul, mul_plain, mul_const, add_const, rotate}` with explicit rescale; `LeveledEvaluator` and rotation-key helpers moved to `ckks/evaluator.rs`; the `pipeline.rs` scalar helpers and the engine's op functions are gone.
- Task `0009-ops-composition`; source: user request 2026-10-13: pipeline.rs is a std-dev grab bag; decompose into composable ops like fhe.rs / tfhe-rs

## 2026-10-14 — ciphertext levels

- Shared modulus chain: `RnsBasis` links to the next level, `Ciphertext` gains `level`/`logq()`/`rescale`/`mod_switch_to`/`match_level`, `logq` is no longer stored, and `align_level`, `decode_wide`, `rescale_ciphertext` and the unimplemented level traits are gone.
- Task `0008-ciphertext-levels`; source: user request 2026-10-13: pipeline.rs is a std-dev grab bag; decompose into composable ops like fhe.rs / tfhe-rs

## 2026-10-14 — module layout

- Moved to the `foo.rs` + `foo/` layout: `crypto`, `encoding`, `keys` now live under `ckks/`; `math` and `rings/backends/rns_ntt` lost their `mod.rs`; import paths and docs updated, no logic change.
- Task `0007-module-layout`; source: user request 2026-10-13: publication cleanup, remove ancient code and golden vectors, foo.rs module layout

## 2026-10-14 — trim dependencies

- Trimmed `Cargo.toml`: removed dhat, serde stack and 8 unused dev-dependencies plus the dhat features; `tracing-subscriber` moved to dev; `tracing` and `crypto-bigint` kept for planned work.
- Task `0006-trim-dependencies`; source: user request 2026-10-13: publication cleanup, remove ancient code and golden vectors, foo.rs module layout

## 2026-10-14 — remove dead code

- Moved dead code and stale artifacts (golden vectors, scripts, uncompiled encoders, `operations.rs`) into `tmp/`; removed commented-out code and `channels_mut`; fixed stale references.
- Task `0005-remove-dead-code`; source: user request 2026-10-13: publication cleanup, remove ancient code and golden vectors, foo.rs module layout

## 2026-10-13 — Sensor std-dev end-to-end (M1 CKKS green + M2 native demo)

- Roadmap M1 + M2 closed through tasks 0001, 0004, 0002, 0003.
- Decisions locked: sensor range contract `[-40, 60]` degC with std envelope `[1, 8]` degC, Chebyshev degree 7, 13 x 61-bit primes at N = 8192, 32 readings (one per 45 min) tiled over all slots, population std. Durable notes live in `specs/overview.md` gotchas.

## 2026-10-13 — sensor demo example

- Native N = 8192 sensor demo: 32 readings, std within 1.1e-2 degC of plaintext in ~25 s; scalars now fill all slots (`encode_constant`), fixing a slot-broadcast bug.
- Task `0003-sensor-demo-example`; source: user request 2026-10-12: revive project, CKKS green + native sensor demo with range normalization

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
