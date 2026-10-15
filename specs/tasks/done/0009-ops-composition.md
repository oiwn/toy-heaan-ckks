---
id: 0009-ops-composition
status: done
scope: [src/ckks/**, src/ckks.rs, examples/**, tests/**]
created: 2026-10-13
source: user request 2026-10-13: pipeline.rs is a std-dev grab bag; decompose into composable ops like fhe.rs / tfhe-rs
depends: [0008-ciphertext-levels]
---
# Task: ops composition

## Plan

- [x] `src/ckks/ops.rs`: `std::ops` Add/Sub/Neg for `&ct + &ct` and `&ct + &pt` (fhe.rs style), replacing `add_ciphertexts`/`sub_ciphertexts`/`neg_ciphertext`.
- [x] Same file, `Ciphertext` methods (decided 2026-10-14): `mul(&other, &rlk)`, `mul_plain(&pt)`, `rotate(&rotk)` moved out of `CkksEngine`; scalars `mul_const(value, scale_bits)` (no rescale) and `add_const(value)`. Rescale stays explicit.
- [x] Remove `mul_plain_scalar`/`mul_plain_slots`/`mul_plain_exact`, `trivial_constant_ct`, `add_constant`/`sub_constant`; Clenshaw rewritten on the new ops.
- [x] Move `LeveledEvaluator` (per-level key cache) and the rotation-key helpers to `src/ckks/evaluator.rs`.
- [x] Update examples and tests; `ops` unit tests; behavior unchanged (demo std error stays ~1e-2 degC); `specs/overview.md` updated.

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features`
- [x] `cargo build --examples`
- [x] `cargo bench --no-run`
- [x] `cargo run --release --example sensor_std_dev`
- [x] `typos`
- [x] `gitleaks dir .`

## Summary

Homomorphic ops moved to `ckks/ops.rs`: `std::ops` Add/Sub/Neg plus `Ciphertext::{mul, mul_plain, mul_const, add_const, rotate}` with explicit rescale; `LeveledEvaluator` and rotation-key helpers moved to `ckks/evaluator.rs`; the `pipeline.rs` scalar helpers and the engine's op functions are gone.

## Review

- `CkksEngine` keeps keygen, encrypt, decrypt, the legacy `mul_ciphertexts` path (0012) and the gadget key switch, now `pub(crate)` so `mul`/`rotate` can call it.
- The `mul_plaintext_and_rescale` logp "lift back" step is gone; with SCALE-bit primes it never fired, and every example prints the same numbers as before (`sensor_std_dev` stages and errors 5.21e-4 / 1.08e-2, `std_dev_8` 8.25e-12, `horner_chain` 7.890e-10).
- Clenshaw seeds use a private `zero_like` trivial ciphertext plus `add_const`; `Ct` alias stays in `pipeline.rs` until 0010.
- Tests: 127 lib (+3 `ops` tests), 8 integration, 4 example, 1 doc; `cargo doc` warnings stay at 4.

## Context

Ops stay concrete on `RnsPoly`; 0012 makes them generic. Keygen, encrypt and decrypt stay on `CkksEngine` for now. `sum_slots` moves in 0010 (algorithms).

## Log

- 2026-10-13 created: ops composition
- 2026-10-13 set source user request 2026-10-13: pipeline.rs is a std-dev grab bag; decompose into composable ops like fhe.rs / tfhe-rs
- 2026-10-13 set depends 0008-ciphertext-levels
- 2026-10-13 scope add src/ckks/** — ops module, evaluator move, update callers
- 2026-10-13 scope add src/ckks.rs — ops module, evaluator move, update callers
- 2026-10-13 scope add examples/** — ops module, evaluator move, update callers
- 2026-10-13 scope add tests/** — ops module, evaluator move, update callers
- 2026-10-14 advance draft → ready
- 2026-10-14 scope approved: src/ckks/**, src/ckks.rs, examples/**, tests/**
- 2026-10-14 advance ready → in-progress/implement
- 2026-10-14 baseline: README.md@f032018d21cb, src/lib.rs@9c6c2bc869bb, src/rings.rs@1c1f26b5f8ab, src/rings/backends/rns_ntt/basis.rs@d575a3994242, src/rings/backends/rns_ntt/errors.rs@a541997fb455, src/rings/backends/rns_ntt/poly.rs@58063d5c53a7, src/rings/traits.rs@2a154fdb666a
- 2026-10-14 advance in-progress/implement → in-progress/verify
- 2026-10-14 advance in-progress/verify → in-progress/review
- 2026-10-14 advance in-progress/review → approval
- 2026-10-14 approved: approved by operator
- 2026-10-14 advance approval → done
