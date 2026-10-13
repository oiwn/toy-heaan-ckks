---
id: 0009-ops-composition
status: draft
scope: [src/ckks/**, src/ckks.rs, examples/**, tests/**]
created: 2026-10-13
source: user request 2026-10-13: pipeline.rs is a std-dev grab bag; decompose into composable ops like fhe.rs / tfhe-rs
depends: [0008-ciphertext-levels]
---
# Task: ops composition

## Plan

- [ ] `src/ckks/ops.rs`: `std::ops` Add/Sub/Neg for ct+ct and ct+pt (fhe.rs style), replacing `add_ciphertexts`/`sub_ciphertexts`/`neg_ciphertext`.
- [ ] Same file: `mul(ct, ct, relin_key)`, `mul_plain(ct, pt)`, `rotate(ct, k, rot_key)` moved out of `CkksEngine`; `sum_slots(ct, window, keys)` = rotate + add.
- [ ] Remove `mul_plain_scalar`/`mul_plain_slots`/`mul_plain_exact`, `trivial_constant_ct`, `add_constant`/`sub_constant`: callers write `mul_plain(ct, &encode_constant(..))` + `rescale`, or `ct + &pt`.
- [ ] Move `LeveledEvaluator` (per-level key cache) and the rotation-key helpers to `src/ckks/evaluator.rs`.
- [ ] Update examples (sensor_std_dev, std_dev_8, rotation_*, encrypt_*) and tests; behavior unchanged (demo std error stays ~1e-2 degC).

## Acceptance

- [ ] `cargo fmt --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features`
- [ ] `cargo build --examples`
- [ ] `cargo run --release --example sensor_std_dev`
- [ ] `typos`
- [ ] `gitleaks dir .`

## Context

On hold: revised after the backend design in `specs/ctx.md` is agreed. Keygen, encrypt and decrypt stay on `CkksEngine` here; whether the engine stays generic over `P` is an audit question.

## Log

- 2026-10-13 created: ops composition
- 2026-10-13 set source user request 2026-10-13: pipeline.rs is a std-dev grab bag; decompose into composable ops like fhe.rs / tfhe-rs
- 2026-10-13 set depends 0008-ciphertext-levels
- 2026-10-13 scope add src/ckks/** — ops module, evaluator move, update callers
- 2026-10-13 scope add src/ckks.rs — ops module, evaluator move, update callers
- 2026-10-13 scope add examples/** — ops module, evaluator move, update callers
- 2026-10-13 scope add tests/** — ops module, evaluator move, update callers
