---
id: 0014-remove-commented-code
status: done
scope: [src/ckks/keys.rs, src/ckks/keys/common.rs, src/math/utils.rs, src/ckks/errors.rs, benches/basic.rs, README.md]
created: 2026-10-14
source: operator 2026-10-14: huge commented-out sections in keys/common.rs
---
# Task: remove commented code

## Plan

- [x] Move the fully commented `src/ckks/keys/common.rs` to `tmp/`; drop `pub mod common;`.
- [x] Move the fully commented `benches/basic.rs` to `tmp/`; drop the `cargo bench --bench basic` README line.
- [x] Delete the commented `RnsBasis` tests in `src/math/utils.rs` and the two commented variants in `src/ckks/errors.rs`.

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features`
- [x] `cargo build --examples`
- [x] `cargo bench --no-run`
- [x] `typos`
- [x] `gitleaks dir .`

## Summary

Removed dead commented-out code from the old `RnsPolyRing` API: `keys/common.rs` and `benches/basic.rs` (both moved to `tmp/`), stale tests in `math/utils.rs` and two `CkksError` variants.

## Review

- `basic` had no `[[bench]]` entry, so Cargo.toml stays untouched; `primes` is the only bench.
- No `/*` block comments remain in `src`, `examples`, `tests` or `benches`. Test counts unchanged (128/8/4/1).

## Log

- 2026-10-14 created: remove commented code
- 2026-10-14 set source operator 2026-10-14: huge commented-out sections in keys/common.rs
- 2026-10-14 scope add src/ckks/keys.rs — remove commented-out code
- 2026-10-14 scope add src/ckks/keys/common.rs — remove commented-out code
- 2026-10-14 scope add src/math/utils.rs — remove commented-out code
- 2026-10-14 scope add src/ckks/errors.rs — remove commented-out code
- 2026-10-14 scope add benches/basic.rs — remove commented-out code
- 2026-10-14 scope add README.md — remove commented-out code
- 2026-10-14 advance draft → ready
- 2026-10-14 scope approved: src/ckks/keys.rs, src/ckks/keys/common.rs, src/math/utils.rs, src/ckks/errors.rs, benches/basic.rs, README.md
- 2026-10-14 advance ready → in-progress/implement
- 2026-10-14 baseline: examples/encrypt_add.rs@818bcae6e49d, examples/encrypt_mul.rs@09945cc61939, examples/horner_chain.rs@d5855b6385e8, examples/rotation_demo.rs@041ec7abf627, examples/rotation_stress.rs@3a088df82cec, examples/sensor_std_dev/contract.rs@8586cd0a931a, examples/sensor_std_dev/main.rs@0ff6dcc8f9a8, examples/sensor_std_dev/std_dev.rs@860ae61969bb, examples/sensor_std_plain.rs@56a00ea5f147, examples/std_dev_8.rs@0740eaaf8827, src/ckks.rs@c17a20cd6292, src/ckks/algorithms.rs@e531593da936, src/ckks/algorithms/chebyshev.rs@09b0f23e0e26, src/ckks/algorithms/sum_slots.rs@1cced47f11e0, src/ckks/engine.rs@786b33d429c9, src/ckks/evaluator.rs@bcda1084cfae, src/ckks/keys/secret_key.rs@bca0326fa4e1, src/ckks/ops.rs@7735139c61b8, src/ckks/pipeline.rs@deleted, src/ckks/types.rs@9e4a55df6049, src/lib.rs@9c6c2bc869bb, src/rings.rs@1c1f26b5f8ab, src/rings/backends/rns_ntt/basis.rs@d575a3994242, src/rings/backends/rns_ntt/errors.rs@a541997fb455, src/rings/backends/rns_ntt/poly.rs@58063d5c53a7, src/rings/traits.rs@2a154fdb666a, tests/integration_mul.rs@d161404eb898
- 2026-10-14 advance in-progress/implement → in-progress/verify
- 2026-10-14 advance in-progress/verify → in-progress/review
- 2026-10-14 advance in-progress/review → approval
- 2026-10-15 approved: approved by operator
- 2026-10-15 advance approval → done
