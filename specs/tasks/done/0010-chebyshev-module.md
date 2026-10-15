---
id: 0010-chebyshev-module
status: done
scope: [src/ckks/**, src/ckks.rs, examples/**, README.md]
created: 2026-10-13
source: user request 2026-10-13: pipeline.rs is a std-dev grab bag; decompose into composable ops like fhe.rs / tfhe-rs
depends: [0009-ops-composition]
---
# Task: chebyshev module

## Plan

- [x] `src/ckks/algorithms.rs` + `algorithms/sum_slots.rs`: `sum_slots(ct, slots, &RotationKeyMap)` (rotate + add) with `required_rotation_offsets` (moved from `evaluator.rs`); `LeveledEvaluator::sum_slots` stays the cached-key wrapper.
- [x] `src/ckks/algorithms/chebyshev.rs`: `eval_chebyshev_encrypted` (seeded by the new `Ciphertext::zero_like` op, no ring access), `fit_chebyshev`, `eval_chebyshev_clenshaw`, the domain maps, and the Clenshaw tests.
- [x] Drop the public `Ct` alias (private aliases where needed); move `src/ckks/pipeline.rs` into `tmp/`.
- [x] Update examples (std_dev_8, sensor_std_dev, sensor_std_plain), README/AGENTS layout lines and `specs/overview.md`.

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

`pipeline.rs` retired: `sum_slots` and the Chebyshev code (encrypted Clenshaw, fit, plaintext reference) now live in `ckks/algorithms/{sum_slots, chebyshev}.rs`, built only on ops, level methods and the evaluator; new `Ciphertext::zero_like` op.

## Review

- Plan said seed Clenshaw with `&two_t - &two_t`; clippy rejects that (`eq_op`), so the trivial `(0, 0)` became a small `Ciphertext::zero_like` op in `ops.rs` (tested there). Algorithms still never call the ring.
- `Ct` is no longer public: `evaluator.rs`, `chebyshev.rs` and the sensor example each keep a private alias (0012 replaces them with the generic type).
- `pipeline.rs` moved whole to `tmp/src/ckks/pipeline.rs` (pre-move content, for reference).
- Numbers identical: `sensor_std_dev` stage table and errors match the baseline; `std_dev_8` 8.25e-12; `sensor_std_plain` and `sqrt_chebyshev_plain` run. Tests: 128 lib (+1 offsets test), 8 integration, 4 example, 1 doc.

## Context

Algorithms use only `ops`, `Ciphertext` methods and the evaluator, never the ring directly (`specs/overview.md` Target architecture).

## Log

- 2026-10-13 created: chebyshev module
- 2026-10-13 set source user request 2026-10-13: pipeline.rs is a std-dev grab bag; decompose into composable ops like fhe.rs / tfhe-rs
- 2026-10-13 set depends 0009-ops-composition
- 2026-10-13 scope add src/ckks/** — chebyshev module, retire pipeline.rs
- 2026-10-13 scope add src/ckks.rs — chebyshev module, retire pipeline.rs
- 2026-10-13 scope add examples/** — chebyshev module, retire pipeline.rs
- 2026-10-14 scope add README.md — layout line names pipeline
- 2026-10-14 advance draft → ready
- 2026-10-14 scope approved: src/ckks/**, src/ckks.rs, examples/**, README.md
- 2026-10-14 advance ready → in-progress/implement
- 2026-10-14 baseline: src/lib.rs@9c6c2bc869bb, src/rings.rs@1c1f26b5f8ab, src/rings/backends/rns_ntt/basis.rs@d575a3994242, src/rings/backends/rns_ntt/errors.rs@a541997fb455, src/rings/backends/rns_ntt/poly.rs@58063d5c53a7, src/rings/traits.rs@2a154fdb666a, tests/integration_mul.rs@d161404eb898
- 2026-10-14 advance in-progress/implement → in-progress/verify
- 2026-10-14 advance in-progress/verify → in-progress/review
- 2026-10-14 advance in-progress/review → approval
- 2026-10-14 approved: approved by operator
- 2026-10-14 advance approval → done
