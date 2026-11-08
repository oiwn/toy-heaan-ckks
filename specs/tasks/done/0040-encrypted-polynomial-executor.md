---
id: 0040-encrypted-polynomial-executor
status: done
scope: [examples/**/*.rs, src/ckks/*.rs, src/ckks/algorithms/**/*.rs, tests/*.rs]
created: 2026-11-02
source: Operator-approved bootstrap split plan 2026-11-02 (from 0037)
depends: [0039-scale-metadata-migration]
---
# Task: encrypted polynomial executor

## Context

Phase two of the scale migration. With actual-prime rescale, values reaching an add by different paths differ in scale by about 2^-41, above the 64*eps tolerance. So the flip, the executor and every cross-path caller land together. Numerical rules are in `specs/research/generic-bootstrap.md`.

Rules: touch only scoped files (extend with `specdev task scope` first). If the same check fails twice, stop, record the failure and a hypothesis in Findings and ask the operator. Never weaken tolerances, disable noise or add `#[ignore]`.

## Plan

Operator decision 2026-11-02: **exact tracking.** Scales are never snapped. Bootstrap and recurrence accept non-power-of-two scales and renormalize only by exact `2^k` lifts. Tests compare `log2(scale)` within a tolerance instead of `==`.

- [x] Flip `Ciphertext::rescale` to `scale / rescale_divisor()` (same `u64 -> f64` conversion as the planner). Run the tests and record the breakage list in Findings before fixing anything.
- [x] `ops.rs`: `mul_const_to_scale(value, target)` multiplies by `value` encoded at `target*q/scale`, then rescales (the target-scale rule).
- [x] `bootstrap/executor.rs`: `ChebyshevPlan::execute(input, multiply)` interprets the plan graph. `multiply` is a relinearize-only callback. After every node, actual level/scale must match the `PlannedValue`; a mismatch is an error naming the node.
- [x] Bootstrap `eval_mod` runs on the executor; double angles stay explicit; `term`/`constant` go away. The inverse transform moves to `level - 10`, so output is level 7.
  - `prepare_input` accepts any scale `s` with `log2 s >= 1`, lifts by `2^round(log2(2^44/s))` (rejects a negative k), and the output scale is `tracked * s'/2^44`.
- [x] `eval_chebyshev_encrypted` returns `CkksResult`. It uses an affine map via `mul_const_to_scale`, then the executor, with `LeveledEvaluator::relin_key` for lazy keys.
- [x] Recurrence accepts `log2` scales in `[43, 45)`, encodes the driver at its own scale, and lifts by `2^k` to within half a bit of `2^44`. That is rounding, not the planned flooring into `(2^43, 2^44]`; see Findings. Examples and integration tests encode cross-path operands at the partner's scale.
- [x] Tests: encrypted vs plaintext Clenshaw for degrees 0/1/2/3/15/31/63 on a mixed 60/61-bit chain; node-metadata mismatch errors; the `mul_const_to_scale` target.
- [x] Bootstrap still meets the 1e-3 gate on 8 seeds x 2 cycles at N=16/64. Record the new demo baseline in `specs/research/generic-bootstrap.md`.

## Findings

Breakage after the rescale flip, before any fix:
- lib: bootstrap `private_stages_16/64` (`term` sub) and `recurrence_comparisons_16/64` (non-2^44 state).
- `tests/bootstrap.rs`: `prepare_input` rejects the non-power-of-two squared input.
- `tests/integration_mul.rs`: `multiply_then_add` and the matched-levels add.
- Examples: `sensor_std_dev` (Clenshaw sub), `horner_chain` (add of beta). `std_dev_8` passed only by luck of a prime near 2^61.

During implementation:
- **Planner drift gap (open).** With 2^61 inputs and 60-bit primes at the top of the chain, each squaring doubles the scale mismatch. At degree 63 the quotient leaves got plaintext scales of 2-8, and quantization destroyed the result; `compile` doesn't catch this. The test puts 61-bit primes under the power tree and keeps one 60-bit divisor. Details are in `specs/research/generic-bootstrap.md`.
- **Recurrence lift rounds.** At N=64 the rescaled result is 2^30 * (1 + about 2^-58), so flooring lifted it only to about 2^43.0. Rounding keeps it within half a bit of 2^44.

## Review

- Scope stays in the four globs; `src/ckks/*.rs` changes are `types.rs` (rescale), `ops.rs` (`mul_const_to_scale`) and `evaluator.rs` (`relin_key`).
- Planned and actual node scales agree because rescale and the planner share the divisor conversion. The executor checks every node and uses fallible encoding.
- Bootstrap output is at level 7, one level saved (`INVERSE_TRANSFORM_DEPTH = 10`). Input scales in [2, 2^44.5) are accepted, and the output scale is exact metadata.
- Accuracy is comparable to the 0039 baseline (worst 1.57e-4, was 2.08e-4; numbers in the research doc). Tests: 150 lib (was 147), all integration/example/doc tests pass, and all examples run.
- Caller-visible changes:
  - `eval_chebyshev_encrypted` returns `CkksResult` and uses `1 + depth` levels.
  - Fresh operands added to rescaled ciphertexts are encoded at their scale.

## Summary

Rescale divides by the actual prime; an encrypted `ChebyshevPlan` executor now drives bootstrap EvalMod and the Chebyshev evaluator, with exact scale tracking through bootstrap, recurrence and the examples.

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features`
- [x] `cargo build --examples`
- [x] `cargo run --release --example bootstrap_demo`
- [x] `cargo run --release --example recurrence_bootstrap_demo`
- [x] `typos`
- [x] `specdev check`
- [x] `git diff --check`

## Log

- 2026-11-02 created: encrypted polynomial executor
- 2026-11-02 set source Operator-approved bootstrap split plan 2026-11-02 (from 0037)
- 2026-11-02 set depends 0039-scale-metadata-migration
- 2026-11-02 scope add src/ckks/algorithms/bootstrap/polynomial.rs — Encrypted graph executor, actual-divisor rescale, migrate cross-path callers
- 2026-11-02 scope add src/ckks/algorithms/bootstrap.rs — Encrypted graph executor, actual-divisor rescale, migrate cross-path callers
- 2026-11-02 scope add src/ckks/algorithms/chebyshev.rs — Encrypted graph executor, actual-divisor rescale, migrate cross-path callers
- 2026-11-02 scope add src/ckks/types.rs — Encrypted graph executor, actual-divisor rescale, migrate cross-path callers
- 2026-11-02 scope add src/ckks/errors.rs — Encrypted graph executor, actual-divisor rescale, migrate cross-path callers
- 2026-11-02 scope add examples/**/*.rs — Encrypted graph executor, actual-divisor rescale, migrate cross-path callers
- 2026-11-02 scope rm src/ckks/algorithms/bootstrap/polynomial.rs — covered by broader glob
- 2026-11-02 scope rm src/ckks/algorithms/bootstrap.rs — covered by broader glob
- 2026-11-02 scope rm src/ckks/algorithms/chebyshev.rs — covered by broader glob
- 2026-11-02 scope rm src/ckks/types.rs — covered by broader glob
- 2026-11-02 scope rm src/ckks/errors.rs — covered by broader glob
- 2026-11-02 scope add src/ckks/*.rs — Actual-divisor rescale touches ops/evaluator/recurrence/executor and integration tests
- 2026-11-02 scope add src/ckks/algorithms/**/*.rs — Actual-divisor rescale touches ops/evaluator/recurrence/executor and integration tests
- 2026-11-02 scope add tests/*.rs — Actual-divisor rescale touches ops/evaluator/recurrence/executor and integration tests
- 2026-11-02 draft approved
- 2026-11-02 draft approved: approved.
- 2026-11-02 advance draft → ready
- 2026-11-02 scope approved: examples/**/*.rs, src/ckks/*.rs, src/ckks/algorithms/**/*.rs, tests/*.rs
- 2026-11-02 advance ready → in-progress/implement
- 2026-11-02 baseline: src/ckks.rs@2615ecc20a4b, src/ckks/encoding/ckks_encoder.rs@8f3e4ed1f6dd
- 2026-11-02 advance in-progress/implement → in-progress/verify
- 2026-11-02 advance in-progress/verify → in-progress/review
- 2026-11-02 advance in-progress/review → approval
- 2026-11-02 approved: approved by operator
- 2026-11-02 advance approval → done
