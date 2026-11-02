---
id: 0038-fallible-ciphertext-arithmetic
status: done
scope: [src/ckks/ops.rs, src/ckks/errors.rs, src/ckks/algorithms/*.rs, examples/**/*.rs, tests/*.rs]
created: 2026-11-02
source: Operator-approved bootstrap split plan 2026-11-02 (from 0037)
---
# Task: fallible ciphertext arithmetic

## Context

First of the bootstrap completion chain 0038-0043 (split out of 0037). Operator decision: every ciphertext add/sub is fallible; the `std::ops` Add/Sub impls for `&Ciphertext` (ct+ct, ct+pt, ct-ct, ct-pt) are removed. `Neg` stays. `logp` is unchanged in this task; numeric results must be identical.

Rules: touch only scoped files (extend with `specdev task scope` first). If the same check fails twice, stop, record the failure and a hypothesis in Findings and ask the operator. Never weaken tolerances, disable noise or add `#[ignore]`.

## Plan

- [x] `src/ckks/errors.rs`: add `LevelMismatch { expected, actual }`. Reuse the existing unused `ScaleMismatch { expected, actual }` (logp as f64 for now) and `RingError::IncompatibleContext` (as `conjugate` does) instead of new variants.
- [x] `src/ckks/ops.rs`: replace the Add/Sub impls with `Ciphertext::{add, sub, add_plain, sub_plain} -> CkksResult<Self>`. Check level, then context (`same_context`), then scale (`logp` equality for now). A level mismatch is an error, never an automatic mod switch. `add_const` stays infallible through a private unchecked helper. Update module docs.
- [x] Migrate library callers: bootstrap.rs, recurrence.rs, sum_slots.rs, chebyshev.rs. Use `?` where the function already returns `CkksResult`; `sum_slots` and `eval_chebyshev_encrypted` keep their documented panic contracts (`expect`).
- [x] Migrate examples (rotation_demo, std_dev_8, encrypt_add, horner_chain, sensor_std_dev/std_dev) and tests/integration_mul.rs.
- [x] Unit tests for each error variant on add/sub/add_plain/sub_plain.
- [x] Confirm `cargo run --release --example bootstrap_demo` max error is identical to the baseline in `specs/research/generic-bootstrap.md`.

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features`
- [x] `cargo build --examples`
- [x] `cargo run --release --example bootstrap_demo`
- [x] `typos`
- [x] `specdev check`
- [x] `git diff --check`

## Review

- `std::ops` Add/Sub for `&Ciphertext` removed; `Neg` kept. `add`/`sub`/`add_plain`/`sub_plain` return `CkksResult` and share one `check_compatible` (level -> context -> logp). Previously a level/basis mismatch only hit a `debug_assert` in `RnsPoly::add_assign`.
- `add_const` stays infallible via private `add_plain_unchecked`. `sum_slots` and `eval_chebyshev_encrypted` keep their panic contracts with explicit `expect` messages.
- Errors: new `LevelMismatch`; reused the formerly unused `ScaleMismatch` and `RingError::IncompatibleContext`.
- Evidence: 144 library tests (3 new error tests covering all four methods), 10 integration, 4 example, 1 doctest pass. Bootstrap demo max errors equal the baseline digit for digit. No ciphertext operator uses remain in src/examples/tests.
- Within scope; no tolerances changed.

## Summary

Ciphertext add/sub are fallible methods with level, context and scale checks; operator impls removed, results unchanged.

## Log

- 2026-11-02 created: fallible ciphertext arithmetic
- 2026-11-02 set source Operator-approved bootstrap split plan 2026-11-02 (from 0037)
- 2026-11-02 scope add src/ckks/ops.rs — Fallible add/sub replace std::ops impls; migrate all callers
- 2026-11-02 scope add src/ckks/errors.rs — Fallible add/sub replace std::ops impls; migrate all callers
- 2026-11-02 scope add src/ckks/algorithms/*.rs — Fallible add/sub replace std::ops impls; migrate all callers
- 2026-11-02 scope add examples/**/*.rs — Fallible add/sub replace std::ops impls; migrate all callers
- 2026-11-02 scope add tests/*.rs — Fallible add/sub replace std::ops impls; migrate all callers
- 2026-11-02 draft approved: approve
- 2026-11-02 advance draft → ready
- 2026-11-02 scope approved: src/ckks/ops.rs, src/ckks/errors.rs, src/ckks/algorithms/*.rs, examples/**/*.rs, tests/*.rs
- 2026-11-02 advance ready → in-progress/implement
- 2026-11-02 baseline: src/ckks.rs@2615ecc20a4b, src/ckks/algorithms/bootstrap/polynomial.rs@c822cc7baf28, src/ckks/algorithms/bootstrap/transforms.rs@022348ef8279, src/ckks/scale.rs@8fac076e0dc8
- 2026-11-02 advance in-progress/implement → in-progress/verify
- 2026-11-02 advance in-progress/verify → in-progress/review
- 2026-11-02 advance in-progress/review → approval
- 2026-11-02 approved: approved by operator
- 2026-11-02 advance approval → done
