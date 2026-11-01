---
id: 0033-bootstrapped-recurrence-demo
status: done
scope: [src/ckks/algorithms/recurrence.rs, src/ckks/algorithms.rs, examples/recurrence_bootstrap_demo.rs, README.md]
created: 2026-11-01
---
# Task: bootstrapped recurrence demo

## Plan

- [x] Add a generic public-key-only recurrence plan and pre-generated level keys; consume one prime per step with internal scale alignment.
- [x] Implement 32-round comparisons: 18 primes with refresh, 18 without (stop safely), and 34 without (finish), at N=16/default and N=64.
- [x] Record levels, modulus bits, refresh events, timing and plaintext checkpoint errors. Validate ranges, scale/depth, immutability and error paths with seeded regressions.

## Context

For encrypted driver u and x0=u, compute x_next=0.5*x+0.125*x^2+0.25*sign*u, alternating sign from +1. Ideal states remain in [-0.625,0.875]; errors contract by at most0.75 per step. Keep scale44, nonzero key noise, explicit bootstrap and the 1e-3 error gate. No total-memory claim. Do not change bootstrap profiles or implement tasks0034/0035.

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features`
- [x] `cargo build --examples`
- [x] `cargo run --release --example recurrence_bootstrap_demo`
- [x] `cargo run --release --example recurrence_bootstrap_demo -- --degree 64`
- [x] `specdev check`
- [x] `typos`
- [x] `git diff --check`

## Findings

At both degrees the short chain finishes32 rounds with8 refreshes or stops before round17 without refresh. The34-prime baseline finishes32 rounds without refresh. Exact modulus bit lengths are1098 and2074. Demo worst checkpoint errors (including refreshes): N=16 about1.50e-5; N=64 about3.89e-5. Setup and diagnostics are excluded from evaluation time. Public evaluation uses nonzero-noise keys and receives no secret or RNG; all scale alignment is inside the library step.

## Review

The weighted square, linear and driver terms are aligned before a single rescale; exact dyadic coefficients preserve the one-prime schedule. Setup validates chain divisors and degree, and evaluation rejects wrong scales/contexts, missing keys and insufficient depth. Seeded tests cover both degrees and all three configurations, range invariants, deterministic evaluation and unchanged driver ciphertext. Release all-features tests passed. Tasks0034/0035 remain drafts with no implementation.

## Summary

Added a 32-round encrypted recurrence comparing a short bootstrapped chain, short-chain exhaustion and a longer leveled baseline.

## Log

- 2026-11-01 created: bootstrapped recurrence demo
- 2026-11-01 scope add src/ckks/algorithms/recurrence.rs — bounded recurrence and baseline comparisons
- 2026-11-01 scope add src/ckks/algorithms.rs — bounded recurrence and baseline comparisons
- 2026-11-01 scope add examples/recurrence_bootstrap_demo.rs — bounded recurrence and baseline comparisons
- 2026-11-01 scope add README.md — bounded recurrence and baseline comparisons
- 2026-11-01 draft approved: Implement the plan.
- 2026-11-01 advance draft → ready
- 2026-11-01 scope approved: src/ckks/algorithms/recurrence.rs, src/ckks/algorithms.rs, examples/recurrence_bootstrap_demo.rs, README.md
- 2026-11-01 advance ready → in-progress/implement
- 2026-11-01 baseline: examples/bootstrap_demo.rs@81389ac629b3, src/ckks/algorithms/bootstrap.rs@41579c204f90, src/ckks/keys.rs@d0b53f26bc7e, src/ckks/keys/conjugation_key.rs@eddfd1c737f7, src/ckks/keys/rotation_key.rs@4a8d2128faf5, src/ckks/ops.rs@d9424d9405fd, src/rings/backends/rns_ntt/poly.rs@5c302b70b49f, src/rings/errors.rs@fea138cd1c12, src/rings/traits.rs@8a54adc1cbc6
- 2026-11-01 advance in-progress/implement → in-progress/verify
- 2026-11-01 advance in-progress/verify → in-progress/review
- 2026-11-01 advance in-progress/review → approval
- 2026-11-01 approved: Ok it works fine, let's to bookkeeping and close the case.
- 2026-11-01 advance approval → done
