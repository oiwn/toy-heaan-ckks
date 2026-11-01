---
id: 0036-bootstrap-depth-optimization
status: done
scope: [src/ckks/algorithms/bootstrap.rs, src/ckks/algorithms/recurrence.rs, examples/recurrence_bootstrap_demo.rs, README.md, specs/research/bootstrap.md]
created: 2026-11-01
---
# Task: bootstrap depth optimization

## Plan

- [x] Fuse half-splitting and EvalMod interval normalization into forward transform diagonals; preserve semantics with stage references and complex regressions.
- [x] Recompute bootstrap schedules/key levels and demonstrate reduced depth without increasing the18-prime chain or weakening accuracy/noise.
- [x] Validate recurrence at level2 using its91-bit raw-product scale and122-bit modulus; refresh only when insufficient headroom remains.
- [x] Update demo counts, before/after evidence and documentation; run full checks at both toy degrees.

## Context

Operator requested fewer bootstraps. First optimize the current modulus budget. Target six output primes and three refreshes for32 recurrence rounds, subject to encrypted validation. The recurrence's weighted product has scale2^91, below the two-prime modulus near2^122; its final scale30 result can be lifted to44 at level1. General working-scale bootstrap multiplication retains its separate headroom policy. Do not disable noise, loosen1e-3 error or reduce the workload.

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

Same18-prime1098-bit modulus; bootstrap spends12 levels instead of14 and returns6 instead of4. Recurrence uses17 initial rounds and five between refreshes, reducing32-round refreshes8 ->3. Demo worst checkpoint error: N=16 about1.29e-5, N=64 about6.67e-5, including refreshes. The unbootstrapped short chain now stops before round18. No precision target, noise parameter or workload was relaxed. This is a depth reduction under the current modulus budget, not a claim of a global optimum.

## Review

Forward diagonals encode1/(2*radius), and real/imaginary extraction uses conjugation and exact polynomial multiplication by +/-i. EvalMod receives already normalized inputs and adds only the phase offset. Inverse plaintext/key levels move from5 to7. Independent coefficient references, scalar EvalMod and full complex regressions passed in release mode. Recurrence-only two-prime headroom is separately covered by boundary slots and both driver signs; bootstrap multiplication keeps its original three-level guard.

## Summary

Reduced the32-round toy recurrence from eight bootstraps to three using the same18-prime budget.

## Log

- 2026-11-01 created: bootstrap depth optimization
- 2026-11-01 scope add src/ckks/algorithms/bootstrap.rs — fuse normalization and use verified recurrence headroom
- 2026-11-01 scope add src/ckks/algorithms/recurrence.rs — fuse normalization and use verified recurrence headroom
- 2026-11-01 scope add examples/recurrence_bootstrap_demo.rs — fuse normalization and use verified recurrence headroom
- 2026-11-01 scope add README.md — fuse normalization and use verified recurrence headroom
- 2026-11-01 draft approved: Let’s do optimisations we need to minimize number of bootstraps
- 2026-11-01 advance draft → ready
- 2026-11-01 scope approved: src/ckks/algorithms/bootstrap.rs, src/ckks/algorithms/recurrence.rs, examples/recurrence_bootstrap_demo.rs, README.md
- 2026-11-01 advance ready → in-progress/implement
- 2026-11-01 baseline: examples/bootstrap_demo.rs@81389ac629b3, src/ckks/algorithms.rs@8898367d0c24, src/ckks/keys.rs@d0b53f26bc7e, src/ckks/keys/conjugation_key.rs@eddfd1c737f7, src/ckks/keys/rotation_key.rs@4a8d2128faf5, src/ckks/ops.rs@d9424d9405fd, src/rings/backends/rns_ntt/poly.rs@5c302b70b49f, src/rings/errors.rs@fea138cd1c12, src/rings/traits.rs@8a54adc1cbc6
- 2026-11-01 scope add specs/research/bootstrap.md — record fused schedule and measured depth improvement
- 2026-11-01 advance in-progress/implement → in-progress/verify
- 2026-11-01 advance in-progress/verify → in-progress/review
- 2026-11-01 advance in-progress/review → approval
- 2026-11-01 approved: Ok it works fine, let's to bookkeeping and close the case.
- 2026-11-01 advance approval → done
