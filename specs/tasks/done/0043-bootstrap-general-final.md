---
id: 0043-bootstrap-general-final
status: done
scope: [src/ckks/algorithms/bootstrap.rs, tests/bootstrap.rs, src/ckks/algorithms/bootstrap/*.rs]
created: 2026-11-02
source: Operator-approved bootstrap split plan 2026-11-02 (from 0037)
depends: [0042-multi-group-transforms]
---
# Task: bootstrap general final

## Context

Replace the reference-profile whitelist with validation derived from the chain and the config. The recurrence/examples migration and the final matrix moved to 0044 (operator split 2026-11-02).

Message ratio (operator request 2026-11-02): the quantity that matters is q0/message_scale, not the absolute scale. 2^44 is 61 - 17: 2^40 (ratio 2^21) missed the N=64 second-cycle gate, and 2^44 passed. See `specs/research/bootstrap.md`. Operator decisions: the message scale stays a power of two, `2^(round(log2 q0) - ratio)`; chains whose primes drift the working scale are rejected by a guard rather than fixed in the planner.

Rules: touch only scoped files (extend with `specdev task scope` first). If the same check fails twice, stop, record the failure and a hypothesis in Findings and ask the operator. Never weaken tolerances, disable noise or add `#[ignore]`.

## Plan

- [x] `BootstrapConfig` takes `message_ratio_bits` (log2(q0/message scale), like Lattigo `LogMessageRatio`; reference 17) instead of an absolute `message_scale`. The schedule derives the message scale from the chain's q0 and exposes it.
- [x] Replace the whitelist (degree, scales, polynomial settings, 18 primes, near-2^61) with derived checks: power-of-two degree, budgets `1..=log2(N/2)`, weight `1..=N`, degree >= 1, ratio >= 1 with message scale below q0, transform headroom, chain length.
- [x] Drift guard: EvalMod and double-angle scales stay within 1 bit of the working scale, otherwise the schedule names the level.
- [x] Tests: invalid configs (ratio, degree, weight, budgets, oversized working scale, alternating 60/61 chain, short chain); accepted configs build (N=32/128, uniform 60-bit chain, degree 63 with 3 double angles, weight N/4); encrypted smoke `general_configs` within 1e-3.

## Findings

- The shortest chain for the reference is 12 primes (11 levels plus the bottom), not 18. The old test asserted that 17 primes fail only because of the exact-18 rule; it now checks the real boundary at 12/11.
- Degree 30 and three double angles were rejected by the whitelist only; they are valid now and moved out of the invalid list.
- A 2^62 working scale on 61-bit primes passes headroom; the drift guard rejects it at the first squaring.
- Smoke errors (cycle0/cycle1): N=32 reference 1.7e-5/5.4e-5; N=16 60-bit chain 2.7e-5/5.1e-5; N=16 degree 63 with 3 double angles 9.5e-6/1.8e-5. Details in `specs/research/generic-bootstrap.md`.

## Review

- Scope: `bootstrap.rs` (validation, docs, private tests), `bootstrap/schedule.rs` (message ratio, drift guard, headroom) and `tests/bootstrap.rs` (smoke test). The planner, transforms, recurrence and examples are untouched.
- The reference demo and both recurrence demos reproduce 0042 bit for bit (message scale still 2^44).
- Approximation accuracy is not validated at construction; the module doc says so. Recurrence still has its own N=16/64 whitelist (0044).

## Summary

Bootstrap takes a message ratio and validates any config against the chain (levels, headroom, scale drift) instead of a fixed reference profile.

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

- 2026-11-02 created: bootstrap general final
- 2026-11-02 set source Operator-approved bootstrap split plan 2026-11-02 (from 0037)
- 2026-11-02 set depends 0042-multi-group-transforms
- 2026-11-02 scope add src/ckks/algorithms/bootstrap.rs — Remove profile restrictions, migrate recurrence/examples, final validation matrix
- 2026-11-02 scope add src/ckks/algorithms/recurrence.rs — Remove profile restrictions, migrate recurrence/examples, final validation matrix
- 2026-11-02 scope add examples/**/*.rs — Remove profile restrictions, migrate recurrence/examples, final validation matrix
- 2026-11-02 scope add tests/bootstrap.rs — Remove profile restrictions, migrate recurrence/examples, final validation matrix
- 2026-11-02 scope add src/ckks/algorithms/bootstrap/*.rs — Config message ratio, drift guard and headroom live in schedule.rs
- 2026-11-02 scope rm src/ckks/algorithms/recurrence.rs — Moved to 0044
- 2026-11-02 scope rm examples/**/*.rs — Moved to 0044
- 2026-11-02 draft approved: approved by operator
- 2026-11-02 advance draft → ready
- 2026-11-02 scope approved: src/ckks/algorithms/bootstrap.rs, tests/bootstrap.rs, src/ckks/algorithms/bootstrap/*.rs
- 2026-11-02 advance ready → in-progress/implement
- 2026-11-02 advance in-progress/implement → in-progress/verify
- 2026-11-02 advance in-progress/verify → in-progress/review
- 2026-11-02 advance in-progress/review → approval
- 2026-11-04 approved: approved by operator
- 2026-11-04 advance approval → done
