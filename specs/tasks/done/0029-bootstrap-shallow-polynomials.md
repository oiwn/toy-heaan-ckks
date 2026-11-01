---
id: 0029-bootstrap-shallow-polynomials
status: done
scope: [src/ckks/algorithms/**, tests/**]
created: 2026-11-01
depends: [0026-bootstrap-automorphisms-keys]
---
# Task: bootstrap shallow polynomials

## Plan

- [x] Add balanced Chebyshev evaluation using cached even/odd identities and public evaluation keys; retain current Clenshaw API and behavior.
- [x] Align levels/scales for branch joins and coefficient sums. For degree31 budget at most five product layers plus affine-map and coefficient-weighting layers; report actual depth/headroom.
- [x] Compare against independent plaintext Clenshaw for constants, odd/even polynomials, degree31 and interval endpoints; test invalid coefficients/domains and insufficient levels.

## Context

Design: `specs/research/bootstrap.md`. Evaluation receives no secret key or RNG. The shallow schedule is a prerequisite for EvalMod, not a rewrite of all existing polynomial callers. Implemented; awaiting final operator review.

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features`
- [x] `cargo build --examples`
- [x] `specdev check`
- [x] `typos`

## Review

Validated with the full all-features suite, release bootstrap regressions (both degrees, eight seeds, six patterns, two refreshes), strict Clippy and example builds. Public bootstrap evaluation has no secret key or RNG. Supported profiles are fixed and validated at construction; general-purpose parameter combinations are not exposed.

## Summary

Implemented bootstrap shallow polynomials.

## Log

- 2026-11-01 created: bootstrap shallow polynomials
- 2026-11-01 scope add src/ckks/algorithms/chebyshev.rs — balanced Chebyshev evaluator
- 2026-11-01 scope add src/ckks/algorithms/sum_slots.rs — balanced Chebyshev evaluator
- 2026-11-01 scope add tests/integration_mul.rs — balanced Chebyshev evaluator
- 2026-11-01 set depends 0026-bootstrap-automorphisms-keys
- 2026-11-01 scope rm src/ckks/algorithms/chebyshev.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm src/ckks/algorithms/sum_slots.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm tests/integration_mul.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope add src/ckks/algorithms/** — declare complete module scope including new files
- 2026-11-01 scope add tests/** — declare complete module scope including new files
- 2026-11-01 draft approved: Implement the plan.
- 2026-11-01 advance draft → ready
- 2026-11-01 scope approved: src/ckks/algorithms/**, tests/**
- 2026-11-01 advance ready → in-progress/implement
- 2026-11-01 baseline: README.md@f33b4f4fc917, examples/bootstrap_demo.rs@38295abf01b4, src/ckks/algorithms.rs@6f2f1e43a81b, src/ckks/keys.rs@d0b53f26bc7e, src/ckks/keys/conjugation_key.rs@eddfd1c737f7, src/ckks/keys/rotation_key.rs@4a8d2128faf5, src/ckks/ops.rs@93d912d83dac, src/rings/backends/rns_ntt/poly.rs@5c302b70b49f, src/rings/errors.rs@fea138cd1c12, src/rings/traits.rs@8a54adc1cbc6
- 2026-11-01 advance in-progress/implement → in-progress/verify
- 2026-11-01 advance in-progress/verify → in-progress/review
- 2026-11-01 advance in-progress/review → approval
- 2026-11-01 approved: ok, looks good.
- 2026-11-01 advance approval → done
