---
id: 0030-bootstrap-eval-mod
status: done
scope: [src/ckks/algorithms/**, src/ckks/algorithms.rs, tests/**]
created: 2026-11-01
depends: [0029-bootstrap-shallow-polynomials]
---
# Task: bootstrap eval mod

## Plan

- [x] Implement continuous cosine Chebyshev approximation, quarter-period shift, four double angles and 1/(2*pi) amplitude restoration; parameterize the domain and degree.
- [x] Use the shallow evaluator with explicit scale/level validation. Start with degree31, wrap bounds5/17 and message ratio near 2^-17; use actual modulus corrections.
- [x] Reproduce the plaintext sweep and compare encrypted real branches at every supported wrap, signs, zero, small-message boundaries and fixed seeds. Measure both approximation and key-switch error after message restoration.
- [x] Confirm the profile against the allocated error budget and depth; document hidden-domain preconditions and invalid observable configuration errors.

## Context

Design and measured plaintext candidate: `specs/research/bootstrap.md`; experiment: `specs/research/bootstrap_experiment.py`. Do not assume a sampled f64 result proves encrypted accuracy. Arcsine correction is deferred. Implemented; awaiting final operator review.

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

Implemented bootstrap eval mod.

## Log

- 2026-11-01 created: bootstrap eval mod
- 2026-11-01 scope add src/ckks/algorithms/chebyshev.rs — bounded EvalMod
- 2026-11-01 scope add src/ckks/algorithms/sum_slots.rs — bounded EvalMod
- 2026-11-01 scope add tests/integration_mul.rs — bounded EvalMod
- 2026-11-01 set depends 0029-bootstrap-shallow-polynomials
- 2026-11-01 scope rm src/ckks/algorithms/chebyshev.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm src/ckks/algorithms/sum_slots.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm tests/integration_mul.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope add src/ckks/algorithms/** — declare complete module scope including new files
- 2026-11-01 scope add src/ckks/algorithms.rs — declare complete module scope including new files
- 2026-11-01 scope add tests/** — declare complete module scope including new files
- 2026-11-01 draft approved: Implement the plan.
- 2026-11-01 advance draft → ready
- 2026-11-01 scope approved: src/ckks/algorithms/**, src/ckks/algorithms.rs, tests/**
- 2026-11-01 advance ready → in-progress/implement
- 2026-11-01 baseline: README.md@f33b4f4fc917, examples/bootstrap_demo.rs@38295abf01b4, src/ckks/keys.rs@d0b53f26bc7e, src/ckks/keys/conjugation_key.rs@eddfd1c737f7, src/ckks/keys/rotation_key.rs@4a8d2128faf5, src/ckks/ops.rs@93d912d83dac, src/rings/backends/rns_ntt/poly.rs@5c302b70b49f, src/rings/errors.rs@fea138cd1c12, src/rings/traits.rs@8a54adc1cbc6
- 2026-11-01 advance in-progress/implement → in-progress/verify
- 2026-11-01 advance in-progress/verify → in-progress/review
- 2026-11-01 advance in-progress/review → approval
- 2026-11-01 approved: ok, looks good.
- 2026-11-01 advance approval → done
