---
id: 0027-bootstrap-mod-raise
status: done
scope: [src/rings/**, tests/**]
created: 2026-11-01
---
# Task: bootstrap mod raise

## Plan

- [x] Add fallible `PolyRing::mod_raise` into a supplied compatible target context. V1 requires bottom-level input and a higher context in its original chain.
- [x] Lift centered signed coefficients into each target prime, return NTT form and preserve target context identity; reject incompatible/same/lower targets.
- [x] Test exact representatives near zero and half-modulus boundaries, source residues preserved, NTT/coefficient agreement and all invalid-input paths.

## Context

Design: `specs/research/bootstrap.md`. ModRaise changes the decryption representative by modulus multiples; it does not refresh a ciphertext on its own. Use exact integers, not f64 reconstruction. Implemented; awaiting final operator review.

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

Implemented bootstrap mod raise.

## Log

- 2026-11-01 created: bootstrap mod raise
- 2026-11-01 scope add src/rings/backends — centered ModRaise backend contract
- 2026-11-01 scope add src/rings/backends.rs — centered ModRaise backend contract
- 2026-11-01 scope add src/rings/errors.rs — centered ModRaise backend contract
- 2026-11-01 scope add src/rings/traits.rs — centered ModRaise backend contract
- 2026-11-01 scope add tests/integration_mul.rs — centered ModRaise backend contract
- 2026-11-01 scope rm src/rings/backends — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm src/rings/backends.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm src/rings/errors.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm src/rings/traits.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm tests/integration_mul.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope add src/rings/** — declare complete module scope including new files
- 2026-11-01 scope add tests/** — declare complete module scope including new files
- 2026-11-01 draft approved: Implement the plan.
- 2026-11-01 advance draft → ready
- 2026-11-01 scope approved: src/rings/**, tests/**
- 2026-11-01 advance ready → in-progress/implement
- 2026-11-01 baseline: README.md@f33b4f4fc917, examples/bootstrap_demo.rs@38295abf01b4, src/ckks/algorithms.rs@6f2f1e43a81b, src/ckks/algorithms/bootstrap.rs@0ad1e98c3a13, src/ckks/keys.rs@d0b53f26bc7e, src/ckks/keys/conjugation_key.rs@eddfd1c737f7, src/ckks/keys/rotation_key.rs@4a8d2128faf5, src/ckks/ops.rs@93d912d83dac
- 2026-11-01 advance in-progress/implement → in-progress/verify
- 2026-11-01 advance in-progress/verify → in-progress/review
- 2026-11-01 advance in-progress/review → approval
- 2026-11-01 approved: ok, looks good.
- 2026-11-01 advance approval → done
