---
id: 0028-bootstrap-dense-transforms
status: done
scope: [src/ckks/algorithms/**, src/ckks/algorithms.rs, tests/**]
created: 2026-11-01
depends: [0026-bootstrap-automorphisms-keys]
---
# Task: bootstrap dense transforms

## Plan

- [x] Derive dense forward/inverse matrices from this encoder's slot ordering; define independent plaintext matrix references and verify normalization.
- [x] Encode diagonals, evaluate rotation/plaintext products, accumulate then rescale once per transform where the schedule allows; split/recombine complex halves with explicit conjugation.
- [x] Fold modulus-ratio and final message-scale restoration into the matrices; certify actual-prime scale corrections and headroom against the research ledger before fixing the profile.
- [x] Test asymmetric full complex packing at N=16/64, basis vectors, zero/signs, independent matrix agreement, round trips, transform norm and consumed levels.

## Context

Use dense diagonals first, not BSGS or FFT factorization. Design and candidate 18-prime ledger: `specs/research/bootstrap.md`. Public evaluation uses only precomputed plaintexts/evaluation keys. Implemented; failure to certify the integer-logp schedule requires a reviewed scale-design follow-up, not silent metadata changes.

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

Implemented bootstrap dense transforms.

## Log

- 2026-11-01 created: bootstrap dense transforms
- 2026-11-01 scope add src/ckks/algorithms/chebyshev.rs — dense bootstrap transforms
- 2026-11-01 scope add src/ckks/algorithms/sum_slots.rs — dense bootstrap transforms
- 2026-11-01 scope add tests/integration_mul.rs — dense bootstrap transforms
- 2026-11-01 set depends 0026-bootstrap-automorphisms-keys
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
