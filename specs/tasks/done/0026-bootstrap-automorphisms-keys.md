---
id: 0026-bootstrap-automorphisms-keys
status: done
scope: [src/rings/traits.rs, src/ckks/keys/**, src/ckks/keys.rs, src/ckks/evaluator.rs, src/ckks/ops.rs, tests/**]
created: 2026-11-01
---
# Task: bootstrap automorphisms keys

## Plan

- [x] Implement true signed slot rotations modulo N/2 and explicit conjugation, with matching evaluation keys; update affected docs and tests.
- [x] Add a public evaluation-key collection with level/context checks and pre-generated relin/rotation/conjugation keys. Setup receives secret key and RNG; evaluation does not.
- [x] Test asymmetric complex vectors, opposite/large/zero offsets, conjugation twice, missing/wrong-level keys and preservation of existing evaluator callers.

## Context

Design and candidate profile: `specs/research/bootstrap.md`. Existing negative rotations include conjugation; correcting them is an intentional behavior change approved by the operator. Expose conjugation separately. No dependency on 0017/0022. Implemented; awaiting final operator review.

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

Implemented bootstrap automorphisms keys.

## Log

- 2026-11-01 created: bootstrap automorphisms keys
- 2026-11-01 scope add src/rings/traits.rs — signed automorphisms and public evaluation keys
- 2026-11-01 scope add src/ckks/keys/key_switching.rs — signed automorphisms and public evaluation keys
- 2026-11-01 scope add src/ckks/keys/public_key.rs — signed automorphisms and public evaluation keys
- 2026-11-01 scope add src/ckks/keys/relin_key.rs — signed automorphisms and public evaluation keys
- 2026-11-01 scope add src/ckks/keys/rotation_key.rs — signed automorphisms and public evaluation keys
- 2026-11-01 scope add src/ckks/keys/secret_key.rs — signed automorphisms and public evaluation keys
- 2026-11-01 scope add src/ckks/ops.rs — signed automorphisms and public evaluation keys
- 2026-11-01 scope add tests/integration_mul.rs — signed automorphisms and public evaluation keys
- 2026-11-01 scope rm src/rings/traits.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm src/ckks/keys/key_switching.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm src/ckks/keys/public_key.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm src/ckks/keys/relin_key.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm src/ckks/keys/rotation_key.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm src/ckks/keys/secret_key.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm src/ckks/ops.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm tests/integration_mul.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope add src/rings/traits.rs — declare complete module scope including new files
- 2026-11-01 scope add src/ckks/keys/** — declare complete module scope including new files
- 2026-11-01 scope add src/ckks/keys.rs — declare complete module scope including new files
- 2026-11-01 scope add src/ckks/evaluator.rs — declare complete module scope including new files
- 2026-11-01 scope add src/ckks/ops.rs — declare complete module scope including new files
- 2026-11-01 scope add tests/** — declare complete module scope including new files
- 2026-11-01 draft approved: Implement the plan.
- 2026-11-01 advance draft → ready
- 2026-11-01 scope approved: src/rings/traits.rs, src/ckks/keys/**, src/ckks/keys.rs, src/ckks/evaluator.rs, src/ckks/ops.rs, tests/**
- 2026-11-01 advance ready → in-progress/implement
- 2026-11-01 baseline: README.md@f33b4f4fc917, examples/bootstrap_demo.rs@38295abf01b4, src/ckks/algorithms.rs@6f2f1e43a81b, src/ckks/algorithms/bootstrap.rs@0ad1e98c3a13, src/rings/backends/rns_ntt/poly.rs@5c302b70b49f, src/rings/errors.rs@fea138cd1c12
- 2026-11-01 advance in-progress/implement → in-progress/verify
- 2026-11-01 advance in-progress/verify → in-progress/review
- 2026-11-01 advance in-progress/review → approval
- 2026-11-01 approved: ok, looks good.
- 2026-11-01 advance approval → done
