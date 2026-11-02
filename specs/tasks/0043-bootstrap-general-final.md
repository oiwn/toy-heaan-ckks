---
id: 0043-bootstrap-general-final
status: draft
scope: [src/ckks/algorithms/bootstrap.rs, src/ckks/algorithms/recurrence.rs, examples/**/*.rs, tests/bootstrap.rs]
created: 2026-11-02
source: Operator-approved bootstrap split plan 2026-11-02 (from 0037)
depends: [0042-multi-group-transforms]
---
# Task: bootstrap general final

## Context

Last task of the chain: remove the remaining profile restrictions and close the general-bootstrap validation matrix in `specs/research/generic-bootstrap.md`. After acceptance, 0037 returns to the file review.

Rules: touch only scoped files (extend with `specdev task scope` first). If the same check fails twice, stop, record the failure and a hypothesis in Findings and ask the operator. Never weaken tolerances, disable noise or add `#[ignore]`.

## Plan

- [ ] Remove the degree/profile whitelist, the exact-18-prime and near-2^61 proximity checks and fixed offsets; replace them with headroom, level and representable-scale validation.
- [ ] Recurrence uses actual scales and derived bootstrap output properties, keeping its formula.
- [ ] Update `bootstrap_demo`, `recurrence_bootstrap_demo` and other affected examples; document the breaking API in overview and the changelog summary.
- [ ] Validation matrix: N=16/32/64/128; varied scales, secret weights, polynomial settings, budgets, chain lengths and mixed widths; original 8 seeds x 2 cycles at N=16/64 within 1e-3.
- [ ] Invalid configuration, wrong context, excessive secret weight, missing keys, insufficient headroom and input immutability tests.
- [ ] Mark generic-bootstrap.md complete.

## Manual checks

- [ ] Operator reviews the full general-bootstrap integration.

## Acceptance

- [ ] `cargo fmt --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features`
- [ ] `cargo build --examples`
- [ ] `cargo run --release --example bootstrap_demo`
- [ ] `cargo run --release --example recurrence_bootstrap_demo`
- [ ] `typos`
- [ ] `specdev check`
- [ ] `git diff --check`

## Log

- 2026-11-02 created: bootstrap general final
- 2026-11-02 set source Operator-approved bootstrap split plan 2026-11-02 (from 0037)
- 2026-11-02 set depends 0042-multi-group-transforms
- 2026-11-02 scope add src/ckks/algorithms/bootstrap.rs — Remove profile restrictions, migrate recurrence/examples, final validation matrix
- 2026-11-02 scope add src/ckks/algorithms/recurrence.rs — Remove profile restrictions, migrate recurrence/examples, final validation matrix
- 2026-11-02 scope add examples/**/*.rs — Remove profile restrictions, migrate recurrence/examples, final validation matrix
- 2026-11-02 scope add tests/bootstrap.rs — Remove profile restrictions, migrate recurrence/examples, final validation matrix
