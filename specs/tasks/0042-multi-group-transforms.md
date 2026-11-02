---
id: 0042-multi-group-transforms
status: draft
scope: [src/ckks/algorithms/bootstrap.rs, src/ckks/algorithms/bootstrap/transforms.rs, tests/bootstrap.rs]
created: 2026-11-02
source: Operator-approved bootstrap split plan 2026-11-02 (from 0037)
depends: [0041-bootstrap-config-schedule]
---
# Task: multi group transforms

## Context

Connects the existing factored transforms (`bootstrap/transforms.rs`: sparse radix-two factors, configurable grouping) to encrypted evaluation at the levels the 0041 schedule assigns, with keys derived per level.

Rules: touch only scoped files (extend with `specdev task scope` first). If the same check fails twice, stop, record the failure and a hypothesis in Findings and ask the operator. Never weaken tolerances, disable noise or add `#[ignore]`.

## Plan

- [ ] Encode each transform group at its scheduled level and scale; store explicit (rotation offset, plaintext) pairs instead of a dense diagonal list.
- [ ] Fuse normalization/restoration into the first/last factor.
- [ ] Schedule derives the rotation/conjugation/relin key set per level from the groups; `generate_keys` follows it.
- [ ] `bootstrap` validates key availability before any evaluation.
- [ ] Tests: encrypted multi-group transforms match the dense reference at N=16/32/64/128 for several budgets; bootstrap passes with budgets > 1; missing-key error.

## Acceptance

- [ ] `cargo fmt --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features`
- [ ] `cargo build --examples`
- [ ] `cargo run --release --example bootstrap_demo`
- [ ] `typos`
- [ ] `specdev check`
- [ ] `git diff --check`

## Log

- 2026-11-02 created: multi group transforms
- 2026-11-02 set source Operator-approved bootstrap split plan 2026-11-02 (from 0037)
- 2026-11-02 set depends 0041-bootstrap-config-schedule
- 2026-11-02 scope add src/ckks/algorithms/bootstrap.rs — Multi-group transforms at scheduled levels with derived keys
- 2026-11-02 scope add src/ckks/algorithms/bootstrap/transforms.rs — Multi-group transforms at scheduled levels with derived keys
- 2026-11-02 scope add tests/bootstrap.rs — Multi-group transforms at scheduled levels with derived keys
