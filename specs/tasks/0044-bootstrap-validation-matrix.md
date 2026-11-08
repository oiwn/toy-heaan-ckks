---
id: 0044-bootstrap-validation-matrix
status: draft
scope: [src/ckks/algorithms/recurrence.rs, examples/**/*.rs, tests/bootstrap.rs, src/ckks/algorithms/bootstrap.rs]
created: 2026-11-02
source: Operator split of 0043, 2026-11-02
depends: [0043-bootstrap-general-final, 0045-bootstrap-ratio-default]
---
# Task: bootstrap validation matrix

## Context

Second half of the old 0043 (split by the operator 2026-11-02). 0043 generalizes bootstrap validation; this task migrates its users, runs the final validation matrix in `specs/research/generic-bootstrap.md` and marks it complete. After acceptance, 0037 returns to the file review.

^^^ wtf is "mgrates it's users"? which one?

Large-degree behaviour (N up to 8192, `examples/bootstrap_large.rs`) is recorded in `specs/research/large-degree-bootstrap.md`: the reference ratio 17 does not hold above N=128, and the error follows q0/message. The matrix uses the degree-dependent reference ratio from 0045.

Rules: touch only scoped files (extend with `specdev task scope` first). If the same check fails twice, stop, record the failure and a hypothesis in Findings and ask the operator. Never weaken tolerances, disable noise or add `#[ignore]`.

## Plan

- [ ] Recurrence drops its N=16/64 and near-2^61 whitelist; it uses actual scales and the bootstrap's derived message scale and output properties, keeping its formula.
- [ ] Update `bootstrap_demo`, `recurrence_bootstrap_demo` and other affected examples; document the breaking API (`message_ratio_bits`) in overview and the changelog summary.
- [ ] Validation matrix: N=16/32/64/128; varied scales, secret weights, polynomial settings, budgets, chain lengths and prime widths; original 8 seeds x 2 cycles at N=16/64 within 1e-3.
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

- 2026-11-02 created: bootstrap validation matrix
- 2026-11-02 set source Operator split of 0043, 2026-11-02
- 2026-11-02 set depends 0043-bootstrap-general-final
- 2026-11-02 scope add src/ckks/algorithms/recurrence.rs — Recurrence/examples migration and the final validation matrix
- 2026-11-02 scope add examples/**/*.rs — Recurrence/examples migration and the final validation matrix
- 2026-11-02 scope add tests/bootstrap.rs — Recurrence/examples migration and the final validation matrix
- 2026-11-02 scope add src/ckks/algorithms/bootstrap.rs — Recurrence/examples migration and the final validation matrix
- 2026-11-03 set depends 0043-bootstrap-general-final, 0045-bootstrap-ratio-default
