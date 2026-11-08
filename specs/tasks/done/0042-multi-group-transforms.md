---
id: 0042-multi-group-transforms
status: done
scope: [src/ckks/algorithms/bootstrap.rs, tests/bootstrap.rs, src/ckks/algorithms/bootstrap/*.rs]
created: 2026-11-02
source: Operator-approved bootstrap split plan 2026-11-02 (from 0037)
depends: [0041-bootstrap-config-schedule]
---
# Task: multi group transforms

## Context

Connects the existing factored transforms (`bootstrap/transforms.rs`: sparse radix-two factors, configurable grouping) to encrypted evaluation at the levels the 0041 schedule assigns, with keys derived per level.

Rules: touch only scoped files (extend with `specdev task scope` first). If the same check fails twice, stop, record the failure and a hypothesis in Findings and ask the operator. Never weaken tolerances, disable noise or add `#[ignore]`.

## Plan

- [x] `bootstrap/linear.rs`: `EncodedTransform` encodes each group at its scheduled level and plaintext scale as explicit (offset, plaintext) pairs for the nonzero diagonals. Evaluation takes a rotation callback; it replaces `BootstrapPlan::transform`/`encode_diagonals`.
- [x] Fuse normalization into the last CtS group (see Findings; the plan said first) and restoration into the last StC group. `StageSchedule.plaintext_scale`: working scale everywhere except the last StC group (message scale).
- [x] Budgets `1..=log2(N/2)` pass validation (the rest of the whitelist stays until 0043). Rotation keys per level follow the groups (already derived per factor since 0041).
- [x] `BootstrapKeys::check` runs before any evaluation and names the missing key and level.
- [x] Tests: encrypted groups match the plaintext factored transform at N=16/32/64/128 for budgets 1/2/log2(slots); bootstrap passes 1e-3 with budgets (2,2), (3,3), (1,3), (3,1); a chain too short is rejected; missing-key error.

## Findings

- **CtS normalization moved to the last group (deviation from the plan).** With `1/(2*radius)` in the first CtS group, multi-group bootstrap missed the gate at N=64: (2,2) reached 3.4e-3 to 6.7e-3, and (3,1) reached 1.4e-2. StC-only grouping was unaffected. Probing per budget showed that only CtS grouping degraded, about 20x at N=16 and 80x at N=64, roughly in proportion to the radius (6 and 18). Every later group adds noise to a signal that is already 2*radius smaller. With the normalization fused into the last group, all budgets stay within 5x of baseline (worst 4.8e-4). With budgets 1/1 the first and last groups coincide, so the reference profile is bit-identical to 0041.
- The standalone transform test runs at scale 2^55 with diagonals at 2^61, as bootstrap does. At 2^40, key-switch noise alone gives about 1e-4.
- A plan's missing-key error reports the first gap. A deeper CtS shifts every later level, so 1/1 keys usually lack a relinearization key first, not a rotation key.

## Review

- Scope: `bootstrap.rs`, `bootstrap/*.rs` (new `linear.rs`, `schedule.rs` gains `plaintext_scale` and `stages_of`) and `tests/bootstrap.rs`.
- The group evaluation keeps the old accumulation order, and the 1/1 demo and recurrence numbers match 0041 exactly.
- The multi-group error table and the depth rule (2*log2(N/2)+9) are in `specs/research/generic-bootstrap.md`.
- Out of scope, unchanged: baby-step/giant-step and hoisted rotations, the whitelist and the message ratio (0043).

## Summary

Bootstrap evaluates factored transforms group by group at scheduled levels, with per-level rotation keys, configurable CtS/StC budgets and an upfront key check.

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

- 2026-11-02 created: multi group transforms
- 2026-11-02 set source Operator-approved bootstrap split plan 2026-11-02 (from 0037)
- 2026-11-02 set depends 0041-bootstrap-config-schedule
- 2026-11-02 scope add src/ckks/algorithms/bootstrap.rs — Multi-group transforms at scheduled levels with derived keys
- 2026-11-02 scope add src/ckks/algorithms/bootstrap/transforms.rs — Multi-group transforms at scheduled levels with derived keys
- 2026-11-02 scope add tests/bootstrap.rs — Multi-group transforms at scheduled levels with derived keys
- 2026-11-02 scope rm src/ckks/algorithms/bootstrap/transforms.rs — covered by glob
- 2026-11-02 scope add src/ckks/algorithms/bootstrap/*.rs — New linear.rs; schedule and transforms changes
- 2026-11-02 draft approved: approved by operator
- 2026-11-02 advance draft → ready
- 2026-11-02 scope approved: src/ckks/algorithms/bootstrap.rs, tests/bootstrap.rs, src/ckks/algorithms/bootstrap/*.rs
- 2026-11-02 advance ready → in-progress/implement
- 2026-11-02 advance in-progress/implement → in-progress/verify
- 2026-11-02 advance in-progress/verify → in-progress/review
- 2026-11-02 advance in-progress/review → approval
- 2026-11-02 approved: approved by operator
- 2026-11-02 advance approval → done
