---
id: 0050-bootstrap-arcsine-correction
status: draft
scope: [src/ckks/algorithms/bootstrap/*.rs, src/ckks/algorithms/bootstrap.rs, tests/bootstrap.rs]
created: 2026-11-03
source: Operator agreed large-degree roadmap, 2026-11-03
depends: [0045-bootstrap-ratio-default]
---
# Task: bootstrap arcsine correction

## Context

Roadmap step 5 (`specs/research/large-degree-bootstrap.md`). Below the optimal ratio the bootstrap error is the cubic sine linearization error (16x per 2 bits; N=64 ratio 4: 2.7e-3). Lattigo composes the sine with a small arcsine polynomial (`Mod1InvDegree` 5 or 7) and reaches ratio 2 with 32 bits at N=2^16. With it the ratio can drop where the sine floor limits it (N<=1024), lowering noise amplification.

Rules: touch only scoped files (extend with `specdev task scope` first). If the same check fails twice, stop, record the failure and a hypothesis in Findings and ask the operator. Never weaken tolerances, disable noise or add `#[ignore]`.

## Plan

- [ ] `BootstrapConfig` gains an optional arcsine degree; the schedule compiles it as a second `ChebyshevPlan` stage after the double angles, with its level cost and scale checks.
- [ ] The plaintext reference and the encrypted EvalMod use the composed function; tests cover degree 0 (off, unchanged results) and 5/7.
- [ ] Sweep the ratio with and without the correction at N=64/1024; refit the 0045 rule if it helps and record the data.

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

- 2026-11-03 created: bootstrap arcsine correction
- 2026-11-03 set source Operator agreed large-degree roadmap, 2026-11-03
- 2026-11-03 set depends 0045-bootstrap-ratio-default
- 2026-11-03 scope add src/ckks/algorithms/bootstrap/*.rs — Arcsine stage in EvalMod
- 2026-11-03 scope add src/ckks/algorithms/bootstrap.rs — Arcsine stage in EvalMod
- 2026-11-03 scope add tests/bootstrap.rs — Arcsine stage in EvalMod
