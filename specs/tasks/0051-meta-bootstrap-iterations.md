---
id: 0051-meta-bootstrap-iterations
status: draft
scope: [src/ckks/algorithms/bootstrap.rs, src/ckks/algorithms/bootstrap/*.rs, tests/bootstrap.rs]
created: 2026-11-03
source: Operator agreed large-degree roadmap, 2026-11-03
depends: [0048-hybrid-key-switching]
---
# Task: meta bootstrap iterations

## Context

Roadmap step 6 (`specs/research/large-degree-bootstrap.md`). Meta-BTS (eprint 2022/1167; Lattigo `IterationsParameters`, OpenFHE `EvalBootstrap(ct, 2, precision)`) bootstraps the scaled error of a first bootstrap and subtracts it, roughly doubling precision for twice the cost and one reserved level. Worth it only after hybrid key switching (0048) raises single-bootstrap precision.

Rules: touch only scoped files (extend with `specdev task scope` first). If the same check fails twice, stop, record the failure and a hypothesis in Findings and ask the operator. Never weaken tolerances, disable noise or add `#[ignore]`.

## Plan

- [ ] Bootstrap takes an iteration count (1 or 2) and a precision in bits; the schedule reserves the extra level and checks headroom for the scaled error.
- [ ] Second iteration: `e = 2^p * (refreshed - input)`, bootstrap `e`, subtract `e * 2^-p`; inputs must be at the same level and scale.
- [ ] Encrypted tests at N=16/64 show the precision gain; record N=1024/8192 results in the research doc.

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

- 2026-11-03 created: meta bootstrap iterations
- 2026-11-03 set source Operator agreed large-degree roadmap, 2026-11-03
- 2026-11-03 set depends 0048-hybrid-key-switching
- 2026-11-03 scope add src/ckks/algorithms/bootstrap.rs — Iterated bootstrap
- 2026-11-03 scope add src/ckks/algorithms/bootstrap/*.rs — Iterated bootstrap
- 2026-11-03 scope add tests/bootstrap.rs — Iterated bootstrap
