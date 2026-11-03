---
id: 0045-bootstrap-ratio-default
status: draft
scope: [src/ckks/algorithms/bootstrap.rs, src/ckks/algorithms/bootstrap/*.rs, tests/bootstrap.rs, examples/bootstrap_large.rs, examples/bootstrap_demo.rs]
created: 2026-11-03
source: Operator agreed large-degree roadmap, 2026-11-03
depends: [0043-bootstrap-general-final]
---
# Task: bootstrap ratio default

## Context

Roadmap step 1 of `specs/research/large-degree-bootstrap.md` (operator agreed 2026-11-03). The reference ratio 17 fails above N=128; the ratio sweep shows the error is U-shaped in the ratio (noise amplification above the optimum, cubic sine error below), with the optimum moving down with N: about 10 at N<=64, 6 at N=1024, 4 at N=8192. No fixed ratio passes everywhere, so `BootstrapConfig::reference(N)` picks it from N, like OpenFHE's fitted correction factor but from our own data. The fit is redone after hybrid key switching (0048).

A smaller ratio means a larger message scale (2^(61 - ratio)): the output level must keep headroom for one product after bootstrap, and recurrence (0044) relies on the current 2^44.

Rules: touch only scoped files (extend with `specdev task scope` first). If the same check fails twice, stop, record the failure and a hypothesis in Findings and ask the operator. Never weaken tolerances, disable noise or add `#[ignore]`.

## Plan

- [ ] Fill the sweep at N=16/32/128/256/512/2048/4096 (`bootstrap_large` gets N=16/32; `bootstrap_demo` for N=16) with weight N/2 up to 64 and record it in the research doc.
- [ ] `BootstrapConfig::reference(N)` derives `message_ratio_bits` from N with a documented rule fitted to the sweep (a step table or a line in log2 N); `bootstrap_large` uses the reference ratio by default.
- [ ] Headroom: the schedule checks that the output level can square the output once; error names the level.
- [ ] Regression numbers in tests/demos follow the new defaults; the 1e-3 gate and noise stay; record old and new errors in Findings.

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

- 2026-11-03 created: bootstrap ratio default
- 2026-11-03 set source Operator agreed large-degree roadmap, 2026-11-03
- 2026-11-03 set depends 0043-bootstrap-general-final
- 2026-11-03 scope add src/ckks/algorithms/bootstrap.rs — Degree-dependent ratio default, output headroom check, regression numbers
- 2026-11-03 scope add src/ckks/algorithms/bootstrap/*.rs — Degree-dependent ratio default, output headroom check, regression numbers
- 2026-11-03 scope add tests/bootstrap.rs — Degree-dependent ratio default, output headroom check, regression numbers
- 2026-11-03 scope add examples/bootstrap_large.rs — Degree-dependent ratio default, output headroom check, regression numbers
- 2026-11-03 scope add examples/bootstrap_demo.rs — Degree-dependent ratio default, output headroom check, regression numbers
