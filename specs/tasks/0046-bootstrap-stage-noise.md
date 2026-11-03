---
id: 0046-bootstrap-stage-noise
status: draft
scope: [tests/bootstrap.rs, src/ckks/algorithms/bootstrap.rs]
created: 2026-11-03
source: Operator agreed large-degree roadmap, 2026-11-03
depends: [0043-bootstrap-general-final]
---
# Task: bootstrap stage noise

## Context

Roadmap step 2 prerequisite in `specs/research/large-degree-bootstrap.md`. Key switching dominates the bootstrap error at 20-bit digits, but with 10-bit digits a second ratio-amplified source remains (N=1024 ratio 17: 3.9e-3) and grows fast with N. Its origin decides how much hybrid key switching (0047/0048) can gain, so it is measured first.

Diagnostics only: library results stay bit-identical. Decryption is allowed in tests and diagnostics, never in evaluation.

Rules: touch only scoped files (extend with `specdev task scope` first). If the same check fails twice, stop, record the failure and a hypothesis in Findings and ask the operator. Never weaken tolerances, disable noise or add `#[ignore]`.

## Plan

- [ ] A test-only stage trace: run the bootstrap stages one by one (ModRaise, each CtS group, conjugate split, EvalMod, each StC group), decrypt after each and compare with a plaintext f64 pipeline of the same stages; report the error in units of the stage's message.
- [ ] Run it at N=256/1024 with 20- and 10-bit digits (temporary, reverted) and ratios 17/8; separate key-switch, rescale-rounding, diagonal-encoding and ModRaise-noise contributions.
- [ ] Record the per-stage table, the named second source and a fix proposal in the research doc; update the roadmap (step 3) and the expected gain of 0048.

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

- 2026-11-03 created: bootstrap stage noise
- 2026-11-03 set source Operator agreed large-degree roadmap, 2026-11-03
- 2026-11-03 set depends 0043-bootstrap-general-final
- 2026-11-03 scope add tests/bootstrap.rs — Test-only per-stage noise trace (stage access may need a cfg(test) hook)
- 2026-11-03 scope add src/ckks/algorithms/bootstrap.rs — Test-only per-stage noise trace (stage access may need a cfg(test) hook)
