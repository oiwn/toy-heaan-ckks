---
id: 0034-nonlinear-smoothing-demo
status: draft
scope: [examples/**, src/ckks/algorithms/**, README.md]
created: 2026-11-01
depends: [0033-bootstrapped-recurrence-demo]
---
# Task: nonlinear smoothing demo

## Plan

- [ ] Design neighbor averaging (left+2*center+right)/4 followed by the bounded quadratic recurrence, with an encrypted driver.
- [ ] Pre-generate rotation/relin keys and derive the scale/level schedule; choose a depth requiring refresh and compare short-chain failure and long-chain completion.
- [ ] Implement a toy example with plaintext checkpoints, refresh/rotation counts, timing and accuracy tests.

## Context

Future variant of0033; draft only. Tests must cover boundary/mixed-sign inputs, correct cyclic neighbors, input immutability, invalid/missing keys and domain preservation. Agree concrete profiles and tolerances before implementation. More representative than independent slots, with additional rotation/key cost. No total-memory claim from prime counts alone.

## Acceptance

- [ ] `cargo fmt --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features`

## Log

- 2026-11-01 created: nonlinear smoothing demo
- 2026-11-01 scope add examples/** — future smoothing demo
- 2026-11-01 scope add src/ckks/algorithms/** — future smoothing demo
- 2026-11-01 scope add README.md — future smoothing demo
- 2026-11-01 set depends 0033-bootstrapped-recurrence-demo
