---
id: 0035-cubic-activation-demo
status: draft
scope: [examples/**, src/ckks/algorithms/**, README.md]
created: 2026-11-01
depends: [0033-bootstrapped-recurrence-demo]
---
# Task: cubic activation demo

## Plan

- [ ] Design repeated activation f(x)=(3*x-x^3)/2 on real slots in [-1,1], deriving the two-multiplication-layer schedule and refresh boundaries.
- [ ] Select inputs/round count that expose useful depth before saturation; compare short-chain failure and long-chain completion against plaintext.
- [ ] Implement the toy demo and test range, saturation, error propagation, scale/depth and invalid-key/context cases.

## Context

Future variant of0033; draft only. This resembles polynomial activation layers but approaches +/-1, making late rounds less informative. Decide concrete profiles and accuracy before implementation; report saturation and avoid presenting convergence as evidence that errors are small. Evaluation must use public keys with nonzero noise, explicit refresh and no caller bit counting.

## Acceptance

- [ ] `cargo fmt --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features`

## Log

- 2026-11-01 created: cubic activation demo
- 2026-11-01 scope add examples/** — future activation demo
- 2026-11-01 scope add src/ckks/algorithms/** — future activation demo
- 2026-11-01 scope add README.md — future activation demo
- 2026-11-01 set depends 0033-bootstrapped-recurrence-demo
