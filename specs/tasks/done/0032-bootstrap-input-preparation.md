---
id: 0032-bootstrap-input-preparation
status: done
scope: [src/ckks/algorithms/bootstrap.rs, examples/bootstrap_demo.rs, README.md]
created: 2026-11-01
source: Operator request: move bootstrap scale and level preparation into the library
depends: [0031-bootstrap-integration]
---
# Task: bootstrap input preparation

## Plan

- [x] Validate input/key contexts and scale; normalize a copy to the configured target and reduce to the bottom level inside bootstrap.
- [x] Remove caller-side scale restoration and level dropping from the demo; document accepted inputs without promising recovered precision.
- [x] Test equivalence to explicit preparation, input immutability, fresh/already-prepared/post-multiply inputs and invalid contexts/scales/keys; retain multi-seed repeated-refresh coverage.

## Context

Follow-up to accepted bootstrap tasks0016 and0026-0031. Preserve explicit multiplication/rescale, method signature and output contract. Accept positive scales up to the configured target; reject larger/zero scales. Skip adjustment for matching scales. Context compatibility must be checked before arithmetic; evaluation receives no secret or RNG.

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features`
- [x] `cargo build --examples`
- [x] `cargo run --release --example bootstrap_demo`
- [x] `specdev check`
- [x] `typos`
- [x] `git diff --check`

## Review

Preparation validates both component contexts, key context, chain membership and supported scales before arithmetic. A matching scale skips adjustment; scale increases multiply the representation and interpretation equally. Exact component differences confirm equivalence to the former manual sequence and input immutability. Higher levels, prepared inputs, invalid scales and wrong contexts are covered at both degrees. Bootstrap internals and output contract are unchanged. The release suite and demo passed with the same demo errors as before.

## Summary

Bootstrap prepares input scale and level internally so callers only multiply, rescale and request refresh.

## Log

- 2026-11-01 created: bootstrap input preparation
- 2026-11-01 scope add src/ckks/algorithms/bootstrap.rs — automatic bootstrap scale and level preparation
- 2026-11-01 scope add examples/bootstrap_demo.rs — automatic bootstrap scale and level preparation
- 2026-11-01 scope add README.md — automatic bootstrap scale and level preparation
- 2026-11-01 draft approved: Implement the plan.
- 2026-11-01 advance draft → ready
- 2026-11-01 scope approved: src/ckks/algorithms/bootstrap.rs, examples/bootstrap_demo.rs, README.md
- 2026-11-01 advance ready → in-progress/implement
- 2026-11-01 baseline: src/ckks/algorithms.rs@6f2f1e43a81b, src/ckks/keys.rs@d0b53f26bc7e, src/ckks/keys/conjugation_key.rs@eddfd1c737f7, src/ckks/keys/rotation_key.rs@4a8d2128faf5, src/ckks/ops.rs@d9424d9405fd, src/rings/backends/rns_ntt/poly.rs@5c302b70b49f, src/rings/errors.rs@fea138cd1c12, src/rings/traits.rs@8a54adc1cbc6
- 2026-11-01 set depends 0031-bootstrap-integration
- 2026-11-01 set source Operator request: move bootstrap scale and level preparation into the library
- 2026-11-01 advance in-progress/implement → in-progress/verify
- 2026-11-01 advance in-progress/verify → in-progress/review
- 2026-11-01 advance in-progress/review → approval
- 2026-11-01 approved: Ok it works fine, let's to bookkeeping and close the case.
- 2026-11-01 advance approval → done
