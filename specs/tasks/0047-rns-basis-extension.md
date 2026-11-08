---
id: 0047-rns-basis-extension
status: draft
scope: [src/rings/traits.rs, src/rings/backends/rns_ntt/*.rs]
created: 2026-11-03
source: Operator agreed large-degree roadmap, 2026-11-03
depends: [0046-bootstrap-stage-noise]
---
# Task: rns basis extension

## Context

Roadmap step 2, first half (`specs/research/large-degree-bootstrap.md`). Hybrid key switching needs an extended basis QP (the ciphertext primes Q plus special primes P), conversion of a digit from its primes to QP (ModUp, OpenFHE `ApproxSwitchCRTBasis`, Lattigo `BasisExtender`), and division by P back to Q (ModDown, OpenFHE `ApproxModDown`, Lattigo `ModDownQPtoQNTT`). This task adds them to the RNS backend and the `PolyRing` trait only; key switching moves in 0048.

Rules: touch only scoped files (extend with `specdev task scope` first). If the same check fails twice, stop, record the failure and a hypothesis in Findings and ask the operator. Never weaken tolerances, disable noise or add `#[ignore]`.

## Plan

- [ ] Basis: a chain can carry special primes P (distinct NTT primes, not used for rescale); contexts expose Q at a level, P, and QP at a level.
- [ ] Fast basis conversion from a subset of primes to another basis (exact up to the known small multiple of the source modulus), and ModDown from QP to Q with rounding.
- [ ] `PolyRing` gains the minimal methods for both, documented with their error bounds; no CKKS code changes.
- [ ] Tests: conversion and ModDown against big-integer references (crypto-bigint) at several levels, including the rounding bound.

## Acceptance

- [ ] `cargo fmt --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features`
- [ ] `cargo build --examples`
- [ ] `typos`
- [ ] `specdev check`
- [ ] `git diff --check`

## Log

- 2026-11-03 created: rns basis extension
- 2026-11-03 set source Operator agreed large-degree roadmap, 2026-11-03
- 2026-11-03 set depends 0046-bootstrap-stage-noise
- 2026-11-03 scope add src/rings/traits.rs — Special primes, fast basis conversion and ModDown in the RNS backend
- 2026-11-03 scope add src/rings/backends/rns_ntt/*.rs — Special primes, fast basis conversion and ModDown in the RNS backend
