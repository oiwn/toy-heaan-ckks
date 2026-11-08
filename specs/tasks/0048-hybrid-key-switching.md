---
id: 0048-hybrid-key-switching
status: draft
scope: [src/ckks/keys/*.rs, src/ckks/keys.rs, src/ckks/ops.rs, src/ckks/engine.rs, src/ckks/algorithms/bootstrap.rs, examples/**/*.rs]
created: 2026-11-03
source: Operator agreed large-degree roadmap, 2026-11-03
depends: [0047-rns-basis-extension]
---
# Task: hybrid key switching

## Context

Roadmap step 2, second half (`specs/research/large-degree-bootstrap.md`). Our BV key switch (20-bit digits, no special modulus) adds about `sqrt(gadget_len * N) * 2^20 * std_dev` per switch, which dominates bootstrap error at large N. OpenFHE (`HYBRID`, 60-bit special primes, dnum 3) and Lattigo (`GadgetProductLazy` + ModDown, LogP 4-6 x 61) divide that term by P, leaving rounding-sized noise. Built on the basis extension from 0047.

Open design points for the draft review: the number of digits (dnum) and P primes as parameters or derived; whether BV stays as an option for comparison; key memory at N=8192 (today about 7.5 GB).

Rules: touch only scoped files (extend with `specdev task scope` first). If the same check fails twice, stop, record the failure and a hypothesis in Findings and ask the operator. Never weaken tolerances, disable noise or add `#[ignore]`.

## Plan

- [ ] `KeySwitchingKey` uses hybrid RNS digits: keys live in QP, `switch` does ModUp per digit, the inner product in QP and one ModDown.
- [ ] Relinearization, rotation, conjugation and bootstrap keys move to it; engines and builders take the special-prime parameters.
- [ ] Key-switch noise test against the expected rounding bound; all existing tests and demos pass.
- [ ] Rerun the large-degree ratio sweep and record the new error and key memory; refit the 0045 ratio rule if the optimum moved.

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

- 2026-11-03 created: hybrid key switching
- 2026-11-03 set source Operator agreed large-degree roadmap, 2026-11-03
- 2026-11-03 set depends 0047-rns-basis-extension
- 2026-11-03 scope add src/ckks/keys/*.rs — Hybrid key-switching keys and their users
- 2026-11-03 scope add src/ckks/keys.rs — Hybrid key-switching keys and their users
- 2026-11-03 scope add src/ckks/ops.rs — Hybrid key-switching keys and their users
- 2026-11-03 scope add src/ckks/engine.rs — Hybrid key-switching keys and their users
- 2026-11-03 scope add src/ckks/builder.rs — Hybrid key-switching keys and their users
- 2026-11-03 scope add src/ckks/algorithms/bootstrap.rs — Hybrid key-switching keys and their users
- 2026-11-03 scope add examples/**/*.rs — Hybrid key-switching keys and their users
- 2026-11-03 scope rm src/ckks/builder.rs — Unused builder; add back only if key parameters need it
