---
id: 0010-chebyshev-module
status: draft
scope: [src/ckks/**, src/ckks.rs, examples/**]
created: 2026-10-13
source: user request 2026-10-13: pipeline.rs is a std-dev grab bag; decompose into composable ops like fhe.rs / tfhe-rs
depends: [0009-ops-composition]
---
# Task: chebyshev module

## Plan

- [ ] `src/ckks/chebyshev.rs`: encrypted Clenshaw (`eval_chebyshev_encrypted`) rewritten on `ops` + `Ciphertext::rescale`/`mod_switch_to`, plus `fit_chebyshev`, plaintext `eval_chebyshev_clenshaw` and the domain maps.
- [ ] Move whatever remains of `src/ckks/pipeline.rs` to its home; move the empty file into `tmp/`.
- [ ] Update examples (sensor_std_dev, sensor_std_plain, sqrt_chebyshev_plain); Chebyshev tests move with the code and stay green.

## Acceptance

- [ ] `cargo fmt --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features`
- [ ] `cargo build --examples`
- [ ] `cargo run --release --example sensor_std_dev`
- [ ] `typos`
- [ ] `gitleaks dir .`

## Context

On hold: revised after the backend design in `specs/ctx.md` is agreed.

## Log

- 2026-10-13 created: chebyshev module
- 2026-10-13 set source user request 2026-10-13: pipeline.rs is a std-dev grab bag; decompose into composable ops like fhe.rs / tfhe-rs
- 2026-10-13 set depends 0009-ops-composition
- 2026-10-13 scope add src/ckks/** — chebyshev module, retire pipeline.rs
- 2026-10-13 scope add src/ckks.rs — chebyshev module, retire pipeline.rs
- 2026-10-13 scope add examples/** — chebyshev module, retire pipeline.rs
