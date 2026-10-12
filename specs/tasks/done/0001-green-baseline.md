---
id: 0001-green-baseline
status: done
scope: [examples/debug_std_dev.rs, examples/rotation_demo.rs, examples/sensor_std_plain.rs, src/crypto/pipeline.rs, src/rings/backends/rns_ntt/basis.rs, examples/std_dev_8.rs, tests/integration_mul.rs]
created: 2026-10-12
source: user request 2026-10-12: revive project, CKKS green + native sensor demo with range normalization
---
# Task: green baseline

## Plan

- [x] Delete `examples/debug_std_dev.rs` (scratch from the Clenshaw debugging).
- [x] Fix clippy errors at `examples/rotation_demo.rs:104` (redundant reference in `println!`) and `examples/sensor_std_plain.rs:217` (binary expression can be simplified).
- [x] Run `cargo fmt` on the files with diffs.
- [x] AGENTS.md Quick Reference: replace the removed examples `rns`, `ntt`, `rns_mod_switch` with real ones (`encrypt_mul`, `rotation_demo`, `std_dev_8`).
- [x] AGENTS.md: `cargo bench --bench end_to_end` -> `cargo bench --bench primes` (`end_to_end` does not exist).

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features -- --skip encrypted_std_dev_end_to_end` (that test is fixed by 0002)

## Summary

Green baseline on rustc 1.99: removed scratch example, fixed 2 clippy errors, applied cargo fmt, fixed stale AGENTS.md commands.

## Context

Checked 2026-10-12 on rustc 1.99: 114/115 lib tests, 7 integration tests, and all 11 examples pass; the red test is the known direct-fit design. The clippy errors in committed examples likely come from lints newer than the old 1.93 pin.

## Review

- Diff matches the plan; no behavior change. `{:?}` -> `{input:?}` and `lower <= 0.0 && lower != 0.0` -> `lower < 0.0` are equivalent.
- `cargo fmt` touched only scoped files; `specdev check` shows no out-of-scope changes.
- Extra checks: `typos` clean; examples `rotation_demo`, `sensor_std_plain`, `std_dev_8` run OK.
- Known red test `encrypted_std_dev_end_to_end` stays skipped until 0002.

## Log

- 2026-10-12 created: green baseline
- 2026-10-12 set source user request 2026-10-12: revive project, CKKS green + native sensor demo with range normalization
- 2026-10-12 scope add examples/debug_std_dev.rs — scratch file from Clenshaw debugging; delete
- 2026-10-12 scope add examples/rotation_demo.rs — clippy errors under toolchain 1.99
- 2026-10-12 scope add examples/sensor_std_plain.rs — clippy errors under toolchain 1.99
- 2026-10-12 scope add src/**/*.rs — cargo fmt touches many files (pipeline.rs, std_dev_8.rs, basis.rs, integration_mul.rs)
- 2026-10-12 scope add examples/*.rs — cargo fmt touches many files (pipeline.rs, std_dev_8.rs, basis.rs, integration_mul.rs)
- 2026-10-12 scope add tests/*.rs — cargo fmt touches many files (pipeline.rs, std_dev_8.rs, basis.rs, integration_mul.rs)
- 2026-10-12 scope rm src/**/*.rs
- 2026-10-12 scope rm examples/*.rs
- 2026-10-12 scope rm tests/*.rs
- 2026-10-12 scope add src/crypto/pipeline.rs — cargo fmt --check diffs
- 2026-10-12 scope add src/rings/backends/rns_ntt/basis.rs — cargo fmt --check diffs
- 2026-10-12 scope add examples/std_dev_8.rs — cargo fmt --check diffs
- 2026-10-12 scope add tests/integration_mul.rs — cargo fmt --check diffs
- 2026-10-12 advance draft → ready
- 2026-10-12 scope approved: examples/debug_std_dev.rs, examples/rotation_demo.rs, examples/sensor_std_plain.rs, src/crypto/pipeline.rs, src/rings/backends/rns_ntt/basis.rs, examples/std_dev_8.rs, tests/integration_mul.rs
- 2026-10-12 advance ready → in-progress/implement
- 2026-10-12 baseline: .codex/hooks.json@f274ace3f299, .gitignore@800d14263080, rust-toolchain.toml@275861857b47, specdev.toml@1ea52da4a446, src/crypto/engine.rs@1a784111cbbc, src/crypto/mod.rs@da20b5b1922b, src/crypto/types.rs@2cf43168b9fe, src/rings/backends/rns_ntt/poly.rs@094f2ffe24e8
- 2026-10-12 advance in-progress/implement → in-progress/verify
- 2026-10-12 advance in-progress/verify → in-progress/review
- 2026-10-12 advance in-progress/review → approval
- 2026-10-12 approved: this one will fix later ? ther rest looks good
- 2026-10-12 advance approval → done
