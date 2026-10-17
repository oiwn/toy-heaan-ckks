---
id: 0017-debug-tooling
status: draft
scope: [src/ckks/**, examples/**, Cargo.toml, README.md]
created: 2026-10-16
source: operator request 2026-10-16
---
# Task: debug tooling

Tooling to see what happens inside a computation: noise, precision, level and scale at each step. Merges the audit's `tracing` candidate.

## Plan

- [ ] Secret-key inspector behind a feature (`debug-oracle`): for a ciphertext and the expected values, report level, scale, measured noise and precision bits per slot.
- [ ] `tracing` spans/events in the ckks ops (mul, relin, rescale, rotate, key switch); a subscriber in the examples (`RUST_LOG`).
- [ ] Use the inspector in one example (`encrypt_mul` or `sensor_std_dev`) to show the noise per step.
- [ ] Document it in `specs/overview.md` and the README.

## Acceptance

- [ ] `cargo fmt --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features`
- [ ] `cargo build --examples`
- [ ] `cargo run --release --example sensor_std_dev`
- [ ] `typos`
- [ ] `gitleaks dir .`

## Log

- 2026-10-16 created: debug tooling
- 2026-10-16 set source operator request 2026-10-16
- 2026-10-16 scope add src/ckks/** — debug inspector and tracing
- 2026-10-16 scope add examples/** — debug inspector and tracing
- 2026-10-16 scope add Cargo.toml — debug inspector and tracing
- 2026-10-16 scope add README.md — document the inspector
