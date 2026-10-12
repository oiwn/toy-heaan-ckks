---
id: 0003-sensor-demo-example
status: ready
scope: [examples/sensor_std_plain.rs, examples/sensor_std_dev/main.rs]
created: 2026-10-12
source: user request 2026-10-12: revive project, CKKS green + native sensor demo with range normalization
depends: [0002-range-contract-std-dev]
---
# Task: sensor demo example

## Plan

- [ ] Extend `examples/sensor_std_dev/main.rs` (from 0002) to the demo parameters: N = 8192, 13 x 61-bit primes, degree 7, contract (-40, 60, 1, 8), 32 diurnal readings (one per 45 min) from a seeded generator with std ~2-3 degC.
- [ ] Print the stages (level, logp per op), the decrypted std, the f64 std, and both errors.
- [ ] Assert std error < 0.1 degC and intermediate variance error < 1e-9.
- [ ] Retarget `examples/sensor_std_plain.rs` to the range contract (plaintext reference).

## Acceptance

- [ ] `cargo fmt --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features`
- [ ] `cargo run --release --example sensor_std_dev`

## Summary

## Log

- 2026-10-12 created: sensor demo example
- 2026-10-12 set source user request 2026-10-12: revive project, CKKS green + native sensor demo with range normalization
- 2026-10-12 set depends 0002-range-contract-std-dev
- 2026-10-12 scope add examples/sensor_std_dev_demo.rs — new acceptance example
- 2026-10-12 scope add examples/sensor_std_plain.rs — retarget plaintext reference to the range contract
- 2026-10-12 advance draft → ready
- 2026-10-12 scope approved: examples/sensor_std_dev_demo.rs, examples/sensor_std_plain.rs
- 2026-10-12 scope rm examples/sensor_std_dev_demo.rs
- 2026-10-12 scope add examples/sensor_std_dev/main.rs — demo lives in the multi-file example created by 0002; extend main to N = 8192
