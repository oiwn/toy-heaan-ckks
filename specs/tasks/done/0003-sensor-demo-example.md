---
id: 0003-sensor-demo-example
status: done
scope: [examples/sensor_std_plain.rs, examples/sensor_std_dev/main.rs, examples/sensor_std_dev/std_dev.rs, src/crypto/pipeline.rs, src/encoding/ckks_encoder.rs, examples/std_dev_8.rs]
created: 2026-10-12
source: user request 2026-10-12: revive project, CKKS green + native sensor demo with range normalization
depends: [0002-range-contract-std-dev]
---
# Task: sensor demo example

## Plan

- [x] `std_dev.rs`: `encrypted_std_dev` also returns a stage trace (stage name, RNS channels, logp per op, plus the encrypted variance) so `main` can print stages and check the variance error.
- [x] `main.rs` data: 32 readings (one per 45 min) from a seeded diurnal curve, 15 + 4*sin(...) degC plus ~0.5 degC jitter (std ~2.9); assert `contract.contains_std` before encryption.
- [x] `main.rs` params: N = 8192 (32 readings tiled over all 4096 slots, rotation offsets 1..16), 13 x 61-bit primes, degree 7, `SensorContract::demo()`.
- [x] Print each stage (name, channels, logp, wall time) plus keygen/encrypt/decrypt times, the encrypted std, the f64 std, and both errors; record the timings in `## Findings`.
- [x] Assert std error < 0.1 degC and intermediate variance error < 1e-2 degC^2 (CKKS key-switch noise at N = 8192; 1e-9 was a toy-params number).
- [x] Scalars fill all slots (SEAL-style): `CkksEncoder::encode_constant`; drop `slots` from `mul_plain_scalar`, `mul_plain_exact`, `trivial_constant_ct`, `eval_chebyshev_encrypted`; regression test with partially filled slots.
- [x] Retarget `examples/sensor_std_plain.rs` to the range contract (plaintext reference).
- [x] Tests stay at toy params (N = 16); only `main` runs at N = 8192.

## Context

Runtime risk: `LeveledEvaluator` regenerates the relin key + 5 rotation keys at each of 11 levels, and encode/decode are O(N^2); the first run may take minutes. If it is painfully slow, the fix (keys generated once and modulus-switched down) is in `specs/ideas.md`; keep it out of 0003 unless needed.

## Findings

- Demo run (release, N = 8192): std 2.878781 vs 2.867966 degC (error 1.08e-2, the degree-7 fit); variance 8.225753 vs 8.225232 degC^2 (error 5.2e-4, CKKS noise). Total 24.8 s: sum readings 11.1 s and sum squares 8.0 s (mostly rotation-key generation), square 2.0 s, Chebyshev 3.65 s, constant muls ~22 ms, encode + encrypt 0.19 s, decrypt + decode 0.36 s.
- Bug 1: `sum_slots` sums windows of `slots` consecutive slots (rotations wrap over all N/2). With 32 of 4096 slots filled only slot 0 held the full sum, so the mean and centering were wrong. Fix: the client tiles the readings over all slots; doc corrected.
- Bug 2: scalar helpers broadcast constants to just `slots` slots, so `0.5 @ 1 bit` encoded to coefficients ~0.008 that round to 0 (std came out 5.21 degC). Fix: scalars fill all slots via `encode_constant` (exact constant polynomial, O(N)); this also cut constant muls from ~175 ms to ~22 ms.

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features`
- [x] `cargo run --release --example sensor_std_dev`

## Review

- Library: `CkksEncoder::encode_constant` (SEAL-style scalar encode) and the scalar helpers in `pipeline.rs` lost their misleading `slots` parameter; `sum_slots` docs describe the window semantics. `examples/std_dev_8.rs` updated for the new signature.
- Example: stage trace in `std_dev.rs`; `main.rs` runs the N = 8192 demo with tiled readings and both budgets asserted; `sensor_std_plain.rs` is now the plaintext reference for the contract (envelope sweep, demo pipeline, out-of-envelope behavior) using the library fitter.
- New tests: `encode_constant` exactness and range errors; Chebyshev with partially filled slots (the bug-2 failure mode).
- All checks green: fmt, clippy, 122 lib + 7 integration + 4 example tests, typos, all examples.

## Summary

Native N = 8192 sensor demo: 32 readings, std within 1.1e-2 degC of plaintext in ~25 s; scalars now fill all slots (`encode_constant`), fixing a slot-broadcast bug.

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
- 2026-10-12 scope add examples/sensor_std_dev/std_dev.rs — return a stage trace (name, channels, logp, variance ct) for per-stage printing and the variance-error check
- 2026-10-12 advance ready → in-progress/implement
- 2026-10-12 scope add src/crypto/pipeline.rs — doc fix: sum_slots only fills every slot when slots == DEGREE/2 or the input is tiled (found at N = 8192)
- 2026-10-12 scope add src/encoding/ckks_encoder.rs — encode_constant: SEAL-style scalar encode (exact constant polynomial, all slots)
- 2026-10-12 scope add examples/std_dev_8.rs — mul_plain_scalar drops its slots argument
- 2026-10-12 advance in-progress/implement → in-progress/verify
- 2026-10-12 advance in-progress/verify → in-progress/review
- 2026-10-12 advance in-progress/review → approval
- 2026-10-13 approved: the demo look solid so it work, right? We probably need to move on.
- 2026-10-13 advance approval → done
