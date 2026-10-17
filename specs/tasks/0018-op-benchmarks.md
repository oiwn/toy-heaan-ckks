---
id: 0018-op-benchmarks
status: draft
scope: [benches/**, Cargo.toml, README.md]
created: 2026-10-16
source: operator request 2026-10-16
---
# Task: op benchmarks

Criterion benchmarks for every operation, the baseline for SIMD (0019) and wasm (0015) numbers.

## Plan

- [ ] Ring benches at N = 2^10..2^14: NTT/inverse NTT, poly add and mul, automorphism, decompose, rescale.
- [ ] Scheme benches: encode/decode, encrypt/decrypt, add, mul + relin, rescale, rotate (key switch), relin and rotation keygen.
- [ ] Shared setup helper (basis, keys) so each bench times one operation.
- [ ] README bench section; record one baseline run in the task Review.

## Acceptance

- [ ] `cargo fmt --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features`
- [ ] `cargo bench --no-run`
- [ ] `typos`

## Log

- 2026-10-16 created: op benchmarks
- 2026-10-16 set source operator request 2026-10-16
- 2026-10-16 scope add benches/** — per-op criterion benches
- 2026-10-16 scope add Cargo.toml — per-op criterion benches
- 2026-10-16 scope add README.md — per-op criterion benches
