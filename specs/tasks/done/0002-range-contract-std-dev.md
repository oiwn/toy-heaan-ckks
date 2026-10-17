---
id: 0002-range-contract-std-dev
status: done
scope: [src/crypto/pipeline.rs, examples/sensor_std_dev/main.rs, examples/sensor_std_dev/contract.rs, examples/sensor_std_dev/std_dev.rs, Cargo.toml]
created: 2026-10-12
source: user request 2026-10-12: revive project, CKKS green + native sensor demo with range normalization
depends: [0004-encoder-wide-staging]
---
# Task: range contract std dev

## Plan

- [x] Library stays generic: replace `fit_sqrt_chebyshev` with `fit_chebyshev(f, lower, upper, degree)` in `pipeline.rs`; drop demo code from the library.
- [x] Remove i64-era workarounds: delete `clamped_scalar_bits`, drop the i64 assert in `mul_plain_exact`, fix stale i64 docs; `decode_wide` delegates to `CkksEncoder::decode`.
- [x] New multi-file example `examples/sensor_std_dev/`: `contract.rs` (`SensorContract`: validation, `normalize`, `variance_interval`, `restore_std`, `contains_std`), `std_dev.rs` (`encrypted_std_dev(ct, slots, &contract, degree, evaluator, rng)`), `main.rs` (toy-params demo run).
- [x] `Cargo.toml`: `[[example]] name = "sensor_std_dev"`, `test = true`, so `cargo test` runs the example's tests.
- [x] Example tests: contract validation, normalize/restore/interval, demo fit |c| < 0.1, end-to-end (N = 16, readings `[12, 14, 17, 20, 21, 19, 16, 13]` degC, degree 7, 13 primes; vs plaintext Chebyshev < 1e-3 degC, vs true std < 0.1 degC).

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features`

## Summary

Sensor std-dev demo moved to `examples/sensor_std_dev/` (`SensorContract` + `encrypted_std_dev`); library keeps generic ops with `fit_chebyshev`; i64-era scalar clamp removed.

## Context

Decisions: `specs/ctx.md` Findings. Target pipeline (inputs already normalized client-side):

```text
ct = Enc([x_n0, ..., x_n{SLOTS-1}])     x_n = (x - lo) / (hi - lo)
  -> mean_ct     = sum_slots(ct) * (1/SLOTS)
  -> centered_ct = ct - mean_ct
  -> variance_ct = sum_slots(centered_ct^2) * (1/SLOTS)
  -> std_n_ct    = eval_cheb(variance_ct)        fit on [(std_min/W)^2, (std_max/W)^2]
client: std = decode(std_n_ct) * (hi - lo)
```

No 1/B level and no encrypted x10 restore, so the pipeline uses 3 + (degree + 1) levels.

## Review

- Revision requested by the operator: demo code (`SensorContract`, `encrypted_std_dev`, sqrt fit) must not live in the library; moved to `examples/sensor_std_dev/`.
- Library: `pipeline.rs` has only generic ops; `fit_sqrt_chebyshev` became `fit_chebyshev(f, lower, upper, degree)`; `clamped_scalar_bits` and the i64 assert in `mul_plain_exact` removed; `decode_wide` delegates to the encoder.
- Example: `contract.rs`, `std_dev.rs`, `main.rs`; `[[example]] test = true` makes `cargo test` run its 4 tests. Demo run (toy params): encrypted std 3.135204 vs plaintext 3.122499 degC, error 1.27e-2 (degree-7 fit error, as predicted).
- All checks green: fmt, clippy, 119 lib + 7 integration + 4 example tests, typos, all 11 examples.

## Log

- 2026-10-12 created: range contract std dev
- 2026-10-12 set source user request 2026-10-12: revive project, CKKS green + native sensor demo with range normalization
- 2026-10-12 set depends 0001-green-baseline
- 2026-10-12 scope add src/crypto/pipeline.rs — encrypted_std_dev and its tests
- 2026-10-12 scope add examples/std_dev_8.rs — may call encrypted_std_dev with the new signature
- 2026-10-12 advance draft → ready
- 2026-10-12 scope approved: src/crypto/pipeline.rs, examples/std_dev_8.rs
- 2026-10-12 set depends 0004-encoder-wide-staging
- 2026-10-12 scope rm examples/std_dev_8.rs
- 2026-10-12 advance ready → in-progress/implement
- 2026-10-12 baseline: .codex/hooks.json@f274ace3f299, .gitignore@800d14263080, _typos.toml@c5a0126833f4, examples/rotation_demo.rs@504a8f668790, examples/sensor_std_plain.rs@9842f3fe4352, examples/std_dev_8.rs@3fc1d2a79141, rust-toolchain.toml@275861857b47, specdev.toml@1ea52da4a446, src/crypto/engine.rs@1a784111cbbc, src/crypto/mod.rs@da20b5b1922b, src/crypto/types.rs@2cf43168b9fe, src/encoding/ckks_encoder.rs@442ab0c719e7, src/encoding/mod.rs@6d00fef3dcad, src/rings/backends/rns_ntt/basis.rs@8bc169c920f4, src/rings/backends/rns_ntt/poly.rs@094f2ffe24e8, tests/integration_mul.rs@8386061ae047
- 2026-10-12 advance in-progress/implement → in-progress/verify
- 2026-10-12 advance in-progress/verify → in-progress/review
- 2026-10-12 advance in-progress/review → approval
- 2026-10-12 advance approval → in-progress/fix (revision)
- 2026-10-12 baseline: .codex/hooks.json@f274ace3f299, .gitignore@800d14263080, _typos.toml@c5a0126833f4, examples/rotation_demo.rs@504a8f668790, examples/sensor_std_plain.rs@9842f3fe4352, examples/std_dev_8.rs@3fc1d2a79141, rust-toolchain.toml@275861857b47, specdev.toml@1ea52da4a446, src/crypto/engine.rs@1a784111cbbc, src/crypto/mod.rs@da20b5b1922b, src/crypto/types.rs@2cf43168b9fe, src/encoding/ckks_encoder.rs@442ab0c719e7, src/encoding/mod.rs@6d00fef3dcad, src/rings/backends/rns_ntt/basis.rs@8bc169c920f4, src/rings/backends/rns_ntt/poly.rs@094f2ffe24e8, tests/integration_mul.rs@8386061ae047
- 2026-10-12 scope add examples/sensor_std_dev/main.rs — operator feedback: split demo code (contract, encrypted std-dev) out of the library into an example
- 2026-10-12 scope add examples/sensor_std_dev/contract.rs — operator feedback: split demo code (contract, encrypted std-dev) out of the library into an example
- 2026-10-12 scope add examples/sensor_std_dev/std_dev.rs — operator feedback: split demo code (contract, encrypted std-dev) out of the library into an example
- 2026-10-12 scope add Cargo.toml — operator feedback: split demo code (contract, encrypted std-dev) out of the library into an example
- 2026-10-12 advance in-progress/fix → in-progress/verify
- 2026-10-12 advance in-progress/verify → in-progress/review
- 2026-10-12 advance in-progress/review → approval
- 2026-10-12 approved: it's ok, so what's left?
- 2026-10-12 advance approval → done
