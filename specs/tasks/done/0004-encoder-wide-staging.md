---
id: 0004-encoder-wide-staging
status: done
scope: [src/encoding/ckks_encoder.rs, src/encoding/mod.rs, _typos.toml]
created: 2026-10-12
source: user request 2026-10-12: encoder i64 staging caps |x| < 4 and saturates silently; widen like SEAL, keep f64 FFT
---
# Task: encoder wide staging

## Plan

- [x] `src/encoding/mod.rs`: minimal `EncodingError` (thiserror) with `CoefficientTooLarge { value, max_bits }`.
- [x] Private `rounded_to_residue(x, prime)`: i128 path below 2^63; mantissa * 2^exp mod p above; negate for x < 0.
- [x] `try_encode` / `try_encode_complex` return `Result`; reject non-finite or `|coeff| >= 2^(total_bits - 1)`; build channels via `RnsPoly::from_channels`.
- [x] `encode` / `encode_complex` wrap `try_*` with `expect` and a `# Panics` doc (call sites unchanged).
- [x] `decode_complex` reconstructs with `RnsBasis::reconstruct_centered_coeff_f64` instead of `to_coeffs()`.
- [x] Tests: large round-trips (10, -250, 1e6 at scale 61); residues bit-exact vs BigUint (powers of two, negatives); Q/2 boundary and NaN/inf errors; `encode` panic message.

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features -- --skip encrypted_std_dev_end_to_end`

## Summary

Encoder stages coefficients in f64 and reduces them straight into RNS (no i64 cap, no silent saturation); `try_encode` errors past Q/2; decode reconstructs in f64.

## Context

The encoder rounds f64 FFT coefficients into i64 (`as i64` saturates silently), capping |x| < 2^(63 - scale_bits), ~4 at scale 61. CKKS only needs |coeff| < Q/2. SEAL keeps a `complex<double>` FFT and splits large coefficients into 64-bit words before RNS decomposition, throwing "encoded values are too large"; we keep the f64 FFT and do the same via the f64 mantissa/exponent. Out of scope (task 0002): removing `clamped_scalar_bits` and the i64 asserts in `pipeline.rs`, folding `decode_wide` into the encoder. Also record whether this alone turns `encrypted_std_dev_end_to_end` green (confirms the saturation diagnosis).

## Findings

- This change alone turns `encrypted_std_dev_end_to_end` green: the old failure was i64 saturation of c0 ~ 6.5, not a q0 wrap.
- Precision is ~2^-52 relative to the largest coefficient (a 1e6 slot costs a 1e-4 slot ~1e-9 absolute); documented on `try_encode`.

## Review

- Diff stays in scope (`src/encoding/ckks_encoder.rs`, `src/encoding/mod.rs`); public `encode`/`encode_complex` signatures unchanged, so no call sites moved.
- `rounded_to_residue` checked bit-exact against BigInt for powers of two, 53-bit mantissas up to 2^150, and negatives.
- Full suite green including the previously red test (120 unit, 7 integration); all 10 examples run.
- `typos` hook: added `_typos.toml` excluding `specs/tasks/done/` (closed logs quote operator approvals verbatim).

## Log

- 2026-10-12 created: encoder wide staging
- 2026-10-12 set source user request 2026-10-12: encoder i64 staging caps |x| < 4 and saturates silently; widen like SEAL, keep f64 FFT
- 2026-10-12 scope add src/encoding/ckks_encoder.rs — f64 -> RNS staging, try_encode, wide decode, EncodingError
- 2026-10-12 scope add src/encoding/mod.rs — f64 -> RNS staging, try_encode, wide decode, EncodingError
- 2026-10-12 advance draft → ready
- 2026-10-12 scope approved: src/encoding/ckks_encoder.rs, src/encoding/mod.rs
- 2026-10-12 advance ready → in-progress/implement
- 2026-10-12 baseline: .codex/hooks.json@f274ace3f299, .gitignore@800d14263080, examples/rotation_demo.rs@504a8f668790, examples/sensor_std_plain.rs@9842f3fe4352, examples/std_dev_8.rs@3fc1d2a79141, rust-toolchain.toml@275861857b47, specdev.toml@1ea52da4a446, src/crypto/engine.rs@1a784111cbbc, src/crypto/mod.rs@da20b5b1922b, src/crypto/pipeline.rs@807165b27582, src/crypto/types.rs@2cf43168b9fe, src/rings/backends/rns_ntt/basis.rs@8bc169c920f4, src/rings/backends/rns_ntt/poly.rs@094f2ffe24e8, tests/integration_mul.rs@8386061ae047
- 2026-10-12 advance in-progress/implement → in-progress/verify
- 2026-10-12 advance in-progress/verify → in-progress/review
- 2026-10-12 advance in-progress/review → approval
- 2026-10-12 scope add _typos.toml — typos hook: skip closed task logs that quote the operator verbatim
- 2026-10-12 approved: let's close
- 2026-10-12 advance approval → done
