---
id: 0039-scale-metadata-migration
status: done
scope: [src/ckks/encoding/ckks_encoder.rs, src/ckks/*.rs, src/ckks/algorithms/**/*.rs, examples/**/*.rs, tests/*.rs]
created: 2026-11-02
source: Operator-approved bootstrap split plan 2026-11-02 (from 0037)
depends: [0038-fallible-ciphertext-arithmetic]
---
# Task: scale metadata migration

## Context

Phase one of the two-phase scale migration (operator decision 2026-11-02): a pure type migration. Rescale still divides by the nominal power of two `2^rescale_bits`, which is exact in f64, so every result stays bit-identical. Switching to the actual prime divisor happens only in 0040, together with the graph executor.

Rules: touch only scoped files (extend with `specdev task scope` first). If the same check fails twice, stop, record the failure and a hypothesis in Findings and ask the operator. Never weaken tolerances, disable noise or add `#[ignore]`.

## Plan

- [x] `src/ckks/types.rs`: `Ciphertext.logp: u32` -> `scale: Scale`; `Plaintext.scale_bits` -> `scale: Scale`. `rescale` divides by `Scale::power_of_two(rescale_bits)` (doc: nominal until 0040); `mod_switch_to` keeps the scale.
- [x] Encoder: `CkksEncoder::new(Scale)` only (operator decision: no bit-count convenience constructors), `scale_factor()` returns `scale.value()`. `CkksParams` and builder take `Scale`.
- [x] `src/ckks/ops.rs`: `mul_const(value, Scale)`; `mul`/`mul_plain` use `Scale::checked_mul` (documented panic on f64 overflow); add/sub compare with `Scale::is_compatible`. `ScaleMismatch` carries log2 values.
- [x] Algorithms: `prepare_input` accepts only `2^k` with integer k in 1..=44 (as today); recurrence requires an exact power-of-two output lift. `BootstrapParameters.scale_bits` stays `u32` until 0041.
- [x] Reuse `src/ckks/scale.rs` (`power_of_two`, `checked_mul`, `checked_div`, `is_compatible`, `log2`); no new scale helpers elsewhere. Diagnostics print `scale.log2()`.
- [x] Migrate every remaining `logp`/`scale_bits` use in src, examples and tests (`grep -rn "logp\|scale_bits"` must return only docs that explain the old name, or nothing). `src/rings/` uses are docs/comments only; extend scope first if code changes there.
- [x] Add encode/decode and add/mul tests with a non-power-of-two scale.
- [x] Confirm bootstrap demo max error is identical to the baseline at both degrees.

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features`
- [x] `cargo build --examples`
- [x] `cargo run --release --example bootstrap_demo`
- [x] `typos`
- [x] `specdev check`
- [x] `git diff --check`

## Review

- Scope: changes stay in `src/ckks/*.rs`, the encoder, `src/ckks/algorithms/**`, `examples/**` and `tests/*.rs`; `src/rings/` is untouched.
- One helper was added to `Scale` itself: `exact_log2()` returns `k` when the scale is exactly `2^k`. Recurrence and bootstrap `prepare_input` both use it for their exact power-of-two rules, so neither algorithm file has its own scale helper.
- Behavior is bit-identical. The bootstrap_demo max errors equal the baseline at N=16/64 for both cycles (8.059113425635853e-6, 1.7666600415874825e-5, 7.9910586165085e-5, 2.0767430746117748e-4). `recurrence_bootstrap_demo` completes 32/32 rounds.
- New tests:
  - non-power-of-two encode/decode (encoder);
  - add, ct x ct mul and ct x pt mul at 1.5 * 2^30 (ops);
  - `exact_log2` (scale);
  - a bootstrap wrong-scale case at 1.5 * 2^44.
  - Test count: 147 lib (was 144), 10 integration, 4 example, 1 doctest.
- Remaining `scale_bits` hits:
  - `BootstrapParameters.scale_bits`, which 0041 replaces;
  - the `rescale_bits` ring method;
  - the `encode_decode` table label, kept by decision.
- Builder setter finding is recorded in `specs/cleanup.md`, not fixed here.

## Summary

Ciphertexts, plaintexts, the encoder and engine params carry a `Scale` instead of bit counts; rescale still divides by the nominal power of two, so results are bit-identical.

## Log

- 2026-11-02 created: scale metadata migration
- 2026-11-02 set source Operator-approved bootstrap split plan 2026-11-02 (from 0037)
- 2026-11-02 set depends 0038-fallible-ciphertext-arithmetic
- 2026-11-02 scope add src/ckks/types.rs — Replace logp/scale_bits with Scale across library, examples and tests
- 2026-11-02 scope add src/ckks/ops.rs — Replace logp/scale_bits with Scale across library, examples and tests
- 2026-11-02 scope add src/ckks/encoding/ckks_encoder.rs — Replace logp/scale_bits with Scale across library, examples and tests
- 2026-11-02 scope add src/ckks/*.rs — Replace logp/scale_bits with Scale across library, examples and tests
- 2026-11-02 scope add src/ckks/algorithms/**/*.rs — Replace logp/scale_bits with Scale across library, examples and tests
- 2026-11-02 scope add examples/**/*.rs — Replace logp/scale_bits with Scale across library, examples and tests
- 2026-11-02 scope add tests/*.rs — Replace logp/scale_bits with Scale across library, examples and tests
- 2026-11-02 scope rm src/ckks/types.rs — covered by src/ckks/*.rs
- 2026-11-02 scope rm src/ckks/ops.rs — covered by src/ckks/*.rs
- 2026-11-02 draft approved: approved by operator
- 2026-11-02 advance draft → ready
- 2026-11-02 scope approved: src/ckks/encoding/ckks_encoder.rs, src/ckks/*.rs, src/ckks/algorithms/**/*.rs, examples/**/*.rs, tests/*.rs
- 2026-11-02 advance ready → in-progress/implement
- 2026-11-02 baseline: src/ckks.rs@2615ecc20a4b
- 2026-11-02 advance in-progress/implement → in-progress/verify
- 2026-11-02 advance in-progress/verify → in-progress/review
- 2026-11-02 advance in-progress/review → approval
- 2026-11-02 approved: approved by operator
- 2026-11-02 advance approval → done
