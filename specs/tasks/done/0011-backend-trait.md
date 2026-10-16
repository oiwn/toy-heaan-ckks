---
id: 0011-backend-trait
status: done
scope: [src/rings/**, src/rings.rs, src/ckks/**, examples/**, tests/**, benches/**, src/lib.rs, README.md]
created: 2026-10-14
source: user request 2026-10-13: pluggable ring backends, scheme over a small trait (overview.md Target architecture)
depends: [0010-chebyshev-module]
---
# Task: backend trait

## Plan

- [x] One trait `PolyRing` in `src/rings/traits.rs` (absorbs `PolySampler`, `PolyAutomorphism`): arithmetic (+ `SubAssign`), sampling (`sample_ternary`), automorphism (`rotate_slots` provided), modulus chain (`level`, `modulus_bits`, `rescale_bits`, `rescale`, `mod_switch_to`), encoding boundary (`from_rounded_f64`, `to_centered_f64`), gadget (`gadget_len`, `decompose`, `gadget_term`); backend-agnostic `RingError` in `src/rings/errors.rs`.
- [x] `RnsPoly` keeps NTT form at rest (OpenFHE/SEAL style) behind a private flag; coefficient form only inside encode/sample, decode, decompose, rescale (last channel) and automorphism. Domain methods, `channels()` and `from_channels` go private.
- [x] `CkksEncoder` generic over the trait; `Ciphertext` level methods generic; engine gadget keygen / key switch built on `decompose` / `gadget_term` (flat key vectors, same draw order).
- [x] Update callers, tests (gadget identity, f64 boundary round trip), `specs/overview.md` (NTT-at-rest gotcha, OpenFHE reference), README/AGENTS lines; numbers bit-identical.

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features`
- [x] `cargo build --examples`
- [x] `cargo bench --no-run`
- [x] `cargo run --release --example sensor_std_dev`
- [x] `typos`
- [x] `gitleaks dir .`

## Summary

One backend trait `PolyRing` (arithmetic, sampling, automorphism, modulus chain, encoding boundary, gadget) replaces the three ring traits; `RnsPoly` rests in NTT form behind a private flag, and the encoder, `Ciphertext` level methods and gadget key switching use only the trait. Output bit-identical, sensor pipeline 24 s -> 9 s.

## Review

- Bit-identical: all 10 examples diffed against pre-change output; only timings differ (sensor 24.1 s -> 9.3 s, keys and ciphertexts no longer round-trip through NTT per mul).
- Deviation: `EncodingError` removed instead of mapped; the encoder returns `RingResult` (it only duplicated `RingError::CoefficientTooLarge`, no users outside the encoder).
- `from_channels`, `channels()` and `mul_assign_naive` removed rather than made private (no callers left), with the two `RnsNttError` variants only they used; the naive product is a test helper. Inherent `RnsPoly::zero`/`from_coeffs(Arc)` gave way to the trait versions (`&Context`).
- The encoder's backend `P` cannot be inferred from the context alone: 3 calls in `examples/encode_decode.rs` (decode-only, no engine) got a `Plaintext<RnsPoly<N>, N>` annotation.
- Tests stay at 128 lib: flag-based tests replaced by NTT-at-rest, mixed-form, chain bits, sub, f64 boundary, gadget identity and rotate-vs-automorphism tests; the BigInt residue test moved into the backend. Leak sweep over `src/ckks`, examples, tests, benches is empty.

## Context

The trait must admit a second backend (M2.6 reference on `crypto-bigint`, one big modulus); nothing RNS-shaped in signatures. Degree stays `const N`. Decided 2026-10-16: NTT at rest with a private flag (convert only where coefficients are needed), one task, keep the name `PolyRing`. Ops, evaluator and algorithms stay on `RnsPoly` until 0012.

## Log

- 2026-10-14 created: backend trait
- 2026-10-14 set source user request 2026-10-13: pluggable ring backends, scheme over a small trait (overview.md Target architecture)
- 2026-10-14 set depends 0010-chebyshev-module
- 2026-10-14 scope add src/rings/** — replace ring traits, hide NTT domain, update callers
- 2026-10-14 scope add src/rings.rs — replace ring traits, hide NTT domain, update callers
- 2026-10-14 scope add src/ckks/** — replace ring traits, hide NTT domain, update callers
- 2026-10-14 scope add examples/** — replace ring traits, hide NTT domain, update callers
- 2026-10-14 scope add tests/** — replace ring traits, hide NTT domain, update callers
- 2026-10-14 scope add benches/** — replace ring traits, hide NTT domain, update callers
- 2026-10-16 scope add src/lib.rs — re-exports and layout line name the old traits
- 2026-10-16 scope add README.md — re-exports and layout line name the old traits
- 2026-10-16 advance draft → ready
- 2026-10-16 scope approved: src/rings/**, src/rings.rs, src/ckks/**, examples/**, tests/**, benches/**, src/lib.rs, README.md
- 2026-10-16 advance ready → in-progress/implement
- 2026-10-16 advance in-progress/implement → in-progress/verify
- 2026-10-16 advance in-progress/verify → in-progress/review
- 2026-10-16 advance in-progress/review → approval
- 2026-10-16 approved: approved by operator
- 2026-10-16 advance approval → done
