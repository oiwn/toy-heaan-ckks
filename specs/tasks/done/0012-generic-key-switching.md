---
id: 0012-generic-key-switching
status: done
scope: [src/ckks/**, src/ckks.rs, src/lib.rs, examples/**, tests/**, README.md]
created: 2026-10-14
source: user request 2026-10-13: pluggable ring backends, scheme over a small trait (overview.md Target architecture)
depends: [0011-backend-trait]
---
# Task: generic key switching

## Plan

- [x] `ckks/keys/key_switching.rs`: generic `KeySwitchingKey<P, N> { a, b }` with `generate(sk, s_old, error_std, rng)` and `switch(d)`, built on `gadget_len`/`decompose`/`gadget_term` (replaces `generate_gadget_key_vectors` / `gadget_key_switch`); `DIGIT_BITS` moves here.
- [x] `RelinearizationKey { ksk }` (s^2) and `RotationKey { ksk, rotation }` wrap it (fhe.rs / OpenFHE `EvalKey` style); legacy single-pair keys (to `tmp/`), `RnsGadget*` types, `mul_ciphertexts` and `generate_relinearization_key` go; engine methods become `generate_relin_key` / `generate_rotation_key`.
- [x] Ops, `LeveledEvaluator`, algorithms, engine and builder generic over `P: PolyRing<N>`; no `impl ... RnsPoly` left in `ckks/` outside tests.
- [x] Drop the `lib.rs` workaround re-exports; update examples (`keys` shows a gadget-key relation), tests, `specs/overview.md`, AGENTS; outputs bit-identical except `keys`.

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

One generic `KeySwitchingKey<P>` (gadget encryption of an old key, `switch(d)`) now backs `RelinearizationKey` and `RotationKey`; the legacy single-pair keys, `RnsGadget*` types and `mul_ciphertexts` are gone, and ops, evaluator, algorithms and engine are generic over `PolyRing`.

## Review

- Bit-identical: 9 examples match the pre-0011 output except timing padding; `keys` changed as planned (relation check now covers every gadget term, max error 6).
- `LeveledEvaluator` builds keys with `RelinearizationKey::generate` / `RotationKey::generate` directly (no throwaway engine); same key error `sqrt(error_variance)` and draw order.
- Legacy `relin_key.rs` / `rotation_key.rs` were copied whole to `tmp/src/ckks/keys/` before the rewrite; their `*Params` / `*Error` types are gone (the `lib.rs` block that re-exported them too).
- Tests: 126 lib (legacy key tests replaced by 2 key-switch residual tests and 2 per-term relation tests), 8 integration, 4 example, 1 doc. Sweeps: no `RnsGadget`/`mul_ciphertexts`/legacy params anywhere; `RnsPoly` in `ckks/` only in tests and one doc example.
- README added to scope for its layout line.

## Log

- 2026-10-14 created: generic key switching
- 2026-10-14 set source user request 2026-10-13: pluggable ring backends, scheme over a small trait (overview.md Target architecture)
- 2026-10-14 set depends 0011-backend-trait
- 2026-10-14 scope add src/ckks/** — generic key switching and scheme, legacy keys to tmp
- 2026-10-14 scope add src/ckks.rs — generic key switching and scheme, legacy keys to tmp
- 2026-10-14 scope add src/lib.rs — generic key switching and scheme, legacy keys to tmp
- 2026-10-14 scope add examples/** — generic key switching and scheme, legacy keys to tmp
- 2026-10-14 scope add tests/** — generic key switching and scheme, legacy keys to tmp
- 2026-10-16 advance draft → ready
- 2026-10-16 scope approved: src/ckks/**, src/ckks.rs, src/lib.rs, examples/**, tests/**
- 2026-10-16 advance ready → in-progress/implement
- 2026-10-16 baseline: README.md@d1e119da32f6, src/rings.rs@c62ccf250a83, src/rings/backends/rns_ntt/errors.rs@031bfad22ba0, src/rings/backends/rns_ntt/poly.rs@c779253984c0, src/rings/errors.rs@3202c7295cf7, src/rings/traits.rs@4272e830a012
- 2026-10-16 scope add README.md — layout line names the key types
- 2026-10-16 advance in-progress/implement → in-progress/verify
- 2026-10-16 advance in-progress/verify → in-progress/review
- 2026-10-16 advance in-progress/review → approval
- 2026-10-16 approved: approved by operator
- 2026-10-16 advance approval → done
