---
id: 0008-ciphertext-levels
status: done
scope: [src/ckks/**, src/rings/**, src/lib.rs, examples/**, tests/**, README.md, src/rings.rs]
created: 2026-10-13
source: user request 2026-10-13: pipeline.rs is a std-dev grab bag; decompose into composable ops like fhe.rs / tfhe-rs
depends: [0007-module-layout]
---
# Task: ciphertext levels

## Plan

- [x] Context chain (fhe.rs): `RnsBasis` links `next: Option<Arc<RnsBasis>>`, built once by `RnsBasis::new`; `next()`, `at_level(n)`; `drop_last` removed. `RnsPoly::rescale` uses `next()` (`rescale_into` removed); `mod_drop_last` -> `mod_switch_to(level)`.
- [x] `Ciphertext`: `logq` field removed, `logq()` derived from the basis; `level()`, `rescale()`, `mod_switch_to(level)`, `match_level(&other)` return `CkksResult` (new `Ring` error wrapping `RnsNttError::InvalidLevel`); `encrypt` drops its `logq` argument.
- [x] `SecretKey::mod_switch_to(level)` replaces `reduce_secret_key`.
- [x] Remove `rescale_ciphertext`, `PolyRescale`, `PolyModSwitch`, the generic `Ciphertext::mod_switch`, `align_level`, `decode_wide` and their re-exports.
- [x] Update callers (pipeline, examples, tests, README); tests: rescales share one basis `Arc`, `mod_switch_to` decrypts to the same slots, level errors.
- [x] Update `specs/overview.md` (ops table, traits, pipeline API, `Arc::ptr_eq` gotcha).

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features`
- [x] `cargo build --examples`
- [x] `cargo bench --no-run`
- [x] `typos`
- [x] `gitleaks dir .`

## Summary

Shared modulus chain: `RnsBasis` links to the next level, `Ciphertext` gains `level`/`logq()`/`rescale`/`mod_switch_to`/`match_level`, `logq` is no longer stored, and `align_level`, `decode_wide`, `rescale_ciphertext` and the unimplemented level traits are gone.

## Review

- One error type instead of two: `RnsNttError::InvalidModDrop` became `InvalidLevel { level, channel_count }`, and `CkksError::Ring` wraps it; no separate `LevelOutOfRange` (it would duplicate the same fact).
- `logp` arithmetic unchanged (subtracts the prime's bit length); `logq()` is the floor-based `total_bits`, so printed `logq` values in examples shift by a few bits. The engine's `logq` equality asserts are gone; `mul_ciphertexts_gadget` asserts equal levels instead.
- Numbers unchanged: `sensor_std_dev` prints the same variance/std errors (5.21e-4, 1.08e-2) as before; `encrypt_add`, `encrypt_mul`, `rotation_demo`, `rotation_stress`, `std_dev_8`, `horner_chain` all pass their own checks.
- Tests: 124 lib (+2), 8 integration (+1 mod-switch test), 4 example, 1 doc. `cargo doc` warnings 6 -> 4 (two fixed in touched files).
- Each chain level clones its prefix of NTT tables (L(L+1)/2 tables); fine at demo sizes, not optimized per the no-optimization rule.

## Context

Terms (HEAAN/fhe.rs): rescale = divide by the last prime and drop it (`logp`, `logq` shrink); mod switch = drop primes without dividing (`logp` stays). Concrete `RnsPoly` methods here; 0011 lifts them into the backend trait.

## Log

- 2026-10-13 created: ciphertext levels
- 2026-10-13 set source user request 2026-10-13: pipeline.rs is a std-dev grab bag; decompose into composable ops like fhe.rs / tfhe-rs
- 2026-10-13 set depends 0007-module-layout
- 2026-10-13 scope add src/ckks/** — Ciphertext level methods, drop traits, update callers
- 2026-10-13 scope add src/rings/** — Ciphertext level methods, drop traits, update callers
- 2026-10-13 scope add src/lib.rs — Ciphertext level methods, drop traits, update callers
- 2026-10-13 scope add examples/** — Ciphertext level methods, drop traits, update callers
- 2026-10-14 scope add tests/** — encrypt signature and level API callers
- 2026-10-14 scope add README.md — encrypt signature and level API callers
- 2026-10-14 advance draft → ready
- 2026-10-14 scope approved: src/ckks/**, src/rings/**, src/lib.rs, examples/**, tests/**, README.md
- 2026-10-14 advance ready → in-progress/implement
- 2026-10-14 advance in-progress/implement → in-progress/verify
- 2026-10-14 advance in-progress/verify → in-progress/review
- 2026-10-14 advance in-progress/review → approval
- 2026-10-14 scope add src/rings.rs — drop PolyRescale/PolyModSwitch re-export
- 2026-10-14 approved: approved by operator
- 2026-10-14 advance approval → done
