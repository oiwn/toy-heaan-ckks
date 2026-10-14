---
id: 0008-ciphertext-levels
status: draft
scope: [src/ckks/**, src/rings/**, src/lib.rs, examples/**]
created: 2026-10-13
source: user request 2026-10-13: pipeline.rs is a std-dev grab bag; decompose into composable ops like fhe.rs / tfhe-rs
depends: [0007-module-layout]
---
# Task: ciphertext levels

## Plan

- [ ] Context chain (fhe.rs): `RnsBasis` gets `next() -> Option<&Arc<RnsBasis>>`, the whole chain built once from the prime list; `rescale_into`/`drop_last` callers use `next()` instead of minting a fresh `Arc`.
- [ ] `RnsPoly::mod_switch_to(basis)`: keep the first channels of a basis further down the same chain (replaces private `reduce_poly_to_basis`).
- [ ] `Ciphertext::rescale()` (moves `CkksEngine::rescale_ciphertext`) and `Ciphertext::mod_switch_to(target: &Ciphertext)` (drops primes, computes `logq` itself; `logp` unchanged).
- [ ] `SecretKey::mod_switch_to(basis)` replaces `reduce_secret_key`.
- [ ] Remove `PolyRescale`, `PolyModSwitch`, the generic `Ciphertext::mod_switch`, `align_level`, `decode_wide` (callers use `CkksEncoder::decode`) and their re-exports.
- [ ] Update callers (pipeline, examples); tests: two rescales of different ciphertexts land on the same basis `Arc`, rescale divides by the last prime, `mod_switch_to` keeps decrypting to the same slots; update the `Arc::ptr_eq` gotcha in `specs/overview.md`.

## Acceptance

- [ ] `cargo fmt --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features`
- [ ] `cargo build --examples`
- [ ] `typos`
- [ ] `gitleaks dir .`

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
