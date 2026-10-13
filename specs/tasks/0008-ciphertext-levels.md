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

- [ ] `RnsPoly::mod_switch_to(basis)`: keep the first channels of a prefix basis (replaces private `reduce_poly_to_basis`).
- [ ] `Ciphertext::rescale()` (moves `CkksEngine::rescale_ciphertext`; shared basis `Arc`) and `Ciphertext::mod_switch_to(target: &Ciphertext)` (drops primes, computes `logq` itself; `logp` unchanged).
- [ ] `SecretKey::mod_switch_to(basis)` replaces `reduce_secret_key`.
- [ ] Remove `PolyRescale`, `PolyModSwitch`, the generic `Ciphertext::mod_switch`, `align_level`, `decode_wide` (callers use `CkksEncoder::decode`) and their re-exports.
- [ ] Update callers (pipeline, examples); tests: rescale divides by the last prime, `mod_switch_to` keeps decrypting to the same slots.

## Acceptance

- [ ] `cargo fmt --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features`
- [ ] `cargo build --examples`
- [ ] `typos`
- [ ] `gitleaks dir .`

## Context

On hold: revised after the backend design in `specs/ctx.md` is agreed (rescale/mod switch become backend trait ops). Terms (HEAAN/fhe.rs): rescale = divide by the last prime and drop it (`logp`, `logq` shrink); mod switch = drop primes without dividing (`logp` stays). Exact RNS methods, no trait while there is one backend.

## Log

- 2026-10-13 created: ciphertext levels
- 2026-10-13 set source user request 2026-10-13: pipeline.rs is a std-dev grab bag; decompose into composable ops like fhe.rs / tfhe-rs
- 2026-10-13 set depends 0007-module-layout
- 2026-10-13 scope add src/ckks/** — Ciphertext level methods, drop traits, update callers
- 2026-10-13 scope add src/rings/** — Ciphertext level methods, drop traits, update callers
- 2026-10-13 scope add src/lib.rs — Ciphertext level methods, drop traits, update callers
- 2026-10-13 scope add examples/** — Ciphertext level methods, drop traits, update callers
