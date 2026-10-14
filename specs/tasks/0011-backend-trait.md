---
id: 0011-backend-trait
status: draft
scope: [src/rings/**, src/rings.rs, src/ckks/**, examples/**, tests/**, benches/**]
created: 2026-10-14
source: user request 2026-10-13: pluggable ring backends, scheme over a small trait (overview.md Target architecture)
depends: [0010-chebyshev-module]
---
# Task: backend trait

## Plan

- [ ] One backend trait in `src/rings/traits.rs` with the six groups from `specs/overview.md` Target architecture (arithmetic, sampling, automorphism, modulus chain, encoding boundary, gadget); replaces `PolyRing`, `PolySampler`, `PolyAutomorphism`.
- [ ] Implement it for `RnsPoly`; the 0008 level methods and the encoder's residue trick (`from_rounded_f64` / `to_centered_f64`) move behind it.
- [ ] NTT vs coefficient domain becomes private to `RnsPoly` (lazy conversion); no domain flag or `to_*_domain` call outside `rings/backends/`.
- [ ] `CkksEncoder` uses only the trait boundary (no `moduli()`/`channels()`).
- [ ] Update callers and the domain-flag gotcha in `specs/overview.md`; behavior unchanged.

## Acceptance

- [ ] `cargo fmt --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features`
- [ ] `cargo build --examples`
- [ ] `cargo run --release --example sensor_std_dev`
- [ ] `typos`
- [ ] `gitleaks dir .`

## Context

The trait must admit a second backend (M2.6 reference on `crypto-bigint`, one big modulus); nothing RNS-shaped in signatures. Degree stays `const N`.

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
