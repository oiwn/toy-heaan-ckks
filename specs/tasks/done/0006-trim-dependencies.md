---
id: 0006-trim-dependencies
status: done
scope: [Cargo.toml]
created: 2026-10-13
source: user request 2026-10-13: publication cleanup, remove ancient code and golden vectors, foo.rs module layout
depends: [0005-remove-dead-code]
---
# Task: trim dependencies

## Plan

- [x] Remove unused `[dependencies]`: `dhat`, `serde`, `serde_with`, `serde_json`. Keep `tracing` (planned logging) and `crypto-bigint` (planned use).
- [x] Move `tracing-subscriber` to `[dev-dependencies]`: only examples install a subscriber.
- [x] Remove unused `[dev-dependencies]`: `axum`, `tower`, `base64`, `reqwest`, `tokio`, `fjall`, `paste`, `proptest`.
- [x] Remove the `dhat-heap` / `dhat-ad-hoc` features.
- [x] Re-grep each removed crate before deleting it; `Cargo.lock` follows. serde comes back with M3 when needed.

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features`
- [x] `cargo build --examples`
- [x] `cargo bench --no-run`
- [x] `cargo tree -e normal --depth 1` lists no removed crate
- [x] `typos`
- [x] `gitleaks dir .`

## Summary

Trimmed `Cargo.toml`: removed dhat, serde stack and 8 unused dev-dependencies plus the dhat features; `tracing-subscriber` moved to dev; `tracing` and `crypto-bigint` kept for planned work.

## Context

`comfy-table` stays while `src/table.rs` exists (used by 3 examples); the audit decides whether it moves to the examples and becomes a dev-dependency.

## Review

- Only `Cargo.toml` (and `Cargo.lock`) changed; matches the Plan. `tracing` and `crypto-bigint` carry a one-line comment on why they stay unused for now.
- Runtime deps now: comfy-table, crypto-bigint, num-bigint, num-complex, num-traits, rand, rand_chacha, rand_distr, rustfft, thiserror, tracing.
- All tests (122 lib + 7 integration + 4 example + 1 doc) and benches build; CI's `--all-features` / `--no-default-features` stay valid with no features defined.
- Noticed, not fixed: `benches/basic.rs` has no `[[bench]] harness = false` entry, unlike `primes`; audit item.

## Log

- 2026-10-13 created: trim dependencies
- 2026-10-13 set source user request 2026-10-13: publication cleanup, remove ancient code and golden vectors, foo.rs module layout
- 2026-10-13 set depends 0005-remove-dead-code
- 2026-10-13 scope add Cargo.toml — remove unused deps and features
- 2026-10-14 advance draft → ready
- 2026-10-14 scope approved: Cargo.toml
- 2026-10-14 advance ready → in-progress/implement
- 2026-10-14 baseline: .gitignore@1259d5434660, README.md@2230332f10de, data/golden_context_light.json@deleted, data/golden_context_small.json@deleted, data/golden_encode_light.json@deleted, scripts/fit_sqrt_chebyshev.sage@deleted, scripts/fit_sqrt_chebyshev.sage.py@deleted, scripts/reference_ckks_encode.py@deleted, src/crypto/mod.rs@fff8c10854fc, src/crypto/operations.rs@deleted, src/encoding/fft.rs@deleted, src/encoding/mod.rs@58cdb459cdb2, src/encoding/textbook.rs@deleted, src/rings/README.md@deleted, src/rings/backends/rns_ntt/poly.rs@fb1b44d207de, src/rings/mod.rs@16a83557c10f
- 2026-10-14 advance in-progress/implement → in-progress/verify
- 2026-10-14 advance in-progress/verify → in-progress/review
- 2026-10-14 advance in-progress/review → approval
- 2026-10-14 approved: let's close 0006 tests passing
- 2026-10-14 advance approval → done
