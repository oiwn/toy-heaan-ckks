---
id: 0006-trim-dependencies
status: draft
scope: [Cargo.toml]
created: 2026-10-13
source: user request 2026-10-13: publication cleanup, remove ancient code and golden vectors, foo.rs module layout
depends: [0005-remove-dead-code]
---
# Task: trim dependencies

## Plan

- [ ] Remove unused `[dependencies]`: `dhat`, `serde`, `serde_with`, `serde_json`. Keep `tracing` (planned logging) and `crypto-bigint` (planned use).
- [ ] Move `tracing-subscriber` to `[dev-dependencies]`: only examples install a subscriber.
- [ ] Remove unused `[dev-dependencies]`: `axum`, `tower`, `base64`, `reqwest`, `tokio`, `fjall`, `paste`, `proptest`.
- [ ] Remove the `dhat-heap` / `dhat-ad-hoc` features.
- [ ] Re-grep each removed crate before deleting it; `Cargo.lock` follows. serde comes back with M3 when needed.

## Acceptance

- [ ] `cargo fmt --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features`
- [ ] `cargo build --examples`
- [ ] `typos`
- [ ] `gitleaks dir .`

## Context

`comfy-table` stays while `src/table.rs` exists (used by 3 examples); the audit decides whether it moves to the examples and becomes a dev-dependency.

## Log

- 2026-10-13 created: trim dependencies
- 2026-10-13 set source user request 2026-10-13: publication cleanup, remove ancient code and golden vectors, foo.rs module layout
- 2026-10-13 set depends 0005-remove-dead-code
- 2026-10-13 scope add Cargo.toml — remove unused deps and features
