---
id: 0012-generic-key-switching
status: draft
scope: [src/ckks/**, src/ckks.rs, src/lib.rs, examples/**, tests/**]
created: 2026-10-14
source: user request 2026-10-13: pluggable ring backends, scheme over a small trait (overview.md Target architecture)
depends: [0011-backend-trait]
---
# Task: generic key switching

## Plan

- [ ] `ckks/keys/key_switching.rs`: one generic `KeySwitchingKey<B>` built on the trait's `decompose`/`gadget_scale` (replaces `generate_gadget_key_matrix` / `gadget_key_switch` in `engine.rs`).
- [ ] `RelinearizationKey` and `RotationKey` wrap it (fhe.rs style); the old gadget-key types and the legacy non-gadget keys move into `tmp/`.
- [ ] `ckks/{ciphertext, ops, evaluator, algorithms}` and `CkksEngine` become generic over the backend; no `impl ... for RnsPoly` left in `ckks/`.
- [ ] Clean `lib.rs` re-exports (drop the workaround block); update examples, tests and `specs/overview.md`.

## Acceptance

- [ ] `cargo fmt --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features`
- [ ] `cargo build --examples`
- [ ] `cargo run --release --example sensor_std_dev`
- [ ] `typos`
- [ ] `gitleaks dir .`

## Log

- 2026-10-14 created: generic key switching
- 2026-10-14 set source user request 2026-10-13: pluggable ring backends, scheme over a small trait (overview.md Target architecture)
- 2026-10-14 set depends 0011-backend-trait
- 2026-10-14 scope add src/ckks/** — generic key switching and scheme, legacy keys to tmp
- 2026-10-14 scope add src/ckks.rs — generic key switching and scheme, legacy keys to tmp
- 2026-10-14 scope add src/lib.rs — generic key switching and scheme, legacy keys to tmp
- 2026-10-14 scope add examples/** — generic key switching and scheme, legacy keys to tmp
- 2026-10-14 scope add tests/** — generic key switching and scheme, legacy keys to tmp
