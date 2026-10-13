---
id: 0005-remove-dead-code
status: draft
scope: [data/**, scripts/**, src/encoding/**, src/crypto/**, src/rings/**, README.md]
created: 2026-10-13
source: user request 2026-10-13: publication cleanup, remove ancient code and golden vectors, foo.rs module layout
---
# Task: remove dead code

## Plan

Files are moved into `tmp/` keeping their relative path (AGENTS.md rule); the operator deletes them for real. Inline code is removed by normal edits.

- [ ] Move golden vectors and scripts into `tmp/`: `data/` (3 JSON files + `ckks.py`, nothing reads them) and `scripts/` (sage sqrt fit, python encode reference; superseded by `fit_chebyshev` and the encoder tests).
- [ ] Move uncompiled encoders `src/encoding/fft.rs`, `src/encoding/textbook.rs` into `tmp/`; remove the commented-out code in `src/encoding/mod.rs` and `src/rings/mod.rs`.
- [ ] Move `src/crypto/operations.rs` (FIXME: duplicate of `CkksEngine`) into `tmp/`; remove its re-exports in `src/crypto/mod.rs`.
- [ ] Remove `RnsPoly::channels_mut` (no users); move `src/rings/README.md` (mojibake, stale) into `tmp/`.
- [ ] Fix references to removed files in README.md, AGENTS.md (Key File Paths, `specs/current_task.md` link) and specs.

## Acceptance

- [ ] `cargo fmt --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features`
- [ ] `cargo build --examples`
- [ ] `typos`
- [ ] `gitleaks dir .`

## Context

Live but questionable code (legacy keys, `table.rs`, examples) is left for the audit step in `specs/ctx.md`; `PolyRescale`/`PolyModSwitch` go to the rescale/mod-switch discussion (task 0008).

## Log

- 2026-10-13 created: remove dead code
- 2026-10-13 set source user request 2026-10-13: publication cleanup, remove ancient code and golden vectors, foo.rs module layout
- 2026-10-13 scope add data/** — delete dead code, golden vectors, scripts; fix references
- 2026-10-13 scope add scripts/** — delete dead code, golden vectors, scripts; fix references
- 2026-10-13 scope add src/encoding/** — delete dead code, golden vectors, scripts; fix references
- 2026-10-13 scope add src/crypto/** — delete dead code, golden vectors, scripts; fix references
- 2026-10-13 scope add src/rings/** — delete dead code, golden vectors, scripts; fix references
- 2026-10-13 scope add README.md — delete dead code, golden vectors, scripts; fix references
