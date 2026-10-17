---
id: 0005-remove-dead-code
status: done
scope: [data/**, scripts/**, src/encoding/**, src/crypto/**, src/rings/**, README.md, .gitignore]
created: 2026-10-13
source: user request 2026-10-13: publication cleanup, remove ancient code and golden vectors, foo.rs module layout
---
# Task: remove dead code

## Plan

Files are moved into `tmp/` keeping their relative path (AGENTS.md rule); the operator deletes them for real. Inline code is removed by normal edits.

- [x] Move golden vectors and scripts into `tmp/`: `data/` (3 JSON files + `ckks.py`, nothing reads them) and `scripts/` (sage sqrt fit, python encode reference; superseded by `fit_chebyshev` and the encoder tests).
- [x] Move uncompiled encoders `src/encoding/fft.rs`, `src/encoding/textbook.rs` into `tmp/`; remove the commented-out code in `src/encoding/mod.rs` and `src/rings/mod.rs`.
- [x] Move `src/crypto/operations.rs` (FIXME: duplicate of `CkksEngine`) into `tmp/`; remove its re-exports in `src/crypto/mod.rs`.
- [x] Remove `RnsPoly::channels_mut` (no users); move `src/rings/README.md` (mojibake, stale) into `tmp/`.
- [x] Fix references to removed files in README.md, AGENTS.md (Key File Paths, `specs/current_task.md` link) and specs.

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features`
- [x] `cargo build --examples`
- [x] `typos`
- [x] `gitleaks dir .`

## Summary

Moved dead code and stale artifacts (golden vectors, scripts, uncompiled encoders, `operations.rs`) into `tmp/`; removed commented-out code and `channels_mut`; fixed stale references.

## Context

Live but questionable code (legacy keys, `table.rs`, examples) is left for the audit step in `specs/ctx.md`; `PolyRescale`/`PolyModSwitch` go to 0008/0011.

## Review

- Diff matches the Plan and stays in scope; also touched AGENTS.md (always allowed) and `specs/ideas.md` (dropped the sage-script idea).
- No behavior change: 122 lib + 7 integration + 4 example + 1 doc tests pass; `rg` finds no references to the moved items outside `tmp/`, `specs/tasks/`, `specs/archive/`.
- Note: `data/ckks.py` was gitignored (untracked); the rest of `data/` and `scripts/` show as deletions in git.

## Log

- 2026-10-13 created: remove dead code
- 2026-10-13 set source user request 2026-10-13: publication cleanup, remove ancient code and golden vectors, foo.rs module layout
- 2026-10-13 scope add data/** — delete dead code, golden vectors, scripts; fix references
- 2026-10-13 scope add scripts/** — delete dead code, golden vectors, scripts; fix references
- 2026-10-13 scope add src/encoding/** — delete dead code, golden vectors, scripts; fix references
- 2026-10-13 scope add src/crypto/** — delete dead code, golden vectors, scripts; fix references
- 2026-10-13 scope add src/rings/** — delete dead code, golden vectors, scripts; fix references
- 2026-10-13 scope add README.md — delete dead code, golden vectors, scripts; fix references
- 2026-10-14 scope add .gitignore — drop data/*.py ignore
- 2026-10-14 advance draft → ready
- 2026-10-14 scope approved: data/**, scripts/**, src/encoding/**, src/crypto/**, src/rings/**, README.md, .gitignore
- 2026-10-14 advance ready → in-progress/implement
- 2026-10-14 advance in-progress/implement → in-progress/verify
- 2026-10-14 advance in-progress/verify → in-progress/review
- 2026-10-14 advance in-progress/review → approval
- 2026-10-14 approved: looks good! let's close it
- 2026-10-14 advance approval → done
