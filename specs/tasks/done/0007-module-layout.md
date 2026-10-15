---
id: 0007-module-layout
status: done
scope: [src/**, examples/**, tests/**, benches/**, README.md]
created: 2026-10-13
source: user request 2026-10-13: publication cleanup, remove ancient code and golden vectors, foo.rs module layout
depends: [0006-trim-dependencies]
---
# Task: module layout

## Plan

Layers from `specs/overview.md` Target architecture; `foo.rs` + `foo/` everywhere, no `mod.rs`.

- [x] `src/math.rs` + `math/{primes, sampling, utils}` (only the `mod.rs` move).
- [x] `src/rings.rs` + `rings/traits.rs` + `rings/backends.rs` + `rings/backends/rns_ntt.rs` + `rns_ntt/{basis, poly, errors}` (keeps `backends/` for pluggable backends).
- [x] `src/ckks.rs` + `ckks/`: move in `crypto/*` (engine, builder, errors, pipeline, types), `encoding/*` -> `ckks/encoding.rs` + `ckks/encoding/{ckks_encoder, special_fft}` (file names kept; renames are an audit item) and `keys/*` -> `ckks/keys.rs` + `ckks/keys/`.
- [x] Update `use` paths in src, examples, tests, benches; `lib.rs` declares `ckks`, `math`, `rings`, `table`.
- [x] Mechanical only: no logic, item, or signature changes; only module paths change (`crypto::`, `encoding::`, `keys::` -> `ckks::...`).
- [x] No `mod.rs` left anywhere; update paths in specs/overview.md, specs/demo.md, README.md and AGENTS.md (Key File Paths).

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features`
- [x] `cargo build --examples`
- [x] `cargo bench --no-run`
- [x] `typos`
- [x] `gitleaks dir .`

## Summary

Moved to the `foo.rs` + `foo/` layout: `crypto`, `encoding`, `keys` now live under `ckks/`; `math` and `rings/backends/rns_ntt` lost their `mod.rs`; import paths and docs updated, no logic change.

## Context

Moves are plain file moves; the operator stages them in git. All moves happen here once, so 0008-0012 only change code inside the new layout.

## Review

- Pure move plus path rewrites (Edit tool); `cargo fmt` reordered the touched `use` blocks. `super::` paths stayed valid since subtrees moved intact. Crate-root key re-exports now point at `ckks::keys` (cleanup is 0012).
- Same test counts as before: 122 lib + 7 integration + 4 example + 1 doc; benches and examples build; no `mod.rs` under `src/`; `rg` finds no `crypto::`/`encoding::`/`keys::` crate paths.
- `cargo doc` shows 6 warnings that predate this task (broken `[X]`, `[rescale]`, `[rescale_into]` links, bare URLs in `basis.rs`, `poly.rs`, `types.rs`, `primes.rs`): audit item.

## Log

- 2026-10-13 created: module layout
- 2026-10-13 set source user request 2026-10-13: publication cleanup, remove ancient code and golden vectors, foo.rs module layout
- 2026-10-13 set depends 0006-trim-dependencies
- 2026-10-13 scope add src/crypto.rs — mod.rs -> <dir>.rs moves
- 2026-10-13 scope add src/encoding.rs — mod.rs -> <dir>.rs moves
- 2026-10-13 scope add src/keys.rs — mod.rs -> <dir>.rs moves
- 2026-10-13 scope add src/math.rs — mod.rs -> <dir>.rs moves
- 2026-10-13 scope add src/rings.rs — mod.rs -> <dir>.rs moves
- 2026-10-13 scope add src/**/mod.rs — mod.rs -> <dir>.rs moves
- 2026-10-13 scope add src/rings/backends.rs — nested mod.rs moves
- 2026-10-13 scope add src/rings/backends/rns_ntt.rs — nested mod.rs moves
- 2026-10-13 scope add src/** — mechanical mod.rs -> <dir>.rs moves across all src modules
- 2026-10-13 scope rm src/crypto.rs
- 2026-10-13 scope rm src/encoding.rs
- 2026-10-13 scope rm src/keys.rs
- 2026-10-13 scope rm src/math.rs
- 2026-10-13 scope rm src/rings.rs
- 2026-10-13 scope rm src/**/mod.rs
- 2026-10-13 scope rm src/rings/backends.rs
- 2026-10-13 scope rm src/rings/backends/rns_ntt.rs
- 2026-10-13 scope add examples/** — import paths change with the module moves
- 2026-10-13 scope add tests/** — import paths change with the module moves
- 2026-10-13 scope add benches/** — import paths change with the module moves
- 2026-10-14 scope add README.md — layout paths in README
- 2026-10-14 advance draft → ready
- 2026-10-14 scope approved: src/**, examples/**, tests/**, benches/**, README.md
- 2026-10-14 advance ready → in-progress/implement
- 2026-10-14 advance in-progress/implement → in-progress/verify
- 2026-10-14 advance in-progress/verify → in-progress/review
- 2026-10-14 advance in-progress/review → approval
- 2026-10-14 approved: tests passing! let's discuss next task.
- 2026-10-14 advance approval → done
