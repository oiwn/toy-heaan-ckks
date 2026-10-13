---
id: 0007-module-layout
status: draft
scope: [src/**, examples/**, tests/**, benches/**]
created: 2026-10-13
source: user request 2026-10-13: publication cleanup, remove ancient code and golden vectors, foo.rs module layout
depends: [0006-trim-dependencies]
---
# Task: module layout

## Plan

Target layout (see `specs/ctx.md` Findings); `foo.rs` + `foo/` everywhere, no `mod.rs`.

- [ ] `src/math.rs` + `math/{primes, sampling, utils}` (only the `mod.rs` move).
- [ ] `src/rings.rs` + `rings/{basis, poly, errors, traits}`: flatten `rings/backends/rns_ntt/` (one backend).
- [ ] `src/ckks.rs` + `ckks/`: move in `crypto/*` (engine, builder, errors, pipeline, types), `encoding/*` (`ckks_encoder` -> `encoding`, `special_fft`) and `keys/*` -> `ckks/keys/`.
- [ ] Update `use` paths in src, examples, tests, benches; `lib.rs` declares `math`, `rings`, `ckks`, `table`.
- [ ] Mechanical only: no logic, item, or signature changes; module paths change (`crypto::` -> `ckks::`, `rings::rns_ntt::` -> `rings::`).
- [ ] No `mod.rs` left anywhere; update paths in specs/overview.md and AGENTS.md (Key File Paths).

## Acceptance

- [ ] `cargo fmt --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features`
- [ ] `cargo build --examples`
- [ ] `typos`
- [ ] `gitleaks dir .`

## Context

On hold: revised after the backend design in `specs/ctx.md` is agreed (keeps `rings/backends/`). Moves are plain file moves; the operator stages them in git. All moves happen here once, so 0008-0010 only change code inside the new layout.

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
