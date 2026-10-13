# Roadmap

Committed future direction, in rough priority order. Promote items here from `ideas.md` when decided; move items into `ctx.md` when they become the active task.

## Next

- **M2.5 publication cleanup** — remove dead code, golden vectors and unused deps; `foo.rs` + `foo/` module layout (`math`, `rings`, `ckks`); one explicit rescale/mod-switch API; dissolve `pipeline.rs` into composable ops (fhe.rs / tfhe-rs style); logging with `tracing`; audit architecture, composition, public API, comments and docs so the code reads as written on purpose. No optimizations. Tasks 0005-0010, then audit tasks.

## Later

- **M3 browser demo** — follows `specs/demo.md` Phases B–D:
  - serde/postcard for `RnsPoly`, `Ciphertext`, `Plaintext`, keys; `RnsBasis` rebuilt from the prime list so polys share one `Arc` after load; round-trip test.
  - wasm build: getrandom wasm backend; `cargo build --target wasm32-unknown-unknown` green.
  - wasm spike with timings, then the Leptos UI.
- **M4 lint debt** — make `cargo clippy -- -D clippy::pedantic -D clippy::nursery` pass (254 errors on 2026-10-12, mostly `must_use`, `const fn`, `# Errors`/`# Panics` docs, casts); work through `specs/cleanup.md`.
