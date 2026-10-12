# Roadmap

Committed future direction, in rough priority order. Promote items here from `ideas.md` when decided; move items into `ctx.md` when they become the active task.

## Next

- **M1 CKKS green** — `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`, and `cargo test --all-features` all pass. Includes widening the encoder beyond its i64 staging (task 0004). Tasks 0001, 0004, 0002.
- **M2 native sensor demo** — sensor range contract (client normalizes readings into [0, 1]) and an N = 8192 example that asserts the std-dev error budget. Tasks 0002, 0003.

## Later

- **M3 browser demo** — follows `specs/demo.md` Phases B–D:
  - serde/postcard for `RnsPoly`, `Ciphertext`, `Plaintext`, keys; `RnsBasis` rebuilt from the prime list so polys share one `Arc` after load; round-trip test.
  - Cargo/wasm cleanup: move `dhat`, `comfy-table`, `tracing-subscriber` to dev-dependencies; getrandom wasm backend; `cargo build --target wasm32-unknown-unknown` green.
  - wasm spike with timings, then the Leptos UI.
- **M4 lint debt** — make `cargo clippy -- -D clippy::pedantic -D clippy::nursery` pass (254 errors on 2026-10-12, mostly `must_use`, `const fn`, `# Errors`/`# Panics` docs, casts); work through `specs/cleanup.md`.
