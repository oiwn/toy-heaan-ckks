---
id: 0015-wasm-feasibility
status: draft
scope: [specs/research/wasm.md, specs/demo.md]
created: 2026-10-16
source: operator request 2026-10-16
---
# Task: wasm feasibility

Can the library run in the browser for a demo? Follows `specs/demo.md` Phase B; go/no-go data, no UI.

## Plan

- [ ] Try `cargo build --target wasm32-unknown-unknown`; ask before adding the target or touching the toolchain.
- [ ] Find what blocks: getrandom wasm backend under rand 0.10, `comfy-table`, `crypto-bigint`, `rustfft`, `tracing`.
- [ ] Plan a minimal wasm-bindgen harness (keygen, encrypt, one mul, one rotation) and record timings and peak memory if it builds.
- [ ] Write `specs/research/wasm.md`: blockers, fixes, timings, go/no-go; replace the estimates in `specs/demo.md`.

## Acceptance

- [ ] `specdev check`
- [ ] `cargo build --target wasm32-unknown-unknown` result recorded in the note

## Log

- 2026-10-16 created: wasm feasibility
- 2026-10-16 set source operator request 2026-10-16
- 2026-10-16 scope add specs/research/wasm.md — research note; replace demo.md estimates
- 2026-10-16 scope add specs/demo.md — research note; replace demo.md estimates
