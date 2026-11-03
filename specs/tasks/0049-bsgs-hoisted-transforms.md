---
id: 0049-bsgs-hoisted-transforms
status: draft
scope: [src/ckks/algorithms/bootstrap/*.rs, src/ckks/algorithms/bootstrap.rs, src/ckks/keys/key_switching.rs]
created: 2026-11-03
source: Operator agreed large-degree roadmap, 2026-11-03
depends: [0048-hybrid-key-switching]
---
# Task: bsgs hoisted transforms

## Context

Roadmap step 4 (`specs/research/large-degree-bootstrap.md`). CtS/StC groups apply every diagonal with its own rotation and key switch. Lattigo (`MultiplyByDiagMatrixBSGS`) and OpenFHE (`EvalFastRotationExt`) split each group into baby and giant steps, decompose the input once, keep baby-step rotations in QP and ModDown once per giant step: fewer key switches and roundings, far fewer rotation keys. At N=8192 this is the main lever on time (40 s per bootstrap) and key memory. Needs the QP basis and hybrid keys from 0048.

Rules: touch only scoped files (extend with `specdev task scope` first). If the same check fails twice, stop, record the failure and a hypothesis in Findings and ask the operator. Never weaken tolerances, disable noise or add `#[ignore]`.

## Plan

- [ ] Hoisted rotation: one digit decomposition of a ciphertext, reused for several rotations, results left in QP.
- [ ] BSGS evaluation of a transform group (baby/giant split chosen per group, giant-step plaintexts pre-rotated), one ModDown per giant step.
- [ ] The schedule derives the reduced rotation-key set; key requirements and memory estimates follow.
- [ ] Bit-for-bit or error-bounded agreement with the current transform; record time, key count and error at N=1024/8192 in the research doc.

## Acceptance

- [ ] `cargo fmt --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features`
- [ ] `cargo build --examples`
- [ ] `cargo run --release --example bootstrap_demo`
- [ ] `typos`
- [ ] `specdev check`
- [ ] `git diff --check`

## Log

- 2026-11-03 created: bsgs hoisted transforms
- 2026-11-03 set source Operator agreed large-degree roadmap, 2026-11-03
- 2026-11-03 set depends 0048-hybrid-key-switching
- 2026-11-03 scope add src/ckks/algorithms/bootstrap/*.rs — Hoisted rotations and BSGS transform groups
- 2026-11-03 scope add src/ckks/algorithms/bootstrap.rs — Hoisted rotations and BSGS transform groups
- 2026-11-03 scope add src/ckks/keys/key_switching.rs — Hoisted rotations and BSGS transform groups
