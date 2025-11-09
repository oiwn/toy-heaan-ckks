# Project Overview
- Educational CKKS codebase centered on an RNS-NTT backend for polynomial arithmetic.
- Goal: implement CKKS exactly as described in the FHE Textbook (<https://fhetextbook.github.io/CKKSScheme.html#ckks-scheme>) so every engine feature maps back to that reference.

# Next Session Prompt
- Finish porting the remaining engine structs onto `RnsNttPoly` (plaintexts, ciphertexts, engine wiring) and replace the schoolbook fallback in `RnsNttPoly::mul_assign` with a true negacyclic NTT multiply so ciphertext add/mul stays in the transform domain end-to-end.

# Current Task Context

## Active Initiative · RNS-NTT Reset
- Goal: strip the project down to a single RNS-NTT backend and rebuild CKKS operations from there, using <https://fhetextbook.github.io/CKKSScheme.html#ckks-scheme> as the behavioral reference.
- Scope: remove (not just ignore) the naive and bigint backends plus any abstractions that only existed to support them; keep room to experiment with multiple RNS parameter sets that share the same code paths.
- Working agreement: keep this file as the day-to-day execution log; AGENTS.md already records meta-process rules, so avoid duplicating them here.

## Execution Priorities (near term)
1. **Purge legacy backends** – ✅ done; keep the codebase single-backend.
2. **Baseline NTT polynomial ops** – ✅ `RnsBasis` now ships NTT tables plus the channel-first `RnsNttPoly` with tests; next up is replacing the temporary schoolbook fallback in `mul_assign` with a real negacyclic NTT multiply.
3. **Engine sanity pass** – rewire the CKKS engine just enough to compile/link using the single backend, targeting plaintext encode/decode and ciphertext addition as the first working operations.

## Guardrails & Notes
- Work entirely in ASCII (no Greek letters or Unicode symbols) and prefer descriptive English names for former Greek variables.
- When pulling formulas or algorithms, cite the exact subsection URL from the FHE textbook.
- Specs live in `specs/overview.md`; use this context file for tactical decisions, open questions, and deviations from the plan.

_Last updated: reset toward single-backend RNS-NTT build._

## Phase Plan (iterates toward working addition)

### Phase 1 · Backend cleanup + compile check
- Remove the naive and bigint modules plus any builder/trait indirection that only serves them (keep the `backends` namespace so alternate RNS variants can drop in later).
- Wire `cargo test` to only exercise the RNS-NTT structures; stubs are fine for anything that used to rely on the deleted backends, as long as the crate compiles.
- Snapshot the previously prototyped RNS components (notably `RnsBasis`, `RnsNttPoly`, and the CRT helpers) and decide what to port verbatim vs. rewrite.

### Phase 2 · NTT-first polynomial operations
- Get `RnsNttPoly` into the codebase with channel-first storage, an `in_ntt_domain` flag, and single-source NTT/INTT routines.
- Restrict all public arithmetic (add, sub, scalar mul) to operate on NTT representations; only expose coefficient views for encoding/decoding paths.
- Build a minimal parameter loader (hard-coded primes for now) so we can instantiate the basis from tests/examples.

### Phase 3 · Engine + addition demo
- Update `CkksEngine`, keys, plaintexts, and ciphertexts to depend solely on the new backend without renaming the public API.
- Ensure encoding→encryption→decryption works for a single level and that ciphertext addition in the NTT domain matches plaintext addition after decode.
- Produce at least one example or integration test that exercises the flow end-to-end (no rescale/relinearization yet).

## Deferred Work (parking lot)
- **Modulus chain + level context** – full textbook rescale/mod-drop logic per <https://fhetextbook.github.io/CKKSScheme.html#rescaling>. Keep notes in the scratch docs but do not implement until addition works.
- **Relinearization + gadget decomposition** – port the `get_gadget_vector`/phase2 tests once multiplication is on deck; needs proper evaluation keys and `s^2` handling.
- **Scale tracking & metadata** – upgrading ciphertexts with `{level_idx, scale_bits}` plus rescale helpers belongs to a later phase once multiplication exists.

## Open Questions
- **Trait surface area** – keep the existing `PolyRing` shims for now so alternate RNS variants (or profiling builds) can slot in without rewriting the engine.
- **Prime selection defaults** – currently undecided; will revisit after Phase 3 proves out the flow.

## Progress Log
- 2025-11-07 — completed the destructive part of Phase 1 by moving the naive/bigint backends, their encoders, benches, and legacy examples/tests into `tmp/`, pruning the builder/export surface down to the RNS path, and making the unit tests/docs reference `RnsPolyRing`.
- 2025-11-07 — queued up Phase 2 (NTT-first polynomial ops) for the next session; no further changes until we resume with the `RnsNttPoly` port.
- 2025-11-08 — unblocked `cargo test` by making the secret-key tests/doc snippets specify their const-degree basis and import `PolyRing`, so we have a clean baseline before wiring up the NTT polynomial work.
- 2025-11-08 — removed the stale `ct_mul_demo` example, taught `RnsBasis` to carry NTT tables, added the channel-first `RnsNttPoly` plus toy parameter loader, and landed regression tests so Phase 2 now has a concrete NTT surface to build on.
- 2025-11-09 — aligned the NTT table/multiply implementation with the textbook algorithm so `RnsNttPoly::mul_assign` stays fully in the transform domain, and rewired the builder, keys, encoder, and examples to instantiate `RnsNttPoly` end-to-end (plaintexts now enter encryption already in NTT form).
