# Current task context

## 2026-02-04 decisions
- Scope: refactor basic RNS math only. Do not touch encoder or crypto code.
- Goal: remove "neuro-slope" (low-quality scaffolding) in RNS math, per
  specs/ckks_add_mul_audit.md.
- Architecture: merge `RnsPolyRing` and `RnsNttPoly` into a single RNS
  polynomial type that can switch between coefficient and NTT domains.
- Basis equality: use strict pointer equality (`Arc::ptr_eq`) for compatibility
  checks.

## Plan of action (spec-first, no code yet)
1. Draft `specs/ctx.md`:
   - define the unified RNS polynomial type and domain switching API
   - list acceptance criteria (tests, behavior, public API changes)
2. Identify all RNS math duplication to remove:
   - `mod_inverse` / `crt_reconstruct` duplicates
   - sampling helpers (uniform/gaussian/ternary) used by RNS backends
   - duplicated add/mul/neg and mod-switch paths
3. Define the migration steps and expected file edits.
4. Review the spec with operator before any code changes.

## Open questions
- Confirm the exact public API surface for the merged RNS type
  (constructor names, domain switching methods, trait impls).
- Decide where shared RNS math helpers live (new module vs reuse `math/`).
