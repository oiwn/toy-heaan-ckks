# Current Task Context: Publication review: architecture, code, docs

State: not started. Roadmap M2.5.

## Plan

- [x] Task 0005-remove-dead-code (dead modules, golden vectors, scripts).
- [x] Task 0006-trim-dependencies.
- [x] Agree the backend design (now `specs/overview.md` Target architecture) and revise 0007-0012 around it.
- [x] Task 0007-module-layout (`foo.rs` + `foo/`; `math`, `rings/backends`, `ckks`).
- [x] Task 0008-ciphertext-levels (context chain, `Ciphertext::rescale`/`mod_switch_to`).
- [x] Task 0009-ops-composition (`ckks/ops.rs`, `ckks/evaluator.rs`).
- [x] Task 0010-chebyshev-module (`ckks/algorithms`, `pipeline.rs` retired).
- [ ] Task 0011-backend-trait (one ring trait, NTT domain private).
- [ ] Task 0012-generic-key-switching (one `KeySwitchingKey`, scheme generic, legacy keys out).
- [ ] Audit what remains; operator marks findings with `^^^`; split accepted ones into tasks (`tracing` logging, errors vs panics, comment/ASCII pass, examples, README + CI).

## Findings

- **Order:** compose first on the concrete `RnsPoly` (0008-0010), then generalize the small, clean ops behind the trait (0011-0012).
- **Rescale today:** `RnsPoly::rescale_into` (exact RNS division) wrapped by `CkksEngine::rescale_ciphertext`. **Mod switch today:** `pipeline::align_level` truncates channels, caller passes `logq`. `PolyRescale`/`PolyModSwitch` have no impls.
- **Leak today:** `CkksEngine<P>` is generic only up to add; rescale, mul_plaintext, gadget keys and key switch are `impl CkksEngine<RnsPoly>` walking `moduli()`/`channels()[i]`; the encoder and `pipeline.rs` read RNS channels too.
- **Decided 2026-10-13:** target architecture accepted; degree stays a compile-time `const N`; reference backend on `crypto-bigint` is M2.6, after this cleanup.
- **Decided 2026-10-14 (0008):** basis links to the next basis; `logq` derived from the basis, not stored; `ct.mod_switch_to(level)` + `ct.match_level(&other)`; level errors are `CkksResult`.
- **Decided 2026-10-14 (0009):** `mul`/`mul_plain`/`rotate` are `Ciphertext` methods; scalars via `mul_const(value, scale_bits)` and `add_const(value)`; rescale always explicit.
- **Audit candidates:** `tracing` spans/events in ckks ops, subscriber in examples; `table.rs` (comfy-table, 3 examples); `crypto/builder.rs`; workaround re-exports in `lib.rs`; overlapping examples (`std_dev_8`, `sqrt_chebyshev_plain`); 6 rustdoc warnings (broken intra-doc links, bare URLs); file renames like `ckks_encoder` -> `encoder`; 179 non-ASCII lines in `src/`.

## Context

Drafts 0005-0012 in `specs/tasks/`, chained by `depends`.

## Next

Plan 0011-backend-trait. Side task 0013 (educational Sage scripts) is a draft, independent of the chain.
