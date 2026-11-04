# Current Task Context: Publication review: architecture, code, docs

State: bootstrap chain 0038-0044 active (0038-0042 done; 0043 awaiting operator acceptance; 0044 draft); 0037 blocked until the chain lands, then the repository review resumes.

## Plan

- [x] Task 0005-remove-dead-code (dead modules, golden vectors, scripts).
- [x] Task 0006-trim-dependencies.
- [x] Agree the backend design (now `specs/overview.md` Target architecture) and revise 0007-0012 around it.
- [x] Task 0007-module-layout (`foo.rs` + `foo/`; `math`, `rings/backends`, `ckks`).
- [x] Task 0008-ciphertext-levels (context chain, `Ciphertext::rescale`/`mod_switch_to`).
- [x] Task 0009-ops-composition (`ckks/ops.rs`, `ckks/evaluator.rs`).
- [x] Task 0010-chebyshev-module (`ckks/algorithms`, `pipeline.rs` retired).
- [x] Task 0011-backend-trait (one ring trait, NTT domain private).
- [x] Task 0012-generic-key-switching (one `KeySwitchingKey`, scheme generic, legacy keys out).
- [ ] Audit what remains; operator marks findings with `^^^`; split accepted ones into tasks (`tracing` logging, errors vs panics, comment/ASCII pass, examples, README + CI).

## Findings

- **Current audit:** task 0037 (cleanup-and-refactoring, blocked on 0038-0043). Authoritative implemented/remaining checklist: `specs/research/generic-bootstrap.md`; inventory/experiment log: `specs/research/file-review.md`. General encrypted bootstrap is incomplete; fixed-profile checks remain. "Toy" means experimental status, not restricted functionality. Keep default metric thresholds.
- **Historical observations below:** "today" describes the original publication audit, not current code. Negative rotations now use inverse offsets without conjugation; encryption/key noise inconsistency is confirmed and needs a separate numerical decision.

- **Order:** compose first on the concrete `RnsPoly` (0008-0010), then generalize the small, clean ops behind the trait (0011-0012).
- **Rescale today:** `RnsPoly::rescale_into` (exact RNS division) wrapped by `CkksEngine::rescale_ciphertext`. **Mod switch today:** `pipeline::align_level` truncates channels, caller passes `logq`. `PolyRescale`/`PolyModSwitch` have no impls.
- **Leak today:** `CkksEngine<P>` is generic only up to add; rescale, mul_plaintext, gadget keys and key switch are `impl CkksEngine<RnsPoly>` walking `moduli()`/`channels()[i]`; the encoder and `pipeline.rs` read RNS channels too.
- **Decided 2026-10-13:** target architecture accepted; degree stays a compile-time `const N`; reference backend on `crypto-bigint` is M2.6, after this cleanup.
- **Decided 2026-10-14 (0008):** basis links to the next basis; `logq` derived from the basis, not stored; `ct.mod_switch_to(level)` + `ct.match_level(&other)`; level errors are `CkksResult`.
- **Decided 2026-10-14 (0009):** `mul`/`mul_plain`/`rotate` are `Ciphertext` methods; scalars via `mul_const(value, scale_bits)` and `add_const(value)`; rescale always explicit.
- **Decided 2026-10-16 (0011):** RNS-NTT keeps NTT form at rest (OpenFHE/SEAL style) behind a private flag, converting only where coefficients are needed; one task for all six trait groups; trait keeps the name `PolyRing`. OpenFHE is the reference (SEAL unmaintained).
- **Audit candidates:** `encrypt` passes `error_variance` (3.2) to `sample_gaussian` as the std_dev (true std_dev 1.79; fixing changes noise, own task); `tracing` spans/events in ckks ops, subscriber in examples; `rotate_slots(-k)` is "rotate by k, then conjugate", not a left rotation (on real data equals `rotate(+k)`; left rotation is `5^(N/2-k)`); `table.rs` (comfy-table, 3 examples); `ckks/builder.rs` (generic since 0012, still unused); overlapping examples (`std_dev_8`, `sqrt_chebyshev_plain`); 2 rustdoc warnings (`basis.rs` intra-doc link, bare URL in `primes.rs`); stale ignored doctests in `keys/{secret,public}_key.rs` (`RnsBasisBuilder`, `RnsNttPoly`); file renames like `ckks_encoder` -> `encoder`; 179 non-ASCII lines in `src/`; README refresh: badges (Tests, Code Quality, Codecov, License; crates.io/docs.rs only after `license`/`description`/`repository` land in Cargo.toml), examples table adds `sensor_std_dev`/`sensor_std_plain`/`sqrt_chebyshev_plain`, Current Focus/What Works name the generic `PolyRing` scheme, slot sums and Chebyshev, snippet moves to a compiled `examples/quickstart.rs` (decided 2026-10-17).

## Context

Active cleanup task0037: review every repository-owned file with the operator, maintain the per-file checklist and agree scoped refactoring tasks. Preserve the publication audit plan and findings above; recheck historical observations against current code. Tasks0034/0035 remain future demo drafts.

## Next

Next: operator accepts 0043; then 0052 or 0045 (baseline: `specs/benchmarks.md`). Research: 0023/0024/0025, `specs/research/large-degree-bootstrap.md`.

Other drafts 0015-0022 (wasm, debug tooling, operation benchmarks, SIMD, GPU, Apple NPU, reference backend) wait for operator review; the `tracing` candidate moved to 0017. Side task 0013 (educational Sage scripts) remains independent.
