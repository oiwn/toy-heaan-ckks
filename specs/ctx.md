# Project Overview
- Single-backend CKKS playground focused on the RNS-NTT pipeline.
- Canonical references: scheme definition <https://fhetextbook.github.io/CKKSScheme.html#ckks-scheme> and encoding/decoding flow <https://fhetextbook.github.io/EncodingandDecoding.html#encoding-and-decoding>.

# Active Tasks
1. **Textbook CKKS encoder/decoder** – replace the RustFFT/`homomorphix` remnants with the textbook algorithm: use Δ = 2^scale_bits, build the conjugate-symmetric slot vector, run the canonical Vandermonde/NTT-friendly transform, pack/unpack across the full CRT product, and keep the implementation co-located with the upcoming `@homomorphix.txt` NTT encode/decode spec.
2. **RNS-NTT reset** – hold the line on “only `RnsNttPoly` everywhere,” removing stray abstractions that assumed multiple backends while keeping space for alternate RNS parameterizations.
3. **Engine sanity** – finish wiring plaintext/ciphertext structs plus `CkksEngine` helpers so encode → encrypt → add → decrypt uses NTT-domain polys end-to-end with no schoolbook fallbacks.
4. **Docs alignment** – keep `specs/overview.md` and `specs/textbook_encode_decode.md` synchronized whenever backend or encode/decode decisions change (single-backend story, cached `Q`, BoxedUint-centered flow, slot-input adapter).

# Constraints & Notes
- ASCII identifiers only; replace Greek symbols with descriptive English terms.
- Cite textbook subsections when copying algorithms.
- `specs/overview.md` holds the architectural spec; this file tracks day-to-day intent.

# Phase Map
- **Phase 1 · Backend cleanup** – done: deleted legacy backends, ensured `cargo test` targets RNS-NTT, and staged the CRT/NTT helpers we’re porting.
- **Phase 2 · NTT-first ops** – `RnsNttPoly` + tables landed; outstanding action is enforcing negacyclic NTT multiply everywhere.
- **Phase 3 · Engine + demo** – rewire engine/key structs, prove encode→encrypt→decrypt and ciphertext addition before attempting rescale/relin.

# Deferred / Parking Lot
- Modulus chains + level tracking per <https://fhetextbook.github.io/CKKSScheme.html#rescaling>.
- Relinearization/gadget decomposition once multiplication stabilizes.
- Formal scale metadata (`{level_idx, scale_bits}`) after rescale work begins.

# Open Questions
- Keep `PolyRing` traits or shrink them now?
- Default prime sets for examples/tests once Phase 3 lands.

## FHE Book Decode/Encode
- **Goal**: implement the canonical CKKS encode/decode flow from <https://fhetextbook.github.io/EncodingandDecoding.html#encoding-and-decoding> for the single RNS-NTT backend. Messages up to `N/2` slots must map to polynomials in `ℤ_Q[X]/(X^N + 1)` using Δ = 2^scale_bits and the textbook embedding/recovery transforms.
- **Math spec**: use `N = 512` (slots = 256) as the baseline for realistic flows, while still leaning on degree 8/16 toy cases inside unit tests to illustrate the math. Build the conjugate-symmetric slot vector, run the Vandermonde/special FFT inverse (`V^{-1}`) to obtain integer coefficients, pack them across all CRT primes (`p_i ≡ 1 mod 1024`), and reverse the process by reconstructing coefficients via CRT before applying `V` and dividing by Δ.
- **Scope**: replace the RustFFT/homomorphix remnants with a `CkksEncoder` that (1) enforces textbook scaling, (2) validates NTT-friendly primes, (3) exposes debug helpers to print coefficients/residues and decoded slots, and (4) targets a default prime triplet for `N = 512` while staying generic over `N`.
- **Tests & tooling**: add encode/decode round trips (deterministic fixtures + randomized property tests), overflow guards (`|Δ · z'| < Q/2`), CRT sanity checks for the selected primes, and a debug example that logs inputs, coefficients, residues, and reconstruction error for degree 512. Keep instrumentation lightweight but always printable from examples/tests.
- **Implementation plan**:
  1. **Parameter plumbing**
     - Introduce `CkksEncodingParams<const N: usize>` that stores `degree`, `slots`, `scale_bits`, `primes: Arc<Vec<u64>>`, cached canonical roots, and derived scaling factor `Δ = 2^{scale_bits}`.
     - Validate invariants: `N.is_power_of_two()`, `slots = N/2`, `scale_bits > 0`, `values.len() ≤ slots`, and every prime satisfies `(prime - 1) % (2 * N) == 0`.
     - Provide helpers to print the active prime set and warn when `Δ` risks overflowing `Q/2`.
  2. **Canonical transforms**
     - Implement a conjugate-symmetric slot builder supporting real-only inputs (auto-wrap into `Complex64`), zero padding, and optional assertions for symmetry.
     - Provide Vandermonde-based `special_idft` (inverse embedding) and `special_dft` (forward embedding) that operate on cached roots; keep APIs generic over `N` with tests down to `N = 8`.
     - Leave TODO hooks for a faster special FFT so we can swap implementations without touching encoder logic.
  3. **Encoding pipeline**
     - `fn encode<const N: usize>(values: &[Complex64], params: &CkksEncodingParams<N>, basis: &Arc<RnsBasis>) -> EncodingResult<RnsNttPlaintext<N>>`.
     - Steps: enforce slot bound, build symmetric vector, multiply by `Δ`, run `special_idft`, assert imag parts ≈ 0, round to `i64`, pack residues via CRT reducers into `RnsNttPoly` in coefficient domain, store slot count in `Plaintext`.
     - Emit debug metrics (max coefficient magnitude vs `Q/2`, first few coefficients) behind a feature flag or optional logger for `examples/encode_debug`.
  4. **Decoding pipeline**
     - `fn decode<const N: usize>(plaintext: &Plaintext<RnsNttPoly<N>, N>, params: &CkksEncodingParams<N>) -> Vec<Complex64>`.
     - Steps: clone poly, `to_coeff_domain`, reconstruct signed coefficients with `convert_from_rns`, run `special_dft`, divide by `Δ`, truncate to `plaintext.slots`, and provide helper wrappers (`decode_real`) that drop imag parts.
     - Include helper to compute `(expected, decoded) -> max_abs_error` for tests/examples.
  5. **Integration & debugging**
     - Remove `RustFftEncoder` usage from examples/tests, replacing with `CkksEncoder` builder tied to new params.
     - Add coefficient/residue pretty-printers (per coefficient show signed value + residues per prime) and domain markers (coeff vs NTT).
     - Ship `examples/encode_debug.rs` (degree 512) that prints inputs, rounded coefficients, per-prime residues, decoded outputs, and error stats; also add a scripted way to dump first `k` coefficients for inspection.
  6. **Testing strategy**
     - Toy-degree (`N = 8/16`) unit tests: verify conjugate symmetry, exact encode/decode, and Vandermonde math steps.
     - Degree-512 tests: deterministic fixtures (e.g., arithmetic progression, sinusoids) with `ε ≤ 1e-6`, plus seeded random property tests for stability.
     - CRT tests: ensure `convert_to_rns`/`convert_from_rns` round-trip for the chosen primes; include overflow guard tests (`|Δ · z'|` hitting boundary).
     - Example smoke tests: run `cargo run --example encode_debug` in CI (or doc instructions) to confirm instrumentation compiles.
  7. **Docs & follow-ups**
     - Update README + specs to describe the textbook encoder, mention degree-512 baseline, and note remaining TODOs (complex slots policy, special FFT optimization).
     - Capture open questions in `specs/ctx.md` so future tasks can branch off (e.g., when to implement complex slot rotations, how to expose instrumentation in release builds).

## Progress Log
- 2025-11-10 — updated `specs/overview.md` to narrate the single RNS-NTT backend and refreshed `specs/textbook_encode_decode.md` to capture the newest encode/decode decisions (cached `Q`/`Q/2` via `crypto_bigint`, centered rounding with `BoxedInt`, and dual real/complex slot inputs); committed to logging future spec tweaks here so context stays current.
- 2025-11-07 — completed the destructive part of Phase 1 by moving the naive/bigint backends, their encoders, benches, and legacy examples/tests into `tmp/`, pruning the builder/export surface down to the RNS path, and making the unit tests/docs reference `RnsPolyRing`.
- 2025-11-07 — queued up Phase 2 (NTT-first polynomial ops) for the next session; no further changes until we resume with the `RnsNttPoly` port.
- 2025-11-08 — unblocked `cargo test` by making the secret-key tests/doc snippets specify their const-degree basis and import `PolyRing`, so we have a clean baseline before wiring up the NTT polynomial work.
- 2025-11-08 — removed the stale `ct_mul_demo` example, taught `RnsBasis` to carry NTT tables, added the channel-first `RnsNttPoly` plus toy parameter loader, and landed regression tests so Phase 2 now has a concrete NTT surface to build on.
- 2025-11-09 — aligned the NTT table/multiply implementation with the textbook algorithm so `RnsNttPoly::mul_assign` stays fully in the transform domain, and rewired the builder, keys, encoder, and examples to instantiate `RnsNttPoly` end-to-end (plaintexts now enter encryption already in NTT form).
- 2025-11-10 — drafted `specs/textbook_encode_decode.md`, added `src/encoding/textbook.rs` + `special_fft.rs` scaffolding, and parked placeholder tests/examples to stage the FHE-book encode/decode work.
