# Textbook Encode/Decode Spec

Authoritative reference for porting the CKKS textbook encoder/decoder into
this repository. Aligns with the "FHE Book Decode/Encode" scope captured in
`specs/ctx.md`.

## References
- FHE Textbook — CKKS scheme (<https://fhetextbook.github.io/CKKSScheme.html#ckks-scheme>)
- FHE Textbook — Encoding/Decoding (<https://fhetextbook.github.io/EncodingandDecoding.html#encoding-and-decoding>)
- Related code from book (<https://raw.githubusercontent.com/fhetextbook/fhetextbook.github.io/refs/heads/main/source%20code/ckks.py>)

## Goals
- Implement the canonical CKKS encode/decode algorithms using the Vandermonde
  (aka special FFT) transform and Δ = 2^scale_bits fixed-point scaling.
- Keep the implementation single-backend (`RnsNttPoly`) while remaining generic
  over the ring degree `N = 2^k`.

## Non-Goals
- Rescaling, relinearization, or ciphertext-level gadget work.
- Supporting legacy backends (there is none ATM); this effort is RNS-NTT only.
- Optimizing the transform beyond a straightforward Vandermonde implementation
  (but the API must leave room for a drop-in faster variant later).


## Parameter Object (`TextbookEncodingParams`)
- `const N: usize` with `N.is_power_of_two()` enforced.
- `slots = N / 2`; reject inputs exceeding this count.
- `scale_bits: u32` (Δ = 2^scale_bits) with helpers for overflow warnings when
  `Δ · |value| ≥ Q/2`.
- `Arc<Vec<u64>> primes` defining the CRT basis, validated via `(p_i - 1) % (2N) == 0`
  using `math::primes::is_ntt_friendly_prime`.
- Cached `Q = ∏p_i` and `Q/2` stored as `crypto_bigint::BoxedUint`, so overflow
  checks and centered-range comparisons never recompute the product.
- Cached canonical 2N-th roots `psi`, `psi_inv`, and Vandermonde twiddle tables
  stored in `Arc<Vec<Complex64>>` for reuse.

  ^^^ looks good to me

## Canonical Transforms
- `special_fft.rs` provides:
  - `build_conjugate_slots(values)` — accepts either `&[Complex64]` or `&[f64]`
    (via an enum/trait adapter) and produces a conjugate-symmetric
    `Vec<Complex64>` of length `N` with zero padding.
  - `special_idft(vec, params)` — Vandermonde-based inverse embedding that maps
    slots to coefficient-domain complex values.
  - `special_dft(coeffs, params)` — forward embedding returning slot values.
- The API exposes plain `Vec<Complex64>` for clarity and houses TODO hooks for
  a future fast special FFT implementation.

  ^^^ right, will need to  process different kind of inputs, i want to be able to use vector of floats (f64) and/or vector of complex numbers

## Encoding Pipeline Outline
1. Validate slot count and wrap real inputs into `Complex64`.
2. Build the conjugate-symmetric slot vector (`N` complex entries).
3. Multiply all slots by Δ.
4. Run `special_idft`, expecting near-real coefficients; assert the imaginary
   parts are below a tolerance configurable per test/example.
5. Round to a signed big integer (prefer `crypto_bigint::BoxedInt`) while
   asserting the centered magnitude `|coeff| < Q/2`; conceptually stay in the
   `[-Q/2, Q/2)` range even though storage will be positive residues.
6. Pack the coefficients into `RnsNttPoly` by reducing each centered integer
   modulo every prime (`(value % p_i + p_i) % p_i`) and then switching to NTT
   domain via `to_ntt_domain`.
7. Return a `Plaintext` carrying the polynomial, `scale_bits`, and slot count.

^^^ Do we need i64 or can use centered representation? or need to go back and forth? I saw in book,
required to go back and forth but i not sure it's valid for RNS-NTT

DO NOT FORGET WE ARE USING RNS-NTT and u64 for coefs it's all we have. But scale is valid question, but before we need to make RNS-NTT working at least for roundtrip.

## Decoding Pipeline Outline
1. Clone polynomial, convert to coefficient domain, and reconstruct integers via
   CRT (`convert_from_rns`) backed by `crypto_bigint::BoxedUint`; subtract `Q`
   whenever the recovered value exceeds `Q/2` to re-center.
2. Interpret the centered coefficients as `Complex64` and run `special_dft`.
3. Divide by Δ and truncate to the stored slot count.
4. Provide helpers for `decode_complex`, `decode_real`, and `max_abs_error`.

## Instrumentation & Debugging
- Feature-gated or opt-in logging that prints:
  - Selected input slots (e.g., first 8 values).
  - Max coefficient magnitude vs. `Q/2`.
  - Per-prime residues for the first `k` coefficients.
  - Decoded outputs plus max absolute error.
- `examples/encode_debug.rs` wires the logging together for `N = 512`.

## Testing & Acceptance Criteria
- Toy degrees (`N = 8/16`) covering:
  - Conjugate symmetry builder sanity (values mirrored correctly).
  - Vandermonde DFT/IDFT round trips vs. analytical expectations.
  - Encode/decode round trips with deterministic fixtures.
- Degree-512 tests:
  - Fixture-based inputs (arith progression, sinusoid) with
    `max_error ≤ 1e-6`.
  - Seeded random property test verifying encode→decode stability and CRT
    reconstruction accuracy.
- CRT reducer tests for the chosen prime set.
- Overflow guard tests forcing `|Δ · value|` near `Q/2`.
- `cargo run --example encode_debug` smoke test documented in README/Justfile.

## Integration Tasks
1. Introduce `TextbookEncodingParams` and `TextbookEncoder` in
   `src/encoding/textbook.rs`.
2. Add `special_fft.rs` with Vandermonde helpers shared by encode/decode.
3. Wire encoder into the new engine surface, phasing out `RustFftEncoder` in
   examples/tests once the implementation is stable.
4. Ship the `encode_debug` example plus targeted tests.
5. Update README/specs to reflect the new canonical encoder.

## Open Questions
- Default prime sets for `N = 512` fixtures; keep this synced with `RnsBasis`.
- Long-term location for transform caching (encoder vs. basis vs. global cache).
- How to expose instrumentation in release binaries without bloating output.
