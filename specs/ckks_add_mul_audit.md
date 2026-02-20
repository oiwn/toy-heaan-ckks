# CKKS Add/Mul Audit (2026-02-03)

## Scope
- Focus: ciphertext addition and multiplication paths (engine + operations) for the RNS-NTT backend.
- Files inspected: `src/crypto/engine.rs`, `src/crypto/operations.rs`,
  `src/crypto/types.rs`, `src/rings/traits.rs`,
  `src/rings/backends/rns/ntt_poly.rs`, `src/rings/backends/rns/poly.rs`,
  `src/keys/relin_key.rs`, `src/keys/public_key.rs`, `src/encoding/fft.rs`,
  `examples/rns.rs`, `examples/ntt.rs`, `specs/overview.md`.

## Summary
- Homomorphic addition is implemented and works for RNS-NTT as long as both
  operands are in the same domain and have matching scale/modulus metadata.
- Homomorphic multiplication is incomplete: rescaling and key switching are
  stubbed or partially implemented, modulus tracking is inconsistent, and no
  backend implements `PolyRescale`.
- "Neuro-slope" (AI-generated low-quality scaffolding) is present in core
  crypto flows (duplicate APIs, placeholders, commented-out tests), and it
  directly obscures correctness gaps.

## Findings (ordered by impact)

### Critical
1. **Rescaling is not implemented for any backend, but required by
   multiplication.**
   - `PolyRescale` is declared with `rescale_assign`, yet no concrete impl exists.
   - `multiply_ciphertexts_kim` calls `rescale_assign` on `P`, so any real use
     will fail to compile for concrete backends.
   - Files: `src/rings/traits.rs`, `src/crypto/operations.rs`.

2. **Key switching step is skipped in multiplication, so the output scale/modulus
   is not properly reduced.**
   - In `multiply_ciphertexts_kim`, the right-shift by `logQ` is commented out
     and explicitly marked TODO. This is a required step in Kim’s algorithm.
   - In `CkksEngine::mul_ciphertexts`, relinearization is applied without any
     modulus management or rescale; `logq` is kept constant while the scale grows.
   - Files: `src/crypto/operations.rs`, `src/crypto/engine.rs`.

3. **`logq` tracking is inconsistent and often wrong at call sites.**
   - `CkksEngine::encrypt` requires a `logq` parameter, but examples pass
     `SCALE_BITS` (precision) instead of modulus bits.
   - `operations::encrypt` hard-codes `logq = 60` as a placeholder.
   - This breaks scale/modulus bookkeeping in addition/multiplication and makes
     comparisons (`logq` equality) unreliable.
   - Files: `src/crypto/engine.rs`, `src/crypto/operations.rs`,
     `examples/rns.rs`, `examples/ntt.rs`.

### High
4. **Error parameter semantics are inconsistent (variance vs std dev).**
   - `CkksParams.error_variance` is passed directly into `sample_gaussian`, which
     expects a standard deviation. The name suggests variance, but usage implies
     std dev.
   - `generate_public_key` ignores `CkksParams` and hard-codes `3.2`.
   - `generate_relinearization_key` uses `sqrt(error_variance)` which would be
     correct only if `error_variance` truly is a variance.
   - This impacts noise growth and therefore multiplication correctness.
   - Files: `src/crypto/engine.rs`, `src/keys/public_key.rs`,
     `src/keys/relin_key.rs`, `src/rings/traits.rs`.

5. **Multiplication does not enforce scale compatibility and has no path to
   align scales.**
   - `CkksEngine::mul_ciphertexts` only checks `logq`; it never checks or aligns
     `logp` before multiplying.
   - `multiply_ciphertexts_kim` requires equal `logp` but has no scale alignment
     helper when they differ.
   - This makes real-world usage fragile when ciphertexts originate from
     different levels.
   - Files: `src/crypto/engine.rs`, `src/crypto/operations.rs`.

### Medium
6. **API duplication introduces inconsistent behavior.**
   - There are two parallel encrypt/decrypt/add/mul APIs: `CkksEngine` methods
     and free functions in `crypto::operations`.
   - The free functions use hard-coded `logq` and different sampling
     parameters (`DEGREE/2` for `u`) than the engine, so behavior diverges.
   - Files: `src/crypto/engine.rs`, `src/crypto/operations.rs`.

7. **Addition and multiplication rely on panics instead of recoverable errors
   in some paths.**
   - `CkksEngine::add_ciphertexts` asserts on `logp/logq` mismatch, while
     `operations::add_ciphertexts` returns an error.
   - This difference is unexpected for callers and complicates robust handling.
   - Files: `src/crypto/engine.rs`, `src/crypto/operations.rs`.

### Low
8. **Test and example coverage for multiplication is effectively absent.**
   - Multiplication tests are commented out in `tests/basic_ckks.rs` and
     `tests/plaintext_mul.rs`.
   - Examples only exercise addition and basic polynomial multiplication, not
     full ciphertext multiplication.
   - Files: `tests/basic_ckks.rs`, `tests/plaintext_mul.rs`, `examples/rns.rs`,
     `examples/ntt.rs`.

## Neuro-slope (AI-generated code) indicators
Based on your definition ("AI-generated low-quality code"), these patterns are
present in the add/mul surface area:

1. **Placeholders that look generated rather than designed**
   - `operations::encrypt` hard-codes `logq = 60` with a TODO note.
   - `multiply_ciphertexts_kim` skips a required step and leaves a TODO.
   - Files: `src/crypto/operations.rs`.

2. **Duplicate APIs with divergent behavior**
   - Engine methods vs free functions differ in parameter handling and error
     strategy; this is typical of scaffolded code that was never reconciled.
   - Files: `src/crypto/engine.rs`, `src/crypto/operations.rs`.

3. **Commented-out tests/examples**
   - Large blocks in `tests/basic_ckks.rs` and `tests/plaintext_mul.rs` are
     commented out, so correctness coverage is effectively absent.
   - File: `tests/basic_ckks.rs`, `tests/plaintext_mul.rs`.

4. **Inconsistent parameter semantics**
   - `error_variance` vs `std_dev` handling is ambiguous and conflicting in
     key generation and encryption. This is a frequent AI-code smell where
     naming and usage drift.
   - Files: `src/crypto/engine.rs`, `src/keys/public_key.rs`,
     `src/keys/relin_key.rs`, `src/rings/traits.rs`.

These items are not proof of AI origin, but they are consistent with
low-quality generated scaffolding and should be prioritized for cleanup.

## Recommendations (toward a full add/mul schema)
1. **Unify the API**: choose `CkksEngine` as the single user-facing path and
   remove or internalize the older `crypto::operations` functions.
2. **Define a single source of truth for scale/modulus**: derive `logq` from
   the modulus chain in the context rather than passing it manually.
3. **Implement rescale for the RNS-NTT backend**:
   - Drop the last prime (modulus chain) and divide coefficients by that prime
     (rounded) to keep scale stable.
   - Update both `logp` and `logq` accordingly.
4. **Implement key switching properly**:
   - Add gadget decomposition parameters and apply them in relinearization.
5. **Add noise instrumentation**:
   - Track an estimated noise bound per ciphertext and log it after add/mul.
6. **Restore and expand tests**:
   - Re-enable multiplication tests with deterministic parameters and assert
     on error thresholds after rescale.

## Open questions
- What did “neuroslope” mean in your request? If it is not noise growth,
  please clarify the intended definition so I can audit that specifically.
- Should we target the Kim multiplication path or a simpler textbook
  relinearization + rescale first for the “full schema” milestone?
