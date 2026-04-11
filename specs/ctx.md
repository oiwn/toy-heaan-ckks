# Current task context

## Session status (2026-04-12) — all passing (101 tests)

### What's working
- Encode/decode: `examples/encode_decode.rs` ✓
- Encrypt/add/decrypt: `examples/encrypt_add.rs` ✓ (errors ~1e-7 at scale 2^30)
- Encrypt/mul/rescale/decrypt: `examples/encrypt_mul.rs` ✓ (errors ~5e-6 at scale 2^30)
- `RnsBasis<N>` + `RnsPoly<N>`: NTT forward/inverse, CRT, `PolyRing` + `PolySampler`
- `mul_assign`: NTT-based negacyclic convolution (pre/post-twist by ψ^j around circular NTT)
- `mul_assign_naive`: schoolbook O(N²), kept as reference
- Keys: `SecretKey`, `PublicKey`, `RelinearizationKey`
- `CkksEngine`: encrypt, decrypt, add_ciphertexts, mul_ciphertexts (schoolbook path)
- Single unified `Plaintext<P,N>` in `crypto::types`; encoder uses it directly

### What's NOT working / missing
- `PolyRescale` not implemented for `RnsPoly` — needed after multiplication
- `PolyModSwitch` not implemented for `RnsPoly` — needed for noise management
- `logq` in `Ciphertext` is a tag only; not yet tied to active RNS channels
- `mul_ciphertexts` in engine uses naive relin (broken for non-trivial params); use `mul_ciphertexts_gadget` + `RnsGadgetRelinKey` instead
- `PolyRescale` trait not wired to RNS (mismatched signature); `rescale`/`rescale_into` are direct methods on `RnsPoly`

## NTT implementation notes

- `NttTable` stores:
  - `forward_roots[i] = ω^i`, `inverse_roots[i] = ω^{-i}` where ω = ψ^2 (N-th root)
  - `twist_factors[j] = ψ^j`, `untwist_factors[j] = ψ^{-j}` where ψ = primitive 2N-th root
- `to_ntt_domain()`: pre-twist (×ψ^j) → forward NTT (circular DFT at ω^k)
- `to_coeff_domain()`: inverse NTT → post-untwist (×ψ^{-j})
- Pointwise multiply in NTT domain = negacyclic convolution in Z[X]/(X^N+1) ✓
- Bug that was fixed: roots were stored as ω^{bit_rev(i)} instead of ω^i, making the butterfly compute a wrong (invertible but non-DFT) transform

## Next goal: FHE multiplication

Needed for encode → encrypt → homomorphic multiply → rescale → decrypt → decode.

### Phase 1 — RNS Rescaling (core missing piece)

**Math:** After multiplication `logp` doubles. Rescaling divides every coefficient by the last
RNS prime `q_L` and drops that channel. In RNS this is done without full CRT reconstruction:

```
rescaled[i][j] = (c[i][j] - (c[L][j] mod q_i)) * inv(q_L, q_i)  mod q_i
```

This is exact (no rounding error) because `c[L][j]` is the exact residue mod `q_L`, so
`c - c_L` is exactly divisible by `q_L`.

**Changes:**
- `src/rings/backends/rns_ntt/poly.rs`: add `pub fn rescale(&self) -> Self`
  - Works in coeff domain (auto-convert from NTT, then restore domain flag)
  - Computes `q_L_inv[i] = mod_inverse(q_last, q_i)` inline
  - Returns new poly with L-1 channels via `new_unchecked` + `basis.drop_last(1)`
- `src/crypto/engine.rs`: add `pub fn rescale_ciphertext(ct: &Ciphertext<P, DEGREE>) -> Ciphertext<P, DEGREE>`
  - Calls `c0.rescale()` and `c1.rescale()`
  - `logp -= floor(log2(q_last))`, `logq -= floor(log2(q_last))`
- Note: skip wiring into the `PolyRescale` trait (its `rescale_assign(scale_factor: f64)`
  signature doesn't match the RNS drop-channel design); implement directly on `RnsPoly`
  like `mod_drop_last` is.

### Phase 2 — Relinearization audit

The current `mul_ciphertexts` uses naive relin:

```
c0_new = d0 + rk.b * d2
c1_new = d1 + rk.a * d2
```

Mathematically correct (`b + a·s ≈ s²`), but relin noise = `d2 · e_rk`. Rough estimate for
N=16, Δ=2^30, three 31-bit primes: relin noise per decoded slot ≈ 1.5e-8 — comparable to
encryption noise. May be acceptable for the toy setting.

**If tests show relin noise is too large**, implement gadget decomposition:
- Generate relin key at extended basis (L regular + K special primes), encrypting `P·s²`
  where `P = prod(special primes)`
- Relinearization: extend `d2` to special-prime channels, multiply vs relin key, scale down
  by P, drop special channels
- This requires refactoring `RelinearizationKey` and `mul_ciphertexts`

**Recommended order:** implement Phase 1 first, run Phase 3 example, measure actual noise,
then decide if gadget decomp is needed.

### Phase 3 — Example `examples/encrypt_mul.rs` ✅

Full pipeline: encode → encrypt → mul_ciphertexts_gadget → rescale_ciphertext → decrypt → decode

Implemented and passing. Key implementation notes:
- `RnsGadgetRelinKey<DEGREE>` in `engine.rs`, one (a_i, b_i) pair per channel
- `rescale_ciphertext` shares one `Arc<RnsBasis>` for c0 and c1 (required for `ptr_eq` checks)
- Decryption after rescale requires reducing `sk` to the same 3-prime basis
- `bits_dropped = bit_length(q_last)` (NOT `floor(log2)`); floor gives factor-of-2 error
- 4 × 31-bit primes, `SCALE_BITS=30`; after rescale `logp=29`, error ≈ 5×10⁻⁶

## Architecture decisions (standing)

- Single `RnsPoly<const DEGREE: usize>` with `in_ntt_domain: bool`
- Single `Plaintext<P, N>` in `crypto::types`; encoder returns it directly
- `logq` = Σ floor(log2(qᵢ)) = bit-size of composite modulus Q
- After rescale: drop last RNS channel, logq -= floor(log2(q_last)), logp unchanged
- Channel i strictly aligned with modulus qᵢ and NTT table i
