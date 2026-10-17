# Toy HEAAN-CKKS Architecture Overview

## Project Goals

Educational implementation of CKKS (Cheon-Kim-Kim-Song) homomorphic encryption, also known as HEAAN (Homomorphic Encryption for Arithmetic of Approximate Numbers). The project now concentrates exclusively on the RNS-NTT backend so every example, test, and spec builds intuition for that concrete execution model.

**⚠️ Educational Purpose Only**: Not optimized for production use.

## Backend Architecture

All polynomial arithmetic is implemented with a single `RnsNttPoly` backend:

- **Ring**: `Z_Q[X]/(X^N + 1)` where `N` is a power of two; coefficients live in an RNS basis of NTT-friendly primes.
- **Residues**: Each prime channel stays in NTT form for fast convolution; CRT reconstruction is available for encode/decode debugging.
- **Transform Support**: Bluestein/Cooley-Tukey style tables are cached per modulus to keep NTT and inverse NTT consistent across the stack.
- Variants for benches/tests may tweak modulus chains or sampling strategies, but the runtime representation is always RNS-NTT.

## CKKS Operations Status

| Operation | Status | Notes |
|-----------|--------|-------|
| Encoding/Decoding | ✅ | Textbook encoder/decoder on RNS-NTT with Vandermonde transforms. |
| Key Generation | ✅ | RLWE keys; gadget relin/rotation keys per level. |
| Encryption/Decryption | ✅ | Matches textbook parameters for toy degrees. |
| Homomorphic Addition | ✅ | Component-wise add. |
| Homomorphic Multiplication | ✅ | `ct.mul(&other, &rlk)` (gadget relin), `mul_plain`, `mul_const`; rescale is explicit. |
| Rescaling | ✅ | `Ciphertext::rescale` divides by the last prime; `mod_switch_to`/`match_level` drop primes. |
| Rotation | ✅ | Gadget rotation keys, slot rotations. |
| Bootstrapping | ❌ | Out of scope for now. |

## Architecture Design

### Backend trait (`src/rings/traits.rs`)

`PolyRing<N>` is the only view `ckks/` has of polynomials:

- **Arithmetic:** `AddAssign`, `SubAssign`, `MulAssign` (negacyclic), `Neg`, `zero`.
- **Sampling:** `sample_uniform`, `sample_gaussian(std_dev)`, `sample_ternary(hamming_weight)`.
- **Automorphism:** `automorphism(k)`; `rotate_slots` is a provided method on top.
- **Modulus chain:** `level`, `modulus_bits`, `rescale_bits`, `rescale`, `mod_switch_to`.
- **Encoding boundary:** `from_rounded_f64`, `to_centered_f64`; `from_coeffs` / `to_coeffs`.
- **Gadget:** `gadget_len`, `decompose(digit_bits)`, `gadget_term(k)`, with `sum_k decompose(x)[k] * t.gadget_term(k) == x * t`.
- Errors are the backend-agnostic `RingError` (`InvalidLevel`, `CoefficientTooLarge`).

### RNS-NTT Implementation

- **Basis Management**: `src/rings/backends/rns_ntt` owns modulus-chain construction, CRT reduction, and per-prime tables built from `src/math/primes.rs`.
- **NTT form at rest**: see the gotcha below; the form is private to `RnsPoly`.
- **Sampling**: Error and ternary samplers operate channel-wise to stay compatible with the NTT layout.
- **Diagnostics**: Helper utilities expose per-prime residues and CRT reconstructions for encode/decode testing.

### Encoding Strategies

- **CkksEncoder**: canonical Vandermonde/special-DFT encode/decode (spec archived in `specs/archive/textbook_encode_decode.md`). Stages coefficients in i64, see gotchas.

### Key Management

#### Secret Keys

- Ternary polynomials: coefficients in {-1, 0, 1}
- RLWE-style representation: `(a, b)` where `b + a·s ≈ 0`
- Backend-agnostic generation

#### Public Keys

- RLWE encryption of secret key
- Standard form: `(a, b)` where `b + a·s ≈ noise`

#### Relinearization and Rotation Keys

- `KeySwitchingKey<P>` (fhetextbook D-3.8, OpenFHE `EvalKey`): `b[k] + a[k]*s ~ s_old.gadget_term(k)`; wrapped by `RelinearizationKey` (`s^2`) and `RotationKey` (`s(X^(5^k))`)
- BV-style base-2^DIGIT_BITS digits; level-bound

## Key Design Decisions

### 1. Educational Focus Over Performance

- Clear, readable implementations preferred over optimization
- Single backend keeps the codebase approachable
- No external crypto libraries for core CKKS math

### 2. Trait-Based Architecture

- Unified API keeps polynomial math decoupled from higher-level CKKS flows
- Engine layer stays agnostic to concrete ring storage
- Leaves the door open for future specialized RNS-NTT variants without churn

### 3. Modulus Chain Planning

- RNS-NTT: Designed for multi-level modulus switching and NTT compatibility using primes from `src/math/primes.rs`.

### 4. Error Management

- CKKS is approximate - small errors expected
- Typical bounds: <1e-6 (good), <1e-3 (acceptable)
- Need better instrumentation to track noise throughout the single backend

## Current Limitations

- No bootstrapping; circuit depth is bounded by the prime chain.
- Encode/decode is O(N^2) (special DFT).
- Encoder i64 staging caps |slot| at ~4 at scale 2^61.
- Demo parameters (N = 8192, 61-bit primes) are below 128-bit security; see `specs/demo.md`.

## Target architecture (agreed 2026-10-13, M2.5)

Experimental library: ring backends are pluggable and the scheme sees the ring only through one trait. Concepts borrowed from fhe.rs and tfhe-rs.

- **Layers:**
  1. `math/`: modular arithmetic, primes, distributions, NTT tables.
  2. `rings/`: the backend trait, the level context chain, `backends/{rns_ntt, reference (M2.6)}`.
  3. `ckks/`, generic over the backend: `params`, `encoding`, `plaintext`, `ciphertext`, `keys/` (secret, public, key_switching, relin, rotation), `ops`, `evaluator` (per-level key cache).
  4. `ckks/algorithms`: `sum_slots`, Chebyshev; built only on `ops`. Applications (sensor std-dev) live in `examples/`.
- **Backend trait** (ring `Z_Q[X]/(X^N + 1)` at a level; degree is a compile-time `const N`):
  - done in 0011 as `PolyRing` (see Backend trait above); the NTT form is private to the backend. Reference for current practice: OpenFHE (SEAL is unmaintained).
- **Context chain** (fhe.rs): every level's context is built once and shared (`next()`); polys at the same level share it, so there are no fresh-basis mismatches.
- **Keys** (fhe.rs): one `KeySwitchingKey`, wrapped by the relin key and the rotation key (done in 0012).
- **Ciphertext** (fhe.rs): level changes are methods (`rescale`, `mod_switch_to`); Add/Sub/Neg are `std::ops` impls.
- **Style** (tfhe-rs): data types separate from algorithms; typed parameter newtypes (`DecompositionBaseLog`, `LevelCount`) where raw integers are ambiguous.

## Ops and algorithms API

`src/ckks/ops.rs` holds the ops: `&a + &b`, `&a - &b`, `-&a`, `&ct + &pt` (`std::ops`), and the `Ciphertext` methods `mul`, `mul_plain`, `mul_const`, `add_const`, `zero_like`, `rotate`; nothing rescales on its own. `src/ckks/evaluator.rs` holds `LeveledEvaluator` (per-level key cache) and `generate_rotation_key_set`. `src/ckks/algorithms/` holds `sum_slots` (+ `required_rotation_offsets`) and `chebyshev` (`eval_chebyshev_encrypted`, the plaintext references `fit_chebyshev` / `eval_chebyshev_clenshaw`, the domain maps); algorithms use only ops, level methods and the evaluator. Application code (the sensor std-dev demo: `SensorContract`, `encrypted_std_dev`) lives in `examples/sensor_std_dev/`; `[[example]] test = true` in Cargo.toml runs its tests under `cargo test`.

Direction and active work: `specs/roadmap.md`, `specs/ctx.md`, `specs/tasks/`.

# Gotchas (durable knowledge)

Verified against fhetextbook D-3.1/D-3.6/D-3.7 during the Phase A work (2026-08-17):

- **`RnsPoly` rests in NTT form** (0011): ring ops are pointwise; coefficient form appears only inside encode/sample, decode, `decompose`, `rescale` (last channel) and `automorphism`. A private flag allows coefficient form internally (future mod raise); mixed inputs convert to NTT. The switch changed no output bits.
- **Operations require `Arc::ptr_eq` bases.** `RnsBasis::new` builds the whole level chain once (`next()`), and rescale / mod switch walk it, so everything at the same level of one chain shares one `Arc`. Two chains built by separate `RnsBasis::new` calls are NOT compatible, even with equal primes. `logq` is derived from the basis (floor-based `total_bits`); only `logp` is stored.
- **Encoding limit is the modulus, not i64** (fixed in task 0004): a scaled coefficient must satisfy |coeff| < Q/2 for the plaintext's basis (checked as `2^(total_bits - 1)`; `try_encode` returns `RingError::CoefficientTooLarge` from the backend's `from_rounded_f64`, `encode` panics). The final result must also fit the modulus left after all rescales.
- **The encoder stages coefficients in f64**: precision is ~2^-52 relative to the largest scaled coefficient, so mixing huge and tiny slots costs the tiny ones absolute accuracy. Before task 0004 the staging was i64 (|slot| < 4 at scale 61, silent saturation); the 1/10 pre-scaling in older examples dates from that era.
- **Gadget keys are level-bound**: a relin/rotation key only works at the channel count it was generated for. `LeveledEvaluator` regenerates and caches per level from the master secret key.
- **Every ct×ct mul needs modulus headroom**: the raw product sits at scale ~2^(2·logp); with 61-bit primes the mul must run at ≥ 3 channels or the product wraps (symptom: exact ±1.0 slot error).
- **Exact division by 2 via modular inverse is invalid in CKKS**: odd plaintext integers shift by Q/2 (garbage ~1e92 in slots). Level-free alternatives that DO work: `ct.mul_const(0.5, 1)` (halving) or `ct.mul_const(1.0, k)` (logp lift); exact integer scalars at 0 bits (BFV-style, D-3.6).
- **Sensor range contract** (decided 2026-10-12): the client normalizes readings `x_n = (x - lo) / W` with `W = hi - lo`, fits sqrt on the std envelope `[(std_min/W)^2, (std_max/W)^2]`, and restores `std = std_n * W` after decrypt. Accuracy depends only on `std_max/std_min` and the degree, not on `W`. The code lives in `examples/sensor_std_dev/`.
- **Scalars fill every slot** (task 0003, SEAL-style): encode a scalar with `CkksEncoder::encode_constant`, the exact constant polynomial `round(value * 2^k)`. Spreading a constant over fewer slots gives a non-constant polynomial whose coefficients can round to 0 (`0.5` at 1 bit became the zero polynomial; std came out 5.21 instead of 2.87 degC).
- **`sum_slots` sums windows, not the whole vector**: slot i gets `x_i + ... + x_{i+slots-1}` with rotations wrapping over all N/2 slots. If the data fills fewer than N/2 slots, tile it (slot i = `x_{i mod slots}`, slots dividing N/2) or only slot 0 holds the full sum.
- **Noise at N = 8192** (61-bit primes, DIGIT_BITS = 20): the std-dev pipeline's variance error is ~5e-8 normalized (5e-4 degC^2 at W = 100), key-switch noise; toy-N errors (~1e-12) do not carry over.
- **DIGIT_BITS noise tradeoff**: key-switch noise grows ~2^DIGIT_BITS·σ; 15→20 raised std_dev_8 variance error 2.5e-13 → 7.1e-12 (32x, one fewer gadget digit per 61-bit channel).

