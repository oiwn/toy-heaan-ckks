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
| Encoding/Decoding | 🚧 Reworking | Textbook encoder/decoder is being ported onto RNS-NTT with Vandermonde transforms. |
| Key Generation | ✅ | RLWE-style keys backed by `RnsNttPoly`. |
| Encryption/Decryption | ✅ | Matches textbook parameters for toy degrees. |
| Homomorphic Addition | ✅ | Component-wise add in NTT domain. |
| Homomorphic Multiplication | ❌ | Requires relinearization + rescale work. |
| Rescaling | ❌ | Not yet implemented for the consolidated backend. |
| Bootstrapping | ❌ | Out of scope for now. |

## Architecture Design

### Core Traits (`src/rings/traits.rs`)
- **PolyRing**: Basic arithmetic operations (add, mul, neg)
- **PolySampler**: Random polynomial generation with error distributions
- **PolyRescale**: Modulus switching and rescaling operations

### RNS-NTT Implementation
- **Basis Management**: `src/rings/backends/rns` owns modulus-chain construction, CRT reduction, and per-prime tables built from `src/math/primes.rs`.
- **NTT Domain Storage**: Polynomials are kept in NTT form by default and lazily converted back when needed (e.g., decoding, CRT inspection).
- **Sampling**: Error and ternary samplers operate channel-wise to stay compatible with the NTT layout.
- **Diagnostics**: Helper utilities expose per-prime residues and CRT reconstructions for encode/decode testing.

### Encoding Strategies

- **TextbookEncoder (in progress)**: Implements the canonical Vandermonde/special FFT encode/decode path using the specification in `specs/textbook_encode_decode.md`. Targets real/complex inputs, Δ-scaling, and centered coefficient reconstruction.
- **RustFftEncoder (legacy)**: Remains temporarily for comparison but will be phased out once the textbook path is stable.

### Key Management

#### Secret Keys
- Ternary polynomials: coefficients in {-1, 0, 1}
- RLWE-style representation: `(a, b)` where `b + a·s ≈ 0`
- Backend-agnostic generation

#### Public Keys
- RLWE encryption of secret key
- Standard form: `(a, b)` where `b + a·s ≈ noise`

#### Relinearization Keys
- For ciphertext multiplication (incomplete)
- Key switching: `(a, b)` with `b + a·s ≈ s²`
- Uses gadget decomposition

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

### System-Wide
- No complete homomorphic multiplication
- No bootstrapping implementation
- Schoolbook multiplication only (O(n²))

### Backend-Specific
- Consolidated backend still lacks rescale/mul implementations.


# RNS-NTT params for toy tests

- `logq` = 60 <<< module bits for plaintext
- `log_small_q` = 30 bit <<< rns module also known as "scale factor"
- `slots` = 4 <<< number of elements in message vector
- `DEGREE` = 8 <<< slots * 2 (polynomial degree)
- `L` = 7 <<< number of modules in RNS
- `K` = 2 <<< additional modules for key switch


## Near-Term Development

1. **Textbook Encode/Decode**
   - Implement `TextbookEncodingParams`, conjugate-slot builders, and Vandermonde transforms.
   - Integrate with `RnsNttPoly` packing and CRT diagnostics.
2. **Homomorphic Multiplication Prep**
   - Finish relinearization/scaling design for the single backend.
   - Add instrumentation for noise growth tracking.
3. **Performance & Tooling**
   - Cache transform tables aggressively.
   - Improve debug examples (`encode_debug`) and document recommended parameter sets.

## Future Goals

- Bootstrapping implementation
- Security parameter analysis
- Production-ready optimizations
- Advanced features (rotation, packing)
