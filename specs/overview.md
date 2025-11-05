# Toy HEAAN-CKKS Architecture Overview

## Project Goals

Educational implementation of CKKS (Cheon-Kim-Kim-Song) homomorphic encryption, also known as HEAAN (Homomorphic Encryption for Arithmetic of Approximate Numbers). Focus on demonstrating internal mechanics through multiple backend implementations.

**⚠️ Educational Purpose Only**: Not optimized for production use.

## Backend Comparison

| Backend | Modulus Support | Arithmetic | Precision | Status |
|---------|----------------|-------------|-----------|---------|
| **Naive** | Single `u64` modulus | Schoolbook O(n²) | Limited by f64 | ✅ Complete |
| **RNS** | Multi-prime CRT | Schoolbook per channel | High | ✅ Complete |
| **BigInt** | `crypto-bigint::U256` | Schoolbook O(n²) | Very High | ✅ Complete |
| **NTT** | Multi-prime + NTT | NTT-optimized | High | 🚧 Experimental |

## CKKS Operations Status

| Operation | Naive | RNS | BigInt | NTT |
|-----------|-------|-----|--------|-----|
| Encoding/Decoding | ✅ | ✅ | ✅ | ✅ |
| Key Generation | ✅ | ✅ | ✅ | ✅ |
| Encryption/Decryption | ✅ | ✅ | ✅ | ✅ |
| Homomorphic Addition | ✅ | ✅ | ✅ | ✅ |
| Homomorphic Multiplication | ❌ | ❌ | ❌ | ❌ |
| Rescaling | ❌ | ❌ | Partial | ❌ |
| Bootstrapping | ❌ | ❌ | ❌ | ❌ |

## Architecture Design

### Core Traits (`src/rings/traits.rs`)
- **PolyRing**: Basic arithmetic operations (add, mul, neg)
- **PolySampler**: Random polynomial generation with error distributions
- **PolyRescale**: Modulus switching and rescaling operations

### Backend Implementations

#### Naive Backend
- **Ring**: `Z_q[X]/(X^N + 1)` with single `u64` modulus
- **Arithmetic**: Schoolbook negacyclic multiplication
- **Use Case**: Simple demonstrations, educational clarity
- **Limitations**: Modulus size limited by 64-bit arithmetic

#### RNS Backend  
- **Ring**: Chinese Remainder Theorem across multiple primes
- **Arithmetic**: Schoolbook per CRT channel
- **Use Case**: High precision, foundation for NTT
- **Advantages**: Native `u64` ops, parallelizable, scalable

#### BigInt Backend
- **Ring**: `crypto-bigint::U256` arithmetic
- **Arithmetic**: Schoolbook negacyclic with large moduli
- **Use Case**: Very large moduli, HEAAN-style operations
- **Features**: Extended modulus flow, high-precision encoder

### Encoding Strategies

#### RustFftEncoder
- Standard FFT-based canonical encoding
- Uses `f64` complex arithmetic
- Symmetric layout, slot packing (N/2 complex numbers)
- Best for: Simple use cases, standard CKKS

#### BigIntEncoder  
- HEAAN-style special FFT with `rotGroup` and `ksi_pows`
- Precomputed tables for performance
- High-precision using Malachite integers
- Best for: HEAAN compatibility, large scale_bits

#### NaiveEncoder
- Simple coefficient-based encoding
- Minimal preprocessing
- Best for: Testing, educational purposes

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
- Multiple backends to explore different approaches
- No external crypto libraries for core operations

### 2. Trait-Based Architecture
- Unified API across all backends
- Backend-agnostic engine layer
- Easy to add new implementations

### 3. Modulus Chain Planning
- RNS: Designed for multi-level modulus switching
- BigInt: Extended modulus `qQ` flow for HEAAN compatibility
- Naive: Single modulus with scale tracking challenges

### 4. Error Management
- CKKS is approximate - small errors expected
- Typical bounds: <1e-6 (good), <1e-3 (acceptable)
- Noise tracking varies by backend

## Current Limitations

### Across All Backends
- No complete homomorphic multiplication
- No bootstrapping implementation
- Schoolbook multiplication only (O(n²))

### Backend-Specific
- **Naive**: Limited modulus size, scale tracking issues
- **RNS**: Complex setup, requires modulus chain planning
- **BigInt**: Performance overhead of big integer ops
- **NTT**: Experimental, not fully optimized

## Near-Term Development

1. **Complete Homomorphic Multiplication**
   - Implement relinearization across all backends
   - Add proper rescaling operations
   - Noise budget management

2. **Modulus Switching**
   - Naive: Scale-aware implementation
   - RNS: Zero-noise ModDrop operations
   - BigInt: Extended modulus flow completion

3. **Performance Optimizations**
   - NTT integration for RNS backend
   - SIMD/AVX optimizations where applicable

## Future Goals

- Bootstrapping implementation
- Security parameter analysis
- Production-ready optimizations
- Advanced features (rotation, packing)