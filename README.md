# Toy HEAAN-CKKS

A Rust educational implementation of the CKKS (Cheon-Kim-Kim-Song) homomorphic
encryption scheme, also known as HEAAN (Homomorphic Encryption for Arithmetic
of Approximate Numbers). This project focuses on demonstrating the internal
mechanics of approximate homomorphic encryption through multiple backend
implementations.

**⚠️ Educational Purpose Only**: This is a minimal, non-optimized implementation
designed for learning and experimentation. It should not be used in production
systems.

## Quick Start

```bash
# Run the main working example
cargo run --example rns

# Try related demos
cargo run --example ntt
cargo run --example rns_mod_switch
```

## Architecture Overview

This implementation provides **multiple polynomial arithmetic backends** to
explore different approaches to CKKS operations:

### Backend Options

- **`NaivePolyRing`** - Simple single-modulus arithmetic using native `u64`
- **`RnsNttPoly`** - Channel-first RNS backend that keeps arithmetic in the NTT domain
- **`BigIntPolyRing`** - Large modulus support using `crypto-bigint::U256`
- **`RnsNttPolyRing`** - NTT-optimized RNS backend (experimental)

### Key Features by Backend

| Feature | Naive | RNS | U256 | NTT |
|---------|-------|-----|------|-----|
| Basic operations (Add, Mul, Neg) | ✅ | ✅ | ✅ | ✅ |
| Sampling functions | ✅ | ✅ | ✅ | ✅ |
| Ciphertext multiplication | ❌ | ❌ | ❌ | ❌ |
| Bootstrapping | ❌ | ❌ | ❌ | ❌ |

### CKKS Operations Status

- [x] **Encoding/Decoding** - Float arrays ↔ polynomial coefficients
- [x] **Key Generation** - Secret keys, public keys, relinearization keys
- [x] **Encryption/Decryption** - Full CKKS encryption scheme
- [x] **Homomorphic Addition** - Encrypted addition with all backends
- [ ] **Homomorphic Multiplication** - Partial implementation, needs relinearization
- [ ] **Rescaling** - Noise management after multiplication
- [ ] **Bootstrapping** - Ciphertext refresh (advanced feature)

## Examples

The crate ships a few runnable demos that cover the current CKKS pipeline:

- `rns`: end-to-end encode → encrypt → add → decrypt flow on the RNS-NTT backend (keys, ciphertexts, and plaintexts stay in the transform domain throughout).
- `ntt`: similar flow with extra logging that shows domain transitions for keys and ciphertext components.
- `rns_mod_switch`: demonstrates dropping primes from the modulus chain (modulus switching) while monitoring decode error.

Run any example with `cargo run --example <name>`; for instance:

```bash
cargo run --example rns
```

### Creating the RNS-NTT Engine

```rust
use toy_heaan_ckks::rings::backends::rns::{RnsBasisBuilder, RnsNttPoly};
use toy_heaan_ckks::CkksEngine;
use std::sync::Arc;

let basis = Arc::new(
    RnsBasisBuilder::new(8)
        .with_prime_bits(vec![17, 19, 23])
        .build()?,
);

let engine = CkksEngine::<RnsNttPoly<8>, 8>::builder()
    .error_variance(3.2)
    .hamming_weight(4)
    .build_rns(basis, 40)?;
```

## Implementation Details

### RNS (Residue Number System)

The RNS backend splits large integers across multiple smaller prime moduli:

```
Coefficient: 12345678901234567890
         ↓
RNS form: [coeff mod p1, coeff mod p2, coeff mod p3, ...]
```

**Benefits:**
- Native `u64` arithmetic throughout
- Better precision with multiple primes
- Parallelizable operations
- Foundation for NTT optimizations

### Error Management

CKKS is an **approximate** scheme - small errors are expected:

```rust
// Typical error bounds
assert!(error < 1e-6);  // Good precision
assert!(error < 1e-3);  // Acceptable for educational cases
```

### Polynomial Representation

All backends use the same abstract interface but different internal representations:

```rust
// All implement PolyRing<DEGREE> trait
trait PolyRing<const DEGREE: usize> {
    fn add_assign(&mut self, other: &Self);
    fn mul_assign(&mut self, other: &Self);
    // ... other operations
}
```

### Encoding Options

Multiple encoding strategies are available:

- RustFftEncoder - FFT-based canonical encoding
- NaiveEncoder - Simple coefficient-based encoding
- BigIntEncoder - Large precision encoding for U256 backend

### Benchmarking

Compare backend performance:
```bash
cargo bench --bench end_to_end
```

This runs comprehensive benchmarks across all backends for encrypt/decrypt
cycles and individual operations.

### Known Issues

- Multiplication requires relinearization (work in progress)
- NTT backend is experimental and not fully optimized
- No modulus switching or noise management
- Schoolbook multiplication only (O(n²) complexity)

### Design Choices

- Focus on clarity over performance
- Educational implementation - not production-ready
- Multiple backends to explore different approaches
- No external crypto libraries for core polynomial arithmetic

### Near Term

- Complete homomorphic multiplication with relinearization
- Implement rescaling for noise management
- Optimize NTT-based multiplication
- Add more comprehensive tests

### Future Goals

- Bootstrapping implementation
- Modulus switching
- SIMD/AVX optimizations
- Security parameter analysis


## References

- [Original HEAAN paper](https://eprint.iacr.org/2016/421.pdf) - Cheon, Kim, Kim, Song
- [Homomorphic Encryption for Arithmetic of Approximate Numbers](https://link.springer.com/chapter/10.1007/978-3-319-70694-8_15)
- [A Guide to Fully Homomorphic Encryption](https://eprint.iacr.org/2015/1192.pdf)
- [SEAL Documentation](https://github.com/Microsoft/SEAL)
- [Lattigo](https://github.com/tuneinsight/lattigo)

## License

[MIT License](LICENSE)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request or
create Issue.

Especially for:

- NTT integration
- Efficient multiplication
- Relinearization and modulus switching
- Clean test coverage for keygen and enc/dec edge cases

This is an educational project. Contributions should maintain the focus on:

- Clarity - Code should be readable and well-documented
- Multiple approaches - Show different ways to implement concepts
- Educational value - Explain the "why" not just the "how"
