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
cargo run --example ckks_sum

# Try different backends
cargo run --example ckks_naive
cargo run --example ckks_rns
cargo run --example ckks_bigint
```

## Architecture Overview

This implementation provides **multiple polynomial arithmetic backends** to
explore different approaches to CKKS operations:

### Backend Options

- **`NaivePolyRing`** - Simple single-modulus arithmetic using native `u64`
- **`RnsPolyRing`** - Multi-prime RNS (Residue Number System) for improved precision
- **`PolyRingU256`** - Large modulus support using `crypto-bigint::U256`
- **`RnsNttPolyRing`** - NTT-optimized RNS backend (experimental)

### Key Features by Backend

| Feature | Naive | RNS | U256 | NTT |
|---------|-------|-----|------|-----|
| Basic operations (Add, Mul) | ✅ | ✅ | ✅ | ✅ |
| Schoolbook multiplication | ✅ | ✅ | ✅ | ✅ |
| Multi-prime arithmetic | ❌ | ✅ | ❌ | ✅ |
| Large modulus support | ❌ | ✅ | ✅ | ✅ |
| NTT optimization | ❌ | ❌ | ❌ | 🔄 |

### CKKS Operations Status

- [x] **Encoding/Decoding** - Float arrays ↔ polynomial coefficients
- [x] **Key Generation** - Secret keys, public keys, relinearization keys
- [x] **Encryption/Decryption** - Full CKKS encryption scheme
- [x] **Homomorphic Addition** - Encrypted addition with all backends
- [ ] **Homomorphic Multiplication** - Partial implementation, needs relinearization
- [ ] **Rescaling** - Noise management after multiplication
- [ ] **Bootstrapping** - Ciphertext refresh (advanced feature)

## Examples

### Basic Encryption/Decryption with RNS Backend

```bash
cargo run --example ckks_naive
```

**Output:**
```
❯ cargo run --example ckks_naive                                                                                                                                          12:03:04 [5/4925]

   Compiling toy-heaan-ckks v0.1.0 (/Users/alexch/code/toy-heaan-ckks)
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 0.26s
     Running `target/debug/examples/ckks_naive`
🔐 CKKS Abstract API Demo with Encryption                                                    
✅ Engine configured with builder pattern

🔑 Generating keys...
✅ Secret and public keys generated: 
SecretKey { poly: NaivePolyRing { coeffs: [1, 0, 1, 0, 0, 1, 0, 741507920154517876], context: 741507920154517877 } }
PublicKey { a: NaivePolyRing { coeffs: [693672716415451047, 302556833492949680, 602596531573667421, 252165852101715923, 147456709156760292, 160315845728327181, 23828265504687498, 17872185
407316598], context: 741507920154517877 }, b: NaivePolyRing { coeffs: [21272487852520571, 1683449651977668, 94896585946528545, 63156790907779483, 590518939257597456, 353033160558854108, 2
49793926592803811, 654396073860657723], context: 741507920154517877 } }

📊 Input data: [1.5, 2.5, 3.5]

🔢 Encoding values...
✅ Values encoded to plaintext: Plaintext { poly: NaivePolyRing { coeffs: [1855425871872, 692078510204, 741507164240273781, 741507640392868089, 481036337152, 741507640392868089, 741507164
240273781, 692078510204], context: 741507920154517877 }, scale: 1099511627776.0 }

🔒 Encrypting plaintext...
✅ Plaintext encrypted to ciphertext: Ciphertext { c0: NaivePolyRing { coeffs: [11589659624651541, 32966279333054012, 246331520150206956, 321734295273184328, 620973310206268142, 411192794
789119066, 439181144981773282, 694361773217369658], context: 741507920154517877 }, c1: NaivePolyRing { coeffs: [418437777927414003, 141321956329697559, 528580068367660287, 371703167324472
070, 106392592596933299, 82037449374961662, 199644371131419166, 596475019354968358], context: 741507920154517877 }, scale: 1099511627776.0 }

🔓 Decrypting ciphertext...
✅ Ciphertext decrypted back to plaintext: Plaintext { poly: NaivePolyRing { coeffs: [1855425871867, 692078510201, 741507164240273777, 741507640392868090, 481036337155, 741507640392868090
, 741507164240273787, 692078510208], context: 741507920154517877 }, scale: 1099511627776.0 }

🔢 Decoding back to floating-point...
✅ Plaintext decoded: [1.5000000000027285, 2.499999999991649, 3.499999999996362, -6.200817637136424e-12]

📊 Results:
  Original: [1.5, 2.5, 3.5]
  Decoded:  [1.5000000000027285, 2.499999999991649, 3.499999999996362]
  Max error: 8.35e-12
🎉 Success! Full CKKS encrypt/decrypt pipeline works!

➕ Testing homomorphic addition...
Second input: [0.5, 1.0, 1.5]
Expected sum: [2.0, 3.5, 5.0]
Computed sum: [2.000000000032742, 3.4999999999793583, 4.999999999994543]
Sum error: 3.27e-11
🎉 Homomorphic addition works perfectly!

💡 Full CKKS pipeline completed:
  1. ✅ Key generation (secret + public keys)
  2. ✅ Encoding (float → plaintext)
  3. ✅ Encryption (plaintext → ciphertext)
  4. ✅ Homomorphic operations (encrypted addition)
  5. ✅ Decryption (ciphertext → plaintext)
  6. ✅ Decoding (plaintext → float)
```

### Creating Different Backends

```rust
// Naive backend - simple single modulus
let engine = CkksEngine::<NaivePolyRing<8>, 8>::builder()
    .error_variance(3.2)
    .hamming_weight(4)
    .build_naive(modulus, scale_bits)?;

// RNS backend - multiple primes for better precision  
let rns_basis = Arc::new(RnsBasisBuilder::new(8)
    .with_prime_bits(vec![17, 19, 23])
    .build()?);
let engine = CkksEngine::<RnsPolyRing<8>, 8>::builder()
    .build_rns(rns_basis, scale_bits)?;

// BigInt backend - large modulus support
let modulus = NonZero::new(U256::from_u128(large_prime))?;
let engine = CkksEngine::<PolyRingU256<8>, 8>::builder()
    .build_bigint_u256(modulus, scale_bits)?;
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
assert!(error < 1e-6);  // Good pre

## Mathematical Background

The CKKS scheme operates on the polynomial ring R = `Z[X]/(X^n + 1)` where n is
a power of 2. It encodes vectors of complex numbers into polynomials such that
when evaluated at certain roots of unity, the polynomial yields the scaled input
values.

## References

- [Original HEAAN paper](https://eprint.iacr.org/2016/421.pdf) - Cheon, Kim, Kim, Song
- [Homomorphic Encryption for Arithmetic of Approximate Numbers](https://link.springer.com/chapter/10.1007/978-3-319-70694-8_15)
- [A Guide to Fully Homomorphic Encryption](https://eprint.iacr.org/2015/1192.pdf)

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
