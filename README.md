# Toy HEAAN-CKKS

A Rust toy implementation of the CKKS (Cheon-Kim-Kim-Song) homomorphic
encryption scheme, also known as HEAAN (Homomorphic Encryption for Arithmetic
of Approximate Numbers). This project is meant for educational and experimental
purposes, focusing on the internal mechanics of approximate homomorphic
encryption.


```bash
> cargo run --example ckks_sum

Values 1: [1.5, 2.5, 3.5, 4.5]
Values 2: [0.5, 1.0, 1.5, 2.0]
Encoded coeffs 1: [3019898880, -178298470, -738197504, 580951654, -738197504, 580951654, -738197504, -178298470]
Encoded coeffs 2: [1275068416, -122703667, -335544320, 256921395, -335544320, 256921395, -335544320, -122703667]
Poly 1: RnsPoly<8>[[3052, 3576, 1486]*x^0, [7005, 2622, 597]*x^1, [6735, 2281, 2420]*x^2, …, [709, 5968, 3776]*x^5, [6735, 2281, 2420]*x^6, [7005, 2622, 597]*x^7]
Poly 2: RnsPoly<8>[[6160, 4665, 1513]*x^0, [250, 1722, 3255]*x^1, [6620, 4110, 1100]*x^2, …, [4931, 3395, 1998]*x^5, [6620, 4110, 1100]*x^6, [250, 1722, 3255]*x^7]
Encoded to polynomials successfully
Encrypted both values successfully
Performed homomorphic addition
Decrypted the sum successfully

=== CKKS Sum Results ===
Expected sum: [2.0, 3.5, 5.0, 6.5]
Computed sum: [1.9999999767169356, 3.4999999777574473, 5.00000000372529, 6.500000027830488]
Maximum error: 2.78e-8
✅ Success! Error within acceptable bounds
```

## Overview

CKKS is a homomorphic encryption scheme for approximate arithmetic over
encrypted real (floating-point) numbers. It’s designed for privacy-preserving
computations in applications such as machine learning and signal processing.

This crate provides a minimal, step-by-step Rust implementation using schoolbook
methods and native u64 arithmetic with an RNS (Residue Number System) basis.

Key features of this implementation:
- [x] Encoding/decoding of real numbers
- [x] Keys generation (secret and public)
- [x] Polynomial ring operations (Add, Mul, Negate)
- [x] Basic homomorphic operation - Addition
- [ ] Relinearization and rescale
- [ ] NTT-based polynomial multiplication
- [ ] Basic homomorphic operation - Multiplication
- [ ] Bootstrapping / Modulus switching

## Examples

Only ckks_sum is currently functional:

```bash
cargo run --example ckks_sum
```

This demonstrates:
- encoding two real-valued vectors into polynomials
- encrypting them with CKKS
- performing ciphertext addition
- decrypting and decoding the result

## Implementation Notes

- Polynomials are stored in RNS format: `Vec<[u64; DEGREE]>` where each
  `[u64; DEGREE]` represents coefficients mod a prime.
- No BigInts are used — all computation remains in native u64.
- Polynomial multiplication is currently schoolbook, no NTT optimizations yet.
- The ciphertext format includes an optional c2 component for future
  multiplication support but is unused.

## Limitations

This is not a secure or production-ready library. It is a research prototype
with the following constraints:

- No NTT or FFT-based polynomial multiplication
- No modulus switching, relinearization, or noise tracking
- No ciphertext multiplication yet
- No parameter validation or error estimation
- No SIMD-style ciphertext rotations or packing operations

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
