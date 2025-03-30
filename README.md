# Toy HEAAN-CKKS

A Rust implementation of the CKKS (Cheon-Kim-Kim-Song) homomorphic encryption
scheme, also known as HEAAN (Homomorphic Encryption for Arithmetic of Approximate
Numbers). This educational implementation demonstrates the core concepts of
approximate homomorphic encryption.

```text
Running `target/debug/examples/full_cycle`
Original values 1: [1.5, 2.5, 3.5, 4.5]
Original values 2: [0.5, 1.0, 1.5, 2.0]
Expected sum: [2.0, 3.5, 5.0, 6.5]
Decrypted result: [2.000000011175871, 3.4999999979189345, 4.99999999627471, 6.4999999834546145]
```

## Overview

CKKS is a homomorphic encryption scheme designed for performing computations
on encrypted real (floating-point) numbers. Unlike other HE schemes that work
with integers or bits, CKKS is specialized for approximate arithmetic, making it
suitable for machine learning and statistical analysis on encrypted data.

Key features of this implementation:
- Polynomial ring operations
- Encoding/decoding of real numbers
- Key generation (secret and public)
- Basic homomorphic operations (addition)
- FFT-based encoding for efficient slot packing

## Getting Started

### Prerequisites

- Rust 1.85+

### Installation

```bash
cargo run --example full_cycle 
```

## Limitations

This is an educational implementation meant for learning purposes:

- Not optimized for performance (no RNS, no bootstrapping)
- Limited security analysis
- Not suitable for production use
- Only basic homomorphic operations implemented
- No parameter selection guidance

## Mathematical Background

The CKKS scheme operates on the polynomial ring R = `Z[X]/(X^n + 1)` where n is
a power of 2. It encodes vectors of complex numbers into polynomials such that
when evaluated at certain roots of unity, the polynomial yields the scaled input
values.

Key components:
1. **Encoding**: Maps a vector of real/complex numbers to a polynomial
2. **Scaling**: Uses fixed-point representation with a scaling factor Δ = 2^scale_bits
3. **Encryption**: Adds randomized noise to ensure security
4. **Homomorphic Operations**: Preserves encrypted arithmetic
5. **Rescaling**: Manages the growth of coefficients and noise

## References

- [Original HEAAN paper](https://eprint.iacr.org/2016/421.pdf) - Cheon, Kim, Kim, Song
- [Homomorphic Encryption for Arithmetic of Approximate Numbers](https://link.springer.com/chapter/10.1007/978-3-319-70694-8_15)
- [A Guide to Fully Homomorphic Encryption](https://eprint.iacr.org/2015/1192.pdf)

## License

[MIT License](LICENSE)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
