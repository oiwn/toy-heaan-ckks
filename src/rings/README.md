# Polynomial Rings

This module provides polynomial ring implementations for CKKS homomorphic encryption.

## Overview

CKKS operates on polynomials in the quotient ring **R_q = Z_q[X]/(X^N + 1)**, where:
- **N**: Polynomial degree (power of 2, e.g., 4096, 8192, 16384)
- **q**: Coefficient modulus (large prime, e.g., 2^50 bits)
- **X^N + 1**: Cyclotomic polynomial defining the reduction relation X^N = -1

## Polynomial Representation

Polynomials are represented as coefficient arrays:
```
p(X) = cÄ + cÅ∑X + cÇ∑X≤ + ... + c_{N-1}∑X^{N-1}
```

All coefficients are reduced modulo q, and all polynomial operations respect the relation X^N = -1.

## Operations

- **Addition**: Coefficient-wise addition modulo q
- **Multiplication**: Convolution with reduction X^N = -1, coefficients modulo q
- **Negation**: Coefficient-wise negation modulo q

## Implementations

This module provides multiple backend implementations optimized for different use cases:
- **NaivePolyRing**: Simple u64-based implementation for small rings
- **BigIntPolyRing**: Arbitrary-precision implementation for large moduli

See `traits.rs` for the trait definitions and `backends/` for concrete implementations.
