# Toy HEAAN-CKKS

Educational Rust implementation of CKKS (HEAAN): homomorphic encryption for
approximate arithmetic on encrypted real and complex numbers. Small, explicit
code; every step of the scheme has a runnable example. Not production
cryptography: parameters are toy-sized and there is no security claim.

## Quick start

```bash
cargo run --example encrypt_add             # encode, encrypt, add, decrypt
cargo run --release --example bootstrap_demo
```

## Examples

`cargo run --example <name>`; add `--release` for bootstrap examples.

| Step | Example | Shows |
|---|---|---|
| Encoding | `encode_decode` | Slots to polynomial and back, scale precision |
| Keys | `keys` | Secret, public and relinearization keys |
| Add | `encrypt_add` | Encrypt, add, decrypt, decode |
| Multiply | `encrypt_mul` | Multiply, relinearize, rescale |
| Rotate | `rotation_demo`, `rotation_stress` | Slot rotation; noise over many rotations |
| Depth | `horner_chain` | Multi-level polynomial over 4096 slots |
| Bootstrap | `bootstrap_demo` | Refresh, multiply, refresh again (N=16/64) |
| Bootstrap | `bootstrap_large` | Large degrees up to N=8192 |
| Bootstrap | `recurrence_bootstrap_demo` | 32 rounds of a recurrence, with and without refresh |
| Application | `sensor_std_dev`, `sensor_std_plain`, `std_dev_8` | Encrypted standard deviation |
| Approximation | `sqrt_chebyshev_plain` | Chebyshev approximation of `sqrt(x)` |
| Arithmetic | `modmul_methods` | u128, Barrett and Shoup modular multiplication |

## Develop

```bash
cargo test
cargo clippy --all-targets -- -D warnings
cargo bench          # see specs/benchmarks.md
```

Design notes: `specs/overview.md`. Backend: RNS-NTT over `Z_Q[X]/(X^N+1)`
(`src/rings/backends/rns_ntt/`), behind the `PolyRing` trait.

## References

[CKKS paper](https://eprint.iacr.org/2016/421.pdf) ·
[FHE Textbook](https://fhetextbook.github.io) ·
[OpenFHE](https://github.com/openfheorg/openfhe-development) ·
[Lattigo](https://github.com/tuneinsight/lattigo)

[MIT License](LICENSE)
