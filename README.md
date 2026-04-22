# Toy HEAAN-CKKS

Educational Rust implementation of CKKS, also known as HEAAN: homomorphic encryption for approximate arithmetic.

This repository is for learning the mechanics of CKKS.  It is intentionally small, explicit, and instrumented through examples.  It is not production cryptography.

## Current Focus

The codebase currently centers on one concrete backend:

- `RnsPoly` / `RnsBasis` in `src/rings/backends/rns_ntt/`
- arithmetic over `Z_Q[X] / (X^N + 1)`
- an RNS modulus chain made from NTT-friendly primes
- plaintext/ciphertext operations that keep polynomial channels in NTT form
  where possible

Older README notes about separate naive, BigInt, and non-NTT backends are stale.
The active path is the RNS-NTT backend used by the examples and tests.

## Build And Test

```bash
cargo build
cargo test
cargo clippy --all-targets -- -D warnings
cargo fmt --check
```

Useful stricter lint pass:

```bash
cargo clippy -- -D clippy::pedantic -D clippy::nursery
```

Benchmarks currently available:

```bash
cargo bench --bench basic
cargo bench --bench primes
```

## Examples

Run examples with:

```bash
cargo run --example <name>
```

Current examples:

| Example | Purpose |
|---------|---------|
| `encode_decode` | CKKS encode/decode behavior and scale precision. |
| `keys` | Secret, public, and relinearization key generation checks. |
| `encrypt_add` | End-to-end encode, encrypt, homomorphic add, decrypt, decode. |
| `encrypt_mul` | Homomorphic multiplication with gadget relinearization and rescale. |
| `rotation_demo` | Slot rotation with gadget rotation keys. |
| `rotation_stress` | Rotation checks across a wider set of offsets. |
| `horner_chain` | Multi-level affine Horner-style chain over many slots. |
| `std_dev_8` | Documented scaffold for encrypted standard deviation over 8 slots. |

Good starting points:

```bash
cargo run --example encode_decode
cargo run --example encrypt_add
cargo run --example encrypt_mul
cargo run --example rotation_demo
```

## What Works

| Area | Status | Notes |
|------|--------|-------|
| RNS-NTT basis generation | Working | Uses generated primes congruent to `1 mod 2N`. |
| CKKS encoding/decoding | Working for demos | Approximate real slots with scale `2^scale_bits`. |
| Secret/public keys | Working | Ternary secret keys and RLWE-style public keys. |
| Encryption/decryption | Working | Used by all end-to-end examples. |
| Homomorphic addition | Working | Component-wise ciphertext addition. |
| Homomorphic multiplication | Working in examples | Uses gadget relinearization, then rescale. |
| Rescale | Working for current examples | Drops one RNS prime and lowers ciphertext scale. |
| Slot rotation | Working in examples | Uses gadget rotation keys. |
| Bootstrapping | Not implemented | Out of scope for now. |

This is still a toy implementation.  The examples are the source of truth for
which parameter sets currently behave well.

## Project Layout

```text
src/crypto/              CKKS engine, operations, ciphertext types
src/encoding/            CKKS encoders and FFT helpers
src/keys/                Secret, public, relinearization, rotation keys
src/math/                Prime generation, sampling, utility math
src/rings/backends/      RNS-NTT polynomial backend
src/rings/traits.rs      Backend traits used by the engine
examples/                Runnable demonstrations
tests/                   Integration tests
specs/                  Design notes and active task context
```

Important specs:

- `specs/overview.md` - architecture notes
- `specs/parameters.md` - parameter planning notes
- `specs/ctx.md` - active task context; keep this updated while developing

## Minimal Engine Shape

Most examples follow this structure:

```rust
use std::sync::Arc;

use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::crypto::engine::{CkksEngine, CkksParams};
use toy_heaan_ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

const N: usize = 16;
const SCALE_BITS: u32 = 30;

let primes = generate_primes(31, 3, N as u64);
let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));
let logq = basis.total_bits();

let engine = CkksEngine::<RnsPoly<N>, N>::new(
    basis.clone(),
    CkksParams {
        error_variance: 3.2,
        hamming_weight: N / 2,
        scale_bits: SCALE_BITS,
    },
);

let encoder = CkksEncoder::<N>::new(SCALE_BITS);
let mut rng = ChaCha20Rng::seed_from_u64(42);
let sk = engine.generate_secret_key(&mut rng).expect("secret key");
let pk = engine.generate_public_key(&sk, &mut rng).expect("public key");

let values = [1.0, 2.0, 3.0, 4.0];
let pt = encoder.encode(&values, basis);
let ct = engine.encrypt(&pt, &pk, logq, &mut rng);
let decoded = encoder.decode(&CkksEngine::<RnsPoly<N>, N>::decrypt(&ct, &sk));
```

For multiplication, see `examples/encrypt_mul.rs`.  For rotations, see
`examples/rotation_demo.rs`.

## Known Limitations

- No bootstrapping.
- No production security parameter selection.
- Small toy parameters are common in examples.
- Noise tracking is mostly manual and example-driven.
- Some naming around error variance vs. standard deviation is still being
  audited.
- The code is designed for education, not side-channel resistance,
  performance, or deployment.

## References

- [HEAAN paper](https://eprint.iacr.org/2016/421.pdf)
- [CKKS paper](https://link.springer.com/chapter/10.1007/978-3-319-70694-8_15)
- [FHE Textbook](https://fhetextbook.github.io)
- [Microsoft SEAL](https://github.com/microsoft/SEAL)
- [Lattigo](https://github.com/tuneinsight/lattigo)

## License

[MIT License](LICENSE)
