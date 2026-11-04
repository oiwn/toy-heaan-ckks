# Toy HEAAN-CKKS

Educational Rust implementation of CKKS, also known as HEAAN: homomorphic encryption for approximate arithmetic.

This repository is for learning the mechanics of CKKS.  It is intentionally small, explicit, and instrumented through examples.  It is not production cryptography.

## Current Focus

For a deeper computation, run:

```bash
cargo run --release --example recurrence_bootstrap_demo
cargo run --release --example recurrence_bootstrap_demo -- --degree 64
```

The encrypted driver `u` initializes `x=u`. Each of 32 rounds computes
`x_next=x/2+x*x/8+sign*u/4`, alternating the sign. Library code handles scale
alignment; key generation is confined to setup. The ideal state stays in
[-0.625,0.875], with a state derivative at most0.75 in magnitude.

Three runs compare the same inputs: 18 primes with bootstrap (32 rounds,
three refreshes), 18 without (stops before round18), and 34 without (32 rounds).
The output reports levels, error at every round, exact initial modulus bit
lengths, refresh events and separate setup/evaluation/diagnostic timings.
Client diagnostics decrypt checkpoints; evaluation receives only public keys.
Accuracy must stay within1e-3. This demonstrates a shorter modulus chain,
not a measured total-memory saving. All parameters remain educational.

Educational complex CKKS bootstrapping is available through
`ckks::algorithms::bootstrap::{BootstrapParameters, BootstrapPlan, BootstrapKeys}`.
Run `cargo run --release --example bootstrap_demo` for encrypted computation,
refresh, multiplication and a second refresh at degrees 16 and 64.

The supported toy profile uses 18 NTT primes within 2^20 of 2^61, input/output
scale 2^44, degree-31 cosine and four double-angle steps. Bootstrap accepts any
level in its context chain and a positive scale no greater than 2^44; it
normalizes a copy and prepares the bottom level internally. It returns level 6
at scale 2^44. All N/2 complex slots must have magnitude at most 1. This is a caller
precondition; encrypted values cannot be inspected during evaluation. Setup
precomputes transform plaintexts and nonzero-noise evaluation keys; bootstrap
evaluation receives neither a secret key nor an RNG. These are educational
parameters with no production security claim.

Client code keeps rescaling explicit, with no bootstrap-specific bit counting:

```rust,ignore
let squared = initial.mul(&initial, &relin).rescale()?;
let refreshed = plan.bootstrap(&squared, &keys)?;
```

Preparation preserves represented values; it cannot recover precision lost
before bootstrap. Inputs already at level 1 and scale 2^44 remain supported.

Bootstrap normalization is fused into the forward transform: the half factor
for real/imaginary splitting and the EvalMod interval scaling consume no
separate levels. Exact multiplication by +/-i also consumes no level.
Bootstrap therefore spends12 primes instead of14 in the same18-prime chain.
The recurrence's91-bit weighted product safely fits two near-2^61 primes,
so its application-specific threshold is two levels. This permits17 initial
rounds and five between refreshes, reducing the32-round demo from eight
refreshes to three without increasing modulus size or weakening accuracy.

Signed rotations now satisfy `output[j] = input[(j + offset) mod slots]`.
Negative offsets invert positive rotations. Complex conjugation has a
separate `ConjugationKey` and `Ciphertext::conjugate` operation.

The codebase currently centers on one concrete backend:

- `RnsPoly` / `RnsBasis` in `src/rings/backends/rns_ntt/`
- arithmetic over `Z_Q[X] / (X^N + 1)`
- an RNS modulus chain made from NTT-friendly primes
- polynomials that rest in NTT form; coefficient form only where an
  operation needs it (encode/decode, key switching, rescale)

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

## Benchmarks

Criterion benches, generic over the ring backend (`benches/common.rs`):

| Bench | Covers |
|---|---|
| `arith` | modular mul/add/sub/pow on 30/50/61-bit primes, batched products |
| `ring` | NTT, add/mul/neg, automorphism, decompose, rescale, mod switch/raise, coefficient conversion; per N, level and prime width |
| `scheme` | encode/decode, encrypt/decrypt, mul + relin, rescale, rotate, conjugate, key generation, bootstrap (plan, keys, refresh) |
| `primes` | primality tests and prime generation |

```bash
cargo bench                          # quick survey, ~130 benches, about 5 minutes
BENCH_FULL=1 cargo bench             # N up to 2^14, 20-prime chains, bootstrap N=1024, budget sweep
cargo bench --bench ring -- ntt      # one group (Criterion filter)
BENCH_FULL=1 cargo bench --bench ring -- 'ring/mul/.*/N4096'  # precise numbers for a subset
```

The default run uses short timing (10 samples, 1 s each, no plots); numbers
are stable to about 1%. `BENCH_FULL=1` restores Criterion's defaults (100
samples, 5 s, HTML plots) and adds the large parameters, so use it with a
filter when you need precise comparisons.

Bench ids read `<group>/<op>/<backend>/N<degree>/<bits>b/L<level>`. To compare
a change, save a baseline first and compare against it:

```bash
cargo bench -- --save-baseline before
# ... change the code ...
cargo bench -- --baseline before
```

Each bench checks its results once before timing (decryption error, NTT round
trip, bootstrap within 1e-3), so a fast but wrong path fails instead of
reporting numbers. With `BENCH_FULL=1`, HTML reports land in
`target/criterion/report/index.html`.
Bootstrap above N=1024 is too slow for Criterion sampling; use
`cargo run --release --example bootstrap_large` instead.

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
src/ckks/                CKKS engine, ops, evaluator, ciphertext types
src/ckks/algorithms/     Slot sums, Chebyshev approximation
src/ckks/encoding/       CKKS encoder and special FFT
src/ckks/keys/           Secret, public, key switching (relin, rotation) keys
src/math/                Prime generation, sampling, utility math
src/rings/backends/      RNS-NTT polynomial backend
src/rings/traits.rs      Backend trait used by the engine
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
use toy_heaan_ckks::ckks::encoding::ckks_encoder::CkksEncoder;
use toy_heaan_ckks::ckks::engine::{CkksEngine, CkksParams};
use toy_heaan_ckks::math::generate_primes;
use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

const N: usize = 16;
const SCALE_BITS: u32 = 30;

let primes = generate_primes(31, 3, N as u64);
let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));

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
let ct = engine.encrypt(&pt, &pk, &mut rng);
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
