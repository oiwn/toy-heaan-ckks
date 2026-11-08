# Benchmarks

Measured results of the Criterion benches (`benches/`, task 0018). Each run records its environment and the change it measures; newer runs go on top. Bench ids read `<group>/<op>/<backend>/N<degree>/<bits>b/L<level>`. 

## How to run

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
cargo bench -- --save-baseline before    # then, after a change:
cargo bench -- --baseline before         # --baseline-lenient also measures new benches
```

The default run uses short timing (10 samples, 1 s each, no plots), stable to about 1%. `BENCH_FULL=1` restores Criterion's defaults (100 samples, 5 s, HTML plots in `target/criterion/report/index.html`) and adds the large parameters; use it with a filter. Each bench checks its result once before timing (decryption error, NTT round trip, bootstrap within 1e-3). Bootstrap above N=1024 is too slow for Criterion; use `cargo run --release --example bootstrap_large`.

## 2026-11-06 — Shoup/Barrett (`post-0052`)

Task 0052: NTT butterflies, twist/untwist and fixed factors (rescale, gadget terms) use Shoup multiplication; pointwise products use Barrett. Plain Shoup, values always fully reduced; results bit-identical to the baseline.

- Machine, toolchain and config as the baseline below.
- Command: `cargo bench -- --baseline-lenient pre-0052`, then `--save-baseline post-0052`. Before = `pre-0052` median, after = this run's median.

### Modular arithmetic (`arith`)

| Op | 30-bit | 50-bit | 61-bit |
|---|---|---|---|
| mul u128 | 1.48 ns | 5.52 ns | 5.52 ns |
| mul barrett | 6.06 ns | 6.13 ns | 6.12 ns |
| mul shoup | 0.89 ns | 0.87 ns | 0.87 ns |
| pow u128 / barrett | 145 / 125 ns | 739 / 214 ns | 967 / 264 ns |
| pointwise mul x4096, u128 / barrett | 6.08 / 4.12 us | 28.7 / 4.15 us | 21.3 / 4.14 us |
| pointwise mul by a constant x4096, u128 / shoup | 6.09 / 1.83 us | 26.3 / 1.83 us | 21.1 / 1.82 us |

The scalar `mul` bench measures latency: one Barrett product is a chain of four dependent multiplies (6 ns), while the u128 remainder on a 30-bit modulus takes the routine's fast path. In a loop the products overlap, and Barrett costs 1.0 ns against 5.2 ns (61-bit) — the case the library uses.

### Before -> after (61-bit primes)

| Bench | Before | After | Speedup |
|---|---|---|---|
| NTT forward N1024 / N4096 | 33.0 / 156 us | 6.10 / 31.9 us | 5.4x / 4.9x |
| NTT inverse N1024 / N4096 | 38.4 / 177 us | 6.03 / 31.9 us | 6.4x / 5.5x |
| ring mul N1024 L3 / N4096 L12 | 16.7 / 267 us | 3.45 / 57.3 us | 4.8x / 4.7x |
| automorphism N4096 L12 | 4.02 ms | 788 us | 5.1x |
| decompose N1024 L12 / N4096 L12 | 18.2 / 87.1 ms | 4.16 / 20.7 ms | 4.4x / 4.2x |
| rescale N4096 L12 | 2.20 ms | 452 us | 4.9x |
| mod_raise N4096 L12 | 2.19 ms | 557 us | 3.9x |
| from_coeffs N4096 L12 | 2.01 ms | 471 us | 4.3x |
| to_coeffs N4096 L3 / L12 | 1.64 / 11.6 ms | 1.11 / 9.50 ms | 1.5x / 1.2x |
| gadget_term N4096 L12 | 8.43 us | 3.37 us | 2.5x |
| encode / decode N2048 | 9.41 / 11.8 ms | 8.84 / 11.1 ms | 1.06x |
| encrypt / decrypt N2048 | 2.09 ms / 93.6 us | 538 / 23.8 us | 3.9x |
| mul + relin N1024 / N2048 | 11.4 / 24.4 ms | 2.73 / 5.71 ms | 4.2x / 4.3x |
| rotate N2048 | 26.5 ms | 6.09 ms | 4.4x |
| relin keygen N2048 | 48.0 ms | 15.1 ms | 3.2x |
| bootstrap N64: plan / keygen / bootstrap | 694 us / 55.0 / 60.7 ms | 429 us / 21.9 / 15.9 ms | 1.6x / 2.5x / 3.8x |
| bootstrap N256: plan / keygen / bootstrap | 9.95 / 532 / 503 ms | 7.35 / 186 / 126 ms | 1.4x / 2.9x / 4.0x |
| generate_primes 61-bit x5 | 164 us | 76.4 us | 2.1x |

50-bit primes now cost the same as 61-bit (mul 3.49 vs 3.45 us, rescale 22.1 vs 22.2 us at N1024 L3); the gap in the baseline came from the u128 remainder. Add, neg and mod_switch are unchanged (no products).

Outside Criterion (release):

| Run | Before | After |
|---|---|---|
| `bootstrap_large --degree 1024`: keygen / one bootstrap / whole run | 3.35 / 2.82 / 9.75 s | 1.11 / 0.68 / 3.14 s |
| `bootstrap_demo` whole run | 1.08 s | 0.65 s |
| `cargo test --all-features`: lib / `tests/bootstrap.rs` | 4.75 / 36.2 s | 2.39 / 15.2 s |

### Observations

- Everything built on the NTT is 4-5x faster; the remaining costs are elsewhere. `to_coeffs` (BigUint CRT) and encode/decode (special FFT in f64, BigUint) barely moved and now dominate their operations.
- `decompose` (20.7 ms at N4096 L12) is still about 360x a multiply: 4L digit polynomials, each one forward NTT per channel, so it grows with L^2. Key switching follows it (rotate about 300x a ring multiply at N1024 L8). Hybrid key switching (0048) cuts the digit count.
- `automorphism` is still a coefficient-form round trip (~14x a multiply); an NTT-domain permutation would remove both transforms.
- Lazy Harvey butterflies would save the conditional subtractions in the NTT; worth it only once the NTT dominates again.

## 2026-11-04 — u128 baseline (`pre-0052`)

Baseline before Shoup/Barrett (0052): every modular product is `(a * b as u128) % q`.

- Machine: Apple M4 Pro, macOS 26 (Darwin 25.6.0); rustc 1.99.1; `bench` profile (release + debuginfo).
- Command: `cargo bench -- --save-baseline pre-0052` (default config: 0.3 s warm-up, 1 s measurement, 10 samples; 132 benches in 5.6 min). A run with Criterion's defaults (24 min) gave the same numbers within about 1%.
- Values are Criterion's median estimate.
- Library state: worktree of 0043/0018 (BV key switching, 20-bit digits).

### Modular arithmetic (`arith`, `modmul`)

| Op | 30-bit | 50-bit | 61-bit |
|---|---|---|---|
| mul | 1.55 ns | 5.58 ns | 5.60 ns |
| add | 0.79 ns | 0.62 ns | 0.62 ns |
| sub | 0.84 ns | 0.63 ns | 0.63 ns |
| pow (exp q-2) | 155 ns | 760 ns | 996 ns |
| pointwise mul x4096 | 6.22 us | 29.1 us | 21.6 us |

### NTT, one 61-bit channel (`ring`, `ntt`)

| N | forward | inverse |
|---|---|---|
| 1024 | 33.0 us | 38.4 us |
| 2048 | 72.0 us | 82.9 us |
| 4096 | 156 us | 177 us |

### Ring operations, 61-bit primes (`ring`)

| Op | N1024 L3 | N1024 L12 | N2048 L3 | N2048 L12 | N4096 L3 | N4096 L12 |
|---|---|---|---|---|---|---|
| add | 0.48 us | 2.05 us | 0.97 us | 4.07 us | 2.04 us | 8.20 us |
| mul | 16.7 us | 66.8 us | 33.2 us | 133 us | 67.2 us | 267 us |
| neg | 0.89 us | 3.50 us | 1.70 us | 6.75 us | 3.31 us | 13.2 us |
| automorphism | 216 us | 862 us | 468 us | 1.87 ms | 1.01 ms | 4.02 ms |
| decompose (20-bit) | 1.22 ms | 18.2 ms | 2.69 ms | 39.8 ms | 5.86 ms | 87.1 ms |
| gadget_term | 2.05 us | 2.99 us | 4.06 us | 5.88 us | 7.27 us | 8.43 us |
| rescale | 118 us | 476 us | 255 us | 1.02 ms | 546 us | 2.20 ms |
| mod_switch (to L1) | 0.11 us | 0.10 us | 0.19 us | 0.19 us | 0.32 us | 0.32 us |
| mod_raise (L1 to L) | 157 us | 472 us | 339 us | 1.02 ms | 725 us | 2.19 ms |
| from_coeffs | 107 us | 432 us | 234 us | 937 us | 505 us | 2.01 ms |
| to_coeffs | 388 us | 2.90 ms | 802 us | 6.00 ms | 1.64 ms | 11.6 ms |

### Prime width, level 3 (`ring`, `prime_width`)

| Op | N1024 50-bit | N1024 61-bit | N4096 50-bit | N4096 61-bit |
|---|---|---|---|---|
| mul | 22.0 us | 16.6 us | 91.0 us | 66.7 us |
| rescale | 144 us | 118 us | 687 us | 551 us |

### Scheme operations, 8-prime 61-bit chain, scale 2^50 (`scheme`)

| Op | N1024 | N2048 |
|---|---|---|
| encode | 2.46 ms | 9.41 ms |
| decode | 3.76 ms | 11.8 ms |
| encrypt | 964 us | 2.09 ms |
| decrypt | 46.9 us | 93.6 us |
| add | 4.89 us | 10.2 us |
| mul + relin | 11.4 ms | 24.4 ms |
| rescale | 636 us | 1.36 ms |
| rotate | 12.4 ms | 26.5 ms |
| conjugate | 12.3 ms | 26.5 ms |
| relin keygen | 22.4 ms | 48.0 ms |
| rotation keygen | 22.9 ms | 49.0 ms |

### Bootstrap (`scheme`, `bootstrap`)

Config as `examples/bootstrap_large.rs`: weight min(N/2, 64), ratio 8 above N=64, chain = depth + 2 primes of 61 bits.

| N | budgets cts/stc | plan | keygen | bootstrap |
|---|---|---|---|---|
| 64 | 3/2 | 694 us | 55.0 ms | 60.7 ms |
| 256 | 4/3 | 9.95 ms | 532 ms | 503 ms |

Outside Criterion (`bootstrap_large`, release, 2026-11-03): N=1024 about 2.7 s per bootstrap; N=8192 plan 16 s, keys 50 s (~7.5 GB), 40 s per bootstrap.

### Primes (`primes`)

| Bench | Time |
|---|---|
| is_prime miller_rabin 2^31-1 | 1.83 us |
| is_prime reference 2^31-1 | 8.12 us |
| is_prime miller_rabin ~2^63 | 12.0 us |
| is_prime reference ~2^63 | 540 ms |
| generate_primes 30-bit x5, N=8192 | 14.3 us |
| generate_primes 61-bit x5, N=8192 | 164 us |

### Observations

- The u128 remainder dominates: a 61-bit mod mul (5.6 ns) costs 3.6x a 30-bit one (1.55 ns), where the hardware product is the same. 50-bit primes are slower than 61-bit for whole operations (mul +33%, rescale +22% at N=1024), so narrower primes buy nothing here. 0052 targets this.
- `decompose` is the largest ring cost: 87 ms at N=4096 L12, about 330x a multiply. It produces 4 x L digit polynomials and moves each to NTT form, so it grows with L^2 (L3 to L12: x15). Key switching inherits it: rotate and mul + relin cost about 12 ms at N=1024 L8, about 275x a ring multiply at that level.
- `automorphism` costs about 15x a multiply (coefficient-form round trip); an NTT-domain permutation would remove the transforms.
- `to_coeffs` costs 3-7x `from_coeffs` because of the CRT reconstruction (big integers above two channels).
- Ring operations scale linearly in N and in L, except `decompose` (quadratic in L) and `mod_switch` (truncation only).
- Bootstrap scales about 8x from N=64 to N=256 and about 5x from N=256 to N=1024; key generation costs about one bootstrap.
