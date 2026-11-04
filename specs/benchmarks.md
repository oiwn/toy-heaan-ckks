# Benchmarks

Measured results of the Criterion benches (`benches/`, task 0018). Each run records its environment and the change it measures; newer runs go on top. Bench ids read `<group>/<op>/<backend>/N<degree>/<bits>b/L<level>`. How to run: README "Benchmarks".

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
