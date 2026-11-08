---
id: 0018-op-benchmarks
status: done
scope: [benches/**, Cargo.toml, README.md, src/math.rs, src/math/*.rs, src/rings/backends/rns_ntt/*.rs]
created: 2026-10-16
source: operator request 2026-10-16; benchmark plan 2026-11-03
---
# Task: op benchmarks

## Context

The benchmark foundation (operator, 2026-11-03): measure the crate's parameters, compare backends and compare arithmetic variants. It is also the baseline for 0052 (Shoup/Barrett modmul, depends on this task), SIMD (0019) and wasm (0015). 0025 researches complete workloads on top of this harness; 0022 registers the reference backend in it.

Decisions (operator, 2026-11-03): Criterion only, with saved baselines (`--save-baseline` / `--baseline`) for before/after comparisons; Valgrind-based tools do not run on Apple Silicon, and CI instruction counts are 0025's question. Arithmetic variants live in a public `math::modular` module, like `is_prime_reference` today; this task moves the existing u128 `mul_mod` there as the baseline.

Limits: bootstrap at large N stays in `examples/bootstrap_large.rs` (one N=8192 run takes 40 s, too slow for Criterion sampling). `DIGIT_BITS` is a const, so digit width becomes a bench parameter only with hybrid key switching (0048).

## Plan

- [x] Entry points, no behaviour change: public `math::modular` (u128 `mul_mod`, `add_mod`, `sub_mod`, `pow_mod`; `primes.rs` uses it) and `NttTable::forward`/`inverse` wrappers over the private twisted kernels (operator, 2026-11-04).
- [x] Shared harness `benches/common.rs` (`autobenches = false`, explicit `[[bench]]` entries): `BENCH_FULL=1` enables the full matrix, the default stays small (operator, 2026-11-04); seeded setup (basis from `generate_primes`, keys, ciphertexts) generic over `P: PolyRing<N>`; one macro instantiates a group per degree and backend, so a new backend is one line. Bench IDs encode backend, N, prime bits and level.
- [x] `benches/arith.rs`: modular mul, reduce and pow over 30/50/61-bit primes, and a batch of N products (pointwise multiply). Only the u128 baseline now; 0052 adds its variants to the same group.
- [x] `benches/ring.rs`: forward/inverse NTT, add, mul, automorphism, decompose, gadget term, rescale, mod switch and mod raise at N = 2^10..2^12 (full: up to 2^14), levels 3/12 (full: 3/20), prime width 50/61.
- [x] `benches/scheme.rs`: encode/decode, encrypt/decrypt, add, mul + relin + rescale, rotate, conjugate, relin/rotation keygen at N = 2^10..2^11 (full: up to 2^13); `BootstrapPlan::new`, key generation and `bootstrap` at N=64/256 (full: 1024, plus budgets sweep), `sample_size(10)`. The ratio changes accuracy, not cost, so it is not a timing parameter. Each bench checks correctness once outside timing.
- [x] README bench section (commands, baselines, `BENCH_FULL` and filters); record one baseline run (machine, profile, key numbers) in Review.

## Findings

- `cargo bench -- --save-baseline` failed: the lib target's libtest harness rejects Criterion flags. Fixed with `[lib] bench = false` (standard Criterion setup).
- Criterion's default timing made the default run 24 min (132 benches x 8 s). `common::config()` uses 0.3 s warm-up, 1 s measurement, 10 samples and no plots by default (5.6 min, key numbers within ~1% of the long run); `BENCH_FULL` restores Criterion's defaults.
- `to_coeffs` returns i64, so conversions are benched on a message-sized polynomial (|c| < 2^40); a uniform 183-bit coefficient panics by design.
- Key-switched results (mul, rotate, conjugate) use a 1e-2 sanity bound: BV noise reaches 1.8e-4 at N=2048 with scale 2^50 (see `specs/research/large-degree-bootstrap.md`). The checks catch wrong paths, not precision.
- Bench config for bootstrap follows `examples/bootstrap_large.rs` (weight min(N/2, 64), ratio 8 above N=64) until 0045 makes the reference degree-dependent.
- First hot spots: `decompose` 87 ms at N=4096 L12 (48 digit polynomials, each converted to NTT form), `automorphism` ~15x a multiply (coefficient-form round trip), mod mul 5.6 ns on 61-bit primes against 1.5 ns on 30-bit (u128 division).

## Review

- Scope: `src/math/modular.rs` (new) and `primes.rs` use; `NttTable::forward`/`inverse` delegating to the existing kernels (now `pub(super)`); `Cargo.toml` bench targets; `benches/{common,arith,ring,scheme}.rs`; README. No library behaviour change: `bootstrap_demo` reproduces N=16 1.235318712267856e-5 / 1.5196209730106734e-5 and N=64 1.2215913830506398e-4 / 2.6291749341777456e-4 exactly; 158 lib tests (two new: add/sub wrap, public NTT round trip and negacyclic product).
- Full matrix setups and checks pass (`BENCH_FULL=1 ... -- no-such-bench`: ring 1 s, scheme 17 s including the N=1024 bootstrap and budget sweep). Full timing runs were not executed.
- Baseline `pre-0052` saved (Apple M4 Pro, rustc 1.99.1, default config); all results and observations in `specs/benchmarks.md`.
- `benches/primes.rs` keeps Criterion's defaults (4 benches, unchanged).

## Summary

Criterion benchmark foundation: backend-generic harness with arithmetic, ring, scheme and bootstrap benches, a quick default run plus `BENCH_FULL`, and a saved u128 baseline for 0052.

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features`
- [x] `cargo bench --no-run`
- [x] `cargo bench --bench arith -- --quick`
- [x] `typos`
- [x] `specdev check`
- [x] `git diff --check`

## Log

- 2026-10-16 created: op benchmarks
- 2026-10-16 set source operator request 2026-10-16
- 2026-10-16 scope add benches/** — per-op criterion benches
- 2026-10-16 scope add Cargo.toml — per-op criterion benches
- 2026-10-16 scope add README.md — per-op criterion benches
- 2026-11-03 set source operator request 2026-10-16; benchmark plan 2026-11-03
- 2026-11-03 scope add src/math/modular.rs — Public math::modular with the existing u128 mul_mod as the bench baseline
- 2026-11-03 scope add src/math.rs — Public math::modular with the existing u128 mul_mod as the bench baseline
- 2026-11-03 scope add src/math/primes.rs — Public math::modular with the existing u128 mul_mod as the bench baseline
- 2026-11-04 scope rm src/math/modular.rs — Covered by src/math/*.rs
- 2026-11-04 scope rm src/math/primes.rs — Covered by src/math/*.rs
- 2026-11-04 scope add src/math/*.rs — math::modular baseline and public NTT entry points
- 2026-11-04 scope add src/rings/backends/rns_ntt/*.rs — math::modular baseline and public NTT entry points
- 2026-11-04 draft approved: approved by operator
- 2026-11-04 advance draft → ready
- 2026-11-04 scope approved: benches/**, Cargo.toml, README.md, src/math.rs, src/math/*.rs, src/rings/backends/rns_ntt/*.rs
- 2026-11-04 advance ready → in-progress/implement
- 2026-11-04 advance in-progress/implement → in-progress/verify
- 2026-11-04 advance in-progress/verify → in-progress/review
- 2026-11-04 advance in-progress/review → approval
- 2026-11-04 approved: approved by operator
- 2026-11-04 advance approval → done
