---
id: 0052-rns-fast-modmul
status: done
scope: [src/rings/backends/rns_ntt/*.rs, src/math/primes.rs, src/math/modular.rs, benches/arith.rs, README.md]
created: 2026-11-03
source: Operator: tests too slow, 2026-11-03
depends: [0018-op-benchmarks]
---
# Task: rns fast modmul

## Context

Tests were slow (operator, 2026-11-03): `tests/bootstrap.rs` took 154 s and the lib tests 30 s unoptimized. A `[profile.test] opt-level = 1` (direct change, same day) brought that to 37 s / 4.9 s. The rest is arithmetic: every modular product is `(a as u128 * b as u128) % q` (`rns_ntt/poly.rs` `mul_mod`, `basis.rs` CRT and `mod_pow`, `math/primes.rs`), and the u128 remainder is a software routine (`__umodti3` / `u128_div_rem`). macOS `sample` on release `complex_bootstrap_64` (19 s) puts about 87% of CPU time there.

OpenFHE and Lattigo avoid the division: Shoup multiplication for a fixed operand with a precomputed companion `floor(w * 2^64 / q)` (NTT twiddles, twist factors, scalar constants), and Barrett or Montgomery reduction for general products. Results here must stay bit-identical; this is a speed change only.

Rules: touch only scoped files (extend with `specdev task scope` first). If the same check fails twice, stop, record the failure and a hypothesis in Findings and ask the operator. Never weaken tolerances, disable noise or add `#[ignore]`.

## Plan

- [x] Record the baseline: Criterion baseline `pre-0052` (saved by 0018; rerun `cargo bench -- --save-baseline pre-0052` on a new machine), `cargo test --all-features` times per binary, release `bootstrap_demo` time and errors, `bootstrap_large --degree 1024` time.
- [x] `NttTable` stores Shoup companions for twiddles, twist factors and the untwist factors with `n_inv` folded in; forward and inverse NTT butterflies use plain Shoup multiplication, values always fully reduced (operator, 2026-11-06: lazy Harvey butterflies only as a follow-up if the numbers ask for it). `NttTable::new` rejects moduli >= 2^62 with a new `RnsNttError` (operator, 2026-11-06).
- [x] Per-prime Barrett constant for general products (pointwise multiply); Shoup for fixed factors (rescale, gadget terms); `pow_mod` and basis helpers use `modular`, duplicates in `rns_ntt` go. Shoup and Barrett live in the public `math::modular` (from 0018), next to the u128 reference, shared by `rns_ntt` and `math/primes.rs`; `benches/arith.rs` gains both variants.
- [x] Property tests against the u128 reference: random operands plus 0, 1 and q-1 for 30-, 50-, 59-, 60-, 61- and 62-bit primes (plus small and even moduli), Shoup, Barrett, wide Barrett inputs and `pow`.
- [x] Bit-identical results: demo errors match the baseline exactly; compare with `cargo bench -- --baseline pre-0052` and record the new run on top of `specs/benchmarks.md` (summary in Findings).

## Findings

- Before (2026-11-06, this machine): lib tests 4.75 s, `tests/bootstrap.rs` 36.2 s; `bootstrap_demo` 1.08 s, errors N=16 1.235318712267856e-5 / 1.5196209730106734e-5, N=64 1.2215913830506398e-4 / 2.6291749341777456e-4; `bootstrap_large --degree 1024` keygen 3.35 s, 2.82 s per bootstrap.
- After: lib tests 2.39 s, `tests/bootstrap.rs` 15.2 s; demo 0.65 s with the same errors bit for bit; N=1024 keygen 1.11 s, 0.68 s per bootstrap. NTT 5-6x, ring mul, rescale, automorphism and decompose 4-5x, key switching 4.3x, bootstrap 3.8-4.0x. Full tables: `specs/benchmarks.md` (2026-11-06).
- A single Barrett product is slower than u128 in the scalar latency bench (6 ns, four dependent multiplies); in loops, which is how the library uses it, it is 1.0 ns against 5.2 ns. Shoup is 0.9 ns scalar, 0.45 ns in a loop.
- `cargo bench -- --baseline` aborts on bench ids missing from the baseline (the new variants); `--baseline-lenient` compares the rest. README notes it.
- `primes.rs` needed no edit: Miller-Rabin's `pow_mod` takes the Barrett path below 2^62 and the u128 reference above (it tests 63-bit candidates); its squaring loop stays on `mul_mod`.
- Next hot spots: `to_coeffs` (BigUint CRT), encode/decode (f64 FFT), `decompose` (L^2 forward NTTs, 0048), automorphism round trip; lazy NTT only after those.

## Review

- Scope: `math/modular.rs` (Shoup, Barrett, `pow_mod_reference`, `MAX_FAST_MODULUS`, 4 property tests); `rns_ntt/errors.rs` (`ModulusTooLarge`); `basis.rs` (`ShoupFactors`, companions and Barrett per table, `n_inv` folded into the untwist, iterative table powers, private `mod_pow`/`sub_mod` removed, size check test); `poly.rs` (kernels on Shoup, `mul_assign` on Barrett, rescale/gadget term on Shoup, private helpers removed); `benches/arith.rs` (variants with a correctness check each); README (`--baseline-lenient`).
- Public pub fields of `NttTable` are unchanged; the reconstruct helper still reduces unchecked residues first, as the removed `sub_mod` did.
- No tolerance changed; 163 lib tests (5 new) and every integration test pass; results bit-identical.

## Summary

Shoup and Barrett modular multiplication replace the u128 remainder in the RNS-NTT backend: NTT 5-6x, ring and key-switching operations 4-5x, bootstrap 4x faster, results bit-identical.

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features`
- [x] `cargo build --examples`
- [x] `cargo bench --no-run`
- [x] `cargo bench -- --baseline-lenient pre-0052`
- [x] `cargo run --release --example bootstrap_demo`
- [x] `typos`
- [x] `specdev check`
- [x] `git diff --check`

## Log

- 2026-11-03 created: rns fast modmul
- 2026-11-03 set source Operator: tests too slow, 2026-11-03
- 2026-11-03 scope add src/rings/backends/rns_ntt/*.rs — Shoup/Barrett modular multiplication replaces u128 division
- 2026-11-03 scope add src/math/primes.rs — Shoup/Barrett modular multiplication replaces u128 division
- 2026-11-03 set depends 0018-op-benchmarks
- 2026-11-03 scope add src/math/modular.rs — Shoup/Barrett in math::modular, compared in the arith bench
- 2026-11-03 scope add benches/arith.rs — Shoup/Barrett in math::modular, compared in the arith bench
- 2026-11-06 draft approved: approved by operator
- 2026-11-06 advance draft → ready
- 2026-11-06 scope approved: src/rings/backends/rns_ntt/*.rs, src/math/primes.rs, src/math/modular.rs, benches/arith.rs
- 2026-11-06 advance ready → in-progress/implement
- 2026-11-06 scope add README.md — Bench section: --baseline-lenient for benches added after the baseline
- 2026-11-06 advance in-progress/implement → in-progress/verify
- 2026-11-06 advance in-progress/verify → in-progress/review
- 2026-11-06 advance in-progress/review → approval
- 2026-11-06 approved: approved by operator
- 2026-11-06 advance approval → done
