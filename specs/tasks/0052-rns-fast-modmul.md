---
id: 0052-rns-fast-modmul
status: draft
scope: [src/rings/backends/rns_ntt/*.rs, src/math/primes.rs, src/math/modular.rs, benches/arith.rs]
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

- [ ] Record the baseline: Criterion baseline `pre-0052` (saved by 0018; rerun `cargo bench -- --save-baseline pre-0052` on a new machine), `cargo test --all-features` times per binary, release `bootstrap_demo` time and errors, `bootstrap_large --degree 1024` time.
- [ ] `NttTable` stores Shoup companions for twiddles, twist/untwist factors and `n_inv`; forward and inverse NTT butterflies use Shoup multiplication (lazy reduction only if it stays simple).
- [ ] Per-prime Barrett constant for general products (pointwise multiply, rescale, gadget terms, CRT); `mod_pow` and basis helpers use it. Shoup and Barrett live in the public `math::modular` (from 0018), next to the u128 reference, shared by `rns_ntt` and `math/primes.rs`; `benches/arith.rs` gains both variants.
- [ ] Property tests against the u128 reference: random operands plus 0 and q-1 for 30-, 59-, 60- and 61-bit primes, Shoup and Barrett paths.
- [ ] Bit-identical results: demo errors match the baseline exactly; compare with `cargo bench -- --baseline pre-0052` and record the new run on top of `specs/benchmarks.md` (summary in Findings).

## Acceptance

- [ ] `cargo fmt --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features`
- [ ] `cargo build --examples`
- [ ] `cargo bench --no-run`
- [ ] `cargo run --release --example bootstrap_demo`
- [ ] `typos`
- [ ] `specdev check`
- [ ] `git diff --check`

## Log

- 2026-11-03 created: rns fast modmul
- 2026-11-03 set source Operator: tests too slow, 2026-11-03
- 2026-11-03 scope add src/rings/backends/rns_ntt/*.rs — Shoup/Barrett modular multiplication replaces u128 division
- 2026-11-03 scope add src/math/primes.rs — Shoup/Barrett modular multiplication replaces u128 division
- 2026-11-03 set depends 0018-op-benchmarks
- 2026-11-03 scope add src/math/modular.rs — Shoup/Barrett in math::modular, compared in the arith bench
- 2026-11-03 scope add benches/arith.rs — Shoup/Barrett in math::modular, compared in the arith bench
