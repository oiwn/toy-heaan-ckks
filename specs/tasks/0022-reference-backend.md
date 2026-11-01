---
id: 0022-reference-backend
status: draft
scope: [src/rings/**, src/rings.rs, src/ckks/**, tests/**, Cargo.toml, README.md]
created: 2026-10-17
source: roadmap M2.6; operator request 2026-10-17
---
# Task: reference backend

Second `PolyRing` backend on big integers (roadmap M2.6): one modulus Q per level, textbook math, a correctness oracle for `rns_ntt` at toy N. `ckks/` needs no changes beyond making tests generic.

## Plan

- [ ] `rings/backends/reference.rs` + `reference/`: level chain built from the same primes as `RnsBasis` (same Q per level), shared `Arc` per level like `next()`.
- [ ] `PolyRing` impl: add/sub/neg mod Q, schoolbook negacyclic mul, automorphism by reindexing, `rescale` = `round(c / q_last)` mod `Q_(l-1)`, `mod_switch_to` = reduce mod `Q_l`, `from_rounded_f64`/`to_centered_f64` with the bound check, sampling mod Q.
- [ ] Gadget: base-2^d digits of the centered coefficient, `gadget_len = ceil(logQ / d)`, `gadget_term(k) = t * 2^(d*k)`; unit test of the gadget invariant.
- [ ] Backend-generic scheme tests run on both backends; cross-check: same inputs decrypt within tolerance, deterministic ops (mul, rescale, automorphism on equal coefficients) agree exactly.
- [ ] Docs: `specs/overview.md` backend section, README, roadmap M2.6 done.

## Open questions

- Big-integer crate: `crypto-bigint` (roadmap, fixed/boxed width, constant-time) or `num-bigint` (simplest); both are dependencies today, the other one goes.
- Generic test scope: all scheme tests, or a core set (encrypt/decrypt, add, mul + relin + rescale, rotate, `sum_slots`, Chebyshev).

## Acceptance

- [ ] `cargo fmt --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features`
- [ ] `cargo build --examples`
- [ ] `cargo bench --no-run`
- [ ] `cargo run --release --example sensor_std_dev`
- [ ] `typos`
- [ ] `gitleaks dir .`

## Log

- 2026-10-17 created: reference backend
- 2026-10-17 set source roadmap M2.6; operator request 2026-10-17
- 2026-10-17 scope add src/rings/** — bigint reference backend and backend-generic scheme tests
- 2026-10-17 scope add src/rings.rs — bigint reference backend and backend-generic scheme tests
- 2026-10-17 scope add src/ckks/** — bigint reference backend and backend-generic scheme tests
- 2026-10-17 scope add tests/** — bigint reference backend and backend-generic scheme tests
- 2026-10-17 scope add Cargo.toml — bigint reference backend and backend-generic scheme tests
- 2026-10-17 scope add README.md — bigint reference backend and backend-generic scheme tests
