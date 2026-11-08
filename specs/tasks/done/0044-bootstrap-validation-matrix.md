---
id: 0044-bootstrap-validation-matrix
status: done
scope: [src/ckks/algorithms/recurrence.rs, examples/**/*.rs, tests/bootstrap.rs, src/ckks/algorithms/bootstrap.rs, README.md]
created: 2026-11-02
source: Operator split of 0043, 2026-11-02
depends: [0043-bootstrap-general-final, 0045-bootstrap-ratio-default]
---
# Task: bootstrap validation matrix

## Context

Split from 0043, which made bootstrap accept general inputs. The callers below still assume the old profile: N=16/64, primes near 2^61, scale 2^44. The matrix uses `reference(N)` from 0045. `bootstrap_large` and `benches/scheme.rs` already use it. After acceptance, 0037 resumes.

Decided 2026-11-08: `RecurrencePlan::new(context, scale)` takes the scale and derives its refresh level from the chain; `recurrence_bootstrap_demo` runs N=16/32/64/128 with chain lengths and expected rounds derived from the plans. README's stale bootstrap rows are updated too.

Rules: touch only scoped files. If the same check fails twice, stop, record the failure and a hypothesis in Findings, and ask. Never weaken tolerances, disable noise or add `#[ignore]`.

## Plan

- [x] `recurrence.rs`: drop the N=16/64 and near-2^61 checks (lines 41-53) and the hard-coded 2^44 (lines 121-125, 162). Use actual scales and bootstrap's derived message scale; keep the formula.
- [x] Update `bootstrap_demo`, `recurrence_bootstrap_demo` and `tests/bootstrap.rs`. Document the breaking `message_ratio_bits` API in overview and the Summary.
- [x] Matrix at N=16/32/64/128, varying scale, weight, polynomial, budgets, chain length and prime width. Plus 8 seeds x 2 cycles at N=16/64, all within 1e-3.
- [x] Mark `specs/research/generic-bootstrap.md` complete.

## Manual checks

- [x] Operator reviews the full general-bootstrap integration. — user: reviewed and approved by operator

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features`
- [x] `cargo build --examples`
- [x] `cargo run --release --example bootstrap_demo`
- [x] `cargo run --release --example recurrence_bootstrap_demo`
- [x] `typos`
- [x] `specdev check`
- [x] `git diff --check`

## Findings

- Recurrence refresh rule: the lowest level L >= 2 whose modulus holds `2*(scale_bits+1)+3` bits plus one, and whose level L-1 holds the lifted result. Every prime from L up must leave the rescaled result no more than half a bit above `scale`. On 61-bit primes at 2^44 to 2^54 this gives L=2, as before.
- Unchanged results: `bootstrap_demo` (N=16 1.235e-5/1.520e-5, N=64 1.222e-4/2.629e-4) and the recurrence demo at N=16/64 (8.592e-6 / 3.289e-5) match the recorded numbers exactly.
- New: recurrence demo N=32 3.73e-5 (18 primes), N=128 1.34e-6 (19 primes, scale 2^51), each 32 rounds with three refreshes. Matrix errors are in `specs/research/generic-bootstrap.md`. All are within 1e-3, at least 17x below it. The 61-bit-bottom/60-bit mixed chain passes, so the drift guard does not reject it.
- Debug-profile `tests/bootstrap.rs` takes 15 s, including N=128.

## Review

The diff matches the plan and scope. Recurrence keeps its formula and inputs are unmodified. The derived checks name the failing level. Tests assert the derived round and refresh counts, and anchor them to the old literals (17 rounds, 32 rounds with 3 refreshes) for the reference chain. The 1e-3 gate is unchanged; no `#[ignore]` and no noise changes. One gap was filled: a secret above the configured weight is now tested at key generation.

## Summary

Recurrence and the bootstrap demos/tests run at any degree and scale the chain supports. Breaking change: `RecurrencePlan::new(context, scale)`, with a refresh level derived from the chain. `recurrence_bootstrap_demo` runs N=16/32/64/128, the final validation matrix passes, and `generic-bootstrap.md` is complete.

## Log

- 2026-11-02 created: bootstrap validation matrix
- 2026-11-02 set source Operator split of 0043, 2026-11-02
- 2026-11-02 set depends 0043-bootstrap-general-final
- 2026-11-02 scope add src/ckks/algorithms/recurrence.rs — Recurrence/examples migration and the final validation matrix
- 2026-11-02 scope add examples/**/*.rs — Recurrence/examples migration and the final validation matrix
- 2026-11-02 scope add tests/bootstrap.rs — Recurrence/examples migration and the final validation matrix
- 2026-11-02 scope add src/ckks/algorithms/bootstrap.rs — Recurrence/examples migration and the final validation matrix
- 2026-11-03 set depends 0043-bootstrap-general-final, 0045-bootstrap-ratio-default
- 2026-11-08 draft approved: approved by operator
- 2026-11-08 scope add README.md — Stale bootstrap status rows in What Works and Known Limitations
- 2026-11-08 advance draft → ready
- 2026-11-08 scope approved: src/ckks/algorithms/recurrence.rs, examples/**/*.rs, tests/bootstrap.rs, src/ckks/algorithms/bootstrap.rs, README.md
- 2026-11-08 advance ready → in-progress/implement
- 2026-11-08 advance in-progress/implement → in-progress/verify
- 2026-11-08 advance in-progress/verify → in-progress/review
- 2026-11-08 advance in-progress/review → approval
- 2026-11-08 approved: approved by operator
- 2026-11-08 advance approval → done
