---
id: 0045-bootstrap-ratio-default
status: done
scope: [src/ckks/algorithms/bootstrap.rs, src/ckks/algorithms/bootstrap/*.rs, tests/bootstrap.rs, examples/bootstrap_large.rs, benches/scheme.rs, README.md]
created: 2026-11-03
source: Operator agreed large-degree roadmap, 2026-11-03
depends: [0043-bootstrap-general-final]
---
# Task: bootstrap ratio default

## Context

Roadmap step 1 of `specs/research/large-degree-bootstrap.md` (operator agreed 2026-11-03). The reference ratio 17 fails above N=128; the ratio sweep shows the error is U-shaped in the ratio (noise amplification above the optimum, cubic sine error below), with the optimum moving down with N: about 10 at N<=64, 6 at N=1024, 4 at N=8192. No fixed ratio passes everywhere, so `BootstrapConfig::reference(N)` picks it from N, like OpenFHE's fitted correction factor but from our own data. The fit is redone after hybrid key switching (0048).

A smaller ratio means a larger message scale (2^(61 - ratio)): the output level must keep headroom for one product after bootstrap, and recurrence (0044) relies on the current 2^44.

Decided 2026-11-08 (operator): keep ratio 17 (2^44) at N<=64, where recurrence and the small-N tests rely on it (worst 2.6e-4); above, ratio = 17 - log2 N clamped at 4 (10 at N=128, 4 at N=8192), a line in log2 N fitted to the sweep. Revised after the sweep: 17 at N=128 left only 1.4x margin (7.0e-4). `reference(N)` also caps the secret weight at 64 and uses five double angles above weight 32 (as the sweep and `benches/scheme.rs`). Stop and ask if no ratio passes at some N.

Rules: touch only scoped files (extend with `specdev task scope` first). If the same check fails twice, stop, record the failure and a hypothesis in Findings and ask the operator. Never weaken tolerances, disable noise or add `#[ignore]`.

## Plan

- [x] Fill the sweep at N=16/32/128/256/512/2048/4096 (`bootstrap_large` gets N=16/32; `bootstrap_demo` for N=16) with weight N/2 up to 64 and record it in the research doc.
- [x] `BootstrapConfig::reference(N)` derives `message_ratio_bits` from N with a documented rule fitted to the sweep (17 at N<=64, else 17 - log2 N, at least 4), weight min(N/2, 64) and five double angles above weight 32; `bootstrap_large` and `benches/scheme.rs` use the reference by default.
- [x] Headroom: the schedule checks that the output level can square the output once; error names the level.
- [x] Regression numbers in tests/demos follow the new defaults; the 1e-3 gate and noise stay; record old and new errors in Findings.

## Findings

- Sweep table and rule: `specs/research/large-degree-bootstrap.md` (Measurements, "Full sweep"). Optimum falls one bit per doubling of N; `17 - log2 N` sits at or next to it from N=128 up.
- Defaults, `bootstrap_large` with no overrides (cycle0 / cycle1): N=256 ratio 9 7.1e-6 / 1.4e-5; N=1024 ratio 7 1.6e-5 / 2.7e-5; N=4096 ratio 5 4.6e-5 / 9.9e-5. Before (ratio 17, same weight): 1.7e-3 / 4.1e-3, 1.7e-2 / 3.2e-2, 1.5e-1 / 3.3e-1, all failing. N=8192 ratio 4 (earlier sweep): 2.0e-4 / 3.2e-4.
- N=128 reference: ratio 17 3.7e-4 / 7.0e-4 (weight 64, five angles) -> ratio 10 2.7e-6 / 7.3e-6.
- N<=64 unchanged by construction (ratio 17, weight N/2 <= 32, four angles): `bootstrap_demo` reproduces N=16 1.235318712267856e-5 / 1.5196209730106734e-5 and N=64 1.2215913830506398e-4 / 2.6291749341777456e-4 exactly; `tests/bootstrap.rs` and recurrence untouched; recurrence demo N=16 8.592e-6, N=64 3.289e-5 (32 rounds, 3 refreshes).
- Headroom: the reference's shortest chain grows from 12 to 13 primes (the output needs a spare level for one squaring). Existing 18-prime chains and `bootstrap_large`/bench chains (depth + 2) are unaffected.
- `benches/scheme.rs` now uses the reference for ratio, weight and angles; N=256 bench moves from ratio 8 to 9 and still passes its check.

## Review

Diff matches the plan and scope. `reference_ratio_bits` is a pure function of N with named constants and a doc comment citing the sweep, OpenFHE's analogue and the 0048 refit; it is re-exported and unit-tested at 16..2^14. The weight cap and angle rule mirror what the bench helper did, which is now gone. The headroom check reuses `Builder::headroom` and names the level; the shortest-chain test asserts the message. No tolerances changed, no noise disabled, no tests ignored. Only N>=128 behaviour changes, where the old default failed or had 1.4x margin. Residual: the reference keeps one-level transforms, impractical above N~1024 (N/2-1 rotation keys); large-N users pass grouped budgets as `bootstrap_large` does. Out of scope here.

## Summary

`BootstrapConfig::reference(N)` picks the message ratio from N (17 up to N=64, then `17 - log2 N`, at least 4), caps the secret weight at 64 with five double angles above weight 32, and the schedule checks the output can be squared once; bootstrap now passes 1e-3 at N=16 to 8192 with the defaults.

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

## Log

- 2026-11-03 created: bootstrap ratio default
- 2026-11-03 set source Operator agreed large-degree roadmap, 2026-11-03
- 2026-11-03 set depends 0043-bootstrap-general-final
- 2026-11-03 scope add src/ckks/algorithms/bootstrap.rs — Degree-dependent ratio default, output headroom check, regression numbers
- 2026-11-03 scope add src/ckks/algorithms/bootstrap/*.rs — Degree-dependent ratio default, output headroom check, regression numbers
- 2026-11-03 scope add tests/bootstrap.rs — Degree-dependent ratio default, output headroom check, regression numbers
- 2026-11-03 scope add examples/bootstrap_large.rs — Degree-dependent ratio default, output headroom check, regression numbers
- 2026-11-03 scope add examples/bootstrap_demo.rs — Degree-dependent ratio default, output headroom check, regression numbers
- 2026-11-08 draft approved: approved by operator
- 2026-11-08 scope add benches/scheme.rs — Drop the bench-local ratio/weight override; the reference now holds it
- 2026-11-08 scope add README.md — Describe the degree-dependent bootstrap default
- 2026-11-08 advance draft → ready
- 2026-11-08 scope approved: src/ckks/algorithms/bootstrap.rs, src/ckks/algorithms/bootstrap/*.rs, tests/bootstrap.rs, examples/bootstrap_large.rs, examples/bootstrap_demo.rs, benches/scheme.rs, README.md
- 2026-11-08 advance ready → in-progress/implement
- 2026-11-08 baseline: Cargo.toml@0c00a118600a, benches/arith.rs@de3344c65714, examples/modmul_methods.rs@687a1da29bf5, src/math/modular.rs@3219ba9fd0cc, src/rings/backends/rns_ntt/basis.rs@87e619629038, src/rings/backends/rns_ntt/errors.rs@e654b9d5f7f0, src/rings/backends/rns_ntt/poly.rs@38851611e556
- 2026-11-08 scope rm examples/bootstrap_demo.rs — Not edited: N<=64 defaults unchanged, demo runs as acceptance only
- 2026-11-08 advance in-progress/implement → in-progress/verify
- 2026-11-08 advance in-progress/verify → in-progress/review
- 2026-11-08 advance in-progress/review → approval
- 2026-11-08 approved: approved by operator
- 2026-11-08 advance approval → done
