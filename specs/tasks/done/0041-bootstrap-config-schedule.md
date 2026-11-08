---
id: 0041-bootstrap-config-schedule
status: done
scope: [src/ckks/algorithms/bootstrap.rs, tests/bootstrap.rs, examples/bootstrap_demo.rs, examples/recurrence_bootstrap_demo.rs, src/ckks/algorithms/recurrence.rs, src/ckks/algorithms/bootstrap/schedule.rs]
created: 2026-11-02
source: Operator-approved bootstrap split plan 2026-11-02 (from 0037)
depends: [0040-encrypted-polynomial-executor]
---
# Task: bootstrap config schedule

## Context

Replaces the fixed toy profile with explicit configuration and a derived schedule. Transforms stay one level per direction in this task; multi-group transforms are 0042. The current numerical settings become a regression case.

Rules: touch only scoped files (extend with `specdev task scope` first). If the same check fails twice, stop, record the failure and a hypothesis in Findings and ask the operator. Never weaken tolerances, disable noise or add `#[ignore]`.

## Plan

Operator decisions 2026-11-02:
- One `message_scale` (input target and inverse-transform assumption) plus `working_scale`. OpenFHE and Lattigo don't have a separate output scale either.
- The slot magnitude bound stays a documented precondition until 0043.
- The profile whitelist stays here; 0043 removes it.

- [x] `bootstrap/schedule.rs`: `BootstrapConfig` (message/working `Scale`, polynomial degree, double angles, CtS/StC level budgets 1/1, max secret weight) replaces `BootstrapParameters`. `BootstrapConfig::reference(ring_degree)` gives the current settings. Radius and interval bound are derived from the secret weight.
- [x] Read the bottom modulus from the backend (`mod_switch_to(1)` + `rescale_divisor`) instead of taking it as a parameter.
- [x] `BootstrapSchedule` with named stages (mod raise, CtS, conjugate split, EvalMod polynomial, double angles, StC): input/output level, exact simulated scale and divisors per stage. The EvalMod `ChebyshevPlan` is compiled once here.
- [x] Replace `INVERSE_TRANSFORM_DEPTH` and the constants 44/61/31/4/16.0 with config/schedule values. Expose `depth()`, `output_level()`, `output_scale(prepared)` and `required_keys()`.
- [x] `generate_keys` follows `required_keys()`; update the examples and the recurrence callers (constructor only).
- [x] Tests: the reference schedule has depth 11 (<= 12) and contiguous stages; keys match the requirements; invalid configs are rejected; output level/scale are asserted from the schedule; 8 seeds x 2 cycles within 1e-3.

## Findings

- Only the keys the schedule needs are generated: relinearization at the polynomial and double-angle levels, not at every level in `3..=18`. That changes the RNG stream, so errors move within noise. Worst demo error is now 2.63e-4 (0040: 1.57e-4; gate 1e-3). Setup is faster.
- The EvalMod private check now feeds the polynomial at its scheduled level and scale (17, about 2^61), not level 16, because the plan is compiled once.

## Review

- Scope: the 6 approved entries; `schedule.rs` is new (about 300 lines).
- The schedule uses the same `Scale` ops as evaluation; the executor's per-node check confirms planned and actual EvalMod scales agree in every bootstrap call.
- API: `BootstrapParameters::toy(q0)` is now `BootstrapConfig::reference(N)`. `BootstrapPlan` gains `config()`, `schedule()` and `bottom_modulus()`.
- The whitelist is unchanged in substance; the residue check on a supplied modulus is gone because nothing is supplied.

## Summary

`BootstrapConfig` and a chain-derived `BootstrapSchedule` replace the toy parameters and literal offsets/constants; keys follow the schedule's requirements.

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

- 2026-11-02 created: bootstrap config schedule
- 2026-11-02 set source Operator-approved bootstrap split plan 2026-11-02 (from 0037)
- 2026-11-02 set depends 0040-encrypted-polynomial-executor
- 2026-11-02 scope add src/ckks/algorithms/bootstrap.rs — Explicit BootstrapConfig and BootstrapSchedule replace toy profile constants
- 2026-11-02 scope add tests/bootstrap.rs — Explicit BootstrapConfig and BootstrapSchedule replace toy profile constants
- 2026-11-02 scope add examples/bootstrap_demo.rs — Explicit BootstrapConfig and BootstrapSchedule replace toy profile constants
- 2026-11-02 scope add examples/recurrence_bootstrap_demo.rs — Explicit BootstrapConfig and BootstrapSchedule replace toy profile constants
- 2026-11-02 scope add src/ckks/algorithms/recurrence.rs — Explicit BootstrapConfig and BootstrapSchedule replace toy profile constants
- 2026-11-02 scope add src/ckks/algorithms/bootstrap/schedule.rs — New config and schedule module
- 2026-11-02 draft approved: approved by operator
- 2026-11-02 advance draft → ready
- 2026-11-02 scope approved: src/ckks/algorithms/bootstrap.rs, tests/bootstrap.rs, examples/bootstrap_demo.rs, examples/recurrence_bootstrap_demo.rs, src/ckks/algorithms/recurrence.rs, src/ckks/algorithms/bootstrap/schedule.rs
- 2026-11-02 advance ready → in-progress/implement
- 2026-11-02 baseline: examples/encode_decode.rs@71d16406d270, examples/encrypt_add.rs@f2b8f286d018, examples/encrypt_mul.rs@95ea46780614, examples/horner_chain.rs@287e92e3251b, examples/rotation_demo.rs@ae35a052aac3, examples/rotation_stress.rs@0f1a817be753, examples/sensor_std_dev/main.rs@a1e4fae743d4, examples/sensor_std_dev/std_dev.rs@acd54fd1f4ac, examples/std_dev_8.rs@1dcb840d3609, src/ckks.rs@2615ecc20a4b, src/ckks/algorithms/bootstrap/executor.rs@04985b2bfb99, src/ckks/algorithms/bootstrap/polynomial.rs@c822cc7baf28, src/ckks/algorithms/bootstrap/transforms.rs@022348ef8279, src/ckks/algorithms/chebyshev.rs@76f9da012353, src/ckks/algorithms/sum_slots.rs@b4954675e46c, src/ckks/builder.rs@ee3fa125c3ae, src/ckks/encoding/ckks_encoder.rs@8f3e4ed1f6dd, src/ckks/engine.rs@829dde234c75, src/ckks/errors.rs@9e945a3b1aaa, src/ckks/evaluator.rs@83a706cb9989, src/ckks/ops.rs@ee512ce46da1, src/ckks/scale.rs@aff870709ee3, src/ckks/types.rs@601b4b28533e, tests/integration_mul.rs@ce478646e741
- 2026-11-02 advance in-progress/implement → in-progress/verify
- 2026-11-02 advance in-progress/verify → in-progress/review
- 2026-11-02 advance in-progress/review → approval
- 2026-11-02 approved: approved by operator
- 2026-11-02 advance approval → done
