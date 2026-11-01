---
id: 0031-bootstrap-integration
status: done
scope: [src/ckks/**, examples/**, tests/**, Cargo.toml, README.md]
created: 2026-11-01
depends: [0027-bootstrap-mod-raise, 0028-bootstrap-dense-transforms, 0030-bootstrap-eval-mod]
---
# Task: bootstrap integration

## Plan

- [x] Compose parameter validation, centered ModRaise, complex transforms, EvalMod and output normalization behind a fallible public bootstrap API; expose setup separately from evaluation.
- [x] Validate the 18-prime candidate ledger at N=16, then N=64; freeze a working profile only after actual-prime scale and measured error checks pass.
- [x] Add a runnable full-packing complex bootstrap example and fixed-seed regression tests: maximum error <=1e-3 for slot magnitude<=1, versus original input and pre-refresh decrypt.
- [x] Demonstrate output at least level4/logp44 (or a reviewed equivalent), another multiplication/rescale with explicit scale restoration, and two refresh cycles; test missing keys, invalid contexts and insufficient depth.
- [x] Document supported parameters, input preconditions, output depth/error, rotation behavior change and educational security limits; update overview and README.

## Context

Design: `specs/research/bootstrap.md`. Depends on ModRaise, dense transforms and EvalMod, with public evaluation-key support transitively required. No decrypt/re-encrypt or secret-key access in evaluation. Implemented; awaiting final operator review. Production-size tuning, performance benchmarking and formal verification are separate work.

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features`
- [x] `cargo build --examples`
- [x] `specdev check`
- [x] `typos`

## Review

Validated with the full all-features suite, release bootstrap regressions (both degrees, eight seeds, six patterns, two refreshes), strict Clippy and example builds. Public bootstrap evaluation has no secret key or RNG. Supported profiles are fixed and validated at construction; general-purpose parameter combinations are not exposed.

## Summary

Implemented bootstrap integration.

## Log

- 2026-11-01 created: bootstrap integration
- 2026-11-01 scope add src/ckks/algorithms — bootstrap composition
- 2026-11-01 scope add src/ckks/algorithms.rs — bootstrap composition
- 2026-11-01 scope add src/ckks/builder.rs — bootstrap composition
- 2026-11-01 scope add src/ckks/encoding — bootstrap composition
- 2026-11-01 scope add src/ckks/encoding.rs — bootstrap composition
- 2026-11-01 scope add src/ckks/engine.rs — bootstrap composition
- 2026-11-01 scope add src/ckks/errors.rs — bootstrap composition
- 2026-11-01 scope add src/ckks/evaluator.rs — bootstrap composition
- 2026-11-01 scope add src/ckks/keys — bootstrap composition
- 2026-11-01 scope add src/ckks/keys.rs — bootstrap composition
- 2026-11-01 scope add src/ckks/ops.rs — bootstrap composition
- 2026-11-01 scope add src/ckks/types.rs — bootstrap composition
- 2026-11-01 scope add examples/encode_decode.rs — bootstrap composition
- 2026-11-01 scope add examples/encrypt_add.rs — bootstrap composition
- 2026-11-01 scope add examples/encrypt_mul.rs — bootstrap composition
- 2026-11-01 scope add examples/horner_chain.rs — bootstrap composition
- 2026-11-01 scope add examples/keys.rs — bootstrap composition
- 2026-11-01 scope add examples/rotation_demo.rs — bootstrap composition
- 2026-11-01 scope add examples/rotation_stress.rs — bootstrap composition
- 2026-11-01 scope add examples/sensor_std_dev — bootstrap composition
- 2026-11-01 scope add examples/sensor_std_plain.rs — bootstrap composition
- 2026-11-01 scope add examples/sqrt_chebyshev_plain.rs — bootstrap composition
- 2026-11-01 scope add examples/std_dev_8.rs — bootstrap composition
- 2026-11-01 scope add tests/integration_mul.rs — bootstrap composition
- 2026-11-01 set depends 0027-bootstrap-mod-raise, 0028-bootstrap-dense-transforms, 0030-bootstrap-eval-mod
- 2026-11-01 scope rm src/ckks/algorithms — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm src/ckks/algorithms.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm src/ckks/builder.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm src/ckks/encoding — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm src/ckks/encoding.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm src/ckks/engine.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm src/ckks/errors.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm src/ckks/evaluator.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm src/ckks/keys — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm src/ckks/keys.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm src/ckks/ops.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm src/ckks/types.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm examples/encode_decode.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm examples/encrypt_add.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm examples/encrypt_mul.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm examples/horner_chain.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm examples/keys.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm examples/rotation_demo.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm examples/rotation_stress.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm examples/sensor_std_dev — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm examples/sensor_std_plain.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm examples/sqrt_chebyshev_plain.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm examples/std_dev_8.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope rm tests/integration_mul.rs — replace shell-expanded paths with intended module scope
- 2026-11-01 scope add src/ckks/** — declare complete module scope including new files
- 2026-11-01 scope add examples/** — declare complete module scope including new files
- 2026-11-01 scope add tests/** — declare complete module scope including new files
- 2026-11-01 scope add Cargo.toml — declare complete module scope including new files
- 2026-11-01 scope add README.md — document bootstrap example and contract
- 2026-11-01 draft approved: Implement the plan.
- 2026-11-01 advance draft → ready
- 2026-11-01 scope approved: src/ckks/**, examples/**, tests/**, Cargo.toml, README.md
- 2026-11-01 advance ready → in-progress/implement
- 2026-11-01 baseline: src/rings/backends/rns_ntt/poly.rs@5c302b70b49f, src/rings/errors.rs@fea138cd1c12, src/rings/traits.rs@8a54adc1cbc6
- 2026-11-01 advance in-progress/implement → in-progress/verify
- 2026-11-01 advance in-progress/verify → in-progress/review
- 2026-11-01 advance in-progress/review → approval
- 2026-11-01 approved: ok, looks good.
- 2026-11-01 advance approval → done
