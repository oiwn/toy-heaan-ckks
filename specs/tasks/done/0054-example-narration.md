---
id: 0054-example-narration
status: done
scope: [examples/*.rs, examples/sensor_std_dev/**, README.md, examples/support/narrate.rs, .agents/skills/educational-examples/SKILL.md]
created: 2026-11-09
source: "Operator decision 2026-11-09 (issue #39): rework all examples to the educational-examples skill"
depends: [0053-scheme-step-examples]
---
# Task: example narration

## Context

Issue #39 and the operator's 2026-11-09 decision: every example follows the `educational-examples` skill (`.agents/skills/educational-examples/SKILL.md`, from 0053). Each gets a narrated intro (why), parameters with reasons (`param`), worked steps (`formula = value`), a colored check against its plain model with a nonzero exit on a miss, and takeaways; colors only through `examples/support/narrate.rs`. Replace Unicode (box drawing, arrows, Greek) with ASCII. Math, gates and reported values stay the same.

Not here: the overlap between `std_dev_8` and `sensor_std_plain`, and moving `recurrence.rs` out of the library. Both are decisions for 0037.

Rules: touch only scoped files. Never weaken tolerances or disable noise.

## Plan

- [x] Basics: `encode_decode`, `keys`, `encrypt_add`, `encrypt_mul`.
- [x] Slots and depth: `rotation_demo`, `rotation_stress`, `horner_chain`.
- [x] Bootstrap: `bootstrap_demo`, `bootstrap_large`, `recurrence_bootstrap_demo`; errors unchanged.
- [x] Applications: `sensor_std_dev`, `sensor_std_plain`, `std_dev_8`, `sqrt_chebyshev_plain`.
- [x] README Examples table: update "Shows" where an example's content changed.

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features`
- [x] `cargo build --examples`
- [x] every debug example runs (loop over `examples/*.rs` minus the bootstrap ones and `modmul_methods`)
- [x] `cargo run --release --example bootstrap_demo`
- [x] `cargo run --release --example recurrence_bootstrap_demo` (`--degree 16/32/64/128`)
- [x] `cargo run --release --example bootstrap_large -- --degree 256` (the default N=8192 needs ~7.5 GB)
- [x] `typos`
- [x] `specdev check`
- [x] `git diff --check`

## Findings

- All 18 examples (plus `sensor_std_dev` unit tests) run and pass; output is ASCII (checked over every run's output). Reported values unchanged: `bootstrap_demo` N=16 1.235e-5 / 1.520e-5, N=64 1.222e-4 / 2.629e-4; recurrence N=16 8.592e-6, N=32 3.732e-5, N=64 3.289e-5, N=128 1.339e-6; `sensor_std_dev` variance 5.09e-4 degC^2, std_dev 1.08e-2 degC.
- New `narrate.rs` helpers: `confirm` (yes/no check), `table` (ASCII table; replaces the library's Unicode `toy_heaan_ckks::table` in examples), `elapsed` (ASCII durations; `{:?}` prints a micro sign). No example uses `toy_heaan_ckks::table` now: candidate for removal with `comfy-table` in 0037.
- Bug fixed: `sqrt_chebyshev_plain` fitted sqrt on [-1, 1]; every value at x < 0 was NaN, so the sweep printed 0 errors. It now fits [0.01, 0.2] (the std-dev demos' variances) with a 1e-3 gate (degree 7: 5.1e-4).
- Gates that were print-only now exit nonzero: `keys`, `encode_decode`, `std_dev_8` (new 1e-6 bound, measured 1.8e-10), `sensor_std_plain` (0.1 degC), `bootstrap_large` (1e-3 per cycle), `rotation_stress` (1e-3). Bounds otherwise unchanged; `encrypt_add` keeps its 3.2-as-std_dev bound (the known `error_variance` quirk).
- `sensor_std_plain` printed the affine slope as `4 / width` (635); the library map `t = (2v - (hi + lo)) / (hi - lo)` has slope `2 / width` = 317.5, now printed. The `~635` in `std_dev.rs`'s test comment may be stale too (left for 0037).
- `horner_chain` takes ~30 s in debug (N=8192); its doc and the README now say `--release`.

## Review

Diff matches the plan and scope (plus `narrate.rs` and the skill, added for the new helpers). Every example now has the skill's shape: intro, parameters with reasons, narrated steps, colored checks with a nonzero exit, takeaways; colors only through `narrate.rs`. Computations are the same code paths with the same seeds; the only behavior changes are the `sqrt_chebyshev_plain` interval fix and print-only gates becoming exits, both listed in Findings. `recurrence_bootstrap_demo` now prints one level trace per run instead of 32 round lines, keeping output near one screen.

## Summary

Every example reworked to the `educational-examples` skill (narrated, colored, self-checking, ASCII only), with new `narrate.rs` helpers (`confirm`, `table`, `elapsed`) and a fix for `sqrt_chebyshev_plain`'s NaN interval.

## Log

- 2026-11-09 created: example narration
- 2026-11-09 set source Operator decision 2026-11-09 (issue #39): rework all examples to the educational-examples skill
- 2026-11-09 set depends 0053-scheme-step-examples
- 2026-11-09 scope add examples/*.rs — every existing example moves to the narrated colored style; README Shows text
- 2026-11-09 scope add examples/sensor_std_dev/** — every existing example moves to the narrated colored style; README Shows text
- 2026-11-09 scope add README.md — every existing example moves to the narrated colored style; README Shows text
- 2026-11-09 draft approved: approved by operator
- 2026-11-09 advance draft → ready
- 2026-11-09 scope approved: examples/*.rs, examples/sensor_std_dev/**, README.md
- 2026-11-09 advance ready → in-progress/implement
- 2026-11-09 baseline: .agents/.gitignore@4c41cd39e6d2, .agents/skills/.gitignore@6737dd17c5ee, .agents/skills/educational-examples/SKILL.md@9e5db3a6e99e, examples/support/narrate.rs@fab9fb541ef4, src/ckks/algorithms/bootstrap.rs@a5e92b067b64
- 2026-11-09 scope add examples/support/narrate.rs — new helpers (confirm, table, elapsed) found while reworking; skill lists them
- 2026-11-09 scope add .agents/skills/educational-examples/SKILL.md — new helpers (confirm, table, elapsed) found while reworking; skill lists them
- 2026-11-09 advance in-progress/implement → in-progress/verify
- 2026-11-09 advance in-progress/verify → in-progress/review
- 2026-11-09 advance in-progress/review → approval
- 2026-11-09 approved: approved by operator
- 2026-11-09 advance approval → done
