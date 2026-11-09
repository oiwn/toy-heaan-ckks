---
id: 0055-encrypted-linear-regression
status: draft
scope: [examples/linear_regression.rs, README.md]
created: 2026-11-09
source: Operator request 2026-11-09: ML example (linear regression) with bootstrap, checked against a plain reference
---
# Task: encrypted linear regression

## Context

An example of a real multi-step computation that needs bootstrapping and is easy to check: train linear regression by gradient descent on encrypted data. Decided 2026-11-09:

- data and model are both encrypted;
- 2 features plus a bias;
- 8 samples at N=16 (default), 32 with `--degree 64`.

Follows the `educational-examples` skill. No library changes: demo code stays in `examples/`.

**Math.**
- Model: `y ~ b + w1*x1 + w2*x2`. The bias is a weight on a constant feature x0 = 1, so the weights are `[b, w1, w2]` and the features `[1, x1, x2]`. One loop handles all three, and the shared scale arithmetic keeps the three weight ciphertexts at identical scales, so they can be added.
- Step: `r = sum_j w_j*x_j - y`, then `w_j <- w_j - lr/n * sum_i r_i*x_ij`.
- Slots: sample i in slot i of each feature and label ciphertext; each weight has its value in every slot (`sum_slots` produces that).
- Data:
  - fixed seed, `x1, x2` in [-1, 1];
  - `y = 0.2 + 0.7*x1 - 0.4*x2 + noise`, so the weights stay in (-1, 1) as bootstrapping requires;
  - `lr` about 0.5, tuned so plain gradient descent gets within about 1e-2 of least squares in 20 iterations.

**Two levels per iteration**, using the scale pattern of `RecurrencePlan::step` (`src/ckks/algorithms/recurrence.rs`):
- bring every term to the product's scale with `mul_const(1, target/own)`;
- rescale once;
- lift back with an exact power of two `mul_const(1, 2^k)`, which costs no level.

The two levels:
- **Residual:** one level.
- **Update:** `sum_slots(r * xs_j) + w_j` lifted, one level. `xs_j = x_j * lr/n` is precomputed once by the server.

Level budget:
- A product at about 2^88 needs to start at level 3 or higher.
- Weights start as Enc(0) at level 18, so the first 8 iterations need no refresh.
- After that, each refresh returns them at level 7, room for 3 iterations: 7 -> 5 -> 3 -> 1.
- 20 iterations = 8 + 4*3, so 4 refreshes of 3 weights each: 12 bootstraps.

**Keys** (client setup, done once; server code never touches `sk`):
- relinearization keys at levels 2..18;
- `generate_rotation_key_set` at levels 1..17;
- `BootstrapPlan::generate_keys` with `reference(N)` on 18 primes.

**Fallbacks.**
- **Scale mismatch on an add:** align with `mul_const_to_scale` (+1 level per iteration, 2 iterations per refresh), and record it in Findings.
- **N=64 bootstrap error** (about 2.6e-4 per refresh): measure it, set the bound once, and never loosen it.
- **Possible extension:** pack `b + i*w1` into one complex ciphertext to halve the bootstraps.

Rules: touch only scoped files. Never weaken tolerances or disable noise. If the same check fails twice, record it and ask.

## Plan

- [ ] Dataset with a fixed seed. Plain f64 gradient descent with the same steps, and least squares by normal equations (3x3 solve). Tune `lr` and the iteration count.
- [ ] Client setup: encrypt features, labels and Enc(0) weights; per-level relinearization and rotation keys; bootstrap plan and keys.
- [ ] Encrypted iteration: residual (one level) and update (one level) with lock-step scales; refresh all three weights when the level drops below 3.
- [ ] Narrated output per the skill:
  - intro: training with only + and *, the bias as a feature, why refresh;
  - parameters with reasons, and a dataset table;
  - a per-iteration trace: level, encrypted against plain weights, difference, `R` on refresh rows;
  - result: encrypted against plain gradient descent (`check`, 1e-3 at N=16; N=64 bound measured once), against least squares (distance, both MSEs), and predictions for 3 new points;
  - takeaways.
- [ ] `--degree 16|64`. README Application row: "Train a 2-feature model on encrypted data, with bootstrapping".

## Acceptance

- [ ] `cargo fmt --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features`
- [ ] `cargo build --examples`
- [ ] `cargo run --release --example linear_regression` (4 refreshes in the trace)
- [ ] `cargo run --release --example linear_regression -- --degree 64`
- [ ] output is ASCII; plain with `NO_COLOR=1`
- [ ] `typos`
- [ ] `specdev check`
- [ ] `git diff --check`

## Log

- 2026-11-09 created: encrypted linear regression
- 2026-11-09 set source Operator request 2026-11-09: ML example (linear regression) with bootstrap, checked against a plain reference
- 2026-11-09 scope add examples/linear_regression.rs — new example and its README row
- 2026-11-09 scope add README.md — new example and its README row
