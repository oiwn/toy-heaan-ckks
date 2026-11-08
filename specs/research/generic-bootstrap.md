# General CKKS bootstrapping implementation

## Status

Incomplete. Ciphertexts carry a `Scale` (0039). Rescale divides by the actual prime, and EvalMod and the Chebyshev evaluator run on the encrypted `ChebyshevPlan` executor (0040). Bootstrap is driven by `BootstrapConfig` and a derived `BootstrapSchedule` (0041) and evaluates multi-group transforms at scheduled levels (0042). Since 0043 the reference whitelist is gone: validation is derived from the config and the chain. Recurrence/examples migration and the final matrix are 0044.

The umbrella `specs/tasks/0037-cleanup-and-refactoring.md` is blocked until the chain 0038-0044 lands. Do not mark the general-bootstrap acceptance complete before 0044.

## Implemented components

- [x] `src/ckks/scale.rs`: validated positive finite f64 Scale/ScaleError, representable power-of-two construction including subnormals, checked multiplication/division, square root and symmetric compatibility tolerance 64*f64::EPSILON. Exported through ckks; CkksError propagates ScaleError. This is not a metadata migration.
- [x] Structural bootstrap cleanup: validation/encoding helpers, explicit type names, documented key indices, end-to-end tests under tests/bootstrap.rs and private checks at the bottom of bootstrap.rs.
- [x] `src/ckks/algorithms/bootstrap/transforms.rs`: sparse radix-two factors, inverse factors, configurable even grouping with extra stages in earlier groups, plaintext application and entry counts. Bit-reversed coefficient halves are shared by both directions; slotwise modular reduction commutes with that permutation.
- [x] `src/ckks/algorithms/bootstrap/linear.rs` (0042): `EncodedTransform` encodes each group at its scheduled level and plaintext scale, as (offset, plaintext) pairs for the nonzero diagonals. It evaluates with a rotation callback. Bootstrap takes CtS/StC budgets `1..=log2(N/2)`, rotation keys follow the groups per level, and `BootstrapKeys::check` names any missing key before evaluation starts.
  - CtS normalization (`1/(2*radius)`) sits in the last CtS group. In the first group, every later group added noise to a 2*radius smaller signal: errors were 20x worse at N=16 and 80x worse at N=64. StC restoration sits in the last group.
- [x] `src/ckks/algorithms/bootstrap/polynomial.rs`: target-scale operation graph, recursive Chebyshev decomposition with linear leaves, actual-divisor scale simulation, conservative noiseless node/coefficient headroom checks, reference evaluation and required relin-level inspection. It does not simulate encoding quantization or incoming noise.
- [x] `src/ckks/algorithms/bootstrap/executor.rs`: `ChebyshevPlan::execute` runs the same graph on ciphertexts and checks level and scale at every node. Bootstrap EvalMod and `eval_chebyshev_encrypted` use it.
  - The planner itself still accepts chains whose prime widths differ from the scale in the power tree (drift doubles with every squaring). Bootstrap guards it since 0043: the schedule rejects EvalMod and double-angle scales more than `MAX_SCALE_DRIFT_BITS` (1 bit) from the working scale. Re-targeting power scales (Lattigo-style) for mixed widths is a later option.
- [x] `src/ckks/algorithms/bootstrap/schedule.rs`: `BootstrapConfig`, plus a `BootstrapSchedule` derived from the chain. The schedule names each stage with its levels, exact simulated scale and divisors; it also holds the compiled EvalMod plan and the required keys, and `generate_keys` follows those. The literal offsets and constants are gone. The bottom modulus is read from the backend.

## Validation in executable bootstrap (0043)

`BootstrapConfig` takes `message_ratio_bits` (reference 17); the message scale is `2^(round(log2 q0) - ratio)`, so the reference stays 2^44. `BootstrapPlan::new` checks:
- Config: power-of-two ring degree >= 8, transform budgets `1..=log2(N/2)`, secret weight `1..=N`, polynomial degree >= 1.
- Chain (in `BootstrapSchedule::derive`): ratio >= 1 leaving a scale >= 2, enough levels (the reference needs 12 primes), transform and double-angle product headroom, and scale drift within 1 bit of the working scale.
- Not checked: approximation accuracy for a given weight/degree/double-angle choice (caller's choice), ternary secret weight (checked at key generation), slot magnitudes at most 1 (caller precondition).

Encrypted smoke (seed 0, asymmetric, cycle0/cycle1 max error, 18 primes, depth 11):
- N=32 reference: 1.7e-5 / 5.4e-5.
- N=16, uniform 60-bit chain, working 2^60, message 2^43: 2.7e-5 / 5.1e-5.
- N=16, degree 63, three double angles: 9.5e-6 / 1.8e-5.

Plans also build at N=128 (reference) and N=64 with weight 16 and ratio 20.

Rejected: alternating 60/61-bit primes and a 2^62 working scale on 61-bit primes (drift), 11-prime chains, ratio 0 or 61, degree 0, weight 0 or N+1, budgets 0 or log2(N/2)+1.

Large degrees (N up to 8192, `examples/bootstrap_large.rs`): the reference ratio 17 fails above N=128; data, hypothesis and fix candidates in `specs/research/large-degree-bootstrap.md`.

## Remaining implementation, in dependency order

Tracked as separate task files (operator decision 2026-11-02); each has its own scope, plan and acceptance.

- [x] `specs/tasks/0038-fallible-ciphertext-arithmetic.md` — remove ciphertext Add/Sub operator impls; every add/sub returns `CkksResult` with level/context/scale checks. Results unchanged.
- [x] `specs/tasks/0039-scale-metadata-migration.md` — `logp`/`scale_bits` -> `Scale` everywhere; rescale still divides by nominal `2^bits`, so results stay bit-identical.
- [x] `specs/tasks/done/0040-encrypted-polynomial-executor.md` — rescale divides by the actual prime; encrypted `ChebyshevPlan` executor; bootstrap EvalMod and the Chebyshev evaluator migrate to it.
- [x] `specs/tasks/done/0041-bootstrap-config-schedule.md` — `BootstrapConfig` + `BootstrapSchedule` replace `toy()` and fixed offsets/constants (one-level transforms).
- [x] `specs/tasks/done/0042-multi-group-transforms.md` — transform groups at scheduled levels, per-level derived keys.
- [x] `specs/tasks/done/0043-bootstrap-general-final.md` — message ratio config; derived validation replaces the profile whitelist.
- [ ] `specs/tasks/0044-bootstrap-validation-matrix.md` — migrate recurrence/examples, final validation matrix below.

Why two phases: with actual-divisor rescale, values reaching an add by different paths differ in scale by up to about 2^-41 (relative), above the 64*eps tolerance. Flipping rescale before the executor exists breaks bootstrap `term`/`eval_mod`, Clenshaw and examples at once.

## Numerical rules already selected

Keep Chebyshev coefficients in their original basis. Split at the largest power-of-two m not exceeding the degree: P=T_m*Q+R, Q[0]=c[m], Q[k]=2*c[m+k] for k>0, R[m-k]-=c[m+k]. Cache power-of-two terms using T_(2m)=2*T_m^2-1.

For divisor q, a linear leaf encodes its coefficient at target_scale*q/input_scale. A product branch requests its quotient at target_scale*q/power_scale and its remainder at the product's final scale/level. Constant leaves need no multiplication. Use one graph for simulation and execution. Reject encoding/product headroom failures, insufficient levels and unrepresentable scales with stage-specific errors.

Keep f64 precision/range, backend-generic PolyRing and compile-time degree. Retain full-complex continuous-cosine bootstrapping. Packing/ring switching, secret encapsulation and other approximation families are separate features. "Toy" means experimental status, not restricted functionality. Preserve noise conventions; resolve the known variance/std-dev inconsistency separately. Slot/noise bounds are caller preconditions, not encrypted facts construction can certify.

## Required final validation

- [ ] Scale integration: non-power-of-two encode/decode, actual-divisor and mixed-width rescaling, mismatches, invalid factors and modulus switching.
- [ ] Encrypted polynomial evaluation and multi-group transforms with varied budgets/divisors; verify actual execution follows the graph.
- [ ] Noisy bootstrap at N=16/64 with original eight seeds/two cycles, plus N=32/128, different scales, secret weights, polynomial settings, budgets, chain lengths and mixed-width chains. Retain the 1e-3 accuracy gate; do not weaken tolerances or disable noise.
- [ ] At most 12 levels consumed for legacy numerical settings with one-level transforms; assert derived output properties rather than a hardcoded six.
- [ ] Invalid configuration, wrong context, excessive secret weight, missing keys, insufficient headroom and input immutability.
- [ ] Full formatting, strict Clippy, tests, example builds, bootstrap/recurrence demos at both degrees, spelling, specdev and diff checks on the completed integration.

## Verification evidence and skill experiment

Baseline for 0038-0043 (2026-11-02, current worktree incl. headroom fix): fmt, strict all-target/all-feature Clippy, cargo test --all-features (141 library, 10 integration, 4 example tests, 1 doctest, 5 ignored doctests) and example builds pass. Release bootstrap_demo max_error, output level 6 / logp 44: N=16 cycle0 8.059113425635853e-6, cycle1 1.7666600415874825e-5; N=64 cycle0 7.9910586165085e-5, cycle1 2.0767430746117748e-4. 0038/0039 reproduced these exactly.

After 0040 (actual-divisor rescale, executor EvalMod), the release bootstrap_demo outputs level 7 with log2(scale) about 44. EvalMod saves one level. Max errors:
- N=16: cycle0 9.045642510505948e-6, cycle1 2.606107558812035e-5.
- N=64: cycle0 7.411339230750804e-5, cycle1 1.5688101781596232e-4.

Recurrence demo: 32/32 rounds with 3 refreshes, max_error 1.13e-5. tests/bootstrap.rs passes 8 seeds x 6 scenarios x 2 cycles at N=16/64 against 1e-3.

After 0041, keys follow the schedule, so the RNG stream differs from 0040. The release bootstrap_demo has depth 11, output level 7 and log2(scale) about 44. Max errors:
- N=16: cycle0 1.235318712267856e-5, cycle1 1.5196209730106734e-5.
- N=64: cycle0 1.2215913830506398e-4, cycle1 2.6291749341777456e-4.

Recurrence demo with refresh: 32/32 rounds with 3 refreshes, max_error 8.59e-6 at N=16 and 3.29e-5 at N=64. 0042 reproduces these numbers bit for bit with budgets 1/1.

Multi-group (0042, seed 0, second-cycle max error over the asymmetric and random_disk cases):

| Budgets (CtS,StC) | Depth | N=16 | N=64 |
|---|---|---|---|
| (1,1) | 11 | 2.4e-5 | 1.1e-4 |
| (2,2) | 13 | 2.5e-5 | 1.9e-4 |
| (3,3) | 15 | 5.8e-5 | 3.7e-4 |
| (1,3) | 13 | 2.2e-5 | 1.0e-4 |
| (3,1) | 13 | 3.7e-5 | 4.8e-4 |

Full radix-two budgets need 2*log2(N/2)+9 levels; at N=64 that is 19 > 17, so the schedule rejects it. Encrypted groups match the plaintext factored transform within 1e-6 at N=16/32/64/128 for budgets 1, 2 and log2(slots).

Compiler tests cover degrees 0/1/2/3/15/31/63 over 101 reference inputs, mixed 60/61-bit bases, node dependencies, multiply/add/rescale metadata and invalid coefficients/levels/headroom. These tests exclude encrypted noise and coefficient quantization.

Transform dense references pass at N=16/32/64/128 for every level budget. N=1024 plaintext round-trip, linear-per-stage storage and invalid-input checks pass. Release bootstrap demo passes both degrees/two cycles, largest observed error about 2.08e-4; N=64 recurrence demo passes. These demos use the restricted existing bootstrap, not the proposed general executor.

Formatting, strict all-target/all-feature Clippy, example builds, spelling, specdev and diff checks passed during component development. Historical Scale/transform-stage full runs had 136/139 library tests; these are superseded by the 141-test run, not additional acceptance.

Code-health: Scale has no hotspots (max cc 5/nesting 1/length 20); transforms have none (cc 8/nesting 3/length 41). The polynomial compiler's linear leaf was extracted after measurement; recursive compilation and graph interpretation remain candidates, not gates. No metric rules proposed. The original decide priority experiment remains UNSURE at confidence 0.091; its full log is in specs/research/file-review.md. No further decide judgments were used.

Reference analysis and sources: specs/research/bootstrap.md (OpenFHE and Lattigo).
