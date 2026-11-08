# Repository file review

## Coverage

Inventory: 121 tracked files at the start of 0037; no relevant untracked files then. New files introduced during this session are included below. `[x]` means an initial source review, not operator acceptance or proof of numerical correctness. Unchecked files remain unreviewed; tests or authorship alone do not complete a review.

Excluded from cleanup edits: generated/build output (`target/`), `tmp/`, `.git/`, installed agent/tool internals and caches. Tracked archived specs, completed task history and the reference PDF remain in the checklist for relevance/link review only; do not rewrite accepted history or third-party reference content. `_index.md` is CLI-owned. Ignored local files are outside the repository-owned inventory.

## Baseline and experiment

- `specdev check`: passed before changes; working tree initially clean.
- Baseline verification (2026-11-02): `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`, `cargo build --examples`, `typos`, `git diff --check` and `specdev check` passed. `cargo test --all-features` passed: 130 library tests, 8 integration tests, 4 example tests and 1 doctest; 5 doctests ignored. Examples were built, not individually run. These are baseline results, not final acceptance of future refactors.
- `specdev stats src/ --format json`: 36 files, 336 functions, no parse errors, 11 flagged functions. Includes inline tests. Complexity p90/p95/max: 4/5/34; nesting: 2/2/5; length: 29/43/266. Default thresholds remain unchanged.
- 2026-11-02 | first cleanup candidate | bootstrap_constructor (0.394, conf 0.091); recurrence_step 0.304; noise_mismatch 0.303 | UNSURE. Below confidence guideline; no useful model signal. Bootstrap first is based on inspection. No new metric rules proposed.

## First review findings

- **Bootstrap:** operator approved structural cleanup first, with end-to-end tests outside the source and small private tests inline. Implemented private validation/diagonal helpers, removed aliases, documented key indices, named end-to-end scenarios and extracted automorphism/stage checks. End-to-end cases retain eight seeds, original RNG call order, both degrees and two refresh cycles. Private stages run separately at seed zero, including asymmetric transform references, preparation at every level and 44/27-bit scales, idempotence and missing-key rejection. Public API and numerical formulas unchanged. Await operator review.
- **Bootstrap contract:** constructor still requires exactly 18 levels; corrected the diagnostic to match. Operator subsequently authorized general implementation. Scheduling is now active unfinished work, not deferred scope; the remaining whitelist and integration gaps are listed in `specs/research/generic-bootstrap.md`.
- 2026-11-02 | code-health bootstrap helper/test extraction | accepted via "Implement the plan."; implementation awaiting review. Source before: max cc34 (test), constructor cc22/115 lines; after: source max cc16, constructor cc6/62 lines. End-to-end runner reduced from combined 266 to 130 lines; remains a hotspot. No threshold changes proposed.
- **Bootstrap verification:** full `cargo test --all-features` passed (130 library, 10 integration, 4 example tests, 1 doctest; 5 ignored doctests). After extending the private preparation assertions, both private-stage tests passed again. Strict all-target/all-feature Clippy, formatting, example builds, spelling, specdev and diff checks passed. Release bootstrap demo passed both degrees/two cycles, largest observed error about 2.0e-4; recurrence demo passed at N=16 and N=64. No tolerance or profile changes. Review found the seeded end-to-end RNG calls retain their original ordering; moved private stages do not consume RNG in that runner.
- **Recurrence:** `step` complexity 17 is largely guards and error propagation. Keep production flow for now; the 127-line test scenario is a later readability candidate. Review checked level consumption, scale restoration, context checks and test scenarios by inspection; no correctness proof claimed.
- **Engine:** encryption uses `error_variance` directly as standard deviation, relin/rotation generation uses its square root, and public-key generation hardcodes 3.2. Resolve intended noise policy explicitly before code changes.
- **Ring trait:** negative rotations now use modular inverse offsets without conjugation; the old ctx rotation finding is superseded. Trait remains the scheme/backend boundary. Parameter/panic contracts deserve a later focused review.
- **Primes:** `is_prime` complexity 11 follows small-number rejection and Miller-Rabin witness loops; retain the algorithm structure. Existing tests cover small ranges against trial division, selected pseudoprimes and near-u64 inputs. Public `is_ntt_friendly_prime` panics lack a `# Panics` section; bare module rustdoc URL remains a documentation candidate. Large-input trial-division cost and congruence-helper boundary arithmetic deserve focused validation, not an assumption that passing selected ranges covers all inputs.

## Per-file checklist

- [ ] `.github/workflows/additional-checks.yml`
- [ ] `.github/workflows/code-quality.yml`
- [ ] `.github/workflows/platform-coverage.yml`
- [ ] `.github/workflows/publish.yml`
- [ ] `.github/workflows/test.yml`
- [ ] `.gitignore`
- [ ] `AGENTS.md`
- [ ] `CHANGELOG.md`
- [ ] `Cargo.toml`
- [ ] `Justfile`
- [ ] `LICENSE`
- [ ] `README.md`
- [ ] `_typos.toml`
- [ ] `benches/primes.rs`
- [ ] `examples/bootstrap_demo.rs`
- [ ] `examples/encode_decode.rs`
- [ ] `examples/encrypt_add.rs`
- [ ] `examples/encrypt_mul.rs`
- [ ] `examples/horner_chain.rs`
- [ ] `examples/keys.rs`
- [ ] `examples/recurrence_bootstrap_demo.rs`
- [ ] `examples/rotation_demo.rs`
- [ ] `examples/rotation_stress.rs`
- [ ] `examples/sensor_std_dev/contract.rs`
- [ ] `examples/sensor_std_dev/main.rs`
- [ ] `examples/sensor_std_dev/std_dev.rs`
- [ ] `examples/sensor_std_plain.rs`
- [ ] `examples/sqrt_chebyshev_plain.rs`
- [ ] `examples/std_dev_8.rs`
- [ ] `prek.toml`
- [ ] `rust-toolchain.toml`
- [ ] `rustfmt.toml`
- [ ] `specdev.toml`
- [ ] `specs/2025-1935.pdf`
- [ ] `specs/archive/ckks_add_mul_audit.md`
- [ ] `specs/archive/refactoring.md`
- [ ] `specs/archive/textbook_encode_decode.md`
- [ ] `specs/cleanup.md`
- [ ] `specs/ctx.md`
- [ ] `specs/demo.md`
- [ ] `specs/ideas.md`
- [ ] `specs/overview.md`
- [ ] `specs/parameters.md`
- [ ] `specs/research/bootstrap.md`
- [ ] `specs/research/bootstrap_experiment.py`
- [ ] `specs/roadmap.md`
- [ ] `specs/tasks/0013-sage-education-scripts.md`
- [ ] `specs/tasks/0015-wasm-feasibility.md`
- [ ] `specs/tasks/0017-debug-tooling.md`
- [ ] `specs/tasks/0018-op-benchmarks.md`
- [ ] `specs/tasks/0019-simd-research.md`
- [ ] `specs/tasks/0020-gpu-research.md`
- [ ] `specs/tasks/0021-npu-research.md`
- [ ] `specs/tasks/0022-reference-backend.md`
- [ ] `specs/tasks/0023-creusot-integration-research.md`
- [ ] `specs/tasks/0024-kani-evaluation.md`
- [ ] `specs/tasks/0025-ckks-backend-benchmark-research.md`
- [ ] `specs/tasks/0034-nonlinear-smoothing-demo.md`
- [ ] `specs/tasks/0035-cubic-activation-demo.md`
- [ ] `specs/tasks/0037-cleanup-and-refactoring.md`
- [ ] `specs/tasks/_index.md`
- [ ] `specs/tasks/done/0001-green-baseline.md`
- [ ] `specs/tasks/done/0002-range-contract-std-dev.md`
- [ ] `specs/tasks/done/0003-sensor-demo-example.md`
- [ ] `specs/tasks/done/0004-encoder-wide-staging.md`
- [ ] `specs/tasks/done/0005-remove-dead-code.md`
- [ ] `specs/tasks/done/0006-trim-dependencies.md`
- [ ] `specs/tasks/done/0007-module-layout.md`
- [ ] `specs/tasks/done/0008-ciphertext-levels.md`
- [ ] `specs/tasks/done/0009-ops-composition.md`
- [ ] `specs/tasks/done/0010-chebyshev-module.md`
- [ ] `specs/tasks/done/0011-backend-trait.md`
- [ ] `specs/tasks/done/0012-generic-key-switching.md`
- [ ] `specs/tasks/done/0014-remove-commented-code.md`
- [ ] `specs/tasks/done/0016-bootstrap-research.md`
- [ ] `specs/tasks/done/0026-bootstrap-automorphisms-keys.md`
- [ ] `specs/tasks/done/0027-bootstrap-mod-raise.md`
- [ ] `specs/tasks/done/0028-bootstrap-dense-transforms.md`
- [ ] `specs/tasks/done/0029-bootstrap-shallow-polynomials.md`
- [ ] `specs/tasks/done/0030-bootstrap-eval-mod.md`
- [ ] `specs/tasks/done/0031-bootstrap-integration.md`
- [ ] `specs/tasks/done/0032-bootstrap-input-preparation.md`
- [ ] `specs/tasks/done/0033-bootstrapped-recurrence-demo.md`
- [ ] `specs/tasks/done/0036-bootstrap-depth-optimization.md`
- [ ] `src/ckks.rs`
- [ ] `src/ckks/algorithms.rs`
- [x] `src/ckks/algorithms/bootstrap.rs`
- [ ] `src/ckks/algorithms/chebyshev.rs`
- [x] `src/ckks/algorithms/recurrence.rs`
- [ ] `src/ckks/algorithms/sum_slots.rs`
- [ ] `src/ckks/builder.rs`
- [ ] `src/ckks/encoding.rs`
- [ ] `src/ckks/encoding/ckks_encoder.rs`
- [ ] `src/ckks/encoding/special_fft.rs`
- [x] `src/ckks/engine.rs`
- [ ] `src/ckks/errors.rs`
- [ ] `src/ckks/evaluator.rs`
- [ ] `src/ckks/keys.rs`
- [ ] `src/ckks/keys/conjugation_key.rs`
- [ ] `src/ckks/keys/key_switching.rs`
- [ ] `src/ckks/keys/public_key.rs`
- [ ] `src/ckks/keys/relin_key.rs`
- [ ] `src/ckks/keys/rotation_key.rs`
- [ ] `src/ckks/keys/secret_key.rs`
- [ ] `src/ckks/ops.rs`
- [ ] `src/ckks/types.rs`
- [ ] `src/lib.rs`
- [ ] `src/math.rs`
- [x] `src/math/primes.rs`
- [ ] `src/math/sampling.rs`
- [ ] `src/math/utils.rs`
- [ ] `src/rings.rs`
- [ ] `src/rings/backends.rs`
- [ ] `src/rings/backends/rns_ntt.rs`
- [ ] `src/rings/backends/rns_ntt/basis.rs`
- [ ] `src/rings/backends/rns_ntt/errors.rs`
- [ ] `src/rings/backends/rns_ntt/poly.rs`
- [ ] `src/rings/errors.rs`
- [x] `src/rings/traits.rs`
- [ ] `src/table.rs`
- [ ] `tests/integration_mul.rs`
- [x] `tests/bootstrap.rs`
- [ ] `specs/research/file-review.md`
- [ ] `specs/research/generic-bootstrap.md`
- [ ] `src/ckks/scale.rs`
- [ ] `src/ckks/algorithms/bootstrap/transforms.rs`
- [ ] `src/ckks/algorithms/bootstrap/polynomial.rs`

