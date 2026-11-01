---
id: 0025-ckks-backend-benchmark-research
status: draft
scope: [specs/research/backend-benchmarks.md]
created: 2026-10-31
---
# Task: ckks backend benchmark research

## Plan

- [ ] Audit 0018, 0022 and the sensor standard-deviation example; compare a packed polynomial plus reduction, variance and standard deviation as representative workloads.
- [ ] Recommend the smallest workload covering encode, encrypt, multiply/relinearize, rescale, rotations/reduction, decrypt and decode; specify a plaintext reference and correctness tolerances.
- [ ] Compare Criterion (already installed) with Gungraun (formerly Iai-Callgrind): wall-clock statistics, instruction/cache measurements, platform support, reporting and CI suitability.
- [ ] Design a backend-generic harness with separate evaluation-only and complete-pipeline timing; specify setup, cloning/allocation, key generation, cache warm-up, teardown and correctness-check boundaries.
- [ ] Write `specs/research/backend-benchmarks.md`: parameter matrix, seeded inputs, precision criteria, result metadata, framework recommendation and a separate implementation-task proposal.

## Context

0018 remains the per-operation benchmark implementation task. This task researches complete computations; it neither adds a duplicate operation suite nor implements the harness. Research can proceed now; comparative measurements need the second backend from 0022.

- Criterion: https://bheisler.github.io/criterion.rs/book/ . Alternative: https://github.com/gungraun/gungraun . Check current host support rather than assuming Valgrind is available on the development machine.
- Use the same mathematical parameters, input values, slot layout, operation schedule and accuracy target across backends; allow backend-specific construction outside timed evaluation. Include toy degrees for the schoolbook reference backend and separate practical-size RNS runs without misleading speed ratios.
- `LeveledEvaluator` lazily generates per-level keys. Define warm evaluation with all required keys cached and fresh ciphertext inputs per iteration; measure cold/setup costs separately. Complete-pipeline results must state whether context and key generation are included.
- Record degree, slots, prime chain, scale, decomposition settings, seed, backend, compiler/profile/features, hardware, sampling settings, elapsed time and output error. Reject incorrect results before accepting performance numbers.
- Compare standard deviation's approximation/depth overhead against variance and a simpler circuit. The report must select a workload and framework with reasons, rather than leaving the implementation to choose.
- Add bootstrap timing and a computation that continues after refresh only after bootstrap is available; do not block the initial benchmark design on it.

## Acceptance

- [ ] `specdev check`
- [ ] `typos`

## Manual checks

- [ ] Report selects a workload and framework, defines comparable measurements and accuracy checks, and proposes implementation acceptance commands.

## Log

- 2026-10-31 created: ckks backend benchmark research
- 2026-10-31 scope add specs/research/backend-benchmarks.md — complete CKKS workload and framework research
