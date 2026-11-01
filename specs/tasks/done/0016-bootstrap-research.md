---
id: 0016-bootstrap-research
status: done
scope: [specs/research/bootstrap.md, specs/research/bootstrap_experiment.py, src/rings/**, src/ckks/**, examples/bootstrap_demo.rs, README.md]
created: 2026-10-16
source: operator request 2026-10-16
---
# Task: bootstrap research

Determine how to implement educational CKKS bootstrapping. Explain the design, review choices with the operator and derive ordered implementation tasks. References: fhetextbook D-3.13, Cheon et al. eprint 2018/153, OpenFHE and Lattigo.

## Plan

- [x] Map the selected pipeline to pinned OpenFHE/Lattigo references and explain its mathematical stages.
- [x] Explain each stage in `specs/research/bootstrap.md`: mathematical purpose, coefficient/slot representation, input/output scale and modulus level, keys, operations, depth cost and approximation error.
- [x] Map stages onto repository APIs and gaps; propose minimal backend-neutral interfaces, parameter/key ownership and stage-level plaintext references.
- [x] Compare educational transform and EvalMod variants; recommend a concrete first parameter set, error budget and usable output depth, with evidence from references or experiments.
- [x] Define acceptance tests and resolve the open questions with the operator; derive ordered implementation tasks, with explicit prerequisites and acceptance commands.

## Context

Detailed stage explanations and pinned references live in `specs/research/bootstrap.md`; use ASCII formulas.

1. **Input preparation and ModRaise.** Normalize the input and lift centered component coefficients. Explain the exposed modulus multiples and distinguish lifting from rescale and downward modulus switching.
2. **CoeffToSlot.** Explain diagonals, rotations, plaintext multiplications, conjugation and real/imaginary splitting. Compare dense and optimized transforms, including their key and level costs.
3. **Approximate modular reduction.** Remove the unwanted modulus multiples within a bounded input domain using the selected periodic-function approximation. Explain scaling, polynomial interval/degree, angle doubling and any correction step; a sine expression alone is not a complete EvalMod design. Account separately for approximation error, incoming ciphertext noise and arithmetic error.
4. **SlotToCoeff and output.** Recombine coefficient halves, restore scale/layout and report usable depth. Demonstrate multiplication and rescale after refresh.

Reference entry points (pin revisions during research):

- OpenFHE: https://github.com/openfheorg/openfhe-development/blob/main/src/pke/lib/scheme/ckksrns/ckksrns-fhe.cpp . Trace setup/precomputation, bootstrap key generation and evaluation into their helpers.
- Lattigo: https://github.com/tuneinsight/lattigo/blob/main/circuits/ckks/bootstrapping/evaluator.go . Follow its transform, modular-reduction and parameter packages. Use it for a second explanation and cross-check; reference variants need not have identical stage ordering.
- Follow the project reference policy: OpenFHE and Lattigo, not Microsoft SEAL.

Public bootstrap evaluation uses pre-generated evaluation keys, without decrypting or re-encrypting. Keep secret-key diagnostics outside evaluation.

## Findings

- **Pipeline:** ModRaise (level 1 to the top basis), CoeffToSlot (homomorphic special-DFT, needs rotations and conjugation), EvalMod/EvalExp (scaled sine `q0/(2pi) * sin(2pi x / q0)` removes the `q0*K` term), SlotToCoeff (inverse transform).
- **Exists:** generic key switching and rotation keys (0012), complex encode/decode, Chebyshev evaluation, explicit rescale and levels, `automorphism`, a private coefficient-form flag in `RnsPoly`.
- **Implemented after research:** ModRaise, explicit conjugation, dense complex transforms, balanced Chebyshev evaluation and EvalMod. Sparse-secret switching remains deferred.
- Signed rotations now use true inverse offsets; conjugation is separate. This corrects the former negative-offset rotation-plus-conjugation behavior.
- Rough size: 4-6 tasks. Main risk: precision and parameter tuning.

## Agreed decisions and remaining research

- Operator selected core `PolyRing::mod_raise`, corrected signed rotations plus explicit conjugation, dense diagonal transforms, a shallow Chebyshev evaluator and continuous cosine with double angles.
- First milestone: N=16 then64, full complex packing, fixed-weight ternary secret N/2, slot magnitude<=1 and maximum absolute error1e-3. No security claim; do not wait for 0017/0022.
- Supported implementation: 18 near-2^61 primes, input/output logp44, internal logp61, degree31/four double angles, output level4. The original logp40/17-prime candidate and reasons for tuning remain in `specs/research/bootstrap.md`.
- Encrypted validation covers N=16/64, seeds0..7, six input patterns, multiplication and two refresh cycles. Actual-prime proximity is checked at setup; independent coefficient-packing and polynomial tests check normalization/depth.

## Validation design

- Compare each transform and modular-reduction stage with an independent plaintext reference; include real and complex slots, zero, signs, range endpoints and inputs outside the supported domain.
- Check ModRaise coefficient representatives, shared context identity, invalid levels, rotation direction/conjugation and required-key failures. Round-trip transforms before composing the full pipeline.
- Measure absolute error and precision before/after bootstrap across fixed seeds. State a numeric error tolerance and output-level/scale contract for the chosen parameters; report remaining multiplication depth.
- Include repeated refresh and computation after refresh. Separate educational toy parameters from any security claims; refreshing the modulus budget cannot recover precision already lost.

## Acceptance

- [x] `specdev check`
- [x] `typos`
- [x] `python3 specs/research/bootstrap_experiment.py` — ran scalar sweeps; degree31/double4 sampled restored errors 6.86268e-8 (N=16) and 1.5451e-7 (N=64).

## Manual checks

- [x] Research note explains all stages and maps them to pinned reference code and repository gaps. — agent: reviewed reference map and implemented stage mapping.
- [x] Recommended design specifies parameters, precision/depth targets and interfaces. — agent: release demo and fixed-seed encrypted tests validate the tuned profile; no open implementation choice remains for the supported toy profiles.

## Review

Research decisions are implemented and validated by the bootstrap regression suite; the candidate profile was tuned without weakening the 1e-3 target or disabling noise. See implementation tasks0026-0031 for code handoff. No production security claim.

## Summary

Designed and validated educational full-complex bootstrap profiles and the implementation sequence.

## Log

- 2026-10-16 created: bootstrap research
- 2026-10-16 set source operator request 2026-10-16
- 2026-10-16 scope add specs/research/bootstrap.md — research note
- 2026-11-01 draft approved: Implement the plan.
- 2026-11-01 advance draft → ready
- 2026-11-01 scope approved: specs/research/bootstrap.md
- 2026-11-01 advance ready → in-progress/implement
- 2026-11-01 scope add specs/research/bootstrap_experiment.py — reproducible plaintext approximation experiment
- 2026-11-01 scope add src/rings/** — operator approved complete bootstrap implementation while research remains active
- 2026-11-01 scope add src/ckks/** — operator approved complete bootstrap implementation while research remains active
- 2026-11-01 scope add examples/bootstrap_demo.rs — operator approved complete bootstrap implementation while research remains active
- 2026-11-01 scope add README.md — document completed bootstrap example and supported profile
- 2026-11-01 scope add tests/** — bootstrap regression coverage
- 2026-11-01 advance in-progress/implement → in-progress/verify
- 2026-11-01 advance in-progress/verify → in-progress/review
- 2026-11-01 advance in-progress/review → approval
- 2026-11-01 scope rm tests/** — regressions live inside scoped source modules; no external test files changed
- 2026-11-01 approved: ok, looks good.
- 2026-11-01 advance approval → done
