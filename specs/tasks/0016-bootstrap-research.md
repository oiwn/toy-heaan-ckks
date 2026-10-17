---
id: 0016-bootstrap-research
status: draft
scope: [specs/research/bootstrap.md]
created: 2026-10-16
source: operator request 2026-10-16
---
# Task: bootstrap research

What would CKKS bootstrapping take in this codebase? Reference: fhetextbook D-3.13 (Modulus Bootstrapping), Cheon et al. eprint 2018/153. Open-ended: the operator makes the design decisions.

## Plan

- [ ] Read D-3.13 and its subsections in full; note anything the findings below miss.
- [ ] Write `specs/research/bootstrap.md`: the pipeline, what exists, what is missing, and the open questions below.
- [ ] Operator answers the open questions; split the result into implementation tasks.

## Findings

- **Pipeline:** ModRaise (level 1 to the top basis), CoeffToSlot (homomorphic special-DFT, needs rotations and conjugation), EvalMod/EvalExp (scaled sine `q0/(2pi) * sin(2pi x / q0)` removes the `q0*K` term), SlotToCoeff (inverse transform).
- **Exists:** generic key switching and rotation keys (0012), complex encode/decode, Chebyshev evaluation, explicit rescale and levels, `automorphism`, a private coefficient-form flag in `RnsPoly`.
- **Missing:** a ModRaise operation; a conjugation key; homomorphic linear transforms; EvalMod (sine plus angle doubling, mul by `i`, real/imag split); a deep parameter set with a sparse secret.
- `rotate_slots(-k)` today is "rotate by k, then conjugate", not a left rotation; on real data it equals `rotate(+k)`.
- Rough size: 4-6 tasks. Main risk: precision and parameter tuning.

## Open questions

- Trait shape for ModRaise (new `PolyRing` method or backend-specific).
- Negative-rotation and conjugation semantics in `rotate_slots`.
- Linear-transform method (dense diagonal, BSGS, FFT-factored).
- EvalMod variant (Taylor + double angle, Chebyshev sine, arcsin correction).
- Target N, chain length, secret sparsity, accepted precision and security.
- Whether it waits for 0017 (debug tooling) and the M2.6 reference backend.

## Acceptance

- [ ] `specdev check`
- [ ] `specs/research/bootstrap.md` covers pipeline, gaps and open questions

## Log

- 2026-10-16 created: bootstrap research
- 2026-10-16 set source operator request 2026-10-16
- 2026-10-16 scope add specs/research/bootstrap.md — research note
