---
id: 0013-sage-education-scripts
status: draft
scope: [scripts/**]
created: 2026-10-14
source: user request 2026-10-14: revive educational sage scripts demonstrating properties like Chebyshev polynomials
---
# Task: sage education scripts

## Plan

Small, self-contained Sage scripts in `scripts/sage/`, each printing a short explanation plus numbers that can be checked against the Rust code. Old scripts in `tmp/scripts/` are the starting point. Topic list is a proposal; the operator picks.

- [ ] `chebyshev.sage`: `T_k` recurrence and `cos` identity, Chebyshev nodes, interpolant vs power basis (conditioning), Clenshaw evaluation, sqrt error on the sensor envelope `[1e-4, 6.4e-3]` per degree (Rust reference: degree 7 -> 0.055 degC at W = 100).
- [ ] `canonical_embedding.sage`: `X^N + 1` roots, slots and conjugate symmetry, encode/decode as the special DFT, why a scalar must fill all slots (constant polynomial).
- [ ] `rotation.sage`: automorphism `X -> X^(5^k)` rotates slots; conjugation via `X -> X^-1`.
- [ ] `rescale.sage`: scale growth after mul, rescale by the last prime, mod switch vs rescale on a toy RNS chain.
- [ ] `scripts/README.md`: what each script shows, how to run (`sage scripts/sage/<name>.sage`, `SAGE_DOT_SAGE` workaround).

## Acceptance

- [ ] `sage scripts/sage/chebyshev.sage` (and every other script) runs without errors
- [ ] Printed reference numbers match the Rust side (Chebyshev coefficients vs `fit_chebyshev`, encode/decode vs `CkksEncoder` at N = 8)
- [ ] `typos`
- [ ] `gitleaks dir .`

## Context

Independent of 0006-0012 (no Rust changes). ASCII-only comments and names (AGENTS.md).

## Log

- 2026-10-14 created: sage education scripts
- 2026-10-14 set source user request 2026-10-14: revive educational sage scripts demonstrating properties like Chebyshev polynomials
- 2026-10-14 scope add scripts/** — educational sage scripts and their README
