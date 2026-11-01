---
id: 0024-kani-evaluation
status: draft
scope: [specs/research/kani.md]
created: 2026-10-31
depends: [0023-creusot-integration-research]
---
# Task: kani evaluation

## Plan

- [ ] Read the official Kani tutorial and support documentation; record installation, host/toolchain requirements and CI costs.
- [ ] Reuse 0023's functions and property definitions; run isolated Kani harnesses with explicit input restrictions and loop bounds.
- [ ] Record commands, versions, runtime, counterexamples, unsupported constructs, unwinding-check results and reachability checks; distinguish failed, incomplete and successful verification.
- [ ] Write `specs/research/kani.md`: compare setup, annotation effort, functional and safety coverage, arithmetic/dependency support, scalability and maintenance with Creusot.
- [ ] Recommend Creusot alone or a complementary Kani role; propose an integration task only if the evidence supports one.

## Context

Research only; depends on 0023's common proof targets and results, not a production Creusot integration. Use disposable harnesses under `tmp/`; retain reproducible snippets and commands in the report.

- Starting source: https://model-checking.github.io/kani/tutorial-first-steps.html . Also consult verification results, Rust support and loop-unwinding documentation.
- State exactly which input domain and loop bounds each result covers. Check for over-constrained or unreachable assertions; do not equate a bounded small-instance result with a general algebraic proof.
- Evaluate whether concrete counterexamples and machine-integer safety checks add enough value to justify a second toolchain. Compare the same properties; explain any properties only one tool can handle.
- Keep arithmetic correctness, CKKS approximation error and cryptographic security as separate claims. Installation failures or unsupported features are evidence to report.

## Acceptance

- [ ] `specdev check`
- [ ] `typos`

## Manual checks

- [ ] Report contains a like-for-like comparison, reproducible results or blockers, explicit proof limits and an evidence-based adoption recommendation.

## Log

- 2026-10-31 created: kani evaluation
- 2026-10-31 scope add specs/research/kani.md — comparison with Creusot
- 2026-10-31 set depends 0023-creusot-integration-research
