---
id: 0023-creusot-integration-research
status: draft
scope: [specs/research/creusot.md]
created: 2026-10-31
---
# Task: creusot integration research

## Plan

- [ ] Read the official Creusot guide; record versions, installation, supported hosts, Rust toolchain requirements and CI costs.
- [ ] Audit compatibility with this crate: edition 2024, const generics, `PolyRing`, mutable borrows, `Arc`, randomness, big integers, floats and external dependencies.
- [ ] Select two existing proof targets: a modular arithmetic helper and a small polynomial/indexing operation; define preconditions, functional postconditions and overflow/index safety properties.
- [ ] Run isolated proof experiments against those targets; record exact commands, contracts, solver/tool versions, results, runtime, unsupported constructs and trusted assumptions.
- [ ] Write `specs/research/creusot.md`: recommend an integration boundary, annotation strategy, proof maintenance and optional CI approach; propose a small first integration task with acceptance commands.

## Context

Research only; no production dependency, API or CI changes. Use disposable experiments under `tmp/` and retain reproducible source snippets and commands in the report. An installation or proof failure is a finding, not a reason to claim success or silently substitute a different algorithm.

- Sources: https://creusot.rs and https://guide.creusot.rs/ . Pin the versions actually evaluated.
- Candidate starting points: modular arithmetic in `src/math/` and coefficient indexing in the RNS-NTT backend. Select actual functions after auditing their assumptions; record any extraction required to verify existing behavior.
- Distinguish proven contracts from assumed dependency behavior, floating-point error bounds, randomness properties and cryptographic security. Do not promise a proof of the entire CKKS scheme.
- Share target functions and property definitions with 0024 so the comparison measures the same work.
- Decision: is incremental integration practical, what can be proved first, and what code/toolchain changes would it require?

## Acceptance

- [ ] `specdev check`
- [ ] `typos`

## Manual checks

- [ ] Report includes reproducible experiments or concrete blockers, a compatibility assessment and a justified first integration proposal.

## Log

- 2026-10-31 created: creusot integration research
- 2026-10-31 scope add specs/research/creusot.md — integration feasibility report
