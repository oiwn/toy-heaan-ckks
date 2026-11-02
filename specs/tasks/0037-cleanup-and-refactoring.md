---
id: 0037-cleanup-and-refactoring
status: blocked
scope: [specs/cleanup.md, specs/research/file-review.md, src/ckks/algorithms/bootstrap.rs, tests/bootstrap.rs, src/ckks/scale.rs, src/ckks.rs, src/ckks/algorithms/bootstrap/transforms.rs, src/ckks/algorithms/bootstrap/polynomial.rs, src/ckks/errors.rs]
created: 2026-11-01
source: Operator request: cleanup and refactoring; carefully check every file
grade: huge
blocked_reason: waiting on bootstrap chain 0038-0043
---
# Task: cleanup and refactoring

## Plan

- [x] Inventory repository-owned files, including tracked and relevant untracked files; record a per-file review checklist in `specs/research/file-review.md`. Exclude generated/build output, tmp and tool internals; record exclusions explicitly.
- [ ] Review every file with the operator: purpose, correctness/invariants, API clarity, module boundaries, duplication/dead code, error handling, numerical assumptions, naming/ASCII, documentation and test coverage. Mark reviewed files and unresolved remarks without treating an unreviewed file as clean.
- [ ] Reconcile existing audit candidates in ctx.md and specs/cleanup.md against current code. Include Rust source, examples, tests, benches, manifests/dependencies, CI/hooks, project configuration and documentation; identify stale or contradictory specs.
- [ ] Agree each proposed refactor with the operator, group coherent changes into small implementation tasks with exact scope/dependencies/checks, then implement through those tasks. Extend scope through specdev before editing code; do not launch a blanket rewrite.
- [ ] Revisit affected files after changes and verify behavior, examples and documentation; resolve or explicitly defer each finding before final operator acceptance.

## Context

Active repository review and scoped refactoring workstream after accepted bootstrap tasks0016 and0026-0033/0036. The inventory exists; most files remain unreviewed. Initial scope covered audit artifacts; operator-approved bootstrap implementation expanded it through the CLI. Do not edit accepted task history or remove unrelated uncommitted work.

Preserve working CKKS behavior and bootstrap accuracy/depth regressions. Treat numerical changes (including the variance/std-dev mismatch), parameter changes and public API changes as explicit decisions with separate validation. Use OpenFHE/Lattigo as implementation references; do not use SEAL. Do not silently discard historical findings or operator remarks. File removal follows AGENTS.md: move candidates to tmp for the operator, never delete them.

## Findings

- General encrypted bootstrap is incomplete. Implemented components, active restrictions, numerical rules and verification evidence live in `specs/research/generic-bootstrap.md`.
- 2026-11-02 operator decision: bootstrap completion moved out of this task into the dependency chain 0038-0043 (fallible arithmetic, Scale type migration with nominal rescale, executor with actual-divisor rescale, config/schedule, multi-group transforms, final matrix). The bootstrap source paths in this task's scope belong to that work now; this task returns to the file review after 0043 is accepted.
- Completed component checks do not establish acceptance of this umbrella task. Leave Plan, Manual checks and final Acceptance open until their full scope is actually satisfied.

## Manual checks

- [ ] Every inventoried file is reviewed or explicitly excluded with a reason; findings have an agreed disposition.
- [ ] Operator accepts the review and implemented refactors, including any explicitly deferred findings.

## Acceptance

- [ ] `cargo fmt --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features`
- [ ] `cargo build --examples`
- [ ] `specdev check`
- [ ] `typos`
- [ ] `git diff --check`

## Log

- 2026-11-01 created: cleanup and refactoring
- 2026-11-01 scope add specs/cleanup.md — file-by-file audit and scoped refactoring backlog
- 2026-11-01 scope add specs/research/file-review.md — file-by-file audit and scoped refactoring backlog
- 2026-11-01 set source Operator request: cleanup and refactoring; carefully check every file
- 2026-11-02 draft approved: Implement the plan.
- 2026-11-02 advance draft → ready
- 2026-11-02 scope approved: specs/cleanup.md, specs/research/file-review.md
- 2026-11-02 advance ready → in-progress/implement
- 2026-11-02 scope add src/ckks/algorithms/bootstrap.rs — Approved bootstrap structural cleanup and end-to-end test extraction
- 2026-11-02 scope add tests/bootstrap.rs — Approved bootstrap structural cleanup and end-to-end test extraction
- 2026-11-02 scope add src/ckks/scale.rs — Approved general bootstrap plan: validated real-valued scale foundation for clean API migration
- 2026-11-02 scope add src/ckks.rs — Approved general bootstrap plan: validated real-valued scale foundation for clean API migration
- 2026-11-02 scope add src/ckks/algorithms/bootstrap/transforms.rs — General bootstrap implementation: sparse canonical FFT factors and configurable level grouping
- 2026-11-02 set grade huge
- 2026-11-02 scope add src/ckks/algorithms/bootstrap/polynomial.rs — Approved bootstrap completion: target-scale polynomial operation graph and scale errors
- 2026-11-02 scope add src/ckks/errors.rs — Approved bootstrap completion: target-scale polynomial operation graph and scale errors
- 2026-11-02 blocked: waiting on bootstrap chain 0038-0043
