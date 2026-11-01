---
id: 0037-cleanup-and-refactoring
status: draft
scope: [specs/cleanup.md, specs/research/file-review.md]
created: 2026-11-01
source: Operator request: cleanup and refactoring; carefully check every file
---
# Task: cleanup and refactoring

## Plan

- [ ] Inventory repository-owned files, including tracked and relevant untracked files; record a per-file review checklist in `specs/research/file-review.md`. Exclude generated/build output, tmp and tool internals; record exclusions explicitly.
- [ ] Review every file with the operator: purpose, correctness/invariants, API clarity, module boundaries, duplication/dead code, error handling, numerical assumptions, naming/ASCII, documentation and test coverage. Mark reviewed files and unresolved remarks without treating an unreviewed file as clean.
- [ ] Reconcile existing audit candidates in ctx.md and specs/cleanup.md against current code. Include Rust source, examples, tests, benches, manifests/dependencies, CI/hooks, project configuration and documentation; identify stale or contradictory specs.
- [ ] Agree each proposed refactor with the operator, group coherent changes into small implementation tasks with exact scope/dependencies/checks, then implement through those tasks. Extend scope through specdev before editing code; do not launch a blanket rewrite.
- [ ] Revisit affected files after changes and verify behavior, examples and documentation; resolve or explicitly defer each finding before final operator acceptance.

## Context

Next workstream after accepted bootstrap tasks0016 and0026-0033/0036. Start with review and a file checklist, not automatic edits. This task coordinates the review and scoped refactoring work; implementation tasks carry code changes. Initial scope covers audit artifacts only. Do not edit accepted task history or remove unrelated uncommitted work.

Preserve working CKKS behavior and bootstrap accuracy/depth regressions. Treat numerical changes (including the variance/std-dev mismatch), parameter changes and public API changes as explicit decisions with separate validation. Use OpenFHE/Lattigo as implementation references; do not use SEAL. Do not silently discard historical findings or operator remarks. File removal follows AGENTS.md: move candidates to tmp for the operator, never delete them.

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
