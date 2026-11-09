---
id: 0053-scheme-step-examples
status: done
scope: [examples/conjugation.rs, examples/mod_switch.rs, examples/key_switch.rs, examples/bootstrap_stages.rs, README.md, src/ckks/algorithms/bootstrap.rs, .agents/skills/educational-examples/SKILL.md, .agents/.gitignore, .agents/skills/.gitignore, examples/support/narrate.rs, AGENTS.md, examples/modmul_methods.rs]
created: 2026-11-08
source: Operator request 2026-11-08: examples for scheme steps not yet covered (conjugation, modulus switching, key switching, bootstrap stages)
grade: huge
---
# Task: scheme step examples

## Context

Educational goal: every CKKS step has its own runnable example. These steps are exercised only inside other demos or tests: conjugation, modulus switching, key switching, and the bootstrap stages. Each example prints inputs, the step's effect (level, scale, slots) and the error against a plaintext f64 reference, and gates on a tolerance like the existing demos.

Library code: `Ciphertext::conjugate` (`ckks/ops.rs`), `mod_switch_to`/`match_level` (`ckks/types.rs`), `KeySwitchingKey::switch` (`keys/key_switching.rs`). Bootstrap building blocks (`bootstrap::{linear, polynomial, executor, transforms}`) are public; the stage sequence in `BootstrapPlan` (`prepare_input`, `transform`, `eval_mod`, `multiply_by_i`) is private.

Decided 2026-11-08: `BootstrapPlan` exposes a public stage-by-stage API (ModRaise, CoeffsToSlots, EvalMod, SlotsToCoeffs) and `bootstrap` itself is composed from those stages, so the sequence lives in one place (also serves 0046's stage trace).

Revision 2026-11-09 (operator, issue #39): examples must be educational, explain their methods and use colored output. Add a committed project skill `.agents/skills/educational-examples/` (the other `.agents/skills` stay ignored), shared helpers `examples/support/narrate.rs` (from `modmul_methods.rs`, the reference style), and rework the four new examples to it. Grade huge: the operator added the skill to this task. Reworking the other examples follows as 0054.

Rules: touch only scoped files. Never weaken tolerances or disable noise.

## Plan

- [x] `conjugation`: encrypt complex slots, conjugate with a `ConjugationKey`, show `z -> conj(z)` and the real/imag split `(z + conj(z)) / 2`.
- [x] `mod_switch`: drop levels with `mod_switch_to` versus `rescale`; show that mod switch keeps the scale and rescale divides it; `match_level` before adding ciphertexts at different levels.
- [x] `key_switch`: build a `KeySwitchingKey` from `s_old` to `s_new`, switch a ciphertext and decrypt under the new key; print the noise the switch adds against `DIGIT_BITS` (a library constant, so the digit width is explained, not varied).
- [x] `BootstrapPlan` public stages: `prepare`, `mod_raise`, `coeffs_to_slots`, `eval_mod`, `slots_to_coeffs` with docs; `bootstrap` calls them; results bit-identical.
- [x] `bootstrap_stages`: at N=16, run the stages one by one; after each, decrypt and print level, scale and error against the plaintext pipeline.
- [x] README Examples table: add the four rows under their steps; keep it one screen.
- [x] Skill `educational-examples` (shape, color meanings, text, correctness, checklist); `.agents` gitignores let only it through; local `.claude/skills` symlink; AGENTS.md pointer.
- [x] `examples/support/narrate.rs` (heading, subheading, formula, value, ok, bad, dim, step, param, say, check); `modmul_methods` uses it.
- [x] Rework the four examples to the skill: narrated intro, parameters with reasons, worked steps, colored checks, takeaways.

## Acceptance

- [x] `cargo fmt --check`
- [x] `cargo clippy --all-targets --all-features -- -D warnings`
- [x] `cargo test --all-features`
- [x] `cargo build --examples`
- [x] `cargo run --example conjugation`
- [x] `cargo run --example mod_switch`
- [x] `cargo run --example key_switch`
- [x] `cargo run --release --example bootstrap_stages`
- [x] `typos`
- [x] `specdev check`
- [x] `git diff --check`

## Findings

- Stage API: `prepare(input)` (lift + level 1; no keys), `mod_raise(prepared)`, `coeffs_to_slots(raised, keys) -> (real, imaginary)`, `eval_mod(half, keys)`, `slots_to_coeffs(real, imaginary, prepared_scale, keys)`. Replaces private `prepare_input`/`bootstrap_prepared`; keyed stages reject keys of another plan. A unit test asserts the staged result equals `bootstrap` coefficient for coefficient; `bootstrap_demo` errors unchanged (N=16 1.235e-5 / 1.520e-5, N=64 1.222e-4 / 2.629e-4).
- `bootstrap_stages` (N=16): prepare 1.2e-11; ModRaise `max|I|=2` (bound 5), exact integers; CtS halves relative 8.4e-10; EvalMod relative 3.0e-4 (values ~2e-5, the sine approximation); output 8.3e-6 at level 7.
- Skill plumbing: `git status` shows only `.agents/.gitignore`, `.agents/skills/.gitignore` and the new `SKILL.md`; installed skills stay ignored. Colors verified in a pty; plain with `NO_COLOR=1` or piped.
- `key_switch`: 12 gadget terms at logq 180, added noise 2^23.8 per coefficient (DIGIT_BITS 20), switched error 3.8e-8. `mod_switch`: errors <= 1.0e-10; the unrescaled product mod-switched to level 1 decodes to garbage (0.81), shown on purpose. `conjugation`: 1.3e-7.

## Review

Diff matches the plan and scope. The library change only re-cuts the existing sequence into public methods; arithmetic is unchanged and the equality test plus identical demo errors confirm it. Stage docs state each stage's slot contents (fixed a factor 2 in the first draft: CtS halves hold `c / (radius * q0)`). Examples follow the `educational-examples` skill: narrated sections, parameters with reasons, colors only through `narrate.rs`, failures shown on purpose in red with the reason, nonzero exit on a missed bound, ASCII only. `modmul_methods` output unchanged (same helper code, moved).

## Summary

Four new narrated, colored step examples (`conjugation`, `mod_switch`, `key_switch`, `bootstrap_stages`), a public stage-by-stage bootstrap API (`BootstrapPlan::{prepare, mod_raise, coeffs_to_slots, eval_mod, slots_to_coeffs}`) that `bootstrap` is composed from, and the committed `educational-examples` skill with shared helpers `examples/support/narrate.rs` (closes issue #39).

## Log

- 2026-11-08 created: scheme step examples
- 2026-11-08 set source Operator request 2026-11-08: examples for scheme steps not yet covered (conjugation, modulus switching, key switching, bootstrap stages)
- 2026-11-08 scope add examples/conjugation.rs — new step examples, README table rows, possible public stage API
- 2026-11-08 scope add examples/mod_switch.rs — new step examples, README table rows, possible public stage API
- 2026-11-08 scope add examples/key_switch.rs — new step examples, README table rows, possible public stage API
- 2026-11-08 scope add examples/bootstrap_stages.rs — new step examples, README table rows, possible public stage API
- 2026-11-08 scope add README.md — new step examples, README table rows, possible public stage API
- 2026-11-08 scope add src/ckks/algorithms/bootstrap.rs — new step examples, README table rows, possible public stage API
- 2026-11-08 draft approved: approved by operator
- 2026-11-08 advance draft → ready
- 2026-11-08 scope approved: examples/conjugation.rs, examples/mod_switch.rs, examples/key_switch.rs, examples/bootstrap_stages.rs, README.md, src/ckks/algorithms/bootstrap.rs
- 2026-11-08 advance ready → in-progress/implement
- 2026-11-09 advance in-progress/implement → in-progress/verify
- 2026-11-09 advance in-progress/verify → in-progress/review
- 2026-11-09 advance in-progress/review → approval
- 2026-11-09 advance approval → in-progress/fix (revision)
- 2026-11-09 scope add .agents/skills/educational-examples/SKILL.md — issue #39: educational-examples skill, shared narration helpers, commit plumbing
- 2026-11-09 scope add .agents/.gitignore — issue #39: educational-examples skill, shared narration helpers, commit plumbing
- 2026-11-09 scope add .agents/skills/.gitignore — issue #39: educational-examples skill, shared narration helpers, commit plumbing
- 2026-11-09 scope add examples/support/narrate.rs — issue #39: educational-examples skill, shared narration helpers, commit plumbing
- 2026-11-09 scope add AGENTS.md — issue #39: educational-examples skill, shared narration helpers, commit plumbing
- 2026-11-09 scope add examples/modmul_methods.rs — issue #39: educational-examples skill, shared narration helpers, commit plumbing
- 2026-11-09 set grade huge
- 2026-11-09 advance in-progress/fix → in-progress/verify
- 2026-11-09 advance in-progress/verify → in-progress/review
- 2026-11-09 advance in-progress/review → approval
- 2026-11-09 approved: approved by operator
- 2026-11-09 advance approval → done
