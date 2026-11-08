# Cleanup (code smells)

Read on demand — not session context.

## Encryption noise policy is inconsistent

Confirmed 2026-11-02: `CkksParams::error_variance` is passed directly to Gaussian sampling during encryption; relin/rotation key generation takes its square root, while public-key generation hardcodes 3.2. The previous "all call sites" observation was incomplete. A rename alone cannot resolve this; agree intended noise semantics and validate numerical behavior separately.

## Bootstrap readability

Structural pass implemented: private validation/diagonal helpers, explicit type names, documented key indices, end-to-end tests under `tests/` and small private stage helpers at the bottom of the source. The level diagnostic now says "exactly 18". Await operator review; details and measurements are in `specs/research/file-review.md`.

## Generic bootstrap parameters

General encrypted bootstrap remains incomplete. Scale, factored transforms and a polynomial compiler exist, but ciphertext metadata, encrypted graph execution and derived bootstrap configuration/keys are not integrated. Current code still enforces N=16/64, scales44/61, degree31, four double angles and18 near-2^61 primes. The target replaces these fixed-profile checks with mathematical/backend validation. "Toy" indicates experimental status, not functional scope. The authoritative remaining-work and verification checklist is `specs/research/generic-bootstrap.md`.

## Engine builder ignores its scale setter

Found in 0039: `CkksEngineBuilder::scale(Scale)` stores a value that `build(context, scale)` never reads, because `build` takes the scale as an argument. The builder is otherwise unused. Either drop the setter or make `build` use it (and require it), in a separate change.

## Recurrence

Do not split `step` solely for its complexity score: guards and error propagation explain much of it. Consider named helpers for the 127-line combined test scenario when reviewing tests.
