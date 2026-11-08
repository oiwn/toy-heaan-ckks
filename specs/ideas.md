# Ideas (uncommitted backlog)

Not-now possibilities. Promote to the active plan in `ctx.md` (or a roadmap) when we decide to do them.

- Repair `examples/sqrt_chebyshev_plain.rs` (negative LOWER → NaN).
- Remove the non-optimized u128 modular math from the library (`math::modular::mul_mod`, `pow_mod_reference`, the remaining `% q` products in `primes.rs` and the CRT) and use Shoup/Barrett everywhere (operator, 2026-11-07). The u128 version stays only as the explanatory reference in `examples/modmul_methods.rs` and as the test oracle. Prerequisite: Miller-Rabin must handle moduli >= 2^62, where Barrett is capped.
- Replace the per-level gadget-key regeneration in `LeveledEvaluator` with keys generated once at the top level and modulus-switched down.
