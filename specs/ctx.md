# Current Task Context: Publication review: architecture, code, docs

State: not started. Roadmap M2.5.

## Plan

- [ ] Task 0005-remove-dead-code (dead modules, golden vectors, scripts).
- [ ] Task 0006-trim-dependencies.
- [ ] Agree the backend design (Findings), move it to `specs/overview.md`, then revise 0007-0010 around it.
- [ ] Task 0007-module-layout (move to the target layout below, `foo.rs` + `foo/`).
- [ ] Task 0008-ciphertext-levels (`Ciphertext::rescale`/`mod_switch_to`, drop dead traits).
- [ ] Task 0009-ops-composition (`ckks/ops.rs`, `ckks/evaluator.rs`, scalar helpers gone).
- [ ] Task 0010-chebyshev-module (`ckks/chebyshev.rs`, `pipeline.rs` retired).
- [ ] Audit what remains: findings below, terse, each with severity and proposed fix; operator marks them with `^^^`.
- [ ] Split accepted findings into task files (expected: legacy keys, engine generic over `P`, `tracing` logging, errors vs panics, comment/ASCII pass, examples, README + CI).

## Findings

- **Rescale today** (divide by last prime, drop it; `logp` and `logq` shrink): inherent `RnsPoly::rescale_into`/`rescale` (exact RNS division) wrapped by `CkksEngine::rescale_ciphertext` (shared basis `Arc`). Trait `PolyRescale::rescale_assign(f64)` is a naive-backend leftover: no impls, only user is dead `operations.rs`.
- **Mod switch today** (drop primes, no division; `logp` stays): `pipeline::align_level(ct, basis, logq)` truncates channels via private `reduce_poly_to_basis`; `reduce_secret_key` does the same for keys; caller passes `logq`. Trait `PolyModSwitch` (textbook `round(c * q'/q)`) has no impls, so generic `Ciphertext::mod_switch` in `types.rs` is uncallable.
- **Reference composition:** fhe.rs = `fhe-math` (zq, rns, ntt, rq) + scheme crate (`bfv/{parameters, plaintext, ciphertext, encoding, keys/, ops}`); tfhe-rs `core_crypto` = `entities/` + `algorithms/` + `commons/`.
- **Design goal (operator):** experimental library, pluggable ring backends, math abstracted behind a few operations; borrow concepts from fhe.rs / tfhe-rs.
- **Leak today:** `CkksEngine<P>` is generic only up to add; rescale, mul_plaintext, gadget keys and key switch are `impl CkksEngine<RnsPoly>` walking `moduli()`/`channels()[i]` and flipping NTT flags; the encoder and `pipeline.rs` also read RNS channels. The domain-flag and `Arc::ptr_eq` gotchas come from this.
- **Proposed backend trait** (one trait on `Z_Q[X]/(X^N + 1)` at a level, name open):
  - arithmetic: `zero`, `add`, `sub`, `neg`, `mul` (negacyclic), `mul_small(i64)`;
  - sampling: `uniform`, `gaussian(std_dev)`, `ternary(hamming_weight)`;
  - automorphism: `automorphism(k)` (`X -> X^k`), rotation builds on it;
  - modulus chain: `level()`, `rescale()` (divide by last prime), `mod_switch_down(level)` (drop primes); replaces `PolyRescale`/`PolyModSwitch`;
  - encoding boundary: `from_rounded_f64(&[f64; N], level)`, `to_centered_f64()`; the RNS residue trick moves into the backend;
  - gadget: `decompose(digit_bits)`, `gadget_scale(i, digit_bits)`; one generic `KeySwitchingKey` for relin and rotation.
  - NTT vs coeff domain is private to the backend.
- **Borrowed:** fhe.rs context chain (every level's `Arc<Context>` built once, `next()`; kills the fresh-`Arc` bug class), one `KeySwitchingKey` wrapped by relin/galois keys, level changes as `Ciphertext` methods, `std::ops` Add/Sub/Neg; tfhe-rs entities vs algorithms, typed parameter newtypes (`DecompositionBaseLog`, `LevelCount`), decomposition as its own module.
- **Layering:** `math/` (mod arithmetic, primes, distributions, NTT tables) -> `rings/` (trait, context chain, `backends/{rns_ntt, reference}`) -> `ckks/` generic over the backend (`params, encoding, plaintext, ciphertext, keys/, ops, evaluator`) -> `ckks/algorithms` (`sum_slots`, Chebyshev, only on ops); std-dev stays in `examples/`. `reference` = single big modulus via `crypto-bigint`, schoolbook mul; same scheme tests on both backends.
- **Open questions:** keep `const N` generics or put the degree in the runtime context (fhe.rs)? Reference backend in M2.5 or later?
- **Audit candidates:** legacy `RelinearizationKey`/`RotationKey` next to the gadget keys in `engine.rs`; logging: `tracing` spans/events in engine and pipeline, subscriber in examples; `table.rs` (comfy-table, used by 3 examples only); `crypto/builder.rs`; workaround re-exports in `lib.rs`; overlapping examples (`std_dev_8`, `sqrt_chebyshev_plain`); 179 non-ASCII lines in `src/`.

## Context

Drafts 0005-0010 in `specs/tasks/`, chained by `depends`; 0007-0010 are on hold until the backend design is agreed.

## Next

Operator reviews the backend design and open questions (`^^^` here) and drafts 0005-0006; then advance 0005 to ready.
