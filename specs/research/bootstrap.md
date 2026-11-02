# Educational CKKS bootstrap design

## Decisions and evidence status

Operator choices: toy first, full complex packing, slot magnitude at most 1, maximum absolute slot error `1e-3`; extend `PolyRing`, fix signed rotations, use a shallow polynomial evaluator, and do not wait for 0017 or 0022. Use dense diagonal transforms and continuous Chebyshev cosine plus double angles. Defer arcsine correction, sparse-secret switching and performance optimization. No security claim for these toy parameters.

The scalar experiment below ran successfully. The implemented profile now uses 18 primes and input/output logp44; the original 17-prime/logp40 proposal below is retained as the derivation, not the supported profile. Encrypted two-cycle validation is implemented at both toy degrees.

## Pinned reference map

- [Lattigo evaluator](https://github.com/tuneinsight/lattigo/blob/7abb8d5dec86c6cc609e29540ddfc0c16864e14b/circuits/ckks/bootstrapping/evaluator.go): `NewEvaluator`, `checkKeys`, `initialize`, `bootstrap`, `ScaleDown`, `ModUp`, `CoeffsToSlots`, `EvalMod`, `SlotsToCoeffs`. In this revision the pipeline is in evaluator.go, not bootstrapping.go. `initialize` folds normalization into transform matrices; `checkKeys` separates key availability from evaluation.
- [Lattigo modular reduction](https://github.com/tuneinsight/lattigo/blob/7abb8d5dec86c6cc609e29540ddfc0c16864e14b/circuits/ckks/mod1/mod1_evaluator.go): `EvaluateAndScaleNew` schedules polynomial scales and double angles; `EvaluateNew` documents modulus-ratio correction. Our plain `2*y*y-1` recurrence keeps amplitude outside the recurrence; Lattigo can fold amplitude into polynomial coefficients instead.
- [OpenFHE bootstrap](https://github.com/openfheorg/openfhe-development/blob/b59e1f291965f14fbee5172cd83a20f5b856204f/src/pke/lib/scheme/ckksrns/ckksrns-fhe.cpp): cross-check `EvalBootstrapSetup`, `EvalBootstrapKeyGen`, `EvalBootstrap`, `GetModRaiseOverflowBound` and the linear-transform helpers. Its optimized key switching, scale handling and secret distributions are not drop-in replacements for this repository.
- [FHE textbook D-3.13](https://fhetextbook.github.io/ModulusBootstrapping.html): mathematical explanation of modulus multiples and approximate periodic reduction. Reference code guides implementation; derive matrices from our encoder, not another library's slot conventions.

## Mathematical stages

Write the bottom modulus as `bottom_modulus`, input scale as `input_scale`, secret weight as `weight`, and centered encoded coefficients as `message_coeff`. At the bottom, decryption is `input_scale * message_coeff + noise` modulo `bottom_modulus`.

1. Center each ciphertext component coefficient in `(-bottom_modulus/2, bottom_modulus/2]`, then reduce those signed integers into every target prime. Decryption at the raised modulus is now `input_scale * message_coeff + noise + bottom_modulus * integer_wrap`. Raising restores storage capacity, not correctness or precision. This operation is neither rescale (division) nor ordinary downward modulus switching.
2. Interpret the lifted ciphertext at working scale, initially `2^61`. Apply the dense transform that exposes coefficient halves as complex slots; split into two real ciphertexts with conjugation. Include `2^61 / bottom_modulus` in the transform. Each branch then contains `integer_wrap + message_ratio * message_coeff + normalized_noise`, with `message_ratio = input_scale / bottom_modulus`.
3. For each real branch, set `x` to that normalized value. Approximate `cos(2*pi*(x-0.25)/2^4)` with a degree-31 Chebyshev polynomial. Four repetitions of `y = 2*y*y-1` yield approximately `sin(2*pi*x)`. Divide by `2*pi`; near integers this approximates the small fractional message. A sine alone introduces cubic small-message error; polynomial approximation, integer coverage and restored scale all matter.
4. Recombine the two real coefficient halves and apply the inverse linear transform. Fold `bottom_modulus / input_scale` into its diagonals. This restores the original slot values and yields the requested output scale using the schedule below.

For centered component lifts, each coefficient of their integer decryption has magnitude at most `(weight+1)*bottom_modulus/2`. Assuming centered noisy message coefficients remain below `bottom_modulus/2`, `abs(integer_wrap) <= ceil((weight+2)/2)` is a conservative bound. Use 5 for N=16/weight=8 and 17 for N=64/weight=32. These are deterministic bounds under the stated message/noise condition, not transplanted statistical bounds from a reference library.

Normalized slots imply bounds on canonical coefficients through the inverse embedding. Confirm the encoder's normalization independently, then derive its coefficient bound; never equate slot bounds with an observed secret-key diagnostic. Hidden input-domain violations cannot be rejected by inspecting a ciphertext.

## Interfaces and ownership

- Add `PolyRing::mod_raise(&self, target: &Self::Context) -> RingResult<Self>`. V1 accepts only level-1 input and a higher target in the same original chain. RNS checks the target's bottom context identity, lifts in coefficient form and returns NTT form. Reject same/lower/unrelated targets. A later backend must implement the same centered-representative contract.
- Signed rotations normalize offsets modulo `N/2` before computing the automorphism exponent. Conjugation uses exponent `2*N-1`, its own key and ciphertext operation. Opposite rotations must compose to identity on complex data.
- `BootstrapParameters` describes degree, scales, levels, polynomial settings and domain bounds. `BootstrapPlan` owns encoded diagonals and an operation/scale schedule. `BootstrapKeys` owns level-specific relin, rotation and conjugation keys. Construction validates required keys and context compatibility.
- Key generation receives the secret key and RNG. `bootstrap(ciphertext, plan, keys)` returns `CkksResult<Ciphertext>` and receives neither. Leave existing `LeveledEvaluator` behavior available for existing callers; do not route public bootstrap evaluation through it.
- Missing keys, wrong context/scale/level and insufficient headroom return errors. Approximation-domain membership is a caller precondition, checked only by plaintext references or test diagnostics.

## Candidate profile and scale ledger

Use N=16 first and N=64 second, all N/2 complex slots, fixed-weight ternary secrets with weight N/2. Keep current gadget width 20 initially; measure its accumulated error before changing it. Use the existing key-noise convention explicitly; the audit's encryption variance/std-dev mismatch must be recorded, not silently fixed here.

Generate 61-bit NTT primes close to `2^61`, input/output scale `2^40`, working scale `2^61`, degree 31, four double angles. Start with 17 primes. This is a budget candidate: actual arithmetic must certify the schedule and prime-ratio error.

- ModRaise: level 1 to 17, no division. Reinterpret logp 40 as 61 only as part of the explicitly normalized bootstrap representation, not as a value-preserving metadata edit.
- CoeffToSlot: plaintext diagonals at 61 bits; multiply, accumulate, then rescale once: level 17 to 16, working logp 61. Fold half factors into diagonals or account for their scale explicitly; do not use modular inverses to divide ciphertexts by 2.
- Polynomial: use balanced identities `T_(2k)=2*T_k^2-1` and `T_(2k+1)=2*T_k*T_(k+1)-T_1`. Cache terms, align branches by dropping levels, and charge up to five product layers for T_31, one affine-map level and one coefficient-weighting level: at most seven levels, ending at level 9. All coefficient terms must have a common scale before summation.
- Double angles: four multiply/relinearize/rescale layers: level 9 to 5, working logp 61.
- Restore factor `1/(2*pi)` in inverse-transform diagonals. Encode these diagonals, including `bottom_modulus/input_scale`, at 40 bits: product logp 101, rescale by a 61-bit prime gives logp 40, level 5 to 4. The large diagonal coefficients require checking raw-product headroom, not just logp.
- Output at level 4 supports another square: logp 80 before rescale, then 19 at level 3. Exact multiplication by `2^21` with a corresponding logp lift preserves decoded values and restores logp 40. Drop to level 1 before the next refresh.

Primes are not exact powers of two. Correct `2^61/bottom_modulus` explicitly and account for every actual rescale divisor. Current integer logp metadata only models power-of-two scales; certify the discrepancy against the error budget. If corrections cannot meet it, scale representation becomes a separately reviewed prerequisite. Do not claim the candidate works from this ledger alone.

## Reproducible numerical evidence

Run `python3 specs/research/bootstrap_experiment.py` (Python standard library only). It fits Chebyshev coefficients with cosine-node interpolation and evaluates via plaintext Clenshaw; sweeps degrees 15/31/63 and double-angle counts 3/4/5. It samples every permitted integer wrap with 101 small-message values across [-1,1], at ratio `2^-21`. This is a sampled floating-point experiment, not a uniform approximation proof or an encrypted simulation.

- Degree31/double4: maximum restored scalar error `6.86268e-8` at N=16 and `1.5451e-7` at N=64.
- Multiplying by N as a conservative reconstruction allowance gives sampled bounds `1.09803e-6` and `9.88861e-6`. Confirm the actual transform norm independently.
- Degree15/double3 at N=16 gives scalar error about 0.848; degree31/double3 at N=64 gives about 0.0332. Extra depth and degree cannot be chosen independently of the integer range.
- Reserve `1e-4` of slot error for periodic approximation/float effects and `9e-4` for incoming error, quantization, transforms and key switching. Validate their sum end to end; these allocations are targets, not measured guarantees.

## Implementation gates

Order tasks 0026-0031: signed automorphisms/public key setup; ModRaise; dense transforms; shallow polynomial evaluation; EvalMod; integration. Polynomial work is independent of transforms; integration needs both. No dependence on 0017/0022. Each task must run Rust formatting, linting, tests and example builds before handoff.

Validate exact centered lifting, context identity and error paths; asymmetric complex rotation/conjugation vectors; independent dense matrix round trips; Chebyshev identities and measured level cost; scalar and encrypted EvalMod at wrap/range boundaries. Integration checks fixed seeds, maximum slot error <=1e-3 versus original input, error versus pre-refresh decrypt, multiplication after refresh and at least two refresh cycles. Diagnostics may decrypt, evaluation may not. Reject invalid observable configuration; document out-of-domain behavior without promising encrypted detection.

## Implemented profile

The ledger below describes the original implementation. Task0036 optimizes it without changing the18-prime/44-bit profile: forward diagonals include division by2*radius, real/imaginary splitting uses addition and exact multiplication by-i, and EvalMod adds only the normalized quarter-period offset. The resulting ledger is18 -> forward17 -> balanced terms12 -> coefficient sum11 -> four double angles7 -> inverse6. Conjugation keys remain at17; inverse rotation keys/plaintexts move to7. Stage references and the full complex seeded regressions validate the fused representation. Bootstrap now spends12 primes rather than14.

Recurrence-specific headroom permits a91-bit weighted raw product at level2 (modulus near2^122), followed by rescale to logp30 and exact normalization to44 at level1. This does not relax the separate61-bit bootstrap multiplication threshold. The32-round18-prime recurrence now refreshes three times, after rounds17/22/27, instead of eight. At N=64 the demo worst checkpoint error was about6.67e-5; the1e-3 target and nonzero noise remain unchanged.

The implementation in `src/ckks/algorithms/bootstrap.rs` pays one extra level to split coefficient halves at working logp61. The supported 18-prime ledger is: raise to18, forward transform17, complex split16, affine map15, balanced terms10, coefficient sum9, double angles5, inverse transform4. Final diagonals encoded at44 bits return logp44. A subsequent square rescales to27 bits; `mul_const(1.0,17)` restores logp44 exactly without changing decoded values.

Logp40 missed the second-refresh 1e-3 gate at N=64 (about1.6e-3). Raising input/output precision to44 preserves a message ratio near2^-17 and reduces key-switch error. The release demo measures second-cycle errors about1.46e-5 (N=16) and9.46e-5 (N=64). Noise stays enabled. Plan construction rejects divisors farther than2^20 from2^61, bounding each power-of-two rescale discrepancy below roughly4.55e-13. Final multi-seed tests validate the accumulated effect rather than treating the scalar approximation as an encrypted proof.

Public setup binds keys and plaintexts to the plan's original context chain. `mod_raise`, `same_context` and `rescale_divisor` extend the ring contract. Only N=16/64, degree31/four double angles and logp44/61 are accepted; wider settings require new validation. The original scale ledger's lower-bound allowance is not evidence for additional supported profiles.

## Generic implementation reference review (2026-11-02)

Operator now requests generic bootstrapping after the structural pass. Findings from upstream main branches, inspected on this date:

- **Lattigo configuration and schedule:** [parameters_literal.go](https://github.com/tuneinsight/lattigo/blob/main/circuits/ckks/bootstrapping/parameters_literal.go) exposes transform factorization/scales, modular-reduction settings, ring size and secret distributions. [parameters.go](https://github.com/tuneinsight/lattigo/blob/main/circuits/ckks/bootstrapping/parameters.go) validates configuration and builds the circuit parameters from residual parameters. Stage levels are composed from stage depths; evaluator construction checks their consistency. Borrow this separation of user settings from derived stage levels, rather than copying its complete API.
- **Lattigo approximation:** [mod1_parameters.go](https://github.com/tuneinsight/lattigo/blob/main/circuits/ckks/mod1/mod1_parameters.go) explicitly models approximation range K, message ratio, degree, double angles and optional inverse approximation. Depth is derived from the selected polynomial construction and steps. Our conservative wrap bound must remain linked to the permitted secret weight; removing the ring-degree whitelist alone is insufficient.
- **Scale scheduling:** [mod1_evaluator.go](https://github.com/tuneinsight/lattigo/blob/main/circuits/ckks/mod1/mod1_evaluator.go) computes the polynomial target scale backwards through double angles using the actual consumed primes. [scale.go](https://github.com/tuneinsight/lattigo/blob/main/core/rlwe/scale.go) tracks scale with a 128-bit-precision real value. Our integer logp and bit-count rescaling cannot directly reproduce this generality. Recommendation: parameterize the existing near-power-of-two algorithm first, with explicit modulus constraints; arbitrary mixed-width scheduling needs a separately scoped scale design. This is a recommendation, not operator approval of that scope.
- **OpenFHE:** [ckksrns-fhe.cpp](https://github.com/openfheorg/openfhe-development/blob/main/src/pke/lib/scheme/ckksrns/ckksrns-fhe.cpp), especially EvalBootstrapSetup, GetBootstrapDepth, GetModDepthInternal and AdjustCiphertext, derives transform placement/depth from budgets and modular-reduction depth. It uses tuned polynomial tables and angle counts by secret distribution, and empirical default correction factors. Defaults and supported constraints remain legitimate; independent copies of stage offsets are the problem. Its helper notes additional normalization/correction levels outside the reported bootstrap depth, so do not transplant its depth accounting into our fused pipeline.
- **Keys:** Lattigo evaluator checkKeys validates required relinearization/Galois keys. Our BV keys are level-bound, so retain per-level key maps but derive their indices from our circuit schedule; copying upstream key storage is inappropriate.

Implementation direction to specify next: use every configured numerical field in construction and evaluation, derive the balanced-recurrence depth and all stage/key levels, validate headroom with checked arithmetic, and test multiple degrees/scales/degrees-of-approximation/angle counts/chain lengths. Keep the existing preset regression and add noisy encrypted validation for newly demonstrated configurations. Scalar approximation sweeps alone do not establish encrypted accuracy; these references do not validate arbitrary parameter combinations.

Operator clarification: build a general CKKS library. "Toy" in the repository name describes experimental status, not a restricted algorithm or API. The earlier staged near-power-of-two recommendation above is not an agreed scope limit. Scale-representation limitations are implementation gaps to address where general scheduling requires it. Replace the degree/profile whitelist with justified mathematical/backend validation; retain existing profiles as regression cases and name convenience presets descriptively. Generality does not imply every parameter combination satisfies an accuracy target.
