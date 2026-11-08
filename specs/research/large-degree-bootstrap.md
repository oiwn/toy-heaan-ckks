# Research: bootstrap error at large ring degree

Status: open (started 2026-11-03; upstream review and first experiments done, roadmap awaits operator review). Goal: keep the bootstrap error below the 1e-3 gate at N up to 8192 without a degree-specific message ratio, or know why that is impossible here. Tool: `examples/bootstrap_large.rs` (`--plan-only` prints the schedule and key memory).

## Problem

With the reference message ratio 2^17 (message scale 2^44 on a 61-bit q0) the error passes at N=16/64 and fails above about N=128. At N=8192 bootstrap only passes with ratio 2^4-2^6, which is not a general setting: a smaller ratio means a larger message scale (2^(61 - ratio)), so less headroom for the products that follow, and more sine linearization error for messages with large coefficients.

## Measurements

Seed 0, one 61-bit prime per level, chain = depth + 2, degree 31, max error over all slots, `cycle0 / cycle1`.

- Error vs N at ratio 17, weight 8, cycle0: N=64 1.6e-4, 128 2.8e-4, 256 1.6e-3, 512 4.0e-3, 1024 1.2e-2. About N^1.5. The same with weight 32 and with budgets 1/1 or 9/9.
- Error vs ratio at N=1024, weight 8: ratio 17 1.2e-2, 14 1.5e-3, 11 1.8e-4, 8 2.3e-5. A factor of 8 for every 3 bits, so the error is proportional to q0/message_scale.
- Weight 64 (radius 34) needs five double angles with degree 31: N=1024 ratio 8 gives 2.9e-4 with four, 3.3e-5 with five (degree 63 with five: 3.1e-5).
- N=8192, weight 64, budgets 6/4, five double angles, 22 primes: plan 16 s, keys 50 s (10 relin, 82 rotation, 1 conjugation; ~7.5 GB), 40 s per bootstrap. Ratio 8: 2.2e-3 / 4.7e-3. Ratio 6: 5.6e-4 (cycle0).
- Not secure: N=8192 allows about 218 modulus bits at 128-bit security; these chains have about 1300.

Ratio sweep (2026-11-03, `--ratio`, input |x| = 0.8 per slot, 20-bit digits):
- N=64, weight 32, four double angles: ratio 17 1.2e-4 / 2.6e-4, 12 4.1e-6 / 1.1e-5, 10 1.2e-6 / 2.2e-6, 8 1.1e-5 / 1.9e-5, 6 1.7e-4 / 3.1e-4, 4 2.7e-3 / 4.9e-3, 2 4.3e-2 / 7.4e-2. Five double angles give the same at ratio 8/6/4.
- N=1024, weight 64, five double angles: ratio 17 1.7e-2 / 3.2e-2, 12 5.2e-4 / 9.1e-4, 10 1.3e-4 / 2.5e-4, 8 3.3e-5 / 6.9e-5, 6 1.1e-5 / 2.4e-5, 4 1.0e-4 / 2.7e-4.
- N=8192, weight 64, five double angles: ratio 12 3.5e-2 / 8.6e-2, 10 8.8e-3 / 2.1e-2, 8 2.2e-3 / 4.7e-3, 6 5.6e-4 / 1.0e-3, 4 2.0e-4 / 3.2e-4.
- The error is U-shaped in the ratio. Above the optimum it doubles per bit (noise amplification); below it, it grows 16x per 2 bits, the cubic sine linearization error. The optimum moves down with N: 10 at N=64, 6 at N=1024, 4 or lower at N=8192, because large-N message coefficients are small. No single ratio passes everywhere: ratio 8 fails at N=8192, ratio 4 fails at N=64. The best N=8192 result is 12 bits (ratio 4), against 26 bits Lattigo reports at N=2^16 with ratio 8.

Full sweep (2026-11-08, task 0045, after 0052; weight min(N/2, 64), five double angles above weight 32, default budgets, `cycle0 / cycle1`; rows N=64/1024/8192 above):

| N | 17 | 14 | 12 | 10 | 8 | 6 | 4 |
|---|---|---|---|---|---|---|---|
| 16 | 1.6e-5 / 3.3e-5 | 1.9e-6 / 3.9e-6 | 5.3e-7 / 1.2e-6 | 1.3e-6 / 3.0e-6 | 2.0e-5 / 4.5e-5 | 3.3e-4 / 7.2e-4 | 5.2e-3 / 1.1e-2 |
| 32 | 4.9e-5 / 1.0e-4 | 6.2e-6 / 1.5e-5 | 1.4e-6 / 2.8e-6 | 1.5e-6 / 3.0e-6 | 2.1e-5 / 4.4e-5 | 3.3e-4 / 7.0e-4 | 5.3e-3 / 1.1e-2 |
| 128 | 3.7e-4 / 7.0e-4 | 4.4e-5 / 9.8e-5 | 1.1e-5 / 2.3e-5 | 2.7e-6 / 7.3e-6 | 3.4e-6 / 7.4e-6 | 4.9e-5 / 1.0e-4 | 7.8e-4 / 1.7e-3 |
| 256 | 1.7e-3 / 4.1e-3 | 2.1e-4 / 5.0e-4 | 5.4e-5 / 1.3e-4 | 1.3e-5 / 3.4e-5 | 5.0e-6 / 1.1e-5 | 4.2e-5 / 7.8e-5 | 6.7e-4 / 1.2e-3 |
| 512 | 4.3e-3 / 8.9e-3 | 5.3e-4 / 1.1e-3 | 1.3e-4 / 2.9e-4 | 3.3e-5 / 7.0e-5 | 8.1e-6 / 1.6e-5 | 1.7e-5 / 3.2e-5 | 2.6e-4 / 5.2e-4 |
| 2048 | 3.9e-2 / 7.8e-2 | 4.8e-3 / 9.5e-3 | 1.2e-3 / 2.5e-3 | 3.1e-4 / 6.5e-4 | 7.5e-5 / 1.4e-4 | 2.2e-5 / 4.2e-5 | 1.0e-4 / 2.6e-4 |
| 4096 | 1.5e-1 / 3.3e-1 | 1.8e-2 / 3.8e-2 | 4.5e-3 / 8.2e-3 | 1.2e-3 / 2.3e-3 | 2.8e-4 / 5.6e-4 | 7.1e-5 / 1.9e-4 | 1.0e-4 / 2.2e-4 |

The optimum falls about one bit per doubling of N: 12 at N=16, 10 at 64-128, 8 at 256-512, 6 at 1024-4096, 4 at 8192. The reference rule (0045): ratio 17 up to N=64 (recurrence and the small demos rely on 2^44; worst 2.6e-4), above it `17 - log2 N`, at least 4. The odd points on that line also pass: N=256 ratio 9 7.1e-6 / 1.4e-5, N=1024 ratio 7 1.6e-5 / 2.7e-5, N=4096 ratio 5 4.6e-5 / 9.9e-5.

Digit width (temporary `DIGIT_BITS` change, reverted; weight 64, five double angles, cycle0):
- N=256: ratio 17 at 10/20/30 bits: 1.1e-4, 1.7e-3, 2.3. Ratio 8: 2.7e-6, 5.0e-6, 4.5e-3.
- N=1024: ratio 17: 3.9e-3, 1.7e-2, 1.2e1. Ratio 8: 7.7e-6, 3.3e-5, 2.3e-2.
- 20 to 30 bits multiplies the error by about 2^10 (x1300, x700), as the key-switch estimate predicts, so key switching dominates at 20 bits. 20 to 10 bits gains only 2-15x, not 2^10: below 20 bits a second source takes over. It is also amplified by the ratio (N=1024: x500 from ratio 8 to 17) and grows fast with N (ratio 17, 10 bits: x34 from N=256 to 1024). Its origin is open.

Transform budgets (N=1024, ratio 17 / 8, cycle0, max and mean):
- cts/stc 1/1 (one level, 511 rotations): max 1.1e-2 / 2.3e-5, mean 5.4e-4 / 1.1e-6.
- 5/4 (default): max 1.7e-2 / 3.3e-5, mean 4.9e-3 / 9.6e-6.
- 9/9 (two rotations per level): max 1.9e-2 / 3.7e-5, mean 7.2e-3 / 1.4e-5.
- The error does not follow the rotation count: the one-level transform has the most rotations and the smallest mean error, with a few outlier slots setting the max. A likely reason (not measured): noise added at a factored level passes through the gain of every later level, while a dense one-level transform averages it with diagonals of size ~1/N.

## Hypothesis

Bootstrap computes `m + q0 * I` in the slots and removes `q0 * I` with the sine. Any additive noise `e` present before EvalMod (in units of the coefficient) comes out as `e / message_scale` in the result, so it is amplified by q0/message_scale compared with ordinary operations. That matches the ratio data exactly.

At 20-bit digits the noise that grows with N is key switching (confirmed by the digit-width runs above). Our BV switch (`src/ckks/keys/key_switching.rs`) sums `gadget_len * N` products of a digit below 2^20 and a key error of std_dev ~1.8, so it adds about `sqrt(gadget_len * N) * 2^20 * 1.8` per switch. A CoeffsToSlots rotation acts on the raised ciphertext before the diagonal multiply and rescale, so this noise sits directly against the 2^61 working scale. At N=8192, level 22 (gadget_len 88): about 2^30.6, i.e. 2^-30 relative to the slot value `I + m * 2^-ratio`, so about 2^-13 in the message at ratio 17, per rotation, before the embedding growth. Hybrid key switching divides this term by a special modulus P (see the upstream comparison).

The second source seen at 10-bit digits is open. Candidates: rescale rounding in CtS, the f64 encoding of the CtS diagonals (`diagonal_scale`) at large N, ModRaise of the input noise, and the gain of later transform levels on earlier noise. It caps what hybrid key switching alone can deliver: at N=1024 ratio 17, 10-bit digits leave 3.9e-3.

## Upstream comparison (2026-11-03)

Read from the OpenFHE and Lattigo main branches (links in References). Neither runs with our ratio or our key switching.

Message ratio (q0 / message at ModRaise):
- Lattigo: `DefaultLogMessageRatio = 8`. `ScaleDown` multiplies the input by an integer and rescales so the message sits at q0/2^8 before `ModUp`; the rescale error is corrected at the end. The code notes the ratio "directly impacts the precision" and recommends the arcsine correction (`Mod1InvDegree` 5 or 7) for small ratios. The precision-optimized sets use ratio 2.
- OpenFHE: the correction factor `cf` is a fitted line, `round(-0.2419 * (2*log2 N + log2 slots) + 19.081)` clamped to [7, 14] (FLEXIBLEAUTO). Full packing: N=1024 gives 12, N=8192 gives 10, N=2^16 gives 8. `EvalBootstrap` scales the input by `2^-(cf - log2(q0/scale))` with `AdjustCiphertext`, raises, and multiplies back after EvalMod. The comment calls it "emulating a larger first modulus to improve precision".
- Ours: 17 (a 2^44 message on a 61-bit q0), chosen at N=64 when ratio 21 missed the gate. Upstream picks the ratio as the balance between noise amplification (grows with the ratio) and sine linearization error (shrinks as the ratio grows); both libraries land at 8-12.

Key switching:
- OpenFHE: `HYBRID` by default, `numLargeDigits` 3 in the bootstrap examples, special primes of 60 bits (`AUXMODSIZE`), enough of them to cover the largest digit plus a margin "to bring hybrid key switching error closer to the modulus switching error" (`rns-cryptoparameters.cpp`). Flow in `keyswitch-hybrid.cpp`: digit split, `ApproxSwitchCRTBasis` to QP, inner product, `ApproxModDown` by P.
- Lattigo: `GadgetProductLazy` computes the inner product mod QP, then `ModDownQPtoQNTT` divides by P; dnum = ceil(#Qi / #Pi). Default bootstrap sets use LogP 4-6 x 61 bits.
- Effect: after the division by P, a key switch adds about the rounding term `sqrt(N * h)` instead of the digit noise. Ours never divides, so every rotation adds the full digit noise.

Linear transforms:
- Both factor CtS/StC into a few levels (Lattigo CtS 4 levels on 56-bit primes, StC 3 on 39-bit) and use BSGS with hoisting: decompose once, rotate baby steps in QP, one ModDown per giant step (Lattigo `MultiplyByDiagMatrixBSGS`, OpenFHE `EvalFastRotationExt`). Fewer key switches and fewer roundings.

Other levers:
- Arcsine correction after the sine (Lattigo `Mod1InvDegree`): allows ratios down to 2.
- Meta-BTS (Lattigo `IterationsParameters`, OpenFHE `EvalBootstrap(ct, 2, precision)`): bootstrap the scaled error and subtract it; about twice the precision for twice the cost and one extra level.
- Sparse-secret encapsulation (Lattigo `EphemeralSecretWeight = 32`, OpenFHE `SPARSE_ENCAPSULATED`): small K for a dense user key. A security feature; it changes K, not the key-switch noise.
- OpenFHE `FLEXIBLEAUTO` scale tracking: we already track exact scales (0040).

Reported precision (upstream claims, not reproduced here): Lattigo N=2^16 default sets 26.6 bits (H=192, scale 2^40), 32.1 bits with ratio 2 and arcsine degree 7, 23.8-29.8 bits dense.

Reduces noise: hybrid key switching, a smaller ratio, arcsine correction, Meta-BTS, fewer key switches. Reduces cost: BSGS, hoisting, dnum choice, sparse encapsulation (smaller K).

## Experiments

- [x] Ratio sweep at N=64/1024/8192 (Measurements).
- [x] Digit width 10/20/30 at N=256/1024: key switching dominates at 20 bits; a second ratio-amplified source remains below it.
- [x] Transform budgets 1/1, 5/4, 9/9 at N=1024: error does not follow the rotation count.
- [ ] Per-stage error: decrypt after ModRaise, each CtS group, EvalMod and each StC group against a plaintext pipeline, at 10-bit digits, to find the second source.
- [ ] Plaintext precision: compare the CtS output with an f64 reference of the same matrix to rule out diagonal encoding loss.

## Proposed roadmap

Operator order (2026-11-03): ratio default first, then hybrid key switching. Operator agreed the roadmap 2026-11-03; drafts: 0045 ratio default (before 0044), 0046 stage noise, 0047 RNS basis extension, 0048 hybrid key switching, 0049 BSGS/hoisting, 0050 arcsine correction, 0051 Meta-BTS. Step 3 gets a task once 0046 names the source.

1. Done in 0045 (2026-11-08): `reference(N)` uses ratio 17 up to N=64, else `17 - log2 N` (at least 4), weight min(N/2, 64), five double angles above weight 32; the output must survive one squaring. Original note: the sweep shows no fixed ratio passes at all degrees (8 fails at N=8192, 4 fails at N=64), so the default must depend on N, like OpenFHE's fit but from our data: about 10 at N<=64, 6 at N=1024, 4 at N=8192 with 20-bit digits. A larger message scale (2^51-2^57) needs a headroom review of the output level and of recurrence (0044); the N=16/64 regression numbers change. Fold into 0044 or a task before it (operator decision). The fit must be redone after step 2.
2. Hybrid key switching: special primes P, a QP basis, basis-extension ModUp and ModDown in the RNS backend, a `PolyRing` extension and a new `KeySwitchingKey`. Likely two tasks: backend basis conversion, then key switching. Expected: removes the dominant 20-bit term (x4 at N=1024 ratio 8, x15 at N=256 ratio 17 by the 10-bit proxy) and cuts key memory. The per-stage experiment should run first, so the second source is known before this lands.
3. The second noise source, once the per-stage experiment names it.
4. BSGS with hoisted rotations in CtS/StC: fewer key switches and roundings, much less time and key memory at large N.
5. Arcsine correction, to let the ratio drop further where the sine floor limits it (N<=1024).
6. Meta-BTS iterations, after 2 and 3.

## References

- OpenFHE: [ckksrns-fhe.cpp](https://github.com/openfheorg/openfhe-development/blob/main/src/pke/lib/scheme/ckksrns/ckksrns-fhe.cpp) (`EvalBootstrapSetup` correction factor, `EvalBootstrap`, `AdjustCiphertext`), [ckksrns-fhe.h](https://github.com/openfheorg/openfhe-development/blob/main/src/pke/include/scheme/ckksrns/ckksrns-fhe.h) (K and degrees per secret), [keyswitch-hybrid.cpp](https://github.com/openfheorg/openfhe-development/blob/main/src/pke/lib/keyswitch/keyswitch-hybrid.cpp), [rns-cryptoparameters.cpp](https://github.com/openfheorg/openfhe-development/blob/main/src/pke/lib/schemerns/rns-cryptoparameters.cpp) (special primes), [examples](https://github.com/openfheorg/openfhe-development/tree/main/src/pke/examples) (simple/advanced/iterative bootstrapping).
- Lattigo: [default_parameters.go](https://github.com/tuneinsight/lattigo/blob/main/circuits/ckks/bootstrapping/default_parameters.go), [parameters_literal.go](https://github.com/tuneinsight/lattigo/blob/main/circuits/ckks/bootstrapping/parameters_literal.go), [evaluator.go](https://github.com/tuneinsight/lattigo/blob/main/circuits/ckks/bootstrapping/evaluator.go) (`ScaleDown`, `ModUp`, iterations), [evaluator_gadget_product.go](https://github.com/tuneinsight/lattigo/blob/main/core/rlwe/evaluator_gadget_product.go), [lintrans_evaluator.go](https://github.com/tuneinsight/lattigo/blob/main/circuits/common/lintrans/lintrans_evaluator.go) (BSGS, hoisted ModDown).
- Papers: eprint [2021/204](https://eprint.iacr.org/2021/204) (hybrid key switching), [2020/1118](https://eprint.iacr.org/2020/1118) (scaling techniques), [2020/1203](https://eprint.iacr.org/2020/1203) (Bossuat et al., non-sparse-key bootstrapping, double hoisting), [2020/552](https://eprint.iacr.org/2020/552) (arcsine / minimax), [2022/024](https://eprint.iacr.org/2022/024) (sparse-secret encapsulation), [2022/1167](https://eprint.iacr.org/2022/1167) (Meta-BTS).
