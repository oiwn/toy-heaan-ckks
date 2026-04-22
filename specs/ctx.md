# Current task context

  Temperature sensor, 16 readings per minute (one per ~4 seconds). Each package is one ciphertext with 16 slots. We compute std-dev of that minute's readings, fully encrypted.

  Input normalization (decided): values pre-scaled to [0, 1] before encrypting.
  This keeps variance ∈ [0, 0.25] and makes the sqrt approximation tractable.

  package = [t0, t1, ..., t15]   (normalized, range [0, 1])

  Vector semantics: all CKKS operations on ct are element-wise (Hadamard).
  ct is a 16-slot vector; mean_ct replicates the same scalar in every slot.
  (ct - mean_ct)² = [(t0-μ)², (t1-μ)², ..., (t15-μ)²] — exactly the per-sample squared deviations.

  Pipeline

  ct ─── sum_slots ──→ sum_ct  (same value in every slot)
  │                       │
  │       mean_ct = sum / 16   [mul_plain(1/16) + rescale]
  │                       │
  │      ct - mean_ct ←──┘
  │           │
  │        sq_diff = (ct - mean)²  [mul_self + rescale]
  │           │
  │        sum_slots(sq_diff)
  │           │
  │        variance = sum / 16     [mul_plain(1/16) + rescale]
  │           │
  │        std_dev = poly_sqrt(variance)  [2 muls + rescales]
  └──────────────────────────────────────────────────────→ result

  Multiplicative depth

  ┌───────────────────────────┬─────────────────────┬─────────────────┐
  │           step            │         op          │ levels consumed │
  ├───────────────────────────┼─────────────────────┼─────────────────┤
  │ mean = sum / 16           │ mul_plain + rescale │ 1               │
  ├───────────────────────────┼─────────────────────┼─────────────────┤
  │ squaring (ct - mean)²     │ mul + rescale       │ 1               │
  ├───────────────────────────┼─────────────────────┼─────────────────┤
  │ variance = sum / 16       │ mul_plain + rescale │ 1               │
  ├───────────────────────────┼─────────────────────┼─────────────────┤
  │ poly_sqrt degree-3 Horner │ 2× mul + rescale    │ 2               │
  ├───────────────────────────┼─────────────────────┼─────────────────┤
  │ rotations in sum_slots    │ level-free          │ 0               │
  └───────────────────────────┴─────────────────────┴─────────────────┘

  Total: 5 prime levels. With 7 primes (5 + 2 headroom) at SCALE=30 / 30-bit primes.

  What needs to be built

  1. neg_ciphertext(ct) and sub_ciphertexts(ct1, ct2)
  Trivial — negate both components, then add. No key needed, no level cost.

  2. mul_plain_scalar(ct, scalar: f64)  [consumes 1 level]
  In CKKS all arithmetic lives in the scaled integer domain (coefficients ≈ Δ × message).
  A scalar α cannot be applied as a bare float multiply. The correct steps are:
    a. Encode α as a constant plaintext polynomial: all coefficients = round(α × Δ)
       (same encoding step as CkksEncoder, but for a uniform scalar)
    b. Multiply each ciphertext component by that polynomial (NTT-domain, component-wise)
       Result scale is Δ² × m × α  →  logp doubles
    c. Rescale: drop one RNS prime, divide by it  →  logp returns to Δ
  No relinearization needed (result stays degree-1).
  Reference: fhetextbook.github.io/CiphertexttoPlaintextMultiplication1.html

  3. sum_slots(ct, rotk_map: &HashMap<i32, RnsGadgetRotationKey>)  [level-free]
  Binary tree: log₂(16) = 4 steps.
  Needs rotation keys for offsets 1, 2, 4, 8.
  step 1: ct += rotate(ct, 1)   → slots [0+1, 1+2, ..., 15+0]
  step 2: ct += rotate(ct, 2)
  step 3: ct += rotate(ct, 4)
  step 4: ct += rotate(ct, 8)
  → every slot now contains the sum of all 16.
  Key management: HashMap<i32, RnsGadgetRotationKey> keyed by the offset value.

  4. eval_poly_horner(ct, coeffs: &[f64], rlk)  [depth = coeffs.len() - 1 levels]
  Generic Horner evaluator: a0 + ct*(a1 + ct*(a2 + ct*a3))
  Each inner step is: mul_ciphertexts_gadget + rescale + mul_plain_scalar (for the coefficient).
  Lives on CkksEngine; coefficients passed from the example.

  5. Polynomial sqrt approximation  [2 levels, uses eval_poly_horner]
  Approach: direct degree-3 Chebyshev fit to sqrt(x) on [0, 0.25].
  - Fits the level budget (2 muls)
  - Accurate for variance bounded away from 0 (realistic sensor data)
  - Singularity of sqrt'(x) at x=0 causes error for near-zero variance; acceptable for demo

  The standard production approach (Panda 2022, eprint.iacr.org/2022/423) uses
  Newton-Raphson for 1/sqrt(x):
      y_{i+1} = 0.5 · y_i · (3 - x · y_i²)   (depth 2 per step)
      sqrt(x) = x · y_final                    (1 more level)
  This is more general and robust but costs 3–4+ levels total.
  For this toy demo, direct Chebyshev on [0, 0.25] is the pragmatic choice.

  Chebyshev coefficients for sqrt on [0, 0.25] (to be fitted offline):
  p(x) ≈ a0 + a1·x + a2·x² + a3·x³   (fit with numpy/scipy before coding)

  6. The demo  (examples/std_dev.rs)
  - N=32, SLOTS=16, SCALE=30, 7 × 30-bit NTT-friendly primes
  - Synthesize 16 sensor readings: sine wave + noise, normalized to [0, 1]
  - Encrypt as one ciphertext
  - Run the full pipeline homomorphically (mean → sq_diff → variance → sqrt)
  - Decrypt, compare to plaintext std-dev
  - Report per-slot errors and final std-dev error
