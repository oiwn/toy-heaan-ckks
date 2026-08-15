# Current task: demo pre-requisites (Phase A)

Make the library ready for the browser std-dev demo. Demo architecture, UI,
and phases B–D live in `specs/demo.md` — this file is only the working
context for the pre-requisite tasks below.

## Target pipeline

The end state of this phase is a library-level `encrypted_std_dev`:

```text
ct = Enc([x0, x1, ..., x_{SLOTS-1}])
  -> sum_ct      = sum_slots(ct)                    -- rotation tree
  -> mean_ct     = sum_ct * (1/SLOTS)                -- mul_plain + rescale
  -> centered_ct = ct - mean_ct                      -- sub (scale-aligned)
  -> sq_diff_ct  = centered_ct * centered_ct          -- mul + rescale
  -> var_sum_ct  = sum_slots(sq_diff_ct)             -- rotation tree
  -> variance_ct = var_sum_ct * (1/SLOTS)            -- mul_plain + rescale
  -> v_norm_ct   = variance_ct * (1/B)               -- mul_plain + rescale
  -> std_dev_ct  = eval_chebyshev(v_norm_ct) * sqrt(B) -- poly eval
```

## Constants

Demo parameters (acceptance example; rationale in demo.md):

```text
N = 8192 (4096 slots, 32 used: 24h of readings, one per 45 min)
primes = 12 x 61-bit, scale = 2^61
DIGIT_BITS = 20 (4 digits/prime)
```

Toy parameters (unit tests): N = 16, 5 primes, as today.

Sqrt approximation (decided 2026-07-04):

```text
Chebyshev fit of sqrt, degree 7, interval [0.01, 1.0], B = 100
std envelope [1, 10] degC; error ~0.107 degC worst (std=1), ~0.01 typical
below interval: graceful high-biased extrapolation; v_norm > 1: explodes
-> envelope is a client-side contract (clamp before encryption)
```

## Library status

Working, verified 2026-07-04 after the base-β gadget key-switch fix
(`DIGIT_BITS` decomposition, keys indexed `[channel][digit]`):
encrypt/decrypt, add/sub, mul + relin + rescale, rotation, CRT wide
reconstruction. Numbers: all tests pass (110 unit + 6 integration);
`horner_chain` max error 7.9e-10; `std_dev_8` variance error 2.5e-13 at
61-bit primes. Plaintext Chebyshev reference: `examples/sensor_std_plain.rs`.

## Pre-requisite tasks

1. **DIGIT_BITS 15 → 20** (`src/crypto/engine.rs`), then re-verify:
   `cargo test`, `horner_chain`, `std_dev_8`.

2. **Serialization.** Serde derives for `RnsPoly`, `Ciphertext`, `Plaintext`,
   secret/public/gadget keys. `RnsBasis` is not serialized wholesale: encode
   the prime list, rebuild tables (and the `Arc`) on deserialize — polys in a
   ciphertext must share one `Arc<RnsBasis>` after load (`Arc::ptr_eq` checks
   in ops). Postcard round-trip test; assert byte sizes ≈ demo.md table
   (ciphertext ~1.6 MB, gadget key ~75 MB at demo params).

3. **Pipeline ops into the library** (`src/crypto/pipeline.rs`), extracted
   from `examples/std_dev_8.rs`:
   - `sum_slots` (power-of-two rotation tree)
   - `mul_plain_scalar` (encode scalar, mul, rescale)
   - `align_level` (generalizes example-local `reduce_ct_to_basis` +
     `reduce_sk`: drop channels to a target basis)
   - wide f64 decode (BigUint mixed-radix path from example `decode_wide`;
     the `[i64; DEGREE]` plaintext API overflows at high levels)
   - `fit_sqrt_chebyshev` (Chebyshev-basis fitter, move from
     `examples/sensor_std_plain.rs`) so coefficients are computed, not
     hardcoded.

4. **Encrypted Chebyshev eval** — Clenshaw in **Chebyshev basis** (never
   convert to power basis: ill-conditioned for small interval lower bound):
   - domain map `t = (2v - (hi+lo))/(hi-lo)`: plaintext affine op
   - recurrence `b_k = 2t·b_{k+1} - b_{k+2} + c_k`: one ct×ct mul per step;
     `b_{k+2}` must be mod-switched down to match `2t·b_{k+1}` after rescale
   - ~7–8 levels for degree 7
   - test: decrypt/decode vs plaintext Clenshaw on the same inputs.

5. **`encrypted_std_dev`** composing tasks 3 + 4 per the target pipeline.

6. **Acceptance example** `examples/std_dev_demo_params.rs`: N=8192,
   12 × 61-bit primes, 32 diurnal readings (std ≈ 2–3 °C). Asserts:
   - std-dev error < 0.11 °C (expect ~1e-2 at std ≈ 2.8)
   - intermediate variance error < 1e-9 (expect ~1e-13)

7. **Cargo/wasm hygiene**: move `dhat`, `comfy-table`, `tracing-subscriber`
   from `[dependencies]` to dev-dependencies (grep lib usage first);
   getrandom wasm backend feature; gate:
   `cargo build --target wasm32-unknown-unknown` green.

## Key files

| Path | Role |
|------|------|
| `specs/demo.md` | Demo spec + phases B–D |
| `src/crypto/engine.rs` | `DIGIT_BITS`, gadget keys, `gadget_key_switch`, mul/rotate |
| `src/crypto/pipeline.rs` | NEW — pipeline ops (task 3), Chebyshev eval (4), `encrypted_std_dev` (5) |
| `src/rings/backends/rns_ntt/poly.rs` | `RnsPoly` — storage, NTT, arithmetic |
| `src/rings/backends/rns_ntt/basis.rs` | `RnsBasis` — primes, tables, CRT |
| `examples/std_dev_8.rs` | Source of the ops to extract (task 3) |
| `examples/sensor_std_plain.rs` | Plaintext reference + Chebyshev fitter |
| `Cargo.toml` | Task 7 |

## Notes

- `errors_variance` on `CkksParams` is actually a std_dev (all call sites
  pass it to `sample_gaussian(std_dev, ...)`). Rename when touching it.
- `src/crypto/operations.rs` is a stale duplicate of `engine.rs`, tagged
  FIXME for removal.
- `PolyRescale` / `PolyModSwitch` traits are declared but unimplemented;
  real rescaling is `RnsPoly::rescale()` / `CkksEngine::rescale_ciphertext()`.

Later (not blocking): repair `sqrt_chebyshev_plain.rs` (negative LOWER →
NaN), retarget `sensor_std_plain.rs` to the B=100 envelope, update
`scripts/fit_sqrt_chebyshev.sage` to emit Chebyshev-basis coefficients.
