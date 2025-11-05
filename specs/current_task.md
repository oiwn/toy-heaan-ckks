# Current Tasks & Specifications

## Task 1: RNS-based Modulus Switching Implementation

Implement RNS-based modulus switching (ModDrop) for `RnsPolyRing` backend and create a demo example similar to `examples/naive_mod_switch.rs` that demonstrates key differences between naive and RNS approaches.

**References:**
- FHE Textbook: https://fhetextbook.github.io/RNSbasedModDroptextsfModDroptextsubscriptRNS.html
- Existing naive implementation: `examples/naive_mod_switch.rs`
- RNS backend: `src/rings/backends/rns/poly.rs`

---

### Background: Modulus Switching in CKKS

Modulus switching transforms a ciphertext from an old modulus `q` to a smaller new modulus `q'` while preserving decryptability with same secret key.

**Before switching (mod q):**
```
c0 + c1·s = Delta·m + e (mod q)
```

**After switching (mod q'):**
```
c0' + c1'·s = Delta'·m + e' + epsilon (mod q')
```

**Key Differences: Naive vs RNS**

| Aspect | Naive ModSwitch | RNS ModDrop |
|--------|----------------|-------------|
| **Formula** | `c' = round(c · ratio)` | Simply drop channels |
| **Noise added** | Bounded rounding error ≤ 1.5 | **Zero** (when q' divides q) |
| **Arithmetic** | Floating-point | Integer only |
| **Scale change** | Yes: `Delta' = Delta · (q'/q)` | No (RNS repr unchanged) |
| **Moduli constraint** | Any pair works | q' must divide q |

### Implementation Plan

**File:** `src/rings/backends/rns/poly.rs`

**Add PolyModSwitch Import:**
```rust
use crate::{PolyModSwitch, PolyRing, PolySampler, math::is_prime};
```

**Implement PolyModSwitch for RnsPolyRing:**
```rust
impl<const DEGREE: usize> PolyModSwitch<DEGREE> for RnsPolyRing<DEGREE> {
    fn mod_switch(&self, new_context: &Self::Context) -> Self {
        let current_count = self.channels();
        let new_count = new_context.channel_count();

        // Validate: new basis must have fewer primes
        assert!(new_count <= current_count, "ModDrop requires new basis <= current");

        // Validate: new basis primes must match prefix of current basis
        for (i, &new_prime) in new_context.primes().iter().enumerate() {
            let current_prime = self.basis.primes()[i];
            assert_eq!(new_prime, current_prime, "Bases must match prefix");
        }

        // RNS ModDrop: simply truncate to first new_count channels
        let new_coefficients = self.coefficients[..new_count].to_vec();

        Self {
            coefficients: new_coefficients,
            basis: new_context.clone(),
        }
    }
}
```

### Demo Example

**New file:** `examples/rns_mod_switch.rs`

**Parameters:**
```rust
const DEGREE: usize = 8;
const SCALE_BITS: u32 = 20;  // Delta = 2^20
const PRIME_BITS: [usize; 3] = [51, 41, 31];
```

**Demo Flow:**
1. Setup RNS modulus chain (3 levels)
2. Encrypt at highest level (Q2)
3. Decrypt at Q2 (baseline)
4. ModDrop: Q2 → Q1 (drop q3)
5. Decrypt at Q1
6. ModDrop: Q1 → Q0 (drop q2)
7. Decrypt at Q0
8. Show zero additional noise property

---

## Task 2: Ciphertext Multiplication for Naive Backend

Implement complete homomorphic multiplication with relinearization and rescaling for the naive u64 backend.

### Fixed Parameters

```rust
const DEGREE: usize = 8;
const SCALE_BITS: u32 = 20;          // Delta = 2^20
const Q_L: u64 = (1u64<<61) - 1;     // top modulus (level 1)
const Q_0: u64 = (1u64<<41) - 9;     // bottom modulus (level 0)
```

### Pipeline (Textbook Order)

1. **Basic multiplication (mod q_l):**
```
D0 = B1 * B2
D1 = B2 * A1 + B1 * A2
D2 = A1 * A2
```

2. **Relinearization (mod q_l):**
- Build `ctα = (D1, D0)`
- Build `ctβ from D2` using decomposition + eval key
- Output: `ctα+β ≈ RLWE_S(Delta^2 * M1*M2) mod q_l`

3. **Rescale (true modulus switch):**
```
A_hat = round(A / w_l) mod q_{l-1}
B_hat = round(B / w_l) mod q_{l-1}
```
Where `w_l = floor(Q_L / Q_0) ≈ 2^20 = Delta`

### RelinKey Structure

```rust
struct RelinKey {
    base_log2: u32,       // b, e.g. 16
    digits: u32,          // l, e.g. 4 for ~61-bit q1
    evk: Vec<(Poly, Poly)>, // length = l
    modulus_level: u8,    // 1 (matches q1)
}
```

### Key Generation

- Base: `beta = 2^b` (recommend b=16)
- Digits: `l = ceil(bitlen(q1)/b)` (for q1 ≈ 61b, l=4)
- For i in 1..l: `g[i] = Q / beta^i mod q1`
- Construct: `evk[i] = RLWE_S,σ( S^2 * g[i] )` at modulus q1

### Decomposition

Centered base-β digits: `D2 ≈ sum_{i=1..l} d2_i * (Q / beta^i) (mod q1)`
- Digits in range `[-β/2, +β/2)`
- Tie-break: halves toward +infinity

### Rescale Implementation (per coeff)

```rust
fn rescale_coeff_centered(x: u64, q1: u64, q0: u64, w: u128) -> u64 {
    let half_q1 = q1 >> 1;
    let xc_i128: i128 = if x > half_q1 { 
        (x as i128) - (q1 as i128) 
    } else { 
        x as i128 
    };

    let t_i128 = if xc_i128 >= 0 {
        (((xc_i128 as u128) + (w >> 1)) / w) as i128
    } else {
        -(((((-xc_i128) as u128) + (w >> 1)) / w) as i128)
    };

    let q0_i128 = q0 as i128;
    let xhat = ((t_i128 % q0_i128) + q0_i128) % q0_i128;
    xhat as u64
}
```

### API Surface

```rust
fn mul_once(ct1, ct2, &relin_key) -> Ciphertext {
  // requires ct1.level==ct2.level==1
  // returns level 0, scale_bits=SCALE_BITS
}
```

### Correctness Tests

- Vector test (4 slots): encode v1, v2; encrypt; mul_once; decrypt; decode
- Compare to v1 ⊙ v2 with tolerance eps = 1e-3
- Print level, q, scale after each stage (mul, relin, rescale)
- Reject calling mul_once when level==0

### Safety Checks

- With `Delta=2^20 and q1≈2^61`, one multiply is safe for small inputs (|m| ~ 1..4)
- Enforce: `q1 / q0 ≈ 2^20` so w_l≈Delta
- Use u128 temporaries for safe rescale division

### Call Site Example

```rust
let ct_alpha = (D1.clone(), D0.clone());
let ct_beta  = relinearize_d2(&D2, &rk);
let ct_mul   = add_ct(ct_alpha, ct_beta);   // scale ~ Delta^2 @ q1
let ct_out   = true_rescale(ct_mul, q1, q0); // scale ~ Delta @ q0
```

---

## Implementation Priority

1. **Complete RNS ModDrop first** - simpler, validates trait design
2. **Implement naive ciphertext multiplication** - more complex, builds core CKKS functionality
3. **Test both implementations** - ensure zero-noise property for RNS, correctness for naive

## Cross-References

- Architecture Overview: `specs/overview.md`
- Agent Guide: `AGENTS.md`
- Existing examples: `examples/naive_mod_switch.rs`, `examples/ckfs_naive.rs`
- Trait definitions: `src/rings/traits.rs`