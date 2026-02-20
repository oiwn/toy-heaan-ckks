# CKKS Parameter Selection Guide

Comprehensive guide for selecting RNS-CKKS parameters based on production
library practices (SEAL, Lattigo, OpenFHE).

---

## Table of Contents

1. [Core Concepts](#core-concepts)
2. [Modulus Chain Structure](#modulus-chain-structure)
3. [Production Library Comparison](#production-library-comparison)
4. [Parameter Selection Rules](#parameter-selection-rules)
5. [Recommended Configurations](#recommended-configurations)
6. [Implementation Guidelines](#implementation-guidelines)
7. [References](#references)

---

## Core Concepts

### There is NO "Plaintext Modulus" in CKKS!

Unlike BFV/BGV schemes, CKKS does not have a separate plaintext modulus.
Instead:

- **Scale Factor**: Controls encoding precision for real/complex numbers
  - Formula: `scale_factor = 2^scale_bits`
  - Typical values: `2^40`, `2^45`, `2^50`
  - Determines fixed-point precision during encoding

- **Ciphertext Modulus (Q)**: Product of RNS primes, controls noise budget
  - Formula: `Q = q_0 * q_1 * q_2 * ... * q_L * q_special`
  - Total size: typically 200-600 bits
  - Each `q_i` is a prime around 30-60 bits

### Scale Factor vs Modulus Relationship

**Critical insight**: Intermediate RNS primes should match the scale factor.

```
If q_i ~= scale_factor, then rescaling is clean:
  ct' = ct / q_i  ~=  ct / scale_factor

If q_i != scale_factor, rescaling introduces error:
  ct' = ct * (q_i / scale_factor) + error
```

**Example precision loss:**
- Scale: `scale_factor = 2^40 = 1,099,511,627,776`
- Prime: `q_i = 1,099,511,627,689` (closest 40-bit NTT prime)
- Ratio error: `q_i / scale_factor ~= 0.999999999921` (~1e-10 per rescale)

The closer `q_i` is to a power of 2, the lower the rescaling error.

---

## Modulus Chain Structure

### Universal Pattern Across All Libraries

All production CKKS implementations follow this pattern:

```
+-----------------------------------------------+
| q_0   q_1   q_2   ...   q_L   q_special        |
| 60b  scale scale ...  scale  60b              |
| LARGE SCALE SCALE ... SCALE  LARGE            |
+-----------------------------------------------+
  '      '                     '        '
 first  intermediate            last   special
 prime   primes                 scale  prime
```

### Role of Each Prime

1. **First Prime (q_0)**: 50-60 bits
   - **Purpose**: Extra noise budget for fresh ciphertexts
   - **Why larger**: Provides highest precision when decrypting
   - **Not used in rescaling**: Preserved throughout computation

2. **Intermediate Primes (q_1 ... q_L)**: Match scale (30-55 bits)
   - **Purpose**: Support rescaling after multiplications
   - **Critical property**: All same size, as close to `2^scale_bits` as possible
   - **Dropped during rescaling**: One prime removed per multiplication

3. **Last/Special Prime (q_special)**: 50-60 bits
   - **Purpose**: Key switching and rotation operations
   - **Why larger**: Should be as large as the largest intermediate prime
   - **Never dropped**: Reserved for special operations

### Rescaling Mechanics

After homomorphic multiplication, the scale becomes `scale_factor^2`:

```
ct_mult = Enc(m_1 * m_2, scale = scale_factor^2)
```

Rescaling divides by one RNS prime to restore the scale:

```
Rescale: ct' = ct / q_i mod Q_remaining

If q_i ~= scale_factor:
  ct' ~= Enc(m_1 * m_2, scale = scale_factor)  (correct)
```

**Modulus chain evolution:**
```
Level L:   Q_L = q_0 * q_1 * q_2 * ... * q_L * q_special
            -> (multiply + rescale by q_L)
Level L-1: Q_{L-1} = q_0 * q_1 * q_2 * ... * q_{L-1} * q_special
            -> (multiply + rescale by q_{L-1})
Level L-2: Q_{L-2} = q_0 * q_1 * q_2 * ... * q_{L-2} * q_special
            ...
```

**Multiplicative depth**: Number of intermediate primes = number of
multiplications supported.

---

## Production Library Comparison

### Microsoft SEAL (C++)

**Typical configuration:**
```cpp
poly_modulus_degree = 8192
coeff_modulus = CoeffModulus::Create(8192, { 60, 40, 40, 60 })
scale = pow(2.0, 40)
```

**Modulus chain:**
```
q_0:       60 bits  - first prime
q_1:       40 bits  - intermediate (match scale)
q_2:       40 bits  - intermediate (match scale)
q_special: 60 bits  - special prime

Total Q ~= 200 bits
Multiplicative depth: 2
```

**Design rationale** (from SEAL examples):
"Choose a 60-bit prime as the first prime in coeff_modulus for highest
precision when decrypting. Choose another 60-bit prime as the last element,
as this will be used as the special prime and should be as large as the
largest of the other primes. Choose the intermediate primes to be close to
each other."

**Constraints:**
- Minimum prime size: 20 bits
- Maximum prime size: 60 bits (hardware limit for 64-bit arithmetic)
- All intermediate primes: same size
- First/last primes: at least 10 bits larger than intermediate

---

### Lattigo (Go)

**Typical configuration:**
```go
LogN = 14  // N = 16384
LogQ = 438 bits total
Modulus chain: 7 primes (30-60 bits each)
DefaultScale = 2^45
```

**Modulus chain:**
```
Example: LogN=14, LogQP=438

Q chain (7 primes):
  q_0:    60 bits
  q_1:    45 bits  - matches DefaultScale
  q_2:    45 bits
  q_3:    45 bits
  q_4:    45 bits
  q_5:    45 bits
  q_6:    60 bits

P chain (2 primes for key switching):
  p_0:    60 bits
  p_1:    60 bits

Total Q ~= 330 bits
Total P ~= 120 bits
Multiplicative depth: 5
```

**Design philosophy** (from Lattigo docs):
"The used moduli are chosen to be of size 30 to 60 bits for the best
performance. The individual size of each of the moduli also has an effect on
the error introduced during the rescaling, since they cannot be powers of 2,
so they should be chosen as NTT primes as close as possible to a power of 2."

**HE Standards-compliant parameters:**
- `{12, 109, 3.2}`: logN=12, logQ=109 (N=4096, Q ~= 128 bits)
- `{13, 218, 3.2}`: logN=13, logQ=218 (N=8192, Q ~= 218 bits)
- `{14, 438, 3.2}`: logN=14, logQ=438 (N=16384, Q ~= 438 bits)
- `{15, 881, 3.2}`: logN=15, logQ=881 (N=32768, Q ~= 881 bits)

---

### OpenFHE (C++)

**Typical configuration:**
```cpp
multDepth = 5
scaleModSize = 50       // bits per intermediate prime
firstModSize = 60       // first prime size
ringDim = 0             // auto-calculate for security
```

**Modulus chain:**
```
Example: depth=5, scaleModSize=50, firstModSize=60

Q = q_0 * q_1 * q_2 * q_3 * q_4 * q_5 * q_special

q_0:       60 bits  - firstModSize (extra precision)
q_1:       50 bits  - scaleModSize
q_2:       50 bits  - scaleModSize
q_3:       50 bits  - scaleModSize
q_4:       50 bits  - scaleModSize
q_5:       50 bits  - scaleModSize
q_special: 60 bits  - auto-added

Total Q ~= 370 bits
Multiplicative depth: 5
```

**Default values** (if not specified):
- `firstModSize = 60`
- `scalingModSize = 50`
- Ring dimension: auto-calculated for 128-bit security

**Constraint:**
- `scalingModSize < 60` when using 64-bit native integers

---

### Comparison Table

| Library | First Prime | Intermediate Primes | Special Prime | Typical Scale | Max Prime Size |
|---------|-------------|---------------------|---------------|---------------|----------------|
| SEAL    | 60 bits     | 40 bits             | 60 bits       | 2^40          | 60 bits        |
| Lattigo | 60 bits     | 45 bits             | 60 bits       | 2^45          | 60 bits        |
| OpenFHE | 60 bits     | 50 bits             | 60 bits       | 2^50          | 60 bits        |

**Universal findings:**
- First prime: always 50-60 bits (extra precision)
- Intermediate primes: match scale factor (30-55 bits)
- Special prime: always 50-60 bits (key switching)
- Max prime size: 60 bits (hardware constraint)
- All libraries use NTT-friendly primes close to `2^k`

---

## Parameter Selection Rules

### Rule 1: First Prime (q_0)

**Size**: 50-60 bits (typically 60)

**Properties:**
- Largest prime in the chain (or equal to special prime)
- Provides maximum noise budget for fresh encryptions
- Never dropped during rescaling

**Selection:**
```rust
let first_prime_bits = 60;
let q_0 = find_closest_ntt_prime(DEGREE, first_prime_bits);
```

---

### Rule 2: Intermediate Primes (q_1 ... q_L)

**Size**: Match scale factor exactly.

**Properties:**
- All must be the same size
- Should be as close to `2^scale_bits` as possible
- Each supports one rescaling operation
- Count determines multiplicative depth

**Selection:**
```rust
let scale_bits = 40;
let mult_depth = 5;

let intermediate_primes: Vec<u64> = (0..mult_depth)
    .map(|_| find_closest_ntt_prime(DEGREE, scale_bits))
    .collect();
```

**Why same size?**
- Maintains consistent precision across levels
- Simplifies error analysis
- Standard practice in all libraries

---

### Rule 3: Special Prime (q_special)

**Size**: 50-60 bits (typically 60)

**Properties:**
- Used for key switching and rotations
- Should be as large as the largest intermediate prime
- Never dropped during computation
- Required even if you do not use rotations

**Selection:**
```rust
let special_prime_bits = 60;
let q_special = find_closest_ntt_prime(DEGREE, special_prime_bits);
```

---

### Rule 4: NTT-Friendly Primes

**Requirement**: All primes must satisfy `p == 1 (mod 2N)`

**Why?**
- Ensures primitive 2N-th root of unity exists modulo p
- Required for Number Theoretic Transform (NTT)
- Enables fast polynomial multiplication in RNS representation

**Selection algorithm:**
```rust
fn find_closest_ntt_prime(degree: usize, target_bits: u32) -> u64 {
    let target = 1u64 << target_bits;
    let modulus = (2 * degree) as u64;

    // Search for prime p where:
    // 1. p == 1 (mod 2N)
    // 2. |p - 2^target_bits| is minimized
    // 3. p is actually prime

    for offset in 0..10_000_000 {
        for candidate in [target + offset, target - offset] {
            if candidate % modulus == 1 && is_prime(candidate) {
                return candidate;
            }
        }
    }
    panic!("Could not find suitable NTT prime near 2^{}", target_bits);
}
```

**Tip**: Pre-compute and cache prime tables for common configurations.

---

### Rule 5: Security Constraints

**128-bit security targets** (from HE standards):

| Ring Dimension (N) | Max Total Modulus (log2 Q) |
|--------------------|----------------------------|
| 2048               | ~54 bits                   |
| 4096               | ~109 bits                  |
| 8192               | ~218 bits                  |
| 16384              | ~438 bits                  |
| 32768              | ~881 bits                  |

**Constraint check:**
```rust
fn validate_security(ring_dim: usize, total_modulus_bits: usize) -> bool {
    let max_bits = match ring_dim {
        2048 => 54,
        4096 => 109,
        8192 => 218,
        16384 => 438,
        32768 => 881,
        _ => panic!("Unsupported ring dimension"),
    };

    total_modulus_bits <= max_bits
}
```

If security check fails: increase N or decrease total Q (reduce depth or scale).

---

### Rule 6: Precision vs Depth Tradeoff

**Total modulus budget (bits)**: `Q_total ~= q_0 + L * scale_bits + q_special`

**Example calculation:**
```
Target: 128-bit security with N=8192 -> max Q ~= 218 bits

Available budget: 218 bits
- First prime: 60 bits
- Special prime: 60 bits
- Remaining: 218 - 120 = 98 bits

For scale = 2^40:
  Max depth L = 98 / 40 = 2 multiplications

For scale = 2^30:
  Max depth L = 98 / 30 = 3 multiplications
```

**Tradeoff:**
- Higher scale -> better precision, lower depth
- Lower scale -> worse precision, higher depth

**Rule of thumb**: Choose scale based on application precision requirements,
then calculate maximum depth.

---

## Recommended Configurations

### Configuration 1: Educational/Toy (Fast)

**Use case**: Learning, testing, small examples

```rust
const DEGREE: usize = 512;
const SCALE_BITS: u32 = 40;
const MULT_DEPTH: usize = 3;

// Modulus chain: [55, 40, 40, 40, 55]
let primes = vec![
    // First prime: ~55 bits
    35184372088833,   // 2^45 + 1, NTT-friendly for N=512

    // Intermediate primes: ~40 bits (match scale)
    1099511627777,    // Close to 2^40, NTT-friendly
    1099511627777,    // Reuse same prime (simplified)
    1099511627777,

    // Special prime: ~55 bits
    35184372088833,   // Same as first
];

// Properties:
// - Fast operations (small N)
// - Sufficient precision (40 bits ~= 10 decimal digits)
// - Supports 2-3 multiplications
// - NOT production-secure (too small N)
```

**Precision**: ~10-12 decimal digits
**Performance**: fast (milliseconds)
**Security**: educational only (N too small)

---

### Configuration 2: Balanced (Recommended for Development)

**Use case**: Development, prototyping, medium-scale applications

```rust
const DEGREE: usize = 8192;
const SCALE_BITS: u32 = 40;
const MULT_DEPTH: usize = 5;

// Modulus chain: [60, 40, 40, 40, 40, 40, 60]
let primes = vec![
    // First prime: 60 bits
    1152921504606584833,   // Close to 2^60, NTT-friendly for N=8192

    // Intermediate primes: 40 bits each
    1099511627777,
    1099511644161,
    1099511660545,
    1099511676929,
    1099511693313,

    // Special prime: 60 bits
    1152921504606584837,
];

// Total Q: 60 + 5*40 + 60 = 320 bits (exceeds 218-bit limit for N=8192)
// ISSUE: Reduce depth or increase N
```

**Fix for 128-bit security:**
```rust
// Option A: Reduce depth
const MULT_DEPTH: usize = 3;  // Total: 60 + 3*40 + 60 = 240 bits (still too much)

// Option B: Reduce scale
const SCALE_BITS: u32 = 30;   // Total: 60 + 5*30 + 60 = 270 bits (still too much)

// Option C: Increase N
const DEGREE: usize = 16384;  // Max Q ~= 438 bits (fits)
```

**Recommended balanced config:**
```rust
const DEGREE: usize = 16384;
const SCALE_BITS: u32 = 40;
const MULT_DEPTH: usize = 5;

// Total Q: 320 bits < 438 bits limit
```

**Precision**: ~10-12 decimal digits
**Performance**: medium (hundreds of ms)
**Security**: 128-bit secure

---

### Configuration 3: High Precision

**Use case**: Applications requiring high numerical precision

```rust
const DEGREE: usize = 16384;
const SCALE_BITS: u32 = 50;
const MULT_DEPTH: usize = 4;

// Modulus chain: [60, 50, 50, 50, 50, 60]
let primes = vec![
    // First prime: 60 bits
    1152921504606584833,

    // Intermediate primes: 50 bits each (match scale)
    1125899906842597,
    1125899906842601,
    1125899906842613,
    1125899906842619,

    // Special prime: 60 bits
    1152921504606584837,
];

// Total Q: 60 + 4*50 + 60 = 320 bits < 438 bits limit
```

**Precision**: ~15 decimal digits
**Performance**: slow (seconds)
**Security**: 128-bit secure

---

### Configuration 4: Deep Computation

**Use case**: Applications with many multiplications (ML, polynomial evaluation)

```rust
const DEGREE: usize = 32768;
const SCALE_BITS: u32 = 40;
const MULT_DEPTH: usize = 15;

// Modulus chain: [60, 40*15, 60]
let primes = vec![
    // First prime: 60 bits
    1152921504606584833,

    // Intermediate primes: 40 bits * 15
    // ... (15 primes close to 2^40)

    // Special prime: 60 bits
    1152921504606584837,
];

// Total Q: 60 + 15*40 + 60 = 720 bits < 881 bits limit
```

**Precision**: ~10-12 decimal digits
**Performance**: very slow (minutes)
**Security**: 128-bit secure
**Depth**: Supports up to 15 multiplications

---

### Comparison Summary

| Config          | N     | Scale | Depth | Total Q | Precision | Speed  | Security       |
|-----------------|-------|-------|-------|---------|-----------|--------|----------------|
| Toy             | 512   | 2^40  | 3     | ~175b   | 10 digits | fast   | no (educational) |
| Balanced        | 16384 | 2^40  | 5     | ~320b   | 10 digits | medium | 128-bit        |
| High Precision  | 16384 | 2^50  | 4     | ~320b   | 15 digits | slow   | 128-bit        |
| Deep            | 32768 | 2^40  | 15    | ~720b   | 10 digits | very slow | 128-bit     |

---

## Implementation Guidelines

### Step 1: Define Requirements

Before selecting parameters, answer:

1. **What precision do you need?**
   - Financial calculations: 12-15 decimal digits -> scale ~= 2^45
   - Scientific computing: 10-12 digits -> scale ~= 2^35
   - Integer-only: 6-8 digits -> scale ~= 2^25

2. **How many multiplications?**
   - Polynomial evaluation degree d: depth ~= log2(d)
   - Neural network L layers: depth ~= L
   - Simple operations (add/mult): depth = 1-2

3. **What's your performance budget?**
   - Real-time (<100ms): N <= 8192
   - Interactive (<1s): N <= 16384
   - Batch processing (seconds): N <= 32768
   - Offline (minutes): N >= 65536

### Step 2: Calculate Parameters

**Given**: Required precision and depth

**Calculate**:
```rust
// 1. Choose scale from precision requirement
let scale_bits = match required_decimal_digits {
    0..=6 => 25,
    7..=9 => 30,
    10..=12 => 40,
    13..=15 => 50,
    _ => 55,
};

// 2. Calculate total modulus needed
let total_modulus_bits =
    60 +                          // first prime
    (mult_depth * scale_bits) +   // intermediate primes
    60;                           // special prime

// 3. Find minimum N for security
let ring_dim = match total_modulus_bits {
    0..=54 => 2048,
    55..=109 => 4096,
    110..=218 => 8192,
    219..=438 => 16384,
    439..=881 => 32768,
    _ => panic!("Modulus too large, reduce depth or scale"),
};

// 4. Adjust if needed for performance
let ring_dim = ring_dim.max(minimum_n_for_performance);
```

### Step 3: Generate Prime Chain

```rust
fn generate_modulus_chain(
    ring_dim: usize,
    scale_bits: u32,
    mult_depth: usize,
) -> Vec<u64> {
    let mut primes = Vec::new();

    // First prime: 60 bits
    primes.push(find_closest_ntt_prime(ring_dim, 60));

    // Intermediate primes: match scale
    for _ in 0..mult_depth {
        primes.push(find_closest_ntt_prime(ring_dim, scale_bits));
    }

    // Special prime: 60 bits
    primes.push(find_closest_ntt_prime(ring_dim, 60));

    primes
}
```

### Step 4: Validate Configuration

```rust
fn validate_parameters(
    ring_dim: usize,
    primes: &[u64],
) -> Result<(), String> {
    // Check all primes are NTT-friendly
    let modulus = (2 * ring_dim) as u64;
    for (i, &prime) in primes.iter().enumerate() {
        if prime % modulus != 1 {
            return Err(format!(
                "Prime {} at index {} is not NTT-friendly",
                prime, i
            ));
        }
    }

    // Check security constraint
    let total_bits: u32 = primes.iter()
        .map(|p| (p.ilog2() + 1))
        .sum();

    let max_bits = get_max_modulus_bits(ring_dim);
    if total_bits > max_bits {
        return Err(format!(
            "Total modulus {} bits exceeds security limit {} bits for N={}",
            total_bits, max_bits, ring_dim
        ));
    }

    Ok(())
}
```

### Step 5: Document Your Choice

Always document parameter choices in code:

```rust
/// CKKS parameters for high-precision financial calculations
///
/// Security: 128-bit (based on HE standards for N=16384)
/// Precision: ~15 decimal digits (scale = 2^50)
/// Depth: 4 multiplications
///
/// Modulus chain:
///   q_0 = 1152921504606584833  (60 bits, first prime)
///   q_1 = 1125899906842597     (50 bits, matches scale)
///   q_2 = 1125899906842601     (50 bits)
///   q_3 = 1125899906842613     (50 bits)
///   q_4 = 1125899906842619     (50 bits)
///   q_s = 1152921504606584837  (60 bits, special prime)
///
/// Total Q ~= 320 bits < 438 bits limit
const PARAMS: CKKSParams = CKKSParams {
    ring_dim: 16384,
    scale_bits: 50,
    primes: &[
        1152921504606584833,
        1125899906842597,
        1125899906842601,
        1125899906842613,
        1125899906842619,
        1152921504606584837,
    ],
};
```

---

## References

### Primary Sources

- Microsoft SEAL CKKS Basics Example:
  https://github.com/microsoft/SEAL/blob/main/native/examples/5_ckks_basics.cpp
- SEAL CoeffModulus Selection - GitHub Issue #128:
  https://github.com/microsoft/SEAL/issues/128
- Lattigo CKKS Documentation v5:
  https://pkg.go.dev/github.com/tuneinsight/lattigo/v5/schemes/ckks
- OpenFHE Advanced Real Numbers Example:
  https://github.com/openfheorg/openfhe-development/blob/main/src/pke/examples/advanced-real-numbers-128.cpp
- OpenFHE RNS-CKKS Parameter Discussion:
  https://openfhe.discourse.group/t/rns-ckks-parameter-setting/134

### Tutorials and Guides

- CKKS Explained Part 5: Rescaling - OpenMined:
  https://openmined.org/blog/ckks-explained-part-5-rescaling/
- Setting SEAL Parameters - Intuit Engineering (Medium):
  https://medium.com/intuit-engineering/data-science-without-seeing-data-how-to-set-microsoft-open-source-seal-parameters-72929b184058
- Introduction to SEAL - TU Graz:
  https://www.isec.tugraz.at/wp-content/uploads/2021/08/seal.pdf

### Research Papers

- Approximate Homomorphic Encryption with Reduced Approximation Error:
  https://eprint.iacr.org/2020/1118.pdf
- Grafting: Decoupled Scale Factors and Modulus in RNS-CKKS:
  https://eprint.iacr.org/2024/1014.pdf
- Introduction to CKKS - Yongsoo Song (Slides):
  https://yongsoosong.github.io/files/slides/intro_to_CKKS.pdf

### Standards

- Homomorphic Encryption Security Standard:
  http://homomorphicencryption.org/standard/

---

**Last updated**: 2025-12-01
**Next review**: When implementing parameter selection utilities
