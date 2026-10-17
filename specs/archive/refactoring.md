# Rings Refactoring Plan

## Current State

```
src/rings/backends/
├── rns/
│   ├── mod.rs
│   ├── poly.rs        # RnsPolyRing (RNS, schoolbook O(n^2) multiplication)
│   ├── ntt_poly.rs    # RnsNttPoly (RNS + NTT)
│   └── params.rs
└── ntt/
    ├── mod.rs
    └── poly.rs        # RnsNttPolyRing (RNS + NTT, DUPLICATE)
```

### Problems

1. **Duplication**: `RnsNttPoly` and `RnsNttPolyRing` are essentially identical
2. **Confusing split**: RNS and NTT are separate modules but NTT already requires RNS
3. **RnsBasis already has NTT tables**: `NttTable` stored per prime (line 283-284 in poly.rs)

## Proposed Architecture

```
src/rings/backends/
└── rns/
    ├── mod.rs
    ├── basis.rs       # RnsBasis (primes + NTT tables)
    ├── poly.rs        # RnsNttPoly - the ONE unified type
    └── params.rs      # Helper functions for toy/test basis
```

### Single Unified Type: RnsNttPoly

```rust
pub struct RnsNttPoly<const DEGREE: usize> {
    coefficients: Vec<[u64; DEGREE]>,  // RNS channels
    in_ntt_domain: bool,                // Domain flag
    basis: Arc<RnsBasis>,               // Shared basis (primes + NTT tables)
}
```

**Capabilities:**
- RNS representation (always)
- `to_ntt_domain()` / `to_coeff_domain()` transforms
- Addition: works in either domain (assert same domain)
- Multiplication: pointwise in NTT domain
- Mod-switching: drop channels

### Files to Delete

- `ntt/` directory entirely
- `rns/poly.rs` (RnsPolyRing with schoolbook) - unless needed for educational purposes

## Basis Compatibility

### Problem

Two polynomials with identical primes but different `Arc` allocations fail `Arc::ptr_eq()`:

```rust
let basis1 = Arc::new(RnsBasis::new(...));
let basis2 = Arc::new(RnsBasis::new(...));  // Same primes, different Arc
// ptr_eq would fail even though bases are compatible
```

### Solutions Considered

| Approach | Pros | Cons |
|----------|------|------|
| Global registry | Guaranteed uniqueness, fast ptr_eq | Global state, threading complexity |
| Compare by value | Simple, no global state | O(n) comparison overhead |
| Sub-basis caching | Fast mod_switch | Can't cache inside Arc (immutable) |

### Chosen Approach: Compare by Value

```rust
impl RnsBasis {
    pub fn is_compatible_with(&self, other: &RnsBasis) -> bool {
        self.degree == other.degree && self.primes == other.primes
    }
}
```

Usage in arithmetic:
```rust
assert!(self.basis.is_compatible_with(&rhs.basis), "basis mismatch");
```

Overhead is negligible compared to polynomial operations.

## Mod-Switching API

### Simple Approach (No Caching)

```rust
impl RnsNttPoly {
    pub fn mod_switch(&self, channels: usize) -> Self {
        assert!(channels <= self.channels(), "cannot add channels");
        // Validate prefix matches (for safety)

        let new_basis = Arc::new(RnsBasis::from_primes(
            &self.basis.primes()[..channels],
            self.basis.degree(),
        ));

        Self {
            coefficients: self.coefficients[..channels].to_vec(),
            in_ntt_domain: self.in_ntt_domain,
            basis: new_basis,
        }
    }
}
```

### Future: Basis Interning

If performance becomes critical, add a global registry:

```rust
lazy_static! {
    static ref BASIS_CACHE: Mutex<HashMap<Vec<u64>, Weak<RnsBasis>>> = ...;
}

impl RnsBasis {
    pub fn intern(primes: Vec<u64>, degree: usize) -> Arc<RnsBasis> {
        // Return cached or create new
    }
}
```

Defer this until needed.

## Migration Steps

1. [ ] Create `basis.rs` with `RnsBasis` (move from poly.rs)
2. [ ] Consolidate `RnsNttPoly` as single type in `poly.rs`
3. [ ] Add `is_compatible_with()` to `RnsBasis`
4. [ ] Update all arithmetic to use value comparison
5. [ ] Delete `ntt/` directory
6. [ ] Delete `RnsPolyRing` (schoolbook version)
7. [ ] Update `mod.rs` exports
8. [ ] Update all examples and tests
9. [ ] Run full test suite

## Open Questions

- Keep schoolbook `RnsPolyRing` for educational comparison?
- Any code currently depending on the `ntt/` module?
