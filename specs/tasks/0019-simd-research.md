---
id: 0019-simd-research
status: draft
scope: [specs/research/simd.md]
created: 2026-10-16
source: operator request 2026-10-16
depends: [0018-op-benchmarks]
---
# Task: simd research

Where could SIMD speed up the RNS-NTT backend, and at what cost to readability? Research only; M2.5 stays free of optimizations.

## Plan

- [ ] Profile with the 0018 benches: which ops dominate (NTT butterflies, pointwise mul, key switch).
- [ ] Survey options: compiler auto-vectorization, `wide`/`pulp` crates, `std::simd` (nightly; toolchain changes need operator approval), Shoup/Barrett per element.
- [ ] Note the hardware limits: no 64x64->128 vector mul in AVX2/NEON, AVX-512 IFMA with 52-bit primes (Intel HEXL), wasm simd128.
- [ ] Write `specs/research/simd.md`: options, expected speedup from the 0018 numbers, readability cost, recommendation.

## Acceptance

- [ ] `specdev check`
- [ ] `specs/research/simd.md` covers profile, options, limits and estimates

## Log

- 2026-10-16 created: simd research
- 2026-10-16 set source operator request 2026-10-16
- 2026-10-16 set depends 0018-op-benchmarks
- 2026-10-16 scope add specs/research/simd.md — research note
