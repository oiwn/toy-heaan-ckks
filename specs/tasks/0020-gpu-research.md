---
id: 0020-gpu-research
status: draft
scope: [specs/research/gpu.md]
created: 2026-10-17
source: operator request 2026-10-17
depends: [0018-op-benchmarks]
---
# Task: gpu research

Could a GPU backend run the RNS-NTT ring ops, and how would it fit the `PolyRing` trait? Research only.

## Plan

- [ ] Read the references; note which ops they move to the GPU and where the bottleneck is (memory bandwidth, per the 100x paper).
- [ ] Survey Rust options that also cover Apple GPUs and the browser: `cubecl` (wgpu/Metal/CUDA/Vulkan), plain `wgpu` compute shaders; note 64-bit integer support per target.
- [ ] Sketch how a GPU backend maps onto `PolyRing` (device-resident polys, batching, transfer cost) and what the trait would have to change.
- [ ] Write `specs/research/gpu.md`: options, fit with the trait, expected gain against the 0018 numbers at toy and demo N.

## References

- Over 100x faster bootstrapping with GPUs (TCHES 2021): https://eprint.iacr.org/2021/508
- TensorFHE (NTT on tensor cores): https://ar5iv.labs.arxiv.org/html/2212.14191
- Cheddar, CUDA CKKS library: https://arxiv.org/pdf/2407.13055
- FIDESlib, open-source CKKS on GPU, OpenFHE-compatible: https://par.nsf.gov/biblio/10612798-fideslib-fully-fledged-open-source-fhe-library-efficient-ckks-gpus
- Efficient NTT on GPU for HE: https://eprint.iacr.org/2021/124.pdf
- Homomorphic Encryption on GPU (survey): https://eprint.iacr.org/2022/1222.pdf
- FHECore, GPU microarchitecture for FHE: https://arxiv.org/pdf/2602.22229
- CubeCL (Rust GPU compute): https://docs.rs/crate/cubecl/latest

## Acceptance

- [ ] `specdev check`
- [ ] `specs/research/gpu.md` covers options, trait fit and estimates

## Log

- 2026-10-17 created: gpu research
- 2026-10-17 set source operator request 2026-10-17
- 2026-10-17 set depends 0018-op-benchmarks
- 2026-10-17 scope add specs/research/gpu.md — research note
