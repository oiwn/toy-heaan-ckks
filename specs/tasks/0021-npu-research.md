---
id: 0021-npu-research
status: draft
scope: [specs/research/npu.md]
created: 2026-10-17
source: operator request 2026-10-17
depends: [0018-op-benchmarks]
---
# Task: npu research

Can Apple's neural hardware (Neural Engine, M5 GPU neural accelerators) speed up NTT or key switching? Research only.

## Plan

- [ ] Find how the hardware is reachable: Neural Engine only through Core ML (no public low-level API, to confirm); M5 per-core neural accelerators through Metal.
- [ ] Study the CROSS idea: high-precision modular arithmetic turned into low-precision (INT8) matrix multiplications; check whether it maps to Apple's supported formats.
- [ ] Write `specs/research/npu.md`: reachability, a possible mapping (NTT as matmul), expected gain against the 0018 numbers, go/no-go.

## References

- CROSS, HE on AI ASICs (TPU), HPCA 2026: https://arxiv.org/pdf/2501.07047
- CROSS project page (Microsoft Research): https://www.microsoft.com/en-us/research/?p=1163453
- TensorFHE (NTT on tensor cores): https://ar5iv.labs.arxiv.org/html/2212.14191
- Apple Silicon domain-specific architectures for general-purpose work: https://arxiv.org/pdf/2511.13450
- Apple M5 neural accelerators (BaseRT): https://arxiv.org/pdf/2607.19438
- Apple M-series for HPC: https://arxiv.org/pdf/2502.05317
- Apple M5 (hardware summary): https://en.wikipedia.org/wiki/Apple_M5

## Acceptance

- [ ] `specdev check`
- [ ] `specs/research/npu.md` covers reachability, mapping and go/no-go

## Log

- 2026-10-17 created: npu research
- 2026-10-17 set source operator request 2026-10-17
- 2026-10-17 set depends 0018-op-benchmarks
- 2026-10-17 scope add specs/research/npu.md — research note
