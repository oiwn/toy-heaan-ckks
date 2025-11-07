# Toy HEAAN-CKKS Agent Guide

**Project:** Educational Rust implementation of CKKS (Cheon-Kim-Kim-Song)
homomorphic encryption scheme with multiple backends.

## Quick Reference

**Build & Test**
```bash
cargo build && cargo test                    # build and run tests
cargo clippy --all-targets -- -D warnings  # lint
cargo run --example ckks_naive              # Naive demo
cargo run --example ckks_rns                # RNS demo
cargo run --example ckks_bigint             # BigInt demo
cargo bench --bench end_to_end              # benchmarks
```

**Pre-commit Setup**
```bash
pip install pre-commit                      # install pre-commit
pre-commit install                          # install git hooks
```

**Code Style**
- rustfmt (edition 2024), max width ~84 chars
- Imports: std -> external -> internal
- Naming: snake_case (functions/vars), PascalCase (types)
- Errors: `thiserror`, return `Result<T, E>`; avoid panics outside tests
- Types: prefer explicit types; use `const DEGREE: usize` generics
- Docs: `///` with `# Panics` and `# Errors` where relevant
- Clippy: fix warnings; allow specific lints only when necessary

**Key File Paths**
- `src/crypto/` - engine, builder, operations, types
- `src/rings/backends/` - Naive, RNS, BigInt, NTT backends
- `src/encoding/` - FFT, BigInt HEAAN-style, naive encoders
- `src/keys/` - secret/public/relinearization keys
- `examples/` - runnable demos per backend
- `tests/` - encoder tables, properties, e2e tests
- `src/rings/traits.rs` - abstract traits: `PolyRing`, `PolySampler`, `PolyRescale`

**Scope**
- Applies to entire repository; follow conventions for any files you touch
- Targets CKKS/HEAAN development with Naive, RNS, and BigInt (U256) backends
- Educational implementation - not production-ready

## Documentation Links
- Project Architecture: `specs/overview.md`
- Current Tasks & Specifications: `specs/current_task.md`
- FHE Textbook: https://fhetextbook.github.io - Comprehensive reference for FHE schemes
- Build and test all examples to verify implementation

**Important: Naming Conventions**
- Use ASCII-only variable names and formulas
- Replace Greek symbols with descriptive English names:
  - `Delta` -> `scale` or `scaling_factor`
  - `sigma` -> `std_dev` (standard deviation)
  - `epsilon` -> `error` or `noise`
  - `q`, `Q` -> `modulus`, `extended_modulus`
  - `alpha`, `beta` -> use descriptive names like `base`, `digit`
- This ensures code readability and avoids encoding issues
