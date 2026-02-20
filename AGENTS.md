# Toy HEAAN-CKKS Agent Guide

**Project:** Educational Rust implementation of CKKS (Cheon-Kim-Kim-Song)
homomorphic encryption scheme with multiple backends.

## Quick Reference

**Build & Test**
```bash
cargo build && cargo test                    # build and run tests
cargo clippy --all-targets -- -D warnings  # lint
cargo clippy -- -D clippy::pedantic -D clippy::nursery  # also lint
cargo run --example rns                     # end-to-end CKKS flow
cargo run --example ntt                     # verbose domain logging
cargo run --example rns_mod_switch          # modulus-switching demo
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
- prefer to use plain ascii instead of unicode for math!

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
- Running task context: `specs/ctx.md` (keep this file up to date and remind the operator to record every active-task decision there)
- FHE Textbook: https://fhetextbook.github.io - Comprehensive reference for FHE schemes
- Build and test all examples to verify implementation

**Specification Rhythm**
- Kick off every dev cycle by drafting or reviewing a spec: start with the rough idea, iterate with the AI until it becomes a crisp PRD with acceptance criteria, and treat requirements/design as living, versioned artifacts that move through branches like code.

**Human-in-the-loop Rule**
- Whenever you're uncertain about the correct approach or hit a decision point with multiple viable options, pause and ask the operator for guidance before proceeding.

**Important: Naming Conventions**
- Use ASCII-only variable names and formulas
- Replace Greek symbols with descriptive English names:
  - `Delta` -> `scale` or `scaling_factor`
  - `sigma` -> `std_dev` (standard deviation)
  - `epsilon` -> `error` or `noise`
  - `q`, `Q` -> `modulus`, `extended_modulus`
  - `alpha`, `beta` -> use descriptive names like `base`, `digit`
- This ensures code readability and avoids encoding issues
