# Current task context

## New session restart (read this first)

### Goal right now
- Rebuild `rings/backends/rns_ntt` from scratch with a cleaner structure.

### Finalized file layout
- `errors.rs`
- `basis.rs` (contains both `RnsBasis` and `NttTable`)
- `poly.rs` (single `RnsPoly` type; includes modular helpers + NTT kernels)
- `mod.rs` (exports only)

### Explicit constraints
- Do not split coeff/ntt into separate polynomial types.
- Do not create `poly_ntt.rs`.
- Do not create `modular.rs`.
- Do not create standalone `ntt.rs`.

### Immediate next implementation steps
- Replace current `rns_ntt/poly.rs` scaffold with the finalized single-type design.
- Extract errors into `rns_ntt/errors.rs`.
- Move basis + table ownership into `rns_ntt/basis.rs`.
- Wire `rns_ntt/mod.rs` exports to the new files.
- Add roundtrip and mod-drop tests before removing old backends.

## Active task: rebuild `rns_ntt` backend (composition first)

### Decisions agreed
- Keep trait compatibility with `src/rings/traits.rs` as the primary contract.
- Implement `PolyRing` and `PolySampler` for coeff-domain `RnsPoly<const DEGREE: usize>`.
- Keep `RnsNttPoly<const DEGREE: usize>` as an auxiliary NTT-domain type.
- Expose explicit conversion methods: `RnsPoly::to_ntt()` and `RnsNttPoly::to_coeff()`.
- Use per-channel storage for clarity and mod-drop ergonomics: `Vec<[u64; DEGREE]>`.
- Keep modulus and NTT table alignment strict: channel `i` <-> modulus `q_i` <-> table `i`.

### Planned module split (finalized)
- `errors.rs`: backend-specific errors.
- `basis.rs`: `RnsBasis<DEGREE>` and `NttTable<DEGREE>` in one file; basis owns moduli and NTT tables.
- `poly.rs`: single `RnsPoly<DEGREE>` type with per-channel storage, domain state (`Coeff`/`Ntt`), arithmetic, modular helpers, and NTT kernels inline.
- `mod.rs`: exports only.

### Explicit non-goals for this phase
- No `poly_ntt.rs`.
- No `modular.rs`.
- No separate `ntt.rs`.
- No split coeff/ntt polynomial types.

### Invariants to enforce
- `channels.len() == basis.moduli.len() == basis.ntt_tables.len()`.
- Every coefficient in channel `i` is reduced modulo `q_i`.
- `DEGREE` is power-of-two for NTT-enabled contexts.
- Mod-drop removes aligned suffixes from moduli, tables, and channels.

### Rollout sequence
- Build basis and polynomial skeletons first.
- Add coeff <-> NTT conversion and roundtrip tests.
- Add trait operation smoke tests on `RnsPoly`.
- Only then start removing old `rns` / `ntt` backends.

### Operator note
- Record every active-task architecture/API decision in this file as it is made.
