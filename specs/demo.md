# Browser FHE demo specification

## Goal & story

A single-page browser demo, fully client-side, telling one story end to end:

> A temperature sensor collected 24 hours of readings. The readings are
> encrypted **in the browser**. The standard deviation is computed **on the
> ciphertext** — the "server" code never sees a plaintext value. The result is
> decrypted and compared against canonical f64 math, digit by digit.

Everything runs in the browser via WASM: key generation, encryption, the homomorphic pipeline, decryption. No server, no network calls.

```text
ct = Enc([t0 .. t31])
  -> sum_slots (rot 1,2,4,8,16)  -> mean = sum/32     -> centered = ct - mean
  -> squared = centered²          -> sum_slots again   -> variance = sum/32
  -> v_norm = variance / B        -> chebyshev sqrt    -> std = result · sqrt(B)
```

## Decisions (agreed 2026-07-05)

| Topic | Decision | Rationale |
|-------|----------|-----------|
| Frontend | Leptos (CSR) + trunk, own crate `demo-web/` | Rust end to end; SPA fits the staged UI |
| Compute placement | 100% client-side WASM, heavy ops in a Web Worker | No server needed; worker keeps UI responsive |
| Fallback server | none (axum rejected) | Lib compiles to wasm32; no blocker found |
| Storage | In-memory SPA state + export/import files | fjall is a disk LSM — cannot run in a browser; IndexedDB deferred |
| Ring dimension | `N = 8192` (4096 slots) | Browser-viable keys/compute; "huge key" story intact |
| Prime chain | 12 × 61-bit NTT primes, scale = 2^61 | Prime size the whole codebase is validated at; 1e-13 noise floor makes the comparison panel shine |
| DIGIT_BITS | bump 15 → 20 (4 digits per 61-bit prime) | ~20% smaller keys and fewer key-switch muls; noise still ≪ scale |
| Security labeling | Honest banner: reduced security | logQ ≈ 732 > 218-bit budget for 128-bit at N=8192 (see below) |
| Data | 32 readings = one per 45 min over 24h | 5 rotation keys only; matches sqrt envelope |
| Sqrt approx | degree-7 Chebyshev on [0.01, 1], B = 100 | Decided earlier — see ctx.md "Sqrt approximation decision" |

### Security banner (display in UI, verbatim-ish)

> **Educational parameters.** RLWE security depends jointly on the ring
> dimension N and the total modulus size logQ. This demo uses N = 8192 with
> logQ ≈ 732 bits — far above the ≈218-bit budget that N = 8192 allows for
> 128-bit security — so its effective security is roughly 40–60 bits.
> True 128-bit security for this 12-level pipeline needs N = 32768, where the
> gadget keys grow to ~2.3 GB: not browser material. The *mechanics* shown
> here are exactly the real ones.

### Size math (drives the UI byte displays)

```text
poly          = l · N · 8 bytes           = 12 · 8192 · 8   ≈ 786 KB
ciphertext    = 2 polys                                     ≈ 1.6 MB
gadget key    = 2 · (l·d) polys = 16·l²·d·N
              = 16 · 144 · 4 · 8192                         ≈ 75 MB
key material  = 1 relin + 5 rotation keys                   ≈ 450 MB RAM
secret key    = 1 poly                                      ≈ 786 KB
public key    = 2 polys                                     ≈ 1.6 MB
```

Level consumption across the pipeline (12 primes):

```text
fresh ct           level 12
mean (·1/32)       11
center (sub)       11   (no level)
square             10
variance (·1/32)   9
normalize (·1/B)   8
chebyshev deg 7    ~1    (7 Clenshaw muls + domain map)
scale (·√B)        folded into last plaintext mul
```

## UI flow (one page, staged top to bottom)

1. **Data panel.** Seeded generator produces 32 readings from a diurnal curve + jitter (std ≈ 2–3 °C, inside the [1, 10] °C envelope); readings editable by hand; regenerate button. Shows canonical f64 mean / variance / std-dev alongside.
2. **Parameters panel.** N, prime list, logQ, scale, slot count, expected key sizes, and the security banner.
3. **Key generation.** Button starts the worker; progress bar per key (sk, pk, relin, rot×5). After: table of actual byte sizes + generation time per key.
4. **Encrypt.** Encode + encrypt in the worker. Display: ciphertext byte size, first-bytes hex preview, logp / logq / level, wall time. Export button (download binary file).
5. **Compute.** Runs the pipeline stage by stage; each stage row shows the operation, level before → after, scale bits, wall time. The intent is that the level-consumption table above becomes *visible live*.
6. **Decrypt & compare.** Decrypt + decode in the worker. Comparison table:

   | | value |
   |---|---|
   | homomorphic std-dev (decrypted) | 2.79052... |
   | canonical f64 std-dev | 2.79054... |
   | absolute error | ~1e-2 max (Chebyshev-dominated) |
   | encrypted-variance error (intermediate) | ~1e-13 (CKKS noise floor) |

   The two error rows are the punchline: crypto noise is invisible; the polynomial approximation is the entire error budget.
7. **Export / import.** Download ciphertext / result as files; upload a previously exported ciphertext and continue from step 5.

## Crate layout

```text
toy-heaan-ckks/            existing lib (root package)
  Cargo.toml               + exclude = ["demo-web"]; deps cleaned for wasm
demo-web/
  Cargo.toml               leptos (csr), wasm-bindgen, web-sys, postcard,
                           serde, path dep on toy-heaan-ckks
  Trunk.toml, index.html
  src/
    main.rs                Leptos app, panels/state
    worker.rs              Web Worker entry: keygen/encrypt/pipeline/decrypt
    messages.rs            worker protocol (serde), progress events
```

Built with `trunk serve` / `trunk build --release` inside `demo-web/`; the root workspace never builds it (excluded), so `cargo test --workspace` stays native-only.

## Gap analysis — missing puzzle pieces

Status of every part the demo needs:

| Piece | Status |
|-------|--------|
| Encrypt/decrypt, add/sub, mul+relin+rescale, rotate | ✅ works (61-bit verified) |
| Base-β gadget key switch | ✅ done 2026-07-04 |
| Chebyshev sqrt coefficients + envelope | ✅ decided, plaintext-verified |
| Encrypted Chebyshev eval (Clenshaw) | ✅ done 2026-08-17 |
| Pipeline ops as library API | ✅ `src/crypto/pipeline.rs` |
| **Serialization (serde) for ct/keys** | ❌ no derives yet (serde already a dep) |
| **wasm32 build of the lib** | ❌ untested; `dhat`/`comfy-table`/`tracing-subscriber` sit in `[dependencies]`, getrandom needs the wasm backend |
| DIGIT_BITS = 20 | ✅ done 2026-08-15 |
| Browser performance numbers | ❓ unmeasured (estimates below) |

Performance estimates to validate in the spike (Phase B):

```text
key-switch  ≈ l·d = 48 digit muls ≈ 1700 channel-NTTs  → ~0.3-1 s wasm
pipeline    ≈ 10 rotations + ~10 muls                   → ~15-60 s wasm
keygen      ≈ 288 key-poly muls + 450 MB alloc          → ~10-30 s wasm
encode/decode (O(N²) special DFT)                       → ~1-2 s each, rare
```

Contingency if too slow: keep gadget keys and key-switch accumulators in NTT domain (saves ~3× NTTs per key-switch). Only implement if the spike says so.

## Work plan

### Phase A — library gaps (native, testable without any browser)

Status 2026-10-12: steps 1, 3, 4 done (see CHANGELOG). Steps 5–6 are now tasks 0002/0003 (roadmap M1/M2); step 2 moved to roadmap M3. The 1/B normalization and the 12-prime count below are superseded by the sensor range contract (2026-10-12): range [-40, 60] degC, std envelope [1, 8] degC, degree 7, 13 x 61-bit primes, std restored on the client (see `specs/overview.md` gotchas).

1. `DIGIT_BITS` 15 → 20; re-run `cargo test`, `horner_chain`, `std_dev_8`.
2. Serde derives for `RnsPoly`, `Ciphertext`, `Plaintext`, key types (`RnsBasis` reconstructed from prime list, not serialized wholesale); postcard round-trip test asserting exact byte sizes from the table above.
3. Extract pipeline ops from `examples/std_dev_8.rs` into the library (`src/crypto/pipeline.rs`): `sum_slots`, `mul_plain_scalar`, `align_level` (ct+sk reduction), wide f64 decode.
4. Encrypted Chebyshev eval via Clenshaw in Chebyshev basis (ctx.md Phase 2): domain map is a plaintext affine op; each recurrence step one ct×ct mul; b-terms mod-switched to stay aligned. Test vs plaintext Clenshaw.
5. `encrypted_std_dev` composing 3+4 (lives in `examples/sensor_std_dev/`, not the library).
6. New native example `std_dev_demo_params.rs`: N=8192, 12 × 61-bit primes, 32 diurnal readings — asserts std-dev error < 0.11 °C (expect ~1e-2) and variance error < 1e-9. This is the acceptance gate for Phase A.

### Phase B — wasm spike (go/no-go data, no UI)

1. Cargo cleanup: move `dhat`, `comfy-table`, `tracing-subscriber` to dev-dependencies (grep first — lib must not use them); add getrandom wasm backend feature; `cargo build --target wasm32-unknown-unknown` green.
2. Minimal wasm-bindgen harness (no Leptos): run keygen / encrypt / one rotation / one mul / full pipeline in the browser console; record timings and peak memory next to the estimates table.
3. Decide on the NTT-domain-keys optimization based on the numbers.

### Phase C — demo app

1. Scaffold `demo-web/` (trunk + Leptos CSR + worker plumbing + messages).
2. Panels 1–6 from the UI flow, wired to the worker with progress events.
3. Export/import (postcard files, panel 7).

### Phase D — polish

Progress instrumentation granularity, security banner copy, README with screenshots, `trunk build --release` artifact check.

## Verification

- **Phase A**: `cargo test --all --workspace` green; `std_dev_demo_params` example prints and asserts the error budget (std < 0.11 °C, variance < 1e-9).
- **Phase B**: wasm32 target builds; browser console benchmark table recorded in this file (replace the estimates).
- **Phase C/D**: manual walkthrough — generate → keygen → encrypt → compute → decrypt shows both error rows within budget; export then import a ciphertext and re-run compute to the same result; page stays responsive during compute.
