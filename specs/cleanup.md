# Cleanup (code smells)

Read on demand — not session context.

## `CkksParams::errors_variance` is actually a std_dev

All call sites pass it to `sample_gaussian(std_dev, ...)`; rename when next touching `src/crypto/engine.rs`.

## `src/crypto/operations.rs` is a stale duplicate of `engine.rs`

Tagged FIXME for removal.

## `PolyRescale` / `PolyModSwitch` traits are declared but unimplemented

The real rescaling lives in `RnsPoly::rescale()` / `CkksEngine::rescale_ciphertext()`.

## `RnsPoly::channels_mut` is dead code

Added for the `div_pow2` experiment (since removed — modular-inverse halving is invalid in CKKS). No users; remove or keep deliberately.
