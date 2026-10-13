# Cleanup (code smells)

Read on demand — not session context.

## `CkksParams::errors_variance` is actually a std_dev

All call sites pass it to `sample_gaussian(std_dev, ...)`; rename when next touching `src/crypto/engine.rs`.

## `PolyRescale` / `PolyModSwitch` traits are declared but unimplemented

The real rescaling lives in `RnsPoly::rescale()` / `CkksEngine::rescale_ciphertext()`.