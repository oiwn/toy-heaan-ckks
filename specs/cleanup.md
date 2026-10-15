# Cleanup (code smells)

Read on demand — not session context.

## `CkksParams::errors_variance` is actually a std_dev

All call sites pass it to `sample_gaussian(std_dev, ...)`; rename when next touching `src/ckks/engine.rs`.