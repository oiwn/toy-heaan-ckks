---
name: educational-examples
description: How to write and change runnable examples in this educational CKKS library - narrated, colored, checked walk-throughs, one per scheme step. Use whenever you add, rewrite or review anything under examples/.
---

# Educational examples

toy-heaan-ckks is an educational library. The examples are the teaching
surface: one runnable example per CKKS step (encode, keys, encrypt, add,
multiply, relinearize, rescale, mod switch, key switch, rotate, conjugate,
bootstrap stages, applications). A reader should learn the step from the
program's output alone, without opening the source.

Reference style: `examples/modmul_methods.rs`.

## Shape of an example

Print these sections, in this order:

1. **Intro** - the problem, and why this step exists in CKKS. Explain why
   before what. Two or three short paragraphs at most.
2. **Parameters** - print every parameter (N, primes and their bits, scale,
   secret weight, seed) with `param(label, value, reason)`: one line each,
   the reason for the value dimmed after it.
3. **Worked steps** - for each operation: a narration line, then
   `formula = value` lines with `step`. Show level, scale and `logq` when
   they change. When the real numbers are unreadable (61-bit primes,
   N=8192), first show a small toy: decimal stand-ins (`q = 97`, `10^4` for
   `2^128`) or N=8.
4. **Check** - compare against a plain f64 model of the same step, print
   the max error against its bound with `check`, and on a miss return `Err`
   from `main` (nonzero exit). Never weaken a bound to pass.
5. **Takeaways** - one short paragraph: what to remember, the limits, and
   which example comes next.

Show failures on purpose when they teach (a wrong key, an overflowing
scale, a level mismatch): print them with `bad` and say why they fail.

The module doc (`//!`) states what the example shows, the run command
(`--release` for bootstrap-sized examples), and the color note: "Colors
are on in a terminal and off when piped or with `NO_COLOR=1`."

## Colors

Use the shared helpers in `examples/support/narrate.rs` only:

```rust
#[path = "support/narrate.rs"]
mod narrate;
use narrate::{bad, check, dim, formula, heading, ok, param, say, step, value};
```

Printing helpers: `say` (a wrapped paragraph), `step` (`formula = value`),
`param` (a parameter and its reason), `check` (max error against a bound,
green or red; returns whether it passed), `confirm` (a yes/no property,
green or red; returns it), `table` (an ASCII table of plain cells),
`elapsed` (an ASCII duration; `{:?}` on a `Duration` prints a non-ASCII
micro sign). Do not
use the library's `toy_heaan_ckks::table` (Unicode box drawing).

Each color has one meaning:

| Helper | Style | Use |
|---|---|---|
| `heading` | bold cyan | section titles |
| `subheading` | bold | sub-sections |
| `formula` | yellow | formulas and their labels |
| `value` | bold | computed numbers |
| `ok` | green | passed checks |
| `bad` | red | failures, wrong results shown on purpose |
| `dim` | dimmed | side notes, references |

- Helpers go through `if_supports_color(Stream::Stdout, ..)`; never print
  raw ANSI codes.
- Pad before coloring (`format!("{x:>8}")`, then color), so escape codes do
  not break table columns.
- Need a new meaning? Add a helper to `narrate.rs` and this table, not a
  one-off color in an example.

## Text

- ASCII only: no Greek letters, emoji, box drawing or arrows (`->` is fine).
  Names follow AGENTS.md (`scale`, not Delta; `std_dev`, not sigma).
- Wrap narration at about 76 columns.
- Name the reference when one applies: fhetextbook.github.io sections,
  OpenFHE or Lattigo. Never Microsoft SEAL.
- Keep the output to about one or two screens; move detail into `dim`
  notes or drop it.

## Correctness

- Fixed RNG seed (`ChaCha20Rng::seed_from_u64`), so output is reproducible.
- Toy parameters that run in seconds in debug; bootstrap examples document
  `--release`.
- Decrypting in the middle of a pipeline is a diagnostic only: label it
  as such ("decrypting here to show the stage").
- Demo-only code stays in `examples/`, never in the library API.
- Avoid `unwrap` in `main`; return `Result<(), Box<dyn std::error::Error>>`.

## Checklist when adding or changing an example

- [ ] Sections in the order above; every number explained or labeled.
- [ ] `check` against a plaintext model, nonzero exit on failure.
- [ ] Row in the README Examples table under its step.
- [ ] `cargo build --examples`, then run it in a terminal and with
      `NO_COLOR=1` (plain output must still read well).
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`.
