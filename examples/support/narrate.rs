//! Narration and color helpers shared by the examples.
//!
//! Include with `#[path = "support/narrate.rs"] mod narrate;`. Each color
//! has one meaning (see `.agents/skills/educational-examples/SKILL.md`).
//! Colors are on in a terminal and off when piped or with `NO_COLOR=1`.
#![allow(dead_code)]

use std::fmt::Display;

use owo_colors::{OwoColorize, Stream};

/// Section title: bold cyan.
pub fn heading(text: &str) -> String {
    let text = format!("== {text} ==");
    format!(
        "{}",
        text.if_supports_color(Stream::Stdout, |t| t.bold().cyan().to_string())
    )
}

/// Sub-title: bold.
pub fn subheading(text: &str) -> String {
    let text = format!("-- {text} --");
    format!("{}", text.if_supports_color(Stream::Stdout, |t| t.bold()))
}

/// A formula: yellow.
pub fn formula(text: &str) -> String {
    format!("{}", text.if_supports_color(Stream::Stdout, |t| t.yellow()))
}

/// A computed number: bold.
pub fn value(x: impl Display) -> String {
    let text = x.to_string();
    format!("{}", text.if_supports_color(Stream::Stdout, |t| t.bold()))
}

/// A passed check: green.
pub fn ok(text: &str) -> String {
    format!("{}", text.if_supports_color(Stream::Stdout, |t| t.green()))
}

/// A failure, or a wrong result shown on purpose: red.
pub fn bad(text: &str) -> String {
    format!("{}", text.if_supports_color(Stream::Stdout, |t| t.red()))
}

/// A side note: dimmed.
pub fn dim(text: &str) -> String {
    format!("{}", text.if_supports_color(Stream::Stdout, |t| t.dimmed()))
}

/// Prints `label = value`, label as a formula, value in bold.
pub fn step(label: &str, x: impl Display) {
    println!("  {} = {}", formula(label), value(x));
}

/// An ASCII duration: `12.34 ms` or `1.23 s` (Debug formatting would
/// print the non-ASCII `us` sign).
pub fn elapsed(duration: std::time::Duration) -> String {
    let seconds = duration.as_secs_f64();
    if seconds < 1.0 {
        format!("{:.2} ms", seconds * 1e3)
    } else {
        format!("{seconds:.2} s")
    }
}

/// Prints a parameter: `label = value  reason`, the reason dimmed.
pub fn param(label: &str, x: impl Display, reason: &str) {
    let x = format!("{x:<14}");
    println!("  {} = {} {}", formula(label), value(x), dim(reason));
}

/// Prints a narration paragraph (already wrapped).
pub fn say(text: &str) {
    println!("{text}");
}

/// Prints an ASCII table: right-aligned columns, `|` separators, a dashed
/// rule under the header. Cells are plain text (color would break the
/// widths).
pub fn table<R: AsRef<[String]>>(headers: &[&str], rows: &[R]) {
    let mut widths: Vec<usize> = headers.iter().map(|h| h.len()).collect();
    for row in rows {
        for (width, cell) in widths.iter_mut().zip(row.as_ref()) {
            *width = (*width).max(cell.len());
        }
    }
    let line = |cells: Vec<&str>| {
        let padded: Vec<String> = cells
            .iter()
            .zip(&widths)
            .map(|(cell, width)| format!("{cell:>width$}"))
            .collect();
        println!("  {}", padded.join(" | "));
    };
    line(headers.to_vec());
    let total = widths.iter().sum::<usize>() + 3 * (widths.len() - 1);
    println!("  {}", "-".repeat(total));
    for row in rows {
        line(row.as_ref().iter().map(String::as_str).collect());
    }
}

/// Prints a yes/no check in green when it holds, red otherwise, and
/// returns it.
pub fn confirm(text: &str, holds: bool) -> bool {
    if holds {
        println!("  {}", ok(&format!("{text}: yes")));
    } else {
        println!("  {}", bad(&format!("{text}: NO")));
    }
    holds
}

/// Prints `name: max error E (bound B)` in green, or red when E > B, and
/// returns whether it passed.
pub fn check(name: &str, error: f64, bound: f64) -> bool {
    let passed = error.is_finite() && error <= bound;
    let line = format!("{name}: max error {error:.3e} (bound {bound:.0e})");
    println!("  {}", if passed { ok(&line) } else { bad(&line) });
    passed
}
