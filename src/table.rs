//! Pretty-printing helpers for CKKS examples.
//!
//! Thin wrapper around `comfy-table` with a dense preset applied by default.
//!
//! # Usage
//! ```ignore
//! use toy_heaan_ckks::table;
//!
//! let mut t = table::new(["slot", "decoded", "expected"]);
//! t.add_row(["0", "7.000001", "7.0"]);
//! println!("{t}");
//! ```

use comfy_table::{Table, presets};

/// Create a dense table with the given column headers.
///
/// Uses `UTF8_FULL_CONDENSED`: unicode box-drawing, no inter-row dividers.
/// Column widths are inferred automatically from content.
pub fn new<H, S>(headers: H) -> Table
where
    H: IntoIterator<Item = S>,
    S: Into<comfy_table::Cell>,
{
    let mut t = Table::new();
    t.load_preset(presets::UTF8_FULL_CONDENSED);
    t.set_header(headers);
    t
}
