//! CKKS algorithms built only on the ops, the level methods and the
//! evaluator: slot sums and Chebyshev approximation.
//!
//! Applications built from these (e.g. the sensor std-dev demo) live in
//! `examples/`.

pub mod chebyshev;
pub mod sum_slots;
