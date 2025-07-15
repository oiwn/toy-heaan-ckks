// mod basis;
// mod display;
pub mod backends;
pub mod traits;
// mod poly_ring;

pub use backends::NaivePolyRing;
pub use traits::{PolyRing, PolySampler};

#[derive(Debug, Clone, Copy)]
pub enum BackendType {
    Naive,
    // Rns,     // Future
    // BigInt,  // Future
}

impl Default for BackendType {
    fn default() -> Self {
        Self::Naive
    }
}
