pub mod encoding;
pub mod keys;
pub mod rings;

pub use keys::{KeyGenParams, PublicKey, SecretKey};
pub use rings::PolyRing;

use crypto_bigint::nlimbs;

// Helper macro to create type aliases with bits
macro_rules! poly_ring_bits {
    ($name:ident, $bits:expr) => {
        pub type $name = PolyRing<{ nlimbs!($bits) }>;
    };
}

// Type aliases using bits instead of limbs
poly_ring_bits!(PolyRing128, 128);
poly_ring_bits!(PolyRing256, 256);
poly_ring_bits!(PolyRing512, 512);
