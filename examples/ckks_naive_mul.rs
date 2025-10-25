use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use toy_heaan_ckks::{
    Ciphertext, CkksEngine, NaivePolyRing, Plaintext, PolyRing,
    crypto::operations::multiply_ciphertexts,
};

// CT Multiplication Constants from spec
const DEGREE: usize = 8;
const SCALE_BITS: u32 = 20; // Δ = 2^20
const Q_L: u64 = (1u64 << 61) - 1; // Top modulus (level 1)
const Q_0: u64 = (1u64 << 41) - 9; // Bottom modulus (level 0)

// Relinearization parameters
const BASE_LOG2: u32 = 16; // β = 2^16
const DIGITS: u32 = 4; // l = ceil(61/16) = 4

type Engine = CkksEngine<NaivePolyRing<DEGREE>, DEGREE>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = ChaCha20Rng::from_seed([42u8; 32]);

    Ok(())
}
