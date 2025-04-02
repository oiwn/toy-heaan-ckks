use crate::PolyRing;
use rand::Rng;

// Helper functions needed for public key generation

/// Generate a polynomial with random coefficients in [0, modulus)
pub fn generate_random_poly<T: Rng>(
    n: usize,
    modulus: u64,
    rng: &mut T,
) -> PolyRing {
    let mut coeffs = Vec::with_capacity(n);

    for _ in 0..n {
        // Generate a random value modulo q
        let random_val = rng.random::<u64>() % modulus;
        coeffs.push(random_val);
    }

    PolyRing::from_coeffs(&coeffs, modulus, 8)
}

/// Generate a polynomial with small random coefficients (error polynomial)
/// using discretized Gaussian distribution
pub fn generate_error_poly<T: Rng>(
    n: usize,
    modulus: u64,
    variance: f64,
    rng: &mut T,
) -> PolyRing {
    // Get standard deviation from variance
    let std_dev = variance.sqrt();
    let mut coeffs = Vec::with_capacity(n);

    for _ in 0..n {
        // Simplified approximation of Gaussian noise
        // Sum three random values and shift to get approximately normal distribution
        let noise =
            (rng.random::<f64>() + rng.random::<f64>() + rng.random::<f64>() - 1.5)
                * std_dev;
        let noise_int = noise.round() as i64;

        // Convert to unsigned and reduce modulo q
        let coeff = if noise_int < 0 {
            let abs_noise = noise_int.unsigned_abs();
            modulus - (abs_noise % modulus)
        } else {
            noise_int as u64 % modulus
        };

        coeffs.push(coeff);
    }

    PolyRing::from_coeffs(&coeffs, modulus, 8)
}

/// Generate a ternary polynomial with coefficients in {-1, 0, 1}
/// sparsity controls what fraction of coefficients will be 0
pub fn generate_ternary_poly<T: Rng>(
    n: usize,
    modulus: u64,
    sparsity: f64,
    rng: &mut T,
) -> PolyRing {
    let mut coeffs = Vec::with_capacity(n);

    for _ in 0..n {
        // First decide if coefficient is zero based on sparsity
        if rng.random::<f64>() < sparsity {
            coeffs.push(0);
        } else {
            // Otherwise, generate -1 or 1 with equal probability
            let coeff = if rng.random::<bool>() {
                1u64
            } else {
                // -1 is represented as (modulus - 1) in modular arithmetic
                modulus - 1
            };
            coeffs.push(coeff);
        }
    }

    PolyRing::from_coeffs(&coeffs, modulus, 8)
}

/// Compute the negative of a polynomial in the ring
pub fn negate_poly(poly: &PolyRing, modulus: u64) -> PolyRing {
    // Create a polynomial with coefficients (q - a_i) for each coefficient a_i
    let mut neg_coeffs = Vec::with_capacity((poly.poly_degree() + 1) as usize);

    for coeff in poly {
        if *coeff == 0 {
            neg_coeffs.push(0);
        } else {
            neg_coeffs.push(modulus - *coeff);
        }
    }

    PolyRing::from_coeffs(&neg_coeffs, modulus, 8)
}
