use crate::PolyRing;
use rand::Rng;
use rand_distr::{Distribution, Normal};

/// Generate a polynomial with uniformly random coefficients
pub fn generate_random_poly<R: Rng>(
    n: usize,
    modulus: u64,
    ring_dim: usize,
    rng: &mut R,
) -> PolyRing {
    let mut coeffs = Vec::with_capacity(n);

    for _ in 0..n {
        // Use gen_range for unbiased sampling
        let random_val = rng.random_range(0..modulus);
        coeffs.push(random_val);
    }

    PolyRing::from_coeffs(&coeffs, modulus, ring_dim)
}

/// Generate a polynomial with coefficients from discrete Gaussian distribution
pub fn generate_error_poly<R: Rng>(
    n: usize,
    modulus: u64,
    variance: f64,
    ring_dim: usize,
    rng: &mut R,
) -> PolyRing {
    let std_dev = variance.sqrt();
    let normal = Normal::new(0.0, std_dev).expect("Invalid standard deviation");
    let mut coeffs = Vec::with_capacity(n);

    for _ in 0..n {
        // Sample from Gaussian distribution and round to nearest integer
        let noise = normal.sample(rng);
        let noise_int = noise.round() as i64;

        // Convert to unsigned, correctly handling negative values with modular arithmetic
        let coeff = if noise_int < 0 {
            modulus - (noise_int.unsigned_abs() % modulus)
        } else {
            noise_int as u64 % modulus
        };

        coeffs.push(coeff);
    }

    PolyRing::from_coeffs(&coeffs, modulus, ring_dim)
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
    let mut neg_coeffs =
        Vec::with_capacity((poly.polynomial_degree() + 1) as usize);

    for coeff in poly {
        if *coeff == 0 {
            neg_coeffs.push(0);
        } else {
            neg_coeffs.push(modulus - *coeff);
        }
    }

    PolyRing::from_coeffs(&neg_coeffs, modulus, 8)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_uniform_distribution() {
        let mut rng = ChaCha20Rng::seed_from_u64(0);
        let modulus = 17;
        let poly_len = 8;
        let samples_count = 100000;

        let mut counts = vec![0usize; modulus as usize];

        for _ in 0..samples_count {
            let poly = generate_random_poly(poly_len, modulus, poly_len, &mut rng);
            for &coeff in poly.into_iter() {
                counts[coeff as usize] += 1;
            }
        }

        let expected = samples_count as f64 * poly_len as f64 / modulus as f64;
        let tolerance_ratio = 1.2;

        for &count in &counts {
            assert!(
                (count as f64 - expected).abs() / expected < tolerance_ratio,
                "Distribution not uniform enough: count {}, expected {}",
                count,
                expected
            );
        }
    }

    #[test]
    fn test_gaussian_distribution() {
        let mut rng = ChaCha20Rng::seed_from_u64(1);
        let modulus = 17;
        let poly_len = 8;
        let variance = 3.0;
        let samples_count = 10000;

        let mut values = Vec::with_capacity(samples_count * poly_len);

        for _ in 0..samples_count {
            let poly = generate_error_poly(
                poly_len, modulus, variance, poly_len, &mut rng,
            );
            for &coeff in poly.into_iter() {
                let val = coeff as i64;
                let val = if val > modulus as i64 / 2 {
                    val - modulus as i64
                } else {
                    val
                };
                values.push(val as f64);
            }
        }

        let mean: f64 = values.iter().copied().sum::<f64>() / values.len() as f64;
        let std_dev = (values.iter().map(|x| (x - mean).powi(2)).sum::<f64>()
            / values.len() as f64)
            .sqrt();

        let expected_std_dev = variance.sqrt();
        let tolerance = 0.2;

        assert!(
            (std_dev - expected_std_dev).abs() < tolerance,
            "Standard deviation deviates too much: actual {}, expected {}",
            std_dev,
            expected_std_dev
        );
    }
}
