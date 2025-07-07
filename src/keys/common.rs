use crate::{RnsBasis, RnsPolyRing};
use rand::{Rng, seq::SliceRandom};
use rand_distr::{Distribution, Normal};
use std::sync::Arc;

// Sample uniformly random integer coefficients in range [0, max_value)
pub fn sample_uniform_u64<const DEGREE: usize, R: Rng + ?Sized>(
    max_value: u64,
    rng: &mut R,
) -> [u64; DEGREE] {
    let mut coeffs = [0u64; DEGREE];
    for i in 0..DEGREE {
        coeffs[i] = rng.random_range(0..max_value);
    }
    coeffs
}

/// Sample Gaussian noise as signed integers, then convert to unsigned mod max_value
pub fn sample_gaussian_u64<const DEGREE: usize, R: Rng + ?Sized>(
    std_dev: f64,
    max_value: u64,
    rng: &mut R,
) -> [u64; DEGREE] {
    let normal = Normal::new(0.0, std_dev).expect("Invalid Gaussian std_dev");
    let mut coeffs = [0u64; DEGREE];

    for i in 0..DEGREE {
        let sample = normal.sample(rng);
        let noise_int = sample.round() as i64;

        // Convert to unsigned: handle negative values properly
        coeffs[i] = if noise_int < 0 {
            let abs_val = noise_int.unsigned_abs() % max_value;
            if abs_val == 0 { 0 } else { max_value - abs_val }
        } else {
            (noise_int as u64) % max_value
        };
    }

    coeffs
}

/// Sample a ternary polynomial (coeffs in \{-1,0,1\}) with given Hamming weight.
pub fn sample_ternary_i64<const DEGREE: usize, R: Rng + ?Sized>(
    hamming_weight: usize,
    rng: &mut R,
) -> [i64; DEGREE] {
    let mut out = [0i64; DEGREE];
    // shuffle indices and pick the first `hamming_weight`
    let mut indices: Vec<usize> = (0..DEGREE).collect();
    indices.shuffle(rng);
    for &idx in indices.iter().take(hamming_weight) {
        out[idx] = if rng.random_bool(0.5) { 1 } else { -1 };
    }
    out
}

/// Convert integer coefficients to RNS representation
pub fn u64_to_rns<const DEGREE: usize>(
    coeffs: &[u64; DEGREE],
    basis: Arc<RnsBasis>,
) -> RnsPolyRing<DEGREE> {
    let mut channels: Vec<[u64; DEGREE]> = Vec::with_capacity(basis.primes().len());

    for &q in basis.primes() {
        let mut arr = [0u64; DEGREE];
        for i in 0..DEGREE {
            arr[i] = coeffs[i] % q;
        }
        channels.push(arr);
    }

    RnsPolyRing {
        coefficients: channels,
        basis,
    }
}

/// Convert signed integer coefficients to RNS representation
pub fn i64_to_rns<const DEGREE: usize>(
    coeffs: &[i64; DEGREE],
    basis: Arc<RnsBasis>,
) -> RnsPolyRing<DEGREE> {
    let mut channels: Vec<[u64; DEGREE]> = Vec::with_capacity(basis.primes().len());

    for &q in basis.primes() {
        let mut arr = [0u64; DEGREE];
        let q_i64 = q as i64;
        for i in 0..DEGREE {
            // Ensure positive residue: ((x % q) + q) % q
            arr[i] = ((coeffs[i] % q_i64 + q_i64) % q_i64) as u64;
        }
        channels.push(arr);
    }

    RnsPolyRing {
        coefficients: channels,
        basis,
    }
}

/// Simplified uniform polynomial sampling: integers first, then convert to RNS
pub fn sample_uniform_poly<const DEGREE: usize, R: Rng + ?Sized>(
    basis: Arc<RnsBasis>,
    rng: &mut R,
) -> RnsPolyRing<DEGREE> {
    // Find the largest prime in the basis to use as max_value
    let max_prime = *basis.primes().iter().max().unwrap_or(&1);

    // Sample uniform integers
    let coeffs = sample_uniform_u64::<DEGREE, _>(max_prime, rng);

    // Convert to RNS
    u64_to_rns(&coeffs, basis)
}

/// Simplified Gaussian polynomial sampling: integers first, then convert to RNS
pub fn sample_gaussian_poly<const DEGREE: usize, R: Rng + ?Sized>(
    basis: Arc<RnsBasis>,
    std_dev: f64,
    rng: &mut R,
) -> RnsPolyRing<DEGREE> {
    // Use the product of all primes as the working modulus for Gaussian sampling
    let modulus_product: u64 = basis.primes().iter().product();

    // Sample Gaussian integers
    let coeffs = sample_gaussian_u64::<DEGREE, _>(std_dev, modulus_product, rng);

    // Convert to RNS
    u64_to_rns(&coeffs, basis)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rings::RnsBasisBuilder;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn test_sample_uniform_integers() {
        const DEGREE: usize = 8;
        let mut rng = StdRng::seed_from_u64(42);

        let coeffs = sample_uniform_u64::<DEGREE, _>(100, &mut rng);

        // All coefficients should be in range [0, 100)
        for &coeff in &coeffs {
            assert!(coeff < 100, "Coefficient {} should be < 100", coeff);
        }

        // Should have DEGREE coefficients
        assert_eq!(coeffs.len(), DEGREE);
    }

    #[test]
    fn test_sample_gaussian_integers() {
        const DEGREE: usize = 16;
        let mut rng = StdRng::seed_from_u64(123);

        let coeffs = sample_gaussian_u64::<DEGREE, _>(3.2, 1000, &mut rng);

        // All coefficients should be in range [0, 1000)
        for &coeff in &coeffs {
            assert!(coeff < 1000, "Coefficient {} should be < 1000", coeff);
        }

        // Most coefficients should be small (within a few standard deviations)
        let small_count = coeffs.iter().filter(|&&x| x < 20 || x > 980).count();
        assert!(
            small_count > DEGREE / 2,
            "Most coefficients should be small (near zero)"
        );
    }

    #[test]
    fn test_integers_to_rns_conversion() {
        const DEGREE: usize = 4;
        let basis = Arc::new(
            RnsBasisBuilder::new(DEGREE)
                .with_custom_primes(vec![17, 19, 23])
                .build()
                .unwrap(),
        );

        let coeffs = [100u64, 200, 300, 400];
        let rns_poly = u64_to_rns(&coeffs, basis.clone());

        // Verify each channel has correct residues
        let expected_residues = [
            [100 % 17, 200 % 17, 300 % 17, 400 % 17], // mod 17
            [100 % 19, 200 % 19, 300 % 19, 400 % 19], // mod 19
            [100 % 23, 200 % 23, 300 % 23, 400 % 23], // mod 23
        ];

        for (channel_idx, expected) in expected_residues.iter().enumerate() {
            assert_eq!(
                rns_poly.coefficients[channel_idx], *expected,
                "Channel {} residues don't match",
                channel_idx
            );
        }
    }

    #[test]
    fn test_signed_integers_to_rns_conversion() {
        const DEGREE: usize = 4;
        let basis = Arc::new(
            RnsBasisBuilder::new(DEGREE)
                .with_custom_primes(vec![17, 19])
                .build()
                .unwrap(),
        );

        let coeffs = [1i64, -1, 0, -5];
        let rns_poly: RnsPolyRing<DEGREE> = i64_to_rns(&coeffs, basis.clone());

        // -1 mod 17 = 16, -1 mod 19 = 18
        // -5 mod 17 = 12, -5 mod 19 = 14
        let expected_17 = [1u64, 16, 0, 12];
        let expected_19 = [1u64, 18, 0, 14];

        assert_eq!(rns_poly.coefficients[0], expected_17);
        assert_eq!(rns_poly.coefficients[1], expected_19);
    }

    #[test]
    fn test_sample_uniform_simple() {
        const DEGREE: usize = 8;
        let basis = Arc::new(
            RnsBasisBuilder::new(DEGREE)
                .with_custom_primes(vec![17, 19, 23])
                .build()
                .unwrap(),
        );

        let mut rng = StdRng::seed_from_u64(42);
        let poly: RnsPolyRing<DEGREE> =
            sample_uniform_poly(basis.clone(), &mut rng);

        // Verify it's a valid RNS polynomial
        assert_eq!(poly.channels(), 3);
        assert_eq!(poly.len(), DEGREE);

        // Verify all residues are in range
        for (channel_idx, &prime) in basis.primes().iter().enumerate() {
            for &residue in &poly.coefficients[channel_idx] {
                assert!(residue < prime, "Residue {} >= prime {}", residue, prime);
            }
        }
    }

    #[test]
    fn test_sample_gaussian_simple() {
        const DEGREE: usize = 8;
        let basis = Arc::new(
            RnsBasisBuilder::new(DEGREE)
                .with_custom_primes(vec![17, 19])
                .build()
                .unwrap(),
        );

        let mut rng = StdRng::seed_from_u64(123);
        let poly: RnsPolyRing<DEGREE> =
            sample_gaussian_poly(basis.clone(), 3.2, &mut rng);

        // Verify it's a valid RNS polynomial
        assert_eq!(poly.channels(), 2);
        assert_eq!(poly.len(), DEGREE);

        // Verify all residues are in range
        for (channel_idx, &prime) in basis.primes().iter().enumerate() {
            for &residue in &poly.coefficients[channel_idx] {
                assert!(residue < prime, "Residue {} >= prime {}", residue, prime);
            }
        }
    }

    #[test]
    fn test_reproducible_sampling() {
        const DEGREE: usize = 4;
        let basis = Arc::new(
            RnsBasisBuilder::new(DEGREE)
                .with_custom_primes(vec![17, 19])
                .build()
                .unwrap(),
        );

        // Same seed should produce same results
        let mut rng1 = StdRng::seed_from_u64(999);
        let mut rng2 = StdRng::seed_from_u64(999);

        let poly1: RnsPolyRing<DEGREE> =
            sample_uniform_poly(basis.clone(), &mut rng1);
        let poly2 = sample_uniform_poly(basis.clone(), &mut rng2);

        assert_eq!(poly1.coefficients, poly2.coefficients);
    }
}
