use super::basis::RnsBasis;
use super::math::crt_reconstruct;
use std::sync::Arc;

/// An RNS-encoded polynomial in ℤ[X]/(X^DEGREE + 1) using const generics.
///
/// Coefficients are stored modulus-by-modulus:
/// - Outer `Vec`: one entry per prime (channel) in the RNS basis.
/// - Inner `[u64; DEGREE]`: the residues of all `DEGREE` slots for that prime.
///
/// Arithmetic (NTT, add, mul, rescale) is performed per channel. The polynomial
/// degree (number of slots) is the compile-time constant `DEGREE`.
#[derive(Debug, Clone)]
pub struct RnsPolyRing<const DEGREE: usize> {
    /// `coefficients[c][i]` = i-th slot residue mod basis.primes[c]
    pub coefficients: Vec<[u64; DEGREE]>,
    /// Shared RNS basis with primes and NTT tables.
    pub basis: Arc<RnsBasis>,
}

impl<const DEGREE: usize> RnsPolyRing<DEGREE> {
    /// Create the zero polynomial: all residues = 0.
    pub fn zero(basis: Arc<RnsBasis>) -> Self {
        let channels = basis.channel_count();
        Self {
            coefficients: vec![[0; DEGREE]; channels],
            basis,
        }
    }

    /// Number of slots (polynomial degree = DEGREE).
    pub fn len(&self) -> usize {
        DEGREE
    }

    /// Number of RNS channels (primes).
    pub fn channels(&self) -> usize {
        self.coefficients.len()
    }

    pub fn from_integer_coeffs(ints: &[i64], basis: Arc<RnsBasis>) -> Self {
        assert_eq!(ints.len(), DEGREE, "Input slice length must match DEGREE");
        let mut channels: Vec<[u64; DEGREE]> =
            Vec::with_capacity(basis.primes().len());
        for &q in basis.primes() {
            let mut arr = [0u64; DEGREE];
            let q_i64 = q as i64;
            for i in 0..DEGREE {
                // ensure non-negative residue
                let v = ((ints[i] % q_i64 + q_i64) % q_i64) as u64;
                arr[i] = v;
            }
            channels.push(arr);
        }
        Self {
            coefficients: channels,
            basis,
        }
    }

    pub fn coefficient_to_u64(&self, position: usize) -> u64 {
        assert!(position < DEGREE);

        let residues: Vec<u64> = self
            .coefficients
            .iter()
            .map(|channel| channel[position])
            .collect();

        crt_reconstruct(&residues, self.basis.primes())
    }

    /// Convert entire polynomial to Vec<u64> for debugging
    pub fn to_u64_coefficients(&self) -> Vec<u64> {
        (0..DEGREE).map(|i| self.coefficient_to_u64(i)).collect()
    }
}

#[cfg(test)]
mod tests {
    use crate::RnsBasisBuilder;
    use crate::rings::NttTables;

    use super::*;
    use std::sync::Arc;

    /// Build a minimal RNS basis with `channels` identical primes.
    fn dummy_basis<const D: usize>(channels: usize) -> Arc<RnsBasis> {
        let primes = vec![17u64; channels];
        let ntt_tables = NttTables::build_ntt_tables_for_primes(&primes).unwrap();
        Arc::new(RnsBasis { primes, ntt_tables })
    }

    #[test]
    fn zero_has_correct_dimensions() {
        const D: usize = 4;
        let basis = dummy_basis::<D>(3);
        let poly = RnsPolyRing::<D>::zero(basis.clone());
        assert_eq!(poly.len(), D);
        assert_eq!(poly.channels(), 3);
        // All residues are zero
        for channel in &poly.coefficients {
            for &res in channel.iter() {
                assert_eq!(res, 0, "Expected zero residue");
            }
        }
        // Check shared basis Arc count
        assert_eq!(Arc::strong_count(&basis), 2);
    }

    #[test]
    fn test_rns_roundtrip_conversion() {
        const DEGREE: usize = 8;

        // Create RNS basis with small coprime primes for testing
        let basis = Arc::new(
            RnsBasisBuilder::new(DEGREE)
                .with_custom_primes(vec![17u64, 19u64, 23u64]) // Product = 7429
                .build()
                .expect("Failed to build test RNS basis"),
        );

        // Test various coefficient vectors
        let test_cases = vec![
            // Case 1: Small positive numbers
            vec![1i64, 2, 3, 4, 5, 6, 7, 8],
            // Case 2: Mix of positive and negative (within product range)
            vec![-1i64, 10, -50, 100, -200, 300, -500, 1000],
            // Case 3: Edge cases including zero
            vec![0i64, 1, -1, 7428, -7428, 3714, -3714, 42],
            // Case 4: Larger numbers that should wrap around modulo product
            vec![8000i64, -8000, 15000, -15000, 7429, -7429, 7430, -7430],
        ];

        for (case_idx, original_coeffs) in test_cases.iter().enumerate() {
            println!("Testing case {}: {:?}", case_idx + 1, original_coeffs);

            // Step 1: Create RNS polynomial from integer coefficients
            let rns_poly: RnsPolyRing<8> =
                RnsPolyRing::from_integer_coeffs(original_coeffs, basis.clone());

            // Step 2: Verify that residues are correct for each prime
            let product: u64 = basis.primes().iter().product();
            for (i, &coeff) in original_coeffs.iter().enumerate() {
                for (channel_idx, &prime) in basis.primes().iter().enumerate() {
                    let expected_residue = ((coeff % prime as i64 + prime as i64)
                        % prime as i64)
                        as u64;
                    let actual_residue = rns_poly.coefficients[channel_idx][i];
                    assert_eq!(
                        actual_residue, expected_residue,
                        "Channel {} position {}: expected residue {}, got {}",
                        channel_idx, i, expected_residue, actual_residue
                    );
                }
            }

            // Step 3: Reconstruct back to u64 coefficients
            let reconstructed = rns_poly.to_u64_coefficients();

            // Step 4: Verify reconstruction (accounting for modular reduction)
            for (i, (&original, &reconstructed)) in
                original_coeffs.iter().zip(reconstructed.iter()).enumerate()
            {
                // Convert original to expected value in [0, product) range
                let expected = ((original % product as i64 + product as i64)
                    % product as i64) as u64;

                assert_eq!(
                    reconstructed, expected,
                    "Position {}: original={}, reconstructed={}, expected={}",
                    i, original, reconstructed, expected
                );
            }

            // Step 5: Test individual coefficient reconstruction
            for i in 0..DEGREE {
                let single_coeff = rns_poly.coefficient_to_u64(i);
                assert_eq!(
                    single_coeff, reconstructed[i],
                    "Single coefficient reconstruction mismatch at position {}",
                    i
                );
            }

            println!(
                "✓ Case {} passed: {:?} -> {:?}",
                case_idx + 1,
                original_coeffs,
                reconstructed
            );
        }

        // Additional test: Verify that different coefficient sets produce different RNS representations
        let coeffs1 = vec![1i64, 2, 3, 4, 5, 6, 7, 8];
        let coeffs2 = vec![1i64, 2, 3, 4, 5, 6, 7, 9]; // Only last coefficient differs

        let poly1: RnsPolyRing<8> =
            RnsPolyRing::from_integer_coeffs(&coeffs1, basis.clone());
        let poly2: RnsPolyRing<8> =
            RnsPolyRing::from_integer_coeffs(&coeffs2, basis.clone());

        // They should be different (at least in the last coefficient)
        let last_idx = DEGREE - 1;
        let different_channels =
            basis.primes().iter().enumerate().any(|(channel_idx, _)| {
                poly1.coefficients[channel_idx][last_idx]
                    != poly2.coefficients[channel_idx][last_idx]
            });

        assert!(
            different_channels,
            "Different input coefficients should produce different RNS representations"
        );

        println!("✓ All roundtrip tests passed!");
    }
}
