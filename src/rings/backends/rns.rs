//! RNS (Residue Number System) backend for CKKS polynomial ring operations
//!
//! This module provides a complete RNS implementation that can handle multiple prime
//! moduli simultaneously for improved performance and precision in CKKS operations.

use crate::{PolyRing, PolySampler, math::is_prime};
use rand::Rng;
use rand_distr::{Distribution, Normal};
use std::{
    fmt,
    ops::{Add, AddAssign, MulAssign, Neg},
    sync::Arc,
};
use thiserror::Error;

// ============================================================================
// Error Types
// ============================================================================

/// Errors for RNS basis construction and operations.
#[derive(Error, Debug)]
pub enum RnsError {
    /// Requested number of primes is zero or invalid.
    #[error("Invalid prime count: {0}")]
    InvalidPrimeCount(usize),

    /// Failed to find a suitable prime for the given constraints.
    #[error("Unable to find prime for bit size {0} and ring degree {1}")]
    PrimeGenerationFailed(usize, usize),

    /// Provided basis does not match expected channel count.
    #[error(
        "RNS basis mismatch: expected {expected_count} primes, got {actual_count}"
    )]
    BasisMismatch {
        expected_count: usize,
        actual_count: usize,
    },

    #[error("Invalid prime: {0}")]
    InvalidPrime(u64),
}

/// Arithmetic operation errors
#[derive(Error, Debug)]
pub enum ArithmeticError {
    #[error(
        "Basis mismatch: cannot operate on polynomials with different RNS bases"
    )]
    BasisMismatch,

    #[error("Channel count mismatch: {left} vs {right}")]
    ChannelMismatch { left: usize, right: usize },
}

/// Result type for RNS operations.
pub type RnsResult<T> = Result<T, RnsError>;
pub type ArithmeticResult<T> = Result<T, ArithmeticError>;

/// Generate primes for RNS basis
pub fn generate_primes(
    bit_size: usize,
    count: usize,
    ring_degree: usize,
) -> RnsResult<Vec<u64>> {
    if count == 0 {
        return Err(RnsError::InvalidPrimeCount(count));
    }

    let min_prime = 1u64 << (bit_size - 1);
    let max_prime = (1u64 << bit_size) - 1;
    let mut primes = Vec::with_capacity(count);
    let mut candidate = max_prime;

    // For NTT support, primes should be of the form k * 2n + 1
    // For now, just find primes in the range
    while primes.len() < count && candidate >= min_prime {
        if is_prime(candidate) {
            // Additional check: ensure 2*ring_degree divides (candidate - 1) for NTT
            if (candidate - 1) % (2 * ring_degree as u64) == 0 {
                primes.push(candidate);
            }
        }
        candidate -= 1;
    }

    if primes.len() < count {
        return Err(RnsError::PrimeGenerationFailed(bit_size, ring_degree));
    }

    Ok(primes)
}

/// Chinese Remainder Theorem reconstruction
pub fn crt_reconstruct(residues: &[u64], primes: &[u64]) -> u64 {
    if residues.len() != primes.len() {
        panic!("Residues and primes length mismatch");
    }

    let product: u64 = primes.iter().product();
    let mut result = 0u128;

    for (&residue, &prime) in residues.iter().zip(primes) {
        let quotient = product / prime;
        let inv = mod_inverse(quotient, prime);
        result +=
            (residue as u128 * quotient as u128 * inv as u128) % product as u128;
    }

    (result % product as u128) as u64
}

/// Modular inverse using extended Euclidean algorithm
pub fn mod_inverse(a: u64, m: u64) -> u64 {
    fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
        if a == 0 {
            (b, 0, 1)
        } else {
            let (g, y, x) = extended_gcd(b % a, a);
            (g, x - (b / a) * y, y)
        }
    }

    let (g, x, _) = extended_gcd(a as i64, m as i64);
    if g != 1 {
        panic!("Modular inverse does not exist");
    }
    ((x % m as i64 + m as i64) % m as i64) as u64
}

// ============================================================================
// RNS Basis
// ============================================================================

/// An RNS basis containing prime moduli and NTT tables.
#[derive(Debug)]
pub struct RnsBasis {
    primes: Vec<u64>,
}

impl RnsBasis {
    pub fn primes(&self) -> &Vec<u64> {
        &self.primes
    }

    pub fn channel_count(&self) -> usize {
        self.primes.len()
    }

    pub fn validate(&self, expected: usize) -> RnsResult<()> {
        let actual = self.primes.len();
        if actual != expected {
            Err(RnsError::BasisMismatch {
                expected_count: expected,
                actual_count: actual,
            })
        } else {
            Ok(())
        }
    }
}

/// Builder for constructing an RNS basis.
pub struct RnsBasisBuilder {
    ring_degree: usize,
    prime_bits: Vec<usize>,
    custom_primes: Option<Vec<u64>>,
}

impl RnsBasisBuilder {
    pub fn new(ring_degree: usize) -> Self {
        Self {
            ring_degree,
            prime_bits: Vec::new(),
            custom_primes: None,
        }
    }

    pub fn with_prime_bits(mut self, bits: Vec<usize>) -> Self {
        self.prime_bits = bits;
        self
    }

    pub fn with_custom_primes(mut self, primes: Vec<u64>) -> Self {
        self.custom_primes = Some(primes);
        self
    }

    pub fn with_uniform_prime_bits(mut self, bits: usize, count: usize) -> Self {
        self.prime_bits = vec![bits; count];
        self
    }

    pub fn build(self) -> RnsResult<RnsBasis> {
        let primes = self.get_or_generate_primes()?;

        Ok(RnsBasis { primes })
    }

    fn get_or_generate_primes(&self) -> RnsResult<Vec<u64>> {
        if let Some(ref primes) = self.custom_primes {
            if primes.is_empty() {
                return Err(RnsError::InvalidPrimeCount(0));
            }
            for &p in primes {
                if !is_prime(p) {
                    return Err(RnsError::InvalidPrime(p));
                }
            }
            Ok(primes.clone())
        } else {
            if self.prime_bits.is_empty() {
                return Err(RnsError::InvalidPrimeCount(0));
            }

            let mut all_primes = Vec::new();
            for &bits in &self.prime_bits {
                let mut primes = generate_primes(bits, 1, self.ring_degree)?;
                all_primes.append(&mut primes);
            }
            Ok(all_primes)
        }
    }
}

// ============================================================================
// RNS Polynomial Ring
// ============================================================================

/// RNS-encoded polynomial ring supporting CKKS operations.
#[derive(Debug, Clone)]
pub struct RnsPolyRing<const DEGREE: usize> {
    /// coefficients[c][i] = i-th coefficient residue mod basis.primes[c]
    pub coefficients: Vec<[u64; DEGREE]>,
    /// Shared RNS basis with primes and NTT tables
    pub basis: Arc<RnsBasis>,
}

impl<const DEGREE: usize> RnsPolyRing<DEGREE> {
    /// Create the zero polynomial
    pub fn zero(basis: Arc<RnsBasis>) -> Self {
        let channels = basis.channel_count();
        Self {
            coefficients: vec![[0; DEGREE]; channels],
            basis,
        }
    }

    /// Number of polynomial coefficients
    pub fn len(&self) -> usize {
        DEGREE
    }

    /// Number of RNS channels (primes)
    pub fn channels(&self) -> usize {
        self.coefficients.len()
    }

    /// Create from signed integer slice
    pub fn from_i64_slice(ints: &[i64], basis: Arc<RnsBasis>) -> Self {
        assert_eq!(ints.len(), DEGREE, "Input slice length must match DEGREE");

        let mut channels: Vec<[u64; DEGREE]> =
            Vec::with_capacity(basis.primes().len());

        for &q in basis.primes() {
            let mut arr = [0u64; DEGREE];
            let q_i64 = q as i64;
            for i in 0..DEGREE {
                // Ensure positive residue: ((x % q) + q) % q
                arr[i] = ((ints[i] % q_i64 + q_i64) % q_i64) as u64;
            }
            channels.push(arr);
        }

        Self {
            coefficients: channels,
            basis,
        }
    }

    /// Reconstruct coefficient at given position using CRT
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

    /// Convert RNS coefficients to signed i64 values via CRT reconstruction
    pub fn to_i64_coefficients(&self) -> Vec<i64> {
        let u64s = self.to_u64_coefficients();
        let product: u64 = self.basis.primes.iter().product();
        let half_mod = product / 2;

        u64s.into_iter()
            .map(|x| {
                if x >= half_mod {
                    (x as i128 - product as i128) as i64
                } else {
                    x as i64
                }
            })
            .collect()
    }
}

// ============================================================================
// Trait Implementations for CKKS Integration
// ============================================================================

impl<const DEGREE: usize> PolyRing<DEGREE> for RnsPolyRing<DEGREE> {
    type Context = Arc<RnsBasis>;

    fn zero(context: &Self::Context) -> Self {
        Self::zero(context.clone())
    }

    fn from_coeffs(coeffs: &[i64], context: &Self::Context) -> Self {
        assert!(coeffs.len() >= DEGREE, "Not enough coefficients provided");
        Self::from_i64_slice(&coeffs[..DEGREE], context.clone())
    }

    fn to_coeffs(&self) -> [i64; DEGREE] {
        let i64_coeffs = self.to_i64_coefficients();
        let mut result = [0i64; DEGREE];
        result.copy_from_slice(&i64_coeffs[..DEGREE]);
        result
    }

    fn context(&self) -> &Self::Context {
        &self.basis
    }
}

// ============================================================================
// Arithmetic Operations
// ============================================================================

impl<const DEGREE: usize> AddAssign<&RnsPolyRing<DEGREE>> for RnsPolyRing<DEGREE> {
    fn add_assign(&mut self, rhs: &RnsPolyRing<DEGREE>) {
        // Check basis compatibility
        assert!(
            std::ptr::eq(self.basis.as_ref(), rhs.basis.as_ref()),
            "Cannot add RnsPolyRing with different bases"
        );

        let channel_count = self.channels();
        assert_eq!(channel_count, rhs.channels(), "Channel count mismatch");

        for channel_idx in 0..channel_count {
            let prime = self.basis.primes()[channel_idx];

            for coeff_idx in 0..DEGREE {
                let sum = self.coefficients[channel_idx][coeff_idx]
                    + rhs.coefficients[channel_idx][coeff_idx];
                self.coefficients[channel_idx][coeff_idx] = sum % prime;
            }
        }
    }
}

impl<const DEGREE: usize> MulAssign<&RnsPolyRing<DEGREE>> for RnsPolyRing<DEGREE> {
    fn mul_assign(&mut self, rhs: &RnsPolyRing<DEGREE>) {
        self.mul_assign_checked(rhs)
            .expect("Multiplication failed due to incompatible bases");
    }
}

impl<const DEGREE: usize> Neg for RnsPolyRing<DEGREE> {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        self.negate_assign();
        self
    }
}

impl<const DEGREE: usize> Neg for &RnsPolyRing<DEGREE> {
    type Output = RnsPolyRing<DEGREE>;

    fn neg(self) -> Self::Output {
        self.negate()
    }
}

impl<const DEGREE: usize> Add<&RnsPolyRing<DEGREE>> for &RnsPolyRing<DEGREE> {
    type Output = RnsPolyRing<DEGREE>;

    fn add(self, rhs: &RnsPolyRing<DEGREE>) -> Self::Output {
        assert!(
            std::ptr::eq(self.basis.as_ref(), rhs.basis.as_ref()),
            "Cannot add RnsPolyRing with different bases"
        );

        let mut result = self.clone();
        result += rhs;
        result
    }
}

impl<const DEGREE: usize> RnsPolyRing<DEGREE> {
    /// Checked in-place addition that returns a Result
    pub fn add_assign_checked(
        &mut self,
        rhs: &RnsPolyRing<DEGREE>,
    ) -> ArithmeticResult<()> {
        if !std::ptr::eq(self.basis.as_ref(), rhs.basis.as_ref()) {
            return Err(ArithmeticError::BasisMismatch);
        }

        let left_channels = self.channels();
        let right_channels = rhs.channels();

        if left_channels != right_channels {
            return Err(ArithmeticError::ChannelMismatch {
                left: left_channels,
                right: right_channels,
            });
        }

        for channel_idx in 0..left_channels {
            let prime = self.basis.primes()[channel_idx];

            for coeff_idx in 0..DEGREE {
                let sum = self.coefficients[channel_idx][coeff_idx]
                    + rhs.coefficients[channel_idx][coeff_idx];
                self.coefficients[channel_idx][coeff_idx] = sum % prime;
            }
        }

        Ok(())
    }

    /// In-place polynomial multiplication in quotient ring Z[X]/(X^n + 1)
    pub fn mul_assign_checked(
        &mut self,
        rhs: &RnsPolyRing<DEGREE>,
    ) -> ArithmeticResult<()> {
        if !std::ptr::eq(self.basis.as_ref(), rhs.basis.as_ref()) {
            return Err(ArithmeticError::BasisMismatch);
        }

        let channel_count = self.channels();
        if channel_count != rhs.channels() {
            return Err(ArithmeticError::ChannelMismatch {
                left: channel_count,
                right: rhs.channels(),
            });
        }

        // Perform multiplication channel by channel
        for channel_idx in 0..channel_count {
            let prime = self.basis.primes()[channel_idx] as u128;

            // Store original coefficients
            let lhs_coeffs = self.coefficients[channel_idx];
            let rhs_coeffs = rhs.coefficients[channel_idx];

            // Initialize result to zero
            for coeff in &mut self.coefficients[channel_idx] {
                *coeff = 0;
            }

            // Schoolbook multiplication with quotient reduction
            for i in 0..DEGREE {
                for j in 0..DEGREE {
                    let prod =
                        (lhs_coeffs[i] as u128 * rhs_coeffs[j] as u128) % prime;
                    let degree_sum = i + j;

                    if degree_sum < DEGREE {
                        // Normal coefficient: result[i+j] += a[i] * b[j]
                        let current =
                            self.coefficients[channel_idx][degree_sum] as u128;
                        self.coefficients[channel_idx][degree_sum] =
                            ((current + prod) % prime) as u64;
                    } else {
                        // Quotient reduction: X^(n+k) = -X^k, so subtract from lower degree
                        let wrapped_idx = degree_sum - DEGREE;
                        let current =
                            self.coefficients[channel_idx][wrapped_idx] as u128;
                        self.coefficients[channel_idx][wrapped_idx] =
                            ((current + prime - prod) % prime) as u64;
                    }
                }
            }
        }

        Ok(())
    }

    /// In-place negation
    pub fn negate_assign(&mut self) {
        let channel_count = self.channels();

        for channel_idx in 0..channel_count {
            let prime = self.basis.primes()[channel_idx];

            for coeff_idx in 0..DEGREE {
                let current = self.coefficients[channel_idx][coeff_idx];

                if current == 0 {
                    self.coefficients[channel_idx][coeff_idx] = 0;
                } else {
                    self.coefficients[channel_idx][coeff_idx] = prime - current;
                }
            }
        }
    }

    /// Returns negated polynomial
    pub fn negate(&self) -> Self {
        let mut result = self.clone();
        result.negate_assign();
        result
    }
}

// ============================================================================
// Sampling Operations
// ============================================================================

/// Sample uniformly random integer coefficients
fn sample_uniform_u64<const DEGREE: usize, R: Rng + ?Sized>(
    max_value: u64,
    rng: &mut R,
) -> [u64; DEGREE] {
    let mut coeffs = [0u64; DEGREE];
    for coeff in &mut coeffs {
        *coeff = rng.random_range(0..max_value);
    }
    coeffs
}

/// Sample Gaussian noise as signed integers, then convert to unsigned
fn sample_gaussian_u64<const DEGREE: usize, R: Rng + ?Sized>(
    std_dev: f64,
    max_value: u64,
    rng: &mut R,
) -> [u64; DEGREE] {
    let normal = Normal::new(0.0, std_dev).expect("Invalid Gaussian std_dev");
    let mut coeffs = [0u64; DEGREE];

    for coeff in &mut coeffs {
        let sample = normal.sample(rng);
        let noise_int = sample.round() as i64;

        *coeff = if noise_int < 0 {
            let abs_val = noise_int.unsigned_abs() % max_value;
            if abs_val == 0 { 0 } else { max_value - abs_val }
        } else {
            (noise_int as u64) % max_value
        };
    }

    coeffs
}

/// Sample a ternary polynomial with given Hamming weight
fn sample_ternary_i64<const DEGREE: usize, R: Rng + ?Sized>(
    hamming_weight: usize,
    rng: &mut R,
) -> [i64; DEGREE] {
    let mut out = [0i64; DEGREE];
    let mut indices: Vec<usize> = (0..DEGREE).collect();

    // Shuffle and pick first hamming_weight indices
    for i in 0..hamming_weight.min(DEGREE) {
        let j = rng.random_range(i..DEGREE);
        indices.swap(i, j);
        out[indices[i]] = if rng.random_bool(0.5) { 1 } else { -1 };
    }

    out
}

/// Convert unsigned integer coefficients to RNS representation
fn u64_to_rns<const DEGREE: usize>(
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
fn i64_to_rns<const DEGREE: usize>(
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

impl<const DEGREE: usize> PolySampler<DEGREE> for RnsPolyRing<DEGREE> {
    fn sample_uniform<R: Rng>(rng: &mut R, context: &Self::Context) -> Self {
        let max_prime = *context.primes().iter().max().unwrap_or(&1);
        let coeffs = sample_uniform_u64::<DEGREE, _>(max_prime, rng);
        u64_to_rns(&coeffs, context.clone())
    }

    fn sample_gaussian<R: Rng>(
        rng: &mut R,
        std_dev: f64,
        context: &Self::Context,
    ) -> Self {
        let modulus_product: u64 = context.primes().iter().product();
        let coeffs =
            sample_gaussian_u64::<DEGREE, _>(std_dev, modulus_product, rng);
        u64_to_rns(&coeffs, context.clone())
    }

    fn sample_tribits<R: Rng>(
        rng: &mut R,
        hamming_weight: usize,
        context: &Self::Context,
    ) -> Self {
        let ternary = sample_ternary_i64::<DEGREE, _>(hamming_weight, rng);
        i64_to_rns(&ternary, context.clone())
    }

    fn sample_noise<R: Rng>(
        rng: &mut R,
        variance: f64,
        context: &Self::Context,
    ) -> Self {
        Self::sample_gaussian(rng, variance.sqrt(), context)
    }
}

// ============================================================================
// Display Implementation
// ============================================================================

impl<const DEGREE: usize> fmt::Display for RnsPolyRing<DEGREE> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let deg = DEGREE;
        let channels = self.basis.channel_count();
        let num = 3; // Number of terms to show at start/end

        write!(f, "RnsPoly<{deg}>[")?;

        if deg <= num * 2 {
            // Show all coefficients
            for i in 0..deg {
                if i > 0 {
                    write!(f, ", ")?;
                }
                write!(f, "[")?;
                for c in 0..channels {
                    if c > 0 {
                        write!(f, ", ")?;
                    }
                    write!(f, "{}", self.coefficients[c][i])?;
                }
                write!(f, "]")?;
                if i > 0 {
                    write!(f, "*x^{i}")?;
                }
            }
        } else {
            // Show first `num` terms
            for i in 0..num {
                if i > 0 {
                    write!(f, ", ")?;
                }
                write!(f, "[")?;
                for c in 0..channels {
                    if c > 0 {
                        write!(f, ", ")?;
                    }
                    write!(f, "{}", self.coefficients[c][i])?;
                }
                write!(f, "]*x^{i}")?;
            }
            write!(f, ", …")?;
            // Show last `num` terms
            for i in (deg - num)..deg {
                write!(f, ", [")?;
                for c in 0..channels {
                    if c > 0 {
                        write!(f, ", ")?;
                    }
                    write!(f, "{}", self.coefficients[c][i])?;
                }
                write!(f, "]*x^{i}")?;
            }
        }

        write!(f, "]")
    }
}
