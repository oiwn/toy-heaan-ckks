//! NTT-based RNS polynomial ring implementation for CKKS
//!
//! This module provides an RNS polynomial ring that tracks domain state
//! and supports efficient polynomial operations through NTT transforms.

use crate::math::{
    gaussian_coefficients, ternary_coefficients, uniform_coefficients,
};
use crate::{PolyRing, PolySampler};
use rand::Rng;
use std::{
    fmt,
    ops::{Add, AddAssign, MulAssign, Neg},
    sync::Arc,
};
use thiserror::Error;

// Re-use existing RNS infrastructure
use super::rns::{RnsBasis, crt_reconstruct, mod_inverse};

// ============================================================================
// NTT Error Types
// ============================================================================

#[derive(Error, Debug)]
pub enum NttError {
    #[error("No primitive root of order {order} exists for prime {prime}")]
    NoPrimitiveRoot { prime: u64, order: usize },

    #[error("Invalid degree {degree}: must be power of 2")]
    InvalidDegree { degree: usize },

    #[error("Prime {prime} not NTT-friendly for degree {degree}")]
    PrimeNotNttFriendly { prime: u64, degree: usize },

    #[error("RNS basis error: {0}")]
    RnsBasisError(String),
}

pub type NttResult<T> = Result<T, NttError>;

// ============================================================================
// NTT Tables
// ============================================================================

/// Precomputed NTT tables for a specific prime and ring degree
#[derive(Debug, Clone)]
pub struct NttTables {
    /// Forward NTT roots: powers of primitive root
    pub forward_roots: Vec<u64>,
    /// Inverse NTT roots: powers of inverse primitive root
    pub inverse_roots: Vec<u64>,
    /// Inverse of n modulo prime (for normalization)
    pub n_inv: u64,
    /// The prime modulus
    pub prime: u64,
}

impl NttTables {
    /// Generate NTT tables for given prime and degree
    /// Requires: prime ≡ 1 (mod 2*degree) for primitive root to exist
    pub fn new(prime: u64, degree: usize) -> NttResult<Self> {
        if !degree.is_power_of_two() {
            return Err(NttError::InvalidDegree { degree });
        }

        // Check if prime is NTT-friendly: (prime-1) % (2*degree) == 0
        if (prime - 1) % (2 * degree as u64) != 0 {
            return Err(NttError::PrimeNotNttFriendly { prime, degree });
        }

        // Find primitive root of order 2*degree
        let primitive_root = find_primitive_root(prime, 2 * degree)?;

        // Precompute forward roots: ω^i for i = 0..degree
        let mut forward_roots = Vec::with_capacity(degree);
        let mut power = 1u64;
        for _ in 0..degree {
            forward_roots.push(power);
            power = (power * primitive_root) % prime;
        }

        // Precompute inverse roots
        let inv_root = mod_inverse(primitive_root, prime);
        let mut inverse_roots = Vec::with_capacity(degree);
        power = 1u64;
        for _ in 0..degree {
            inverse_roots.push(power);
            power = (power * inv_root) % prime;
        }

        // Inverse of degree for normalization
        let n_inv = mod_inverse(degree as u64, prime);

        Ok(NttTables {
            forward_roots,
            inverse_roots,
            n_inv,
            prime,
        })
    }
}

// ============================================================================
// Extended RNS Basis with NTT Support
// ============================================================================

/// Extended RNS basis that includes NTT tables for each prime
#[derive(Debug)]
pub struct RnsNttBasis {
    /// Base RNS functionality
    pub base: RnsBasis,
    /// NTT tables for each prime
    pub ntt_tables: Vec<NttTables>,
}

impl RnsNttBasis {
    /// Create NTT-enabled RNS basis from existing RnsBasis
    pub fn new(base: RnsBasis, ring_degree: usize) -> NttResult<Self> {
        let mut ntt_tables = Vec::with_capacity(base.channel_count());

        for &prime in base.primes() {
            let table = NttTables::new(prime, ring_degree)?;
            ntt_tables.push(table);
        }

        Ok(RnsNttBasis { base, ntt_tables })
    }

    pub fn primes(&self) -> &Vec<u64> {
        self.base.primes()
    }

    pub fn channel_count(&self) -> usize {
        self.base.channel_count()
    }
}

// ============================================================================
// RNS NTT Polynomial Ring
// ============================================================================

/// RNS polynomial ring with NTT domain tracking
#[derive(Debug, Clone)]
pub struct RnsNttPolyRing<const DEGREE: usize> {
    /// coefficients[c][i] = i-th coefficient residue mod basis.primes[c]
    pub coefficients: Vec<[u64; DEGREE]>,
    /// Shared NTT-enabled RNS basis
    pub basis: Arc<RnsNttBasis>,
    /// Current domain: true = NTT form, false = coefficient form
    pub ntt_form: bool,
}

impl<const DEGREE: usize> RnsNttPolyRing<DEGREE> {
    /// Create the zero polynomial in coefficient form
    pub fn zero(basis: Arc<RnsNttBasis>) -> Self {
        let channels = basis.channel_count();
        Self {
            coefficients: vec![[0; DEGREE]; channels],
            basis,
            ntt_form: false,
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

    /// Create from signed integer slice (always in coefficient form)
    pub fn from_i64_slice(ints: &[i64], basis: Arc<RnsNttBasis>) -> Self {
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
            ntt_form: false,
        }
    }

    /// Convert to NTT form (in-place)
    pub fn to_ntt_form(&mut self) {
        if self.ntt_form {
            return; // Already in NTT form
        }

        for (channel_idx, channel) in self.coefficients.iter_mut().enumerate() {
            let table = &self.basis.ntt_tables[channel_idx];
            forward_ntt(channel, table);
        }

        self.ntt_form = true;
    }

    /// Convert to coefficient form (in-place)
    pub fn to_coeff_form(&mut self) {
        if !self.ntt_form {
            return; // Already in coefficient form
        }

        for (channel_idx, channel) in self.coefficients.iter_mut().enumerate() {
            let table = &self.basis.ntt_tables[channel_idx];
            inverse_ntt(channel, table);
        }

        self.ntt_form = false;
    }

    /// Check if polynomial is in NTT form
    pub fn is_ntt_form(&self) -> bool {
        self.ntt_form
    }

    /// Reconstruct coefficient at given position using CRT (requires coefficient form)
    pub fn coefficient_to_u64(&self, position: usize) -> u64 {
        assert!(position < DEGREE);
        assert!(
            !self.ntt_form,
            "Must be in coefficient form for CRT reconstruction"
        );

        let residues: Vec<u64> = self
            .coefficients
            .iter()
            .map(|channel| channel[position])
            .collect();

        crt_reconstruct(&residues, self.basis.primes())
    }

    /// Convert entire polynomial to Vec<u64> for debugging (requires coefficient form)
    pub fn to_u64_coefficients(&self) -> Vec<u64> {
        assert!(
            !self.ntt_form,
            "Must be in coefficient form for reconstruction"
        );
        (0..DEGREE).map(|i| self.coefficient_to_u64(i)).collect()
    }

    /// Convert RNS coefficients to signed i64 values via CRT reconstruction
    pub fn to_i64_coefficients(&self) -> Vec<i64> {
        let u64s = self.to_u64_coefficients();
        let product: u64 = self.basis.primes().iter().product();
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

    /// In-place negation
    pub fn negate_assign(&mut self) {
        for (channel_idx, channel) in self.coefficients.iter_mut().enumerate() {
            let prime = self.basis.primes()[channel_idx];
            for coeff in channel.iter_mut() {
                *coeff = if *coeff == 0 { 0 } else { prime - *coeff };
            }
        }
    }

    /// Create negated copy
    pub fn negate(&self) -> Self {
        let mut result = self.clone();
        result.negate_assign();
        result
    }
}

// ============================================================================
// Trait Implementations for CKKS Integration
// ============================================================================

impl<const DEGREE: usize> PolyRing<DEGREE> for RnsNttPolyRing<DEGREE> {
    type Context = Arc<RnsNttBasis>;

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
// Arithmetic Operations (work in both domains)
// ============================================================================

impl<const DEGREE: usize> AddAssign<&RnsNttPolyRing<DEGREE>>
    for RnsNttPolyRing<DEGREE>
{
    fn add_assign(&mut self, rhs: &RnsNttPolyRing<DEGREE>) {
        // Check basis compatibility
        assert!(
            std::ptr::eq(self.basis.as_ref(), rhs.basis.as_ref()),
            "Cannot add RnsNttPolyRing with different bases"
        );

        // Both must be in same domain
        assert_eq!(
            self.ntt_form, rhs.ntt_form,
            "Both polynomials must be in the same domain for addition"
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

impl<const DEGREE: usize> MulAssign<&RnsNttPolyRing<DEGREE>>
    for RnsNttPolyRing<DEGREE>
{
    fn mul_assign(&mut self, rhs: &RnsNttPolyRing<DEGREE>) {
        // Check basis compatibility
        assert!(
            std::ptr::eq(self.basis.as_ref(), rhs.basis.as_ref()),
            "Cannot multiply RnsNttPolyRing with different bases"
        );

        // For now, use naive multiplication (TODO: use NTT multiplication)
        let channel_count = self.channels();
        assert_eq!(channel_count, rhs.channels(), "Channel count mismatch");

        // Ensure both are in coefficient form for naive multiplication
        let mut lhs_copy = self.clone();
        let mut rhs_copy = rhs.clone();
        lhs_copy.to_coeff_form();
        rhs_copy.to_coeff_form();

        // Perform multiplication channel by channel
        for channel_idx in 0..channel_count {
            let prime = self.basis.primes()[channel_idx] as u128;

            // Store original coefficients
            let lhs_coeffs = lhs_copy.coefficients[channel_idx];
            let rhs_coeffs = rhs_copy.coefficients[channel_idx];

            // Compute polynomial multiplication in Z[X]/(X^DEGREE + 1)
            let mut result = [0u64; DEGREE];

            for i in 0..DEGREE {
                for j in 0..DEGREE {
                    let coeff_product =
                        (lhs_coeffs[i] as u128 * rhs_coeffs[j] as u128) % prime;

                    if i + j < DEGREE {
                        // Normal term
                        result[i + j] = ((result[i + j] as u128 + coeff_product)
                            % prime) as u64;
                    } else {
                        // Wrap around with negation: X^n = -1
                        let wrapped_idx = (i + j) - DEGREE;
                        result[wrapped_idx] =
                            ((result[wrapped_idx] as u128 + prime - coeff_product)
                                % prime) as u64;
                    }
                }
            }

            self.coefficients[channel_idx] = result;
        }

        // Result is in coefficient form
        self.ntt_form = false;
    }
}

impl<const DEGREE: usize> Neg for RnsNttPolyRing<DEGREE> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let mut result = self;
        result.negate_assign();
        result
    }
}

impl<const DEGREE: usize> Add<&RnsNttPolyRing<DEGREE>> for &RnsNttPolyRing<DEGREE> {
    type Output = RnsNttPolyRing<DEGREE>;

    fn add(self, rhs: &RnsNttPolyRing<DEGREE>) -> Self::Output {
        let mut result = self.clone();
        result += rhs;
        result
    }
}

// ============================================================================
// Sampling Implementation
// ============================================================================

impl<const DEGREE: usize> PolySampler<DEGREE> for RnsNttPolyRing<DEGREE> {
    fn sample_uniform<R: Rng>(rng: &mut R, context: &Self::Context) -> Self {
        // Find the largest prime in the basis to use as max_value
        let max_prime = *context.primes().iter().max().unwrap_or(&1);

        // Sample uniform integers and convert to RNS
        let coeffs = uniform_coefficients::<DEGREE, _>(max_prime, rng);

        // Convert u64 array to i64 for from_i64_slice
        let i64_coeffs: [i64; DEGREE] = coeffs.map(|x| x as i64);
        Self::from_i64_slice(&i64_coeffs, context.clone())
    }

    fn sample_gaussian<R: Rng>(
        rng: &mut R,
        std_dev: f64,
        context: &Self::Context,
    ) -> Self {
        // Use the product of all primes as the working modulus for Gaussian sampling
        let modulus_product: u64 = context.primes().iter().product();

        // Sample Gaussian integers
        let coeffs =
            gaussian_coefficients::<DEGREE, _>(std_dev, modulus_product, rng);

        // Convert u64 array to i64 for from_i64_slice
        let i64_coeffs: [i64; DEGREE] = coeffs.map(|x| x as i64);
        Self::from_i64_slice(&i64_coeffs, context.clone())
    }

    fn sample_tribits<R: Rng>(
        rng: &mut R,
        hamming_weight: usize,
        context: &Self::Context,
    ) -> Self {
        let ternary = ternary_coefficients::<DEGREE, _>(hamming_weight, rng);
        Self::from_i64_slice(&ternary, context.clone())
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
// NTT Algorithms (placeholder for now)
// ============================================================================

/// Forward NTT: coefficient domain -> evaluation domain
pub fn forward_ntt(coeffs: &mut [u64], tables: &NttTables) {
    let n = coeffs.len();
    assert!(n.is_power_of_two(), "NTT size must be power of 2");
    assert_eq!(n, tables.forward_roots.len(), "Table size mismatch");

    // Bit-reverse permutation
    bit_reverse_permute(coeffs);

    // Cooley-Tukey NTT algorithm
    cooley_tukey_ntt(coeffs, &tables.forward_roots, tables.prime);
}

/// Inverse NTT: evaluation domain -> coefficient domain
pub fn inverse_ntt(evals: &mut [u64], tables: &NttTables) {
    let n = evals.len();
    assert!(n.is_power_of_two(), "NTT size must be power of 2");
    assert_eq!(n, tables.inverse_roots.len(), "Table size mismatch");

    // Bit-reverse permutation
    bit_reverse_permute(evals);

    // Gentleman-Sande NTT algorithm (for inverse)
    gentleman_sande_ntt(evals, &tables.inverse_roots, tables.prime);

    // Normalize by 1/n
    for coeff in evals.iter_mut() {
        *coeff = (*coeff * tables.n_inv) % tables.prime;
    }
}

// ============================================================================
// Helper Functions
// ============================================================================

/// Find primitive root of given order modulo prime (simplified version)
fn find_primitive_root(prime: u64, order: usize) -> NttResult<u64> {
    // Check if primitive root can exist
    if (prime - 1) % (order as u64) != 0 {
        return Err(NttError::NoPrimitiveRoot { prime, order });
    }

    // Simple brute force search (can be optimized)
    for candidate in 2..prime {
        if is_primitive_root(candidate, prime, order) {
            return Ok(candidate);
        }
    }

    Err(NttError::NoPrimitiveRoot { prime, order })
}

/// Check if g is a primitive root of order `order` modulo `prime`
fn is_primitive_root(g: u64, prime: u64, order: usize) -> bool {
    // Check if g^order ≡ 1 (mod prime)
    if mod_exp(g, order as u64, prime) != 1 {
        return false;
    }

    // Check that g^(order/p) ≢ 1 for all prime divisors p of order
    // For simplicity, just check a few common cases
    if order % 2 == 0 && mod_exp(g, (order / 2) as u64, prime) == 1 {
        return false;
    }

    true
}

/// Modular exponentiation: base^exp mod m
fn mod_exp(mut base: u64, mut exp: u64, m: u64) -> u64 {
    let mut result = 1;
    base %= m;

    while exp > 0 {
        if exp & 1 == 1 {
            result = (result * base) % m;
        }
        exp >>= 1;
        base = (base * base) % m;
    }

    result
}

/// Bit-reverse permutation for NTT input ordering
fn bit_reverse_permute(data: &mut [u64]) {
    let n = data.len();
    let log_n = n.trailing_zeros() as usize;

    for i in 0..n {
        let j = bit_reverse(i, log_n);
        if i < j {
            data.swap(i, j);
        }
    }
}

/// Reverse bits in integer (for bit-reverse permutation)
fn bit_reverse(mut x: usize, log_n: usize) -> usize {
    let mut result = 0;
    for _ in 0..log_n {
        result = (result << 1) | (x & 1);
        x >>= 1;
    }
    result
}

/// Cooley-Tukey NTT algorithm (decimation-in-time)
fn cooley_tukey_ntt(data: &mut [u64], roots: &[u64], prime: u64) {
    let n = data.len();
    let mut len = 2;

    while len <= n {
        let step = n / len;
        for i in (0..n).step_by(len) {
            for j in 0..len / 2 {
                let u = data[i + j];
                let v = (data[i + j + len / 2] * roots[step * j]) % prime;
                data[i + j] = (u + v) % prime;
                data[i + j + len / 2] = (u + prime - v) % prime;
            }
        }
        len <<= 1;
    }
}

/// Gentleman-Sande NTT algorithm (decimation-in-frequency, for inverse)
fn gentleman_sande_ntt(data: &mut [u64], inv_roots: &[u64], prime: u64) {
    let n = data.len();
    let mut len = n;

    while len >= 2 {
        let step = n / len;
        for i in (0..n).step_by(len) {
            for j in 0..len / 2 {
                let u = data[i + j];
                let v = data[i + j + len / 2];
                data[i + j] = (u + v) % prime;
                data[i + j + len / 2] =
                    ((u + prime - v) * inv_roots[step * j]) % prime;
            }
        }
        len >>= 1;
    }
}

// ============================================================================
// Builder Integration
// ============================================================================

/// Extend CkksEngineBuilder to support NTT backend
pub trait NttEngineBuilder<const DEGREE: usize> {
    fn build_rns_ntt(
        self,
        ntt_basis: Arc<RnsNttBasis>,
        scale_bits: u32,
    ) -> Result<crate::CkksEngine<RnsNttPolyRing<DEGREE>, DEGREE>, crate::CkksError>;
}

// Note: This would be implemented on CkksEngineBuilder in the main codebase
// For now, we'll add this as a utility function for the example

/// Create CKKS engine with NTT backend (helper function)
pub fn create_ntt_engine<const DEGREE: usize>(
    ntt_basis: Arc<RnsNttBasis>,
    scale_bits: u32,
    error_variance: f64,
    hamming_weight: usize,
) -> crate::CkksEngine<RnsNttPolyRing<DEGREE>, DEGREE> {
    let params = crate::crypto::engine::CkksParams {
        error_variance,
        hamming_weight,
        scale_bits,
    };

    crate::CkksEngine::new(ntt_basis, params)
}

// ============================================================================
// Display Implementation
// ============================================================================

impl<const DEGREE: usize> fmt::Display for RnsNttPolyRing<DEGREE> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let domain = if self.ntt_form { "NTT" } else { "COEFF" };
        write!(f, "RnsNttPolyRing[{}][{}] = [", DEGREE, domain)?;

        // Show coefficient reconstruction only if in coefficient form and small degree
        if !self.ntt_form && DEGREE <= 8 {
            let coeffs = self.to_u64_coefficients();
            for (i, &coeff) in coeffs.iter().enumerate() {
                if i > 0 {
                    write!(f, ", ")?;
                }
                write!(f, "{}", coeff)?;
            }
        } else {
            write!(f, "..{} channels..", self.channels())?;
        }

        write!(f, "]")
    }
}
