use crate::{
    Ciphertext, PolyRescale, PolyRing, PolySampler, math::ternary_coefficients,
};
use crypto_bigint::{NonZero, U256, Zero};
use rand::Rng;
use rand_distr::{Distribution, Normal};
use std::ops::{AddAssign, MulAssign, Neg};

pub use super::display::CompactPolyDisplay;

/// BigInt context that supports Kim's HEAAN modulus requirements
///
/// Kim's HEAAN uses multiple moduli for different operations:
/// - q: Base modulus (2^logq) for regular ciphertext operations
/// - Q: Extended modulus (2^logQ) for key switching operations
/// - qQ: Product modulus (q * Q) for relinearization
/// - QQ: Extended key modulus (Q^2) for key generation
#[derive(Debug, Clone, PartialEq)]
pub struct BigIntContext {
    pub q: NonZero<U256>, // Base ciphertext modulus (2^logq)
    pub logq: u32,        // Base modulus bits
    pub log_q: u32,       // Extended modulus bits (logQ for key switching)
}

impl BigIntContext {
    pub fn new(q: NonZero<U256>, logq: u32, log_q: u32) -> Self {
        Self { q, logq, log_q }
    }

    /// Create from Kim's standard parameters: q ≈ 2^logq (closest odd), Q = 2^logQ
    pub fn from_kim_params(logq: u32, log_q: u32) -> Result<Self, String> {
        if logq >= 256 || log_q >= 256 {
            return Err("Modulus bits too large for U256".to_string());
        }

        // Use 2^logq - 1 to ensure odd modulus (safer for cryptographic operations)
        let q_val = (U256::ONE << logq) - U256::ONE;
        println!("Debug: logq={}, q_val={:?} (2^{} - 1)", logq, q_val, logq);

        if q_val.is_zero().into() {
            return Err(format!("Base modulus q is zero for 2^{} - 1", logq));
        }

        let q = NonZero::new(q_val).expect("q_val should not be zero");

        Ok(Self { q, logq, log_q })
    }

    /// Extended modulus Q = 2^logQ for key switching
    pub fn q_extended(&self) -> NonZero<U256> {
        let q_val = U256::ONE << self.log_q;
        NonZero::new(q_val).expect("Extended Q modulus should not be zero")
    }

    /// Get Q as U256 (for multiplication operations)
    pub fn q_extended_u256(&self) -> U256 {
        U256::ONE << self.log_q
    }

    /// Combined modulus qQ = q * Q for relinearization key operations
    /// This is used in Step 3 of Kim's multiplication algorithm
    pub fn q_times_q(&self) -> U256 {
        let q = self.q.get();
        let q_ext = self.q_extended_u256();
        // For U256, we need to be careful about overflow
        // In practice, we might need to use a larger type or handle this differently
        q.saturating_mul(&q_ext)
    }

    /// Extended key modulus QQ = Q^2 for key generation
    /// Used in Kim's key generation (addMultKey, addConjKey, etc.)
    pub fn qq_modulus(&self) -> U256 {
        let q_ext = self.q_extended_u256();
        q_ext.saturating_mul(&q_ext)
    }

    /// Right shift by logQ bits (key switching scaling operation)
    /// This is the critical Step 4 in Kim's multiplication algorithm
    pub fn right_shift_by_log_q(&self, value: U256) -> U256 {
        value >> self.log_q
    }

    /// Left shift by logQ bits (for key generation scaling)
    pub fn left_shift_by_log_q(&self, value: U256) -> U256 {
        value << self.log_q
    }

    /// Get modulus for specific domain
    #[allow(dead_code)]
    pub fn get_modulus(&self, domain: ModulusDomain) -> NonZero<U256> {
        match domain {
            ModulusDomain::Base => self.q,
            ModulusDomain::Extended => self.q_extended(),
            ModulusDomain::Combined => NonZero::new(self.q_times_q())
                .expect("Combined modulus qQ should not be zero"),
            ModulusDomain::ExtendedKey => NonZero::new(self.qq_modulus())
                .expect("Extended key modulus QQ should not be zero"),
        }
    }
}

/// Domain types for different modulus operations in Kim's HEAAN
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ModulusDomain {
    Base,        // q = 2^logq (regular operations)
    Extended,    // Q = 2^logQ (key switching)
    Combined,    // qQ = q * Q (relinearization)
    ExtendedKey, // QQ = Q^2 (key generation)
}

#[derive(Debug, Clone, PartialEq)]
pub struct BigIntPolyRing<const DEGREE: usize> {
    pub coeffs: [U256; DEGREE], // Heap-allocated for large polynomials
    context: BigIntContext,
}

impl<const DEGREE: usize> BigIntPolyRing<DEGREE> {
    pub fn with_context(context: BigIntContext) -> Self {
        Self {
            coeffs: [U256::ZERO; DEGREE],
            context,
        }
    }

    /// Backward compatibility - create with simple modulus
    pub fn with_modulus(modulus: NonZero<U256>) -> Self {
        // Default values: logq = 50, logQ = 20
        let context = BigIntContext::new(modulus, 50, 20);
        Self::with_context(context)
    }

    pub fn coefficients(&self) -> [U256; DEGREE] {
        self.coeffs
    }

    pub fn modulus(&self) -> NonZero<U256> {
        self.context.q
    }

    pub fn big_int_context(&self) -> &BigIntContext {
        &self.context
    }

    /// Create polynomial from U256 coefficients directly
    pub fn from_u256_coeffs(coeffs: &[U256], context: &BigIntContext) -> Self {
        let mut result_coeffs = [U256::ZERO; DEGREE];
        for (i, &coeff) in coeffs.iter().enumerate().take(DEGREE) {
            result_coeffs[i] = coeff % context.q.get(); // Ensure coefficient is in range
        }

        Self {
            coeffs: result_coeffs,
            context: context.clone(),
        }
    }

    // ============================================================================
    // Kim's HEAAN Modulus Switching Operations
    // ============================================================================

    /// Switch to extended modulus domain (qQ) for key switching operations
    /// This is used in Step 3 of Kim's multiplication algorithm
    pub fn to_extended_domain(&self) -> Self {
        let q_q_modulus = self.context.q_times_q();
        let mut extended_coeffs = [U256::ZERO; DEGREE];

        for (i, &coeff) in self.coeffs.iter().enumerate() {
            // Convert coefficient to extended domain by reducing mod qQ
            extended_coeffs[i] = coeff % q_q_modulus;
        }

        Self {
            coeffs: extended_coeffs,
            context: self.context.clone(),
        }
    }

    /// Scale down by Q after key switching (Step 4 of Kim's multiplication)
    /// This is the critical missing operation identified in BIGINT_MUL.md
    pub fn scale_down_by_q(&mut self) {
        for coeff in self.coeffs.iter_mut() {
            *coeff = self.context.right_shift_by_log_q(*coeff);
        }
    }

    /// Scale down by Q and return new polynomial (immutable version)
    pub fn scaled_down_by_q(&self) -> Self {
        let mut scaled_coeffs = [U256::ZERO; DEGREE];

        for (i, &coeff) in self.coeffs.iter().enumerate() {
            scaled_coeffs[i] = self.context.right_shift_by_log_q(coeff);
        }

        Self {
            coeffs: scaled_coeffs,
            context: self.context.clone(),
        }
    }

    /// Scale up by Q for key generation (left shift by logQ bits)
    pub fn scale_up_by_q(&mut self) {
        for coeff in self.coeffs.iter_mut() {
            *coeff = self.context.left_shift_by_log_q(*coeff);
        }
    }

    /// Reduce coefficients modulo specific domain
    pub fn reduce_mod(&mut self, domain: ModulusDomain) {
        let modulus = self.context.get_modulus(domain).get();

        for coeff in self.coeffs.iter_mut() {
            *coeff = *coeff % modulus;
        }
    }

    /// Create copy with coefficients reduced modulo specific domain
    pub fn with_mod_reduction(&self, domain: ModulusDomain) -> Self {
        let modulus = self.context.get_modulus(domain).get();
        let mut reduced_coeffs = [U256::ZERO; DEGREE];

        for (i, &coeff) in self.coeffs.iter().enumerate() {
            reduced_coeffs[i] = coeff % modulus;
        }

        Self {
            coeffs: reduced_coeffs,
            context: self.context.clone(),
        }
    }

    /// Multiply in extended domain (qQ) - used for key switching
    pub fn mul_in_extended_domain(&mut self, other: &Self) {
        let q_q_modulus = self.context.q_times_q();

        for (i, coeff) in self.coeffs.iter_mut().enumerate() {
            let product = coeff.saturating_mul(&other.coeffs[i]);
            *coeff = product % q_q_modulus;
        }
    }

    /// Create polynomial operating in specific modulus domain
    pub fn in_domain(&self, domain: ModulusDomain) -> Self {
        self.with_mod_reduction(domain)
    }
}

impl<const DEGREE: usize> PolyRing<DEGREE> for BigIntPolyRing<DEGREE> {
    type Context = BigIntContext;

    fn zero(context: &Self::Context) -> Self {
        Self {
            coeffs: [U256::ZERO; DEGREE],
            context: context.clone(),
        }
    }

    fn from_coeffs(coeffs: &[i64], context: &Self::Context) -> Self {
        assert!(
            coeffs.len() >= DEGREE,
            "Insufficient number of coefficients"
        );

        let mut poly_coeffs = [U256::ZERO; DEGREE];
        for (i, &coeff) in coeffs.iter().take(DEGREE).enumerate() {
            poly_coeffs[i] = if coeff >= 0 {
                U256::from(coeff as u64)
            } else {
                // For negative coefficients: modulus - |coeff|
                let abs_coeff = U256::from((-coeff) as u64);
                context.q.wrapping_sub(&abs_coeff)
            };
        }

        Self {
            coeffs: poly_coeffs,
            context: context.clone(),
        }
    }

    fn context(&self) -> &Self::Context {
        &self.context
    }

    fn to_coeffs(&self) -> [i64; DEGREE] {
        let mut signed_coeffs = [0i64; DEGREE];
        let half = self.context.q.wrapping_shr(1); // modulus / 2

        for (i, &c) in self.coeffs.iter().enumerate() {
            signed_coeffs[i] = if c < half {
                // Convert to u64 first, then to i64 (assuming it fits)
                // For toy implementation, we assume coefficients fit in i64
                let as_u64 = c.as_words()[0]; // Get lowest 64 bits
                as_u64 as i64
            } else {
                // Negative representation: -(modulus - c)
                let diff = self.context.q.wrapping_sub(&c);
                let as_u64 = diff.as_words()[0];
                -(as_u64 as i64)
            };
        }

        signed_coeffs
    }
}

impl<const DEGREE: usize> AddAssign<&Self> for BigIntPolyRing<DEGREE> {
    fn add_assign(&mut self, rhs: &Self) {
        assert_eq!(
            self.context.q, rhs.context.q,
            "Cannot add polynomials with different moduli"
        );

        for i in 0..DEGREE {
            // Use crypto-bigint modular addition
            self.coeffs[i] =
                self.coeffs[i].add_mod(&rhs.coeffs[i], &self.context.q);
        }
    }
}

impl<const DEGREE: usize> MulAssign<&Self> for BigIntPolyRing<DEGREE> {
    fn mul_assign(&mut self, rhs: &Self) {
        assert_eq!(
            self.context.q, rhs.context.q,
            "Cannot multiply polynomials with different moduli"
        );

        // Schoolbook polynomial multiplication in Z[X]/(X^DEGREE + 1)
        let mut result = [U256::ZERO; DEGREE];

        for i in 0..DEGREE {
            for j in 0..DEGREE {
                let coeff_product =
                    self.coeffs[i].mul_mod(&rhs.coeffs[j], &self.context.q);

                if i + j < DEGREE {
                    // Normal coefficient
                    result[i + j] =
                        result[i + j].add_mod(&coeff_product, &self.context.q);
                } else {
                    // Wrap around with negation due to X^DEGREE = -1
                    let wrapped_idx = (i + j) - DEGREE;
                    result[wrapped_idx] = result[wrapped_idx]
                        .sub_mod(&coeff_product, &self.context.q);
                }
            }
        }

        self.coeffs = result;
    }
}

impl<const DEGREE: usize> Neg for BigIntPolyRing<DEGREE> {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        for coeff in &mut self.coeffs {
            *coeff = self.context.q.wrapping_sub(coeff);
        }
        self
    }
}

impl<const DEGREE: usize> PolySampler<DEGREE> for BigIntPolyRing<DEGREE> {
    fn sample_uniform<R: Rng>(context: &Self::Context, rng: &mut R) -> Self {
        // Generate uniform coefficients using full U256 range
        let mut coeffs = [U256::ZERO; DEGREE];

        for coeff in &mut coeffs {
            // Build a full U256 from four u64 values for proper uniform distribution
            let words = [
                rng.random::<u64>(),
                rng.random::<u64>(),
                rng.random::<u64>(),
                rng.random::<u64>(),
            ];
            let random_u256 = U256::from_words(words);
            *coeff = random_u256.rem(&context.q);
        }

        Self {
            coeffs,
            context: context.clone(),
        }
    }

    fn sample_gaussian<R: Rng>(
        std_dev: f64,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self {
        // Sample centered Gaussian integers, map into Z_q using full U256 arithmetic
        let q = context.q.get();
        let normal = Normal::new(0.0, std_dev).expect("invalid Gaussian std_dev");
        let mut coeffs = [U256::ZERO; DEGREE];

        for c in &mut coeffs {
            let sample = normal.sample(rng).round() as i128; // centered integer
            if sample == 0 {
                *c = U256::ZERO;
                continue;
            }
            if sample > 0 {
                // Reduce positive sample modulo q
                let v = U256::from_u128(sample as u128);
                *c = v.rem(&context.q);
            } else {
                // Map negative sample to q - (|v| mod q)
                let v_abs = U256::from_u128((-sample) as u128).rem(&context.q);
                *c = if v_abs.is_zero().into() {
                    U256::ZERO
                } else {
                    q.wrapping_sub(&v_abs)
                };
            }
        }

        Self {
            coeffs,
            context: context.clone(),
        }
    }

    fn sample_tribits<R: Rng>(
        hamming_weight: usize,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self {
        let ternary = ternary_coefficients::<DEGREE, R>(hamming_weight, rng);
        let mut coeffs = [U256::ZERO; DEGREE];

        for (i, &val) in ternary.iter().enumerate() {
            coeffs[i] = if val == -1 {
                context.q.wrapping_sub(&U256::ONE) // modulus - 1
            } else {
                U256::from(val as u64)
            };
        }

        Self {
            coeffs,
            context: context.clone(),
        }
    }

    fn sample_noise<R: Rng>(
        variance: f64,
        context: &Self::Context,
        rng: &mut R,
    ) -> Self {
        Self::sample_gaussian(variance.sqrt(), context, rng)
    }
}

impl<const DEGREE: usize> PolyRescale<DEGREE> for BigIntPolyRing<DEGREE> {
    fn rescale_assign(&mut self, scale_factor: f64) {
        // Extract scale_bits from scale_factor (assuming scale_factor = 2^k)
        let k = (scale_factor.log2().round()) as u32;
        let modulus = self.context.q.get();

        for coeff in &mut self.coeffs {
            *coeff = rescale_coeff_pow2_u256(*coeff, modulus, k);
        }
    }
}

/// Round-to-nearest division by 2^k in Z_q, sign-aware, mapping back to [0, q).
fn rescale_coeff_pow2_u256(c: U256, q: U256, k: u32) -> U256 {
    let half = q >> 1;
    // centered lift to Z: x in (-(q/2), q/2]
    let (neg, x_abs) = if c > half {
        (true, q.wrapping_sub(&c))
    } else {
        (false, c)
    };

    // y = round(x_abs / 2^k) = (x_abs + 2^(k-1)) >> k
    let round = U256::ONE << (k - 1);
    let y = x_abs.saturating_add(&round) >> k;

    // map back to [0, q)
    if neg {
        if y.is_zero().into() {
            U256::ZERO
        } else {
            q.wrapping_sub(&y)
        }
    } else {
        y
    }
}

pub fn rescale_ciphertext_u256_inplace<const DEGREE: usize>(
    ct: &mut Ciphertext<BigIntPolyRing<DEGREE>, DEGREE>,
    k: u32,
) {
    let q = ct.c0.context.q.get();

    for i in 0..DEGREE {
        let c0 = ct.c0.coeffs[i];
        let c1 = ct.c1.coeffs[i];

        let r0 = rescale_coeff_pow2_u256(c0, q, k);
        let r1 = rescale_coeff_pow2_u256(c1, q, k);

        ct.c0.coeffs[i] = r0;
        ct.c1.coeffs[i] = r1;
    }

    ct.logp -= k; // Rescaling reduces precision
}
