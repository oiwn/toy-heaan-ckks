//! Homomorphic operations on ciphertexts.
//!
//! Neg is a `std::ops` impl (`-&ct`); everything else is a `Ciphertext`
//! method. Additions and subtractions are fallible (`ct.add(&other)?`): the
//! operands must share level, ring context and a compatible scale
//! ([`Scale::is_compatible`]). Nothing rescales or mod-switches on its own:
//! multiplications multiply the scales, and the caller decides when to call
//! [`Ciphertext::rescale`] or [`Ciphertext::match_level`].

use std::ops::Neg;

use crate::PolyRing;
use crate::ckks::Scale;
use crate::ckks::encoding::ckks_encoder::CkksEncoder;
use crate::ckks::errors::{CkksError, CkksResult};
use crate::ckks::keys::{RelinearizationKey, RotationKey};
use crate::ckks::types::{Ciphertext, Plaintext};
use crate::rings::RingError;

// ── Additive ops (any ring backend) ──────────────────────────────────────────

impl<P: PolyRing<DEGREE>, const DEGREE: usize> Neg for &Ciphertext<P, DEGREE> {
    type Output = Ciphertext<P, DEGREE>;

    fn neg(self) -> Self::Output {
        Ciphertext {
            c0: -self.c0.clone(),
            c1: -self.c1.clone(),
            scale: self.scale,
        }
    }
}

impl<P: PolyRing<DEGREE>, const DEGREE: usize> Ciphertext<P, DEGREE> {
    /// Checks that a right-hand operand at `poly`'s level and context with
    /// encoding `scale` can be combined with `self`.
    fn check_compatible(&self, poly: &P, scale: Scale) -> CkksResult<()> {
        if self.level() != poly.level() {
            return Err(CkksError::LevelMismatch {
                expected: self.level(),
                actual: poly.level(),
            });
        }
        if !self.c0.same_context(&self.c1) || !self.c0.same_context(poly) {
            return Err(RingError::IncompatibleContext.into());
        }
        if !self.scale.is_compatible(scale) {
            return Err(CkksError::ScaleMismatch {
                expected: self.scale.log2(),
                actual: scale.log2(),
            });
        }
        Ok(())
    }

    /// Component-wise sum `(c0 + c0', c1 + c1')`.
    ///
    /// # Errors
    /// Returns [`CkksError::LevelMismatch`], a ring
    /// [`RingError::IncompatibleContext`] or [`CkksError::ScaleMismatch`] if
    /// the operands differ in level, context or scale. Levels are never
    /// aligned implicitly; use [`Ciphertext::match_level`] first.
    pub fn add(&self, other: &Self) -> CkksResult<Self> {
        self.check_compatible(&other.c0, other.scale)?;
        if !other.c0.same_context(&other.c1) {
            return Err(RingError::IncompatibleContext.into());
        }
        let mut c0 = self.c0.clone();
        c0 += &other.c0;
        let mut c1 = self.c1.clone();
        c1 += &other.c1;
        Ok(Self {
            c0,
            c1,
            scale: self.scale,
        })
    }

    /// Component-wise difference `self - other`.
    ///
    /// # Errors
    /// Same conditions as [`Ciphertext::add`].
    pub fn sub(&self, other: &Self) -> CkksResult<Self> {
        self.add(&-other)
    }

    /// Adds a known plaintext: only `c0` changes, since `(c0 + m', c1)`
    /// decrypts to `m + m'`.
    ///
    /// # Errors
    /// Same conditions as [`Ciphertext::add`]: `pt` must be encoded at the
    /// ciphertext's scale in its basis.
    pub fn add_plain(&self, pt: &Plaintext<P, DEGREE>) -> CkksResult<Self> {
        self.check_compatible(&pt.poly, pt.scale)?;
        Ok(self.add_plain_unchecked(pt))
    }

    /// Subtracts a known plaintext: `(c0 - m', c1)`.
    ///
    /// # Errors
    /// Same conditions as [`Ciphertext::add_plain`].
    pub fn sub_plain(&self, pt: &Plaintext<P, DEGREE>) -> CkksResult<Self> {
        self.check_compatible(&pt.poly, pt.scale)?;
        let mut c0 = -pt.poly.clone();
        c0 += &self.c0;
        Ok(Self {
            c0,
            c1: self.c1.clone(),
            scale: self.scale,
        })
    }

    /// [`Ciphertext::add_plain`] for plaintexts the caller encoded at this
    /// ciphertext's scale and context.
    fn add_plain_unchecked(&self, pt: &Plaintext<P, DEGREE>) -> Self {
        let mut c0 = self.c0.clone();
        c0 += &pt.poly;
        Self {
            c0,
            c1: self.c1.clone(),
            scale: self.scale,
        }
    }
}

// ── Multiplicative ops and rotation ──────────────────────────────────────────

impl<P: PolyRing<DEGREE>, const DEGREE: usize> Ciphertext<P, DEGREE> {
    /// Conjugate all slots using a key generated at this level.
    ///
    /// # Errors
    /// Returns an error if the key belongs to a different level or context.
    pub fn conjugate(
        &self,
        key: &crate::ckks::keys::ConjugationKey<P, DEGREE>,
    ) -> crate::ckks::errors::CkksResult<Self> {
        if !self.c0.same_context(&self.c1)
            || !key.ksk.a.first().is_some_and(|a| self.c0.same_context(a))
        {
            return Err(crate::rings::RingError::IncompatibleContext.into());
        }
        let mut c0 = self.c0.automorphism((2 * DEGREE - 1) as u64);
        let (extra, c1) = key
            .ksk
            .switch(&self.c1.automorphism((2 * DEGREE - 1) as u64));
        c0 += &extra;
        Ok(Self {
            c0,
            c1,
            scale: self.scale,
        })
    }
    /// Ciphertext x ciphertext multiplication with gadget relinearization.
    ///
    /// Computes the tensor product `(c0 + c1*s)(c0' + c1'*s) = d0 + d1*s +
    /// d2*s^2`, then key-switches the `d2*s^2` term back to `s` with the
    /// relinearization key, yielding a valid degree-1 ciphertext.
    ///
    /// The output scale is `self.scale * other.scale`; call
    /// [`Ciphertext::rescale`] afterwards to bring it back down.
    ///
    /// # Panics
    /// Panics if the operands sit at different levels, or if the product
    /// scale overflows `f64` (beyond `2^1023`).
    pub fn mul(&self, other: &Self, rlk: &RelinearizationKey<P, DEGREE>) -> Self {
        assert_eq!(
            self.level(),
            other.level(),
            "level mismatch in multiplication"
        );

        // ── Tensor product ────────────────────────────────────────────────────
        let mut d0 = self.c0.clone();
        d0 *= &other.c0; // c0 * c0'

        let mut d1a = self.c0.clone();
        d1a *= &other.c1; // c0 * c1'
        let mut d1b = self.c1.clone();
        d1b *= &other.c0; // c1 * c0'
        let mut d1 = d1a;
        d1 += &d1b;

        let mut d2 = self.c1.clone();
        d2 *= &other.c1; // c1 * c1'

        // ── Relinearization: replace d2*s^2 ──────────────────────────────────
        let (relin_c0, relin_c1) = rlk.ksk.switch(&d2);
        d0 += &relin_c0;
        d1 += &relin_c1;

        Self {
            c0: d0,
            c1: d1,
            scale: self
                .scale
                .checked_mul(other.scale)
                .expect("mul: product scale overflows f64"),
        }
    }

    /// Ciphertext x plaintext multiplication: `(c0, c1) -> (c0*p, c1*p)`.
    ///
    /// Since `p` is known, the product still decrypts under the same key — no
    /// `s^2` term appears, so no relinearization is needed. `pt` must be
    /// encoded in the ciphertext's basis; the scale is multiplied by
    /// `pt.scale`.
    ///
    /// # Panics
    /// Panics if the product scale overflows `f64` (beyond `2^1023`).
    pub fn mul_plain(&self, pt: &Plaintext<P, DEGREE>) -> Self {
        let mut c0 = self.c0.clone();
        c0 *= &pt.poly;
        let mut c1 = self.c1.clone();
        c1 *= &pt.poly;
        Self {
            c0,
            c1,
            scale: self
                .scale
                .checked_mul(pt.scale)
                .expect("mul_plain: product scale overflows f64"),
        }
    }

    /// Multiplies every slot by `value`, encoded at `scale`; the ciphertext
    /// scale is multiplied by `scale` and no prime is consumed.
    ///
    /// With `scale = self.scale` follow it by a rescale (plain scalar
    /// multiply). Small exact encodings need no rescale: `(0.5, 2^1)` is the
    /// integer 1 (halves the slots, doubles the scale) and `(1.0, 2^k)` is
    /// `2^k` (a level-free scale lift).
    ///
    /// # Panics
    /// Panics if `value * scale` does not fit the basis modulus, or if the
    /// product scale overflows `f64`.
    pub fn mul_const(&self, value: f64, scale: Scale) -> Self {
        let pt = CkksEncoder::<DEGREE>::new(scale)
            .encode_constant(value, self.c0.context().clone());
        self.mul_plain(&pt)
    }

    /// Multiplies every slot by `value` and rescales so the result lands at
    /// `target` scale (up to `f64` roundoff, see [`Scale::is_compatible`]).
    ///
    /// `value` is encoded at `target * q / self.scale`, where `q` is the
    /// prime this rescale divides by. Use it to bring a value onto another
    /// ciphertext's scale before adding the two.
    ///
    /// # Errors
    /// Returns an error at level 1, if a scale is unrepresentable, or if the
    /// encoded constant does not fit the modulus.
    pub fn mul_const_to_scale(
        &self,
        value: f64,
        target: Scale,
    ) -> CkksResult<Self> {
        let divisor = Scale::new(self.c0.rescale_divisor() as f64)?;
        let plaintext_scale =
            target.checked_mul(divisor)?.checked_div(self.scale)?;
        let pt = CkksEncoder::<DEGREE>::new(plaintext_scale)
            .try_encode_constant(value, self.c0.context().clone())?;
        let scale = self.scale.checked_mul(plaintext_scale)?;
        let mut c0 = self.c0.clone();
        c0 *= &pt.poly;
        let mut c1 = self.c1.clone();
        c1 *= &pt.poly;
        Self { c0, c1, scale }.rescale()
    }

    /// Adds `value` to every slot (encoded at `self.scale`); consumes no
    /// level.
    ///
    /// # Panics
    /// Panics if `value * scale` does not fit the basis modulus.
    pub fn add_const(&self, value: f64) -> Self {
        let pt = CkksEncoder::<DEGREE>::new(self.scale)
            .encode_constant(value, self.c0.context().clone());
        // Encoded at `self.scale` in `self`'s context: always compatible.
        self.add_plain_unchecked(&pt)
    }

    /// The trivial ciphertext `(0, 0)` at this level and scale: decrypts to
    /// zero under any key. Seed for sums of constants, e.g.
    /// `ct.zero_like().add_const(c)`.
    pub fn zero_like(&self) -> Self {
        let zero = P::zero(self.c0.context());
        Self {
            c0: zero.clone(),
            c1: zero,
            scale: self.scale,
        }
    }

    /// Rotates the slots by `rotk.rotation` positions.
    ///
    /// 1. Apply the automorphism `X -> X^(5^k)` to both components. This
    ///    rotates the slots but changes the decryption key from `s` to `s_k`.
    /// 2. Key-switch `sigma_k(c1)` back to `s` with the rotation key.
    ///
    /// The output keeps scale and level: rotation consumes no level.
    pub fn rotate(&self, rotk: &RotationKey<P, DEGREE>) -> Self {
        let mut c0 = self.c0.rotate_slots(rotk.rotation);
        let c1_rot = self.c1.rotate_slots(rotk.rotation);

        let (ks_c0, ks_c1) = rotk.ksk.switch(&c1_rot);
        c0 += &ks_c0;

        Self {
            c0,
            c1: ks_c1,
            scale: self.scale,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    use super::*;
    use crate::ckks::engine::{CkksEngine, CkksParams};
    use crate::ckks::keys::SecretKey;
    use crate::math::generate_primes;
    use crate::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

    const N: usize = 16;
    const SCALE: i32 = 30;

    type Ct = Ciphertext<RnsPoly<N>, N>;

    struct Setup {
        encoder: CkksEncoder<N>,
        sk: SecretKey<RnsPoly<N>, N>,
        rlk: RelinearizationKey<RnsPoly<N>, N>,
        a: Ct,
        b: Ct,
    }

    const A: [f64; 4] = [0.5, -0.25, 0.75, 0.1];
    const B: [f64; 4] = [0.2, 0.4, -0.6, 0.3];

    fn pow2(exponent: i32) -> Scale {
        Scale::power_of_two(exponent).expect("valid scale")
    }

    fn setup() -> Setup {
        setup_at(pow2(SCALE))
    }

    /// Three 30-bit primes; `a` and `b` encrypted at `scale`.
    fn setup_at(scale: Scale) -> Setup {
        let primes = generate_primes(SCALE as usize, 3, N as u64);
        let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid primes"));
        let params = CkksParams {
            error_variance: 3.2,
            hamming_weight: N / 2,
            scale,
        };
        let engine = CkksEngine::<RnsPoly<N>, N>::new(basis.clone(), params);
        let mut rng = ChaCha20Rng::seed_from_u64(9);
        let sk = engine.generate_secret_key(&mut rng).expect("sk");
        let pk = engine.generate_public_key(&sk, &mut rng).expect("pk");
        let rlk = engine.generate_relin_key(&sk, &mut rng);
        let encoder = CkksEncoder::<N>::new(scale);
        let a = engine.encrypt(&encoder.encode(&A, basis.clone()), &pk, &mut rng);
        let b = engine.encrypt(&encoder.encode(&B, basis), &pk, &mut rng);
        Setup {
            encoder,
            sk,
            rlk,
            a,
            b,
        }
    }

    impl Setup {
        fn decrypt(&self, ct: &Ct) -> Vec<f64> {
            let sk = self.sk.mod_switch_to(ct.level()).expect("sk level");
            let pt = CkksEngine::<RnsPoly<N>, N>::decrypt(ct, &sk);
            self.encoder.decode(&pt).into_iter().take(A.len()).collect()
        }
    }

    fn assert_close(got: &[f64], want: impl Iterator<Item = f64>, tol: f64) {
        for (slot, (g, w)) in got.iter().zip(want).enumerate() {
            assert!((g - w).abs() < tol, "slot {slot}: got {g}, want {w}");
        }
    }

    #[test]
    fn add_sub_neg_match_plaintext() {
        let s = setup();
        let pairs = || A.iter().zip(B.iter());
        assert_close(
            &s.decrypt(&s.a.add(&s.b).expect("add")),
            pairs().map(|(x, y)| x + y),
            1e-6,
        );
        assert_close(
            &s.decrypt(&s.a.sub(&s.b).expect("sub")),
            pairs().map(|(x, y)| x - y),
            1e-6,
        );
        assert_close(&s.decrypt(&-&s.a), A.iter().map(|x| -x), 1e-6);
    }

    #[test]
    fn plaintext_add_and_sub() {
        let s = setup();
        let pt = s.encoder.encode(&B, s.a.c0.context().clone());
        let pairs = || A.iter().zip(B.iter());
        let sum = s.a.add_plain(&pt).expect("add_plain");
        assert_close(&s.decrypt(&sum), pairs().map(|(x, y)| x + y), 1e-6);
        let diff = s.a.sub_plain(&pt).expect("sub_plain");
        assert_close(&s.decrypt(&diff), pairs().map(|(x, y)| x - y), 1e-6);
    }

    /// Runs all four additive ops of `s.a` against a ciphertext and a
    /// plaintext that are both incompatible in the same way.
    fn additive_errors(
        s: &Setup,
        ct: &Ct,
        pt: &Plaintext<RnsPoly<N>, N>,
    ) -> [CkksError; 4] {
        [
            s.a.add(ct).expect_err("add"),
            s.a.sub(ct).expect_err("sub"),
            s.a.add_plain(pt).expect_err("add_plain"),
            s.a.sub_plain(pt).expect_err("sub_plain"),
        ]
    }

    #[test]
    fn scale_mismatch_is_an_error() {
        let s = setup();
        let ct = s.b.mul_const(1.0, pow2(1));
        let pt = CkksEncoder::<N>::new(pow2(SCALE + 1))
            .encode(&B, s.a.c0.context().clone());
        for err in additive_errors(&s, &ct, &pt) {
            assert!(matches!(err, CkksError::ScaleMismatch { .. }), "{err}");
        }
    }

    #[test]
    fn level_mismatch_is_an_error() {
        let s = setup();
        let ct = s.b.mod_switch_to(2).expect("mod switch");
        let pt = s.encoder.encode(&B, ct.c0.context().clone());
        for err in additive_errors(&s, &ct, &pt) {
            assert!(
                matches!(
                    err,
                    CkksError::LevelMismatch {
                        expected: 3,
                        actual: 2
                    }
                ),
                "{err}"
            );
        }
    }

    #[test]
    fn context_mismatch_is_an_error() {
        let s = setup();
        // Same primes and level, but a different basis instance.
        let primes = s.a.c0.basis().moduli().to_vec();
        let other = Arc::new(RnsBasis::<N>::new(primes).expect("valid primes"));
        let zero = RnsPoly::<N>::zero(&other);
        let ct = Ct {
            c0: zero.clone(),
            c1: zero,
            scale: pow2(SCALE),
        };
        let pt = s.encoder.encode(&B, other);
        for err in additive_errors(&s, &ct, &pt) {
            assert!(
                matches!(err, CkksError::Ring(RingError::IncompatibleContext)),
                "{err}"
            );
        }
    }

    #[test]
    fn scalar_ops() {
        let s = setup();

        let scaled = s.a.mul_const(3.0, pow2(SCALE)).rescale().expect("rescale");
        assert_eq!(scaled.level(), 2);
        assert_close(&s.decrypt(&scaled), A.iter().map(|x| 3.0 * x), 1e-5);

        let shifted = s.a.add_const(-0.5);
        assert_eq!(shifted.level(), s.a.level());
        assert_close(&s.decrypt(&shifted), A.iter().map(|x| x - 0.5), 1e-6);

        let constant = s.a.zero_like().add_const(0.25);
        assert_close(&s.decrypt(&constant), A.iter().map(|_| 0.25), 1e-6);

        // (0.5 at 1 bit) is the integer 1: halves the slots, no rescale.
        let halved = s.a.mul_const(0.5, pow2(1));
        assert_eq!(
            (halved.level(), halved.scale),
            (s.a.level(), pow2(SCALE + 1))
        );
        assert_close(&s.decrypt(&halved), A.iter().map(|x| x / 2.0), 1e-6);
    }

    /// Scales need not be powers of two: add and multiply at 1.5 * 2^30.
    #[test]
    fn non_power_of_two_scale() {
        let scale = Scale::new(1.5 * 2f64.powi(SCALE)).expect("valid scale");
        let s = setup_at(scale);
        let pairs = || A.iter().zip(B.iter());

        let sum = s.a.add(&s.b).expect("add");
        assert_eq!(sum.scale, scale);
        assert_close(&s.decrypt(&sum), pairs().map(|(x, y)| x + y), 1e-6);

        // Decrypt the product before rescaling: the decoder divides by the
        // recorded `scale^2`, which is exact here.
        let product = s.a.mul(&s.b, &s.rlk);
        let squared = scale.checked_mul(scale).expect("product scale");
        assert_eq!(product.scale, squared);
        assert_close(&s.decrypt(&product), pairs().map(|(x, y)| x * y), 1e-6);

        let pt = s.encoder.encode(&B, s.a.c0.context().clone());
        let plain_product = s.a.mul_plain(&pt);
        assert_eq!(plain_product.scale, squared);
        assert_close(
            &s.decrypt(&plain_product),
            pairs().map(|(x, y)| x * y),
            1e-6,
        );
    }

    /// `mul_const_to_scale` lands on another ciphertext's scale, so the two
    /// can be added after a rescale by a prime that is not a power of two.
    #[test]
    fn mul_const_to_scale_meets_target() {
        let s = setup();
        // Start off-scale: 2^31 instead of the partner's 2^30.
        let lifted = s.a.mul_const(1.0, pow2(1));
        let tripled = lifted.mul_const_to_scale(3.0, s.b.scale).expect("to scale");
        assert_eq!(tripled.level(), 2);
        assert!(tripled.scale.is_compatible(s.b.scale));

        let partner = s.b.mod_switch_to(2).expect("mod switch");
        let sum = tripled.add(&partner).expect("compatible scales");
        let want = A.iter().zip(B.iter()).map(|(x, y)| 3.0 * x + y);
        assert_close(&s.decrypt(&sum), want, 1e-5);

        let bottom = sum.mod_switch_to(1).expect("mod switch");
        assert!(bottom.mul_const_to_scale(1.0, s.b.scale).is_err());
    }
}
