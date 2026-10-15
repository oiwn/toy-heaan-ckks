//! Homomorphic operations on ciphertexts.
//!
//! Add, Sub and Neg are `std::ops` impls (`&a + &b`, `&ct - &pt`, `-&ct`);
//! everything else is a `Ciphertext` method. Nothing rescales on its own:
//! multiplications raise `logp`, and the caller decides when to call
//! [`Ciphertext::rescale`].

use std::ops::{Add, Neg, Sub};

use crate::PolyRing;
use crate::ckks::encoding::ckks_encoder::CkksEncoder;
use crate::ckks::engine::{CkksEngine, RnsGadgetRelinKey, RnsGadgetRotationKey};
use crate::ckks::types::{Ciphertext, Plaintext};
use crate::rings::backends::rns_ntt::RnsPoly;
use crate::rings::traits::PolyAutomorphism;

// ── Additive ops (any ring backend) ──────────────────────────────────────────

impl<P: PolyRing<DEGREE>, const DEGREE: usize> Add for &Ciphertext<P, DEGREE> {
    type Output = Ciphertext<P, DEGREE>;

    /// Component-wise sum; both operands must share `logp` and level.
    fn add(self, rhs: Self) -> Self::Output {
        assert_eq!(self.logp, rhs.logp, "logp mismatch in addition");
        let mut c0 = self.c0.clone();
        c0 += &rhs.c0;
        let mut c1 = self.c1.clone();
        c1 += &rhs.c1;
        Ciphertext {
            c0,
            c1,
            logp: self.logp,
        }
    }
}

impl<P: PolyRing<DEGREE>, const DEGREE: usize> Neg for &Ciphertext<P, DEGREE> {
    type Output = Ciphertext<P, DEGREE>;

    fn neg(self) -> Self::Output {
        Ciphertext {
            c0: -self.c0.clone(),
            c1: -self.c1.clone(),
            logp: self.logp,
        }
    }
}

impl<P: PolyRing<DEGREE>, const DEGREE: usize> Sub for &Ciphertext<P, DEGREE> {
    type Output = Ciphertext<P, DEGREE>;

    fn sub(self, rhs: Self) -> Self::Output {
        assert_eq!(self.logp, rhs.logp, "logp mismatch in subtraction");
        self + &(-rhs)
    }
}

impl<P: PolyRing<DEGREE>, const DEGREE: usize> Add<&Plaintext<P, DEGREE>>
    for &Ciphertext<P, DEGREE>
{
    type Output = Ciphertext<P, DEGREE>;

    /// Adds a known plaintext: only `c0` changes, since `(c0 + m', c1)`
    /// decrypts to `m + m'`. The plaintext must be encoded at `logp` in the
    /// ciphertext's basis.
    fn add(self, rhs: &Plaintext<P, DEGREE>) -> Self::Output {
        assert_eq!(self.logp, rhs.scale_bits, "logp mismatch in plaintext add");
        let mut c0 = self.c0.clone();
        c0 += &rhs.poly;
        Ciphertext {
            c0,
            c1: self.c1.clone(),
            logp: self.logp,
        }
    }
}

impl<P: PolyRing<DEGREE>, const DEGREE: usize> Sub<&Plaintext<P, DEGREE>>
    for &Ciphertext<P, DEGREE>
{
    type Output = Ciphertext<P, DEGREE>;

    fn sub(self, rhs: &Plaintext<P, DEGREE>) -> Self::Output {
        assert_eq!(self.logp, rhs.scale_bits, "logp mismatch in plaintext sub");
        let mut c0 = -rhs.poly.clone();
        c0 += &self.c0;
        Ciphertext {
            c0,
            c1: self.c1.clone(),
            logp: self.logp,
        }
    }
}

// ── Multiplicative ops and rotation (RNS-NTT backend) ────────────────────────

impl<const DEGREE: usize> Ciphertext<RnsPoly<DEGREE>, DEGREE> {
    /// Ciphertext x ciphertext multiplication with gadget relinearization.
    ///
    /// Computes the tensor product `(c0 + c1*s)(c0' + c1'*s) = d0 + d1*s +
    /// d2*s^2`, then replaces the `d2*s^2` term with
    /// `sum_i sum_j digit_j(alpha_i(d2)) * rlk[i][j]`, yielding a valid
    /// degree-1 ciphertext.
    ///
    /// The output has `logp = self.logp + other.logp`; call
    /// [`Ciphertext::rescale`] afterwards to bring it back down.
    ///
    /// # Panics
    /// Panics if the operands sit at different levels.
    pub fn mul(&self, other: &Self, rlk: &RnsGadgetRelinKey<DEGREE>) -> Self {
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
        d2.to_coeff_domain(); // the key switch reads channels[i]

        // ── Gadget relinearization: replace d2*s^2 ───────────────────────────
        let (relin_c0, relin_c1) =
            CkksEngine::<RnsPoly<DEGREE>, DEGREE>::gadget_key_switch(
                &d2, &rlk.a, &rlk.b,
            );
        d0 += &relin_c0;
        d1 += &relin_c1;

        Self {
            c0: d0,
            c1: d1,
            logp: self.logp + other.logp,
        }
    }

    /// Ciphertext x plaintext multiplication: `(c0, c1) -> (c0*p, c1*p)`.
    ///
    /// Since `p` is known, the product still decrypts under the same key — no
    /// `s^2` term appears, so no relinearization is needed. `pt` must be
    /// encoded in the ciphertext's basis; `logp` grows by `pt.scale_bits`.
    pub fn mul_plain(&self, pt: &Plaintext<RnsPoly<DEGREE>, DEGREE>) -> Self {
        let mut c0 = self.c0.clone();
        c0 *= &pt.poly;
        let mut c1 = self.c1.clone();
        c1 *= &pt.poly;
        Self {
            c0,
            c1,
            logp: self.logp + pt.scale_bits,
        }
    }

    /// Multiplies every slot by `value`, encoded at `scale_bits` precision;
    /// `logp` grows by `scale_bits` and no prime is consumed.
    ///
    /// With `scale_bits = logp` follow it by a rescale (plain scalar
    /// multiply). Small exact encodings need no rescale: `(0.5, 1)` is the
    /// integer 1 (halves the slots, `logp + 1`) and `(1.0, k)` is `2^k` (a
    /// level-free `logp` lift).
    ///
    /// # Panics
    /// Panics if `scale_bits` is 0 or `value * 2^scale_bits` does not fit the
    /// basis modulus.
    pub fn mul_const(&self, value: f64, scale_bits: u32) -> Self {
        let pt = CkksEncoder::<DEGREE>::new(scale_bits)
            .encode_constant(value, self.c0.basis().clone());
        self.mul_plain(&pt)
    }

    /// Adds `value` to every slot (encoded at `logp`); consumes no level.
    ///
    /// # Panics
    /// Panics if `value * 2^logp` does not fit the basis modulus.
    pub fn add_const(&self, value: f64) -> Self {
        let pt = CkksEncoder::<DEGREE>::new(self.logp)
            .encode_constant(value, self.c0.basis().clone());
        self + &pt
    }

    /// The trivial ciphertext `(0, 0)` at this level and `logp`: decrypts to
    /// zero under any key. Seed for sums of constants, e.g.
    /// `ct.zero_like().add_const(c)`.
    pub fn zero_like(&self) -> Self {
        let zero = RnsPoly::zero(self.c0.basis().clone());
        Self {
            c0: zero.clone(),
            c1: zero,
            logp: self.logp,
        }
    }

    /// Rotates the slots by `rotk.rotation` positions.
    ///
    /// 1. Apply the automorphism `X -> X^(5^k)` to both components. This
    ///    rotates the slots but changes the decryption key from `s` to `s_k`.
    /// 2. Key-switch `sigma_k(c1)` back to `s` with the gadget rotation key.
    ///
    /// The output keeps `logp` and level: rotation consumes no level.
    pub fn rotate(&self, rotk: &RnsGadgetRotationKey<DEGREE>) -> Self {
        let mut c0 = self.c0.rotate_slots(rotk.rotation);
        let mut c1_rot = self.c1.rotate_slots(rotk.rotation);
        c1_rot.to_coeff_domain(); // already coeff; explicit for clarity

        let (ks_c0, ks_c1) =
            CkksEngine::<RnsPoly<DEGREE>, DEGREE>::gadget_key_switch(
                &c1_rot, &rotk.a, &rotk.b,
            );
        c0 += &ks_c0;

        Self {
            c0,
            c1: ks_c1,
            logp: self.logp,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    use super::*;
    use crate::ckks::engine::CkksParams;
    use crate::ckks::keys::SecretKey;
    use crate::math::generate_primes;
    use crate::rings::backends::rns_ntt::RnsBasis;

    const N: usize = 16;
    const SCALE: u32 = 30;

    type Ct = Ciphertext<RnsPoly<N>, N>;

    struct Setup {
        encoder: CkksEncoder<N>,
        sk: SecretKey<RnsPoly<N>, N>,
        a: Ct,
        b: Ct,
    }

    const A: [f64; 4] = [0.5, -0.25, 0.75, 0.1];
    const B: [f64; 4] = [0.2, 0.4, -0.6, 0.3];

    fn setup() -> Setup {
        let primes = generate_primes(SCALE as usize, 3, N as u64);
        let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid primes"));
        let params = CkksParams {
            error_variance: 3.2,
            hamming_weight: N / 2,
            scale_bits: SCALE,
        };
        let engine = CkksEngine::<RnsPoly<N>, N>::new(basis.clone(), params);
        let mut rng = ChaCha20Rng::seed_from_u64(9);
        let sk = engine.generate_secret_key(&mut rng).expect("sk");
        let pk = engine.generate_public_key(&sk, &mut rng).expect("pk");
        let encoder = CkksEncoder::<N>::new(SCALE);
        let a = engine.encrypt(&encoder.encode(&A, basis.clone()), &pk, &mut rng);
        let b = engine.encrypt(&encoder.encode(&B, basis), &pk, &mut rng);
        Setup { encoder, sk, a, b }
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
            &s.decrypt(&(&s.a + &s.b)),
            pairs().map(|(x, y)| x + y),
            1e-6,
        );
        assert_close(
            &s.decrypt(&(&s.a - &s.b)),
            pairs().map(|(x, y)| x - y),
            1e-6,
        );
        assert_close(&s.decrypt(&-&s.a), A.iter().map(|x| -x), 1e-6);
    }

    #[test]
    fn plaintext_add_and_sub() {
        let s = setup();
        let pt = s.encoder.encode(&B, s.a.c0.basis().clone());
        let pairs = || A.iter().zip(B.iter());
        assert_close(&s.decrypt(&(&s.a + &pt)), pairs().map(|(x, y)| x + y), 1e-6);
        assert_close(&s.decrypt(&(&s.a - &pt)), pairs().map(|(x, y)| x - y), 1e-6);
    }

    #[test]
    fn scalar_ops() {
        let s = setup();

        let scaled = s.a.mul_const(3.0, SCALE).rescale().expect("rescale");
        assert_eq!(scaled.level(), 2);
        assert_close(&s.decrypt(&scaled), A.iter().map(|x| 3.0 * x), 1e-5);

        let shifted = s.a.add_const(-0.5);
        assert_eq!(shifted.level(), s.a.level());
        assert_close(&s.decrypt(&shifted), A.iter().map(|x| x - 0.5), 1e-6);

        let constant = s.a.zero_like().add_const(0.25);
        assert_close(&s.decrypt(&constant), A.iter().map(|_| 0.25), 1e-6);

        // (0.5 at 1 bit) is the integer 1: halves the slots, no rescale.
        let halved = s.a.mul_const(0.5, 1);
        assert_eq!((halved.level(), halved.logp), (s.a.level(), SCALE + 1));
        assert_close(&s.decrypt(&halved), A.iter().map(|x| x / 2.0), 1e-6);
    }
}
