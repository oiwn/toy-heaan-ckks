//! Educational, fully packed complex CKKS bootstrap.
//!
//! Toy profiles only. Evaluation uses public keys; secret-key diagnostics
//! belong to callers. Hidden slot magnitudes must be at most one.

use std::collections::HashMap;

use num_complex::Complex64;
use rand::Rng;

use crate::PolyRing;
use crate::ckks::encoding::ckks_encoder::CkksEncoder;
use crate::ckks::errors::{CkksError, CkksResult};
use crate::ckks::keys::{KeySwitchingKey, RelinearizationKey, SecretKey};
use crate::ckks::types::{Ciphertext, Plaintext};

use super::chebyshev::fit_chebyshev;

fn invalid(message: &str) -> CkksError {
    CkksError::InvalidParameter {
        message: message.into(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ckks::engine::{CkksEngine, CkksParams};
    use crate::ckks::keys::{ConjugationKey, RotationKey};
    use crate::math::generate_primes;
    use crate::rings::backends::rns_ntt::{RnsBasis, RnsPoly};
    use rand::{RngExt, SeedableRng};
    use rand_chacha::ChaCha20Rng;
    use std::sync::Arc;

    fn assert_same<const N: usize>(
        left: &Ct<RnsPoly<N>, N>,
        right: &Ct<RnsPoly<N>, N>,
    ) {
        assert_eq!((left.level(), left.logp), (right.level(), right.logp));
        for (a, b) in [(&left.c0, &right.c0), (&left.c1, &right.c1)] {
            let mut difference = a.clone();
            difference -= b;
            assert_eq!(difference.to_coeffs(), [0; N]);
        }
    }

    fn suite<const N: usize>() {
        let primes = generate_primes(61, 18, N as u64);
        let basis = Arc::new(RnsBasis::<N>::new(primes.clone()).unwrap());
        let plan = BootstrapPlan::<RnsPoly<N>, N>::new(
            basis.clone(),
            BootstrapParameters::toy(primes[0]),
        )
        .unwrap();
        let encoder = CkksEncoder::new(44);
        let engine = CkksEngine::new(
            basis.clone(),
            CkksParams {
                error_variance: 3.2,
                hamming_weight: N / 2,
                scale_bits: 44,
            },
        );
        let other = Arc::new(RnsBasis::<N>::new(primes.clone()).unwrap());
        let bottom = basis.at_level(1).unwrap();
        let mut coeffs = [0i64; N];
        coeffs[0] = (primes[0] / 2) as i64;
        coeffs[1] = -coeffs[0];
        coeffs[2] = -1;
        let poly = RnsPoly::from_coeffs(&coeffs, &bottom);
        assert_eq!(poly.mod_raise(&basis).unwrap().to_coeffs(), coeffs);
        assert!(poly.mod_raise(&other).is_err());
        assert!(poly.mod_raise(&bottom).is_err());
        assert!(RnsPoly::zero(&basis).mod_raise(&basis).is_err());
        for seed in 0..8 {
            let mut rng = ChaCha20Rng::seed_from_u64(seed);
            let secret = engine.generate_secret_key(&mut rng).unwrap();
            let public = engine.generate_public_key(&secret, &mut rng).unwrap();
            let mut keys = plan
                .generate_keys(&secret, 3.2f64.sqrt(), &mut rng)
                .unwrap();
            if seed == 0 {
                let ctx = basis.at_level(16).unwrap();
                let enc = CkksEncoder::new(61);
                let pt: Plaintext<RnsPoly<N>, N> =
                    enc.encode(&vec![0.13; N / 2], ctx.clone());
                let t = Ciphertext {
                    c0: pt.poly,
                    c1: RnsPoly::zero(&ctx),
                    logp: 61,
                };
                let mut cache = HashMap::from([
                    (0, t.zero_like().add_const(1.0)),
                    (1, t.clone()),
                ]);
                let term = plan.term(31, &mut cache, &keys).unwrap();
                assert_eq!(term.level(), 11);
                let decoded = enc.decode(&CkksEngine::decrypt(
                    &term,
                    &secret.mod_switch_to(term.level()).unwrap(),
                ));
                for value in decoded {
                    assert!((value - (31.0 * 0.13f64.acos()).cos()).abs() < 1e-9);
                }
                let bound = (N / 2 + 2) / 2;
                for base in (0..2 * bound + 1).step_by(N / 2) {
                    let values: Vec<_> = (0..N / 2)
                        .map(|j| {
                            let integer =
                                (base + j).min(2 * bound) as f64 - bound as f64;
                            integer
                                + if j % 2 == 0 {
                                    2f64.powi(-17)
                                } else {
                                    -2f64.powi(-17)
                                }
                        })
                        .collect();
                    let pt: Plaintext<RnsPoly<N>, N> = enc.encode(
                        &values.iter().map(|x| x / plan.radius).collect::<Vec<_>>(),
                        ctx.clone(),
                    );
                    let ct = Ciphertext {
                        c0: pt.poly,
                        c1: RnsPoly::zero(&ctx),
                        logp: 61,
                    };
                    let reduced = plan.eval_mod(&ct, &keys).unwrap();
                    let decoded = enc.decode(&CkksEngine::decrypt(
                        &reduced,
                        &secret.mod_switch_to(reduced.level()).unwrap(),
                    ));
                    for (j, value) in decoded.into_iter().enumerate() {
                        let restored =
                            value / (std::f64::consts::TAU * 2f64.powi(-17));
                        assert!(
                            (restored - if j % 2 == 0 { 1.0 } else { -1.0 }).abs()
                                < 1e-4
                        );
                    }
                }
            }
            for case in 0..6 {
                let values: Vec<_> = (0..N / 2)
                    .map(|j| match case {
                        0 => Complex64::default(),
                        1 => Complex64::new(0.5, -0.25),
                        2 => Complex64::from_polar(0.8, j as f64 * 0.7),
                        3 => Complex64::from_polar(1.0, j as f64 * 0.9),
                        4 => Complex64::from_polar(
                            rng.random::<f64>().sqrt(),
                            rng.random::<f64>() * std::f64::consts::TAU,
                        ),
                        _ => {
                            Complex64::new(if j % 2 == 0 { 1.0 } else { -1.0 }, 0.0)
                        }
                    })
                    .collect();
                let initial = engine.encrypt(
                    &encoder.encode_complex(&values, basis.clone()),
                    &public,
                    &mut rng,
                );
                if case == 2 {
                    if seed == 0 {
                        let pt: Plaintext<RnsPoly<N>, N> =
                            encoder.encode_complex(&values, basis.clone());
                        let coefficients = pt.poly.to_centered_f64();
                        let trivial = Ciphertext {
                            c0: pt.poly,
                            c1: RnsPoly::zero(&basis),
                            logp: 44,
                        };
                        let transformed =
                            plan.transform(&trivial, &plan.forward, &keys).unwrap();
                        let decoded = encoder.decode_complex(&CkksEngine::decrypt(
                            &transformed,
                            &secret.mod_switch_to(transformed.level()).unwrap(),
                        ));
                        for j in 0..N / 2 {
                            let want = Complex64::new(
                                coefficients[j],
                                coefficients[j + N / 2],
                            ) / 2f64.powi(44)
                                * 2f64.powi(61)
                                / primes[0] as f64;
                            let want = want / (2.0 * plan.radius);
                            assert!((decoded[j] - want).norm() < 1e-9);
                        }
                    }
                    for offset in [-9, -1, 0, 1, 9] {
                        let key = RotationKey::generate(
                            &secret,
                            offset,
                            3.2f64.sqrt(),
                            &mut rng,
                        );
                        let actual = encoder.decode_complex(&CkksEngine::decrypt(
                            &initial.rotate(&key),
                            &secret,
                        ));
                        for (j, a) in actual.iter().enumerate() {
                            assert!(
                                (*a - values[(j as i32 + offset)
                                    .rem_euclid((N / 2) as i32)
                                    as usize])
                                    .norm()
                                    < 1e-3,
                                "N={N} offset={offset} slot={j} actual={a}"
                            );
                        }
                    }
                    let key =
                        ConjugationKey::generate(&secret, 3.2f64.sqrt(), &mut rng);
                    let actual = encoder.decode_complex(&CkksEngine::decrypt(
                        &initial.conjugate(&key).unwrap(),
                        &secret,
                    ));
                    for (a, v) in actual.iter().zip(&values) {
                        assert!((*a - v.conj()).norm() < 1e-3);
                    }
                }
                let relin =
                    RelinearizationKey::generate(&secret, 3.2f64.sqrt(), &mut rng);
                let mut squared = initial.mul(&initial, &relin).rescale().unwrap();
                let mut expected: Vec<_> = values.iter().map(|x| x * x).collect();
                for scale in [0, 45, u32::MAX] {
                    let mut wrong_scale = squared.clone();
                    wrong_scale.logp = scale;
                    assert!(plan.bootstrap(&wrong_scale, &keys).is_err());
                }
                if seed == 0 && (case == 0 || case == 2) {
                    let snapshot = initial.clone();
                    let fresh = plan.bootstrap(&initial, &keys).unwrap();
                    let decoded = encoder.decode_complex(&CkksEngine::decrypt(
                        &fresh,
                        &secret.mod_switch_to(fresh.level()).unwrap(),
                    ));
                    for (a, b) in decoded.iter().zip(&values) {
                        assert!((*a - *b).norm() < 1e-3);
                    }
                    assert_same(&snapshot, &initial);
                    for level in 1..=18 {
                        let at_level = initial.mod_switch_to(level).unwrap();
                        let prepared =
                            plan.prepare_input(&at_level, &keys).unwrap();
                        assert_same(&prepared, &initial.mod_switch_to(1).unwrap());
                    }
                    let mut wrong_context = squared.clone();
                    wrong_context.c0 =
                        RnsPoly::zero(&other.at_level(squared.level()).unwrap());
                    wrong_context.c1 = wrong_context.c0.clone();
                    assert!(plan.bootstrap(&wrong_context, &keys).is_err());
                    wrong_context.c1 = squared.c1.clone();
                    assert!(plan.bootstrap(&wrong_context, &keys).is_err());
                    let other_plan = BootstrapPlan::<RnsPoly<N>, N>::new(
                        other.clone(),
                        BootstrapParameters::toy(primes[0]),
                    )
                    .unwrap();
                    assert!(other_plan.bootstrap(&squared, &keys).is_err());
                }
                if case == 0 {
                    let key = keys.relin.remove(&9).unwrap();
                    assert!(plan.bootstrap(&squared, &keys).is_err());
                    keys.relin.insert(9, key);
                }
                for cycle in 0..2 {
                    let snapshot = squared.clone();
                    let explicit =
                        squared.mul_const(1.0, 17).mod_switch_to(1).unwrap();
                    let prepared = plan.prepare_input(&squared, &keys).unwrap();
                    assert_same(&prepared, &explicit);
                    assert_same(
                        &plan.prepare_input(&prepared, &keys).unwrap(),
                        &prepared,
                    );
                    let refreshed = plan.bootstrap(&squared, &keys).unwrap();
                    assert_same(&squared, &snapshot);
                    if seed == 0 && case == 0 && cycle == 0 {
                        assert_same(
                            &refreshed,
                            &plan.bootstrap(&explicit, &keys).unwrap(),
                        );
                    }
                    assert_eq!((refreshed.level(), refreshed.logp), (6, 44));
                    let sk = secret.mod_switch_to(refreshed.level()).unwrap();
                    let actual = encoder
                        .decode_complex(&CkksEngine::decrypt(&refreshed, &sk));
                    let before = encoder.decode_complex(&CkksEngine::decrypt(
                        &squared,
                        &secret.mod_switch_to(squared.level()).unwrap(),
                    ));
                    for ((a, want), prior) in
                        actual.iter().zip(&expected).zip(before)
                    {
                        assert!(
                            (*a - *want).norm() < 1e-3,
                            "N={N} seed={seed} case={case} cycle={cycle} got={a} expected={want}"
                        );
                        assert!((*a - prior).norm() < 1e-3);
                    }
                    let key =
                        RelinearizationKey::generate(&sk, 3.2f64.sqrt(), &mut rng);
                    squared = refreshed.mul(&refreshed, &key).rescale().unwrap();
                    for x in &mut expected {
                        *x *= *x;
                    }
                }
            }
        }
    }

    #[test]
    fn complex_bootstrap_16() {
        suite::<16>();
    }
    #[test]
    fn complex_bootstrap_64() {
        suite::<64>();
    }
}

/// Parameters for the supported educational profiles.
#[derive(Clone, Debug)]
pub struct BootstrapParameters {
    /// Exact bottom prime, used for normalization.
    pub bottom_modulus: u64,
    /// Input and output log2 scale.
    pub scale_bits: u32,
    /// Working log2 scale and prime bit length.
    pub working_bits: u32,
    /// Cosine polynomial degree.
    pub polynomial_degree: usize,
    /// Number of double-angle steps.
    pub double_angles: u32,
}

impl BootstrapParameters {
    /// Candidate profile using near-2^61 NTT primes.
    pub fn toy(bottom_modulus: u64) -> Self {
        Self {
            bottom_modulus,
            scale_bits: 44,
            working_bits: 61,
            polynomial_degree: 31,
            double_angles: 4,
        }
    }
}

type Ct<P, const N: usize> = Ciphertext<P, N>;
type Matrix = Vec<Vec<Complex64>>;

/// Public precomputation; contains no secret key.
pub struct BootstrapPlan<P: PolyRing<N>, const N: usize> {
    parameters: BootstrapParameters,
    context: P::Context,
    forward: Vec<Plaintext<P, N>>,
    inverse: Vec<Plaintext<P, N>>,
    coefficients: Vec<f64>,
    radius: f64,
    level: usize,
}

/// Evaluation keys bound to one plan context. Contains no secret key.
pub struct BootstrapKeys<P: PolyRing<N>, const N: usize> {
    context: P::Context,
    relin: HashMap<usize, RelinearizationKey<P, N>>,
    rotations: HashMap<(usize, usize), KeySwitchingKey<P, N>>,
    conjugation: HashMap<usize, KeySwitchingKey<P, N>>,
}

impl<P: PolyRing<N>, const N: usize> BootstrapPlan<P, N> {
    /// Construct and validate a toy profile and its transform matrices.
    ///
    /// # Errors
    /// Rejects unsupported degrees, scales, primes or insufficient depth.
    pub fn new(
        context: P::Context,
        parameters: BootstrapParameters,
    ) -> CkksResult<Self> {
        if ![16, 64].contains(&N)
            || parameters.scale_bits != 44
            || parameters.working_bits != 61
            || parameters.polynomial_degree != 31
            || parameters.double_angles != 4
        {
            return Err(invalid("unsupported educational bootstrap profile"));
        }
        let zero = P::zero(&context);
        let level = zero.level();
        if level != 18 || zero.rescale_bits() != 61 {
            return Err(invalid("bootstrap requires at least 18 near-2^61 primes"));
        }
        let bottom = zero.mod_switch_to(1)?;
        // Check the supplied modulus exactly through its residue, without
        // decoding a ciphertext or reconstructing large integers as f64.
        let modulus = parameters.bottom_modulus;
        let mut modulus_coefficients = [0i64; N];
        modulus_coefficients[0] = modulus as i64;
        if !(1u64 << 60..1u64 << 61).contains(&modulus)
            || modulus > i64::MAX as u64
            || P::from_coeffs(&modulus_coefficients, bottom.context()).to_coeffs()
                [0]
                != 0
        {
            return Err(invalid("incorrect bottom modulus"));
        }
        for current in 2..=level {
            if zero
                .mod_switch_to(current)?
                .rescale_divisor()
                .abs_diff(1u64 << 61)
                > (1u64 << 20)
            {
                return Err(invalid("working primes must be within 2^20 of 2^61"));
            }
        }
        if modulus.abs_diff(1u64 << 61) > (1u64 << 20) {
            return Err(invalid("bottom prime must be near 2^61"));
        }
        let slots = N / 2;
        let mut inverse = vec![vec![Complex64::default(); slots]; slots];
        let mut exponent = 1usize;
        for row in &mut inverse {
            for (k, entry) in row.iter_mut().enumerate() {
                *entry = Complex64::from_polar(
                    1.0,
                    std::f64::consts::PI * (exponent * k) as f64 / N as f64,
                );
            }
            exponent = exponent * 5 % (2 * N);
        }
        let forward: Matrix = (0..slots)
            .map(|j| {
                (0..slots)
                    .map(|k| inverse[k][j].conj() / slots as f64)
                    .collect()
            })
            .collect();
        let radius = ((N / 2 + 2) / 2 + 1) as f64;
        let coefficients = fit_chebyshev(
            |t| (2.0 * std::f64::consts::PI * radius * t / 16.0).cos(),
            -1.0,
            1.0,
            31,
        );
        let encode_diagonals = |matrix: &Matrix,
                                factor: f64,
                                bits: u32,
                                level: usize|
         -> CkksResult<Vec<Plaintext<P, N>>> {
            let ctx = zero.mod_switch_to(level)?.context().clone();
            (0..slots)
                .map(|offset| {
                    let diagonal: Vec<_> = (0..slots)
                        .map(|j| matrix[j][(j + offset) % slots] * factor)
                        .collect();
                    Ok(CkksEncoder::new(bits)
                        .try_encode_complex(&diagonal, ctx.clone())?)
                })
                .collect()
        };
        let forward = encode_diagonals(
            &forward,
            2f64.powi(61) / modulus as f64 / (2.0 * radius),
            61,
            level,
        )?;
        let inverse = encode_diagonals(
            &inverse,
            modulus as f64 / 2f64.powi(44) / (2.0 * std::f64::consts::PI),
            44,
            level - 11,
        )?;
        Ok(Self {
            parameters,
            context,
            forward,
            inverse,
            coefficients,
            radius,
            level,
        })
    }

    /// Generate level-specific evaluation keys during client setup.
    ///
    /// # Errors
    /// Rejects wrong contexts and invalid noise parameters.
    pub fn generate_keys<R: Rng>(
        &self,
        secret: &SecretKey<P, N>,
        error_std: f64,
        rng: &mut R,
    ) -> CkksResult<BootstrapKeys<P, N>> {
        if !secret.poly.same_context(&P::zero(&self.context))
            || !error_std.is_finite()
            || error_std <= 0.0
        {
            return Err(invalid("invalid bootstrap secret context or key noise"));
        }
        let mut keys = BootstrapKeys {
            context: self.context.clone(),
            relin: HashMap::new(),
            rotations: HashMap::new(),
            conjugation: HashMap::new(),
        };
        let coefficients = secret.poly.to_coeffs();
        if coefficients.iter().any(|c| ![-1, 0, 1].contains(c))
            || coefficients.iter().filter(|&&c| c != 0).count() > N / 2
        {
            return Err(invalid(
                "bootstrap requires a ternary secret of weight at most N/2",
            ));
        }
        for level in 3..=self.level {
            let sk = secret.mod_switch_to(level)?;
            keys.relin
                .insert(level, RelinearizationKey::generate(&sk, error_std, rng));
            if level == self.level || level == self.level - 11 {
                for offset in 1..N / 2 {
                    let rotated = sk.poly.rotate_slots(offset as i32);
                    keys.rotations.insert(
                        (level, offset),
                        KeySwitchingKey::generate(&sk, &rotated, error_std, rng),
                    );
                }
            }
            if level == self.level - 1 {
                let conjugated = sk.poly.automorphism((2 * N - 1) as u64);
                keys.conjugation.insert(
                    level,
                    KeySwitchingKey::generate(&sk, &conjugated, error_std, rng),
                );
            }
        }
        Ok(keys)
    }

    fn constant(
        &self,
        ct: &Ct<P, N>,
        value: Complex64,
        bits: u32,
    ) -> CkksResult<Plaintext<P, N>> {
        Ok(CkksEncoder::new(bits)
            .try_encode_complex(&vec![value; N / 2], ct.c0.context().clone())?)
    }

    fn switch(
        &self,
        ct: &Ct<P, N>,
        exponent: u64,
        key: &KeySwitchingKey<P, N>,
    ) -> Ct<P, N> {
        let mut c0 = ct.c0.automorphism(exponent);
        let (extra, c1) = key.switch(&ct.c1.automorphism(exponent));
        c0 += &extra;
        Ciphertext {
            c0,
            c1,
            logp: ct.logp,
        }
    }

    fn transform(
        &self,
        ct: &Ct<P, N>,
        diagonals: &[Plaintext<P, N>],
        keys: &BootstrapKeys<P, N>,
    ) -> CkksResult<Ct<P, N>> {
        let mut sum = ct.zero_like();
        sum.logp += diagonals[0].scale_bits;
        for (offset, diagonal) in diagonals.iter().enumerate() {
            let rotated = if offset == 0 {
                ct.clone()
            } else {
                let key = keys
                    .rotations
                    .get(&(ct.level(), offset))
                    .ok_or_else(|| invalid("missing rotation key"))?;
                let mut exponent = 1u64;
                for _ in 0..offset {
                    exponent = exponent * 5 % (2 * N) as u64;
                }
                self.switch(ct, exponent, key)
            };
            sum = &sum + &rotated.mul_plain(diagonal);
        }
        sum.rescale()
    }

    fn multiply(
        &self,
        a: &Ct<P, N>,
        b: &Ct<P, N>,
        keys: &BootstrapKeys<P, N>,
    ) -> CkksResult<Ct<P, N>> {
        let level = a.level().min(b.level());
        if level < 3 {
            return Err(invalid("insufficient multiplication headroom"));
        }
        let key = keys
            .relin
            .get(&level)
            .ok_or_else(|| invalid("missing relinearization key"))?;
        a.mod_switch_to(level)?
            .mul(&b.mod_switch_to(level)?, key)
            .rescale()
    }

    fn term(
        &self,
        index: usize,
        cache: &mut HashMap<usize, Ct<P, N>>,
        keys: &BootstrapKeys<P, N>,
    ) -> CkksResult<Ct<P, N>> {
        if let Some(term) = cache.get(&index) {
            return Ok(term.clone());
        }
        let left = self.term(index / 2, cache, keys)?;
        let right = self.term(index.div_ceil(2), cache, keys)?;
        let product = self.multiply(&left, &right, keys)?;
        let doubled = &product + &product;
        let result = if index.is_multiple_of(2) {
            doubled.add_const(-1.0)
        } else {
            &doubled - &cache[&1].mod_switch_to(doubled.level())?
        };
        cache.insert(index, result.clone());
        Ok(result)
    }

    fn eval_mod(
        &self,
        ct: &Ct<P, N>,
        keys: &BootstrapKeys<P, N>,
    ) -> CkksResult<Ct<P, N>> {
        // Forward diagonals already apply the interval normalization.
        let t = ct.add_const(-0.25 / self.radius);
        let mut cache =
            HashMap::from([(0, t.zero_like().add_const(1.0)), (1, t.clone())]);
        for index in 2..=31 {
            self.term(index, &mut cache, keys)?;
        }
        let level = cache
            .values()
            .map(Ciphertext::level)
            .min()
            .unwrap_or(t.level());
        let base = t.mod_switch_to(level)?;
        let mut sum = base.zero_like();
        sum.logp += 61;
        for (index, coefficient) in self.coefficients.iter().enumerate() {
            let term = cache[&index].mod_switch_to(level)?;
            sum = &sum
                + &term.mul_plain(&self.constant(
                    &term,
                    Complex64::new(*coefficient, 0.0),
                    61,
                )?);
        }
        let mut result = sum.rescale()?;
        for _ in 0..4 {
            let square = self.multiply(&result, &result, keys)?;
            result = (&square + &square).add_const(-1.0);
        }
        Ok(result)
    }

    /// Refresh a fully packed complex ciphertext, preparing its scale and level.
    ///
    /// Accepts any level in this plan's chain and a positive scale at most
    /// the configured input scale. Preparation preserves represented values
    /// but cannot recover lost precision. The input is not modified.
    ///
    /// # Errors
    /// Rejects incompatible contexts, unsupported scales and missing keys.
    /// Slot magnitudes <=1 are a hidden-data precondition.
    pub fn bootstrap(
        &self,
        input: &Ct<P, N>,
        keys: &BootstrapKeys<P, N>,
    ) -> CkksResult<Ct<P, N>> {
        let prepared = self.prepare_input(input, keys)?;
        self.bootstrap_prepared(&prepared, keys)
    }

    fn prepare_input(
        &self,
        input: &Ct<P, N>,
        keys: &BootstrapKeys<P, N>,
    ) -> CkksResult<Ct<P, N>> {
        if input.logp == 0
            || input.logp > self.parameters.scale_bits
            || !input.c0.same_context(&input.c1)
            || !P::zero(&keys.context).same_context(&P::zero(&self.context))
        {
            return Err(invalid("invalid bootstrap ciphertext or keys"));
        }
        let expected = P::zero(&self.context).mod_switch_to(input.level())?;
        if !input.c0.same_context(&expected) {
            return Err(invalid("ciphertext belongs to a different context chain"));
        }
        let increase = self.parameters.scale_bits - input.logp;
        let normalized = if increase == 0 {
            input.clone()
        } else {
            // Scale both the integer representation and its interpretation
            // by the same exact power of two, without changing slot values.
            input.mul_const(1.0, increase)
        };
        normalized.mod_switch_to(1)
    }

    fn bootstrap_prepared(
        &self,
        input: &Ct<P, N>,
        keys: &BootstrapKeys<P, N>,
    ) -> CkksResult<Ct<P, N>> {
        let raised = Ciphertext {
            c0: input.c0.mod_raise(&self.context)?,
            c1: input.c1.mod_raise(&self.context)?,
            logp: 61,
        };
        let packed = self.transform(&raised, &self.forward, keys)?;
        let conjugation = keys
            .conjugation
            .get(&packed.level())
            .ok_or_else(|| invalid("missing conjugation key"))?;
        let conjugate = self.switch(&packed, (2 * N - 1) as u64, conjugation);
        // Half factors were fused into the forward diagonals.
        let real = &packed + &conjugate;
        let imaginary = Self::multiply_by_i(&(&packed - &conjugate), true);
        let real = self.eval_mod(&real, keys)?;
        let imaginary = self.eval_mod(&imaginary, keys)?;
        let combined = &real + &Self::multiply_by_i(&imaginary, false);
        self.transform(&combined, &self.inverse, keys)
    }

    fn multiply_by_i(input: &Ct<P, N>, negative: bool) -> Ct<P, N> {
        // X^(N/2) evaluates to i at every root indexed by 5^j.
        let mut coefficients = vec![0i64; N];
        coefficients[N / 2] = if negative { -1 } else { 1 };
        let unit = P::from_coeffs(&coefficients, input.c0.context());
        let mut c0 = input.c0.clone();
        c0 *= &unit;
        let mut c1 = input.c1.clone();
        c1 *= &unit;
        Ciphertext {
            c0,
            c1,
            logp: input.logp,
        }
    }
}
