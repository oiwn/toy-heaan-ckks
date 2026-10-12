//! Encrypted standard deviation over range-normalized sensor readings,
//! composed from the library's generic pipeline ops.

use rand::Rng;
use toy_heaan_ckks::crypto::CkksEngine;
use toy_heaan_ckks::crypto::pipeline::{
    Ct, LeveledEvaluator, align_level, eval_chebyshev_encrypted, fit_chebyshev,
    mul_plain_scalar,
};
use toy_heaan_ckks::rings::backends::rns_ntt::RnsPoly;

use crate::contract::SensorContract;

/// Full encrypted standard deviation over packed, normalized slots:
///
/// ```text
/// sum      = sum_slots(ct)                     -- rotation tree
/// mean     = sum * (1/slots)                   -- mul_plain + rescale
/// centered = ct - mean                         -- sub (level-aligned)
/// sq       = centered * centered               -- ct x ct + rescale
/// variance = sum_slots(sq) * (1/slots)         -- rotation tree + mul_plain
/// std_n    = Chebyshev-sqrt(variance)          -- encrypted Clenshaw
/// ```
///
/// `ct` must hold readings normalized with [`SensorContract::normalize`].
/// The result is the std in normalized units; decrypt and map it back with
/// [`SensorContract::restore_std`]. The sqrt is fitted on
/// [`SensorContract::variance_interval`].
///
/// Levels consumed: 3 + degree + 1 (mean, square, variance, affine map,
/// `degree - 1` recurrence muls, final mul). With 61-bit primes the
/// ciphertext must start with at least `degree + 6` RNS channels (the final
/// Clenshaw mul needs >= 3).
///
/// # Panics
/// Panics if the ciphertext runs out of channels.
pub fn encrypted_std_dev<R: Rng, const DEGREE: usize>(
    ct: &Ct<DEGREE>,
    slots: usize,
    contract: &SensorContract,
    degree: usize,
    evaluator: &mut LeveledEvaluator<DEGREE>,
    rng: &mut R,
) -> Ct<DEGREE> {
    let (var_lo, var_hi) = contract.variance_interval();
    let coeffs = fit_chebyshev(f64::sqrt, var_lo, var_hi, degree);

    // mean = sum * 1/slots
    let sum_ct = evaluator.sum_slots(ct, slots, rng);
    let mean_ct = mul_plain_scalar(&sum_ct, 1.0 / slots as f64, ct.logp, slots);

    // centered = input - mean, at the mean's level
    let ct_aligned = align_level(ct, mean_ct.c0.basis().clone(), mean_ct.logq);
    let centered = CkksEngine::<RnsPoly<DEGREE>, DEGREE>::sub_ciphertexts(
        &ct_aligned,
        &mean_ct,
    );

    // variance = sum of squared deviations * 1/slots
    let sq_diff = evaluator.mul(&centered, &centered, rng);
    let var_sum = evaluator.sum_slots(&sq_diff, slots, rng);
    let variance =
        mul_plain_scalar(&var_sum, 1.0 / slots as f64, sq_diff.logp, slots);

    // std_n = sqrt(variance) via encrypted Clenshaw on the contract envelope
    eval_chebyshev_encrypted(
        &variance, &coeffs, var_lo, var_hi, slots, evaluator, rng,
    )
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;
    use toy_heaan_ckks::crypto::engine::CkksParams;
    use toy_heaan_ckks::crypto::pipeline::{
        decode_wide, eval_chebyshev_clenshaw, to_chebyshev_domain,
    };
    use toy_heaan_ckks::encoding::ckks_encoder::CkksEncoder;
    use toy_heaan_ckks::math::generate_primes;
    use toy_heaan_ckks::rings::backends::rns_ntt::RnsBasis;

    use super::*;

    /// Full encrypted std-dev at toy params with the demo contract:
    /// N = 16, 8 diurnal-style readings, degree 7. Levels: 3 + 8 = 11, and
    /// the final mul needs >= 3 channels -> 13 primes, 2 left at the end.
    /// The fit interval is narrow ([1e-4, 6.4e-3]), so the affine map
    /// multiplies by ~635 — no longer an encoder problem.
    #[test]
    fn encrypted_std_dev_end_to_end() {
        const N: usize = 16;
        const SLOTS: usize = 8;
        const SCALE: u32 = 61;

        let contract = SensorContract::demo();
        let degree = 7;

        let readings: [f64; SLOTS] =
            [12.0, 14.0, 17.0, 20.0, 21.0, 19.0, 16.0, 13.0];
        let mean = readings.iter().sum::<f64>() / SLOTS as f64;
        let std_dev = (readings.iter().map(|x| (x - mean).powi(2)).sum::<f64>()
            / SLOTS as f64)
            .sqrt();
        assert!(
            contract.contains_std(std_dev),
            "test data std {std_dev} must honor the contract envelope"
        );
        let normalized: Vec<f64> =
            readings.iter().map(|&x| contract.normalize(x)).collect();

        let primes = generate_primes(SCALE as usize, degree + 6, N as u64);
        let basis = Arc::new(RnsBasis::<N>::new(primes).expect("valid NTT primes"));
        let logq = basis.total_bits();

        let params = CkksParams {
            error_variance: 3.2,
            hamming_weight: N / 2,
            scale_bits: SCALE,
        };
        let engine =
            CkksEngine::<RnsPoly<N>, N>::new(basis.clone(), params.clone());
        let mut rng = ChaCha20Rng::seed_from_u64(4321);

        let sk = engine.generate_secret_key(&mut rng).expect("sk");
        let pk = engine.generate_public_key(&sk, &mut rng).expect("pk");

        let encoder = CkksEncoder::<N>::new(SCALE);
        let pt = encoder.encode(&normalized, basis);
        let ct = engine.encrypt(&pt, &pk, logq, &mut rng);

        let mut evaluator = LeveledEvaluator::new(params, sk);
        let result = encrypted_std_dev(
            &ct,
            SLOTS,
            &contract,
            degree,
            &mut evaluator,
            &mut rng,
        );
        assert_eq!(
            result.c0.basis().channel_count(),
            2,
            "13 primes - 11 levels = 2 channels left"
        );

        let decoded = decode_wide(&evaluator.decrypt(&result));

        // Encryption fidelity: match the plaintext Chebyshev evaluation.
        let (var_lo, var_hi) = contract.variance_interval();
        let coeffs = fit_chebyshev(f64::sqrt, var_lo, var_hi, degree);
        let variance_normalized = (std_dev / contract.width()).powi(2);
        let cheb_plain = contract.restore_std(eval_chebyshev_clenshaw(
            to_chebyshev_domain(variance_normalized, var_lo, var_hi),
            &coeffs,
        ));

        for (slot, &got) in decoded.iter().enumerate() {
            let got = contract.restore_std(got);
            assert!(
                (got - cheb_plain).abs() < 1e-3,
                "slot {slot}: encrypted {got} vs plaintext-chebyshev {cheb_plain}"
            );
            // End-to-end: within the degree-7 approximation budget.
            assert!(
                (got - std_dev).abs() < 0.1,
                "slot {slot}: encrypted std {got} vs true std {std_dev}"
            );
        }
    }
}
