//! Encrypted standard deviation over range-normalized sensor readings,
//! composed from the library's generic pipeline ops.

use std::time::{Duration, Instant};

use rand::Rng;
use toy_heaan_ckks::ckks::algorithms::chebyshev::{
    eval_chebyshev_encrypted, fit_chebyshev,
};
use toy_heaan_ckks::ckks::evaluator::LeveledEvaluator;
use toy_heaan_ckks::ckks::types::Ciphertext;
use toy_heaan_ckks::rings::backends::rns_ntt::RnsPoly;

use crate::contract::SensorContract;

type Ct<const DEGREE: usize> = Ciphertext<RnsPoly<DEGREE>, DEGREE>;

/// One pipeline stage as observed after it ran.
#[derive(Clone, Debug)]
pub struct Stage {
    pub name: &'static str,
    /// RNS channels left in the stage's output ciphertext.
    pub channels: usize,
    /// Scale bits of the stage's output ciphertext.
    pub logp: u32,
    pub elapsed: Duration,
}

/// Result of [`encrypted_std_dev`] plus what the demo wants to show.
pub struct StdDevRun<const DEGREE: usize> {
    /// Std in normalized units (restore with `SensorContract::restore_std`).
    pub std_normalized: Ct<DEGREE>,
    /// Intermediate encrypted variance in normalized units.
    pub variance: Ct<DEGREE>,
    pub stages: Vec<Stage>,
}

/// Records `ct`'s level and scale under `name`, timed from `started`.
fn record<const DEGREE: usize>(
    stages: &mut Vec<Stage>,
    name: &'static str,
    ct: &Ct<DEGREE>,
    started: Instant,
) {
    stages.push(Stage {
        name,
        channels: ct.c0.basis().channel_count(),
        logp: ct.logp,
        elapsed: started.elapsed(),
    });
}

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
/// `ct` must hold readings normalized with [`SensorContract::normalize`],
/// tiled over all `DEGREE / 2` slots (slot `i` = reading `i mod slots`) so
/// the rotation tree puts the full sum, and hence the mean, in every slot.
/// The result holds the std in normalized units (decrypt and map it back
/// with [`SensorContract::restore_std`]), the intermediate variance, and a
/// per-stage trace. The sqrt is fitted on
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
    evaluator: &mut LeveledEvaluator<RnsPoly<DEGREE>, DEGREE>,
    rng: &mut R,
) -> StdDevRun<DEGREE> {
    let (var_lo, var_hi) = contract.variance_interval();
    let coeffs = fit_chebyshev(f64::sqrt, var_lo, var_hi, degree);
    let mut stages = Vec::new();

    // mean = sum * 1/slots
    let started = Instant::now();
    let sum_ct = evaluator.sum_slots(ct, slots, rng);
    record(&mut stages, "sum readings", &sum_ct, started);

    let started = Instant::now();
    let mean_ct = sum_ct
        .mul_const(1.0 / slots as f64, ct.logp)
        .rescale()
        .expect("mean: rescale");
    record(&mut stages, "mean (x 1/n)", &mean_ct, started);

    // centered = input - mean, at the mean's level
    let started = Instant::now();
    let ct_aligned = ct
        .match_level(&mean_ct)
        .expect("input sits above the mean's level");
    let centered = &ct_aligned - &mean_ct;
    record(&mut stages, "center", &centered, started);

    // variance = sum of squared deviations * 1/slots
    let started = Instant::now();
    let sq_diff = evaluator.mul(&centered, &centered, rng);
    record(&mut stages, "square", &sq_diff, started);

    let started = Instant::now();
    let var_sum = evaluator.sum_slots(&sq_diff, slots, rng);
    record(&mut stages, "sum squares", &var_sum, started);

    let started = Instant::now();
    let variance = var_sum
        .mul_const(1.0 / slots as f64, sq_diff.logp)
        .rescale()
        .expect("variance: rescale");
    record(&mut stages, "variance (x 1/n)", &variance, started);

    // std_n = sqrt(variance) via encrypted Clenshaw on the contract envelope
    let started = Instant::now();
    let std_normalized = eval_chebyshev_encrypted(
        &variance, &coeffs, var_lo, var_hi, evaluator, rng,
    );
    record(&mut stages, "chebyshev sqrt", &std_normalized, started);

    StdDevRun {
        std_normalized,
        variance,
        stages,
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;
    use toy_heaan_ckks::ckks::CkksEngine;
    use toy_heaan_ckks::ckks::algorithms::chebyshev::{
        eval_chebyshev_clenshaw, to_chebyshev_domain,
    };
    use toy_heaan_ckks::ckks::encoding::ckks_encoder::CkksEncoder;
    use toy_heaan_ckks::ckks::engine::CkksParams;
    use toy_heaan_ckks::math::generate_primes;
    use toy_heaan_ckks::rings::backends::rns_ntt::{RnsBasis, RnsPoly};

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
        let ct = engine.encrypt(&pt, &pk, &mut rng);

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
            result.std_normalized.c0.basis().channel_count(),
            2,
            "13 primes - 11 levels = 2 channels left"
        );
        assert_eq!(result.stages.len(), 7, "one trace entry per stage");

        // Intermediate variance: CKKS noise only, restored to degC^2.
        let variance = std_dev * std_dev;
        let variance_got = encoder.decode(&evaluator.decrypt(&result.variance))[0]
            * contract.width().powi(2);
        assert!(
            (variance_got - variance).abs() < 1e-6,
            "encrypted variance {variance_got} vs plaintext {variance}"
        );

        let decoded = encoder.decode(&evaluator.decrypt(&result.std_normalized));

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
