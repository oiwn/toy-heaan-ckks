//! Client-side sensor contract: reading range plus std envelope.

use toy_heaan_ckks::ckks::{CkksError, CkksResult};

/// Client-side contract for the sensor std-dev pipeline.
///
/// Readings are clamped to `[lo, hi]` and sent normalized as
/// `x_n = (x - lo) / width` (so every slot is in `[0, 1]`). The client also
/// promises that the plaintext std lies in `[std_min, std_max]`; the sqrt
/// approximation is fitted on exactly that envelope, in normalized units.
/// Because sqrt is homogeneous, the approximation error in reading units
/// depends only on `std_max / std_min` and the degree, not on `width`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SensorContract {
    pub lo: f64,
    pub hi: f64,
    pub std_min: f64,
    pub std_max: f64,
}

impl SensorContract {
    /// Validates and builds a contract.
    ///
    /// # Errors
    /// Returns [`CkksError::InvalidParameter`] if any bound is not finite,
    /// `lo >= hi`, `std_min <= 0`, `std_max <= std_min`, or
    /// `std_max > (hi - lo) / 2` (no data inside `[lo, hi]` can have a
    /// larger std).
    pub fn new(lo: f64, hi: f64, std_min: f64, std_max: f64) -> CkksResult<Self> {
        let invalid = |message: &str| {
            Err(CkksError::InvalidParameter {
                message: format!("SensorContract: {message}"),
            })
        };
        if ![lo, hi, std_min, std_max]
            .iter()
            .all(|value| value.is_finite())
        {
            return invalid("all bounds must be finite");
        }
        if lo >= hi {
            return invalid("lo must be below hi");
        }
        if std_min <= 0.0 {
            return invalid("std_min must be positive");
        }
        if std_max <= std_min {
            return invalid("std_max must exceed std_min");
        }
        if std_max > (hi - lo) / 2.0 {
            return invalid("std_max cannot exceed (hi - lo) / 2");
        }
        Ok(Self {
            lo,
            hi,
            std_min,
            std_max,
        })
    }

    /// The demo contract: outdoor readings in [-40, 60] degC, daily std in
    /// [1, 8] degC.
    #[must_use]
    pub fn demo() -> Self {
        Self::new(-40.0, 60.0, 1.0, 8.0).expect("demo contract is valid")
    }

    /// Width of the reading range, `hi - lo`.
    #[must_use]
    pub fn width(&self) -> f64 {
        self.hi - self.lo
    }

    /// Clamps a reading to `[lo, hi]` and maps it to `[0, 1]`.
    #[must_use]
    pub fn normalize(&self, reading: f64) -> f64 {
        (reading.clamp(self.lo, self.hi) - self.lo) / self.width()
    }

    /// Normalized variance interval the sqrt is fitted on:
    /// `((std_min / width)^2, (std_max / width)^2)`.
    #[must_use]
    pub fn variance_interval(&self) -> (f64, f64) {
        let width = self.width();
        (
            (self.std_min / width).powi(2),
            (self.std_max / width).powi(2),
        )
    }

    /// Maps a normalized std back to reading units.
    #[must_use]
    pub fn restore_std(&self, std_normalized: f64) -> f64 {
        std_normalized * self.width()
    }

    /// Whether a std (in reading units) lies inside the envelope.
    #[must_use]
    pub fn contains_std(&self, std_dev: f64) -> bool {
        (self.std_min..=self.std_max).contains(&std_dev)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use toy_heaan_ckks::ckks::algorithms::chebyshev::fit_chebyshev;

    #[test]
    fn sensor_contract_validates() {
        assert!(SensorContract::new(-40.0, 60.0, 1.0, 8.0).is_ok());
        assert!(
            SensorContract::new(60.0, -40.0, 1.0, 8.0).is_err(),
            "lo >= hi"
        );
        assert!(SensorContract::new(0.0, 0.0, 1.0, 8.0).is_err(), "lo == hi");
        assert!(
            SensorContract::new(-40.0, 60.0, 0.0, 8.0).is_err(),
            "std_min 0"
        );
        assert!(
            SensorContract::new(-40.0, 60.0, 8.0, 8.0).is_err(),
            "max <= min"
        );
        assert!(
            SensorContract::new(0.0, 10.0, 1.0, 6.0).is_err(),
            "max > W/2"
        );
        assert!(
            SensorContract::new(f64::NAN, 1.0, 0.1, 0.2).is_err(),
            "NaN lo"
        );
    }

    #[test]
    fn sensor_contract_normalize_and_restore() {
        let contract = SensorContract::demo();
        assert_eq!(contract.width(), 100.0);
        assert_eq!(contract.normalize(-50.0), 0.0, "clamped to lo");
        assert_eq!(contract.normalize(70.0), 1.0, "clamped to hi");
        assert!((contract.normalize(10.0) - 0.5).abs() < 1e-15);
        assert!((contract.restore_std(3.0 / 100.0) - 3.0).abs() < 1e-12);

        let (var_lo, var_hi) = contract.variance_interval();
        assert!((var_lo - 1e-4).abs() < 1e-18);
        assert!((var_hi - 6.4e-3).abs() < 1e-16);

        assert!(contract.contains_std(3.0));
        assert!(!contract.contains_std(0.5));
        assert!(!contract.contains_std(9.0));
    }

    /// Small fitted coefficients keep every encoded constant far from any
    /// modulus limit.
    #[test]
    fn demo_fit_coefficients_are_small() {
        let (var_lo, var_hi) = SensorContract::demo().variance_interval();
        for coeff in fit_chebyshev(f64::sqrt, var_lo, var_hi, 7) {
            assert!(coeff.abs() < 0.1, "coefficient {coeff} too large");
        }
    }
}
