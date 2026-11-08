use super::{CkksEngine, CkksError};
use crate::PolyRing;
use crate::ckks::Scale;
use crate::ckks::engine::CkksParams;

pub struct CkksEngineBuilder<const DEGREE: usize> {
    error_variance: Option<f64>,
    hamming_weight: Option<usize>,
    scale: Option<Scale>,
}

impl<const DEGREE: usize> Default for CkksEngineBuilder<DEGREE> {
    fn default() -> Self {
        Self::new()
    }
}

impl<const DEGREE: usize> CkksEngineBuilder<DEGREE> {
    pub fn new() -> Self {
        Self {
            error_variance: None,
            hamming_weight: None,
            scale: None,
        }
    }

    pub fn error_variance(mut self, variance: f64) -> Self {
        self.error_variance = Some(variance);
        self
    }

    pub fn hamming_weight(mut self, weight: usize) -> Self {
        self.hamming_weight = Some(weight);
        self
    }

    pub fn scale(mut self, scale: Scale) -> Self {
        self.scale = Some(scale);
        self
    }

    pub fn build<P: PolyRing<DEGREE>>(
        self,
        context: P::Context,
        scale: Scale,
    ) -> Result<CkksEngine<P, DEGREE>, CkksError> {
        let params = CkksParams {
            error_variance: self.error_variance.unwrap_or(3.2),
            hamming_weight: self.hamming_weight.unwrap_or(DEGREE / 2),
            scale,
        };

        Ok(CkksEngine::new(context, params))
    }
}
