use super::{CkksEngine, CkksError};
use crate::crypto::engine::CkksParams;
use crate::rings::backends::rns::{RnsBasis, RnsNttPoly};
use std::sync::Arc;

pub struct CkksEngineBuilder<const DEGREE: usize> {
    error_variance: Option<f64>,
    hamming_weight: Option<usize>,
    scale_bits: Option<u32>,
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
            scale_bits: None,
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

    pub fn scale_bits(mut self, scale_bits: u32) -> Self {
        self.scale_bits = Some(scale_bits);
        self
    }

    pub fn build_rns(
        self,
        rns_basis: Arc<RnsBasis>,
        scale_bits: u32,
    ) -> Result<CkksEngine<RnsNttPoly<DEGREE>, DEGREE>, CkksError> {
        let params = CkksParams {
            error_variance: self.error_variance.unwrap_or(3.2),
            hamming_weight: self.hamming_weight.unwrap_or(DEGREE / 2),
            scale_bits,
        };

        Ok(CkksEngine::<RnsNttPoly<DEGREE>, DEGREE>::new(
            rns_basis, params,
        ))
    }
}
