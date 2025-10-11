use super::{CkksEngine, CkksError};
use crate::crypto::engine::CkksParams;
use crate::rings::NaivePolyRing;
use crate::rings::backends::bigint::BigIntContext;
use crate::rings::backends::bigint::BigIntPolyRing;
use crate::rings::backends::rns::{RnsBasis, RnsPolyRing};
use std::sync::Arc;

pub enum CkksEngineVariant<const DEGREE: usize> {
    Naive(CkksEngine<NaivePolyRing<DEGREE>, DEGREE>),
    BigIntU256(CkksEngine<BigIntPolyRing<DEGREE>, DEGREE>),
    RNS(CkksEngine<RnsPolyRing<DEGREE>, DEGREE>),
}

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

    pub fn build_naive(
        self,
        modulus: u64,
        scale_bits: u32,
    ) -> Result<CkksEngine<NaivePolyRing<DEGREE>, DEGREE>, CkksError> {
        let params = CkksParams {
            error_variance: self.error_variance.unwrap_or(3.2),
            hamming_weight: self.hamming_weight.unwrap_or(DEGREE / 2),
            scale_bits,
        };

        Ok(CkksEngine::<NaivePolyRing<DEGREE>, DEGREE>::new(
            modulus, params,
        ))
    }

    pub fn build_bigint_u256(
        self,
        modulus: crypto_bigint::NonZero<crypto_bigint::U256>,
        scale_bits: u32,
    ) -> Result<CkksEngine<BigIntPolyRing<DEGREE>, DEGREE>, CkksError> {
        let params = CkksParams {
            error_variance: self.error_variance.unwrap_or(3.2),
            hamming_weight: self.hamming_weight.unwrap_or(DEGREE / 2),
            scale_bits,
        };

        // Create BigIntContext with Kim's HEAAN parameters for proper multiplication
        // Use the same parameters as examples/ckks_full.rs for consistency
        let context = BigIntContext::new(modulus, 50, 20);

        Ok(CkksEngine::<BigIntPolyRing<DEGREE>, DEGREE>::new(
            context, params,
        ))
    }

    pub fn build_rns(
        self,
        rns_basis: Arc<RnsBasis>,
        scale_bits: u32,
    ) -> Result<CkksEngine<RnsPolyRing<DEGREE>, DEGREE>, CkksError> {
        let params = CkksParams {
            error_variance: self.error_variance.unwrap_or(3.2),
            hamming_weight: self.hamming_weight.unwrap_or(DEGREE / 2),
            scale_bits,
        };

        Ok(CkksEngine::<RnsPolyRing<DEGREE>, DEGREE>::new(
            rns_basis, params,
        ))
    }
}
