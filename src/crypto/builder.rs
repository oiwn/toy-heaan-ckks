use super::{CkksEngine, CkksError};
use crate::crypto::engine::CkksParams;
use crate::encoding::{Encoder, EncoderType, RustFftEncoder};
use crate::rings::NaivePolyRing;
use crate::rings::backends::PolyRingU256;

pub enum CkksEngineVariant<const DEGREE: usize> {
    Naive(CkksEngine<NaivePolyRing<DEGREE>, DEGREE>),
    BigIntU256(CkksEngine<PolyRingU256<DEGREE>, DEGREE>),
}

pub struct CkksEngineBuilder<const DEGREE: usize> {
    encoder_type: Option<EncoderType>,
    error_variance: Option<f64>,
    hamming_weight: Option<usize>,
    scale_bits: Option<u32>,
}

impl<const DEGREE: usize> CkksEngineBuilder<DEGREE> {
    pub fn new() -> Self {
        Self {
            encoder_type: None,
            error_variance: None,
            hamming_weight: None,
            scale_bits: None,
        }
    }

    pub fn encoder(mut self, encoder_type: EncoderType) -> Self {
        self.encoder_type = Some(encoder_type);
        self
    }

    pub fn error_variance(mut self, variance: f64) -> Self {
        self.error_variance = Some(variance);
        self
    }

    pub fn hamming_weight(mut self, weight: usize) -> Self {
        self.hamming_weight = Some(weight);
        self
    }

    pub fn scale_bits(mut self, bits: u32) -> Self {
        self.scale_bits = Some(bits);
        self
    }

    pub fn build_naive(
        self,
        modulus: u64,
    ) -> Result<CkksEngine<NaivePolyRing<DEGREE>, DEGREE>, CkksError> {
        let encoder: Box<dyn Encoder<DEGREE>> = match self.encoder_type {
            Some(EncoderType::RustFft) => Box::new(RustFftEncoder::<DEGREE>::new(
                self.scale_bits.expect("Use .scale_bits(...)"),
            )?),
            _ => panic!("Only RustFft encoder supported currently"),
        };

        let params = CkksParams {
            error_variance: self.error_variance.unwrap_or(3.2),
            hamming_weight: self.hamming_weight.unwrap_or(DEGREE / 2),
            scale_bits: self.scale_bits.unwrap_or(40),
        };

        Ok(CkksEngine::<NaivePolyRing<DEGREE>, DEGREE>::new(
            modulus, encoder, params,
        ))
    }

    pub fn build_bigint_u256(
        self,
        modulus: crypto_bigint::NonZero<crypto_bigint::U256>,
    ) -> Result<CkksEngine<PolyRingU256<DEGREE>, DEGREE>, CkksError> {
        let encoder: Box<dyn Encoder<DEGREE>> = match self.encoder_type {
            Some(EncoderType::RustFft) => Box::new(RustFftEncoder::<DEGREE>::new(
                self.scale_bits.expect("Use .scale_bits(...)"),
            )?),
            _ => panic!("Only RustFft encoder supported currently"),
        };

        let params = CkksParams {
            error_variance: self.error_variance.unwrap_or(3.2),
            hamming_weight: self.hamming_weight.unwrap_or(DEGREE / 2),
            scale_bits: self.scale_bits.unwrap_or(40),
        };

        Ok(CkksEngine::<PolyRingU256<DEGREE>, DEGREE>::new(
            modulus, encoder, params,
        ))
    }
}
