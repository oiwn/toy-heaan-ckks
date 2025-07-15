use super::{CkksEngine, CkksResult};
use crate::crypto::engine::CkksParams;
use crate::encoding::{Encoder, EncoderType, RustFftEncoder};
use crate::rings::{BackendType, NaivePolyRing};

pub struct CkksEngineBuilder<const DEGREE: usize> {
    encoder_type: Option<EncoderType>,
    backend_type: Option<BackendType>,
    error_variance: Option<f64>,
    hamming_weight: Option<usize>,
    scale_bits: Option<u32>,
}

impl<const DEGREE: usize> CkksEngineBuilder<DEGREE> {
    pub fn new() -> Self {
        Self {
            encoder_type: None,
            backend_type: None,
            error_variance: None,
            hamming_weight: None,
            scale_bits: None,
        }
    }

    pub fn encoder(mut self, encoder_type: EncoderType) -> Self {
        self.encoder_type = Some(encoder_type);
        self
    }

    pub fn backend(mut self, backend: BackendType) -> Self {
        self.backend_type = Some(backend);
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

    pub fn build(self) -> CkksResult<CkksEngine<NaivePolyRing<DEGREE>, DEGREE>> {
        let encoder: Box<dyn Encoder<DEGREE>> = match self.encoder_type {
            Some(EncoderType::RustFft) => Box::new(RustFftEncoder::<DEGREE>::new(
                self.scale_bits.expect("Use .scale_bits(...)"),
            )?),
            _ => panic!(),
        };

        let backend_type = self.backend_type.unwrap_or(BackendType::Naive(23));

        match backend_type {
            BackendType::Naive(modulus) => {
                let params = CkksParams {
                    error_variance: self.error_variance.unwrap_or(3.2),
                    hamming_weight: self.hamming_weight.unwrap_or(DEGREE / 2),
                    scale_bits: self.scale_bits.unwrap_or(40),
                };

                Ok(CkksEngine::<NaivePolyRing<DEGREE>, DEGREE>::new(
                    modulus, encoder, params,
                ))
            } // _ => Err(CkksError::InvalidParameter {
              //     message: "Only Naive backend supported currently".to_string(),
              // }),
        }
    }
}
