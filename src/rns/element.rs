//! Core RnsElement type and basic operations

use super::{RnsBasis, RnsError, RnsResult};
use std::sync::Arc;

/// An integer represented in Residue Number System
#[derive(Debug, Clone)]
pub struct RnsElement {
    residues: Vec<u64>,
    basis: Arc<RnsBasis>,
}

impl RnsElement {
    pub fn new(residues: Vec<u64>, basis: Arc<RnsBasis>) -> RnsResult<Self> {
        if residues.len() != basis.prime_count() {
            return Err(RnsError::IncompatibleBases);
        }
        Ok(Self { residues, basis })
    }

    pub fn from_u64(value: u64, basis: Arc<RnsBasis>) -> Self {
        let residues = basis.primes().iter().map(|&prime| value % prime).collect();

        Self { residues, basis }
    }

    pub fn residues(&self) -> &[u64] {
        &self.residues
    }

    pub fn basis(&self) -> &Arc<RnsBasis> {
        &self.basis
    }

    pub fn prime_count(&self) -> usize {
        self.residues.len()
    }
}
