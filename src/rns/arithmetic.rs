use super::{RnsElement, RnsError, RnsResult};

pub trait RnsOps {
    fn add(&self, other: &Self) -> RnsResult<Self>;
    fn sub(&self, other: &Self) -> RnsResult<Self>;
    fn neg(&self) -> Self;
    fn mul_scalar(&self, scalar: u64) -> Self;
}

impl RnsOps for RnsElement {
    fn add(&self, other: &Self) -> RnsResult<Self> {
        ensure_compatible_basis(self, other)?;

        let residues = self
            .residues
            .iter()
            .zip(&other.residues)
            .zip(self.basis.primes())
            .map(|((&a, &b), &p)| (a + b) % p)
            .collect();

        Ok(RnsElement::new(residues, self.basis.clone())?)
    }

    fn sub(&self, other: &Self) -> RnsResult<Self> {
        ensure_compatible_basis(self, other)?;

        let residues = self
            .residues
            .iter()
            .zip(&other.residues)
            .zip(self.basis.primes())
            .map(|((&a, &b), &p)| (a + p - b) % p) // Add p to avoid underflow
            .collect();

        Ok(RnsElement::new(residues, self.basis.clone())?)
    }

    fn neg(&self) -> Self {
        let residues = self
            .residues
            .iter()
            .zip(self.basis.primes())
            .map(|(&a, &p)| if a == 0 { 0 } else { p - a })
            .collect();

        RnsElement::new(residues, self.basis.clone()).unwrap()
    }

    fn mul_scalar(&self, scalar: u64) -> Self {
        let residues = self
            .residues
            .iter()
            .zip(self.basis.primes())
            .map(|(&a, &p)| (a * (scalar % p)) % p)
            .collect();

        RnsElement::new(residues, self.basis.clone()).unwrap()
    }
}

// Helper function
fn ensure_compatible_basis(a: &RnsElement, b: &RnsElement) -> RnsResult<()> {
    if a.basis.primes() != b.basis.primes() {
        return Err(RnsError::IncompatibleBases);
    }
    Ok(())
}
