use super::{PolyRing, PolySampler};

#[derive(Clone, Debug)]
pub struct RnsNttPoly<const DEGREE: usize> {
    pub(crate) coefficients: Vec<[u64; DEGREE]>,
    pub(crate) in_ntt_domain: bool,
    pub(crate) basis: Arc<RnsBasis>,
}
