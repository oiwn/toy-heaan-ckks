use crate::PolyRing;

pub struct Plaintext<P, const DEGREE: usize>
where
    P: PolyRing<DEGREE>,
{
    pub poly: P,
    pub scale: f64,
}
