use crate::PolyRing;

#[derive(Debug)]
pub struct Plaintext<P, const DEGREE: usize>
where
    P: PolyRing<DEGREE>,
{
    pub poly: P,
    pub scale_bits: u32,
}

#[derive(Debug)]
pub struct Ciphertext<P, const DEGREE: usize>
where
    P: PolyRing<DEGREE>,
{
    pub c0: P,
    pub c1: P,
    pub scale_bits: u32,
}
