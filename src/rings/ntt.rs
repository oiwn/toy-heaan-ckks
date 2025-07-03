use super::RnsResult;

#[derive(Debug)]
pub struct NttTables {
    /// Forward NTT twiddle factors per prime: roots[i][k] = ω^k mod primes[i]
    roots: Vec<Vec<u64>>,
    /// Inverse NTT twiddles per prime.
    inv_roots: Vec<Vec<u64>>,
    /// Modular inverses of the polynomial degree for each prime.
    inv_degree: Vec<u64>,
}

impl NttTables {
    pub fn build_ntt_tables_for_primes(primes: &[u64]) -> RnsResult<NttTables> {
        let capacity = primes.len();
        let mut tables = NttTables {
            roots: Vec::with_capacity(capacity),
            inv_roots: Vec::with_capacity(capacity),
            inv_degree: Vec::with_capacity(capacity),
        };

        for &p in primes {
            // Self::build_tables_for_prime(p, &mut tables)?;
        }

        Ok(tables)
    }
}
