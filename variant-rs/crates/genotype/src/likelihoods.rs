//! Pieces of `AlleleLikelihoods` used after the likelihood calculation:
//! best-allele search (with tie-breaking priorities), marginalisation and the informative test.

/// `AlleleLikelihoods.LOG_10_INFORMATIVE_THRESHOLD`
pub const LOG_10_INFORMATIVE_THRESHOLD: f64 = 0.2;

#[derive(Clone, Copy, Debug)]
pub struct BestAllele {
    pub allele: Option<usize>,
    pub likelihood: f64,
    pub second: Option<usize>,
    pub second_likelihood: f64,
    pub confidence: f64,
}

impl BestAllele {
    #[inline]
    pub fn is_informative(&self) -> bool {
        self.confidence > LOG_10_INFORMATIVE_THRESHOLD
    }
}

/// `AlleleLikelihoods.searchBestAllele(sample, evidence, canBeReference = true, priorities)`.
/// `lks[a]` is the likelihood of the read under allele `a`.
pub fn best_allele(lks: &[f64], priorities: Option<&[f64]>) -> BestAllele {
    best_allele_t(lks, priorities, LOG_10_INFORMATIVE_THRESHOLD)
}

/// [`best_allele`] with the tie threshold `getInformativeThreshold()` (0.2, or 0.2 * ln 10 after
/// `switchToNaturalLog`).
pub fn best_allele_t(lks: &[f64], priorities: Option<&[f64]>, threshold: f64) -> BestAllele {
    let n = lks.len();
    if n == 0 {
        return BestAllele {
            allele: None,
            likelihood: f64::NEG_INFINITY,
            second: None,
            second_likelihood: f64::NEG_INFINITY,
            confidence: 0.0,
        };
    }
    let mut best = 0usize;
    let mut second = 0usize;
    let mut best_lk = lks[0];
    let mut second_lk = f64::NEG_INFINITY;
    for (a, &c) in lks.iter().enumerate().skip(1) {
        if c > best_lk {
            second = best;
            best = a;
            second_lk = best_lk;
            best_lk = c;
        } else if c > second_lk {
            second = a;
            second_lk = c;
        }
    }
    if let Some(pr) = priorities {
        if best_lk - second_lk < threshold {
            let mut best_pr = pr[best];
            let mut second_pr = pr[second];
            for a in 0..n {
                let c = lks[a];
                if a == best || best_lk - c > threshold {
                    continue;
                }
                let cp = pr[a];
                if cp > best_pr {
                    second = best;
                    best = a;
                    second_pr = best_pr;
                    best_pr = cp;
                } else if cp > second_pr {
                    second = a;
                    second_pr = cp;
                }
            }
        }
    }
    let likelihood = lks[best];
    let second_likelihood = if second != best { lks[second] } else { f64::NEG_INFINITY };
    let confidence = if likelihood == second_likelihood { 0.0 } else { likelihood - second_likelihood };
    BestAllele { allele: Some(best), likelihood, second: Some(second), second_likelihood, confidence }
}

/// Default tie-breaking of `bestAllelesBreakingTies()`: reference allele gets priority 1.
pub fn ref_priorities(n_alleles: usize, ref_index: Option<usize>) -> Vec<f64> {
    (0..n_alleles).map(|a| if Some(a) == ref_index { 1.0 } else { 0.0 }).collect()
}
