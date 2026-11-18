//! Mutect2 numerics (GATK / commons-math 3.5 semantics; `FastMath.log`/`exp` as libm, checked against GATK by
//! `somatic/javacheck/Llr.java`).

/// `NaturalLogUtils.PHRED_TO_LOG_ERROR_PROB_FACTOR` = -Math.log(10)/10
const PHRED_TO_LOG_ERROR_PROB_FACTOR: f64 = -std::f64::consts::LN_10 / 10.0;

/// commons-math `FastMath.exp` (table-driven; differs from libm in the last bit for some inputs).
pub struct FastExp {
    int_a: Vec<f64>,
    int_b: Vec<f64>,
    frac_a: Vec<f64>,
    frac_b: Vec<f64>,
}

impl FastExp {
    pub fn load() -> FastExp {
        let b = include_bytes!("fastmath_exp.bin");
        let mut o = 0usize;
        let mut next = || {
            let n = u32::from_le_bytes(b[o..o + 4].try_into().unwrap()) as usize;
            o += 4;
            let v: Vec<f64> = (0..n).map(|i| f64::from_le_bytes(b[o + 8 * i..o + 8 * i + 8].try_into().unwrap())).collect();
            o += 8 * n;
            v
        };
        let (int_a, int_b, frac_a, frac_b) = (next(), next(), next(), next());
        FastExp { int_a, int_b, frac_a, frac_b }
    }

    pub fn exp(&self, x: f64) -> f64 {
        const EXP_INT_TABLE_MAX_INDEX: i32 = 750;
        if x.is_nan() {
            return f64::NAN;
        }
        let mut int_val = x as i32;
        if x < 0.0 {
            if x < -746.0 {
                return 0.0;
            }
            if int_val < -709 {
                return self.exp(x + 40.19140625) / 285040095144011776.0;
            }
            if int_val == -709 {
                return self.exp(x + 1.494140625) / 4.455505956692756620;
            }
            int_val -= 1;
        } else if int_val > 709 {
            return f64::INFINITY;
        }
        let ia = self.int_a[(EXP_INT_TABLE_MAX_INDEX + int_val) as usize];
        let ib = self.int_b[(EXP_INT_TABLE_MAX_INDEX + int_val) as usize];
        let int_frac = ((x - int_val as f64) * 1024.0) as i32;
        let fa = self.frac_a[int_frac as usize];
        let fb = self.frac_b[int_frac as usize];
        let epsilon = x - (int_val as f64 + int_frac as f64 / 1024.0);
        let mut z = 0.04168701738764507;
        z = z * epsilon + 0.1666666505023083;
        z = z * epsilon + 0.5000000000042687;
        z = z * epsilon + 1.0;
        z = z * epsilon + -3.940510424527919E-20;
        let temp_a = ia * fa;
        let temp_b = ia * fb + ib * fa + ib * fb;
        let temp_c = temp_b + temp_a;
        temp_c * z + temp_b + temp_a
    }
}

pub struct M2Math {
    pub fexp: FastExp,
    digamma: Vec<f64>,
    qual_to_error_prob: [f64; 256],
    qual_to_log_prob: [f64; 256],
}

/// commons-math `Gamma.digamma` (recursive up to 49, then the asymptotic series).
pub fn digamma(x: f64) -> f64 {
    const GAMMA: f64 = 0.577215664901532860606512090082;
    if x > 0.0 && x <= 1e-5 {
        return -GAMMA - 1.0 / x;
    }
    if x >= 49.0 {
        let inv = 1.0 / (x * x);
        return crate::fastmath::log(x) - 0.5 / x - inv * ((1.0 / 12.0) + inv * (1.0 / 120.0 - inv / 252.0));
    }
    digamma(x + 1.0) - 1.0 / x
}

/// `NaturalLogUtils.log1mexp`
pub fn log1mexp(a: f64) -> f64 {
    if a > 0.0 {
        return f64::NAN;
    }
    if a == 0.0 {
        return f64::NEG_INFINITY;
    }
    if a < crate::jmath::log(0.5) { crate::jmath::log1p(-crate::jmath::exp(a)) } else { crate::jmath::log(-crate::jmath::expm1(a)) }
}

/// commons-math `CombinatoricsUtils.binomialCoefficient` (exact for n < 67).
fn binomial_exact(n: i64, k: i64) -> i64 {
    let k = k.min(n - k);
    let mut r: u128 = 1;
    for j in 1..=k {
        r = r * (n - k + j) as u128 / j as u128;
    }
    r as i64
}

/// commons-math `CombinatoricsUtils.binomialCoefficientDouble`
fn binomial_double(n: i64, k: i64) -> f64 {
    if n == k || k == 0 {
        return 1.0;
    }
    if k == 1 || k == n - 1 {
        return n as f64;
    }
    if k > n / 2 {
        return binomial_double(n, n - k);
    }
    if n < 67 {
        return binomial_exact(n, k) as f64;
    }
    let mut r = 1.0f64;
    for i in 1..=k {
        r *= (n - k + i) as f64 / i as f64;
    }
    (r + 0.5).floor()
}

/// commons-math `CombinatoricsUtils.binomialCoefficientLog`
pub fn binomial_coefficient_log(n: i64, k: i64) -> f64 {
    if n == k || k == 0 {
        return 0.0;
    }
    if k == 1 || k == n - 1 {
        return crate::fastmath::log(n as f64);
    }
    if n < 67 {
        return crate::fastmath::log(binomial_exact(n, k) as f64);
    }
    if n < 1030 {
        return crate::fastmath::log(binomial_double(n, k));
    }
    if k > n / 2 {
        return binomial_coefficient_log(n, n - k);
    }
    let mut s = 0.0;
    for i in n - k + 1..=n {
        s += crate::fastmath::log(i as f64);
    }
    for i in 2..=k {
        s -= crate::fastmath::log(i as f64);
    }
    s
}

/// Java `Math.round(double)` (as a double): floor(x + 1/2) computed exactly; NaN -> 0.
pub fn java_round(x: f64) -> f64 {
    if x.is_nan() {
        return 0.0;
    }
    if x.is_infinite() {
        return if x > 0.0 { i64::MAX as f64 } else { i64::MIN as f64 };
    }
    let f = x.floor();
    if x - f >= 0.5 { f + 1.0 } else { f }
}

/// `QualityUtils.errorProbToQual(errorRate)` (max qual `MAX_SAM_QUAL_SCORE` = 93), as a Java byte.
pub fn error_prob_to_qual(error_rate: f64) -> i8 {
    let d = java_round(-10.0 * crate::jmath::log10(error_rate));
    let q = d as i32; // Java (int) of a double saturates, as Rust `as` does
    q.min(93).max(1) as i8
}

/// `MathUtils.fastBernoulliEntropy`
#[inline]
pub fn fast_bernoulli_entropy(p: f64) -> f64 {
    let product = p * (1.0 - p);
    product * (11.0 + 33.0 * product) / (2.0 + 20.0 * product)
}

impl M2Math {
    pub fn new() -> Self {
        let mut qe = [0f64; 256];
        let mut ql = [0f64; 256];
        for i in 0..256 {
            qe[i] = crate::jmath::pow(10.0, i as f64 / -10.0);
            ql[i] = log1mexp(i as f64 * PHRED_TO_LOG_ERROR_PROB_FACTOR);
        }
        M2Math { fexp: FastExp::load(), digamma: (0..4096).map(|i| digamma(i as f64)).collect(), qual_to_error_prob: qe, qual_to_log_prob: ql }
    }

    #[inline]
    pub fn digamma_int(&self, i: usize) -> f64 {
        if i < self.digamma.len() { self.digamma[i] } else { digamma(i as f64) }
    }

    /// `Mutect2Engine.logLikelihoodRatio(nRef, altQuals, 1)` with the flat Beta(1,1) prior. Quals are Java bytes.
    pub fn log_likelihood_ratio(&self, n_ref: i64, alt_quals: &[i8]) -> f64 {
        self.log_likelihood_ratio_rep(n_ref, alt_quals, 1)
    }

    /// `Mutect2Engine.logLikelihoodRatio(refCount, altCount, errorProbability)`
    pub fn log_likelihood_ratio_err(&self, ref_count: i64, alt_count: i64, error_probability: f64) -> f64 {
        let q = error_prob_to_qual(error_probability);
        self.log_likelihood_ratio_rep(ref_count, &[q], alt_count)
    }

    /// `Mutect2Engine.logLikelihoodRatio(nRef, altQuals, repeatFactor)` (flat prior).
    pub fn log_likelihood_ratio_rep(&self, n_ref: i64, alt_quals: &[i8], repeat_factor: i64) -> f64 {
        let n_alt = repeat_factor * alt_quals.len() as i64;
        let n = n_ref + n_alt;
        let f_tilde_ratio = self.fexp.exp(self.digamma_int((n_ref + 1) as usize) - self.digamma_int((n_alt + 1) as usize));
        let mut read_sum = 0.0;
        for &q in alt_quals {
            let qi = (q as u8) as usize;
            let eps = self.qual_to_error_prob[qi];
            let z_bar_alt = (1.0 - eps) / (1.0 - eps + eps * f_tilde_ratio);
            let log_eps = qi as f64 * PHRED_TO_LOG_ERROR_PROB_FACTOR;
            let log_one_minus_eps = self.qual_to_log_prob[qi];
            read_sum += z_bar_alt * (log_one_minus_eps - log_eps) + fast_bernoulli_entropy(z_bar_alt);
        }
        let beta_entropy = -crate::jmath::log((n + 1) as f64) - binomial_coefficient_log(n, n_alt);
        beta_entropy + read_sum * repeat_factor as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    /// GATK values from `somatic/javacheck/Llr.java` (run it to regenerate llr.txt).
    #[test]
    fn llr_matches_gatk() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../somatic/javacheck/llr.txt");
        let Ok(text) = std::fs::read_to_string(path) else { return };
        let m = M2Math::new();
        let (mut n, mut bad) = (0, 0);
        for l in text.lines() {
            let f: Vec<&str> = l.split(' ').collect();
            let n_ref: i64 = f[0].parse().unwrap();
            let q: Vec<i8> = if f[1] == "-" { vec![] } else { f[1].split(',').map(|x| x.parse().unwrap()).collect() };
            let want = f64::from_bits(u64::from_str_radix(f[2], 16).unwrap());
            let got = m.log_likelihood_ratio(n_ref, &q);
            n += 1;
            if got.to_bits() != want.to_bits() {
                if bad < 10 {
                    eprintln!("nRef {n_ref} nAlt {} got {got:e} want {want:e}", q.len());
                }
                bad += 1;
            }
        }
        assert_eq!(bad, 0, "{bad} of {n} differ");
    }
}
