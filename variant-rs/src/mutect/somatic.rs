//! `SomaticLikelihoodsEngine` and the numeric helpers of `SomaticGenotypingEngine`, with Java's
//! evaluation order (plain left-to-right sums, `DoubleStream.sum` compensated summation where GATK uses
//! streams). `Gamma.logGamma` / `digamma` are the exact commons-math ports; `Math.exp` / `Math.log` use
//! libm (HotSpot's intrinsics can differ in the last bit).
//!
//! Matrices are `[allele][evidence]` natural-log likelihoods.

use hc_assembly::fastmath;
use hc_assembly::m2math::digamma;

pub const CONVERGENCE_THRESHOLD: f64 = 0.001;
const NEGLIGIBLE_RESPONSIBILITY: f64 = 1.0e-10;
/// `NaturalLogUtils.LOG_ONE_HALF = Math.log(0.5)`
pub fn log_one_half() -> f64 {
    hc_assembly::hotspot::log(0.5)
}

/// `MathUtils.sum(double[])`
#[inline]
pub fn sum(v: &[f64]) -> f64 {
    let mut s = 0.0;
    for &x in v {
        s += x;
    }
    s
}

/// `DoubleStream.sum()` (Kahan summation as in `Collectors.sumWithCompensation` / `computeFinalSum`, JDK 17).
pub fn stream_sum(v: impl IntoIterator<Item = f64>) -> f64 {
    let (mut hi, mut comp, mut simple) = (0.0f64, 0.0f64, 0.0f64);
    for d in v {
        let tmp = d - comp;
        let velvel = hi + tmp;
        comp = (velvel - hi) - tmp;
        hi = velvel;
        simple += d;
    }
    let tmp = hi - comp;
    if tmp.is_nan() && simple.is_infinite() {
        simple
    } else {
        tmp
    }
}

/// `MathUtils.maxElementIndex`: first index of the maximum.
#[inline]
fn max_element_index(v: &[f64]) -> usize {
    let mut m = 0;
    for i in 1..v.len() {
        if v[i] > v[m] {
            m = i;
        }
    }
    m
}

/// `NaturalLogUtils.logSumExp`
pub fn log_sum_exp(v: &[f64]) -> f64 {
    let mi = max_element_index(v);
    let max = v[mi];
    if max == f64::NEG_INFINITY {
        return max;
    }
    let mut s = 1.0;
    for (i, &c) in v.iter().enumerate() {
        if i == mi || c == f64::NEG_INFINITY {
            continue;
        }
        s += hc_assembly::hotspot::exp(c - max);
    }
    assert!(!(s.is_nan() || s == f64::INFINITY), "logValues must be non-infinite and non-NAN");
    max + if s != 1.0 { hc_assembly::hotspot::log(s) } else { 0.0 }
}

/// `NaturalLogUtils.normalizeFromLogToLinearSpace` (in place)
pub fn normalize_from_log_to_linear(v: &mut [f64]) {
    let ls = log_sum_exp(v);
    for x in v.iter_mut() {
        *x -= ls;
    }
    for x in v.iter_mut() {
        *x = hc_assembly::hotspot::exp(*x);
    }
}

/// `NaturalLogUtils.posteriors(logPriors, logLikelihoods)`
fn posteriors(log_priors: &[f64], lks: &[f64]) -> Vec<f64> {
    let mut v: Vec<f64> = log_priors.iter().zip(lks).map(|(a, b)| a + b).collect();
    normalize_from_log_to_linear(&mut v);
    v
}

/// `Dirichlet.effectiveLogMultinomialWeights`
fn effective_log_weights(alpha: &[f64]) -> Vec<f64> {
    let ds = digamma(sum(alpha));
    alpha.iter().map(|&a| digamma(a) - ds).collect()
}

#[inline]
fn column(m: &[Vec<f64>], r: usize, out: &mut Vec<f64>) {
    out.clear();
    out.extend(m.iter().map(|row| row[r]));
}

/// `SomaticLikelihoodsEngine.getEffectiveCounts(logLikelihoods, dirichletPrior, null)` via
/// `MathUtils.sumArrayFunction` (first column, then += the others). Needs at least one read.
fn effective_counts(m: &[Vec<f64>], n_reads: usize, alpha: &[f64]) -> Vec<f64> {
    let w = effective_log_weights(alpha);
    let mut col = Vec::with_capacity(m.len());
    column(m, 0, &mut col);
    let mut res = posteriors(&w, &col);
    for r in 1..n_reads {
        column(m, r, &mut col);
        let p = posteriors(&w, &col);
        for (x, y) in res.iter_mut().zip(p) {
            *x += y;
        }
    }
    res
}

/// `SomaticLikelihoodsEngine.alleleFractionsPosterior(logLikelihoods, priorPseudocounts)`
pub fn allele_fractions_posterior(m: &[Vec<f64>], n_reads: usize, prior: &[f64]) -> Vec<f64> {
    let mut post = vec![1.0f64; m.len()];
    loop {
        let counts = effective_counts(m, n_reads, &post);
        let new: Vec<f64> = counts.iter().zip(prior).map(|(a, b)| a + b).collect();
        // MathArrays.distance1
        let mut d = 0.0;
        for (a, b) in post.iter().zip(&new) {
            d += (a - b).abs();
        }
        let converged = d / sum(&new) < CONVERGENCE_THRESHOLD;
        post = new;
        if converged {
            return post;
        }
    }
}

/// `SomaticLikelihoodsEngine.logDirichletNormalization`
pub fn log_dirichlet_normalization(params: &[f64]) -> f64 {
    let num = fastmath::log_gamma(sum(params));
    let lg: Vec<f64> = params.iter().map(|&a| fastmath::log_gamma(a)).collect();
    num - sum(&lg)
}

#[inline]
fn x_log_x(x: f64) -> f64 {
    if x < 1e-8 { 0.0 } else { x * hc_assembly::hotspot::log(x) }
}

/// `SomaticLikelihoodsEngine.logEvidence(logLikelihoods, priorPseudocounts)`
pub fn log_evidence(m: &[Vec<f64>], n_reads: usize, prior: &[f64]) -> f64 {
    let post = allele_fractions_posterior(m, n_reads, prior);
    let prior_c = log_dirichlet_normalization(prior);
    let post_c = -log_dirichlet_normalization(&post);
    let log_af = effective_log_weights(&post);
    let mut lae = 0.0;
    let mut col = Vec::with_capacity(m.len());
    for r in 0..n_reads {
        column(m, r, &mut col);
        let resp = posteriors(&log_af, &col);
        let entropy = stream_sum(resp.iter().map(|&x| x_log_x(x)));
        let mut lc = 0.0;
        for n in 0..col.len() {
            lc += if resp[n] < NEGLIGIBLE_RESPONSIBILITY { 0.0 } else { col[n] * resp[n] };
        }
        lae += lc - entropy;
    }
    prior_c + post_c + lae
}

/// `SomaticGenotypingEngine.somaticLogOdds` with pseudocounts 1 (minAF 0): for each allele other than the
/// reference (index 0), logEvidence(all alleles) - logEvidence(without it). Entry 0 is unused.
pub fn somatic_log_odds(m: &[Vec<f64>], n_reads: usize) -> Vec<f64> {
    let n = m.len();
    let mut out = vec![0.0; n];
    let all = if n_reads == 0 { 0.0 } else { log_evidence(m, n_reads, &vec![1.0; n]) };
    for a in 1..n {
        let sub: Vec<Vec<f64>> = (0..n).filter(|&k| k != a).map(|k| m[k].clone()).collect();
        let without = if n_reads == 0 { 0.0 } else { log_evidence(&sub, n_reads, &vec![1.0; n - 1]) };
        out[a] = all - without;
    }
    out
}

/// `SomaticGenotypingEngine.diploidAltLogOdds` (reference = index 0). Entry 0 is unused.
pub fn diploid_alt_log_odds(m: &[Vec<f64>], n_reads: usize) -> Vec<f64> {
    let n = m.len();
    let mut hom_ref = 0.0;
    for r in 0..n_reads {
        hom_ref += m[0][r];
    }
    let half = log_one_half();
    let mut out = vec![0.0; n];
    for a in 1..n {
        let mut het = 0.0;
        for r in 0..n_reads {
            het += log_sum_exp(&[m[0][r], m[a][r]]) + half;
        }
        out[a] = hom_ref - het;
    }
    out
}

/// `SomaticGenotypingEngine.getEffectiveCounts` (AD before annotation; unused for output) and the AF:
/// `normalizeSumToOne(alleleFractionsPosterior(flat prior))`.
pub fn allele_fractions_mean(m: &[Vec<f64>], n_reads: usize) -> Vec<f64> {
    let n = m.len();
    let post = if n_reads == 0 { vec![1.0; n] } else { allele_fractions_posterior(m, n_reads, &vec![1.0; n]) };
    let s = sum(&post);
    post.iter().map(|&x| x / s).collect()
}

/// Java `String.hashCode` (UTF-16 code units; names are ASCII).
pub fn java_string_hash(s: &str) -> i32 {
    let mut h: i32 = 0;
    for c in s.encode_utf16() {
        h = h.wrapping_mul(31).wrapping_add(c as i32);
    }
    h
}

/// Iteration order of a `java.util.HashMap` (default capacity, load factor 0.75) filled with `keys` in
/// order (duplicates collapse to the first insertion): returns the distinct keys' first-insertion indices
/// in `values()` order. Bins are assumed not to be treeified (>= 8 colliding keys in one bin).
pub fn java_hashmap_order(keys: &[&str]) -> Vec<usize> {
    let mut first: Vec<(usize, u32)> = Vec::new(); // (insertion index, spread hash)
    let mut seen: std::collections::HashMap<&str, ()> = std::collections::HashMap::new();
    for (i, k) in keys.iter().enumerate() {
        if seen.insert(k, ()).is_none() {
            let h = java_string_hash(k) as u32;
            first.push((i, h ^ (h >> 16)));
        }
    }
    let n = first.len();
    let mut cap: u64 = 16;
    while n as u64 > cap * 3 / 4 {
        cap *= 2;
    }
    let mask = (cap - 1) as u32;
    // stable sort by bucket keeps insertion order within a bucket (resize splits preserve it)
    let mut order: Vec<(u32, usize)> = first.iter().map(|&(i, h)| (h & mask, i)).collect();
    order.sort_by_key(|&(b, _)| b);
    order.into_iter().map(|(_, i)| i).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn hashmap_order_small() {
        // "a"=97, "b"=98: buckets 1 and 2; "q" = 113 -> bucket 1 (113 & 15 = 1) after "a"
        assert_eq!(java_hashmap_order(&["b", "q", "a", "b"]), vec![1, 2, 0]);
        assert_eq!(java_string_hash("hello"), 99162322);
    }
}
