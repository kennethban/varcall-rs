//! Diploid genotype likelihoods (`GenotypeLikelihoodCalculator`), the allele frequency
//! calculator (`AlleleFrequencyCalculator`), allele subsetting (`AlleleSubsettingUtils`) and
//! genotype assignment (`GATKVariantContextUtils.makeGenotypeCall`).

use crate::javamath::java_round;
use crate::mathutils::{self, approximate_log10_sum_log10, log10_sum_log10, log10_sum_log10_2, LOG10_2};

pub const PLOIDY: usize = 2;

/// Number of diploid genotypes for `n` alleles.
#[inline]
pub fn genotype_count(n: usize) -> usize {
    n * (n + 1) / 2
}

/// Diploid genotype index of alleles (a, b), canonical VCF order.
#[inline]
pub fn genotype_index(a: usize, b: usize) -> usize {
    let (i, j) = if a <= b { (a, b) } else { (b, a) };
    j * (j + 1) / 2 + i
}

/// Alleles (i <= j) of each diploid genotype index, in index order.
pub fn genotype_alleles(n: usize) -> Vec<(usize, usize)> {
    let mut v = Vec::with_capacity(genotype_count(n));
    for j in 0..n {
        for i in 0..=j {
            v.push((i, j));
        }
    }
    v
}

/// `GenotypeLikelihoodCalculator.computeLog10GenotypeLikelihoods` for ploidy 2.
/// `lk[a][r]` are allele log10 likelihoods of the (retained) reads.
pub fn log10_genotype_likelihoods(lk: &[Vec<f64>], read_count: usize) -> Vec<f64> {
    let n = lk.len();
    let table = mathutils::jacobian();
    let mut out = Vec::with_capacity(genotype_count(n));
    let log10_ploidy_total = read_count as f64 * LOG10_2;
    for (i, j) in genotype_alleles(n) {
        if i == j {
            out.push(mathutils::sum(&lk[i][..read_count]));
        } else {
            // log10Count1 = log10Count2 = Math.log10(1) = 0.0
            let (l1, l2) = (&lk[i], &lk[j]);
            let mut s = 0.0f64;
            for r in 0..read_count {
                s += approximate_log10_sum_log10(l1[r] + 0.0, l2[r] + 0.0, table);
            }
            out.push(s - log10_ploidy_total);
        }
    }
    out
}

/// `GenotypeLikelihoods.fromLog10Likelihoods(gls).getAsPLs()`
pub fn gls_to_pls(gls: &[f64]) -> Vec<i32> {
    let mut max = f64::NEG_INFINITY;
    for &g in gls {
        max = max.max(g);
    }
    gls.iter().map(|&g| java_round((-10.0 * (g - max)).min(2147483647.0)) as i32).collect()
}

/// `GenotypeLikelihoods.getGQLog10FromLikelihoods`
pub fn gq_log10_from_likelihoods(idx: usize, lks: &[f64]) -> f64 {
    let mut qual = f64::NEG_INFINITY;
    for (i, &l) in lks.iter().enumerate() {
        if i == idx {
            continue;
        }
        if l >= qual {
            qual = l;
        }
    }
    qual = lks[idx] - qual;
    if qual < 0.0 {
        // htsjdk GeneralUtils.normalizeFromLog10 (not reached for best genotypes)
        let max = lks.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let lin: Vec<f64> = lks.iter().map(|l| hc_assembly::jmath::pow(10.0, l - max)).collect();
        let s: f64 = lin.iter().sum();
        hc_assembly::jmath::log10(1.0 - lin[idx] / s)
    } else {
        -1.0 * qual
    }
}

/// Result of `AlleleFrequencyCalculator.calculate` for one diploid sample.
#[derive(Debug, Clone)]
pub struct AfResult {
    /// MLE allele counts of the alt alleles (index 0 = allele 1)
    pub mle_alt_counts: Vec<i32>,
    pub log10_p_no_variant: f64,
    /// log10 posterior of each allele being absent; index = allele index (entry 0 unused)
    pub log10_p_ref_by_allele: Vec<f64>,
}

pub struct AfCalculator {
    ref_pseudo: f64,
    snp_pseudo: f64,
    indel_pseudo: f64,
}

impl Default for AfCalculator {
    fn default() -> Self {
        // heterozygosity 0.001, indel heterozygosity 1.25E-4, heterozygosity stdev 0.01
        let snp_het = 0.001f64;
        let ref_pseudo = snp_het / hc_assembly::jmath::pow(0.01, 2.0);
        AfCalculator { ref_pseudo, snp_pseudo: snp_het * ref_pseudo, indel_pseudo: 1.25E-4 * ref_pseudo }
    }
}

/// `GenotypeAlleleCounts.log10CombinationCount()` for diploid genotypes.
#[inline]
fn log10_comb(i: usize, j: usize) -> f64 {
    if i == j {
        0.0
    } else {
        LOG10_2
    }
}

fn log10_normalized_posteriors(gls: &[f64], log10_af: &[f64], gts: &[(usize, usize)]) -> Vec<f64> {
    let mut post = Vec::with_capacity(gts.len());
    for (g, &(i, j)) in gts.iter().enumerate() {
        let s = if i == j { 0.0 + 2.0 * log10_af[i] } else { 0.0 + 1.0 * log10_af[i] + 1.0 * log10_af[j] };
        post.push(log10_comb(i, j) + gls[g] + s);
    }
    mathutils::normalize_log10_in_place(&mut post);
    post
}

impl AfCalculator {
    /// `AlleleFrequencyCalculator.calculate(numAlleles, alleles, genotypes, ploidy=2, refLength)` for a single
    /// sample whose genotype carries `pls`. `allele_lens[a]` is each allele's length, `ref_len` the
    /// reference allele's, `span_del_index` the index of `*` if present.
    pub fn calculate(&self, pls: &[i32], allele_lens: &[usize], span_del_index: Option<usize>) -> AfResult {
        self.calculate_multi(&[pls], allele_lens, span_del_index)
    }

    /// `AlleleFrequencyCalculator.calculate` over several diploid samples (`pls[s]`: sample s's PLs, in the
    /// VariantContext's genotype order; every sample is usable, as all carry likelihoods). One sample gives
    /// exactly [`AfCalculator::calculate`].
    pub fn calculate_multi(&self, pls: &[&[i32]], allele_lens: &[usize], span_del_index: Option<usize>) -> AfResult {
        let n = allele_lens.len();
        let ref_len = allele_lens[0];
        let gls: Vec<Vec<f64>> = pls.iter().map(|p| p.iter().map(|&p| p as f64 / -10.0).collect()).collect();
        let gts = genotype_alleles(n);
        let prior: Vec<f64> = (0..n)
            .map(|a| {
                if a == 0 {
                    self.ref_pseudo
                } else if allele_lens[a] == ref_len {
                    self.snp_pseudo
                } else {
                    self.indel_pseudo
                }
            })
            .collect();
        let mut allele_counts = vec![0.0f64; n];
        let flat = -hc_assembly::jmath::log10(n as f64);
        let mut log10_af = vec![flat; n];
        let log10_count = |c: usize| if c == 1 { 0.0 } else { LOG10_2 };
        loop {
            // effectiveAlleleCounts
            let mut log10_res = vec![f64::NEG_INFINITY; n];
            for g_s in &gls {
                let post = log10_normalized_posteriors(g_s, &log10_af, &gts);
                for (g, &(i, j)) in gts.iter().enumerate() {
                    if i == j {
                        log10_res[i] = log10_sum_log10_2(log10_res[i], post[g] + log10_count(2));
                    } else {
                        log10_res[i] = log10_sum_log10_2(log10_res[i], post[g] + log10_count(1));
                        log10_res[j] = log10_sum_log10_2(log10_res[j], post[g] + log10_count(1));
                    }
                }
            }
            let new_counts: Vec<f64> = log10_res.iter().map(|&x| 10f64.powf(x)).collect();
            let mut max_diff = f64::NEG_INFINITY;
            for a in 0..n {
                max_diff = max_diff.max((allele_counts[a] - new_counts[a]).abs());
            }
            allele_counts = new_counts;
            let posterior_pseudo: Vec<f64> = (0..n).map(|a| prior[a] + allele_counts[a]).collect();
            let s = mathutils::sum(&posterior_pseudo);
            log10_af = posterior_pseudo.iter().map(|&x| (x / s).log10()).collect();
            if !(max_diff > 0.1) {
                break;
            }
        }

        let mut log10_p_of_zero = vec![0.0f64; n];
        let mut log10_p_no_variant = 0.0f64;
        let biallelic = n == 2 && span_del_index.is_none();
        for g_s in &gls {
            let post = log10_normalized_posteriors(g_s, &log10_af, &gts);
            log10_p_no_variant += match span_del_index {
                None => post[0],
                Some(sd) => {
                    let idx = [genotype_index(0, 0), genotype_index(0, sd), genotype_index(sd, sd)];
                    let v: Vec<f64> = idx.iter().map(|&i| post[i]).collect();
                    (0.0f64).min(log10_sum_log10(&v))
                }
            };
            if !biallelic {
                let mut buffers: Vec<Vec<f64>> = vec![Vec::new(); n];
                for (g, &(i, j)) in gts.iter().enumerate() {
                    for (a, b) in buffers.iter_mut().enumerate() {
                        if a != i && a != j {
                            b.push(post[g]);
                        }
                    }
                }
                for a in 0..n {
                    let x = log10_sum_log10(&buffers[a]).min(0.0);
                    log10_p_of_zero[a] += x;
                }
            }
        }
        if biallelic {
            log10_p_of_zero[1] = log10_p_no_variant;
        }
        let mle_alt_counts = allele_counts[1..].iter().map(|&x| java_round(x) as i32).collect();
        AfResult { mle_alt_counts, log10_p_no_variant, log10_p_ref_by_allele: log10_p_of_zero }
    }
}

/// `AlleleSubsettingUtils.subsettedPLIndices` for ploidy 2: for each new genotype, the old index.
/// `keep[k]` = old index of the k-th kept allele.
pub fn subsetted_pl_indices(keep: &[usize]) -> Vec<usize> {
    genotype_alleles(keep.len()).into_iter().map(|(i, j)| genotype_index(keep[i], keep[j])).collect()
}

/// Output of `AlleleSubsettingUtils.subsetAlleles` + `makeGenotypeCall(USE_PLS_TO_ASSIGN)` for one sample.
#[derive(Clone, Debug)]
pub struct CalledGenotype {
    /// genotype allele indices (into the kept allele list); None = no-call
    pub alleles: Option<(usize, usize)>,
    pub gq: Option<i32>,
    pub pl: Vec<i32>,
}

/// Subsets the PLs of the original genotype (over `n_orig` alleles) to `keep` and assigns the
/// genotype from the subsetted PLs (`GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN`).
pub fn subset_and_call(orig_pls: &[i32], keep: &[usize]) -> CalledGenotype {
    let orig_gls: Vec<f64> = orig_pls.iter().map(|&p| p as f64 / -10.0).collect();
    let idx = subsetted_pl_indices(keep);
    let mut new_lks: Vec<f64> = idx.iter().map(|&i| orig_gls[i]).collect();
    // scaleLogSpaceArrayForNumericalStability
    let max = new_lks[mathutils::max_element_index(&new_lks)];
    for x in new_lks.iter_mut() {
        *x -= max;
    }
    let pl = gls_to_pls(&new_lks);
    let best = mathutils::max_element_index(&new_lks);
    let gq = gq_log10_from_likelihoods(best, &new_lks);
    let gts = genotype_alleles(keep.len());
    let alleles = if best == 0 && gq > -0.1 { None } else { Some(gts[best]) };
    let gq_int = if keep.len() > 1 { if gq == 1.0 { None } else { Some(java_round(gq * -10.0) as i32) } } else { None };
    CalledGenotype { alleles, gq: gq_int, pl }
}
