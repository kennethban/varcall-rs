//! `SomaticClusteringModel` (FilterMutectCalls): allele-fraction clusters of somatic variants (a background
//! beta-binomial, a high-AF beta-binomial and up to 5 binomial clusters), learned by EM, plus the somatic
//! priors per indel length.

use super::fmc_math::{beta_binomial_log_probability, binomial_probability, jexp, jlog, log1mexp};
use super::somatic::{log_dirichlet_normalization, log_sum_exp, normalize_from_log_to_linear, stream_sum};
use hc_assembly::m2math::digamma;
use hc_genotype::javamath::format_fixed;

const MAX_INDEL_SIZE_IN_PRIOR_MAP: i32 = 10;
const NUM_INITIALIZATION_QUANTILES: usize = 50;
const MIN_QUANTILE_INDEX_FOR_MAKING_CLUSTER: usize = 5; // (int) (0.1 * 50)
const INITIAL_HIGH_AF_WEIGHT: f64 = 0.01;
const MAX_FRACTION_OF_BACKGROUND_TO_SPLIT_OFF: f64 = 0.9;
const REGULARIZING_PSEUDOCOUNT: f64 = 1.0;
const NUM_ITERATIONS: usize = 5;
const MAX_BINOMIAL_CLUSTERS: usize = 5;
const OBVIOUS_ARTIFACT_PROBABILITY_THRESHOLD: f64 = 0.9;

#[derive(Clone, Debug)]
pub struct Datum {
    pub tumor_log_odds: f64,
    pub artifact_prob: f64,
    pub non_sequencing_error_prob: f64,
    pub alt_count: i32,
    pub total_count: i32,
    pub indel_length: i32,
}

impl Datum {
    pub fn new(tumor_log_odds: f64, artifact_prob: f64, non_somatic_prob: f64, alt_count: i32, total_count: i32, indel_length: i32) -> Datum {
        Datum {
            tumor_log_odds,
            artifact_prob,
            non_sequencing_error_prob: 1.0 - (1.0 - artifact_prob) * (1.0 - non_somatic_prob),
            alt_count,
            total_count,
            indel_length,
        }
    }
}

#[derive(Clone, Copy, Debug)]
enum Cluster {
    /// BetaBinomialCluster(alpha, beta)
    BetaBinomial(f64, f64),
    /// BinomialCluster: a narrow beta (alpha, beta)
    Binomial(f64, f64),
}

fn g(a: f64, b: f64) -> f64 {
    log_dirichlet_normalization(&[a, b])
}

/// `BetaBinomialCluster.correctedLogLikelihood(datum, shape)`
fn corrected_log_likelihood(d: &Datum, alpha: f64, beta: f64) -> f64 {
    let alt = d.alt_count as f64;
    let rf = (d.total_count - d.alt_count) as f64;
    d.tumor_log_odds + (g(alpha, beta) - g(alpha + alt, beta + rf) - g(1.0, 1.0) + g(1.0 + alt, 1.0 + rf))
}

/// `BinomialCluster.getFuzzyBinomial(mean)`
fn fuzzy_binomial(unbounded_mean: f64) -> (f64, f64) {
    const STD_DEV_OVER_MEAN: f64 = 0.01;
    let mean = unbounded_mean.min(1.0 - STD_DEV_OVER_MEAN);
    let alpha_plus_beta = ((1.0 - mean) / (mean * (STD_DEV_OVER_MEAN * STD_DEV_OVER_MEAN))) - 1.0;
    let alpha = mean * alpha_plus_beta;
    (alpha, alpha_plus_beta - alpha)
}

impl Cluster {
    fn shape(&self) -> (f64, f64) {
        match *self {
            Cluster::BetaBinomial(a, b) | Cluster::Binomial(a, b) => (a, b),
        }
    }
    fn corrected_log_likelihood(&self, d: &Datum) -> f64 {
        let (a, b) = self.shape();
        corrected_log_likelihood(d, a, b)
    }
    fn log_likelihood(&self, total: i32, alt: i32) -> f64 {
        let (a, b) = self.shape();
        beta_binomial_log_probability(a, b, total, alt)
    }
    fn learn(&mut self, data: &[Datum], resp: &[f64]) {
        match self {
            Cluster::BetaBinomial(alpha, beta) => {
                const RATE: f64 = 0.01;
                let (mut a, mut b) = (*alpha, *beta);
                for _epoch in 0..10 {
                    for (n, d) in data.iter().enumerate() {
                        let alt = d.alt_count as f64;
                        let rf = (d.total_count - d.alt_count) as f64;
                        let dg_total = digamma(d.total_count as f64 + a + b);
                        let dg_ab = digamma(a + b);
                        let ga = digamma(a + alt) - dg_total - digamma(a) + dg_ab;
                        let gb = digamma(b + rf) - dg_total - digamma(b) + dg_ab;
                        a = (a + RATE * ga * resp[n]).max(1.0);
                        b = (b + RATE * gb * resp[n]).max(0.5);
                    }
                }
                *alpha = a;
                *beta = b;
            }
            Cluster::Binomial(alpha, beta) => {
                let alt = stream_sum(data.iter().enumerate().map(|(n, d)| d.alt_count as f64 * resp[n])) + 0.0001;
                let total = stream_sum(data.iter().enumerate().map(|(n, d)| d.total_count as f64 * resp[n])) + 0.0001;
                let (a, b) = fuzzy_binomial(alt / total);
                *alpha = a;
                *beta = b;
            }
        }
    }
    fn describe(&self) -> String {
        match *self {
            Cluster::BetaBinomial(a, b) => format!("alpha = {}, beta = {}", format_fixed(a, 2), format_fixed(b, 2)),
            Cluster::Binomial(a, b) => format!("mean = {}", format_fixed(a / (a + b), 3)),
        }
    }
}

/// Iteration order of a `java.util.HashMap<Integer, _>` filled with `keys` (distinct) in order.
fn java_int_hashmap_order(keys: &[i32]) -> Vec<usize> {
    let n = keys.len();
    let mut cap: u64 = 16;
    while n as u64 > cap * 3 / 4 {
        cap *= 2;
    }
    let mask = (cap - 1) as u32;
    let mut order: Vec<(u32, usize)> = keys
        .iter()
        .enumerate()
        .map(|(i, &k)| {
            let h = k as u32;
            ((h ^ (h >> 16)) & mask, i)
        })
        .collect();
    order.sort_by_key(|&(b, _)| b);
    order.into_iter().map(|(_, i)| i).collect()
}

pub struct SomaticClusteringModel {
    initialized: bool,
    /// `Int2DoubleArrayMap` (insertion order)
    log_variant_priors: Vec<(i32, f64)>,
    log_variant_vs_artifact_prior: f64,
    callable_sites: Option<f64>,
    log_cluster_weights: Vec<f64>,
    clusters: Vec<Cluster>,
    data: Vec<Datum>,
    obvious_artifact_count: i32,
}

impl SomaticClusteringModel {
    pub fn new(log_snv_prior: f64, log_indel_prior: f64, log_artifact_prior: f64, callable: Option<f64>) -> Self {
        let mut pri: Vec<(i32, f64)> = (-MAX_INDEL_SIZE_IN_PRIOR_MAP..=MAX_INDEL_SIZE_IN_PRIOR_MAP).map(|n| (n, log_indel_prior)).collect();
        pri.iter_mut().find(|e| e.0 == 0).unwrap().1 = log_snv_prior;
        let callable_sites = match callable {
            Some(c) if c < 1.0 => {
                eprintln!("WARN: No callable sites found in Mutect stats.  Running without the full somatic clustering model.");
                None
            }
            c => c,
        };
        SomaticClusteringModel {
            initialized: false,
            log_variant_priors: pri,
            log_variant_vs_artifact_prior: log_artifact_prior,
            callable_sites,
            log_cluster_weights: vec![crate::jmath::log1p(INITIAL_HIGH_AF_WEIGHT), jlog(INITIAL_HIGH_AF_WEIGHT)],
            clusters: vec![Cluster::BetaBinomial(1.0, 1.0), Cluster::BetaBinomial(10.0, 1.0)],
            data: Vec::new(),
            obvious_artifact_count: 0,
        }
    }

    /// `record(tumorADs, tumorLogOdds, artifactProbabilities, nonSomaticProbabilities, vc)`
    /// `record`: symbolic alts (`symbolic[i]`) add nothing, and their AD is left out of the total. `artifact` and
    /// `non_somatic` have no entries for symbolic alts and are indexed by alt index as GATK does (its symbolic alt,
    /// <NON_REF>, comes last).
    pub fn record(&mut self, tumor_ads: &[i32], tumor_log_odds: &[f64], artifact: &[f64], non_somatic: &[f64], indel_lengths: &[i32], symbolic: &[bool]) {
        let total: i32 = tumor_ads.iter().enumerate().filter(|&(a, _)| a == 0 || !symbolic[a - 1]).map(|(_, &x)| x).sum();
        for i in 0..tumor_log_odds.len() {
            if symbolic[i] {
                continue;
            }
            if artifact[i] > OBVIOUS_ARTIFACT_PROBABILITY_THRESHOLD {
                self.obvious_artifact_count += 1;
                continue;
            } else if non_somatic[i] > OBVIOUS_ARTIFACT_PROBABILITY_THRESHOLD {
                continue;
            }
            self.data.push(Datum::new(tumor_log_odds[i], artifact[i], non_somatic[i], tumor_ads[i + 1], total, indel_lengths[i]));
        }
    }

    /// `getLogPriorOfSomaticVariant(indelLength)` (adds the minimum prior for lengths not in the map)
    pub fn log_prior_of_somatic_variant(&mut self, indel_length: i32) -> f64 {
        let v = match self.log_variant_priors.iter().find(|e| e.0 == indel_length) {
            Some(e) => e.1,
            None => {
                let m = self.log_variant_priors.iter().map(|e| e.1).fold(f64::INFINITY, f64::min);
                self.log_variant_priors.push((indel_length, m));
                m
            }
        };
        v + if indel_length == 0 { -jlog(3.0) } else { 0.0 }
    }

    pub fn log_prior_of_variant_versus_artifact(&self) -> f64 {
        self.log_variant_vs_artifact_prior
    }

    /// `Mutect2FilteringEngine.posteriorProbabilityOfError(logOddsOfRealVersusError, logPriorOfReal)`
    pub fn posterior_probability_of_error(log_odds: f64, log_prior: f64) -> f64 {
        let mut v = [log_odds + log_prior, log1mexp(log_prior)];
        normalize_from_log_to_linear(&mut v);
        v[1]
    }

    pub fn probability_of_sequencing_error(&mut self, d: &Datum) -> f64 {
        let lls: Vec<f64> = (0..self.clusters.len()).map(|c| self.log_cluster_weights[c] + self.clusters[c].corrected_log_likelihood(d)).collect();
        let variant_ll = log_sum_exp(&lls);
        let prior = self.log_prior_of_somatic_variant(d.indel_length);
        Self::posterior_probability_of_error(variant_ll, prior)
    }

    fn probability_of_somatic_variant(&mut self, d: &Datum) -> f64 {
        let seq = self.probability_of_sequencing_error(d);
        (1.0 - d.artifact_prob) * (1.0 - d.non_sequencing_error_prob) * (1.0 - seq)
    }

    fn cluster_log_likelihoods(&self, total: i32, alt: i32) -> Vec<f64> {
        (0..self.clusters.len()).map(|c| self.log_cluster_weights[c] + self.clusters[c].log_likelihood(total, alt)).collect()
    }

    pub fn log_likelihood_given_somatic(&self, total: i32, alt: i32) -> f64 {
        log_sum_exp(&self.cluster_log_likelihoods(total, alt))
    }

    fn background_prob_given_somatic(&self, total: i32, alt: i32) -> f64 {
        let mut v = self.cluster_log_likelihoods(total, alt);
        normalize_from_log_to_linear(&mut v);
        v[0]
    }

    fn data_somatic_probs(&mut self) -> Vec<f64> {
        let data = std::mem::take(&mut self.data);
        let v = data.iter().map(|d| self.probability_of_somatic_variant(d)).collect();
        self.data = data;
        v
    }

    fn initialize_clusters(&mut self) {
        assert!(!self.initialized);
        let somatic_probs = self.data_somatic_probs();
        let mut previous_bic = f64::NEG_INFINITY;
        for _cluster in 0..MAX_BINOMIAL_CLUSTERS {
            let old_weights = self.log_cluster_weights.clone();
            let background_given_somatic: Vec<f64> = self.data.iter().map(|d| self.background_prob_given_somatic(d.total_count, d.alt_count)).collect();
            let background: Vec<f64> = somatic_probs.iter().zip(&background_given_somatic).map(|(a, b)| a * b).collect();
            let quantiles = self.allele_fraction_quantiles();
            let resp = self.quantile_background_responsibilities(&quantiles, &background);
            let peaks = peaks_and_masses(&quantiles, &resp);
            if peaks.is_empty() {
                break;
            }
            // sorted by mass, descending (stable), first
            let mut biggest = peaks[0];
            for &p in &peaks[1..] {
                if p.1 > biggest.1 {
                    biggest = p;
                }
            }
            if std::env::var_os("VARIANT_RS_FMC_DUMP").is_some() {
                println!("peak {} {} npeaks {}", super::fmc_math::java_double_to_string(biggest.0), super::fmc_math::java_double_to_string(biggest.1), peaks.len());
            }
            if biggest.0 < quantiles[MIN_QUANTILE_INDEX_FOR_MAKING_CLUSTER.min(quantiles.len() - 1)] {
                break;
            }
            let total_mass = stream_sum(peaks.iter().map(|p| p.1));
            let fraction = MAX_FRACTION_OF_BACKGROUND_TO_SPLIT_OFF.min(biggest.1 / total_mass);
            let new_cluster_log_weight = jlog(fraction) + self.log_cluster_weights[0];
            let new_background_weight = crate::jmath::log1p(fraction) + self.log_cluster_weights[0];
            let (a, b) = fuzzy_binomial(biggest.0);
            self.clusters.push(Cluster::Binomial(a, b));
            self.log_cluster_weights.push(new_cluster_log_weight);
            self.log_cluster_weights[0] = new_background_weight;
            for _ in 0..NUM_ITERATIONS {
                self.em_iteration(false);
                if std::env::var_os("VARIANT_RS_FMC_DUMP").is_some() {
                    println!("init em{}", self.debug_state());
                }
            }
            let lls: Vec<f64> = self.data.iter().map(|d| self.log_likelihood_given_somatic(d.total_count, d.alt_count)).collect();
            let mut weighted = 0.0;
            for (p, l) in somatic_probs.iter().zip(&lls) {
                weighted += p * l;
            }
            let mut effective = 0.0;
            for &p in &somatic_probs {
                effective += p;
            }
            let num_parameters = 2.0 * self.clusters.len() as f64;
            let bic = weighted - num_parameters * jlog(effective);
            if std::env::var_os("VARIANT_RS_FMC_DUMP").is_some() {
                println!("bic {}", super::fmc_math::java_double_to_string(bic));
            }
            if bic < previous_bic {
                self.clusters.pop();
                self.log_cluster_weights = old_weights;
                break;
            }
            previous_bic = bic;
        }
        self.initialized = true;
    }

    fn allele_fraction_quantiles(&mut self) -> Vec<f64> {
        let data = std::mem::take(&mut self.data);
        let mut afs: Vec<(f64, f64)> = data.iter().map(|d| (d.alt_count as f64 / d.total_count as f64, self.probability_of_somatic_variant(d))).collect();
        self.data = data;
        afs.sort_by(|a, b| crate::jmath::java_cmp(a.0, b.0));
        let total = stream_sum(afs.iter().map(|p| p.1));
        let mut cumulative = 0.0;
        let step = total / NUM_INITIALIZATION_QUANTILES as f64;
        let mut quantile_prob = step;
        let mut out: Vec<f64> = Vec::new();
        for p in &afs {
            cumulative += p.1;
            if cumulative > quantile_prob {
                out.push(p.0);
                while cumulative > quantile_prob {
                    quantile_prob += step;
                }
            }
        }
        // distinct (Double.equals)
        let mut d: Vec<f64> = Vec::new();
        for x in out {
            if !d.iter().any(|y| y.to_bits() == x.to_bits()) {
                d.push(x);
            }
        }
        d
    }

    fn quantile_background_responsibilities(&self, quantiles: &[f64], background: &[f64]) -> Vec<f64> {
        let mut total = vec![0.0; quantiles.len()];
        for (n, d) in self.data.iter().enumerate() {
            let bp = background[n];
            for (q, &f) in quantiles.iter().enumerate() {
                let x = binomial_probability(d.total_count, d.alt_count, f);
                total[q] += x * bp * (d.total_count + 1) as f64;
            }
        }
        total
    }

    pub fn debug_data(&self) -> String {
        let d = super::fmc_math::java_double_to_string;
        let mut s = String::new();
        for x in &self.data {
            s.push_str(&format!("datum {} {} {} {} {} {}\n", d(x.tumor_log_odds), d(x.artifact_prob), d(x.non_sequencing_error_prob), x.alt_count, x.total_count, x.indel_length));
        }
        s.push_str(&format!("obvious {}", self.obvious_artifact_count));
        s
    }

    pub fn learn_and_clear(&mut self) {
        if !self.initialized {
            if std::env::var_os("VARIANT_RS_FMC_DUMP").is_some() {
                let data = self.data.clone();
                for x in &data {
                    let sp = self.probability_of_somatic_variant(x);
                    let seq = self.probability_of_sequencing_error(x);
                    println!("sprob {} seq {}", super::fmc_math::java_double_to_string(sp), super::fmc_math::java_double_to_string(seq));
                }
            }
            if std::env::var_os("VARIANT_RS_FMC_DUMP").is_some() {
                let q = self.allele_fraction_quantiles();
                println!("quantiles {}", q.iter().map(|&x| super::fmc_math::java_double_to_string(x)).collect::<Vec<_>>().join(" "));
            }
            self.initialize_clusters();
            if std::env::var_os("VARIANT_RS_FMC_DUMP").is_some() {
                println!("after init{}", self.debug_state());
            }
        }
        for _ in 0..NUM_ITERATIONS {
            self.em_iteration(true);
        }
        self.data.clear();
        self.obvious_artifact_count = 0;
    }

    fn em_iteration(&mut self, update_priors: bool) {
        // variantCountsByIndelLength: HashMap<Integer, MutableDouble>, keys -10..10 then putIfAbsent
        let mut count_keys: Vec<i32> = (-MAX_INDEL_SIZE_IN_PRIOR_MAP..=MAX_INDEL_SIZE_IN_PRIOR_MAP).collect();
        let mut counts: Vec<f64> = vec![0.0; count_keys.len()];
        let nc = self.clusters.len();
        let mut responsibilities: Vec<Vec<f64>> = Vec::with_capacity(self.data.len());
        let mut totals = vec![0.0; nc];
        let data = std::mem::take(&mut self.data);
        for d in &data {
            let somatic = self.probability_of_somatic_variant(d);
            let k = match count_keys.iter().position(|&k| k == d.indel_length) {
                Some(k) => k,
                None => {
                    count_keys.push(d.indel_length);
                    counts.push(0.0);
                    counts.len() - 1
                }
            };
            counts[k] += somatic;
            let mut lls = self.cluster_log_likelihoods(d.total_count, d.alt_count);
            normalize_from_log_to_linear(&mut lls);
            let r: Vec<f64> = lls.iter().map(|x| somatic * x).collect();
            for c in 0..nc {
                totals[c] += r[c];
            }
            responsibilities.push(r);
        }
        self.data = data;
        for t in totals.iter_mut() {
            *t += REGULARIZING_PSEUDOCOUNT;
        }
        let mut s = 0.0;
        for &t in &totals {
            s += t;
        }
        self.log_cluster_weights = totals.iter().map(|t| jlog(t / s)).collect();
        let technical_artifact_count = self.obvious_artifact_count as f64 + stream_sum(self.data.iter().map(|d| d.artifact_prob));
        let order = java_int_hashmap_order(&count_keys);
        let variant_count = stream_sum(order.iter().map(|&i| counts[i]));
        if update_priors {
            self.log_variant_vs_artifact_prior =
                jlog((variant_count + REGULARIZING_PSEUDOCOUNT) / (variant_count + technical_artifact_count + REGULARIZING_PSEUDOCOUNT * 2.0));
            if let Some(callable) = self.callable_sites {
                for n in -MAX_INDEL_SIZE_IN_PRIOR_MAP..=MAX_INDEL_SIZE_IN_PRIOR_MAP {
                    let c = count_keys.iter().position(|&k| k == n).map_or(0.0, |i| counts[i]);
                    let empirical = c / callable;
                    let v = jlog(empirical.max(if n == 0 { 1.0e-8 } else { 1.0e-9 }));
                    self.log_variant_priors.iter_mut().find(|e| e.0 == n).unwrap().1 = v;
                }
            }
        }
        for c in 0..nc {
            let col: Vec<f64> = responsibilities.iter().map(|r| r[c]).collect();
            let mut cl = self.clusters[c];
            cl.learn(&self.data, &col);
            self.clusters[c] = cl;
        }
    }

    /// Learned state in the format of `somatic/javacheck/FmcCheck.java` (debugging).
    pub fn debug_state(&self) -> String {
        let d = super::fmc_math::java_double_to_string;
        let mut s = format!(" artifactPrior={} weights=", d(self.log_variant_vs_artifact_prior));
        for &w in &self.log_cluster_weights {
            s.push_str(&format!("{},", d(w)));
        }
        s.push_str(" clusters=");
        for c in &self.clusters {
            let (a, b) = c.shape();
            s.push_str(&format!("{}/{},", d(a), d(b)));
        }
        s.push_str(" priors=");
        for n in -MAX_INDEL_SIZE_IN_PRIOR_MAP..=MAX_INDEL_SIZE_IN_PRIOR_MAP {
            s.push_str(&format!("{n}:{},", d(self.log_variant_priors.iter().find(|e| e.0 == n).unwrap().1)));
        }
        s
    }

    /// `clusteringMetadata()`: (key, value) pairs for the filtering stats
    pub fn metadata(&self) -> Vec<(String, String)> {
        let mut out = Vec::new();
        for n in -MAX_INDEL_SIZE_IN_PRIOR_MAP..=MAX_INDEL_SIZE_IN_PRIOR_MAP {
            let lp = self.log_variant_priors.iter().find(|e| e.0 == n).unwrap().1;
            let ty = if n == 0 { "SNV".to_string() } else { format!("{} of length {}", if n < 0 { "deletion" } else { "insertion" }, n.abs()) };
            out.push((format!("Ln prior of {ty}"), super::fmc_math::java_double_to_string(lp)));
        }
        let w = |c: usize| format_fixed(jexp(self.log_cluster_weights[c]), 4);
        out.push(("Background beta-binomial cluster".into(), format!("weight = {}, {}", w(0), self.clusters[0].describe())));
        out.push(("High-AF beta-binomial cluster".into(), format!("weight = {}, {}", w(1), self.clusters[1].describe())));
        let mut rest: Vec<usize> = (2..self.clusters.len()).collect();
        rest.sort_by(|&a, &b| crate::jmath::java_cmp(-self.log_cluster_weights[a], -self.log_cluster_weights[b]));
        for c in rest {
            out.push(("Binomial cluster".into(), format!("weight = {}, {}", w(c), self.clusters[c].describe())));
        }
        out
    }
}

/// `calculatePeaksAndMasses`
fn peaks_and_masses(q: &[f64], resp: &[f64]) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    let mut mass = 0.0;
    let mut peak = 0.0;
    let mut peak_resp = 0.0;
    let n = q.len();
    for i in 0..n {
        let left_r = if i == 0 { 0.0 } else { resp[i - 1] };
        let r = resp[i];
        let right_r = if i == n - 1 { 0.0 } else { resp[i + 1] };
        let left_af = if i == 0 { 0.0 } else { q[i - 1] };
        let af = q[i];
        mass += (af - left_af) * (left_r + r) / 2.0;
        if r > peak_resp {
            peak = af;
            peak_resp = r;
        }
        let lc = crate::jmath::java_cmp(r, left_r);
        let rc = crate::jmath::java_cmp(r, right_r);
        use std::cmp::Ordering::*;
        let local_min = (lc == Less && rc != Greater) || (lc != Greater && rc == Less);
        if (local_min && i > 0) || i == n - 1 {
            out.push((peak, mass));
            mass = 0.0;
            peak = af;
            peak_resp = r;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn int_hashmap_order() {
        // keys -10..10: Integer hash spread, capacity 32
        let keys: Vec<i32> = (-10..=10).collect();
        let o = java_int_hashmap_order(&keys);
        // bucket 0 holds -1 (0xffffffff ^ 0xffff = 0xffff0000) and then 0, in insertion order
        assert_eq!(keys[o[0]], -1);
        assert_eq!(keys[o[1]], 0);
        assert_eq!(o.len(), 21);
    }
}
