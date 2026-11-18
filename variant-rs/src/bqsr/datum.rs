//! `RecalDatum` arithmetic, `QualQuantizer`, and the Java numerics they depend on.

use crate::jmath;
use hc_genotype::javamath::{binomial_log_probability, java_round};

pub const MULTIPLIER: f64 = 100000.0;
const SMOOTHING: i64 = 1;
const MAX_RECALIBRATED_Q_SCORE: i32 = 93;
const MAX_REASONABLE_Q_SCORE: usize = 60;
const MAX_GATK_USABLE_Q_SCORE: usize = 40;

/// `QualityUtils.qualToErrorProb(double)` = `Math.pow(10, q / -10)` (HotSpot's pow: glibc's differs
/// in the last bit for some non-integer q, e.g. the combined reported qualities of GatherBQSRReports).
#[inline]
pub fn qual_to_error_prob(q: f64) -> f64 {
    jmath::pow(10.0, q / -10.0)
}

/// `Math.ulp(double)` for finite x.
#[inline]
fn java_ulp(x: f64) -> f64 {
    let a = x.abs();
    if a == f64::MAX {
        return 2f64.powi(971);
    }
    f64::from_bits(a.to_bits() + 1) - a
}

/// `MathUtils.roundToNDecimalPlaces`.
pub fn round_to_n_decimal_places(x: f64, n: i32) -> f64 {
    let mult = jmath::pow(10.0, n as f64); // Math.pow(10, n)
    java_round((x + java_ulp(x)) * mult) as f64 / mult
}

/// `NormalDistribution(0, 0.5).logDensity(i)` for i in 0..=40 (`RecalDatum.logPriorCache`).
fn log_prior_cache() -> &'static [f64; MAX_GATK_USABLE_Q_SCORE + 1] {
    use std::sync::OnceLock;
    static C: OnceLock<[f64; MAX_GATK_USABLE_Q_SCORE + 1]> = OnceLock::new();
    C.get_or_init(|| {
        let k = hc_assembly::fastmath::log(0.5) + 0.5 * hc_assembly::fastmath::log(2.0 * std::f64::consts::PI); // FastMath.log
        let mut t = [0.0; MAX_GATK_USABLE_Q_SCORE + 1];
        for (i, v) in t.iter_mut().enumerate() {
            let x1 = i as f64 / 0.5;
            *v = -0.5 * x1 * x1 - k;
        }
        t
    })
}

fn error_probs() -> &'static [f64; MAX_REASONABLE_Q_SCORE + 1] {
    use std::sync::OnceLock;
    static C: OnceLock<[f64; MAX_REASONABLE_Q_SCORE + 1]> = OnceLock::new();
    C.get_or_init(|| {
        let mut t = [0.0; MAX_REASONABLE_Q_SCORE + 1];
        for (i, v) in t.iter_mut().enumerate() {
            *v = qual_to_error_prob(i as f64);
        }
        t
    })
}

/// `RecalDatum.bayesianEstimateOfEmpiricalQuality`.
pub fn bayesian_estimate(n_obs: i64, n_err: i64, prior: f64) -> i32 {
    let lp = log_prior_cache();
    let ep = error_probs();
    let mut best = 0usize;
    let mut best_v = f64::NEG_INFINITY;
    for q in 0..=MAX_REASONABLE_Q_SCORE {
        // getLogPrior: (int) truncation of the difference
        let diff = (((q as f64) - prior) as i32).unsigned_abs().min(MAX_GATK_USABLE_Q_SCORE as u32) as usize;
        let log_prior = lp[diff];
        // getLogBinomialLikelihood
        let ll = if n_obs == 0 {
            0.0
        } else {
            let (mut n, mut k) = (n_obs, n_err);
            const MAX_OBS: i64 = i32::MAX as i64 - 1;
            if n > MAX_OBS {
                let fraction = MAX_OBS as f64 / n as f64;
                k = java_round(k as f64 * fraction);
                n = MAX_OBS;
            }
            let l = binomial_log_probability(n as i32, k as i32, ep[q]);
            if l.is_infinite() || l.is_nan() {
                -f64::MAX
            } else {
                l
            }
        };
        let v = log_prior + ll;
        if q == 0 || v > best_v {
            best = q;
            best_v = v;
        }
    }
    best as i32
}

/// One table cell as GATK keeps it: observations, mismatches (stored x 100000), reported quality.
#[derive(Clone, Copy, Debug)]
pub struct Datum {
    pub obs: i64,
    pub mism: f64,
    pub reported: f64,
}

impl Datum {
    pub fn from_counts(obs: u64, err: u64, qual: u8) -> Datum {
        // RecalDatum(1, isError, qual) followed by increment(1, isError) per further observation:
        // integer-valued doubles, exact in any order
        Datum { obs: obs as i64, mism: err as f64 * MULTIPLIER, reported: qual as f64 }
    }
    #[inline]
    pub fn num_mismatches(&self) -> f64 {
        self.mism / MULTIPLIER
    }
    fn expected_errors(&self) -> f64 {
        self.obs as f64 * qual_to_error_prob(self.reported)
    }
    /// `RecalDatum.combine`.
    pub fn combine(&mut self, o: &Datum) {
        let expected = self.expected_errors() + o.expected_errors();
        self.obs += o.obs;
        self.mism += o.num_mismatches() * MULTIPLIER;
        self.reported = -10.0 * jmath::log10(expected / self.obs as f64);
    }
    /// `BaseRecalibrationEngine.roundTableValues` for this datum.
    pub fn round(&mut self) {
        let m = round_to_n_decimal_places(self.num_mismatches(), 2);
        self.mism = m * MULTIPLIER;
        self.reported = round_to_n_decimal_places(self.reported, 4);
    }
    /// `RecalDatum.calcEmpiricalQuality(prior)`.
    pub fn empirical_quality(&self, prior: f64) -> i32 {
        let mismatches = (self.num_mismatches() + 0.5) as i64 + SMOOTHING;
        let observations = self.obs + SMOOTHING + SMOOTHING;
        bayesian_estimate(observations, mismatches, prior).min(MAX_RECALIBRATED_Q_SCORE)
    }
}

/// `QualityUtils.errorProbToQual(p)` (max 93).
fn error_prob_to_qual(p: f64) -> u8 {
    let d = java_round(-10.0 * jmath::log10(p)) as f64;
    let q = d as i32; // (int) of a double: saturating
    q.min(93).max(1) as u8
}

#[derive(Clone)]
struct QualInterval {
    q_start: usize,
    q_end: usize,
    n_obs: i64,
    n_err: i64,
    fixed_qual: i32,
    subs: Vec<usize>, // indices into the arena
}

/// `QualQuantizer(nObservationsPerQual, nLevels, minInterestingQual).getOriginalToQuantizedMap()`.
pub fn quantize(hist: &[i64], n_levels: usize, min_interesting_qual: usize) -> Vec<u8> {
    let mut arena: Vec<QualInterval> = Vec::new();
    let qtab: Vec<f64> = (0..256).map(|i| qual_to_error_prob(i as f64)).collect();
    // ordered set by q_start (TreeSet); intervals never overlap
    let mut set: Vec<usize> = Vec::new();
    for (q, &n) in hist.iter().enumerate() {
        let err = qtab[q];
        let n_err = (n as f64 * err).floor() as i32 as i64; // (int) Math.floor(...)
        arena.push(QualInterval { q_start: q, q_end: q, n_obs: n, n_err, fixed_qual: q as i32, subs: vec![] });
        set.push(arena.len() - 1);
    }
    fn error_rate(a: &[QualInterval], qtab: &[f64], i: usize) -> f64 {
        let x = &a[i];
        if x.fixed_qual != -1 {
            qtab[(x.fixed_qual as usize) & 0xff]
        } else if x.n_obs == 0 {
            0.0
        } else {
            (x.n_err + 1) as f64 / (1.0 * (x.n_obs + 1) as f64)
        }
    }
    fn calc_penalty(a: &[QualInterval], qtab: &[f64], i: usize, global: f64, min_q: usize) -> f64 {
        if global == 0.0 {
            return 0.0;
        }
        let x = &a[i];
        if x.subs.is_empty() {
            if x.q_end <= min_q {
                0.0
            } else {
                (jmath::log10(error_rate(a, qtab, i)) - jmath::log10(global)).abs() * x.n_obs as f64
            }
        } else {
            let mut sum = 0.0;
            for &s in &x.subs {
                sum += calc_penalty(a, qtab, s, global, min_q);
            }
            sum
        }
    }
    while set.len() > n_levels {
        let mut best: Option<(usize, f64)> = None; // (position of left in set, penalty)
        let mut merged_best: Option<QualInterval> = None;
        for p in 0..set.len() - 1 {
            let (l, r) = (set[p], set[p + 1]);
            assert_eq!(arena[l].q_end + 1, arena[r].q_start);
            let m = QualInterval {
                q_start: arena[l].q_start,
                q_end: arena[r].q_end,
                n_obs: arena[l].n_obs + arena[r].n_obs,
                n_err: arena[l].n_err + arena[r].n_err,
                fixed_qual: -1,
                subs: vec![l, r],
            };
            arena.push(m);
            let mi = arena.len() - 1;
            let global = error_rate(&arena, &qtab, mi);
            let pen = calc_penalty(&arena, &qtab, mi, global, min_interesting_qual);
            let m = arena.pop().unwrap();
            if best.map_or(true, |(_, bp)| pen < bp) {
                best = Some((p, pen));
                merged_best = Some(m);
            }
        }
        let (p, _) = best.unwrap();
        arena.push(merged_best.unwrap());
        let mi = arena.len() - 1;
        set.remove(p + 1);
        set[p] = mi;
    }
    let mut map = vec![0u8; hist.len()];
    for &i in &set {
        let q = if arena[i].fixed_qual != -1 { arena[i].fixed_qual as u8 } else { error_prob_to_qual(error_rate(&arena, &qtab, i)) };
        for x in arena[i].q_start..=arena[i].q_end {
            map[x] = q;
        }
    }
    map
}
