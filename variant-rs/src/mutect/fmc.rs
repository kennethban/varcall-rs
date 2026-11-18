//! `variant-rs filtermutectcalls`: port of GATK 4.7 FilterMutectCalls (default arguments, plus contamination
//! tables / estimate and tumor segmentation tables; no orientation priors or panel of normals). Four passes:
//! two learning passes (strand-artifact parameters, filtered haplotypes, the somatic clustering model and
//! the threshold), one for the threshold only, and the output pass (`Mutect2FilteringEngine`).

use super::cluster::{Datum, SomaticClusteringModel};
use super::fmc_math::*;
use super::somatic::{log_sum_exp, stream_sum};
use std::collections::HashMap;
use std::io::BufRead;

const USAGE: &str = "usage: variant-rs filtermutectcalls [-R ref.fa] -V calls.vcf -O filtered.vcf [--stats calls.vcf.stats] [--filtering-stats out.tsv] \
[--contamination-table t.table ...] [--contamination-estimate 0.0] [--tumor-segmentation segments.table ...] \
[--mitochondria-mode] [--min-allele-fraction 0.0] [--max-alt-allele-count 1] [--f-score-beta 1.0]";
const FMC_HEADER_LINES: &str = include_str!("../../tables/fmc_header_lines.txt");

const LN_10: f64 = std::f64::consts::LN_10;
const EPSILON: f64 = 1.0e-10;
const MIN_REPORTABLE_ERROR_PROBABILITY: f64 = 0.1;
const SITE: &str = "SITE";

// M2FiltersArgumentCollection defaults
const INITIAL_THRESHOLD: f64 = 0.1;
const MAX_EVENTS_IN_REGION: i32 = 3;
const MAX_EVENTS_IN_HAPLOTYPE: i32 = 2;
const UNIQUE_ALT_READ_COUNT: i32 = 0;
const MIN_MEDIAN_MAPPING_QUALITY: i32 = 30;
const LONG_INDEL_LENGTH: i32 = 5;
const MIN_MEDIAN_BASE_QUALITY: i32 = 20;
const MAX_MEDIAN_FRAGMENT_LENGTH_DIFFERENCE: i32 = 10000;
const MIN_MEDIAN_READ_POSITION: i32 = 1;
const MAX_N_RATIO: f64 = f64::INFINITY;
const MIN_READS_ON_EACH_STRAND: i32 = 0;
const CONTAMINATION_ESTIMATE: f64 = 0.0;
const NORMAL_P_VALUE_THRESHOLD: f64 = 0.001;
const MIN_SLIPPAGE_LENGTH: i32 = 8;
const SLIPPAGE_RATE: f64 = 0.1;
const MAX_INTRA_HAPLOTYPE_DISTANCE: i64 = 100;

// ------------------------------------------------------------------------------------------------
// VCF records

struct Genotype {
    ad: Vec<i32>,
    gt: Option<String>,
    f1r2: Option<Vec<i32>>,
    f2r1: Option<Vec<i32>>,
    af: Option<String>,
    pgt: Option<String>,
    pid: Option<String>,
}

struct Rec {
    /// CHROM .. QUAL
    head: Vec<String>,
    start: i64,
    alleles: Vec<Vec<u8>>,
    info: Vec<(String, Option<String>)>,
    /// FORMAT and sample columns, unchanged
    tail: String,
    genotypes: Vec<Genotype>,
}

impl Rec {
    fn n_alt(&self) -> usize {
        self.alleles.len() - 1
    }
    /// Per alt: htsjdk `Allele.isSymbolic` (<NON_REF>, <DEL>, breakends)
    fn symbolic(&self) -> Vec<bool> {
        self.alleles[1..]
            .iter()
            .map(|a| {
                a.len() > 1
                    && (a[0] == b'<' || a[a.len() - 1] == b'>' || a.contains(&b'[') || a.contains(&b']') || a[0] == b'.' || a[a.len() - 1] == b'.')
            })
            .collect()
    }
    fn get(&self, key: &str) -> Option<&str> {
        self.info.iter().find(|(k, _)| k == key).map(|(_, v)| v.as_deref().unwrap_or(""))
    }
    fn has(&self, key: &str) -> bool {
        self.info.iter().any(|(k, _)| k == key)
    }
    /// `getAttributeAsIntList(key, default)`
    fn int_list(&self, key: &str, default: i32) -> Vec<i32> {
        match self.get(key) {
            None => Vec::new(),
            Some(v) => v.split(',').map(|s| if s == "." || s.is_empty() { default } else { s.trim().parse().unwrap() }).collect(),
        }
    }
    fn int(&self, key: &str, default: i32) -> i32 {
        self.get(key).map_or(default, |v| if v == "." { default } else { v.trim().parse().unwrap() })
    }
    /// `VariantContextGetters.getAttributeAsDoubleArray(vc, key)` (`.` -> default)
    fn double_list(&self, key: &str, default: f64) -> Option<Vec<f64>> {
        self.get(key).map(|v| v.split(',').map(|s| if s == "." { default } else { s.trim().parse().unwrap() }).collect())
    }
    /// `SomaticClusteringModel.indelLength`: htsjdk `Allele.length()` is 0 for a symbolic allele
    fn indel_length(&self, alt: usize) -> i32 {
        let l = if self.symbolic()[alt] { 0 } else { self.alleles[alt + 1].len() as i32 };
        l - self.alleles[0].len() as i32
    }
    /// `Mutect2FilteringEngine.getTumorLogOdds`: TLOD in natural log
    fn tumor_log_odds(&self) -> Option<Vec<f64>> {
        self.double_list("TLOD", f64::NAN).map(|v| v.into_iter().map(|x| x * LN_10).collect())
    }
}

fn parse_record(line: &str, n_samples: usize) -> Rec {
    let cols: Vec<&str> = line.splitn(9, '\t').collect();
    assert!(cols.len() >= 8, "malformed VCF line: {line}");
    let head: Vec<String> = cols[..6].iter().map(|s| s.to_string()).collect();
    let start: i64 = cols[1].parse().unwrap();
    let mut alleles = vec![cols[3].as_bytes().to_vec()];
    alleles.extend(cols[4].split(',').map(|a| a.as_bytes().to_vec()));
    let info: Vec<(String, Option<String>)> = if cols[7] == "." {
        Vec::new()
    } else {
        cols[7]
            .split(';')
            .map(|kv| match kv.split_once('=') {
                Some((k, v)) => (k.to_string(), Some(v.to_string())),
                None => (kv.to_string(), None),
            })
            .collect()
    };
    let tail = cols.get(8).map_or(String::new(), |s| s.to_string());
    let mut genotypes = Vec::with_capacity(n_samples);
    if !tail.is_empty() {
        let mut it = tail.split('\t');
        let keys: Vec<&str> = it.next().unwrap().split(':').collect();
        for s in it {
            let vals: Vec<&str> = s.split(':').collect();
            let get = |k: &str| keys.iter().position(|x| *x == k).and_then(|i| vals.get(i)).filter(|v| **v != ".").map(|v| v.to_string());
            let ad = get("AD").map_or(Vec::new(), |v| v.split(',').map(|x| x.parse().unwrap_or(0)).collect());
            let ints = |k: &str| get(k).map(|v| v.split(',').map(|x| x.parse().unwrap_or(0)).collect::<Vec<i32>>());
            genotypes.push(Genotype { ad, gt: get("GT"), f1r2: ints("F1R2"), f2r1: ints("F2R1"), af: get("AF"), pgt: get("PGT"), pid: get("PID") });
        }
    }
    Rec { head, start, alleles, info, tail, genotypes }
}

// ------------------------------------------------------------------------------------------------
// filters

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum ErrorType {
    Artifact,
    NonSomatic,
    Sequencing,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum F {
    TumorEvidence,
    BaseQuality,
    MappingQuality,
    DuplicatedAltRead,
    StrandArtifact,
    Contamination,
    StrictStrandBias,
    ReadPosition,
    MinAlleleFraction,
    NormalArtifact,
    NRatio,
    PanelOfNormals,
    ReadOrientation,
    ClusteredEvents,
    Multiallelic,
    FragmentLength,
    PolymeraseSlippage,
    FilteredHaplotype,
    Germline,
}

/// `buildFiltersList` order (ReadOrientation is inserted after PanelOfNormals with --ob-priors)
const FILTERS: [F; 18] = [
    F::TumorEvidence,
    F::BaseQuality,
    F::MappingQuality,
    F::DuplicatedAltRead,
    F::StrandArtifact,
    F::Contamination,
    F::StrictStrandBias,
    F::ReadPosition,
    F::MinAlleleFraction,
    F::NormalArtifact,
    F::NRatio,
    F::PanelOfNormals,
    F::ClusteredEvents,
    F::Multiallelic,
    F::FragmentLength,
    F::PolymeraseSlippage,
    F::FilteredHaplotype,
    F::Germline,
];

impl F {
    fn name(self) -> &'static str {
        match self {
            F::TumorEvidence => "weak_evidence",
            F::BaseQuality => "base_qual",
            F::MappingQuality => "map_qual",
            F::DuplicatedAltRead => "duplicate",
            F::StrandArtifact => "strand_bias",
            F::Contamination => "contamination",
            F::StrictStrandBias => "strict_strand",
            F::ReadPosition => "position",
            F::MinAlleleFraction => "low_allele_frac",
            F::NormalArtifact => "normal_artifact",
            F::NRatio => "n_ratio",
            F::PanelOfNormals => "panel_of_normals",
            F::ReadOrientation => "orientation",
            F::ClusteredEvents => "clustered_events",
            F::Multiallelic => "multiallelic",
            F::FragmentLength => "fragment",
            F::PolymeraseSlippage => "slippage",
            F::FilteredHaplotype => "haplotype",
            F::Germline => "germline",
        }
    }
    fn error_type(self) -> ErrorType {
        match self {
            F::TumorEvidence => ErrorType::Sequencing,
            F::Contamination | F::Germline => ErrorType::NonSomatic,
            _ => ErrorType::Artifact,
        }
    }
    /// `Mutect2VariantFilter` (site-level) vs `Mutect2AlleleFilter`
    fn is_variant_filter(self) -> bool {
        matches!(
            self,
            F::NormalArtifact | F::NRatio | F::PanelOfNormals | F::ReadOrientation | F::ClusteredEvents | F::Multiallelic | F::FragmentLength | F::PolymeraseSlippage | F::FilteredHaplotype | F::Germline
        )
    }
    fn required_info(self) -> &'static [&'static str] {
        match self {
            F::TumorEvidence | F::Multiallelic => &["TLOD"],
            F::BaseQuality => &["MBQ"],
            F::MappingQuality => &["MMQ"],
            F::DuplicatedAltRead => &["UNIQ_ALT_READ_COUNT"],
            F::StrandArtifact | F::StrictStrandBias => &["AS_SB_TABLE"],
            F::Contamination => &["POPAF"],
            F::ReadPosition => &["MPOS"],
            F::NormalArtifact => &["NALOD", "TLOD"],
            F::NRatio => &["NCount"],
            F::ClusteredEvents => &["ECNT", "ECNTH"],
            F::FragmentLength => &["MFRL"],
            F::PolymeraseSlippage => &["RPA", "RU"],
            F::Germline => &["TLOD", "POPAF"],
            F::MinAlleleFraction | F::PanelOfNormals | F::FilteredHaplotype | F::ReadOrientation => &[],
        }
    }
    /// `phredScaledPosteriorAnnotationName` (only used for variant filters)
    fn qual_annotation(self) -> Option<&'static str> {
        match self {
            F::PolymeraseSlippage => Some("STRQ"),
            F::Germline => Some("GERMQ"),
            F::TumorEvidence => Some("SEQQ"),
            F::StrandArtifact => Some("STRANDQ"),
            F::Contamination => Some("CONTQ"),
            F::ReadOrientation => Some("ROQ"),
            _ => None,
        }
    }
}

#[inline]
fn round_finite(p: f64) -> f64 {
    // Math.max(Math.min(p, 1.0), 0.0) (NaN stays NaN)
    if p.is_nan() { p } else { p.min(1.0).max(0.0) }
}

/// `Double.compare`-based maximum of a stream (NaN is the largest)
fn java_max(it: impl Iterator<Item = f64>) -> Option<f64> {
    // BinaryOperator.maxBy(Double::compare): (a, b) -> compare(a, b) >= 0 ? a : b
    let cmp = crate::jmath::java_cmp;
    it.fold(None, |m, x| match m {
        None => Some(x),
        Some(m) => Some(if cmp(m, x) != std::cmp::Ordering::Less { m } else { x }),
    })
}

#[derive(Clone, Copy)]
struct EStep {
    fwd_resp: f64,
    rev_resp: f64,
    fwd_count: i32,
    rev_count: i32,
    fwd_alt: i32,
    rev_alt: i32,
}

impl EStep {
    fn artifact_probability(&self) -> f64 {
        self.fwd_resp + self.rev_resp
    }
}

struct ErrorProbs {
    by_filter: Vec<(F, Vec<f64>)>,
    artifact: Vec<f64>,
    non_somatic: Vec<f64>,
    combined: Vec<f64>,
}

#[derive(Default)]
struct OutputStats {
    pass: i64,
    tps: f64,
    fps: f64,
    fns: f64,
    filter_fps: Vec<f64>,
    filter_fns: Vec<f64>,
}

struct Engine {
    model: SomaticClusteringModel,
    threshold: f64,
    threshold_probs: Vec<f64>,
    is_normal: Vec<bool>,
    filters: Vec<F>,
    // ReadOrientationFilter: per sample artifact priors, and the reference
    ob_priors: Vec<Option<HashMap<String, super::orientation::Prior>>>,
    fasta: Option<crate::reference::Fasta>,
    // ContaminationFilter: per sample (table value or --contamination-estimate)
    contamination: Vec<f64>,
    // GermlineFilter: per sample tumor segments (contig, start, end, minor allele fraction)
    segments: Vec<Vec<(String, i64, i64, f64)>>,
    // StrandArtifactFilter
    alpha_strand: f64,
    beta_strand: f64,
    strand_artifact_prior: f64,
    esteps: Vec<EStep>,
    // FilteredHaplotypeFilter
    accumulating_phased: HashMap<String, Vec<(i64, f64)>>,
    phased: HashMap<String, Vec<(i64, f64)>>,
    stats: OutputStats,
    // --min-allele-fraction, --max-alt-allele-count, --f-score-beta
    min_af: f64,
    max_alt: usize,
    f_score_beta: f64,
}

impl Engine {
    fn sum_ads(&self, r: &Rec, tumor: bool, normal: bool) -> Vec<i32> {
        let mut ads = vec![0i32; r.alleles.len()];
        for (s, g) in r.genotypes.iter().enumerate() {
            if (tumor && !self.is_normal[s]) || (normal && self.is_normal[s]) {
                for n in 0..ads.len() {
                    ads[n] += g.ad[n];
                }
            }
        }
        ads
    }

    fn posterior_of_error(&mut self, r: &Rec, log_odds: f64, alt: usize) -> f64 {
        let prior = self.model.log_prior_of_somatic_variant(r.indel_length(alt));
        SomaticClusteringModel::posterior_probability_of_error(log_odds, prior)
    }

    /// One filter's `errorProbabilities(vc)`; None = not applied (empty list).
    fn filter_probs(&mut self, f: F, r: &Rec) -> Option<Vec<f64>> {
        let nalt = r.n_alt();
        let present = f.required_info().iter().all(|k| r.has(k));
        if f.is_variant_filter() {
            let p = if present { self.variant_probability(f, r) } else { 0.0 };
            return Some(vec![round_finite(p); nalt]);
        }
        if !present {
            return None;
        }
        let v = match f {
            F::TumorEvidence => {
                let tlods = r.tumor_log_odds().unwrap();
                let ads = self.sum_ads(r, true, false);
                let total: i32 = ads.iter().sum();
                (0..tlods.len())
                    .map(|i| {
                        let d = Datum::new(tlods[i], 0.0, 0.0, ads[i + 1], total, r.indel_length(i));
                        self.model.probability_of_sequencing_error(&d)
                    })
                    .collect()
            }
            F::BaseQuality => r.int_list("MBQ", 0).iter().skip(1).map(|&q| if (q as f64) < MIN_MEDIAN_BASE_QUALITY as f64 { 1.0 } else { 0.0 }).collect(),
            F::MappingQuality => {
                let lens = indel_lengths_if_indel(r);
                let mut mq = r.int_list("MMQ", 0);
                let ref_q = mq.remove(0);
                for i in 0..mq.len() {
                    if let Some(l) = &lens {
                        if l[i] >= LONG_INDEL_LENGTH {
                            mq[i] = ref_q;
                        }
                    }
                }
                mq.iter().map(|&q| if q < MIN_MEDIAN_MAPPING_QUALITY { 1.0 } else { 0.0 }).collect()
            }
            F::DuplicatedAltRead => r.int_list("UNIQ_ALT_READ_COUNT", 1).iter().map(|&c| if c <= UNIQUE_ALT_READ_COUNT { 1.0 } else { 0.0 }).collect(),
            F::StrandArtifact => {
                let es = self.strand_artifact_probabilities(r);
                if es.is_empty() {
                    return None;
                }
                es.iter().map(|e| e.fwd_resp + e.rev_resp).collect()
            }
            F::Contamination => self.contamination(r),
            F::StrictStrandBias => {
                let sbs = sbs_for_alleles(r);
                if MIN_READS_ON_EACH_STRAND == 0 || sbs.len() <= 1 {
                    return None;
                }
                sbs[1..].iter().map(|s| if s.iter().any(|&x| x == 0) { 1.0 } else { 0.0 }).collect()
            }
            F::ReadPosition => r.int_list("MPOS", 0).iter().map(|&p| if p > -1 && (p as f64) < MIN_MEDIAN_READ_POSITION as f64 { 1.0 } else { 0.0 }).collect(),
            F::MinAlleleFraction => {
                // max tumor AF per alt (default 1.0) < minAf
                let mut maxes: Vec<Option<f64>> = vec![None; nalt];
                for (s, g) in r.genotypes.iter().enumerate() {
                    if self.is_normal[s] {
                        continue;
                    }
                    if let Some(af) = &g.af {
                        for (i, x) in af.split(',').enumerate().take(nalt) {
                            let x: f64 = if x == "." { 1.0 } else { x.parse().unwrap() };
                            maxes[i] = Some(maxes[i].map_or(x, |m: f64| if x > m { x } else { m }));
                        }
                    }
                }
                maxes.iter().map(|m| if m.unwrap_or(1.0) < self.min_af { 1.0 } else { 0.0 }).collect()
            }
            _ => unreachable!(),
        };
        let v: Vec<f64> = v.into_iter().map(round_finite).collect();
        if v.is_empty() { None } else { Some(v) }
    }

    fn variant_probability(&mut self, f: F, r: &Rec) -> f64 {
        match f {
            F::NormalArtifact => self.normal_artifact(r),
            F::NRatio => {
                let ads = self.sum_ads(r, true, true);
                let alt = ads.iter().sum::<i32>() - ads[0];
                if alt == 0 {
                    0.0
                } else if r.int("NCount", 0) as f64 / alt as f64 >= MAX_N_RATIO {
                    1.0
                } else {
                    0.0
                }
            }
            F::ReadOrientation => self.read_orientation(r),
            F::PanelOfNormals => {
                if r.has("PON") {
                    1.0
                } else {
                    0.0
                }
            }
            F::ClusteredEvents => {
                let h = r.int_list("ECNTH", 0).into_iter().max().unwrap();
                if h > MAX_EVENTS_IN_HAPLOTYPE || r.int("ECNT", 0) > MAX_EVENTS_IN_REGION {
                    1.0
                } else {
                    0.0
                }
            }
            F::Multiallelic => {
                let n = r.tumor_log_odds().unwrap().iter().filter(|&&l| l > 5.0).count();
                if n > self.max_alt {
                    1.0
                } else {
                    0.0
                }
            }
            F::FragmentLength => {
                let fl = r.int_list("MFRL", 0);
                if ((fl[1] - fl[0]).abs() as f64) > MAX_MEDIAN_FRAGMENT_LENGTH_DIFFERENCE as f64 {
                    1.0
                } else {
                    0.0
                }
            }
            F::PolymeraseSlippage => self.slippage(r),
            F::FilteredHaplotype => self.filtered_haplotype(r),
            F::Germline => self.germline(r),
            _ => unreachable!(),
        }
    }

    fn normal_artifact(&mut self, r: &Rec) -> f64 {
        let tlods = r.tumor_log_odds().unwrap();
        let idx = max_element_index(&tlods);
        let t_ads = self.sum_ads(r, true, false);
        let t_depth: i32 = t_ads.iter().sum();
        let t_alt = t_ads[idx + 1];
        let n_ads = self.sum_ads(r, false, true);
        let n_depth: i32 = n_ads.iter().sum();
        let n_alt = n_ads[idx + 1];
        let t_af = t_alt as f64 / t_depth as f64;
        let n_af = if n_depth == 0 { 0.0 } else { n_alt as f64 / n_depth as f64 };
        if n_af < 0.1 * t_af {
            return 0.0;
        }
        let neg: Vec<f64> = r.double_list("NALOD", f64::NAN).unwrap().into_iter().map(|x| -(x * LN_10)).collect();
        let prob = SomaticClusteringModel::posterior_probability_of_error(neg[idx], self.model.log_prior_of_variant_versus_artifact());
        let ref_bq = r.int_list("MBQ", 30)[0];
        let p_err = crate::jmath::pow(10.0, ref_bq as f64 / -10.0);
        let p_value = 1.0 - binomial_cdf(n_depth, p_err, n_alt - 1);
        if p_value < NORMAL_P_VALUE_THRESHOLD { 1.0 } else { prob }
    }

    fn slippage(&mut self, r: &Rec) -> f64 {
        let rpa = r.int_list("RPA", 0);
        if rpa.len() < 2 {
            return 0.0;
        }
        let ru = r.get("RU").unwrap_or("");
        let ref_str_bases = ru.len() as i32 * rpa[0];
        let slips = rpa[0] - rpa[1];
        if ref_str_bases >= MIN_SLIPPAGE_LENGTH && slips.abs() == 1 {
            let ads = self.sum_ads(r, true, false);
            if ads.len() < 2 {
                return 0.0;
            }
            let depth: i32 = ads.iter().sum();
            let alt = depth - ads[0];
            let log_somatic = self.model.log_likelihood_given_somatic(depth, alt);
            let lik = regularized_beta(SLIPPAGE_RATE, ads[1] as f64 + 1.0, ads[0] as f64 + 1.0);
            let log_odds = log_somatic - jlog(lik);
            self.posterior_of_error(r, log_odds, 0)
        } else {
            0.0
        }
    }

    fn filtered_haplotype(&self, r: &Rec) -> f64 {
        // tumor genotype with the largest AF (first on ties)
        let mut best: Option<(usize, f64)> = None;
        for (s, g) in r.genotypes.iter().enumerate() {
            if self.is_normal[s] {
                continue;
            }
            let afs: Vec<f64> = g.af.as_ref().map_or(vec![0.0], |v| v.split(',').map(|x| if x == "." { 0.0 } else { x.parse().unwrap() }).collect());
            let m = afs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            if best.map_or(true, |(_, b)| m > b) {
                best = Some((s, m));
            }
        }
        let g = &r.genotypes[best.unwrap().0];
        let ps = match phasing_string(g) {
            Some(p) => p,
            None => return 0.0,
        };
        match self.phased.get(&ps) {
            None => 0.0,
            Some(v) => v.iter().filter(|(pos, _)| (pos - r.start).abs() <= MAX_INTRA_HAPLOTYPE_DISTANCE).map(|p| p.1).fold(None, |m: Option<f64>, x| Some(m.map_or(x, |m| m.max(x)))).unwrap_or(0.0),
        }
    }

    fn weighted_average_of_tumor_afs(&self, r: &Rec) -> Vec<f64> {
        let mut total = 0.0;
        let mut afs = vec![0.0; r.n_alt()];
        for (s, g) in r.genotypes.iter().enumerate() {
            if self.is_normal[s] {
                continue;
            }
            let w: f64 = g.ad.iter().map(|&x| x as i64).sum::<i64>() as f64;
            total += w;
            let sample: Vec<f64> = g.af.as_ref().map_or(vec![0.0], |v| v.split(',').map(|x| if x == "." { 0.0 } else { x.parse().unwrap() }).collect());
            for (i, x) in sample.iter().enumerate().take(afs.len()) {
                afs[i] += w * x;
            }
        }
        let inv = 1.0 / total;
        afs.iter().map(|x| x * inv).collect()
    }

    fn germline(&mut self, r: &Rec) -> f64 {
        let tlods = r.tumor_log_odds().unwrap();
        let idx = max_element_index(&tlods);
        let normal_lods: Option<Vec<f64>> = r.double_list("NLOD", f64::NAN).map(|v| v.into_iter().map(|x| x * LN_10).collect());
        let neg_log10_af = r.double_list("POPAF", f64::NAN).unwrap();
        let pop_af = crate::jmath::pow(10.0, -neg_log10_af[idx]);
        if pop_af < EPSILON {
            return 0.0;
        } else if pop_af > 1.0 - EPSILON {
            return 1.0;
        }
        let counts = self.sum_ads(r, true, false);
        let total: i32 = counts.iter().sum();
        if total == 0 {
            return 0.0;
        }
        let alt = counts[idx + 1];
        let alt_af = self.weighted_average_of_tumor_afs(r)[idx];
        // computeMinorAlleleFraction: the overlapping tumor segment's MAF, else 0.5
        let (contig, start) = (r.head[0].as_str(), r.head[1].parse::<i64>().unwrap());
        let end = start + r.alleles[0].len() as i64 - 1;
        let mut wsum = 0.0;
        for (s, g) in r.genotypes.iter().enumerate() {
            if !self.is_normal[s] {
                let maf = self.segments[s].iter().find(|x| x.0 == contig && x.1 <= end && start <= x.2).map_or(0.5, |x| x.3);
                wsum += maf * g.ad.iter().map(|&x| x as i64).sum::<i64>() as f64;
            }
        }
        let maf = wsum / counts.iter().map(|&x| x as i64).sum::<i64>() as f64;
        let log_germline = jlog(0.5) + log_sum_exp(&[binomial_log_probability(total, alt, maf), binomial_log_probability(total, alt, 1.0 - maf)]);
        let log_somatic = self.model.log_likelihood_given_somatic(total, alt);
        let het_vs_somatic = log_germline - log_somatic;
        let hom_alt_vs_somatic = if alt_af < 0.9 { f64::NEG_INFINITY } else { 0.0 };
        let normal_lod = normal_lods.map_or(0.0, |v| v[idx]);
        let log_prior_somatic = self.model.log_prior_of_somatic_variant(r.indel_length(idx));
        germline_probability(-normal_lod, het_vs_somatic, hom_alt_vs_somatic, pop_af, log_prior_somatic)
    }

    /// `ReadOrientationFilter.calculateErrorProbability`
    fn read_orientation(&mut self, r: &Rec) -> f64 {
        let len = r.alleles[0].len();
        let snp_or_mnp = r.alleles.iter().all(|a| a.len() == len && !a.contains(&b'*') && !a.starts_with(b"<")) && r.alleles[1..].iter().any(|a| *a != r.alleles[0]);
        if !snp_or_mnp {
            return 0.0;
        }
        let mut v: Vec<(i32, f64)> = Vec::new();
        for (s, g) in r.genotypes.iter().enumerate() {
            if self.is_normal[s] {
                continue;
            }
            let p = self.orientation_artifact_probability(r, s, g);
            let alt = g.ad.iter().sum::<i32>() - g.ad.first().copied().unwrap_or(0);
            v.push((alt, p));
        }
        weighted_median_posterior(v)
    }

    fn orientation_artifact_probability(&self, r: &Rec, s: usize, g: &Genotype) -> f64 {
        let hom_ref = g.gt.as_ref().is_some_and(|gt| gt.split(['/', '|']).all(|a| a == "0"));
        if hom_ref {
            return 0.0;
        }
        let Some(priors) = &self.ob_priors[s] else { return 0.0 };
        let tlods = r.double_list("TLOD", -1.0).unwrap_or_default();
        let idx = max_element_index(&tlods);
        let alt = &r.alleles[idx + 1];
        let fasta = self.fasta.as_ref().expect("the read orientation filter needs the reference (-R)");
        let contig = r.head[0].as_str();
        let mut best: Option<f64> = None;
        for (n, &b) in alt.iter().enumerate() {
            let pos = r.start + n as i64;
            let p = (|| -> f64 {
                let (seq, clen) = fasta.fetch(contig, pos - 1, pos + 1);
                if pos - 1 < 1 || pos + 1 > clen || seq.len() != 3 || seq.contains(&b'N') {
                    return 0.0;
                }
                let base = |c: u8| match c {
                    b'A' => Some(0usize),
                    b'C' => Some(1),
                    b'G' => Some(2),
                    b'T' => Some(3),
                    _ => None,
                };
                let (Some(ref_allele), Some(alt_base)) = (base(seq[1]), base(b.to_ascii_uppercase())) else { return 0.0 };
                let (Some(f1r2), Some(f2r1)) = (&g.f1r2, &g.f2r1) else { return 0.0 };
                let ref_count = f1r2[0] + f2r1[0];
                let alt_f1r2 = f1r2[idx + 1];
                let alt_count = alt_f1r2 + f2r1[idx + 1];
                let Some(prior) = priors.get(std::str::from_utf8(&seq).unwrap()) else { return 0.0 };
                let post = super::orientation::compute_responsibilities(ref_allele, alt_base, alt_count, alt_f1r2, ref_count + alt_count, &prior.pi, true);
                super::cm_linalg::jmax(post[alt_base], post[4 + alt_base])
            })();
            best = Some(best.map_or(p, |m| super::cm_linalg::jmax(m, p)));
        }
        best.unwrap_or(0.0)
    }

    fn contamination(&mut self, r: &Rec) -> Vec<f64> {
        let nalt = r.n_alt();
        let mut per_allele: Vec<Vec<(i32, f64)>> = vec![Vec::new(); nalt];
        for (s, g) in r.genotypes.iter().enumerate() {
            if self.is_normal[s] {
                continue;
            }
            let c = self.contamination[s].min(1.0 - EPSILON).max(0.0);
            let total: i32 = g.ad.iter().sum();
            let alt_ads = &g.ad[1..];
            let neg = r.double_list("POPAF", f64::INFINITY).unwrap_or(vec![f64::INFINITY, f64::INFINITY]);
            let afs: Vec<f64> = neg.iter().map(|&x| crate::jmath::pow(10.0, -x)).collect();
            let log_somatic: Vec<f64> = alt_ads.iter().map(|&a| self.model.log_likelihood_given_somatic(total, a)).collect();
            let mut log_odds = vec![0.0; afs.len()];
            for i in 0..afs.len() {
                let single = 2.0 * afs[i] * (1.0 - afs[i]) * binomial_probability(total, alt_ads[i], c / 2.0) + (afs[i] * afs[i]) * binomial_probability(total, alt_ads[i], c);
                let many = binomial_probability(total, alt_ads[i], c * afs[i]);
                let log_contaminant = jlog(single.max(many));
                log_odds[i] = log_somatic[i] - log_contaminant;
            }
            for i in 0..afs.len() {
                let p = self.posterior_of_error(r, log_odds[i], i);
                per_allele[i].push((alt_ads[i], p));
            }
        }
        per_allele.into_iter().map(|v| if v.is_empty() { f64::NAN } else { weighted_median_posterior(v) }).collect()
    }

    /// `StrandArtifactFilter.calculateArtifactProbabilities`
    fn strand_artifact_probabilities(&self, r: &Rec) -> Vec<EStep> {
        let mut sbs = sbs_for_alleles(r);
        if sbs.len() <= 1 {
            return Vec::new();
        }
        // symbolic alts' counts removed (removeDataForSymbolicAlleles: by allele index, the reference included)
        let sym = r.symbolic();
        if sym.iter().any(|&x| x) {
            sbs = sbs.into_iter().enumerate().filter(|&(a, _)| a == 0 || !sym.get(a - 1).copied().unwrap_or(false)).map(|(_, x)| x).collect();
        }
        let ref_len = r.alleles[0].len() as i32;
        let total_fwd: i32 = sbs.iter().map(|s| s[0]).sum();
        let total_rev: i32 = sbs.iter().map(|s| s[1]).sum();
        sbs[1..]
            .iter()
            .enumerate()
            .map(|(i, sb)| {
                let alt_len = if sym.get(i).copied().unwrap_or(false) { 0 } else { r.alleles[i + 1].len() as i32 };
                let indel = (ref_len - alt_len).abs();
                if sb.iter().sum::<i32>() == 0 || indel > 4 {
                    EStep { fwd_resp: 0.0, rev_resp: 0.0, fwd_count: total_fwd, rev_count: total_rev, fwd_alt: sb[0], rev_alt: sb[1] }
                } else {
                    self.strand_artifact_probability(total_fwd, total_rev, sb[0], sb[1], indel)
                }
            })
            .collect()
    }

    fn strand_artifact_probability(&self, fc: i32, rc: i32, fa: i32, ra: i32, indel: i32) -> EStep {
        let art = |n: i32, k: i32| beta_binomial_log_probability(self.alpha_strand, self.beta_strand, n, k);
        let non_art = |n: i32, k: i32| {
            let beta_seq = if indel == 0 { 1000.0 } else if indel < 3 { 5000.0 } else { 50000.0 };
            beta_binomial_log_probability(1.0, beta_seq, n, k)
        };
        let fwd_ll = art(fc, fa) + non_art(rc, ra);
        let rev_ll = art(rc, ra) + non_art(fc, fa);
        let none_ll = binomial_coefficient_log(fc as i64, fa as i64) + binomial_coefficient_log(rc as i64, ra as i64)
            - binomial_coefficient_log((fc + rc) as i64, (fa + ra) as i64)
            + beta_binomial_log_probability(1.0, 1.0, fc + rc, fa + ra);
        let fwd_prior = jlog(self.strand_artifact_prior / 2.0);
        let rev_prior = jlog(self.strand_artifact_prior / 2.0);
        let none_prior = jlog(1.0 - self.strand_artifact_prior);
        let l10e = std::f64::consts::LOG10_E;
        let mut v = [(fwd_ll + fwd_prior) * l10e, (rev_ll + rev_prior) * l10e, (none_ll + none_prior) * l10e];
        normalize_log10_to_linear(&mut v);
        EStep { fwd_resp: v[0], rev_resp: v[1], fwd_count: fc, rev_count: rc, fwd_alt: fa, rev_alt: ra }
    }

    fn learn_strand(&mut self) {
        let potential: Vec<EStep> = self.esteps.iter().filter(|e| e.artifact_probability() > 0.1).cloned().collect();
        let total_artifacts = stream_sum(potential.iter().map(|e| e.artifact_probability()));
        let total_non_artifacts = stream_sum(self.esteps.iter().map(|e| 1.0 - e.artifact_probability()));
        self.strand_artifact_prior = (total_artifacts + 1.0) / (total_artifacts + 1.0 + total_non_artifacts + 1000.0);
        let alt_count = stream_sum(potential.iter().map(|e| e.fwd_resp * e.fwd_alt as f64 + e.rev_resp * e.rev_alt as f64));
        let depth = stream_sum(potential.iter().map(|e| e.fwd_resp * e.fwd_count as f64 + e.rev_resp * e.rev_count as f64));
        let mean = (alt_count + 1.0) / (depth + 1.0 + 20.0);
        let objective = |alpha: f64| {
            let beta = (1.0 / mean - 1.0) * alpha;
            stream_sum(potential.iter().map(|e| {
                e.fwd_resp * beta_binomial_log_probability(alpha, beta, e.fwd_count, e.fwd_alt)
                    + e.rev_resp * beta_binomial_log_probability(alpha, beta, e.rev_count, e.rev_alt)
            }))
        };
        self.alpha_strand = brent_max(objective, 0.01, 100.0, 1.0, 0.01, 0.01, 100);
        self.beta_strand = (1.0 / mean - 1.0) * self.alpha_strand;
        self.esteps.clear();
    }

    fn error_probabilities(&mut self, r: &Rec) -> ErrorProbs {
        let mut by_filter: Vec<(F, Vec<f64>)> = Vec::new();
        for f in self.filters.clone() {
            if let Some(v) = self.filter_probs(f, r) {
                by_filter.push((f, v));
            }
        }
        // ErrorProbabilities: data of symbolic alts (GVCF <NON_REF>) removed from every filter's list (after dropping
        // the filters that were not applied, so a list may become empty)
        let sym = r.symbolic();
        if sym.iter().any(|&s| s) {
            for (_, v) in by_filter.iter_mut() {
                *v = v.iter().enumerate().filter(|&(i, _)| !sym.get(i).copied().unwrap_or(false)).map(|(_, &p)| p).collect();
            }
        }
        let nalt = r.n_alt() - sym.iter().filter(|&&s| s).count();
        // per error type (LinkedHashMap order of first appearance), max over its filters per allele
        let mut types: Vec<(ErrorType, Vec<f64>)> = Vec::new();
        for (f, v) in &by_filter {
            let t = f.error_type();
            match types.iter_mut().find(|(tt, _)| *tt == t) {
                Some((_, m)) => {
                    for i in 0..nalt {
                        m[i] = java_max([m[i], v[i]].into_iter()).unwrap();
                    }
                }
                None => types.push((t, v.clone())),
            }
        }
        let get = |t: ErrorType| types.iter().find(|(tt, _)| *tt == t).map(|(_, v)| v.clone());
        // HashMap<ErrorType, _>.values(): enum identity hashes; declaration order used
        let mut combined = vec![1.0f64; nalt];
        for t in [ErrorType::Artifact, ErrorType::NonSomatic, ErrorType::Sequencing] {
            if let Some(v) = get(t) {
                for i in 0..nalt {
                    combined[i] *= 1.0 - v[i];
                }
            }
        }
        let combined = combined.into_iter().map(|p| round_finite(1.0 - p)).collect();
        ErrorProbs { artifact: get(ErrorType::Artifact).unwrap(), non_somatic: get(ErrorType::NonSomatic).unwrap(), by_filter, combined }
    }

    fn accumulate(&mut self, r: &Rec) {
        let ep = self.error_probabilities(r);
        // StrandArtifactFilter.accumulateDataForLearning
        if r.has("AS_SB_TABLE") {
            let es = self.strand_artifact_probabilities(r);
            self.esteps.extend(es);
        }
        // FilteredHaplotypeFilter.accumulateDataForLearning
        let germline_prob = java_max(ep.by_filter.iter().filter(|(f, _)| *f == F::Germline).flat_map(|(_, v)| v.iter().cloned())).unwrap_or(0.0);
        let ignore_normal_artifact = germline_prob > 0.25;
        let artifact_prob = java_max(
            ep.by_filter
                .iter()
                .filter(|(f, _)| f.error_type() != ErrorType::NonSomatic)
                .filter(|(f, _)| !(ignore_normal_artifact && *f == F::NormalArtifact))
                .filter(|(f, _)| *f != F::FilteredHaplotype)
                .flat_map(|(_, v)| v.iter().cloned()),
        )
        .unwrap_or(0.0);
        for (s, g) in r.genotypes.iter().enumerate() {
            if self.is_normal[s] {
                continue;
            }
            if let Some(ps) = phasing_string(g) {
                self.accumulating_phased.entry(ps).or_default().push((r.start, artifact_prob));
            }
        }
        let tumor_ads = self.sum_ads(r, true, false);
        let tlods = r.tumor_log_odds().unwrap();
        let lens: Vec<i32> = (0..r.n_alt()).map(|i| r.indel_length(i)).collect();
        self.model.record(&tumor_ads, &tlods, &ep.artifact, &ep.non_somatic, &lens, &r.symbolic());
        self.threshold_probs.extend(ep.combined);
    }

    fn learn_parameters(&mut self) {
        self.learn_strand();
        self.phased = std::mem::take(&mut self.accumulating_phased);
        self.model.learn_and_clear();
        self.learn_threshold();
    }

    fn learn_threshold(&mut self) {
        self.threshold = optimal_f_score_threshold(&mut self.threshold_probs, self.f_score_beta);
        self.threshold_probs.clear();
        self.stats = OutputStats::default();
    }

    /// `applyFiltersAndAccumulateOutputStats`: (FILTER column, extra INFO attributes)
    fn apply(&mut self, r: &Rec) -> (String, Vec<(String, Option<String>)>) {
        let ep = self.error_probabilities(r);
        self.record_stats(&ep, self.threshold - EPSILON);
        let error_threshold = (1.0 - EPSILON).min(EPSILON.max(self.threshold));
        let mut site: Vec<(String, f64)> = Vec::new(); // LinkedHashMap
        let put = |site: &mut Vec<(String, f64)>, k: &str, v: f64| match site.iter_mut().find(|(kk, _)| kk == k) {
            Some(e) => e.1 = v,
            None => site.push((k.to_string(), v)),
        };
        let allele_status: Vec<Vec<&str>> = ep
            .by_filter
            .iter()
            .filter(|(f, v)| !f.is_variant_filter() && !v.is_empty())
            .map(|(f, v)| v.iter().map(|&p| if p > error_threshold { f.name() } else { SITE }).collect())
            .collect();
        // non-symbolic alts (the filters' lists have no entries for symbolic ones)
        let nalt = if allele_status.is_empty() { 0 } else { allele_status[0].len() };
        let distinct_by_allele: Vec<Vec<&str>> = (0..nalt)
            .map(|i| {
                let mut d: Vec<&str> = Vec::new();
                for st in &allele_status {
                    if !d.contains(&st[i]) {
                        d.push(st[i]);
                    }
                }
                if d.len() > 1 && d.contains(&SITE) {
                    d.retain(|x| *x != SITE);
                }
                if d.is_empty() {
                    d.push(SITE);
                }
                d
            })
            .collect();
        let sym = r.symbolic();
        let mut merged = distinct_by_allele.iter().map(|d| d.join(","));
        let as_status = sym.iter().map(|&s| if s { SITE.to_string() } else { merged.next().expect("one filter status per non-symbolic alt") }).collect::<Vec<_>>().join("|");
        let mut extra: Vec<(String, Option<String>)> = vec![("AS_FilterStatus".into(), Some(as_status))];
        for st in &allele_status {
            if !st.is_empty() && st.iter().all(|x| *x == st[0]) && st[0] != SITE {
                put(&mut site, st[0], 1.0);
            }
        }
        for (f, v) in ep.by_filter.iter().filter(|(f, v)| f.is_variant_filter() && !v.is_empty()) {
            let p = v[0];
            if let Some(ann) = f.qual_annotation() {
                if f.required_info().iter().all(|k| r.has(k)) {
                    extra.push((ann.into(), Some(error_prob_to_qual(p).to_string())));
                }
            }
            if p > error_threshold {
                put(&mut site, f.name(), p);
            }
        }
        if site.is_empty() && !distinct_by_allele.iter().all(|d| d.is_empty()) {
            // GATK removes the symbolic alts' indices again (a no-op when they come last, as <NON_REF> does)
            let non_symbolic = distinct_by_allele.iter().enumerate().filter(|&(i, _)| !sym.get(i).copied().unwrap_or(false)).map(|(_, d)| d);
            if !non_symbolic.into_iter().any(|d| d.contains(&SITE)) {
                put(&mut site, "FAIL", 1.0);
            }
        }
        let max_p = site.iter().map(|e| e.1).fold(None, |m: Option<f64>, x| Some(m.map_or(x, |m| m.max(x)))).unwrap_or(1.0);
        let mut filters: Vec<String> = site.iter().filter(|e| e.1 >= max_p.min(MIN_REPORTABLE_ERROR_PROBABILITY)).map(|e| e.0.clone()).collect();
        filters.sort();
        filters.dedup();
        (if filters.is_empty() { "PASS".into() } else { filters.join(";") }, extra)
    }

    /// `FilteringOutputStats.recordCall`
    fn record_stats(&mut self, ep: &ErrorProbs, threshold: f64) {
        let st = &mut self.stats;
        st.filter_fps.resize(self.filters.len(), 0.0);
        st.filter_fns.resize(self.filters.len(), 0.0);
        let filtered: Vec<bool> = ep.combined.iter().map(|&p| p > threshold).collect();
        for &p in &ep.combined {
            if p > threshold {
                st.fns += 1.0 - p;
            } else {
                st.pass += 1;
                st.fps += p;
                st.tps += 1.0 - p;
            }
        }
        for i in 0..ep.combined.len() {
            for (f, v) in &ep.by_filter {
                let k = self.filters.iter().position(|x| x == f).unwrap();
                let p = v[i];
                if p > EPSILON && p > threshold - EPSILON {
                    st.filter_fns[k] += 1.0 - ep.combined[i];
                } else if !filtered[i] {
                    st.filter_fps[k] += p;
                }
            }
        }
    }

    fn filtering_stats(&self) -> String {
        let st = &self.stats;
        let mut out = String::new();
        let r3 = |x: f64| java_double_to_string(round_to_n_decimal_places(x, 3));
        let r2 = |x: f64| java_double_to_string(round_to_n_decimal_places(x, 2));
        for (k, v) in self.model.metadata() {
            out.push_str(&format!("#<METADATA>{k}={v}\n"));
        }
        out.push_str(&format!("#<METADATA>threshold={}\n", r3(self.threshold)));
        out.push_str(&format!("#<METADATA>fdr={}\n", r3(st.fps / st.pass as f64)));
        out.push_str(&format!("#<METADATA>sensitivity={}\n", r3(st.tps / (st.tps + st.fns))));
        out.push_str("filter\tFP\tFDR\tFN\tFNR\n");
        let total_true = st.tps + st.fns;
        for (k, f) in self.filters.iter().enumerate() {
            let (fp, fnc) = (st.filter_fps.get(k).copied().unwrap_or(0.0), st.filter_fns.get(k).copied().unwrap_or(0.0));
            if fp > 0.0 || fnc > 0.0 {
                out.push_str(&format!("{}\t{}\t{}\t{}\t{}\n", f.name(), r2(fp), r2(fp / st.pass as f64), r2(fnc), r2(fnc / total_true)));
            }
        }
        out
    }
}

/// `GermlineFilter.germlineProbability`
fn germline_probability(normal_log_odds: f64, het_vs_somatic: f64, hom_alt_vs_somatic: f64, pop_af: f64, log_prior_somatic: f64) -> f64 {
    let log_prior_not_somatic = log1mexp(log_prior_somatic);
    let prior_het = jlog(2.0 * pop_af * (1.0 - pop_af));
    let prior_hom_alt = jlog(pop_af * pop_af);
    let prior_not_germline = jlog((1.0 - pop_af) * (1.0 - pop_af));
    let p_het = prior_het + het_vs_somatic + normal_log_odds + log_prior_not_somatic;
    let p_hom = prior_hom_alt + hom_alt_vs_somatic + normal_log_odds + log_prior_not_somatic;
    let p_germline = log_sum_exp(&[p_het, p_hom]);
    let p_somatic = prior_not_germline + log_prior_somatic;
    let lse = log_sum_exp(&[p_germline, p_somatic]);
    jexp(p_germline - lse)
}

/// `Mutect2Filter.weightedMedianPosteriorProbability`
fn weighted_median_posterior(mut v: Vec<(i32, f64)>) -> f64 {
    let total: i32 = v.iter().map(|p| p.0).sum();
    v.sort_by(|a, b| crate::jmath::java_cmp(a.1, b.1));
    let mut cum = 0;
    for (d, p) in v {
        cum += d;
        if cum * 2 >= total {
            return p;
        }
    }
    0.0
}

/// `ThresholdCalculator.calculateThresholdBasedOnOptimalFScore`
fn optimal_f_score_threshold(posteriors: &mut [f64], beta: f64) -> f64 {
    posteriors.sort_by(|a, b| crate::jmath::java_cmp(*a, *b));
    let expected_tp = stream_sum(posteriors.iter().map(|p| 1.0 - p));
    let (mut tp, mut fp, mut fnn) = (0.0, 0.0, expected_tp);
    let mut optimal: i64 = -1;
    let mut optimal_f = 0.0;
    let n = posteriors.len();
    let b2 = beta * beta;
    for i in 0..n {
        tp += 1.0 - posteriors[i];
        fp += posteriors[i];
        fnn -= 1.0 - posteriors[i];
        let f = (1.0 + b2) * tp / ((1.0 + b2) * tp + b2 * fnn + fp);
        if f >= optimal_f {
            optimal = i as i64;
            optimal_f = f;
        }
    }
    if optimal == -1 {
        0.0
    } else if optimal as usize == n - 1 {
        1.0
    } else {
        posteriors[optimal as usize]
    }
}

fn max_element_index(v: &[f64]) -> usize {
    let mut m = 0;
    for i in 1..v.len() {
        if v[i] > v[m] {
            m = i;
        }
    }
    m
}

/// `StrandBiasUtils.getSBsForAlleles`
fn sbs_for_alleles(r: &Rec) -> Vec<Vec<i32>> {
    match r.get("AS_SB_TABLE") {
        None | Some("") => Vec::new(),
        Some(s) => s.split('|').map(|a| a.split(',').map(|x| x.trim().parse().unwrap()).collect()).collect(),
    }
}

/// `vc.getIndelLengths()`: alt - ref lengths when the type is INDEL or MIXED, else None
fn indel_lengths_if_indel(r: &Rec) -> Option<Vec<i32>> {
    #[derive(PartialEq, Clone, Copy)]
    enum T {
        Snp,
        Mnp,
        Indel,
    }
    let rl = r.alleles[0].len();
    let types: Vec<T> = r.alleles[1..].iter().map(|a| if a.len() == rl { if rl == 1 { T::Snp } else { T::Mnp } } else { T::Indel }).collect();
    let mixed = types.iter().any(|t| *t != types[0]);
    if mixed || types[0] == T::Indel {
        Some((0..r.n_alt()).map(|i| r.indel_length(i)).collect())
    } else {
        None
    }
}

fn phasing_string(g: &Genotype) -> Option<String> {
    match (&g.pgt, &g.pid) {
        (Some(a), Some(b)) => Some(format!("{a}{b}")),
        _ => None,
    }
}

// ------------------------------------------------------------------------------------------------
// driver

pub fn run(argv: &[String]) {
    let t0 = std::time::Instant::now();
    let (mut input, mut output, mut stats_path, mut fstats_path) = (None::<String>, None::<String>, None::<String>, None::<String>);
    let (mut ob_prior_paths, mut reference) = (Vec::<String>::new(), None::<String>);
    let (mut contamination_tables, mut segmentation_tables, mut contamination_estimate) = (Vec::<String>::new(), Vec::<String>::new(), CONTAMINATION_ESTIMATE);
    let (mut mito, mut min_af, mut max_alt, mut f_score_beta) = (false, 0.0f64, 1usize, 1.0f64);
    let mut i = 0;
    while i < argv.len() {
        let v = || argv.get(i + 1).cloned().unwrap_or_else(|| panic!("{USAGE}"));
        if argv[i] == "--mitochondria-mode" {
            mito = true;
            i += 1;
            continue;
        }
        match argv[i].as_str() {
            "--min-allele-fraction" => min_af = v().parse().unwrap(),
            "--max-alt-allele-count" => max_alt = v().parse().unwrap(),
            "--f-score-beta" => f_score_beta = v().parse().unwrap(),
            "-V" | "--variant" => input = Some(v()),
            "-O" | "--output" => output = Some(v()),
            "--stats" => stats_path = Some(v()),
            "--filtering-stats" => fstats_path = Some(v()),
            "--contamination-table" => contamination_tables.push(v()),
            "--contamination-estimate" => contamination_estimate = v().parse().unwrap(),
            "--tumor-segmentation" => segmentation_tables.push(v()),
            "--ob-priors" | "--orientation-bias-artifact-priors" => ob_prior_paths.push(v()),
            "-R" | "--reference" => reference = Some(v()),
            a => panic!("unknown argument {a}\n{USAGE}"),
        }
        i += 2;
    }
    let input = input.unwrap_or_else(|| panic!("{USAGE}"));
    let output = output.unwrap_or_else(|| panic!("{USAGE}"));
    let stats_path = stats_path.unwrap_or_else(|| format!("{input}.stats"));
    let fstats_path = fstats_path.unwrap_or_else(|| format!("{output}.filteringStats.tsv"));

    // input VCF
    // plain or bgzipped (by content, as htsjdk)
    let mut header: Vec<String> = Vec::new();
    let mut samples: Vec<String> = Vec::new();
    let mut recs: Vec<Rec> = Vec::new();
    for line in super::contam::open_text(&input).lines() {
        let line = line.unwrap();
        if let Some(h) = line.strip_prefix("#CHROM") {
            samples = h.split('\t').skip(9).map(|s| s.to_string()).collect();
        } else if line.starts_with("##") {
            header.push(line);
        } else if !line.is_empty() {
            recs.push(parse_record(&line, samples.len()));
        }
    }
    let normals: Vec<String> = header.iter().filter_map(|l| l.strip_prefix("##normal_sample=")).map(|s| s.to_string()).collect();
    let is_normal: Vec<bool> = samples.iter().map(|s| normals.contains(s)).collect();
    // ContaminationFilter: first record of each table, by sample
    // (ContaminationFilter and GermlineFilter collect them with Collectors.toMap: a sample twice is an error)
    let mut by_sample: HashMap<String, f64> = HashMap::new();
    for t in &contamination_tables {
        let (s, c, _) = super::contam::read_contamination(t)[0].clone();
        if let Some(c0) = by_sample.insert(s.clone(), c) {
            panic!("java.lang.IllegalStateException: Duplicate key {s} (attempted merging values {c0} and {c}): --contamination-table {t}");
        }
    }
    let contamination: Vec<f64> = samples.iter().map(|s| by_sample.get(s).copied().unwrap_or(contamination_estimate)).collect();
    let seg_tables: Vec<(String, Vec<(String, i64, i64, f64)>)> = segmentation_tables.iter().map(|t| super::contam::read_segments(t)).collect();
    for (i, (s, _)) in seg_tables.iter().enumerate() {
        if seg_tables[..i].iter().any(|(s0, _)| s0 == s) {
            panic!("java.lang.IllegalStateException: Duplicate key {s}: --tumor-segmentation {}", segmentation_tables[i]);
        }
    }
    let prior_sets = super::orientation::read_priors(&ob_prior_paths);
    let ob_priors: Vec<Option<HashMap<String, super::orientation::Prior>>> =
        samples.iter().map(|s| prior_sets.iter().rev().find(|(n, _)| n == s).map(|(_, p)| p.clone())).collect();
    let mut filters: Vec<F> = FILTERS.to_vec();
    if !ob_prior_paths.is_empty() {
        let k = filters.iter().position(|f| *f == F::PanelOfNormals).unwrap();
        filters.insert(k + 1, F::ReadOrientation);
    }
    if mito {
        // buildFiltersList: no ClusteredEvents, Multiallelic, FragmentLength, PolymeraseSlippage, FilteredHaplotype, Germline
        filters.retain(|f| {
            !matches!(f, F::ClusteredEvents | F::Multiallelic | F::FragmentLength | F::PolymeraseSlippage | F::FilteredHaplotype | F::Germline)
        });
    }
    let fasta = reference.as_ref().map(|p| crate::reference::Fasta::open(p));
    let segments: Vec<Vec<(String, i64, i64, f64)>> =
        samples.iter().map(|s| seg_tables.iter().rev().find(|(n, _)| n == s).map(|(_, v)| v.clone()).unwrap_or_default()).collect();

    // Mutect stats: callable sites
    let stats_text = std::fs::read_to_string(&stats_path).unwrap_or_else(|e| panic!("Mutect stats table {stats_path} not found ({e})"));
    let callable = stats_text.lines().filter_map(|l| l.split_once('\t')).find(|(k, _)| *k == "callable").map(|(_, v)| v.parse::<f64>().unwrap());

    // getLogSnvPrior / getLogIndelPrior: mitochondria defaults
    let log_snv_prior = if mito { -2.5 * LN_10 } else { -6.0 * LN_10 };
    let log_indel_prior = if mito { -3.75 * LN_10 } else { -7.0 * LN_10 };
    let log_artifact_prior = -1.0 * LN_10;
    let mut eng = Engine {
        model: SomaticClusteringModel::new(log_snv_prior, log_indel_prior, log_artifact_prior, callable),
        threshold: INITIAL_THRESHOLD,
        threshold_probs: Vec::new(),
        is_normal,
        filters,
        ob_priors,
        fasta,
        contamination,
        segments,
        alpha_strand: 1.0,
        beta_strand: 20.0,
        strand_artifact_prior: 0.001,
        esteps: Vec::new(),
        accumulating_phased: HashMap::new(),
        phased: HashMap::new(),
        stats: OutputStats::default(),
        min_af,
        max_alt,
        f_score_beta,
    };
    let has_alt = |r: &Rec| r.alleles[1..].iter().any(|a| a.as_slice() != b"<NON_REF>");
    // passes 0, 1: learn parameters; pass 2: threshold
    for pass in 0..3 {
        for r in recs.iter().filter(|r| has_alt(r)) {
            eng.accumulate(r);
        }
        if pass == 0 && std::env::var_os("VARIANT_RS_FMC_DUMP").is_some() {
            println!("{}", eng.model.debug_data());
        }
        if pass < 2 {
            eng.learn_parameters();
        } else {
            eng.learn_threshold();
        }
        if std::env::var_os("VARIANT_RS_FMC_DUMP").is_some() {
            let d = java_double_to_string;
            println!(
                "after pass {pass} threshold={}{} strand={}/{}/{}",
                d(eng.threshold),
                eng.model.debug_state(),
                d(eng.alpha_strand),
                d(eng.beta_strand),
                d(eng.strand_artifact_prior)
            );
        }
    }
    // pass 3: output
    // `.vcf.gz`: bgzip + .tbi, as GATK (--create-output-variant-index true)
    let cmdline = argv.join(" ");
    let level = crate::output::block_compressed(&output).then_some(5);
    let mut out = crate::output::OutFile::create(&output, &fmc_header(&header, &samples, &format!("filtermutectcalls {cmdline}"), &crate::now_string()), level)
        .unwrap_or_else(|e| panic!("{e}"));
    let mut body: Vec<u8> = Vec::new();
    let mut n_pass = 0;
    for r in &recs {
        let (filter, extra) = eng.apply(r);
        if filter == "PASS" {
            n_pass += 1;
        }
        let mut info = r.info.clone();
        for (k, v) in extra {
            match info.iter_mut().find(|(kk, _)| *kk == k) {
                Some(e) => e.1 = v,
                None => info.push((k, v)),
            }
        }
        info.sort_by(|a, b| a.0.cmp(&b.0));
        let info_s = info.iter().map(|(k, v)| match v {
            Some(v) => format!("{k}={v}"),
            None => k.clone(),
        });
        let mut line = r.head.join("\t");
        line.push('\t');
        line.push_str(&filter);
        line.push('\t');
        line.push_str(&info_s.collect::<Vec<_>>().join(";"));
        if !r.tail.is_empty() {
            line.push('\t');
            line.push_str(&r.tail);
        }
        line.push('\n');
        body.extend_from_slice(line.as_bytes());
    }
    out.add_text(&body).unwrap();
    out.finish().unwrap();
    std::fs::write(&fstats_path, eng.filtering_stats()).unwrap_or_else(|e| panic!("{fstats_path}: {e}"));
    eprintln!("filtermutectcalls: {} records, {} PASS, threshold {:.3}, {:.2}s", recs.len(), n_pass, eng.threshold, t0.elapsed().as_secs_f64());
}

/// Output header: the input lines without `##filtering_status`, plus FilterMutectCalls' FILTER / INFO lines
/// and command line, `##fileformat` first and the rest sorted by text (contig lines kept in input order).
fn fmc_header(input: &[String], samples: &[String], cmdline: &str, date: &str) -> String {
    let mut lines: Vec<String> = input.iter().filter(|l| !l.starts_with("##fileformat=") && !l.starts_with("##filtering_status=")).cloned().collect();
    for l in FMC_HEADER_LINES.lines() {
        if !lines.iter().any(|x| x == l) {
            lines.push(l.to_string());
        }
    }
    lines.push(format!(
        "##variant-rsCommandLine=<ID=FilterMutectCalls,CommandLine=\"{}\",Version=\"{}\",Date=\"{}\">",
        cmdline.replace('"', "\\\""),
        env!("CARGO_PKG_VERSION"),
        date
    ));
    let key = |l: &String| -> String {
        if l.starts_with("##contig=<") {
            "##contig=<".to_string()
        } else if let Some(rest) = l.strip_prefix("##variant-rsCommandLine").or_else(|| l.strip_prefix("##ard-rsCommandLine")) {
            // ours (and ard-rs, the name before 2026-11-15) sort where GATK's command lines do
            format!("##GATKCommandLine{rest}")
        } else {
            l.clone()
        }
    };
    lines.sort_by_key(|l| key(l));
    let mut out = String::from("##fileformat=VCFv4.2\n");
    for l in lines {
        out.push_str(&l);
        out.push('\n');
    }
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for s in samples {
        out.push('\t');
        out.push_str(s);
    }
    out.push('\n');
    out
}
