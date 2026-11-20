//! VQSR: GATK 4.7's `VariantRecalibrator` (`variant-rs vqsr`) and `ApplyVQSR` (`variant-rs applyvqsr`), with
//! GATK's outputs byte for byte.
//!
//! VariantRecalibrator is sequential in GATK and draws every random number from one `java.util.Random` seeded with
//! 47382911: 400 doubles at start-up, annotation jitter while reading the input, normalization of missing values,
//! training-set shuffling, k-means initialization, marginalization of missing annotations. All of these happen here
//! in GATK's order on the same generator; the per-datum arithmetic around them runs in parallel. Not supported
//! (refused): allele-specific mode, `--aggregate`, model input / output, scattered tranches, R plots, intervals.

pub mod applyvqsr;
mod gmm;
mod jama;
mod tranche;
pub mod vcfio;

use crate::intervals::Dict;
use hc_assembly::jmath;
use hc_genotype::javamath::JavaRandom;
use rayon::prelude::*;
use std::collections::HashMap;
use vcfio::{Lines, VType};

/// `VariantRecalibratorArgumentCollection.Mode`
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Mode {
    Snp,
    Indel,
    Both,
}

impl Mode {
    pub fn parse(s: &str) -> Mode {
        match s {
            "SNP" => Mode::Snp,
            "INDEL" => Mode::Indel,
            "BOTH" => Mode::Both,
            _ => crate::gatk_cli::refuse(&format!("-mode {s}: expected SNP, INDEL or BOTH")),
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            Mode::Snp => "SNP",
            Mode::Indel => "INDEL",
            Mode::Both => "BOTH",
        }
    }
    /// `VariantDataManager.checkVariationClass(evalVC, mode)`
    pub fn accepts(self, t: VType) -> bool {
        match self {
            Mode::Snp => matches!(t, VType::Snp | VType::Mnp),
            Mode::Indel => matches!(t, VType::Indel | VType::Mixed | VType::Symbolic),
            Mode::Both => true,
        }
    }
}

/// `VariantDatum`
#[derive(Clone)]
pub struct Datum {
    pub ann: Vec<f64>,
    pub is_null: Vec<bool>,
    pub is_known: bool,
    pub lod: f64,
    pub at_truth: bool,
    pub at_training: bool,
    pub at_anti_training: bool,
    pub is_transition: bool,
    pub is_snp: bool,
    pub failing_std: bool,
    pub prior: f64,
    pub tid: usize,
    pub start: i64,
    pub end: i64,
    pub worst: i32,
    /// (variant type, for the training-set class check)
    vtype: VType,
}

/// A `--resource:name,known=..,training=..,truth=..,bad=..,prior=..` input (`TrainingSet`).
#[derive(Clone, Debug)]
pub struct Resource {
    pub name: String,
    pub path: String,
    pub known: bool,
    pub training: bool,
    pub truth: bool,
    pub anti_training: bool,
    pub prior: f64,
}

impl Resource {
    /// From GATK's tagged argument: `--resource:<name>[,key=value...]` and the file.
    pub fn parse(tag: &str, path: &str) -> Resource {
        let mut parts = tag.split(',');
        let name = parts.next().unwrap_or("").to_string();
        let mut r = Resource { name, path: path.to_string(), known: false, training: false, truth: false, anti_training: false, prior: 0.0 };
        for kv in parts {
            let (k, v) = kv.split_once('=').unwrap_or_else(|| crate::gatk_cli::refuse(&format!("--resource:{tag}: {kv} is not key=value")));
            match k {
                "known" => r.known = v == "true",
                "training" => r.training = v == "true",
                "truth" => r.truth = v == "true",
                "bad" => r.anti_training = v == "true",
                "consensus" => {
                    if v == "true" {
                        crate::gatk_cli::refuse(&format!("--resource:{tag}: consensus tracks are not supported"));
                    }
                }
                "prior" => r.prior = vcfio::parse_vcf_double(v).filter(|x| !x.is_nan() || v == "NaN").unwrap_or_else(|| crate::gatk_cli::refuse(&format!("Malformed floating point valueprior"))),
                _ => {} // other tags (e.g. a file type) do not change the result
            }
        }
        r
    }
}

pub struct RecalArgs {
    pub input: String,
    pub output: String,
    pub tranches_file: String,
    pub reference: Option<String>,
    pub resources: Vec<Resource>,
    pub annotations: Vec<String>,
    pub mode: Mode,
    pub tranches: Vec<f64>,
    pub ignore_filters: Vec<String>,
    pub ignore_all_filters: bool,
    pub trust_all_polymorphic: bool,
    pub max_gaussians: usize,
    pub max_negative_gaussians: usize,
    pub std_threshold: f64,
    pub max_training: usize,
    pub min_bad: usize,
    pub bad_lod_cutoff: f64,
    pub mq_cap: i32,
    pub mq_jitter: f64,
    pub settings: gmm::Settings,
    pub compression_level: i32,
}

const USAGE: &str = "usage: variant-rs vqsr -V in.vcf.gz -O out.recal.vcf.gz --tranches-file out.tranches -mode SNP|INDEL|BOTH \
-an ANN [-an ANN ...] --resource:name,known=B,training=B,truth=B,prior=Q file [...] [-R ref.fa] [-tranche T ...] \
[--trust-all-polymorphic] [--max-gaussians 8] [--max-negative-gaussians 2] [--ignore-filter F] [--ignore-all-filters] (GATK's names)";

fn need<'a>(argv: &'a [String], i: usize) -> &'a str {
    argv.get(i + 1).map(|s| s.as_str()).unwrap_or_else(|| crate::gatk_cli::refuse(&format!("{} needs a value\n{USAGE}", argv[i])))
}

fn num<T: std::str::FromStr>(argv: &[String], i: usize) -> T {
    need(argv, i).parse::<T>().unwrap_or_else(|_| crate::gatk_cli::refuse(&format!("{}: not a number: {}", argv[i], need(argv, i))))
}

fn bool_val(v: &str) -> bool {
    match v.to_ascii_lowercase().as_str() {
        "true" | "t" | "yes" | "y" => true,
        "false" | "f" | "no" | "n" => false,
        _ => crate::gatk_cli::refuse(&format!("not a boolean: {v}")),
    }
}

/// A boolean option given as a flag or with a value
fn flag(argv: &[String], i: &mut usize) -> bool {
    match argv.get(*i + 1).map(|s| s.to_ascii_lowercase()) {
        Some(v) if ["true", "false", "t", "f", "yes", "no", "y", "n"].contains(&v.as_str()) => {
            *i += 1;
            bool_val(&v)
        }
        _ => true,
    }
}

fn parse_args(argv: &[String]) -> RecalArgs {
    let mut a = RecalArgs {
        input: String::new(),
        output: String::new(),
        tranches_file: String::new(),
        reference: None,
        resources: vec![],
        annotations: vec![],
        mode: Mode::Snp,
        tranches: vec![],
        ignore_filters: vec![],
        ignore_all_filters: false,
        trust_all_polymorphic: false,
        max_gaussians: 8,
        max_negative_gaussians: 2,
        std_threshold: 10.0,
        max_training: 2_500_000,
        min_bad: 1000,
        bad_lod_cutoff: -5.0,
        mq_cap: 0,
        mq_jitter: 0.05,
        settings: gmm::Settings { max_iterations: 150, kmeans_iterations: 100, shrinkage: 1.0, dirichlet: 0.001, prior_counts: 20.0 },
        compression_level: 6,
    };
    let mut mode_given = false;
    let mut i = 0;
    while i < argv.len() {
        let arg = argv[i].as_str();
        if let Some(tag) = arg.strip_prefix("--resource:").or_else(|| arg.strip_prefix("-resource:")) {
            a.resources.push(Resource::parse(tag, need(argv, i)));
            i += 2;
            continue;
        }
        match arg {
            "-V" | "--variant" => {
                if !a.input.is_empty() {
                    crate::gatk_cli::refuse("VariantRecalibrator with several -V inputs is not supported");
                }
                a.input = need(argv, i).to_string();
                i += 1;
            }
            "-O" | "--output" => {
                a.output = need(argv, i).to_string();
                i += 1;
            }
            "--tranches-file" => {
                a.tranches_file = need(argv, i).to_string();
                i += 1;
            }
            "-R" | "--reference" => {
                a.reference = Some(need(argv, i).to_string());
                i += 1;
            }
            "-an" | "--use-annotation" => {
                a.annotations.push(need(argv, i).to_string());
                i += 1;
            }
            "-mode" | "--mode" => {
                a.mode = Mode::parse(need(argv, i));
                mode_given = true;
                i += 1;
            }
            "-tranche" | "--truth-sensitivity-tranche" => {
                let t = vcfio::parse_vcf_double(need(argv, i)).unwrap_or_else(|| crate::gatk_cli::refuse(&format!("-tranche {}: not a number", need(argv, i))));
                a.tranches.push(t);
                i += 1;
            }
            "--ignore-filter" => {
                a.ignore_filters.push(need(argv, i).to_string());
                i += 1;
            }
            "--ignore-all-filters" => a.ignore_all_filters = flag(argv, &mut i),
            "--trust-all-polymorphic" => a.trust_all_polymorphic = flag(argv, &mut i),
            "--max-gaussians" => {
                a.max_gaussians = num(argv, i);
                i += 1;
            }
            "--max-negative-gaussians" => {
                a.max_negative_gaussians = num(argv, i);
                i += 1;
            }
            "--max-iterations" => {
                a.settings.max_iterations = num(argv, i);
                i += 1;
            }
            "--k-means-iterations" => {
                a.settings.kmeans_iterations = num(argv, i);
                i += 1;
            }
            "-std" | "--standard-deviation-threshold" => {
                a.std_threshold = num(argv, i);
                i += 1;
            }
            "--shrinkage" => {
                a.settings.shrinkage = num(argv, i);
                i += 1;
            }
            "--dirichlet" => {
                a.settings.dirichlet = num(argv, i);
                i += 1;
            }
            "--prior-counts" => {
                a.settings.prior_counts = num(argv, i);
                i += 1;
            }
            "--maximum-training-variants" => {
                a.max_training = num(argv, i);
                i += 1;
            }
            "--minimum-bad-variants" => {
                a.min_bad = num(argv, i);
                i += 1;
            }
            "-bad-lod-cutoff" | "--bad-lod-score-cutoff" => {
                a.bad_lod_cutoff = num(argv, i);
                i += 1;
            }
            "-mq-cap" | "--mq-cap-for-logit-jitter-transform" => {
                a.mq_cap = num(argv, i);
                i += 1;
            }
            "--mq-jitter" => {
                a.mq_jitter = num(argv, i);
                i += 1;
            }
            "-titv" | "--target-titv" => i += 1, // for plots only
            "--compression-level" => {
                a.compression_level = num(argv, i);
                i += 1;
            }
            x => crate::gatk_cli::refuse(&format!("VariantRecalibrator: {x} is not supported\n{USAGE}")),
        }
        i += 1;
    }
    if a.input.is_empty() || a.output.is_empty() || a.tranches_file.is_empty() || a.annotations.is_empty() || !mode_given {
        crate::gatk_cli::refuse(USAGE);
    }
    if a.tranches.is_empty() {
        a.tranches = vec![100.0, 99.9, 99.0, 90.0];
    }
    if !a.resources.iter().any(|r| r.training) {
        crate::gatk_cli::refuse("No training set found! Please provide sets of known polymorphic loci marked with the training=true feature input tag. For example, -resource hapmap,VCF,known=false,training=true,truth=true,prior=12.0 hapmapFile.vcf");
    }
    if !a.resources.iter().any(|r| r.truth) {
        crate::gatk_cli::refuse("No truth set found! Please provide sets of known polymorphic loci marked with the truth=true feature input tag. For example, -resource hapmap,VCF,known=false,training=true,truth=true,prior=12.0 hapmapFile.vcf");
    }
    if a.settings.kmeans_iterations == 0 || a.max_gaussians == 0 || a.max_negative_gaussians == 0 {
        crate::gatk_cli::refuse("--max-gaussians, --max-negative-gaussians and --k-means-iterations must be positive");
    }
    // `VariantDataManager`: duplicate annotations are used once
    let mut seen = std::collections::HashSet::new();
    a.annotations.retain(|x| seen.insert(x.clone()));
    if a.annotations.len() > 16 {
        crate::gatk_cli::refuse("more than 16 annotations are not supported");
    }
    a
}

/// The sequence dictionary: the reference's when given, else the input VCF's `##contig` lines.
fn dictionary(reference: &Option<String>, header: &[String]) -> Dict {
    if let Some(r) = reference {
        return Dict::from_reference(r);
    }
    let (mut names, mut lens) = (vec![], vec![]);
    for l in header {
        if let Some(rest) = l.strip_prefix("##contig=<") {
            let mut id = None;
            let mut len = 0i64;
            for kv in rest.trim_end_matches('>').split(',') {
                if let Some(v) = kv.strip_prefix("ID=") {
                    id = Some(v.to_string());
                } else if let Some(v) = kv.strip_prefix("length=") {
                    len = v.parse().unwrap_or(0);
                }
            }
            if let Some(id) = id {
                names.push(id);
                lens.push(len);
            }
        }
    }
    if names.is_empty() {
        crate::gatk_cli::refuse("the input VCF has no ##contig lines: give the reference with -R");
    }
    Dict::new(names, lens)
}

/// `logitTransform`
fn logit(x: f64, xmin: f64, xmax: f64) -> f64 {
    jmath::log((x - xmin) / (xmax - x))
}

/// `MathUtils.compareDoubles(a, b, epsilon) == 0`
fn near(a: f64, b: f64, eps: f64) -> bool {
    (a - b).abs() < eps
}

/// `VariantDataManager.decodeAnnotation` with jitter (as `addDatum` always calls it).
fn decode_annotation(key: &str, info: &[u8], a: &RecalArgs, rng: &mut JavaRandom) -> Result<f64, String> {
    const PRECISION: f64 = 0.01;
    const SAFETY_OFFSET: f64 = 0.01;
    const LOG_OF_TWO: f64 = 0.6931472;
    let mut value = vcfio::info_double(info, key)?;
    if value.is_infinite() {
        value = f64::NAN;
    }
    let k = key.to_ascii_lowercase();
    if k == "haplotypescore" && near(value, 0.0, PRECISION) {
        value += 0.01 * rng.next_gaussian();
    }
    if (k == "fs" || k == "as_filterstatus") && near(value, 0.0, PRECISION) {
        value += 0.01 * rng.next_gaussian();
    }
    if k == "inbreedingcoeff" && near(value, 0.0, PRECISION) {
        value += 0.01 * rng.next_gaussian();
    }
    if (k == "sor" || k == "as_sor") && near(value, LOG_OF_TWO, PRECISION) {
        value += 0.01 * rng.next_gaussian();
    }
    if k == "mq" {
        if a.mq_cap > 0 {
            let cap = a.mq_cap as f64;
            value = logit(value, -SAFETY_OFFSET, cap + SAFETY_OFFSET);
            if near(value, logit(cap, -SAFETY_OFFSET, cap + SAFETY_OFFSET), PRECISION) {
                value += a.mq_jitter * rng.next_gaussian();
            }
        } else if near(value, a.mq_cap as f64, PRECISION) {
            value += a.mq_jitter * rng.next_gaussian();
        }
    }
    if k == "as_mq" {
        value += a.mq_jitter * rng.next_gaussian();
    }
    Ok(value)
}

/// One training-set record at a datum's start, as `parseTrainingSets` needs it.
#[derive(Clone, Copy)]
struct TrainRec {
    filtered: bool,
    vtype: VType,
    /// hasGenotypes() && !isPolymorphicInSamples() (only computed without --trust-all-polymorphic)
    monomorphic: bool,
}

/// The records of one resource starting at each key (dictionary index, start), keys sorted.
fn scan_resource(r: &Resource, dict: &Dict, keys: &[(usize, i64)], need_genotypes: bool) -> Vec<Vec<TrainRec>> {
    let mut out: Vec<Vec<TrainRec>> = vec![Vec::new(); keys.len()];
    let mut lines = Lines::open(&r.path);
    let mut line = Vec::new();
    let mut k = 0usize;
    let mut last: (usize, i64) = (0, 0);
    let mut cur_tid: Option<(Vec<u8>, Option<usize>)> = None;
    while k < keys.len() && lines.next_line(&mut line) {
        if line.is_empty() || line[0] == b'#' {
            continue;
        }
        let tab1 = memchr::memchr(b'\t', &line).unwrap_or_else(|| panic!("{}: malformed VCF line", r.path));
        let chrom = &line[..tab1];
        let tid = match &cur_tid {
            Some((c, t)) if c.as_slice() == chrom => *t,
            _ => {
                let t = dict.index.get(std::str::from_utf8(chrom).unwrap_or("")).copied();
                cur_tid = Some((chrom.to_vec(), t));
                t
            }
        };
        let Some(tid) = tid else { continue };
        let rest = &line[tab1 + 1..];
        let tab2 = memchr::memchr(b'\t', rest).unwrap_or(rest.len());
        let pos: i64 = std::str::from_utf8(&rest[..tab2]).ok().and_then(|s| s.parse().ok()).unwrap_or_else(|| panic!("{}: bad POS", r.path));
        let key = (tid, pos);
        if key < last {
            crate::gatk_cli::refuse(&format!("{}: records are not sorted in the order of the sequence dictionary (at {}:{pos})", r.path, String::from_utf8_lossy(chrom)));
        }
        last = key;
        while k < keys.len() && keys[k] < key {
            k += 1;
        }
        if k == keys.len() || keys[k] != key {
            continue;
        }
        // the sample columns only when polymorphism is checked
        let cols: Vec<&[u8]> = if need_genotypes { line.split(|&c| c == b'\t').collect() } else { line.splitn(9, |&c| c == b'\t').collect() };
        if cols.len() < 8 {
            panic!("{}: VCF line with fewer than 8 columns", r.path);
        }
        let vtype = vcfio::variant_type(cols[3], cols[4]);
        let filtered = vcfio::filters(cols[6]).is_some();
        let monomorphic = need_genotypes && vtype != VType::NoVariation && !filtered && !vcfio::polymorphic_in_samples(&cols);
        out[k].push(TrainRec { filtered, vtype, monomorphic });
    }
    out
}

/// `VariantDataManager.checkVariationClass(evalVC, trainVC)`
fn class_matches(eval: VType, train: VType) -> bool {
    match train {
        VType::Snp | VType::Mnp => Mode::Snp.accepts(eval),
        VType::Indel | VType::Mixed | VType::Symbolic => Mode::Indel.accepts(eval),
        VType::NoVariation => false,
    }
}

/// `VariantDataManager` after reading: means, standard deviations and annotation order.
struct Normalized {
    keys: Vec<String>,
    mean: Vec<f64>,
    std: Vec<f64>,
}

/// `VariantDataManager.normalizeData(true, null)`
fn normalize(data: &mut [Datum], keys: &[String], std_threshold: f64, rng: &mut JavaRandom) -> Result<Normalized, String> {
    let n = keys.len();
    let mean_of = |data: &[Datum], i: usize, training: bool| -> f64 {
        let (mut sum, mut cnt) = (0.0, 0usize);
        for d in data {
            if training == d.at_training && !d.is_null[i] {
                sum += d.ann[i];
                cnt += 1;
            }
        }
        sum / cnt as f64
    };
    let mut mean = vec![0.0; n];
    let mut std = vec![0.0; n];
    let mut zero_variance = false;
    for i in 0..n {
        let m = mean_of(data, i, true);
        let (mut sum, mut cnt) = (0.0, 0usize);
        for d in data.iter() {
            if d.at_training && !d.is_null[i] {
                sum += (d.ann[i] - m) * (d.ann[i] - m);
                cnt += 1;
            }
        }
        let s = (sum / cnt as f64).sqrt();
        if m.is_nan() {
            return Err(format!("Values for {} annotation not detected for ANY training variant in the input callset. VariantAnnotator may be used to add these annotations.", keys[i]));
        }
        zero_variance = zero_variance || s < 1e-5;
        mean[i] = m;
        std[i] = s;
        eprintln!("[vqsr] {}: \t mean = {}\t standard deviation = {}", keys[i], hc_genotype::javamath::format_fixed(m, 2), hc_genotype::javamath::format_fixed(s, 2));
        for d in data.iter_mut() {
            d.ann[i] = if d.is_null[i] { 0.1 * rng.next_gaussian() } else { (d.ann[i] - m) / s };
        }
    }
    if zero_variance {
        return Err("Found annotations with zero variance. They must be excluded before proceeding.".into());
    }
    data.par_iter_mut().for_each(|d| d.failing_std = d.ann.iter().any(|v| v.abs() > std_threshold));
    // calculateSortOrder: by -|training mean - mean of the other (normalized) data|, stable
    let mut order: Vec<(f64, usize)> = (0..n).map(|i| (-1.0 * (mean[i] - mean_of(data, i, false)).abs(), i)).collect();
    order.sort_by(|a, b| hc_assembly::jmath::java_cmp(a.0, b.0));
    let order: Vec<usize> = order.into_iter().map(|x| x.1).collect();
    let keys: Vec<String> = order.iter().map(|&i| keys[i].clone()).collect();
    let mean: Vec<f64> = order.iter().map(|&i| mean[i]).collect();
    let std: Vec<f64> = order.iter().map(|&i| std[i]).collect();
    data.par_iter_mut().for_each(|d| {
        d.ann = order.iter().map(|&i| d.ann[i]).collect();
        d.is_null = order.iter().map(|&i| d.is_null[i]).collect();
    });
    eprintln!("[vqsr] Annotation order is: [{}]", keys.join(", "));
    Ok(Normalized { keys, mean, std })
}

/// `VariantDataManager.getTrainingData`: indices of training data within the standard-deviation threshold
/// (shuffled and cut to `max` when there are more).
fn training_data(data: &[Datum], a: &RecalArgs, rng: &mut JavaRandom) -> Vec<usize> {
    let mut t: Vec<usize> = (0..data.len()).filter(|&i| data[i].at_training && !data[i].failing_std).collect();
    eprintln!("[vqsr] Training with {} variants after standard deviation thresholding.", t.len());
    if t.len() < a.min_bad {
        eprintln!("[vqsr] WARNING: Training with very few variant sites! Please check the model reporting PDF to ensure the quality of the model is reliable.");
    } else if t.len() > a.max_training {
        eprintln!("[vqsr] WARNING: Very large training set detected. Downsampling to {} training variants.", a.max_training);
        rng.shuffle(&mut t);
        t.truncate(a.max_training);
    }
    t
}

/// `VariantDataManager.selectWorstVariants`: marks and returns the data with lod below the cutoff.
fn select_worst(data: &mut [Datum], cutoff: f64) -> Vec<usize> {
    let mut out = Vec::new();
    for (i, d) in data.iter_mut().enumerate() {
        if !d.failing_std && !d.lod.is_infinite() && d.lod < cutoff {
            d.at_anti_training = true;
            out.push(i);
        }
    }
    eprintln!("[vqsr] Selected worst {} scoring variants --> variants with LOD <= {}.", out.len(), hc_genotype::javamath::format_fixed(cutoff, 4));
    out
}

/// `VariantRecalibratorEngine.calculateWorstPerformingAnnotation`
fn worst_annotations(data: &mut [Datum], good: &gmm::Gmm, bad: &gmm::Gmm) {
    data.par_iter_mut().for_each(|d| {
        let mut worst = -1i32;
        let mut min_prob = f64::MAX;
        for i in 0..d.ann.len() {
            if let (Some(g), Some(b)) = (good.evaluate_in_one_dimension(d, i), bad.evaluate_in_one_dimension(d, i)) {
                let p = g - b;
                if p < min_prob {
                    min_prob = p;
                    worst = i as i32;
                }
            }
        }
        d.worst = worst;
    });
}

/// The input's header lines and data, with the datums GATK's `addDatum` makes (annotations decoded and jittered
/// in input order).
fn read_input(a: &RecalArgs, rng: &mut JavaRandom) -> (Vec<String>, Dict, Vec<Datum>) {
    let mut lines = Lines::open(&a.input);
    let mut line = Vec::new();
    let mut header = Vec::new();
    loop {
        if !lines.next_line(&mut line) {
            crate::gatk_cli::refuse(&format!("{}: no #CHROM line", a.input));
        }
        let l = String::from_utf8_lossy(&line).into_owned();
        if l.starts_with("##") {
            header.push(l);
        } else {
            break;
        }
    }
    let dict = dictionary(&a.reference, &header);
    let ignore: std::collections::HashSet<&[u8]> = a.ignore_filters.iter().map(|s| s.as_bytes()).collect();
    let mut data = Vec::new();
    let mut cur: Option<(Vec<u8>, usize)> = None;
    while lines.next_line(&mut line) {
        if line.is_empty() || line[0] == b'#' {
            continue;
        }
        let cols: Vec<&[u8]> = line.splitn(9, |&c| c == b'\t').collect();
        if cols.len() < 8 {
            crate::gatk_cli::refuse(&format!("{}: VCF line with fewer than 8 columns", a.input));
        }
        if let Some(f) = vcfio::filters(cols[6]) {
            if !a.ignore_all_filters && !f.iter().all(|x| ignore.contains(x)) {
                continue;
            }
        }
        let vtype = vcfio::variant_type(cols[3], cols[4]);
        if !a.mode.accepts(vtype) {
            continue;
        }
        let tid = match &cur {
            Some((c, t)) if c.as_slice() == cols[0] => *t,
            _ => {
                let name = String::from_utf8_lossy(cols[0]).into_owned();
                let t = *dict.index.get(&name).unwrap_or_else(|| crate::gatk_cli::refuse(&format!("{}: contig {name} is not in the sequence dictionary", a.input)));
                cur = Some((cols[0].to_vec(), t));
                t
            }
        };
        let start: i64 = std::str::from_utf8(cols[1]).ok().and_then(|s| s.parse().ok()).unwrap_or_else(|| crate::gatk_cli::refuse(&format!("{}: bad POS", a.input)));
        let mut ann = Vec::with_capacity(a.annotations.len());
        let mut is_null = Vec::with_capacity(a.annotations.len());
        for key in &a.annotations {
            let v = decode_annotation(key, cols[7], a, rng).unwrap_or_else(|e| crate::gatk_cli::refuse(&format!("{}:{start}: {e}", String::from_utf8_lossy(cols[0]))));
            is_null.push(v.is_nan());
            ann.push(v);
        }
        // htsjdk: END from INFO when present, else POS + len(REF) - 1
        let end = match vcfio::info_double(cols[7], "END") {
            Ok(e) if !e.is_nan() => e as i64,
            _ => start + cols[3].len() as i64 - 1,
        };
        let biallelic = !cols[4].contains(&b',');
        let is_snp = vtype == VType::Snp && biallelic;
        let is_transition = is_snp && {
            let (r, x) = (cols[3][0].to_ascii_uppercase(), cols[4][0].to_ascii_uppercase());
            matches!((r, x), (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C'))
        };
        data.push(Datum {
            ann,
            is_null,
            is_known: false,
            lod: 0.0,
            at_truth: false,
            at_training: false,
            at_anti_training: false,
            is_transition,
            is_snp,
            failing_std: false,
            prior: 2.0,
            tid,
            start,
            end,
            worst: -1,
            vtype,
        });
    }
    (header, dict, data)
}

/// `parseTrainingSets` for every datum, then the prior as a log10 odds.
fn label_training(data: &mut [Datum], a: &RecalArgs, dict: &Dict) {
    let mut keys: Vec<(usize, i64)> = data.iter().map(|d| (d.tid, d.start)).collect();
    if keys.windows(2).any(|w| w[1] < w[0]) {
        crate::gatk_cli::refuse(&format!("{}: records are not sorted in the order of the sequence dictionary", a.input));
    }
    keys.dedup();
    let need_gt = !a.trust_all_polymorphic;
    for r in &a.resources {
        eprintln!("[vqsr] Found {} track: \tKnown = {} \tTraining = {} \tTruth = {} \tPrior = Q{}", r.name, r.known, r.training, r.truth, hc_genotype::javamath::format_fixed(r.prior, 1));
    }
    let found: Vec<Vec<Vec<TrainRec>>> = a.resources.par_iter().map(|r| scan_resource(r, dict, &keys, need_gt)).collect();
    let key_index: HashMap<(usize, i64), usize> = keys.iter().enumerate().map(|(i, &k)| (k, i)).collect();
    data.par_iter_mut().for_each(|d| {
        let k = key_index[&(d.tid, d.start)];
        for (r, recs) in a.resources.iter().zip(&found) {
            for t in &recs[k] {
                let valid = !t.filtered && t.vtype != VType::NoVariation && class_matches(d.vtype, t.vtype) && (a.trust_all_polymorphic || !t.monomorphic);
                if valid {
                    d.is_known = d.is_known || r.known;
                    d.at_truth = d.at_truth || r.truth;
                    d.at_training = d.at_training || r.training;
                    d.prior = d.prior.max(r.prior);
                }
                d.at_anti_training = d.at_anti_training || r.anti_training;
            }
        }
        // QualityUtils.qualToProb, then the log10 odds
        assert!(d.prior >= 0.0, "qual must be >= 0.0 but got {}", d.prior);
        let p = 1.0 - jmath::pow(10.0, d.prior / -10.0);
        d.prior = jmath::log10(p) - jmath::log10(1.0 - p);
    });
}

fn fail(msg: &str) -> ! {
    eprintln!("variant-rs vqsr: {msg}");
    std::process::exit(3);
}

pub fn run(argv: &[String]) {
    let a = parse_args(argv);
    let t0 = std::time::Instant::now();
    let mut rng = JavaRandom::new(hc_genotype::javamath::GATK_RANDOM_SEED);
    // onTraversalStart: 2 * REPLICATE (200) doubles for `replicate`, unused but drawn
    for _ in 0..400 {
        rng.next_double();
    }
    let (in_header, dict, mut data) = read_input(&a, &mut rng);
    eprintln!("[vqsr] {} {} records in {:.1}s", data.len(), a.mode.name(), t0.elapsed().as_secs_f64());
    label_training(&mut data, &a, &dict);
    eprintln!("[vqsr] training sets matched in {:.1}s", t0.elapsed().as_secs_f64());
    let norm = normalize(&mut data, &a.annotations, a.std_threshold, &mut rng).unwrap_or_else(|e| fail(&e));
    let positive = training_data(&data, &a, &mut rng);
    let good = gmm::generate_model(&data, &positive, a.max_gaussians, &a.settings, &mut rng).unwrap_or_else(|e| fail(&e));
    let mut good = good;
    gmm::evaluate_data(&mut data, &mut good, false, &mut rng);
    if good.failed_to_converge {
        fail("Positive training model failed to converge.  One or more annotations (usually MQ) may have insufficient variance.  Please consider lowering the maximum number of Gaussians allowed for use in the model (via --max-gaussians 4, for example).");
    }
    let negative = select_worst(&mut data, a.bad_lod_cutoff);
    let mut bad = gmm::generate_model(&data, &negative, a.max_negative_gaussians.min(a.max_gaussians), &a.settings, &mut rng).unwrap_or_else(|e| fail(&e));
    if bad.failed_to_converge {
        fail("NaN LOD value assigned. Clustering with this few variants and these annotations is unsafe. Please consider raising the number of variants used to train the negative model (via --minimum-bad-variants 5000, for example).");
    }
    gmm::evaluate_data(&mut data, &mut bad, true, &mut rng);
    worst_annotations(&mut data, &good, &bad);
    eprintln!("[vqsr] models in {:.1}s", t0.elapsed().as_secs_f64());
    // tranches (sorts the data by lod), then the recal file in locus order
    let text = tranche::tranches_text(&mut data, &a.tranches, a.mode).unwrap_or_else(|e| fail(&e));
    std::fs::write(&a.tranches_file, text).unwrap_or_else(|e| fail(&format!("{}: {e}", a.tranches_file)));
    data.sort_by(|x, y| (x.tid, x.start, x.end).cmp(&(y.tid, y.start, y.end)));
    write_recal(&a, &in_header, &dict, &norm, &data);
    let _ = &norm.mean;
    let _ = &norm.std;
    eprintln!("[vqsr] done in {:.1}s", t0.elapsed().as_secs_f64());
}

/// `writeOutRecalibrationTable`: one `N` / `<VQSR>` record per datum, in locus order.
fn write_recal(a: &RecalArgs, _in_header: &[String], dict: &Dict, norm: &Normalized, data: &[Datum]) {
    let header = recal_header(a, dict);
    let level = if crate::output::block_compressed(&a.output) { Some(a.compression_level) } else { None };
    let mut of = crate::output::OutFile::create(&a.output, &header, level).unwrap_or_else(|e| fail(&format!("{}: {e}", a.output)));
    let mut buf = Vec::with_capacity(8 << 20);
    for d in data {
        use std::io::Write;
        let culprit = if d.worst >= 0 { norm.keys[d.worst as usize].as_str() } else { "NULL" };
        write!(buf, "{}\t{}\t.\tN\t<VQSR>\t.\t.\tEND={}", dict.names[d.tid], d.start, d.end).unwrap();
        if d.at_anti_training {
            buf.extend_from_slice(b";NEGATIVE_TRAIN_SITE");
        }
        if d.at_training {
            buf.extend_from_slice(b";POSITIVE_TRAIN_SITE");
        }
        writeln!(buf, ";VQSLOD={};culprit={}", hc_genotype::javamath::format_fixed(d.lod, 4), culprit).unwrap();
        if buf.len() >= 8 << 20 {
            of.add_text(&buf).unwrap();
            buf.clear();
        }
    }
    of.add_text(&buf).unwrap();
    of.finish().unwrap_or_else(|e| fail(&format!("{}: {e}", a.output)));
}

/// The recal VCF header: the VQSR lines, `##source`, the command line, sorted as htsjdk writes them (our command line
/// where GATK's `##GATKCommandLine` sorts), then the dictionary's contigs.
fn recal_header(a: &RecalArgs, dict: &Dict) -> String {
    let cmdline: Vec<String> = std::env::args().collect();
    let mut lines = vec![
        "##FILTER=<ID=PASS,Description=\"Site contains at least one allele that passes filters\">".to_string(),
        "##source=VariantRecalibrator".to_string(),
        format!(
            "##variant-rsCommandLine=<ID=VariantRecalibrator,CommandLine=\"{}\",Version=\"{}\",Date=\"{}\">",
            cmdline.join(" ").replace('"', "\\\""),
            env!("CARGO_PKG_VERSION"),
            crate::now_string()
        ),
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">".to_string(),
        "##INFO=<ID=NEGATIVE_TRAIN_SITE,Number=0,Type=Flag,Description=\"This variant was used to build the negative training set of bad variants\">".to_string(),
        "##INFO=<ID=POSITIVE_TRAIN_SITE,Number=0,Type=Flag,Description=\"This variant was used to build the positive training set of good variants\">".to_string(),
        "##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=\"Log odds of being a true variant versus being false under the trained gaussian mixture model\">".to_string(),
        "##INFO=<ID=culprit,Number=1,Type=String,Description=\"The annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out\">".to_string(),
    ];
    // with a reference, GATK's `VcfUtils.updateHeaderContigLines`: the fasta's file name as each contig's assembly,
    // and `##reference=` with that name up to its last '.'
    let fasta_name = a.reference.as_deref().map(|r| std::path::Path::new(r).file_name().map(|f| f.to_string_lossy().into_owned()).unwrap_or_default());
    for (n, l) in dict.names.iter().zip(&dict.lens) {
        match &fasta_name {
            Some(f) => lines.push(format!("##contig=<ID={n},length={l},assembly={f}>")),
            None => lines.push(format!("##contig=<ID={n},length={l}>")),
        }
    }
    if let Some(f) = &fasta_name {
        lines.push(format!("##reference={}", f.rsplit_once('.').map(|(s, _)| s).unwrap_or(f)));
    }
    // htsjdk: lines sorted as text, the contig lines as one block (in dictionary order) under "##contig=<"
    let key = |l: &String| -> String {
        if l.starts_with("##contig=<") {
            "##contig=<".into()
        } else {
            l.strip_prefix("##variant-rsCommandLine").map(|r| format!("##GATKCommandLine{r}")).unwrap_or_else(|| l.clone())
        }
    };
    lines.sort_by_key(key);
    let _ = a;
    let mut out = String::from("##fileformat=VCFv4.2\n");
    for l in lines {
        out.push_str(&l);
        out.push('\n');
    }
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    out
}
