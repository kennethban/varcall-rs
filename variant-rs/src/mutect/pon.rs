//! CreateSomaticPanelOfNormals (`variant-rs createsomaticpanelofnormals`, GATK 4.7): the sites where at least
//! `--min-sample-count` normals show an artifact, read from a GenomicsDB workspace of Mutect2 calls on the normals
//! (GATK's `mutect2_pon.wdl`), with the fraction of normals showing it and a beta distribution fitted to their allele
//! fractions.
//!
//! * The records are what a walker receives from the whole workspace (`gdb_merge::whole_workspace_records`: every
//!   array over its whole range, arrays in `ChrArrayFolderComparator` order), as htsjdk parses them: a sample's AD
//!   is absent when the text is `.` or any element is (`decodeInts`).
//! * A site with no ALT, or only `*`, is skipped. At a biallelic site a sample counts when it has alt reads and its
//!   germline probability (het / hom-alt priors from the germline resource's AF, against an artifact beta-binomial
//!   prior) is below `--max-germline-probability`; at a multiallelic site every sample counts (GATK's own shortcut).
//! * FRACTION = counting samples / samples; BETA = (alpha, beta) maximizing the beta-binomial likelihood of the
//!   counting samples' (alt, ref) depths with alpha / beta fixed to the empirical ratio (Brent, `OptimizationUtils`).
//! * Output: sites only, `BETA` and `FRACTION`, a `##normal_sample` line per sample, the workspace's contigs.

use super::fmc_math::{beta_binomial_log_probability, binomial_probability, brent_max, jexp};
use crate::gdb_store::Workspace;
use crate::vqsr::vcfio::{self, Lines};
use std::collections::VecDeque;

const ARTIFACT_PRIOR: f64 = 0.001;
const ARTIFACT_ALPHA: f64 = 1.0;
const ARTIFACT_BETA: f64 = 7.0;
const NEGLIGIBLE_ALLELE_FREQUENCY: f64 = 1.0e-8;

fn refuse(why: &str) -> ! {
    crate::gatk_cli::refuse(why)
}

struct Args {
    workspace: String,
    output: String,
    reference: String,
    germline: Option<String>,
    min_sample_count: usize,
    max_germline_probability: f64,
    compression_level: i32,
}

const USAGE: &str = "usage: variant-rs createsomaticpanelofnormals -R ref.fa -V gendb://ws -O pon.vcf.gz [--germline-resource af.vcf.gz] \
[--min-sample-count 2] [--max-germline-probability 0.5] (GATK's names)";

fn parse_args(argv: &[String]) -> Args {
    let mut a = Args {
        workspace: String::new(),
        output: String::new(),
        reference: String::new(),
        germline: None,
        min_sample_count: 2,
        max_germline_probability: 0.5,
        compression_level: 6,
    };
    let val = |i: usize| -> String { argv.get(i + 1).cloned().unwrap_or_else(|| refuse(&format!("{} needs a value\n{USAGE}", argv[i]))) };
    let mut i = 0;
    while i < argv.len() {
        match argv[i].as_str() {
            "-V" | "--variant" => {
                let v = val(i);
                a.workspace = v.strip_prefix("gendb://").unwrap_or_else(|| refuse("CreateSomaticPanelOfNormals: -V must be a GenomicsDB workspace (gendb://)")).to_string();
                i += 1;
            }
            "-O" | "--output" => {
                a.output = val(i);
                i += 1;
            }
            "-R" | "--reference" => {
                a.reference = val(i);
                i += 1;
            }
            "--germline-resource" => {
                a.germline = Some(val(i));
                i += 1;
            }
            "--min-sample-count" => {
                a.min_sample_count = val(i).parse().unwrap_or_else(|_| refuse("--min-sample-count: not an integer"));
                i += 1;
            }
            "--max-germline-probability" => {
                a.max_germline_probability = vcfio::parse_vcf_double(&val(i)).unwrap_or_else(|| refuse("--max-germline-probability: not a number"));
                i += 1;
            }
            "--compression-level" => {
                a.compression_level = val(i).parse().unwrap_or_else(|_| refuse("--compression-level: not a number"));
                i += 1;
            }
            x => refuse(&format!("CreateSomaticPanelOfNormals: {x} is not supported\n{USAGE}")),
        }
        i += 1;
    }
    if a.workspace.is_empty() || a.output.is_empty() || a.reference.is_empty() {
        refuse(USAGE);
    }
    a
}

// ---------------------------------------------------------------- the germline resource

/// The germline resource read along the (sorted) sites of one contig: `FeatureContext.getValues(resource,
/// [start, end])`'s first record, i.e. the first in file order that overlaps.
struct GermlineStream {
    path: String,
    contig: String,
    lines: Option<Lines>,
    line: Vec<u8>,
    /// records read so far that may still overlap: (start, end, AF sum)
    window: VecDeque<(i64, i64, f64)>,
    peek: Option<(i64, i64, f64)>,
}

impl GermlineStream {
    fn new(path: &str) -> GermlineStream {
        GermlineStream { path: path.to_string(), contig: String::new(), lines: None, line: Vec::new(), window: VecDeque::new(), peek: None }
    }

    fn next_record(&mut self) -> Option<(i64, i64, f64)> {
        let lines = self.lines.as_mut()?;
        while lines.next_line(&mut self.line) {
            if self.line.is_empty() || self.line[0] == b'#' {
                continue;
            }
            let text = String::from_utf8_lossy(&self.line);
            let c: Vec<&str> = text.split('\t').collect();
            if c.len() < 8 || c[0] != self.contig {
                continue;
            }
            let start: i64 = c[1].parse().unwrap_or_else(|_| refuse(&format!("{}: bad POS {}", self.path, c[1])));
            let end = match vcfio::info_double(c[7].as_bytes(), "END") {
                Ok(e) if !e.is_nan() => e as i64,
                _ => start + c[3].len() as i64 - 1,
            };
            // Mutect2Engine.getAttributeAsDoubleList(vc, AF, 0.0), summed by a DoubleStream ('.' is 0)
            let afs: Vec<f64> = c[7]
                .split(';')
                .find_map(|kv| kv.strip_prefix("AF="))
                .map(|v| v.split(',').map(|x| if x == "." { 0.0 } else { vcfio::parse_vcf_double(x).unwrap_or_else(|| refuse(&format!("{}: AF {x}", self.path))) }).collect())
                .unwrap_or_default();
            return Some((start, end, super::somatic::stream_sum(afs)));
        }
        None
    }

    /// The AF sum of the first record overlapping [start, end] on `contig` (0 if none). Queries come in increasing
    /// start order within a contig.
    fn af(&mut self, contig: &str, start: i64, end: i64) -> f64 {
        if contig != self.contig {
            self.contig = contig.to_string();
            self.window.clear();
            self.peek = None;
            let tbi = format!("{}.tbi", self.path);
            self.lines = if std::path::Path::new(&tbi).exists() {
                vcfio::tbi_seek(&tbi, contig, 1).map(|v| Lines::open_at(&self.path, v))
            } else {
                Some(Lines::open(&self.path))
            };
        }
        loop {
            let r = match self.peek.take() {
                Some(r) => r,
                None => match self.next_record() {
                    Some(r) => r,
                    None => break,
                },
            };
            if r.0 > end {
                self.peek = Some(r);
                break;
            }
            self.window.push_back(r);
        }
        self.window.retain(|r| r.1 >= start);
        self.window.iter().find(|r| r.0 <= end && r.1 >= start).map_or(0.0, |r| r.2)
    }
}

// ---------------------------------------------------------------- the model

/// `germlineProbability(alleleFrequency, altCount, totalCount)`
fn germline_probability(af: f64, alt: i32, total: i32) -> f64 {
    if af < NEGLIGIBLE_ALLELE_FREQUENCY || af > 1.0 {
        return 0.0;
    }
    let het_prior = af * (1.0 - af) * 2.0;
    let hom_prior = af * af;
    let het = binomial_probability(total, alt, 0.5);
    let hom = binomial_probability(total, alt, 0.98);
    let artifact = jexp(beta_binomial_log_probability(ARTIFACT_ALPHA, ARTIFACT_BETA, total, alt));
    let probs = [het_prior * het + hom_prior * hom, ARTIFACT_PRIOR * artifact];
    let sum = probs[0] + probs[1];
    if sum < 0.0 { 0.0 } else { probs[0] / sum }
}

/// `fitBeta`: (alpha, beta)
fn fit_beta(counts: &[(i32, i32)]) -> (f64, f64) {
    let total_alt: i32 = counts.iter().map(|c| c.0).fold(0i32, |a, b| a.wrapping_add(b));
    let total_ref: i32 = counts.iter().map(|c| c.1).fold(0i32, |a, b| a.wrapping_add(b));
    let min = total_alt.min(total_ref);
    let base_alpha = (total_alt as f64 + 1.0) / (min + 1) as f64;
    let base_beta = (total_ref as f64 + 1.0) / (min + 1) as f64;
    let ll = |s: f64| -> f64 {
        let (alpha, beta) = (base_alpha * s, base_beta * s);
        super::somatic::stream_sum(counts.iter().map(|&(k, r)| beta_binomial_log_probability(alpha, beta, k + r, k)))
    };
    let scale = brent_max(ll, 0.01, 100.0, 1.0, 0.01, 0.1, 100);
    (base_alpha * scale, base_beta * scale)
}

/// A genotype's AD as htsjdk decodes it (`decodeInts`): None for '.', or when any element is not an integer.
fn ad_of(format: &[&str], sample: &str) -> Option<Vec<i32>> {
    let k = format.iter().position(|f| *f == "AD")?;
    let v = sample.split(':').nth(k)?;
    v.split(',').map(|x| x.parse::<i32>().ok()).collect()
}

pub fn run(argv: &[String]) {
    let a = parse_args(argv);
    let t0 = std::time::Instant::now();
    let ws = Workspace::open(std::path::Path::new(&a.workspace));
    let n_samples = ws.samples.len();
    let fasta = crate::reference::Fasta::open(&a.reference);
    let fmt = hc_genotype::javamath::format_vcf_double;
    // header: the tool lines, BETA / FRACTION, a ##normal_sample line per sample, the workspace's contigs
    let mut extra: Vec<String> = vec![
        "##INFO=<ID=FRACTION,Number=1,Type=Float,Description=\"Fraction of samples exhibiting artifact\">".into(),
        "##INFO=<ID=BETA,Number=2,Type=Float,Description=\"Beta distribution parameters to fit artifact allele fractions\">".into(),
        "##source=CreateSomaticPanelOfNormals".into(),
        crate::vcfheader::command_line("CreateSomaticPanelOfNormals"),
    ];
    extra.extend(ws.samples.iter().map(|s| format!("##normal_sample={s}")));
    extra.extend(ws.header.lines().filter(|l| l.starts_with("##contig=<")).map(|l| l.to_string()));
    let hdr = crate::vcfheader::write_header(&[], extra, &[]);
    let level = if crate::output::block_compressed(&a.output) { Some(a.compression_level) } else { None };
    let mut out = crate::output::OutFile::create(&a.output, &hdr, level).unwrap_or_else(|e| refuse(&format!("{}: {e}", a.output)));
    let mut germline = a.germline.as_deref().map(GermlineStream::new);
    let mut buf = String::new();
    let (mut n_in, mut n_out) = (0u64, 0u64);
    crate::gdb_merge::whole_workspace_records(&ws, &fasta, true, |line| {
        n_in += 1;
        let c: Vec<&str> = line.split('\t').collect();
        let alts: Vec<&str> = if c[4] == "." { vec![] } else { c[4].split(',').collect() };
        if alts.is_empty() || (alts.len() == 1 && alts[0] == "*") {
            return;
        }
        let pos: i64 = c[1].parse().unwrap();
        let end = pos + c[3].len() as i64 - 1;
        let af = germline.as_mut().map_or(0.0, |g| g.af(c[0], pos, end));
        let format: Vec<&str> = c.get(8).map_or(vec![], |f| f.split(':').collect());
        let ads: Vec<Option<Vec<i32>>> = c.iter().skip(9).map(|s| ad_of(&format, s)).collect();
        let alt_count = |ad: &Option<Vec<i32>>| -> i32 {
            ad.as_ref().map_or(0, |v| (v.iter().map(|&x| x as i64).sum::<i64>() - v[0] as i64) as i32)
        };
        let variant: Vec<&Option<Vec<i32>>> = if alts.len() > 1 {
            ads.iter().collect()
        } else {
            ads.iter()
                .filter(|ad| {
                    let k = alt_count(ad);
                    if k == 0 {
                        return false;
                    }
                    let total = ad.as_ref().unwrap().iter().map(|&x| x as i64).sum::<i64>() as i32;
                    germline_probability(af, k, total) < a.max_germline_probability
                })
                .collect()
        };
        if variant.len() < a.min_sample_count {
            return;
        }
        let fraction = variant.len() as f64 / n_samples as f64;
        let counts: Vec<(i32, i32)> = variant.iter().filter_map(|ad| ad.as_ref().map(|v| (alt_count(ad), v[0]))).collect();
        let (alpha, beta) = fit_beta(&counts);
        let alleles: Vec<String> = alts.iter().map(|x| crate::mergevcfs::allele(x)).collect();
        buf.push_str(&format!(
            "{}\t{}\t.\t{}\t{}\t.\t.\tBETA={},{};FRACTION={}\n",
            c[0],
            pos,
            crate::mergevcfs::allele(c[3]),
            alleles.join(","),
            fmt(alpha),
            fmt(beta),
            fmt(fraction)
        ));
        n_out += 1;
    });
    out.add_text(buf.as_bytes()).unwrap();
    out.finish().unwrap_or_else(|e| refuse(&format!("{}: {e}", a.output)));
    eprintln!("[createsomaticpanelofnormals] {n_samples} samples, {n_in} sites, {n_out} in the panel, in {:.1}s", t0.elapsed().as_secs_f64());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ad() {
        let f = ["GT", "AD", "AF"];
        assert_eq!(ad_of(&f, "0/1:10,5:0.3"), Some(vec![10, 5]));
        assert_eq!(ad_of(&f, ".:.:."), None);
        assert_eq!(ad_of(&f, "./.:10,.,3:."), None);
        assert_eq!(ad_of(&["GT"], "0/1"), None);
    }
}
