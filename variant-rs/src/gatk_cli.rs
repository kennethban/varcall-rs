//! GATK-compatible command line: `variant-rs <Tool> [GATK arguments]`, as `gatk <Tool> ...` would be run.
//!
//! Every GATK 4.7 argument of a supported tool (tables generated from `gatk <Tool> --help`, see
//! `gatkargs/gen.py`) is in one of four classes:
//! - mapped: translated to the variant-rs tool's own argument;
//! - ignored: has no effect on the output records (temporary directories, logging, threads, index creation,
//!   cloud buffers, ...);
//! - default-only: accepted when its value equals GATK's default (so scripts that spell out defaults still run),
//!   refused otherwise, because it would change the output;
//! - unknown to GATK: refused.
//! variant-rs' own options (`--threads`, `--pairhmm`, `--shard`, ...) are accepted as well. With
//! `GATK_FALLBACK=/path/to/gatk`, a refused command (or an unsupported tool) runs the real GATK instead, with the
//! original arguments.
//!
//! Grammar: `--name value`, `--name=value`, short names (`-I x`), repeated list arguments, booleans with or
//! without a value (`--flag`, `--flag true`), Picard's `NAME=value` for MarkDuplicates, and the `gatk` launcher's
//! own `--java-options "..."` anywhere.

use std::collections::HashMap;

struct GatkArg {
    long: &'static str,
    short: Option<&'static str>,
    typ: &'static str,
    /// None: required or no default listed
    default: Option<&'static str>,
    multi: bool,
    section: &'static str,
}

fn table(tsv: &'static str) -> Vec<GatkArg> {
    tsv.lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .map(|l| {
            let f: Vec<&'static str> = l.split('\t').collect();
            GatkArg {
                long: f[0],
                short: (f[1] != "-").then_some(f[1]),
                typ: f[2],
                default: (f[3] != "-").then_some(f[3]),
                multi: f[4] == "1",
                section: f[5],
            }
        })
        .collect()
}

/// How a GATK argument reaches the variant-rs tool.
struct Mapped {
    gatk: &'static str,
    target: &'static str,
    /// how often it may be given
    max: usize,
}

struct Tool {
    name: &'static str,
    tsv: &'static str,
    /// variant-rs argv prefix (the old lowercase subcommand; empty for HaplotypeCaller)
    sub: &'static [&'static str],
    map: &'static [Mapped],
    /// tool-specific GATK arguments with no effect on the output (besides `COMMON_IGNORED`)
    ignore: &'static [&'static str],
    /// variant-rs options: (name, takes a value)
    extras: &'static [(&'static str, bool)],
    /// whether the tool takes `--threads` (defaulted to all CPUs: outputs do not depend on it)
    threads: bool,
    picard: bool,
    /// what `-L` may be
    intervals: Iv,
}

#[derive(PartialEq)]
enum Iv {
    /// interval lists, BED files and `contig[:start-end]` strings (not VCF / GFF feature files)
    Lists,
    /// a sites VCF, or what `Lists` takes (GetPileupSummaries)
    SitesVcf,
}

const fn m(gatk: &'static str, target: &'static str, max: usize) -> Mapped {
    Mapped { gatk, target, max }
}

const MANY: usize = usize::MAX;

/// No effect on the output records: logging, progress, temporary-file and memory knobs, JDK codecs, cloud
/// buffers, index creation (variant-rs always writes the index).
const COMMON_IGNORED: &[&str] = &[
    "--verbosity",
    "--QUIET",
    "--seconds-between-progress-updates",
    "--use-jdk-deflater",
    "--use-jdk-inflater",
    "--cloud-prefetch-buffer",
    "--cloud-index-prefetch-buffer",
    "--gcs-max-retries",
    "--gcs-project-for-requester-pays",
    "--disable-bam-index-caching",
    "--create-output-variant-index",
    "--create-output-bam-index",
    "--max-variants-per-shard",
    "--VERBOSITY",
    "--USE_JDK_DEFLATER",
    "--USE_JDK_INFLATER",
    "--MAX_RECORDS_IN_RAM",
    "--CREATE_INDEX",
    "--MAX_FILE_HANDLES_FOR_READ_ENDS_MAP",
    "--MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP",
    "--SORTING_COLLECTION_SIZE_RATIO",
    "--GA4GH_CLIENT_SECRETS",
];

const TOOLS: &[Tool] = &[
    Tool {
        name: "HaplotypeCaller",
        tsv: include_str!("../gatkargs/HaplotypeCaller.tsv"),
        sub: &[],
        map: &[
            m("--input", "-I", MANY),
            m("--reference", "-R", 1),
            m("--intervals", "-L", MANY),
            m("--interval-padding", "--interval-padding", 1),
            m("--exclude-intervals", "-XL", MANY),
            m("--interval-exclusion-padding", "--interval-exclusion-padding", 1),
            m("--output", "-O", 1),
            m("--emit-ref-confidence", "-ERC", 1),
            m("--dbsnp", "--dbsnp", 1),
            m("--gvcf-gq-bands", "-GQB", MANY),
            m("--contamination-fraction-to-filter", "-contamination", 1),
            m("--annotation-group", "-G", MANY),
            m("--max-reads-per-alignment-start", "--max-reads-per-alignment-start", 1),
            m("--assembly-region-out", "--assembly-region-out", 1),
            m("--tmp-dir", "--tmp-dir", 1),
        ],
        // Smith-Waterman implementations give identical results (GKL reproduces the Java aligner); the
        // native PairHMM thread count does not change results.
        ignore: &["--native-pair-hmm-threads", "--smith-waterman"],
        extras: &[
            ("--threads", true),
            ("--pairhmm", true),
            ("--contig-workers", true),
            ("--max-workers-per-shard", true),
            ("--shard", true),
            ("--shard-output-dir", true),
            ("--compression-level", true),
            ("--bqsr-recal-file", true),
            ("-bqsr", true),
            ("--gvcf-raw", true),
            ("--vcf-header", true),
            ("--dump-activity", true),
        ],
        threads: true,
        picard: false,
        intervals: Iv::Lists,
    },
    Tool {
        name: "BaseRecalibrator",
        tsv: include_str!("../gatkargs/BaseRecalibrator.tsv"),
        sub: &["bqsr"],
        map: &[
            m("--input", "-I", 1),
            m("--reference", "-R", 1),
            m("--intervals", "-L", MANY),
            m("--output", "-O", 1),
            m("--known-sites", "--known-sites", MANY),
            m("--use-original-qualities", "--use-original-qualities", 1),
        ],
        ignore: &["--tmp-dir"],
        extras: &[("--threads", true), ("--chunk", true)],
        threads: true,
        picard: false,
        intervals: Iv::Lists,
    },
    Tool {
        name: "ApplyBQSR",
        tsv: include_str!("../gatkargs/ApplyBQSR.tsv"),
        sub: &["applybqsr"],
        map: &[
            m("--input", "-I", 1),
            m("--reference", "-R", 1),
            m("--intervals", "-L", MANY),
            m("--output", "-O", 1),
            m("--bqsr-recal-file", "--bqsr-recal-file", 1),
            m("--use-original-qualities", "--use-original-qualities", 1),
            m("--static-quantized-quals", "--static-quantized-quals", MANY),
            m("--round-down-quantized", "--round-down-quantized", 1),
            m("--create-output-bam-md5", "--create-output-bam-md5", 1),
        ],
        ignore: &["--tmp-dir"],
        extras: &[("--threads", true), ("--chunk", true), ("--compression-level", true)],
        threads: true,
        picard: false,
        intervals: Iv::Lists,
    },
    Tool {
        name: "MarkDuplicates",
        tsv: include_str!("../gatkargs/MarkDuplicates.tsv"),
        sub: &["markdup"],
        map: &[
            m("--INPUT", "-I", MANY),
            m("--OUTPUT", "-O", 1),
            m("--METRICS_FILE", "-M", 1),
            m("--REFERENCE_SEQUENCE", "--reference", 1),
            m("--TMP_DIR", "--tmp-dir", MANY),
            m("--COMPRESSION_LEVEL", "--compression-level", 1),
            m("--ASSUME_SORT_ORDER", "--ASSUME_SORT_ORDER", 1),
            m("--OPTICAL_DUPLICATE_PIXEL_DISTANCE", "--OPTICAL_DUPLICATE_PIXEL_DISTANCE", 1),
            m("--CLEAR_DT", "--CLEAR_DT", 1),
            m("--ADD_PG_TAG_TO_READS", "--ADD_PG_TAG_TO_READS", 1),
            m("--CREATE_INDEX", "--CREATE_INDEX", 1),
            m("--CREATE_MD5_FILE", "--CREATE_MD5_FILE", 1),
        ],
        // VALIDATION_STRINGENCY only decides whether Picard stops on malformed records (variant-rs does not
        // validate); on inputs Picard accepts, the output does not depend on it
        ignore: &["--VALIDATION_STRINGENCY"],
        extras: &[("--threads", true), ("--chunk", true), ("--bucket-bp", true)],
        threads: true,
        picard: true,
        intervals: Iv::Lists,
    },
    Tool {
        name: "GatherBQSRReports",
        tsv: include_str!("../gatkargs/GatherBQSRReports.tsv"),
        sub: &["gatherbqsrreports"],
        map: &[m("--input", "-I", MANY), m("--output", "-O", 1)],
        ignore: &["--tmp-dir"],
        extras: &[],
        threads: false,
        picard: false,
        intervals: Iv::Lists,
    },
    Tool {
        name: "MergeMutectStats",
        tsv: include_str!("../gatkargs/MergeMutectStats.tsv"),
        sub: &["mergemutectstats"],
        map: &[m("--stats", "--stats", MANY), m("--output", "-O", 1)],
        ignore: &["--tmp-dir"],
        extras: &[],
        threads: false,
        picard: false,
        intervals: Iv::Lists,
    },
    Tool {
        name: "GatherPileupSummaries",
        tsv: include_str!("../gatkargs/GatherPileupSummaries.tsv"),
        sub: &["gatherpileupsummaries"],
        map: &[m("--I", "-I", MANY), m("--O", "-O", 1), m("--sequence-dictionary", "--sequence-dictionary", 1)],
        ignore: &["--tmp-dir"],
        extras: &[],
        threads: false,
        picard: false,
        intervals: Iv::Lists,
    },
    Tool {
        name: "GatherBamFiles",
        tsv: include_str!("../gatkargs/GatherBamFiles.tsv"),
        sub: &["gatherbamfiles"],
        map: &[
            m("--INPUT", "-I", MANY),
            m("--OUTPUT", "-O", 1),
            m("--CREATE_INDEX", "--CREATE_INDEX", 1),
            m("--CREATE_MD5_FILE", "--CREATE_MD5_FILE", 1),
            m("--COMPRESSION_LEVEL", "--compression-level", 1),
        ],
        // the inputs' blocks are copied: no record is read, so validation and the temporary directory play no part
        ignore: &["--VALIDATION_STRINGENCY", "--TMP_DIR", "--MAX_RECORDS_IN_RAM"],
        extras: &[],
        threads: false,
        picard: true,
        intervals: Iv::Lists,
    },
    Tool {
        name: "MergeVcfs",
        tsv: include_str!("../gatkargs/MergeVcfs.tsv"),
        sub: &["mergevcfs"],
        map: &[
            m("--INPUT", "-I", MANY),
            m("--OUTPUT", "-O", 1),
            m("--CREATE_INDEX", "--CREATE_INDEX", 1),
            m("--COMPRESSION_LEVEL", "--compression-level", 1),
            m("--SEQUENCE_DICTIONARY", "--SEQUENCE_DICTIONARY", 1),
        ],
        // records are streamed (no sorting collection), no reads are validated
        ignore: &["--VALIDATION_STRINGENCY", "--TMP_DIR", "--MAX_RECORDS_IN_RAM"],
        extras: &[],
        threads: false,
        picard: true,
        intervals: Iv::Lists,
    },
    Tool {
        name: "GenotypeGVCFs",
        tsv: include_str!("../gatkargs/GenotypeGVCFs.tsv"),
        sub: &["genotypegvcfs"],
        map: &[
            m("--variant", "-V", 1),
            m("--reference", "-R", 1),
            m("--output", "-O", 1),
            m("--dbsnp", "-D", 1),
            m("--annotation-group", "-G", MANY),
            m("--intervals", "-L", MANY),
            m("--only-output-calls-starting-in-intervals", "--only-output-calls-starting-in-intervals", 1),
            m("--merge-input-intervals", "--merge-input-intervals", 1),
        ],
        ignore: &["--tmp-dir"],
        extras: &[],
        threads: false,
        picard: false,
        intervals: Iv::Lists,
    },
    Tool {
        name: "ReblockGVCF",
        tsv: include_str!("../gatkargs/ReblockGVCF.tsv"),
        sub: &["reblockgvcf"],
        map: &[
            m("--variant", "-V", 1),
            m("--reference", "-R", 1),
            m("--output", "-O", 1),
            m("--do-qual-score-approximation", "--do-qual-score-approximation", 1),
            m("--floor-blocks", "--floor-blocks", 1),
            m("--gvcf-gq-bands", "-GQB", MANY),
        ],
        ignore: &["--tmp-dir"],
        extras: &[],
        threads: false,
        picard: false,
        intervals: Iv::Lists,
    },
    Tool {
        name: "SortSam",
        tsv: include_str!("../gatkargs/SortSam.tsv"),
        sub: &["sortsam"],
        map: &[
            m("--INPUT", "-I", 1),
            m("--OUTPUT", "-O", 1),
            m("--SORT_ORDER", "--SORT_ORDER", 1),
            m("--CREATE_INDEX", "--CREATE_INDEX", 1),
            m("--CREATE_MD5_FILE", "--CREATE_MD5_FILE", 1),
            m("--TMP_DIR", "--tmp-dir", 1),
            m("--COMPRESSION_LEVEL", "--compression-level", 1),
        ],
        // MAX_RECORDS_IN_RAM: the spill size of a stable sort, which cannot change the order; VALIDATION_STRINGENCY:
        // as MarkDuplicates
        ignore: &["--MAX_RECORDS_IN_RAM", "--VALIDATION_STRINGENCY"],
        extras: &[("--threads", true), ("--sort-buffer-mb", true)],
        threads: true,
        picard: true,
        intervals: Iv::Lists,
    },
    // The Mutect2 family. Mutect2 and GetPileupSummaries use all cores (env VARIANT_RS_THREADS overrides); outputs do
    // not depend on it.
    Tool {
        name: "Mutect2",
        tsv: include_str!("../gatkargs/Mutect2.tsv"),
        sub: &["mutect2"],
        map: &[
            m("--input", "-I", MANY),
            m("--reference", "-R", 1),
            m("--intervals", "-L", MANY),
            m("--output", "-O", 1),
            m("--normal-sample", "--normal-sample", MANY),
            m("--germline-resource", "--germline-resource", 1),
            m("--panel-of-normals", "--panel-of-normals", 1),
            m("--f1r2-tar-gz", "--f1r2-tar-gz", 1),
            m("--max-reads-per-alignment-start", "--max-reads-per-alignment-start", 1),
            m("--assembly-region-out", "--assembly-region-out", 1),
            m("--genotype-germline-sites", "--genotype-germline-sites", 1),
            m("--genotype-pon-sites", "--genotype-pon-sites", 1),
            m("--af-of-alleles-not-in-resource", "--af-of-alleles-not-in-resource", 1),
            m("--mitochondria-mode", "--mitochondria-mode", 1),
            m("--annotation", "--annotation", MANY),
            m("--read-filter", "--read-filter", MANY),
            m("--max-mnp-distance", "--max-mnp-distance", 1),
            m("--bam-output", "--bam-output", 1),
        ],
        // --tumor-sample is deprecated and unused since GATK 4.1 (every non-normal sample is a tumor); see
        // HaplotypeCaller for the Smith-Waterman and PairHMM thread arguments.
        ignore: &["--tmp-dir", "--native-pair-hmm-threads", "--smith-waterman", "--tumor-sample"],
        extras: &[("--pairhmm", true), ("--dump-activity", true), ("--dump-s2", true)],
        threads: false,
        picard: false,
        intervals: Iv::Lists,
    },
    Tool {
        name: "FilterMutectCalls",
        tsv: include_str!("../gatkargs/FilterMutectCalls.tsv"),
        sub: &["filtermutectcalls"],
        map: &[
            m("--variant", "--variant", 1),
            m("--reference", "--reference", 1),
            m("--output", "--output", 1),
            m("--stats", "--stats", 1),
            m("--filtering-stats", "--filtering-stats", 1),
            m("--contamination-table", "--contamination-table", MANY),
            m("--contamination-estimate", "--contamination-estimate", 1),
            m("--tumor-segmentation", "--tumor-segmentation", MANY),
            m("--orientation-bias-artifact-priors", "--orientation-bias-artifact-priors", MANY),
            m("--mitochondria-mode", "--mitochondria-mode", 1),
            m("--min-allele-fraction", "--min-allele-fraction", 1),
            m("--max-alt-allele-count", "--max-alt-allele-count", 1),
            m("--f-score-beta", "--f-score-beta", 1),
        ],
        ignore: &["--tmp-dir"],
        extras: &[],
        threads: false,
        picard: false,
        intervals: Iv::Lists,
    },
    Tool {
        name: "GetPileupSummaries",
        tsv: include_str!("../gatkargs/GetPileupSummaries.tsv"),
        sub: &["getpileupsummaries"],
        map: &[
            m("--input", "--input", 1),
            m("--variant", "--variant", 1),
            m("--intervals", "--intervals", 1),
            m("--output", "--output", 1),
            m("--reference", "--reference", 1),
            m("--minimum-population-allele-frequency", "--minimum-population-allele-frequency", 1),
            m("--maximum-population-allele-frequency", "--maximum-population-allele-frequency", 1),
        ],
        ignore: &["--tmp-dir"],
        extras: &[],
        threads: false,
        picard: false,
        intervals: Iv::SitesVcf,
    },
    Tool {
        name: "CalculateContamination",
        tsv: include_str!("../gatkargs/CalculateContamination.tsv"),
        sub: &["calculatecontamination"],
        map: &[
            m("--input", "--input", 1),
            m("--output", "--output", 1),
            m("--matched-normal", "--matched-normal", 1),
            m("--tumor-segmentation", "--tumor-segmentation", 1),
            m("--low-coverage-ratio-threshold", "--low-coverage-ratio-threshold", 1),
            m("--high-coverage-ratio-threshold", "--high-coverage-ratio-threshold", 1),
        ],
        ignore: &["--tmp-dir"],
        extras: &[],
        threads: false,
        picard: false,
        intervals: Iv::Lists,
    },
    Tool {
        name: "LearnReadOrientationModel",
        tsv: include_str!("../gatkargs/LearnReadOrientationModel.tsv"),
        sub: &["learnreadorientationmodel"],
        map: &[m("--input", "--input", MANY), m("--output", "--output", 1)],
        ignore: &["--tmp-dir"],
        extras: &[],
        threads: false,
        picard: false,
        intervals: Iv::Lists,
    },
];

/// Whether variant-rs reads this `-L` value: interval files and `contig[:start-end]` strings, except feature files
/// (VCF, GFF), which only GetPileupSummaries reads (a sites VCF, and nothing else).
fn interval_value_ok(kind: &Iv, v: &str) -> bool {
    let vcf = v.ends_with(".vcf") || v.ends_with(".vcf.gz");
    match kind {
        Iv::Lists => !vcf && ![".bcf", ".gff", ".gtf"].iter().any(|e| v.ends_with(e)),
        Iv::SitesVcf => vcf || interval_value_ok(&Iv::Lists, v),
    }
}

/// The supported GATK tool names.
pub fn tool_names() -> Vec<&'static str> {
    TOOLS.iter().map(|t| t.name).collect()
}

/// Whether `name` looks like a GATK tool name (GATK and Picard tools are CamelCase; variant-rs' own subcommands are
/// lowercase).
pub fn is_tool_name(name: &str) -> bool {
    name.chars().next().is_some_and(|c| c.is_ascii_uppercase()) && name.chars().all(|c| c.is_ascii_alphanumeric())
}

fn same_value(typ: &str, a: &str, b: &str) -> bool {
    match typ {
        "Integer" | "Long" | "Double" | "Float" | "Byte" | "int" | "double" => match (a.parse::<f64>(), b.parse::<f64>()) {
            (Ok(x), Ok(y)) => x == y,
            _ => a == b,
        },
        "Boolean" | "boolean" => a.eq_ignore_ascii_case(b),
        _ => a == b,
    }
}

/// A default list as printed by GATK: `null`, `[]` or `[a, b]`; scalars as a one-element list.
fn default_list(d: &str) -> Vec<String> {
    match d {
        "null" | "[]" => Vec::new(),
        _ if d.starts_with('[') && d.ends_with(']') => d[1..d.len() - 1].split(", ").map(|s| s.to_string()).collect(),
        _ => vec![d.to_string()],
    }
}

/// Splits the arguments into (name, optional value) pairs following GATK's grammar.
fn tokenize(args: &[String], picard: bool, is_bool: &dyn Fn(&str) -> bool) -> Result<Vec<(String, Option<String>)>, String> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        if a == "--java-options" {
            // the gatk launcher's JVM options
            i += 2;
            continue;
        }
        if let Some(rest) = a.strip_prefix("--java-options=") {
            let _ = rest;
            i += 1;
            continue;
        }
        if !a.starts_with('-') {
            // Picard's legacy NAME=value
            if picard {
                if let Some((k, v)) = a.split_once('=') {
                    out.push((format!("--{k}"), Some(v.to_string())));
                    i += 1;
                    continue;
                }
            }
            return Err(format!("unexpected positional argument `{a}`"));
        }
        if let Some((k, v)) = a.split_once('=').filter(|(k, _)| k.starts_with("--")) {
            out.push((k.to_string(), Some(v.to_string())));
            i += 1;
            continue;
        }
        let next = args.get(i + 1);
        if is_bool(a) {
            match next {
                Some(v) if v.eq_ignore_ascii_case("true") || v.eq_ignore_ascii_case("false") => {
                    out.push((a.clone(), Some(v.clone())));
                    i += 2;
                }
                _ => {
                    out.push((a.clone(), Some("true".into())));
                    i += 1;
                }
            }
            continue;
        }
        match next {
            Some(v) => {
                out.push((a.clone(), Some(v.clone())));
                i += 2;
            }
            None => return Err(format!("{a} needs a value")),
        }
    }
    Ok(out)
}

/// The result of translating a GATK command line.
#[derive(Debug, PartialEq)]
pub enum Translation {
    /// run variant-rs with these arguments (the tool's subcommand first)
    Run(Vec<String>),
    /// GATK's own `--help` / `--version` style request, already answered (text to print)
    Info(String),
}

/// Translates `gatk <tool> args...` into variant-rs arguments, or explains why it cannot.
pub fn translate(tool_name: &str, args: &[String], cpus: usize) -> Result<Translation, String> {
    let tool = TOOLS
        .iter()
        .find(|t| t.name == tool_name)
        .ok_or_else(|| format!("{tool_name} is not supported by variant-rs (supported: {})", tool_names().join(", ")))?;
    let gatk = table(tool.tsv);
    let by_name: HashMap<&str, &GatkArg> = gatk.iter().flat_map(|a| std::iter::once((a.long, a)).chain(a.short.map(|s| (s, a)))).collect();
    // GATK's parser takes either prefix for either name (`--ob-priors`, `-normal-sample`); Picard's legacy
    // NAME=value takes short names too (I=, O=, M=)
    let lookup = |n: &str| -> Option<&GatkArg> {
        by_name.get(n).copied().or_else(|| {
            let bare = n.trim_start_matches('-');
            [format!("-{bare}"), format!("--{bare}")].iter().find_map(|x| by_name.get(x.as_str()).copied())
        })
    };
    let extra = |n: &str| tool.extras.iter().find(|(e, _)| *e == n);
    let is_bool = |n: &str| match lookup(n) {
        Some(a) => a.typ == "Boolean" || a.typ == "boolean",
        None => extra(n).is_some_and(|(_, takes)| !takes),
    };
    let toks = tokenize(args, tool.picard, &is_bool)?;
    let mut out: Vec<String> = tool.sub.iter().map(|s| s.to_string()).collect();
    let mut counts: HashMap<&str, usize> = HashMap::new();
    // values of default-only list arguments, compared with the default list at the end
    let mut lists: HashMap<&str, Vec<String>> = HashMap::new();
    let mut problems: Vec<String> = Vec::new();
    let mut threads_given = false;
    for (name, val) in &toks {
        let val = val.clone().unwrap_or_default();
        if let Some((e, takes)) = extra(name) {
            threads_given |= *e == "--threads";
            out.push(e.to_string());
            if *takes {
                out.push(val);
            }
            continue;
        }
        let Some(a) = lookup(name) else {
            problems.push(format!("{name}: not a {} argument", tool.name));
            continue;
        };
        if a.long == "--help" || a.long == "--version" {
            if same_value("Boolean", &val, "true") {
                return Ok(Translation::Info(help_text(tool, &gatk)));
            }
            continue;
        }
        if let Some(mp) = tool.map.iter().find(|mp| mp.gatk == a.long) {
            let c = counts.entry(a.long).or_insert(0);
            *c += 1;
            if *c > mp.max {
                problems.push(format!("{}: given {} times; variant-rs takes it {} time(s)", a.long, c, mp.max));
                continue;
            }
            if (mp.gatk == "--intervals" || mp.gatk == "--exclude-intervals") && !interval_value_ok(&tool.intervals, &val) {
                let ok = match tool.intervals {
                    Iv::Lists => "use .interval_list, .bed, .intervals/.list or contig:start-end",
                    Iv::SitesVcf => "use a sites VCF, .interval_list, .bed, .intervals/.list or contig:start-end",
                };
                problems.push(format!("{} {val}: interval files in this format are not supported ({ok})", a.long));
                continue;
            }
            if tool.name == "Mutect2" && mp.gatk == "--annotation" && !["StrandBiasBySample", "OriginalAlignment"].contains(&val.as_str()) {
                problems.push(format!("{} {val}: only StrandBiasBySample (a default) and OriginalAlignment are supported", a.long));
                continue;
            }
            if mp.gatk == "--read-filter" && !["MateOnSameContigOrNoMappedMateReadFilter", "MateUnmappedAndUnmappedReadFilter"].contains(&val.as_str()) {
                problems.push(format!("{} {val}: only MateOnSameContigOrNoMappedMateReadFilter and MateUnmappedAndUnmappedReadFilter are supported", a.long));
                continue;
            }
            if mp.gatk == "--SORT_ORDER" && val != "coordinate" {
                problems.push(format!("{} {val}: only coordinate is supported", a.long));
                continue;
            }
            if tool.name == "GenotypeGVCFs" && mp.gatk == "--annotation-group" && !["StandardAnnotation", "AS_StandardAnnotation"].contains(&val.as_str()) {
                problems.push(format!("{} {val}: only StandardAnnotation (GenotypeGVCFs' default) and AS_StandardAnnotation are supported", a.long));
                continue;
            }
            if mp.gatk == "--annotation-group" && !["StandardAnnotation", "StandardHCAnnotation", "AS_StandardAnnotation"].contains(&val.as_str()) {
                problems.push(format!("{} {val}: only StandardAnnotation, StandardHCAnnotation (HaplotypeCaller's defaults) and AS_StandardAnnotation are supported", a.long));
                continue;
            }
            // a GATK boolean becomes the variant-rs flag when true, nothing when false (Picard booleans keep their
            // value: markdup takes `--CLEAR_DT false`)
            if (a.typ == "Boolean" || a.typ == "boolean") && !tool.picard {
                if same_value("Boolean", &val, "true") {
                    out.push(mp.target.to_string());
                }
                continue;
            }
            out.push(mp.target.to_string());
            out.push(val);
            continue;
        }
        if COMMON_IGNORED.contains(&a.long) || tool.ignore.contains(&a.long) {
            continue;
        }
        // default-only
        match a.default {
            Some(_) if a.multi => lists.entry(a.long).or_default().push(val),
            Some(d) if same_value(a.typ, &val, d) => {}
            Some(d) => problems.push(format!("{} {val}: only GATK's default ({d}) is supported; other values change the output", a.long)),
            None => problems.push(format!("{} ({} argument): not supported", a.long, a.section)),
        }
    }
    // --mitochondria-mode overrides these unless they are given; given (even at GATK's default), GATK keeps them
    let mito = toks.iter().any(|(n, v)| lookup(n).is_some_and(|a| a.long == "--mitochondria-mode") && same_value("Boolean", v.as_deref().unwrap_or("true"), "true"));
    if mito {
        for (n, _) in &toks {
            if let Some(a) = lookup(n) {
                if ["--tumor-lod-to-emit", "--initial-tumor-lod", "--recover-all-dangling-branches", "--pruning-lod-threshold"].contains(&a.long) {
                    problems.push(format!("{} with --mitochondria-mode: not supported (the mode's value would be replaced)", a.long));
                }
            }
        }
    }
    for (long, vals) in &lists {
        let a = by_name[long];
        let d = default_list(a.default.unwrap());
        let same = vals.len() == d.len() && vals.iter().zip(&d).all(|(x, y)| same_value(a.typ, x, y));
        if !same {
            problems.push(format!("{long} {}: only GATK's default ({}) is supported; other values change the output", vals.join(" "), a.default.unwrap()));
        }
    }
    for a in gatk.iter().filter(|a| a.section == "required") {
        if let Some(mp) = tool.map.iter().find(|mp| mp.gatk == a.long) {
            if !counts.contains_key(mp.gatk) {
                problems.push(format!("{}{}: required", a.long, a.short.map(|s| format!(" ({s})")).unwrap_or_default()));
            }
        }
    }
    if !problems.is_empty() {
        return Err(format!("{} arguments variant-rs cannot run exactly as GATK would:\n  {}", tool.name, problems.join("\n  ")));
    }
    if tool.threads && !threads_given {
        out.push("--threads".into());
        out.push(cpus.to_string());
    }
    Ok(Translation::Run(out))
}

fn help_text(tool: &Tool, gatk: &[GatkArg]) -> String {
    let mut s = format!("variant-rs {} (GATK 4.7 compatible; outputs byte-identical to GATK)\n\nArguments mapped to variant-rs:\n", tool.name);
    for mp in tool.map {
        let a = gatk.iter().find(|a| a.long == mp.gatk);
        let short = a.and_then(|a| a.short).map(|s| format!(",{s}")).unwrap_or_default();
        s.push_str(&format!("  {}{short}\n", mp.gatk));
    }
    s.push_str("\nIgnored (no effect on the output):\n  ");
    // GATK's arguments of this tool not mapped above, each once
    let mut ign: Vec<&str> = Vec::new();
    for n in COMMON_IGNORED.iter().chain(tool.ignore).copied() {
        if gatk.iter().any(|a| a.long == n) && !tool.map.iter().any(|mp| mp.gatk == n) && !ign.contains(&n) {
            ign.push(n);
        }
    }
    s.push_str(&ign.join(" "));
    s.push_str("\n\nEvery other GATK argument is accepted only at its default value.\n\nvariant-rs options:\n  ");
    let ex: Vec<&str> = tool.extras.iter().map(|(e, _)| *e).collect();
    s.push_str(&ex.join(" "));
    s.push_str("\n\nSet GATK_FALLBACK=/path/to/gatk to run GATK for commands variant-rs refuses.\n");
    s
}

/// `variant-rs <Tool> args...`: translates, or falls back to GATK (`GATK_FALLBACK`), or exits with the reason.
/// Returns the variant-rs arguments to run.
pub fn dispatch(argv: &[String]) -> Vec<String> {
    let tool = &argv[1];
    let args = &argv[2..];
    let cpus = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    // VARIANT_RS_DRY_RUN=1: only say what would happen (argument level; refusals that depend on the inputs' contents
    // are found when running)
    if std::env::var_os("VARIANT_RS_DRY_RUN").is_some_and(|v| !v.is_empty() && v != "0") {
        match translate(tool, args, cpus) {
            Ok(Translation::Run(v)) => println!("RUN variant-rs {}", v.join(" ")),
            Ok(Translation::Info(_)) => println!("INFO {tool}"),
            Err(why) => println!("REFUSE {why}"),
        }
        std::process::exit(0);
    }
    match translate(tool, args, cpus) {
        Ok(Translation::Run(v)) => v,
        Ok(Translation::Info(text)) => {
            print!("{text}");
            std::process::exit(0);
        }
        Err(why) => refuse_argv(&why),
    }
}

/// The process arguments without the `gatk` launcher's `--java-options` in front of the tool name.
pub fn strip_launcher_options(mut argv: Vec<String>) -> Vec<String> {
    while argv.len() > 2 && (argv[1] == "--java-options" || argv[1].starts_with("--java-options=")) {
        let n = if argv[1] == "--java-options" { 2 } else { 1 };
        argv.drain(1..1 + n);
    }
    argv
}

/// The GATK release whose outputs variant-rs reproduces (what `gatk --version` reports).
pub const GATK_VERSION: &str = "4.7.0.0";

/// Runs GATK (`GATK_FALLBACK`) with the process' original arguments.
pub fn run_gatk_instead() -> ! {
    refuse_argv("a request for GATK itself")
}

/// A refusal found after translation (from the inputs' contents): with a GATK command line and `GATK_FALLBACK`,
/// runs GATK with the original arguments; otherwise exits with the reason.
pub fn refuse(why: &str) -> ! {
    let argv = strip_launcher_options(std::env::args().collect());
    if argv.len() > 1 && is_tool_name(&argv[1]) {
        refuse_argv(why)
    }
    eprintln!("variant-rs: {why}");
    std::process::exit(2);
}

/// Runs GATK (`GATK_FALLBACK`) with the process' original arguments (launcher options included), or exits.
fn refuse_argv(why: &str) -> ! {
    let argv: Vec<String> = std::env::args().collect();
    match std::env::var_os("GATK_FALLBACK") {
        Some(gatk) => {
            eprintln!("variant-rs: {}\nvariant-rs: running GATK instead ({})", why, gatk.to_string_lossy());
            use std::os::unix::process::CommandExt;
            let err = std::process::Command::new(&gatk).args(&argv[1..]).exec();
            eprintln!("variant-rs: could not run {}: {err}", gatk.to_string_lossy());
            std::process::exit(127);
        }
        None => {
            eprintln!("variant-rs: {why}\n(set GATK_FALLBACK=/path/to/gatk to run GATK for such commands)");
            std::process::exit(2);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn s(v: &[&str]) -> Vec<String> {
        v.iter().map(|x| x.to_string()).collect()
    }
    fn run(tool: &str, a: &[&str]) -> Result<Vec<String>, String> {
        match translate(tool, &s(a), 8)? {
            Translation::Run(v) => Ok(v),
            Translation::Info(_) => Err("info".into()),
        }
    }

    #[test]
    fn haplotypecaller_forms() {
        let v = run("HaplotypeCaller", &["-R", "r.fa", "--input=x.bam", "-L", "chr1", "--intervals", "t.interval_list", "-O", "o.g.vcf.gz", "-ERC", "GVCF"]).unwrap();
        assert_eq!(v, s(&["-R", "r.fa", "-I", "x.bam", "-L", "chr1", "-L", "t.interval_list", "-O", "o.g.vcf.gz", "-ERC", "GVCF", "--threads", "8"]));
    }

    #[test]
    fn ignored_defaults_and_launcher_options() {
        let v = run(
            "HaplotypeCaller",
            &["--java-options", "-Xmx4g", "-R", "r.fa", "-I", "x.bam", "-O", "o.vcf.gz", "--native-pair-hmm-threads", "2", "-OVI", "--tmp-dir", "/t",
              "-ploidy", "2", "--min-base-quality-score", "10", "--standard-min-confidence-threshold-for-calling", "30.0", "--kmer-size", "10", "--kmer-size", "25",
              "--sites-only-vcf-output", "false", "--threads", "4"],
        )
        .unwrap();
        assert_eq!(v, s(&["-R", "r.fa", "-I", "x.bam", "-O", "o.vcf.gz", "--tmp-dir", "/t", "--threads", "4"]));
    }

    #[test]
    fn refusals() {
        let e = run("HaplotypeCaller", &["-R", "r.fa", "-I", "x.bam", "-O", "o.vcf.gz", "-ploidy", "1", "-G", "AlleleFraction", "--kmer-size", "21", "-L", "s.vcf.gz", "--bogus", "1"]).unwrap_err();
        for w in ["--sample-ploidy 1", "--annotation-group AlleleFraction", "--kmer-size 21", "s.vcf.gz", "--bogus"] {
            assert!(e.contains(w), "{w} missing in: {e}");
        }
        let v = run("HaplotypeCaller", &["-R", "r.fa", "-I", "a.bam", "-I", "b.bam", "-O", "o.vcf"]).unwrap();
        assert_eq!(v[..8], s(&["-R", "r.fa", "-I", "a.bam", "-I", "b.bam", "-O", "o.vcf"]));
        let e = run("HaplotypeCaller", &["-R", "r.fa", "-R", "s.fa", "-I", "a.bam", "-O", "o.vcf"]).unwrap_err();
        assert!(e.contains("given 2 times"));
        let e = run("HaplotypeCaller", &["-R", "r.fa", "-O", "o.vcf"]).unwrap_err();
        assert!(e.contains("--input (-I): required"));
        assert!(translate("GenotypeGVCFs", &[], 1).is_err());
    }

    #[test]
    fn booleans() {
        let v = run("HaplotypeCaller", &["-R", "r.fa", "-I", "x.bam", "-O", "o.vcf", "--sites-only-vcf-output", "--threads", "2"]);
        assert!(v.unwrap_err().contains("--sites-only-vcf-output true"));
        let v = run("HaplotypeCaller", &["-R", "r.fa", "-I", "x.bam", "-O", "o.vcf", "--create-output-variant-index", "false", "--threads", "2"]).unwrap();
        assert_eq!(v, s(&["-R", "r.fa", "-I", "x.bam", "-O", "o.vcf", "--threads", "2"]));
    }

    #[test]
    fn markduplicates_picard_syntax() {
        let v = run("MarkDuplicates", &["I=in.bam", "O=out.bam", "M=m.txt", "CREATE_INDEX=true", "--TMP_DIR", "/t"]).unwrap();
        assert_eq!(v, s(&["markdup", "-I", "in.bam", "-O", "out.bam", "-M", "m.txt", "--CREATE_INDEX", "true", "--tmp-dir", "/t", "--threads", "8"]));
        let v = run("MarkDuplicates", &["-I", "in.bam", "-O", "out.bam", "-M", "m.txt", "--CREATE_INDEX", "true", "--REMOVE_DUPLICATES", "true"]);
        assert!(v.unwrap_err().contains("REMOVE_DUPLICATES true"));
        // WARP's command line (tasks/wdl/BamProcessing.wdl)
        let v = run("MarkDuplicates", &["INPUT=a.bam", "INPUT=b.bam", "OUTPUT=o.bam", "METRICS_FILE=m.txt", "VALIDATION_STRINGENCY=SILENT",
            "OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500", "ASSUME_SORT_ORDER=queryname", "CLEAR_DT=false", "ADD_PG_TAG_TO_READS=false"]).unwrap();
        // BamProcessing.wdl's GatherBqsrReports
        let g = run("GatherBQSRReports", &["-I", "a.table", "-I", "b.table", "-O", "g.table"]).unwrap();
        assert_eq!(g, s(&["gatherbqsrreports", "-I", "a.table", "-I", "b.table", "-O", "g.table"]));
        // BamProcessing.wdl's GatherSortedBamFiles
        let gb = run("GatherBamFiles", &["INPUT=a.bam", "INPUT=b.bam", "OUTPUT=o.bam", "CREATE_INDEX=true", "CREATE_MD5_FILE=true"]).unwrap();
        assert_eq!(gb, s(&["gatherbamfiles", "-I", "a.bam", "-I", "b.bam", "-O", "o.bam", "--CREATE_INDEX", "true", "--CREATE_MD5_FILE", "true"]));
        // GermlineVariantDiscovery.wdl's MergeVCFs
        let mv = run("MergeVcfs", &["INPUT=a.g.vcf.gz", "INPUT=b.g.vcf.gz", "OUTPUT=o.g.vcf.gz"]).unwrap();
        assert_eq!(mv, s(&["mergevcfs", "-I", "a.g.vcf.gz", "-I", "b.g.vcf.gz", "-O", "o.g.vcf.gz"]));
        // GermlineVariantDiscovery.wdl's Reblock
        let rb = run("ReblockGVCF", &["-R", "r.fa", "-V", "x.g.vcf.gz", "-do-qual-approx", "--floor-blocks", "-GQB", "20", "-GQB", "30", "-GQB", "40", "-O", "x.rb.g.vcf.gz"]).unwrap();
        assert_eq!(rb, s(&["reblockgvcf", "-R", "r.fa", "-V", "x.g.vcf.gz", "--do-qual-score-approximation", "--floor-blocks", "-GQB", "20", "-GQB", "30", "-GQB", "40", "-O", "x.rb.g.vcf.gz"]));
        assert!(run("ReblockGVCF", &["-R", "r.fa", "-V", "x.g.vcf.gz", "--drop-low-quals", "-O", "y.g.vcf.gz"]).is_err());
        // JointGenotypingTasks.wdl's GenotypeGVCFs (with a GVCF file in place of the GenomicsDB workspace)
        let gg = run("GenotypeGVCFs", &["-R", "r.fa", "-O", "o.vcf.gz", "-D", "d.vcf.gz", "-G", "StandardAnnotation", "-G", "AS_StandardAnnotation",
            "--only-output-calls-starting-in-intervals", "-V", "x.g.vcf.gz", "-L", "chr8", "--merge-input-intervals"]).unwrap();
        assert_eq!(gg, s(&["genotypegvcfs", "-R", "r.fa", "-O", "o.vcf.gz", "-D", "d.vcf.gz", "-G", "StandardAnnotation", "-G", "AS_StandardAnnotation",
            "--only-output-calls-starting-in-intervals", "-V", "x.g.vcf.gz", "-L", "chr8", "--merge-input-intervals"]));
        // BamProcessing.wdl's SortSam
        let ss = run("SortSam", &["INPUT=x.bam", "OUTPUT=y.bam", "SORT_ORDER=coordinate", "CREATE_INDEX=true", "CREATE_MD5_FILE=true", "MAX_RECORDS_IN_RAM=300000"]).unwrap();
        assert_eq!(ss, s(&["sortsam", "-I", "x.bam", "-O", "y.bam", "--SORT_ORDER", "coordinate", "--CREATE_INDEX", "true", "--CREATE_MD5_FILE", "true", "--threads", "8"]));
        assert!(run("SortSam", &["I=x.bam", "O=y.bam", "SO=queryname"]).unwrap_err().contains("only coordinate"));
        assert_eq!(v, s(&["markdup", "-I", "a.bam", "-I", "b.bam", "-O", "o.bam", "-M", "m.txt", "--OPTICAL_DUPLICATE_PIXEL_DISTANCE", "2500",
            "--ASSUME_SORT_ORDER", "queryname", "--CLEAR_DT", "false", "--ADD_PG_TAG_TO_READS", "false", "--threads", "8"]));
    }

    #[test]
    fn bqsr() {
        let v = run("BaseRecalibrator", &["-R", "r.fa", "-I", "x.bam", "--known-sites", "a.vcf.gz", "--known-sites", "b.vcf.gz", "-O", "t.table"]).unwrap();
        assert_eq!(v, s(&["bqsr", "-R", "r.fa", "-I", "x.bam", "--known-sites", "a.vcf.gz", "--known-sites", "b.vcf.gz", "-O", "t.table", "--threads", "8"]));
        let v = run("ApplyBQSR", &["-R", "r.fa", "-I", "x.bam", "-bqsr", "t.table", "-O", "y.bam", "--preserve-qscores-less-than", "10"]);
        assert!(v.unwrap_err().contains("--preserve-qscores-less-than 10"));
        // WARP's command line: booleans become flags (or nothing when false)
        let v = run("ApplyBQSR", &["--create-output-bam-md5", "--add-output-sam-program-record", "-R", "r.fa", "-I", "x.bam", "--use-original-qualities",
            "-O", "y.bam", "-bqsr", "t.table", "--static-quantized-quals", "10", "--static-quantized-quals", "20", "--round-down-quantized", "false"]).unwrap();
        assert_eq!(v, s(&["applybqsr", "--create-output-bam-md5", "-R", "r.fa", "-I", "x.bam", "--use-original-qualities", "-O", "y.bam",
            "--bqsr-recal-file", "t.table", "--static-quantized-quals", "10", "--static-quantized-quals", "20", "--threads", "8"]));
    }

    #[test]
    fn mutect2_family() {
        let v = run("Mutect2", &["-R", "r.fa", "-I", "t.bam", "-I", "n.bam", "--normal", "N", "-L", "20", "--germline-resource", "g.vcf", "-pon", "p.vcf",
            "--f1r2-tar-gz", "f.tar.gz", "-O", "o.vcf", "--tumor-sample", "T", "--native-pair-hmm-threads", "4", "--normal-lod", "2.2"]).unwrap();
        assert_eq!(v, s(&["mutect2", "-R", "r.fa", "-I", "t.bam", "-I", "n.bam", "--normal-sample", "N", "-L", "20", "--germline-resource", "g.vcf",
            "--panel-of-normals", "p.vcf", "--f1r2-tar-gz", "f.tar.gz", "-O", "o.vcf"]));
        let e = run("Mutect2", &["-R", "r.fa", "-I", "t.bam", "-O", "o.vcf", "--mitochondria-mode", "--max-population-af", "0.02", "-emit-lod", "3.0"]).unwrap_err();
        assert!(e.contains("--max-population-af 0.02") && e.contains("--tumor-lod-to-emit with --mitochondria-mode"), "{e}");
        // AlignAndCall.wdl's mitochondria command
        let v = run("Mutect2", &["-R", "r.fa", "-I", "m.bam", "-L", "chrM", "--read-filter", "MateOnSameContigOrNoMappedMateReadFilter",
            "--read-filter", "MateUnmappedAndUnmappedReadFilter", "-O", "o.vcf", "--bam-output", "bamout.bam", "--annotation", "StrandBiasBySample",
            "--mitochondria-mode", "--max-reads-per-alignment-start", "75", "--max-mnp-distance", "0"]).unwrap();
        assert_eq!(v, s(&["mutect2", "-R", "r.fa", "-I", "m.bam", "-L", "chrM", "--read-filter", "MateOnSameContigOrNoMappedMateReadFilter",
            "--read-filter", "MateUnmappedAndUnmappedReadFilter", "-O", "o.vcf", "--bam-output", "bamout.bam", "--annotation", "StrandBiasBySample",
            "--mitochondria-mode", "--max-reads-per-alignment-start", "75", "--max-mnp-distance", "0"]));
        assert!(run("Mutect2", &["-R", "r.fa", "-I", "m.bam", "-O", "o.vcf", "-RF", "ProperlyPairedReadFilter"]).unwrap_err().contains("only MateOnSameContig"));
        let v = run("FilterMutectCalls", &["-R", "r.fa", "-V", "m.vcf", "--ob-priors", "p.tar.gz", "--contamination-table", "c.table", "-O", "f.vcf"]).unwrap();
        assert_eq!(v, s(&["filtermutectcalls", "--reference", "r.fa", "--variant", "m.vcf", "--orientation-bias-artifact-priors", "p.tar.gz",
            "--contamination-table", "c.table", "--output", "f.vcf"]));
        assert!(run("FilterMutectCalls", &["-V", "m.vcf", "-O", "f.vcf"]).unwrap_err().contains("--reference (-R): required"));
        // AlignAndCall.wdl's Filter task
        let v = run("FilterMutectCalls", &["-V", "m.vcf", "-R", "r.fa", "-O", "f.vcf", "--stats", "m.stats", "--max-alt-allele-count", "4",
            "--mitochondria-mode", "--min-allele-fraction", "0", "--f-score-beta", "1.5", "--contamination-estimate", "0.0"]).unwrap();
        assert_eq!(v, s(&["filtermutectcalls", "--variant", "m.vcf", "--reference", "r.fa", "--output", "f.vcf", "--stats", "m.stats",
            "--max-alt-allele-count", "4", "--mitochondria-mode", "--min-allele-fraction", "0", "--f-score-beta", "1.5", "--contamination-estimate", "0.0"]));
        let v = run("GetPileupSummaries", &["-I", "t.bam", "-V", "s.vcf.gz", "-L", "s.vcf.gz", "-O", "p.table", "-min-af", "0.02"]).unwrap();
        assert_eq!(v, s(&["getpileupsummaries", "--input", "t.bam", "--variant", "s.vcf.gz", "--intervals", "s.vcf.gz", "--output", "p.table",
            "--minimum-population-allele-frequency", "0.02"]));
        // interval lists / BED (sarek's scatter) are taken as well; other feature files are not
        assert_eq!(
            run("GetPileupSummaries", &["-I", "t.bam", "-V", "s.vcf.gz", "-L", "x.interval_list", "-O", "p.table"]).unwrap(),
            s(&["getpileupsummaries", "--input", "t.bam", "--variant", "s.vcf.gz", "--intervals", "x.interval_list", "--output", "p.table"])
        );
        assert!(run("GetPileupSummaries", &["-I", "t.bam", "-V", "s.vcf.gz", "-L", "x.gff", "-O", "p.table"]).unwrap_err().contains("sites VCF"));
        let v = run("CalculateContamination", &["-I", "t.table", "-matched", "n.table", "-segments", "s.table", "-O", "c.table"]).unwrap();
        assert_eq!(v, s(&["calculatecontamination", "--input", "t.table", "--matched-normal", "n.table", "--tumor-segmentation", "s.table", "--output", "c.table"]));
        let v = run("LearnReadOrientationModel", &["-I", "a.tar.gz", "-I", "b.tar.gz", "-O", "p.tar.gz", "--num-em-iterations", "20"]).unwrap();
        assert_eq!(v, s(&["learnreadorientationmodel", "--input", "a.tar.gz", "--input", "b.tar.gz", "--output", "p.tar.gz"]));
    }
}
