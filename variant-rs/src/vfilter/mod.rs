//! VariantFiltration (`variant-rs variantfiltration`, GATK 4.7) and Picard MakeSitesOnlyVcf (`variant-rs
//! makesitesonlyvcf`), as WARP's joint calling runs them (`ExcessHet > 54.69`, an exome target mask, then sites
//! only), its DRAGEN mode (`QUAL < 10.4139`) and GATK's hard-filtering examples (several expressions, `||`).
//!
//! VariantFiltration: each record's INFO is decoded as htsjdk's `parseInfo` does (a value is a String, a list of
//! Strings, or a flag), the mask and the expressions (`jexl.rs`) are applied, and the record is re-encoded as htsjdk's
//! writer does: FILTER the old filters plus the new ones, sorted (PASS when none), INFO keys sorted, QUAL and alleles
//! as htsjdk prints them, the genotype columns as read. The header gains the FILTER lines, `##source` and the
//! command line; samples keep their order.
//! Refused: genotype filters (`--genotype-filter-expression`, `--set-filtered-genotype-to-no-call`), SNP clusters
//! (`--cluster-window-size`), allele-specific filters, intervals (`-L` / `-XL`).

pub mod jexl;

use crate::vqsr::vcfio::{self, Lines};
use crate::vcfinfo::{self, Info};
use jexl::{EvalErr, Expr, Val};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::Write;

const USAGE: &str = "usage: variant-rs variantfiltration -V in.vcf.gz -O out.vcf.gz [--filter-expression EXP --filter-name NAME]... \
[--missing-values-evaluate-as-failing] [--invert-filter-expression] [--invalidate-previous-filters] [--mask file \
[--mask-name N] [--mask-description D] [--filter-not-in-mask] [--mask-extension N]] (GATK's names)";

fn refuse(why: &str) -> ! {
    crate::gatk_cli::refuse(why)
}

/// GATK fails on this input too (an exception it does not catch, or a user error): exit with its message.
fn fail(why: &str) -> ! {
    eprintln!("variant-rs VariantFiltration: {why}");
    std::process::exit(1)
}

// ---------------------------------------------------------------- records

/// One record as htsjdk's `VCFEncoder` writes it: CHROM POS ID, REF / ALT as `Allele` prints them, QUAL, `filter`,
/// INFO (`vcfinfo`; a key the header does not declare fails the run), then the genotype columns as read (none when
/// `sites_only`).
fn encode(cols: &[&str], filter: &str, info: &HashMap<&str, Info>, declared: &HashMap<String, bool>, sites_only: bool, out: &mut Vec<u8>) {
    use crate::mergevcfs::{allele, format_qual};
    write!(out, "{}\t{}\t{}\t{}\t", cols[0], cols[1], cols[2], allele(cols[3])).unwrap();
    if cols[4] == "." {
        out.push(b'.');
    } else {
        out.extend_from_slice(cols[4].split(',').map(allele).collect::<Vec<_>>().join(",").as_bytes());
    }
    write!(out, "\t{}\t{}\t", format_qual(cols[5]), filter).unwrap();
    vcfinfo::write(info, declared, &format!("{}:{}", cols[0], cols[1]), out).unwrap_or_else(|e| fail(&e));
    if !sites_only {
        for c in &cols[8.min(cols.len())..] {
            out.push(b'\t');
            out.extend_from_slice(c.as_bytes());
        }
    }
    out.push(b'\n');
}

/// The input's header: `##` lines (standard lines repaired as htsjdk reads them), sample names, INFO flag types.
struct InputHeader {
    meta: Vec<String>,
    samples: Vec<String>,
    /// INFO ID -> whether the header declares it a Flag
    flags: HashMap<String, bool>,
}

fn read_header(lines: &mut Lines, path: &str) -> InputHeader {
    let mut line = Vec::new();
    let mut meta = Vec::new();
    loop {
        if !lines.next_line(&mut line) {
            fail(&format!("{path}: no #CHROM line"));
        }
        let l = String::from_utf8_lossy(&line).into_owned();
        if l.starts_with("##") {
            meta.push(crate::vcfheader::repair_standard_line(l));
        } else {
            let samples = l.split('\t').skip(9).map(|s| s.to_string()).collect();
            let flags = vcfinfo::declared(&meta);
            return InputHeader { meta, samples, flags };
        }
    }
}

// ---------------------------------------------------------------- the mask

/// The `--mask` features by contig: intervals sorted by start, with the running maximum of their ends.
struct Mask {
    by_contig: HashMap<String, (Vec<(i64, i64)>, Vec<i64>)>,
}

impl Mask {
    /// A BED file (0-based starts), a Picard interval list, or a VCF (each record from POS to its end), as GATK's
    /// BEDCodec, IntervalListCodec and VCFCodec give them.
    fn load(path: &str) -> Mask {
        let mut v: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
        let bed = path.ends_with(".bed") || path.ends_with(".bed.gz");
        let ilist = path.ends_with(".interval_list") || path.ends_with(".interval_list.gz");
        let vcf = [".vcf", ".vcf.gz", ".vcf.bgz"].iter().any(|e| path.ends_with(e));
        if !(bed || ilist || vcf) {
            refuse(&format!("--mask {path}: only .bed, .interval_list and VCF masks are supported"));
        }
        // GATK queries the mask by interval, which needs an index (IndexFeatureFile's .idx, or a .tbi)
        if ![".idx", ".tbi"].iter().any(|x| std::path::Path::new(&format!("{path}{x}")).exists()) {
            fail(&format!(
                "Input {path} must support random access to enable queries by interval. If it's a file, please index it using the bundled tool IndexFeatureFile"
            ));
        }
        let mut lines = Lines::open(path);
        let mut line = Vec::new();
        let num = |s: &str| -> i64 { s.parse().unwrap_or_else(|_| fail(&format!("{path}: bad position {s}"))) };
        while lines.next_line(&mut line) {
            let l = String::from_utf8_lossy(&line);
            if l.is_empty() || l.starts_with('#') || l.starts_with('@') || l.starts_with("track") || l.starts_with("browser") {
                continue;
            }
            let f: Vec<&str> = l.split('\t').collect();
            if f.len() < 3 {
                fail(&format!("{path}: line with fewer than 3 columns: {l}"));
            }
            let (s, e) = if bed {
                (num(f[1]) + 1, num(f[2]))
            } else if ilist {
                (num(f[1]), num(f[2]))
            } else {
                let s = num(f[1]);
                (s, record_end(s, f[3], f.get(7).copied().unwrap_or(".")))
            };
            v.entry(f[0].to_string()).or_default().push((s, e));
        }
        let by_contig = v
            .into_iter()
            .map(|(c, mut iv)| {
                iv.sort();
                let mut m = i64::MIN;
                let max_end = iv.iter().map(|&(_, e)| { m = m.max(e); m }).collect();
                (c, (iv, max_end))
            })
            .collect();
        Mask { by_contig }
    }

    /// Whether any feature overlaps [start, end].
    fn overlaps(&self, contig: &str, start: i64, end: i64) -> bool {
        let Some((iv, max_end)) = self.by_contig.get(contig) else { return false };
        let n = iv.partition_point(|&(s, _)| s <= end);
        n > 0 && max_end[n - 1] >= start
    }
}

/// A record's end: its END attribute, else POS + |REF| - 1 (as htsjdk's VCF codec sets it).
fn record_end(start: i64, refa: &str, info: &str) -> i64 {
    match vcfio::info_double(info.as_bytes(), "END") {
        Ok(e) if !e.is_nan() => e as i64,
        _ => start + refa.len() as i64 - 1,
    }
}

// ---------------------------------------------------------------- VariantFiltration

struct Args {
    input: String,
    output: String,
    names: Vec<String>,
    exprs: Vec<String>,
    missing_fails: bool,
    invert: bool,
    invalidate: bool,
    mask: Option<String>,
    mask_name: String,
    mask_description: Option<String>,
    not_in_mask: bool,
    mask_extension: i64,
    compression_level: i32,
}

fn parse_args(argv: &[String]) -> Args {
    let mut a = Args {
        input: String::new(),
        output: String::new(),
        names: vec![],
        exprs: vec![],
        missing_fails: false,
        invert: false,
        invalidate: false,
        mask: None,
        mask_name: "Mask".into(),
        mask_description: None,
        not_in_mask: false,
        mask_extension: 0,
        compression_level: 6,
    };
    let val = |i: usize| -> String { argv.get(i + 1).cloned().unwrap_or_else(|| refuse(&format!("{} needs a value\n{USAGE}", argv[i]))) };
    let int = |i: usize| -> i64 { val(i).parse().unwrap_or_else(|_| refuse(&format!("{}: not an integer", argv[i]))) };
    let flag = |i: &mut usize| -> bool {
        match argv.get(*i + 1).map(|s| s.to_ascii_lowercase()) {
            Some(v) if ["true", "false"].contains(&v.as_str()) => {
                *i += 1;
                v == "true"
            }
            _ => true,
        }
    };
    let (mut cluster_size, mut cluster_window) = (3, 0);
    let mut i = 0;
    while i < argv.len() {
        match argv[i].as_str() {
            "-V" | "--variant" => {
                a.input = val(i);
                i += 1;
            }
            "-O" | "--output" => {
                a.output = val(i);
                i += 1;
            }
            "-filter" | "--filter-expression" => {
                a.exprs.push(val(i));
                i += 1;
            }
            "--filter-name" => {
                a.names.push(val(i));
                i += 1;
            }
            "--missing-values-evaluate-as-failing" => a.missing_fails = flag(&mut i),
            "-invfilter" | "--invert-filter-expression" => a.invert = flag(&mut i),
            "--invalidate-previous-filters" => a.invalidate = flag(&mut i),
            "--mask" => {
                a.mask = Some(val(i));
                i += 1;
            }
            "--mask-name" => {
                a.mask_name = val(i);
                i += 1;
            }
            "--mask-description" => {
                a.mask_description = Some(val(i));
                i += 1;
            }
            "--filter-not-in-mask" => a.not_in_mask = flag(&mut i),
            "--mask-extension" => {
                a.mask_extension = int(i);
                i += 1;
            }
            "-cluster" | "--cluster-size" => {
                cluster_size = int(i);
                i += 1;
            }
            "-window" | "--cluster-window-size" => {
                cluster_window = int(i);
                i += 1;
            }
            "--compression-level" => {
                a.compression_level = val(i).parse().unwrap_or_else(|_| refuse("--compression-level: not a number"));
                i += 1;
            }
            // the output header comes from the input's; the reference does not change it
            "-R" | "--reference" => i += 1,
            x => refuse(&format!("VariantFiltration: {x} is not supported\n{USAGE}")),
        }
        i += 1;
    }
    if a.input.is_empty() || a.output.is_empty() {
        refuse(USAGE);
    }
    // onTraversalStart's checks, in its order
    if cluster_size <= 1 {
        fail("Argument cluster-size has a bad value: values lower than 2 are not allowed");
    }
    if cluster_window > 0 {
        refuse("VariantFiltration --cluster-window-size (SNP clusters) is not supported");
    }
    if a.mask_extension < 0 {
        fail("Argument mask-extension has a bad value: negative values are not allowed");
    }
    if a.not_in_mask && a.mask.is_none() {
        fail("Argument filter-not-in-mask has a bad value: argument not allowed if mask argument is not provided");
    }
    if a.mask_description.is_some() && a.mask.is_none() {
        fail("Argument mask-description has a bad value: argument not allowed if mask argument is not provided");
    }
    if a.names.len() != a.exprs.len() {
        fail(&format!("Inconsistent number of provided filter names and expressions: names={:?} exps={:?}", a.names, a.exprs));
    }
    a
}

/// The variables of a record (`VariantJEXLContext`): CHROM, POS, TYPE, QUAL, N_ALLELES, FILTER ("1" / "0"), then
/// the INFO attributes, then the record's filter names ("1").
struct VcContext<'a> {
    cols: &'a [&'a str],
    info: &'a HashMap<&'a str, Info>,
    filters: &'a BTreeSet<String>,
}

impl jexl::Context for VcContext<'_> {
    fn get(&self, name: &str) -> Option<Val> {
        let c = self.cols;
        Some(match name {
            "CHROM" => Val::Str(c[0].to_string()),
            "POS" => Val::Int(c[1].parse().unwrap_or(0)),
            "TYPE" => Val::Str(type_name(vcfio::variant_type(c[3].as_bytes(), c[4].as_bytes())).into()),
            "QUAL" => Val::Double(qual(c[5])),
            "N_ALLELES" => Val::Int(1 + if c[4] == "." { 0 } else { c[4].split(',').count() as i64 }),
            "FILTER" => Val::Str(if self.filters.is_empty() { "0" } else { "1" }.into()),
            _ => match self.info.get(name) {
                Some(Info::Str(s)) => Val::Str(s.clone()),
                Some(Info::List(l)) => Val::List(l.clone()),
                Some(Info::Flag) => Val::Bool(true),
                None if self.filters.contains(name) => Val::Str("1".into()),
                None => return None,
            },
        })
    }
}

fn type_name(t: vcfio::VType) -> &'static str {
    use crate::vqsr::vcfio::VType::*;
    match t {
        NoVariation => "NO_VARIATION",
        Snp => "SNP",
        Mnp => "MNP",
        Indel => "INDEL",
        Symbolic => "SYMBOLIC",
        Mixed => "MIXED",
    }
}

/// `-10 * getLog10PError()` after `parseQual`: '.' and the VCF 3 missing value (-1) are log10PError 1.0.
fn qual(q: &str) -> f64 {
    if q == "." {
        return -10.0;
    }
    let v = vcfio::parse_vcf_double(q).unwrap_or_else(|| fail(&format!("bad QUAL {q}")));
    if v < 0.0 && (v - -1.0).abs() < 5.0e-5 {
        return -10.0;
    }
    -10.0 * (v / -10.0)
}

pub fn run(argv: &[String]) {
    let a = parse_args(argv);
    let t0 = std::time::Instant::now();
    // initializeMatchExps: a name given twice keeps its last expression
    let mut by_name: BTreeMap<String, String> = BTreeMap::new();
    for (n, e) in a.names.iter().zip(&a.exprs) {
        by_name.insert(n.clone(), e.clone());
    }
    let exps: Vec<(String, String, Expr)> = by_name
        .into_iter()
        .map(|(n, e)| {
            let x = jexl::parse(&e).unwrap_or_else(|why| refuse(&why));
            (n, e, x)
        })
        .collect();
    let mask = a.mask.as_ref().map(|m| Mask::load(m));

    let mut lines = Lines::open(&a.input);
    let h = read_header(&mut lines, &a.input);
    crate::vcfheader::require_sorted_samples(&h.samples, &a.input);
    let mut extra: Vec<String> = Vec::new();
    let inverse = |d: &str| if a.invert { format!("Inverse of: {d}") } else { d.to_string() };
    for (n, e, _) in &exps {
        extra.push(crate::vcfheader::filter_line(n, &inverse(e)));
    }
    if mask.is_some() {
        let d = match (&a.mask_description, a.not_in_mask) {
            (Some(d), _) => d.clone(),
            (None, true) => "Doesn't overlap a user-input mask".into(),
            (None, false) => "Overlaps a user-input mask".into(),
        };
        extra.push(crate::vcfheader::filter_line(&a.mask_name, &d));
    }
    extra.push("##source=VariantFiltration".into());
    extra.push(crate::vcfheader::command_line("VariantFiltration"));
    let hdr = crate::vcfheader::write_header(&h.meta, extra, &h.samples);
    let level = if crate::output::block_compressed(&a.output) { Some(a.compression_level) } else { None };
    let mut of = crate::output::OutFile::create(&a.output, &hdr, level).unwrap_or_else(|e| refuse(&format!("{}: {e}", a.output)));

    let mut line = Vec::new();
    let mut buf: Vec<u8> = Vec::with_capacity(8 << 20);
    let (mut n_in, mut n_filtered) = (0u64, 0u64);
    while lines.next_line(&mut line) {
        if line.is_empty() || line[0] == b'#' {
            continue;
        }
        n_in += 1;
        let text = std::str::from_utf8(&line).unwrap_or_else(|_| fail(&format!("{}: not UTF-8", a.input)));
        let cols: Vec<&str> = text.split('\t').collect();
        if cols.len() < 8 {
            fail(&format!("{}: VCF line with fewer than 8 columns", a.input));
        }
        let info = vcfinfo::parse(cols[7], &h.flags);
        // the record's filters (PASS and '.' are none), cleared by --invalidate-previous-filters, then the mask
        let mut filters: BTreeSet<String> = BTreeSet::new();
        if !a.invalidate && cols[6] != "." && cols[6] != "PASS" {
            filters.extend(cols[6].split(';').map(|s| s.to_string()));
        }
        if let Some(m) = &mask {
            if !filters.contains(&a.mask_name) {
                let start: i64 = cols[1].parse().unwrap_or_else(|_| fail(&format!("{}: bad POS {}", a.input, cols[1])));
                let end = record_end(start, cols[3], cols[7]);
                let q = (start - a.mask_extension).max(1);
                if m.overlaps(cols[0], q, end + a.mask_extension) != a.not_in_mask {
                    filters.insert(a.mask_name.clone());
                }
            }
        }
        let mut out = filters.clone();
        let ctx = VcContext { cols: &cols, info: &info, filters: &filters };
        for (n, _, x) in &exps {
            let m = match jexl::evaluate(x, &ctx) {
                Ok(b) => b,
                Err(EvalErr::Missing) => a.missing_fails,
                Err(EvalErr::Fatal(why)) => fail(&format!("{why} (at {}:{})", cols[0], cols[1])),
            };
            if m != a.invert {
                out.insert(n.clone());
            }
        }
        let f = if out.is_empty() {
            if a.invalidate { ".".to_string() } else { "PASS".to_string() }
        } else {
            n_filtered += 1;
            out.into_iter().collect::<Vec<_>>().join(";")
        };
        encode(&cols, &f, &info, &h.flags, false, &mut buf);
        if buf.len() >= 8 << 20 {
            of.add_text(&buf).unwrap();
            buf.clear();
        }
    }
    of.add_text(&buf).unwrap();
    of.finish().unwrap_or_else(|e| refuse(&format!("{}: {e}", a.output)));
    eprintln!("[variantfiltration] {n_in} records, {n_filtered} filtered, in {:.1}s", t0.elapsed().as_secs_f64());
}

// ---------------------------------------------------------------- MakeSitesOnlyVcf

/// Picard MakeSitesOnlyVcf: the input's header lines without samples (nothing added), each record re-encoded
/// without its genotype columns.
pub fn sites_only(argv: &[String]) {
    let (mut input, mut output, mut level) = (String::new(), String::new(), 2);
    let mut i = 0;
    while i < argv.len() {
        let v = || argv.get(i + 1).cloned().unwrap_or_else(|| refuse(&format!("{} needs a value", argv[i])));
        match argv[i].as_str() {
            "-I" | "--INPUT" => input = v(),
            "-O" | "--OUTPUT" => output = v(),
            "--compression-level" => level = v().parse().unwrap_or_else(|_| refuse("--compression-level: not a number")),
            x => refuse(&format!("MakeSitesOnlyVcf: {x} is not supported")),
        }
        i += 2;
    }
    if input.is_empty() || output.is_empty() {
        refuse("usage: variant-rs makesitesonlyvcf -I in.vcf.gz -O out.vcf.gz");
    }
    let t0 = std::time::Instant::now();
    let mut lines = Lines::open(&input);
    let h = read_header(&mut lines, &input);
    let hdr = crate::vcfheader::write_header(&h.meta, vec![], &[]);
    let lvl = if crate::output::block_compressed(&output) { Some(level) } else { None };
    let mut of = crate::output::OutFile::create(&output, &hdr, lvl).unwrap_or_else(|e| refuse(&format!("{output}: {e}")));
    let mut line = Vec::new();
    let mut buf: Vec<u8> = Vec::with_capacity(8 << 20);
    let mut n = 0u64;
    while lines.next_line(&mut line) {
        if line.is_empty() || line[0] == b'#' {
            continue;
        }
        n += 1;
        let text = std::str::from_utf8(&line).unwrap_or_else(|_| fail(&format!("{input}: not UTF-8")));
        let cols: Vec<&str> = text.split('\t').collect();
        if cols.len() < 8 {
            fail(&format!("{input}: VCF line with fewer than 8 columns"));
        }
        let info = vcfinfo::parse(cols[7], &h.flags);
        let f = if cols[6] == "." || cols[6] == "PASS" {
            cols[6].to_string()
        } else {
            cols[6].split(';').collect::<BTreeSet<_>>().into_iter().collect::<Vec<_>>().join(";")
        };
        encode(&cols, &f, &info, &h.flags, true, &mut buf);
        if buf.len() >= 8 << 20 {
            of.add_text(&buf).unwrap();
            buf.clear();
        }
    }
    of.add_text(&buf).unwrap();
    of.finish().unwrap_or_else(|e| refuse(&format!("{output}: {e}")));
    eprintln!("[makesitesonlyvcf] {n} records in {:.1}s", t0.elapsed().as_secs_f64());
}
