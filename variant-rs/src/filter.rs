//! Hard filtering (`variant-rs filter`): GATK's Best Practices single-sample hard filters, i.e. the
//! output of `SelectVariants -select-type SNP` / `-select-type INDEL` -> `VariantFiltration` (one
//! `--filter-name` per expression) -> `MergeVcfs`, in one pass.
//!
//! * Record type as htsjdk's `VariantContext.getType()`: per ALT allele SNP (same length 1), MNP (same
//!   length > 1), INDEL (different length), SYMBOLIC; MIXED if the ALT alleles differ. Only SNP and
//!   INDEL records survive the two selections (MNP, MIXED and SYMBOLIC records are dropped, as in
//!   the GATK workflow; `--keep-other` passes them through unfiltered instead).
//! * Each expression `ANN op value` is evaluated like JEXL on htsjdk's VariantContext: INFO values
//!   compared as doubles, QUAL from the QUAL column (missing = -10), a missing annotation never fails.
//! * FILTER: existing filters plus the failing names, sorted (htsjdk VCFEncoder); PASS if none.
//! * Header: as htsjdk writes the merged header: `##fileformat` first, the other lines sorted as
//!   strings, contig lines in their own order; new `##FILTER` lines with the expression as description.

use crate::output::OutFile;
use std::io::{BufRead, BufReader, Read};

#[derive(Clone, Copy, PartialEq)]
enum Op {
    Lt,
    Gt,
}

struct Filter {
    name: &'static str,
    ann: &'static str,
    op: Op,
    value: f64,
    expr: &'static str,
}

const SNP_FILTERS: &[Filter] = &[
    Filter { name: "QD2", ann: "QD", op: Op::Lt, value: 2.0, expr: "QD < 2.0" },
    Filter { name: "QUAL30", ann: "QUAL", op: Op::Lt, value: 30.0, expr: "QUAL < 30.0" },
    Filter { name: "SOR3", ann: "SOR", op: Op::Gt, value: 3.0, expr: "SOR > 3.0" },
    Filter { name: "FS60", ann: "FS", op: Op::Gt, value: 60.0, expr: "FS > 60.0" },
    Filter { name: "MQ40", ann: "MQ", op: Op::Lt, value: 40.0, expr: "MQ < 40.0" },
    Filter { name: "MQRankSum-12.5", ann: "MQRankSum", op: Op::Lt, value: -12.5, expr: "MQRankSum < -12.5" },
    Filter { name: "ReadPosRankSum-8", ann: "ReadPosRankSum", op: Op::Lt, value: -8.0, expr: "ReadPosRankSum < -8.0" },
];

const INDEL_FILTERS: &[Filter] = &[
    Filter { name: "QD2", ann: "QD", op: Op::Lt, value: 2.0, expr: "QD < 2.0" },
    Filter { name: "QUAL30", ann: "QUAL", op: Op::Lt, value: 30.0, expr: "QUAL < 30.0" },
    Filter { name: "FS200", ann: "FS", op: Op::Gt, value: 200.0, expr: "FS > 200.0" },
    Filter { name: "ReadPosRankSum-20", ann: "ReadPosRankSum", op: Op::Lt, value: -20.0, expr: "ReadPosRankSum < -20.0" },
];

#[derive(PartialEq, Clone, Copy, Debug)]
enum VType {
    NoVariation,
    Snp,
    Mnp,
    Indel,
    Symbolic,
    Mixed,
}

fn is_symbolic(a: &[u8]) -> bool {
    a.first() == Some(&b'<') || a.contains(&b'[') || a.contains(&b']') || (a.len() > 1 && (a[0] == b'.' || a[a.len() - 1] == b'.'))
}

/// `VariantContext.getType()`
fn variant_type(refa: &[u8], alts: &[u8]) -> VType {
    if alts == b"." || alts.is_empty() {
        return VType::NoVariation;
    }
    let mut t: Option<VType> = None;
    for a in alts.split(|&c| c == b',') {
        let bt = if is_symbolic(a) {
            VType::Symbolic
        } else if a.len() == refa.len() {
            if a.len() == 1 { VType::Snp } else { VType::Mnp }
        } else {
            VType::Indel
        };
        match t {
            None => t = Some(bt),
            Some(x) if x != bt => return VType::Mixed,
            _ => {}
        }
    }
    t.unwrap_or(VType::NoVariation)
}

/// JEXL `toDouble` of an attribute string (Java `Double.parseDouble`; "" is 0).
fn jexl_double(s: &str) -> Option<f64> {
    if s.is_empty() {
        return Some(0.0);
    }
    match s {
        "NaN" => Some(f64::NAN),
        "Infinity" | "+Infinity" => Some(f64::INFINITY),
        "-Infinity" => Some(f64::NEG_INFINITY),
        _ => s.trim_end_matches(['d', 'D', 'f', 'F']).parse::<f64>().ok(),
    }
}

/// Header lines in htsjdk's merged-header order.
fn sorted_header(meta: Vec<String>) -> Vec<String> {
    let mut fileformat = None;
    let mut rest: Vec<(usize, String)> = Vec::new();
    let mut seen = std::collections::HashSet::new();
    for (i, l) in meta.into_iter().enumerate() {
        if l.starts_with("##fileformat=") {
            fileformat = Some(l);
        } else if seen.insert(l.clone()) {
            rest.push((i, l));
        }
    }
    // our command lines (variant-rs, or ard-rs before 2026-11-15) sort where GATK's command lines do
    let key = |l: &str| -> String {
        for p in ["##variant-rs-filterCommandLine", "##variant-rsCommandLine", "##ard-rs-filterCommandLine", "##ard-rsCommandLine"] {
            if let Some(rest) = l.strip_prefix(p) {
                return format!("GATKCommandLine{rest}");
            }
        }
        l[2..].to_string()
    };
    rest.sort_by(|a, b| {
        let ca = a.1.starts_with("##contig=");
        let cb = b.1.starts_with("##contig=");
        if ca && cb {
            a.0.cmp(&b.0)
        } else {
            key(&a.1).cmp(&key(&b.1))
        }
    });
    let mut out: Vec<String> = fileformat.into_iter().collect();
    out.extend(rest.into_iter().map(|x| x.1));
    out
}

pub fn run(argv: &[String]) {
    const USAGE: &str = "usage: variant-rs filter -V in.vcf[.gz] -O out.vcf.gz [--keep-other] [--compression-level N (6)]";
    let (mut input, mut out, mut keep_other, mut level) = (String::new(), String::new(), false, 6);
    let mut i = 0;
    while i < argv.len() {
        match argv[i].as_str() {
            "-V" | "--variant" => {
                input = argv[i + 1].clone();
                i += 2;
            }
            "-O" | "--output" => {
                out = argv[i + 1].clone();
                i += 2;
            }
            "--keep-other" => {
                keep_other = true;
                i += 1;
            }
            "--compression-level" => {
                level = argv[i + 1].parse().expect("--compression-level");
                i += 2;
            }
            x => panic!("unknown argument {x}\n{USAGE}"),
        }
    }
    assert!(!input.is_empty() && !out.is_empty(), "{USAGE}");
    let t0 = std::time::Instant::now();
    let f = std::fs::File::open(&input).unwrap_or_else(|e| panic!("{input}: {e}"));
    let rd: Box<dyn Read> = if input.ends_with(".gz") { Box::new(flate2::read::MultiGzDecoder::new(f)) } else { Box::new(f) };
    let mut rd = BufReader::with_capacity(4 << 20, rd);
    // header
    let mut meta: Vec<String> = Vec::new();
    let chrom_line;
    loop {
        let mut l = String::new();
        assert!(rd.read_line(&mut l).unwrap() > 0, "{input}: no #CHROM line");
        let l = l.trim_end_matches(['\n', '\r']).to_string();
        if l.starts_with("##") {
            meta.push(l);
        } else {
            chrom_line = l;
            break;
        }
    }
    for f in SNP_FILTERS.iter().chain(INDEL_FILTERS) {
        meta.push(format!("##FILTER=<ID={},Description=\"{}\">", f.name, f.expr));
    }
    let cmdline = std::env::args().collect::<Vec<_>>().join(" ");
    meta.push(format!("##variant-rs-filterCommandLine=<ID=variant-rs filter,CommandLine=\"{}\",Version=\"{}\">", cmdline.replace('"', "\\\""), env!("CARGO_PKG_VERSION")));
    let mut header = sorted_header(meta).join("\n");
    header.push('\n');
    header.push_str(&chrom_line);
    header.push('\n');
    let level_opt = if crate::output::block_compressed(&out) { Some(level) } else { None };
    let mut of = OutFile::create(&out, &header, level_opt).unwrap_or_else(|e| panic!("{out}: {e}"));

    let (mut n_in, mut n_out, mut n_dropped, mut n_pass) = (0u64, 0u64, 0u64, 0u64);
    let mut buf: Vec<u8> = Vec::with_capacity(8 << 20);
    let mut line = Vec::new();
    loop {
        line.clear();
        if rd.read_until(b'\n', &mut line).unwrap() == 0 {
            break;
        }
        if line.last() == Some(&b'\n') {
            line.pop();
        }
        if line.is_empty() {
            continue;
        }
        n_in += 1;
        let cols: Vec<&[u8]> = line.splitn(9, |&c| c == b'\t').collect();
        let vt = variant_type(cols[3], cols[4]);
        let set: &[Filter] = match vt {
            VType::Snp => SNP_FILTERS,
            VType::Indel => INDEL_FILTERS,
            _ if keep_other => &[],
            _ => {
                n_dropped += 1;
                continue;
            }
        };
        let old_filter = std::str::from_utf8(cols[6]).unwrap();
        let mut filters: Vec<&str> = Vec::new();
        if old_filter != "." && old_filter != "PASS" {
            filters.extend(old_filter.split(';'));
        }
        let info = std::str::from_utf8(cols[7]).unwrap();
        let qual = std::str::from_utf8(cols[5]).unwrap();
        for f in set {
            let v = if f.ann == "QUAL" {
                Some(if qual == "." { -10.0 } else { qual.parse::<f64>().unwrap() })
            } else {
                info.split(';').find_map(|kv| kv.strip_prefix(f.ann).and_then(|r| r.strip_prefix('='))).and_then(jexl_double)
            };
            if let Some(v) = v {
                let fail = match f.op {
                    Op::Lt => v < f.value,
                    Op::Gt => v > f.value,
                };
                if fail && !filters.contains(&f.name) {
                    filters.push(f.name);
                }
            }
        }
        let new_filter = if filters.is_empty() {
            // VariantFiltration sets filters on every record: unfiltered becomes PASS
            if set.is_empty() { old_filter.to_string() } else { "PASS".to_string() }
        } else {
            filters.sort_unstable();
            filters.dedup();
            filters.join(";")
        };
        if new_filter == "PASS" {
            n_pass += 1;
        }
        for (k, c) in cols.iter().enumerate() {
            if k > 0 {
                buf.push(b'\t');
            }
            if k == 6 {
                buf.extend_from_slice(new_filter.as_bytes());
            } else {
                buf.extend_from_slice(c);
            }
        }
        buf.push(b'\n');
        n_out += 1;
        if buf.len() >= 8 << 20 {
            of.add_text(&buf).unwrap();
            buf.clear();
        }
    }
    of.add_text(&buf).unwrap();
    of.finish().unwrap();
    eprintln!(
        "records: {n_in} in, {n_out} out ({n_pass} PASS), {n_dropped} dropped (not SNP/INDEL{}) in {:.1}s",
        if keep_other { "" } else { "; --keep-other keeps them" },
        t0.elapsed().as_secs_f64()
    );
}
