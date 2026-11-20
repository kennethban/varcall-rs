//! ApplyVQSR (`variant-rs applyvqsr`): GATK 4.7's site-specific filtering from a recal file and tranches. Each
//! input record of the mode's class that is not filtered gets the matching recal record's VQSLOD (re-printed as
//! htsjdk prints a Double: `formatVCFDouble`), culprit and training labels, and a FILTER from the tranches (or the
//! `--lod-score-cutoff`). Every record is re-encoded as htsjdk's writer does (INFO keys sorted, QUAL, alleles,
//! FILTER). The recal file is read alongside the input; both are in dictionary order.
//! Not supported (refused): allele-specific mode (`-AS`).

use super::tranche::{read_tranches, Tranche};
use super::vcfio::{self, Lines};
use super::Mode;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::Write;

const USAGE: &str = "usage: variant-rs applyvqsr -V in.vcf.gz --recal-file x.recal.vcf.gz [--tranches-file x.tranches \
--truth-sensitivity-filter-level T | --lod-score-cutoff L] -mode SNP|INDEL|BOTH -O out.vcf.gz [--ignore-filter F] \
[--ignore-all-filters] [--exclude-filtered] (GATK's names)";

fn refuse(why: &str) -> ! {
    crate::gatk_cli::refuse(why)
}

struct Args {
    input: String,
    recal: String,
    tranches_file: Option<String>,
    output: String,
    ts_filter_level: Option<f64>,
    lod_cutoff: Option<f64>,
    mode: Mode,
    ignore_filters: Vec<String>,
    ignore_all_filters: bool,
    exclude_filtered: bool,
    compression_level: i32,
}

fn parse_args(argv: &[String]) -> Args {
    let mut a = Args {
        input: String::new(),
        recal: String::new(),
        tranches_file: None,
        output: String::new(),
        ts_filter_level: None,
        lod_cutoff: None,
        mode: Mode::Snp,
        ignore_filters: vec![],
        ignore_all_filters: false,
        exclude_filtered: false,
        compression_level: 6,
    };
    let val = |i: usize| -> String { argv.get(i + 1).cloned().unwrap_or_else(|| refuse(&format!("{} needs a value\n{USAGE}", argv[i]))) };
    let dbl = |i: usize| -> f64 { vcfio::parse_vcf_double(&val(i)).unwrap_or_else(|| refuse(&format!("{}: not a number", argv[i]))) };
    let flag = |i: &mut usize| -> bool {
        match argv.get(*i + 1).map(|s| s.to_ascii_lowercase()) {
            Some(v) if ["true", "false"].contains(&v.as_str()) => {
                *i += 1;
                v == "true"
            }
            _ => true,
        }
    };
    let mut i = 0;
    while i < argv.len() {
        match argv[i].as_str() {
            "-V" | "--variant" => {
                if !a.input.is_empty() {
                    refuse("ApplyVQSR with several -V inputs is not supported");
                }
                a.input = val(i);
                i += 1;
            }
            "--recal-file" => {
                a.recal = val(i);
                i += 1;
            }
            "--tranches-file" => {
                a.tranches_file = Some(val(i));
                i += 1;
            }
            "-O" | "--output" => {
                a.output = val(i);
                i += 1;
            }
            "-ts-filter-level" | "--truth-sensitivity-filter-level" => {
                a.ts_filter_level = Some(dbl(i));
                i += 1;
            }
            "--lod-score-cutoff" => {
                a.lod_cutoff = Some(dbl(i));
                i += 1;
            }
            "-mode" | "--mode" => {
                a.mode = Mode::parse(&val(i));
                i += 1;
            }
            "--ignore-filter" => {
                a.ignore_filters.push(val(i));
                i += 1;
            }
            "--ignore-all-filters" => a.ignore_all_filters = flag(&mut i),
            "--exclude-filtered" => a.exclude_filtered = flag(&mut i),
            "--compression-level" => {
                a.compression_level = val(i).parse().unwrap_or_else(|_| refuse("--compression-level: not a number"));
                i += 1;
            }
            // the output header comes from the input's; the reference does not change it
            "-R" | "--reference" => i += 1,
            "-AS" | "--use-allele-specific-annotations" => {
                if flag(&mut i) {
                    refuse("ApplyVQSR -AS (allele-specific mode) is not supported");
                }
            }
            x => refuse(&format!("ApplyVQSR: {x} is not supported\n{USAGE}")),
        }
        i += 1;
    }
    if a.input.is_empty() || a.recal.is_empty() || a.output.is_empty() {
        refuse(USAGE);
    }
    if a.ts_filter_level.is_some() && a.lod_cutoff.is_some() {
        refuse("Arguments --truth-sensitivity-filter-level and --lod-score-cutoff are mutually exclusive. Please only specify one option.");
    }
    a
}

/// The filter for a lod (`generateFilterString`); `tranches` ordered from worst (highest sensitivity) to best.
fn filter_for(lod: f64, tranches: &[Tranche], cutoff: Option<f64>) -> String {
    if let Some(c) = cutoff {
        return if lod < c { "LOW_VQSLOD".into() } else { "PASS".into() };
    }
    for i in (0..tranches.len()).rev() {
        if lod >= tranches[i].min_vqslod {
            return if i == tranches.len() - 1 { "PASS".into() } else { tranches[i].filter_name.clone() };
        }
    }
    format!("{}+", tranches[0].filter_name)
}

/// One recal record: end, VQSLOD text, culprit, labels.
struct Recal {
    end: i64,
    lod: String,
    culprit: Option<String>,
    positive: bool,
    negative: bool,
}

/// The recal file, read alongside the input: records at a (contig, start), contigs in the recal file's order.
struct RecalReader {
    lines: Lines,
    rank: HashMap<String, usize>,
    pending: Option<(usize, i64, Recal)>,
    path: String,
}

impl RecalReader {
    fn open(path: &str) -> RecalReader {
        let mut lines = Lines::open(path);
        let mut line = Vec::new();
        let mut rank = HashMap::new();
        loop {
            if !lines.next_line(&mut line) {
                refuse(&format!("{path}: no #CHROM line"));
            }
            let l = String::from_utf8_lossy(&line);
            if let Some(rest) = l.strip_prefix("##contig=<ID=") {
                let id = rest.split([',', '>']).next().unwrap_or("").to_string();
                let n = rank.len();
                rank.entry(id).or_insert(n);
            }
            if !l.starts_with("##") {
                break;
            }
        }
        let mut r = RecalReader { lines, rank, pending: None, path: path.to_string() };
        r.advance();
        r
    }

    fn advance(&mut self) {
        let mut line = Vec::new();
        self.pending = None;
        while self.lines.next_line(&mut line) {
            if line.is_empty() || line[0] == b'#' {
                continue;
            }
            let cols: Vec<&[u8]> = line.splitn(9, |&c| c == b'\t').collect();
            if cols.len() < 8 {
                refuse(&format!("{}: VCF line with fewer than 8 columns", self.path));
            }
            let chrom = String::from_utf8_lossy(cols[0]);
            let Some(&tid) = self.rank.get(chrom.as_ref()) else { refuse(&format!("{}: contig {chrom} has no ##contig line", self.path)) };
            let start: i64 = std::str::from_utf8(cols[1]).ok().and_then(|s| s.parse().ok()).unwrap_or_else(|| refuse(&format!("{}: bad POS", self.path)));
            let mut rec = Recal { end: start + cols[3].len() as i64 - 1, lod: String::new(), culprit: None, positive: false, negative: false };
            let mut lod = None;
            for kv in cols[7].split(|&c| c == b';') {
                let (k, v) = match memchr::memchr(b'=', kv) {
                    Some(i) => (&kv[..i], Some(&kv[i + 1..])),
                    None => (kv, None),
                };
                let v = v.map(|v| String::from_utf8_lossy(v).into_owned());
                match k {
                    b"END" => rec.end = v.and_then(|s| s.parse().ok()).unwrap_or_else(|| refuse(&format!("{}: the END value in the INFO field is not valid", self.path))),
                    b"VQSLOD" => lod = v,
                    b"culprit" => rec.culprit = v,
                    b"POSITIVE_TRAIN_SITE" => rec.positive = true,
                    b"NEGATIVE_TRAIN_SITE" => rec.negative = true,
                    _ => {}
                }
            }
            rec.lod = lod.unwrap_or_default();
            self.pending = Some((tid, start, rec));
            return;
        }
    }

    /// The recal records at (tid, start); earlier ones are skipped.
    fn at(&mut self, tid: usize, start: i64) -> Vec<Recal> {
        let mut out = Vec::new();
        while let Some((t, s, _)) = &self.pending {
            if (*t, *s) < (tid, start) {
                self.advance();
            } else if (*t, *s) == (tid, start) {
                let (_, _, r) = self.pending.take().unwrap();
                out.push(r);
                self.advance();
            } else {
                break;
            }
        }
        out
    }
}

/// htsjdk's header for the output: the input's lines, the VQSR lines and this run's filters, `##source` and the
/// command line; sorted as text with the contigs as one block (input order), samples sorted.
fn header(input_meta: &[String], extra: Vec<String>, samples: &[String]) -> String {
    let mut lines: Vec<String> = Vec::new();
    let mut seen = HashSet::new();
    for l in input_meta.iter().cloned().chain(extra) {
        if l.starts_with("##fileformat=") {
            continue;
        }
        if seen.insert(l.clone()) {
            lines.push(l);
        }
    }
    let key = |l: &String| -> String {
        if l.starts_with("##contig=<") {
            "##contig=<".into()
        } else {
            l.strip_prefix("##variant-rsCommandLine").map(|r| format!("##GATKCommandLine{r}")).unwrap_or_else(|| l.clone())
        }
    };
    lines.sort_by_key(key);
    let mut out = String::from("##fileformat=VCFv4.2\n");
    for l in lines {
        out.push_str(&l);
        out.push('\n');
    }
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    if !samples.is_empty() {
        out.push_str("\tFORMAT");
        for s in samples {
            out.push('\t');
            out.push_str(s);
        }
    }
    out.push('\n');
    out
}

fn filter_line(id: &str, description: &str) -> String {
    format!("##FILTER=<ID={id},Description=\"{}\">", description.replace('"', "\\\""))
}

/// One record re-encoded as htsjdk writes it: `filter` replaces FILTER when Some, `add` adds or replaces INFO
/// keys (a None value is a flag), `order` permutes the sample columns.
fn reencode(cols: &[&str], filter: Option<&str>, add: &[(&str, Option<String>)], order: &Option<Vec<usize>>, out: &mut Vec<u8>) {
    use crate::mergevcfs::{allele, format_qual};
    write!(out, "{}\t{}\t{}\t{}\t", cols[0], cols[1], cols[2], allele(cols[3])).unwrap();
    if cols[4] == "." {
        out.push(b'.');
    } else {
        out.extend_from_slice(cols[4].split(',').map(allele).collect::<Vec<_>>().join(",").as_bytes());
    }
    write!(out, "\t{}\t", format_qual(cols[5])).unwrap();
    match filter {
        Some(f) => out.extend_from_slice(f.as_bytes()),
        None => {
            if cols[6] == "." || cols[6] == "PASS" {
                out.extend_from_slice(cols[6].as_bytes());
            } else {
                let mut v: Vec<&str> = cols[6].split(';').collect();
                v.sort_unstable();
                v.dedup();
                out.extend_from_slice(v.join(";").as_bytes());
            }
        }
    }
    out.push(b'\t');
    let mut m: BTreeMap<&str, Option<&str>> = BTreeMap::new();
    if cols[7] != "." && !cols[7].is_empty() {
        for kv in cols[7].split(';').filter(|s| !s.is_empty()) {
            match kv.split_once('=') {
                Some((k, v)) if !v.is_empty() => m.insert(k, Some(v)),
                Some((k, _)) => m.insert(k, None),
                None => m.insert(kv, None),
            };
        }
    }
    for (k, v) in add {
        m.insert(k, v.as_deref());
    }
    if m.is_empty() {
        out.push(b'.');
    } else {
        let s: Vec<String> = m.iter().map(|(k, v)| match v {
            Some(v) => format!("{k}={v}"),
            None => k.to_string(),
        }).collect();
        out.extend_from_slice(s.join(";").as_bytes());
    }
    if cols.len() > 8 {
        out.push(b'\t');
        out.extend_from_slice(cols[8].as_bytes());
        match order {
            None => {
                for s in &cols[9..] {
                    out.push(b'\t');
                    out.extend_from_slice(s.as_bytes());
                }
            }
            Some(o) => {
                for &k in o {
                    out.push(b'\t');
                    out.extend_from_slice(cols[9 + k].as_bytes());
                }
            }
        }
    }
    out.push(b'\n');
}

pub fn run(argv: &[String]) {
    let a = parse_args(argv);
    let t0 = std::time::Instant::now();
    // tranches at or above the filter level, worst first
    let mut tranches: Vec<Tranche> = Vec::new();
    let mut extra: Vec<String> = vec![
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">".into(),
        "##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=\"Log odds of being a true variant versus being false under the trained gaussian mixture model\">".into(),
        "##INFO=<ID=culprit,Number=1,Type=String,Description=\"The annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out\">".into(),
        "##INFO=<ID=POSITIVE_TRAIN_SITE,Number=0,Type=Flag,Description=\"This variant was used to build the positive training set of good variants\">".into(),
        "##INFO=<ID=NEGATIVE_TRAIN_SITE,Number=0,Type=Flag,Description=\"This variant was used to build the negative training set of bad variants\">".into(),
        "##FILTER=<ID=PASS,Description=\"Site contains at least one allele that passes filters\">".into(),
    ];
    let mut cutoff = a.lod_cutoff;
    let d2s = crate::mutect::fmc_math::java_double_to_string;
    if let Some(ts) = a.ts_filter_level {
        let path = a.tranches_file.clone().unwrap_or_else(|| refuse("--truth-sensitivity-filter-level needs --tranches-file"));
        tranches = read_tranches(&path).unwrap_or_else(|e| refuse(&e)).into_iter().filter(|t| t.target >= ts).collect();
        tranches.reverse();
        if tranches.is_empty() {
            refuse(&format!("No tranches were found in the file or were above the truth sensitivity filter level {}", d2s(ts)));
        }
        for i in 0..tranches.len().saturating_sub(1) {
            let t = &tranches[i];
            extra.push(filter_line(
                &t.filter_name,
                &format!("Truth sensitivity tranche level for {} model at VQS Lod: {} <= x < {}", t.model, d2s(t.min_vqslod), d2s(tranches[i + 1].min_vqslod)),
            ));
        }
        extra.push(filter_line(&format!("{}+", tranches[0].filter_name), &format!("Truth sensitivity tranche level for {} model at VQS Lod < {}", tranches[0].model, d2s(tranches[0].min_vqslod))));
    } else {
        let c = cutoff.unwrap_or(0.0);
        cutoff = Some(c);
        extra.push(filter_line("LOW_VQSLOD", &format!("VQSLOD < {}", d2s(c))));
    }
    extra.push("##source=ApplyVQSR".into());
    let cmdline: Vec<String> = std::env::args().collect();
    extra.push(format!(
        "##variant-rsCommandLine=<ID=ApplyVQSR,CommandLine=\"{}\",Version=\"{}\",Date=\"{}\">",
        cmdline.join(" ").replace('"', "\\\""),
        env!("CARGO_PKG_VERSION"),
        crate::now_string()
    ));

    let mut lines = Lines::open(&a.input);
    let mut line = Vec::new();
    let mut meta = Vec::new();
    let chrom_line;
    loop {
        if !lines.next_line(&mut line) {
            refuse(&format!("{}: no #CHROM line", a.input));
        }
        let l = String::from_utf8_lossy(&line).into_owned();
        if l.starts_with("##") {
            meta.push(crate::vcfheader::repair_standard_line(l));
        } else {
            chrom_line = l;
            break;
        }
    }
    // samples in sorted order (htsjdk writes the header's TreeSet order)
    let samples: Vec<String> = chrom_line.split('\t').skip(9).map(|s| s.to_string()).collect();
    let mut sorted = samples.clone();
    sorted.sort();
    let order: Option<Vec<usize>> = if sorted == samples { None } else { Some(sorted.iter().map(|s| samples.iter().position(|x| x == s).unwrap()).collect()) };
    let hdr = header(&meta, extra, &sorted);
    let level = if crate::output::block_compressed(&a.output) { Some(a.compression_level) } else { None };
    let mut of = crate::output::OutFile::create(&a.output, &hdr, level).unwrap_or_else(|e| refuse(&format!("{}: {e}", a.output)));

    let mut recal = RecalReader::open(&a.recal);
    let ignore: HashSet<&str> = a.ignore_filters.iter().map(|s| s.as_str()).collect();
    let mut buf: Vec<u8> = Vec::with_capacity(8 << 20);
    let (mut n_in, mut n_eval, mut n_pass) = (0u64, 0u64, 0u64);
    while lines.next_line(&mut line) {
        if line.is_empty() || line[0] == b'#' {
            continue;
        }
        n_in += 1;
        let text = std::str::from_utf8(&line).unwrap_or_else(|_| refuse(&format!("{}: not UTF-8", a.input)));
        let cols: Vec<&str> = text.split('\t').collect();
        if cols.len() < 8 {
            refuse(&format!("{}: VCF line with fewer than 8 columns", a.input));
        }
        let vtype = vcfio::variant_type(cols[3].as_bytes(), cols[4].as_bytes());
        let evaluate = a.mode.accepts(vtype);
        let not_filtered = a.ignore_all_filters
            || vcfio::filters(cols[6].as_bytes()).is_none()
            || (!ignore.is_empty() && cols[6].split(';').all(|f| ignore.contains(f)));
        if !(evaluate && not_filtered) {
            reencode(&cols, None, &[], &order, &mut buf);
        } else {
            n_eval += 1;
            let start: i64 = cols[1].parse().unwrap_or_else(|_| refuse(&format!("{}: bad POS", a.input)));
            let end = match vcfio::info_double(cols[7].as_bytes(), "END") {
                Ok(e) if !e.is_nan() => e as i64,
                _ => start + cols[3].len() as i64 - 1,
            };
            let tid = *recal.rank.get(cols[0]).unwrap_or(&usize::MAX);
            let recs = if tid == usize::MAX { vec![] } else { recal.at(tid, start) };
            let Some(r) = recs.into_iter().find(|r| r.end == end) else {
                refuse(&format!(
                    "Encountered input variant which isn't found in the input recal file. Please make sure VariantRecalibrator and ApplyVQSR were run on the same set of input variants. First seen at: {}:{start}",
                    cols[0]
                ))
            };
            if r.lod.is_empty() {
                refuse(&format!("Encountered a malformed record in the input recal file. There is no lod for the record at: {}:{start}", cols[0]));
            }
            let lod: f64 = vcfio::parse_vcf_double(&r.lod)
                .filter(|_| !r.lod.to_ascii_lowercase().contains("inf") || r.lod.contains("Infinity"))
                .unwrap_or_else(|| refuse(&format!("Encountered a malformed record in the input recal file. The lod is unreadable for the record at: {}:{start}", cols[0])));
            let mut add: Vec<(&str, Option<String>)> = vec![("VQSLOD", Some(hc_genotype::javamath::format_vcf_double(lod))), ("culprit", r.culprit.clone())];
            if r.positive {
                add.push(("POSITIVE_TRAIN_SITE", None));
            }
            if r.negative {
                add.push(("NEGATIVE_TRAIN_SITE", None));
            }
            let f = filter_for(lod, &tranches, if tranches.is_empty() { cutoff } else { None });
            if f == "PASS" {
                n_pass += 1;
            }
            if a.exclude_filtered && f != "PASS" && f != "." {
                continue;
            }
            reencode(&cols, Some(&f), &add, &order, &mut buf);
        }
        if buf.len() >= 8 << 20 {
            of.add_text(&buf).unwrap();
            buf.clear();
        }
    }
    of.add_text(&buf).unwrap();
    of.finish().unwrap_or_else(|e| refuse(&format!("{}: {e}", a.output)));
    eprintln!("[applyvqsr] {n_in} records, {n_eval} {} records recalibrated, {n_pass} PASS, in {:.1}s", a.mode.name(), t0.elapsed().as_secs_f64());
}
