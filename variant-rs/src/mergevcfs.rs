//! MergeVcfs (`variant-rs mergevcfs`): Picard `MergeVcfs`. Inputs with the same samples, merged by
//! `VariantContextComparator` (contig index in the first header, then start) through htsjdk's `MergingIterator`,
//! i.e. a `java.util.PriorityQueue` of the inputs keyed by their next record: equal positions in several inputs come
//! out in the order that heap gives (simulated here operation for operation). Header: `VCFUtils.smartMergeHeaders`
//! over each input's lines in sorted order, first line per key winning (`KEY-ID` for INFO / FILTER / FORMAT / contig /
//! ALT / PEDIGREE / META / SAMPLE lines, the bare key for the others, so one `##GATKCommandLine` survives), written
//! sorted as htsjdk's VCF writer does. Records are decoded and re-encoded by htsjdk (`VCFEncoder`): genotype columns as
//! read, INFO keys sorted (`KEY=` -> `KEY`), QUAL `%.2f` without `.00`, FILTER values sorted, alleles upper-cased;
//! that re-encoding is done here.
//! Refused (rather than guessed): INFO / FORMAT / FILTER header lines htsjdk would write differently, and
//! INFO / FORMAT definitions that conflict between inputs (Picard promotes them).

use crate::output::OutFile;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};

const USAGE: &str = "usage: variant-rs mergevcfs -I a.vcf.gz [-I b.vcf.gz ...] -O out.vcf.gz [--CREATE_INDEX true] [--compression-level 2] [--SEQUENCE_DICTIONARY ref.dict]";

fn bool_arg(v: &str) -> bool {
    match v.to_ascii_lowercase().as_str() {
        "true" | "t" | "yes" | "y" => true,
        "false" | "f" | "no" | "n" => false,
        _ => panic!("not a boolean: {v}"),
    }
}

fn refuse(why: String) -> ! {
    crate::gatk_cli::refuse(&why)
}

pub(crate) fn open_lines(path: &str) -> Box<dyn Iterator<Item = String>> {
    let f = std::fs::File::open(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let r: Box<dyn std::io::Read> = if path.ends_with(".gz") || path.ends_with(".bgz") { Box::new(flate2::read::MultiGzDecoder::new(f)) } else { Box::new(f) };
    let p = path.to_string();
    Box::new(BufReader::with_capacity(1 << 20, r).lines().map(move |l| l.unwrap_or_else(|e| panic!("{p}: {e}"))))
}

/// A `##` header line as htsjdk keys it in `smartMergeHeaders`.
fn merge_key(l: &str) -> String {
    let body = &l[2..];
    let key = body.split('=').next().unwrap_or(body);
    const ID_LINES: [&str; 8] = ["INFO", "FILTER", "FORMAT", "contig", "ALT", "PEDIGREE", "META", "SAMPLE"];
    if ID_LINES.contains(&key) {
        if let Some(rest) = body.strip_prefix(&format!("{key}=<ID=")) {
            let id = rest.split([',', '>']).next().unwrap_or("");
            return format!("{key}-{id}");
        }
    }
    key.to_string()
}

/// INFO / FORMAT / FILTER lines must be in the form htsjdk writes them (`toStringEncoding`), as they are re-written
/// from their parsed fields: (ID, Number, Type) for INFO / FORMAT.
fn check_canonical(l: &str, path: &str) -> Option<(String, String)> {
    let bad = || refuse(format!("{path}: header line {l}: not in the form htsjdk writes; not supported"));
    for kind in ["INFO", "FORMAT"] {
        if let Some(rest) = l.strip_prefix(&format!("##{kind}=<ID=")) {
            let (_id, rest) = rest.split_once(",Number=").unwrap_or_else(|| bad());
            let (num, rest) = rest.split_once(",Type=").unwrap_or_else(|| bad());
            let (typ, rest) = rest.split_once(",Description=\"").unwrap_or_else(|| bad());
            // Description up to the closing quote, then optional Source / Version, then '>'
            if !rest.ends_with('>') || !(rest.ends_with("\">") || rest.contains("\",Source=") || rest.contains("\",Version=")) {
                bad();
            }
            if !["Integer", "Float", "Flag", "Character", "String"].contains(&typ) {
                bad();
            }
            return Some((num.to_string(), typ.to_string()));
        }
    }
    if let Some(rest) = l.strip_prefix("##FILTER=<ID=") {
        if !(rest.contains(",Description=\"") && rest.ends_with("\">")) {
            bad();
        }
    }
    None
}

/// QUAL as htsjdk writes a decoded one: `parseQual` keeps log10PError = QUAL / -10.0 ('.' and the VCF 3 missing value
/// -1, within 5e-5, are missing: NO_LOG10_PERROR = 1.0, which also a QUAL of -10 gives), the writer prints '.' for
/// 1.0, else `formatQualValue(log10PError * -10 + 0.0)`: `%.2f`, a trailing `.00` dropped.
pub(crate) fn format_qual(q: &str) -> String {
    if q == "." {
        return q.to_string();
    }
    let v: f64 = q.parse().unwrap_or_else(|_| panic!("QUAL {q}"));
    if v < 0.0 && (v - -1.0).abs() < 5.0e-5 {
        return ".".to_string();
    }
    let log10_p_error = v / -10.0;
    if log10_p_error == 1.0 {
        return ".".to_string();
    }
    let s = hc_genotype::javamath::format_fixed(log10_p_error * -10.0 + 0.0, 2);
    s.strip_suffix(".00").map(|x| x.to_string()).unwrap_or(s)
}

/// `Allele.getDisplayString` after `Allele.create`: bases upper-cased, symbolic alleles as they are.
pub(crate) fn allele(a: &str) -> String {
    let b = a.as_bytes();
    let symbolic = a.len() > 1 && (b[0] == b'<' || b[b.len() - 1] == b'>' || a.contains('[') || a.contains(']') || b[0] == b'.' || b[b.len() - 1] == b'.');
    if symbolic || a == "*" || a == "." {
        a.to_string()
    } else {
        a.to_ascii_uppercase()
    }
}

/// One record line re-encoded as htsjdk's `VCFEncoder` writes it after `VCFCodec` decoded it.
fn reencode(line: &str, declared: &HashMap<String, bool>) -> String {
    let mut f: Vec<&str> = line.splitn(9, '\t').collect();
    assert!(f.len() >= 8, "VCF record with fewer than 8 columns: {line}");
    let mut out = String::with_capacity(line.len() + 8);
    out.push_str(f[0]);
    out.push('\t');
    out.push_str(f[1]);
    out.push('\t');
    out.push_str(f[2]);
    out.push('\t');
    out.push_str(&allele(f[3]));
    out.push('\t');
    if f[4] == "." {
        out.push('.');
    } else {
        out.push_str(&f[4].split(',').map(allele).collect::<Vec<_>>().join(","));
    }
    out.push('\t');
    out.push_str(&format_qual(f[5]));
    out.push('\t');
    // FILTER: sorted when filtered; PASS / '.' as they are
    if f[6] == "." || f[6] == "PASS" {
        out.push_str(f[6]);
    } else {
        let mut v: Vec<&str> = f[6].split(';').collect();
        v.sort_unstable();
        v.dedup();
        out.push_str(&v.join(";"));
    }
    out.push('\t');
    // INFO as htsjdk decodes and writes it (keys sorted; an undeclared key fails, as in Picard)
    let mut info = Vec::new();
    let at = format!("{}:{}", f[0], f[1]);
    crate::vcfinfo::write(&crate::vcfinfo::parse(f[7], declared), declared, &at, &mut info).unwrap_or_else(|e| {
        eprintln!("variant-rs MergeVcfs: {e}");
        std::process::exit(1)
    });
    out.push_str(std::str::from_utf8(&info).unwrap());
    if f.len() == 9 {
        out.push('\t');
        out.push_str(f.remove(8));
    }
    out
}

/// A record waiting in the merge: its comparator keys and line.
struct Head {
    contig: usize,
    start: i64,
    line: String,
}

/// One input: its remaining record lines.
struct Input {
    lines: Box<dyn Iterator<Item = String>>,
    path: String,
}

/// `VariantContextComparator`
fn cmp_heads(a: &Head, b: &Head) -> Ordering {
    a.contig.cmp(&b.contig).then(a.start.cmp(&b.start))
}

/// `java.util.PriorityQueue` (JDK siftUp / siftDown) of input indices keyed by their head record.
struct JavaPq {
    q: Vec<usize>,
}

impl JavaPq {
    fn offer(&mut self, x: usize, heads: &[Option<Head>]) {
        let key = |i: usize| heads[i].as_ref().unwrap();
        let mut k = self.q.len();
        self.q.push(x);
        while k > 0 {
            let parent = (k - 1) >> 1;
            let e = self.q[parent];
            if cmp_heads(key(x), key(e)) != Ordering::Less {
                break;
            }
            self.q[k] = e;
            k = parent;
        }
        self.q[k] = x;
    }
    fn poll(&mut self, heads: &[Option<Head>]) -> Option<usize> {
        let key = |i: usize| heads[i].as_ref().unwrap();
        if self.q.is_empty() {
            return None;
        }
        let result = self.q[0];
        let x = self.q.pop().unwrap();
        let n = self.q.len();
        if n > 0 {
            let mut k = 0;
            let half = n >> 1;
            while k < half {
                let mut child = 2 * k + 1;
                let mut c = self.q[child];
                let right = child + 1;
                if right < n && cmp_heads(key(c), key(self.q[right])) == Ordering::Greater {
                    child = right;
                    c = self.q[child];
                }
                if cmp_heads(key(x), key(c)) != Ordering::Greater {
                    break;
                }
                self.q[k] = c;
                k = child;
            }
            self.q[k] = x;
        }
        Some(result)
    }
}

pub fn run(argv: &[String]) {
    let t0 = std::time::Instant::now();
    let (mut paths, mut out_path, mut index, mut level) = (Vec::<String>::new(), String::new(), true, 2i32);
    let mut seq_dict: Option<String> = None;
    let mut i = 0;
    while i < argv.len() {
        let v = || argv.get(i + 1).cloned().unwrap_or_else(|| panic!("{} needs a value\n{USAGE}", argv[i]));
        match argv[i].as_str() {
            "-I" | "--INPUT" => paths.push(v()),
            "-O" | "--OUTPUT" => out_path = v(),
            "--CREATE_INDEX" => index = bool_arg(&v()),
            "--compression-level" => level = v().parse().expect("--compression-level"),
            "--SEQUENCE_DICTIONARY" => seq_dict = Some(v()),
            x => panic!("unknown argument {x}\n{USAGE}"),
        }
        i += 2;
    }
    assert!(!paths.is_empty() && !out_path.is_empty(), "{USAGE}");
    // Picard's INPUT list files (IOUtil.unrollPaths with FileExtensions.VCF_LIST)
    let paths = crate::gatherbam::unroll_inputs(&paths, &[".vcf", ".vcf.gz", ".vcf.bgz", ".bcf"]);
    // SEQUENCE_DICTIONARY: Picard gives it to inputs without ##contig lines (refused below) and to the writer, which
    // only puts it in the index (the .tbi still lists the records' contigs): the output does not change. Picard opens
    // it with SamReaderFactory (a .dict / SAM / BAM / CRAM header).
    if let Some(d) = &seq_dict {
        if ![".dict", ".sam", ".bam", ".cram"].iter().any(|e| d.ends_with(e)) {
            refuse(format!("SEQUENCE_DICTIONARY {d}: only a .dict / .sam / .bam / .cram file is supported"));
        }
        if let Err(e) = std::fs::File::open(d) {
            eprintln!("SEQUENCE_DICTIONARY {d}: {e}");
            std::process::exit(1);
        }
    }
    let gz = crate::output::block_compressed(&out_path);
    if index && !gz {
        refuse(format!("{out_path}: an uncompressed VCF with CREATE_INDEX=true needs a Tribble .idx index; not supported (use .vcf.gz or CREATE_INDEX=false)"));
    }

    // headers
    let mut inputs: Vec<Input> = Vec::new();
    let mut merged: Vec<(String, String)> = Vec::new(); // (key, line), first per key
    let mut defs: HashMap<String, (String, String)> = HashMap::new();
    let mut samples_line: Option<String> = None;
    let mut contig_index: HashMap<String, usize> = HashMap::new();
    for p in &paths {
        let mut lines = open_lines(p);
        let mut hdr: Vec<String> = Vec::new();
        let chrom = loop {
            let l = lines.next().unwrap_or_else(|| panic!("{p}: no #CHROM line"));
            if l.starts_with("##") {
                hdr.push(crate::vcfheader::repair_standard_line(l));
            } else if l.starts_with("#CHROM") {
                break l;
            } else {
                panic!("{p}: record before the #CHROM line");
            }
        };
        // the sample columns must match the first input's
        let samples: Vec<&str> = chrom.split('\t').skip(9).collect();
        crate::vcfheader::require_sorted_samples(&samples.iter().map(|s| s.to_string()).collect::<Vec<_>>(), p);
        match &samples_line {
            None => samples_line = Some(chrom.clone()),
            Some(s) => {
                if s.split('\t').skip(9).collect::<Vec<_>>() != samples {
                    panic!("Input path {p} has sample entries that don't match the other files.");
                }
            }
        }
        let contigs: Vec<String> = hdr.iter().filter_map(|l| l.strip_prefix("##contig=<ID=")).map(|r| r.split([',', '>']).next().unwrap().to_string()).collect();
        if contigs.is_empty() {
            refuse(format!("{p}: no ##contig lines (Picard fills them from SEQUENCE_DICTIONARY): not supported"));
        }
        if contig_index.is_empty() {
            contig_index = contigs.iter().enumerate().map(|(k, c)| (c.clone(), k)).collect();
        } else if contigs.iter().enumerate().any(|(k, c)| contig_index.get(c) != Some(&k)) {
            panic!("The contig entries in input path {p} are not compatible with the others.");
        }
        // getMetaDataInSortedOrder: by text, contig lines by their index (in place of the first contig line)
        let mut sorted: Vec<&String> = hdr.iter().filter(|l| !l.starts_with("##fileformat=")).collect();
        let sort_key = |l: &String| if l.starts_with("##contig=<") { "##contig=<".to_string() } else { l.to_string() };
        sorted.sort_by_key(|l| sort_key(l));
        for l in sorted {
            let def = check_canonical(l, p);
            let k = merge_key(l);
            match merged.iter().find(|(mk, _)| *mk == k) {
                None => {
                    if let Some(d) = def {
                        defs.insert(k.clone(), d);
                    }
                    merged.push((k, l.clone()));
                }
                Some((_, first)) if first == l => {}
                Some(_) => {
                    if let (Some(d), Some(d0)) = (def, defs.get(&k)) {
                        if &d != d0 {
                            refuse(format!("{p}: header line {l} conflicts with another input's definition (Picard promotes it); not supported"));
                        }
                    }
                }
            }
        }
        inputs.push(Input { lines, path: p.clone() });
    }
    // the writer: ##fileformat first, then the merged lines sorted (contig lines in their order)
    let declared = crate::vcfinfo::declared(&merged.iter().map(|(_, l)| l.clone()).collect::<Vec<_>>());
    let mut lines: Vec<&String> = merged.iter().map(|(_, l)| l).collect();
    let sort_key = |l: &String| if l.starts_with("##contig=<") { "##contig=<".to_string() } else { l.to_string() };
    lines.sort_by_key(|l| sort_key(l));
    let mut header = String::from("##fileformat=VCFv4.2\n");
    for l in lines {
        header.push_str(l);
        header.push('\n');
    }
    header.push_str(samples_line.as_deref().unwrap());
    header.push('\n');
    let mut out = OutFile::create(&out_path, &header, if gz { Some(level) } else { None }).unwrap_or_else(|e| panic!("{e}"));

    // merge
    let next_head = |inp: &mut Input| -> Option<Head> {
        let line = inp.lines.next()?;
        let mut it = line.splitn(3, '\t');
        let c = it.next().unwrap();
        let start: i64 = it.next().and_then(|s| s.parse().ok()).unwrap_or_else(|| panic!("{}: bad record {line}", inp.path));
        let contig = *contig_index.get(c).unwrap_or_else(|| panic!("{}: contig {c} is not in the header", inp.path));
        Some(Head { contig, start, line })
    };
    let mut heads: Vec<Option<Head>> = inputs.iter_mut().map(|inp| next_head(inp)).collect();
    let mut pq = JavaPq { q: Vec::new() };
    for k in 0..heads.len() {
        if heads[k].is_some() {
            pq.offer(k, &heads);
        }
    }
    let mut buf = String::with_capacity(8 << 20);
    let mut last: Option<(usize, i64)> = None;
    let mut n = 0u64;
    while let Some(k) = pq.poll(&heads) {
        let h = heads[k].take().unwrap();
        if let Some(l) = last {
            if (h.contig, h.start) < l {
                panic!("The elements of the input Iterators are not sorted according to the comparator htsjdk.variant.variantcontext.VariantContextComparator");
            }
        }
        last = Some((h.contig, h.start));
        buf.push_str(&reencode(&h.line, &declared));
        buf.push('\n');
        n += 1;
        if buf.len() >= 4 << 20 {
            out.add_text(buf.as_bytes()).unwrap();
            buf.clear();
        }
        heads[k] = next_head(&mut inputs[k]);
        if heads[k].is_some() {
            pq.offer(k, &heads);
        }
    }
    out.add_text(buf.as_bytes()).unwrap();
    if !index && gz {
        // no .tbi
        out.finish().unwrap();
        let _ = std::fs::remove_file(format!("{out_path}.tbi"));
    } else {
        out.finish().unwrap();
    }
    eprintln!("mergevcfs: {} inputs, {n} records -> {out_path} ({:.1}s)", paths.len(), t0.elapsed().as_secs_f64());
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn encoder() {
        let d: HashMap<String, bool> = [("Z", false), ("A", true), ("K", false)].iter().map(|(k, v)| (k.to_string(), *v)).collect();
        assert_eq!(reencode("c\t5\t.\tac\tG,<NON_REF>\t30.5\tb;a\tZ=1;A;K=\tGT\t0/1", &d), "c\t5\t.\tAC\tG,<NON_REF>\t30.50\ta;b\tA;K=.;Z=1\tGT\t0/1");
        assert_eq!(reencode("c\t5\t.\tA\t.\t30.00\t.\t.", &d), "c\t5\t.\tA\t.\t30\t.\t.");
        assert_eq!(merge_key("##GATKCommandLine=<ID=HaplotypeCaller,CommandLine=\"x\">"), "GATKCommandLine");
        assert_eq!(merge_key("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">"), "INFO-DP");
        assert_eq!(merge_key("##contig=<ID=chr1,length=5>"), "contig-chr1");
    }
}
