//! `variant-rs getpileupsummaries` and `variant-rs calculatecontamination`: ports of GATK 4.7 GetPileupSummaries
//! (LocusWalker over the common-SNP sites; default read filters, no downsampling) and CalculateContamination
//! (ContaminationModel with the kernel segmenter, tumor-only or with a matched normal, optional
//! `--tumor-segmentation` output). Tables are written as GATK's TableWriter does (Double.toString values).

use super::cm_linalg::{java_compare, jmax, jmin, mean, percentile, persistence_minima, svd};
use super::fmc_math::{binomial_probability, brent_max_pv, fexp, java_double_to_string, jlog};
use super::somatic::{java_string_hash, stream_sum};
use crate::active::{good_cigar, read_length, ref_length, OP_D, OP_EQ, OP_I, OP_M, OP_N, OP_S, OP_X};
use crate::bam::{self, Rec};
use crate::input::Input;
use std::cmp::Ordering;
use std::io::{BufRead, BufReader, Write};
use std::sync::atomic::AtomicU64;
use std::sync::Arc;


// ------------------------------------------------------------------------------------------------
// PileupSummary tables

#[derive(Clone, Debug)]
pub struct Site {
    pub contig: String,
    pub pos: i64,
    pub refc: i32,
    pub altc: i32,
    pub other: i32,
    pub total: i32,
    pub af: f64,
}

impl Site {
    fn alt_fraction(&self) -> f64 {
        if self.total == 0 { 0.0 } else { self.altc as f64 / self.total as f64 }
    }
    fn ref_frequency(&self) -> f64 {
        1.0 - self.af
    }
}

/// A text file, gzip-compressed or not (by its content, as htsjdk).
pub(crate) fn open_text(path: &str) -> Box<dyn BufRead> {
    let mut magic = [0u8; 2];
    let gz = std::fs::File::open(path).and_then(|mut f| std::io::Read::read_exact(&mut f, &mut magic)).is_ok() && magic == [0x1f, 0x8b];
    let f = std::fs::File::open(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    if gz {
        Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(f)))
    } else {
        Box::new(BufReader::new(f))
    }
}

/// GATK TableReader: `#` comment lines (metadata `#<METADATA>key=value`), then a header line, then rows.
fn read_table(path: &str) -> (Vec<(String, String)>, Vec<String>, Vec<Vec<String>>) {
    let mut meta = Vec::new();
    let mut header: Option<Vec<String>> = None;
    let mut rows = Vec::new();
    for line in open_text(path).lines() {
        let line = line.unwrap();
        if let Some(c) = line.strip_prefix('#') {
            if let Some(kv) = c.strip_prefix("<METADATA>") {
                if let Some((k, v)) = kv.split_once('=') {
                    meta.push((k.to_string(), v.to_string()));
                }
            }
            continue;
        }
        if line.is_empty() {
            continue;
        }
        let f: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        if header.is_none() {
            header = Some(f);
        } else {
            rows.push(f);
        }
    }
    (meta, header.unwrap_or_default(), rows)
}

fn col(header: &[String], name: &str, path: &str) -> usize {
    header.iter().position(|h| h == name).unwrap_or_else(|| panic!("{path}: missing column {name}"))
}

pub(crate) fn read_pileups(path: &str) -> (String, Vec<Site>) {
    let (meta, h, rows) = read_table(path);
    let sample = meta.iter().find(|(k, _)| k == "SAMPLE").map(|(_, v)| v.clone()).unwrap_or_default();
    let (c, p, r, a, o, f) = (
        col(&h, "contig", path),
        col(&h, "position", path),
        col(&h, "ref_count", path),
        col(&h, "alt_count", path),
        col(&h, "other_alt_count", path),
        col(&h, "allele_frequency", path),
    );
    let sites = rows
        .iter()
        .map(|row| {
            let (refc, altc, other): (i32, i32, i32) = (row[r].parse().unwrap(), row[a].parse().unwrap(), row[o].parse().unwrap());
            Site {
                contig: row[c].clone(),
                pos: row[p].parse().unwrap(),
                refc,
                altc,
                other,
                total: refc + altc + other,
                af: row[f].parse().unwrap(),
            }
        })
        .collect();
    (sample, sites)
}

pub(crate) fn write_pileups(path: &str, sample: &str, sites: &[Site]) {
    let mut w = std::io::BufWriter::new(std::fs::File::create(path).unwrap_or_else(|e| panic!("{path}: {e}")));
    writeln!(w, "#<METADATA>SAMPLE={sample}").unwrap();
    writeln!(w, "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency").unwrap();
    for s in sites {
        writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}", s.contig, s.pos, s.refc, s.altc, s.other, java_double_to_string(s.af)).unwrap();
    }
}

/// `MinorAlleleFractionRecord.readFromFile`: (sample, [(contig, start, end, maf)])
pub fn read_segments(path: &str) -> (String, Vec<(String, i64, i64, f64)>) {
    let (meta, h, rows) = read_table(path);
    let sample = meta.iter().find(|(k, _)| k == "SAMPLE").map(|(_, v)| v.clone()).unwrap_or_default();
    let (c, s, e, m) = (col(&h, "contig", path), col(&h, "start", path), col(&h, "end", path), col(&h, "minor_allele_fraction", path));
    (sample, rows.iter().map(|r| (r[c].clone(), r[s].parse().unwrap(), r[e].parse().unwrap(), r[m].parse().unwrap())).collect())
}

/// `ContaminationRecord.readFromFile`: [(sample, contamination, error)]
pub fn read_contamination(path: &str) -> Vec<(String, f64, f64)> {
    let (_, h, rows) = read_table(path);
    let (s, c, e) = (col(&h, "sample", path), col(&h, "contamination", path), col(&h, "error", path));
    rows.iter().map(|r| (r[s].clone(), r[c].parse().unwrap(), r[e].parse().unwrap())).collect()
}

// ------------------------------------------------------------------------------------------------
// GetPileupSummaries

const GPS_USAGE: &str = "usage: variant-rs getpileupsummaries -I reads.bam -V common.vcf[.gz] -L common.vcf[.gz]|intervals.bed|... -O pileups.table \
[-R ref.fa (CRAM)] [--minimum-population-allele-frequency 0.01] [--maximum-population-allele-frequency 0.2]";

/// One resource VCF record: start, end (inclusive), and Some((ref base, alt base, AF)) when it is a biallelic SNP.
struct VcfRec {
    start: i64,
    end: i64,
    snp: Option<(u8, u8, Option<f64>)>,
}

/// One data line of a sites VCF -> (contig, record).
fn parse_sites_line(line: &str) -> Option<(&str, VcfRec)> {
    if line.starts_with('#') || line.is_empty() {
        return None;
    }
    let f: Vec<&str> = line.splitn(9, '\t').collect();
    let (chrom, pos, refa, alt, info) = (f[0], f[1].parse::<i64>().unwrap(), f[3], f[4], f.get(7).copied().unwrap_or("."));
    let alts: Vec<&str> = if alt == "." { vec![] } else { alt.split(',').collect() };
    let end = pos + refa.len() as i64 - 1;
    let symbolic = |a: &str| a.starts_with('<') || a.contains('[') || a.contains(']') || a == "*" || a == ".";
    let snp = if alts.len() == 1 && refa.len() == 1 && alts[0].len() == 1 && !symbolic(alts[0]) && !refa.eq_ignore_ascii_case(alts[0]) {
        let af = info.split(';').find_map(|kv| kv.strip_prefix("AF=")).map(|v| v.split(',').next().unwrap().parse::<f64>().unwrap());
        Some((refa.as_bytes()[0], alts[0].as_bytes()[0], af))
    } else {
        None
    };
    Some((chrom, VcfRec { start: pos, end, snp }))
}

/// CHROM, POS and end (POS + |REF| - 1) of a sites VCF data line, without parsing the rest.
fn peek_sites_line(line: &str) -> Option<(&str, i64, i64)> {
    let b = line.as_bytes();
    if b.is_empty() || b[0] == b'#' {
        return None;
    }
    let mut it = memchr::memchr_iter(b'\t', b);
    let (t0, t1, t2, t3) = (it.next()?, it.next()?, it.next()?, it.next()?);
    let pos: i64 = line[t0 + 1..t1].parse().ok()?;
    Some((&line[..t0], pos, pos + (t3 - t2 - 1) as i64 - 1))
}

/// Whether the VCF header declares INFO/AF (reads the header only).
fn vcf_has_af(path: &str) -> bool {
    for line in open_text(path).lines() {
        let line = line.unwrap();
        if !line.starts_with('#') {
            break;
        }
        if line.starts_with("##INFO=<ID=AF,") {
            return true;
        }
    }
    false
}

/// Every record of a sites VCF, streamed in file order.
fn for_each_sites_rec(path: &str, mut f: impl FnMut(&str, VcfRec)) {
    for line in open_text(path).lines() {
        let line = line.unwrap();
        if let Some((c, r)) = parse_sites_line(&line) {
            f(c, r);
        }
    }
}

/// The sites GetPileupSummaries counts, chosen from the resource records in file order: the first record
/// overlapping a traversal locus decides it; biallelic SNPs with AF in range are kept.
struct SiteSel {
    min_af: f64,
    max_af: f64,
    /// the first record (in traversal order) without AF, for GATK's one warning
    first_noaf: Option<String>,
    saw_with: bool,
    saw_without: bool,
    /// per contig: pos, ref idx, alt idx, AF
    sites: Vec<Vec<(i64, usize, usize, f64)>>,
    max_end: Vec<i64>,
}

impl SiteSel {
    fn new(n: usize, min_af: f64, max_af: f64) -> SiteSel {
        SiteSel { min_af, max_af, first_noaf: None, saw_with: false, saw_without: false, sites: vec![Vec::new(); n], max_end: vec![i64::MIN; n] }
    }

    /// A record that starts outside the traversal: it only extends the span that decides "first record" later.
    fn note_end(&mut self, tid: usize, end: i64) {
        self.max_end[tid] = self.max_end[tid].max(end);
    }

    /// `ivs`: the contig's merged traversal intervals (None: every record is a locus).
    fn feed(&mut self, tid: usize, contig: &str, r: &VcfRec, ivs: Option<&[(i64, i64)]>) {
        let first = self.max_end[tid] < r.start;
        self.max_end[tid] = self.max_end[tid].max(r.end);
        if !first {
            return;
        }
        if let Some(ivs) = ivs {
            let k = ivs.partition_point(|iv| iv.1 < r.start);
            if !(k < ivs.len() && ivs[k].0 <= r.start) {
                return;
            }
        }
        if let Some((rb, ab, af)) = r.snp {
            match af {
                None => {
                    if self.first_noaf.is_none() {
                        self.first_noaf = Some(format!("{contig}:{}", r.start));
                    }
                    self.saw_without = true;
                }
                Some(af) => {
                    self.saw_with = true;
                    if self.min_af < af && af < self.max_af {
                        let (ri, ai) = (base_index(rb), base_index(ab));
                        self.sites[tid].push((r.start, ri.unwrap_or(usize::MAX), ai.unwrap_or(usize::MAX), af));
                    }
                }
            }
        }
    }
}

/// `BaseUtils.simpleBaseToBaseIndex`
fn base_index(b: u8) -> Option<usize> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// GetPileupSummaries.getDefaultReadFilters
fn passes_gps_filters(r: &Rec, ref_lens: &[i64]) -> bool {
    let unmapped = r.flag & 0x4 != 0 || r.ref_id < 0 || r.start < 1;
    r.mapq >= 50
        && r.mapq != 255
        && r.mapq != 0
        && !unmapped
        && r.flag & 0x100 == 0
        && r.flag & 0x800 == 0
        && r.flag & 0x400 == 0
        && r.flag & 0x200 == 0
        && ref_length(&r.cigar) > 0
        && (r.flag & 0x1 == 0 || r.flag & 0x8 != 0 || r.next_ref_id == r.ref_id)
        && good_cigar(&r.cigar)
        && r.has_rg
        && !r.seq.is_empty()
        && !r.qual_missing
        && r.qual.len() == r.seq.len()
        && read_length(&r.cigar) == r.seq.len()
        && !r.cigar.iter().any(|e| e.0 == OP_N)
        && (r.ref_id as usize) < ref_lens.len()
        && r.start <= ref_lens[r.ref_id as usize]
}

const DEFAULT_ADAPTOR_SIZE: i32 = 100;

/// `ReadUtils.getAdaptorBoundary` (None = CANNOT_COMPUTE_ADAPTOR_BOUNDARY).
fn adaptor_boundary(r: &Rec) -> Option<i64> {
    let paired = r.flag & 0x1 != 0;
    let unmapped = r.flag & 0x4 != 0 || r.ref_id < 0 || r.start < 1;
    if r.tlen == 0 || !paired || unmapped || r.flag & 0x8 != 0 {
        return None;
    }
    let rev = r.flag & 0x10 != 0;
    if rev == (r.flag & 0x20 != 0) {
        return None;
    }
    let end = r.start + ref_length(&r.cigar) - 1;
    let well = if rev { end > r.mate_start } else { r.start <= r.mate_start + r.tlen as i64 };
    if !well {
        return None;
    }
    Some(if rev { r.mate_start - 1 } else { r.start + (r.tlen as i64).abs() })
}

/// Merges 1-based inclusive intervals (overlapping or abutting, as IntervalMergingRule.ALL).
fn merge_intervals(mut v: Vec<(i64, i64)>) -> Vec<(i64, i64)> {
    v.sort();
    let mut out: Vec<(i64, i64)> = Vec::new();
    for (s, e) in v {
        match out.last_mut() {
            Some(l) if s <= l.1 + 1 => l.1 = l.1.max(e),
            _ => out.push((s, e)),
        }
    }
    out
}

pub fn run_gps(argv: &[String]) {
    let t0 = std::time::Instant::now();
    let (mut input, mut vcf, mut intervals, mut output, mut reference) = (None::<String>, None::<String>, None::<String>, None::<String>, String::new());
    let (mut min_af, mut max_af) = (0.01f64, 0.2f64);
    let mut i = 0;
    while i < argv.len() {
        let v = || argv.get(i + 1).cloned().unwrap_or_else(|| panic!("{GPS_USAGE}"));
        match argv[i].as_str() {
            "-I" | "--input" => input = Some(v()),
            "-V" | "--variant" => vcf = Some(v()),
            "-L" | "--intervals" => intervals = Some(v()),
            "-O" | "--output" => output = Some(v()),
            "-R" | "--reference" => reference = v(),
            "-min-af" | "--minimum-population-allele-frequency" => min_af = v().parse().unwrap(),
            "-max-af" | "--maximum-population-allele-frequency" => max_af = v().parse().unwrap(),
            a => panic!("unknown argument {a}\n{GPS_USAGE}"),
        }
        i += 2;
    }
    let (input, vcf, output) = (input.expect(GPS_USAGE), vcf.expect(GPS_USAGE), output.expect(GPS_USAGE));
    let intervals = intervals.expect("GetPileupSummaries requires intervals (-L)");

    if !vcf_has_af(&vcf) {
        panic!("Population vcf does not have an allele frequency (AF) info field in its header.");
    }
    let is_vcf = |p: &str| [".vcf", ".vcf.gz", ".vcf.bgz"].iter().any(|e| p.ends_with(e));
    // -L is the resource itself: every record is a locus
    let walk_all = is_vcf(&intervals) && intervals == vcf;
    let (inp, header) = Input::open(&input, &reference);
    // ReadUtils.getSamplesFromHeader is a TreeSet: the smallest sample name
    let sample = header
        .text
        .lines()
        .filter(|l| l.starts_with("@RG"))
        .filter_map(|l| l.split('\t').find_map(|f| f.strip_prefix("SM:")).map(|s| s.to_string()))
        .min()
        .expect("no sample in the read groups");
    let tid_of = |c: &str| {
        header
            .ref_names
            .iter()
            .position(|n| n == c)
            .unwrap_or_else(|| panic!("Badly formed genome unclippedLoc: Contig {c} given as location, but this contig isn't present in the sequence dictionary"))
    };
    let n_ref = header.ref_names.len();
    let mut sel = SiteSel::new(n_ref, min_af, max_af);
    if walk_all {
        // the whole resource, streamed (it can be a genome-wide gnomAD: never held in memory)
        let mut last: Option<(String, usize)> = None;
        for_each_sites_rec(&vcf, |c, r| {
            if last.as_ref().map_or(true, |l| l.0 != c) {
                last = Some((c.to_string(), tid_of(c)));
            }
            sel.feed(last.as_ref().unwrap().1, c, &r, None);
        });
    } else {
        // traversal intervals per contig: another sites VCF's records, or an interval list / BED / contig:start-end
        let mut trav: Vec<Vec<(i64, i64)>> = vec![Vec::new(); n_ref];
        if is_vcf(&intervals) {
            for_each_sites_rec(&intervals, |c, r| trav[tid_of(c)].push((r.start, r.end)));
        } else {
            let dict = crate::intervals::Dict::new(header.ref_names.clone(), header.ref_lens.clone());
            for iv in crate::intervals::parse_l_arg(&intervals, &dict).unwrap_or_else(|e| panic!("{e}")) {
                trav[iv.tid].push((iv.start, iv.end));
            }
        }
        let trav: Vec<Vec<(i64, i64)>> = trav.into_iter().map(merge_intervals).collect();
        if let Some(tl) = super::germline::TabixLines::open(&vcf) {
            // through the index, as GATK queries it: only the blocks over the intervals are read. Contigs in
            // parallel (each its own selection state); a line starting outside the intervals is only peeked at
            // (it can still decide "first record" for a later one through its end), the rest parsed.
            let todo: Vec<usize> = (0..n_ref).filter(|&t| !trav[t].is_empty()).collect();
            let next = std::sync::atomic::AtomicUsize::new(0);
            let parts: std::sync::Mutex<Vec<(usize, SiteSel)>> = std::sync::Mutex::new(Vec::new());
            let n_threads: usize = std::env::var("VARIANT_RS_THREADS").ok().and_then(|v| v.parse().ok()).unwrap_or_else(|| std::thread::available_parallelism().map_or(4, |n| n.get()));
            std::thread::scope(|scope| {
                for _ in 0..n_threads.clamp(1, todo.len().max(1)) {
                    scope.spawn(|| loop {
                        let k = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        let Some(&tid) = todo.get(k) else { break };
                        let (contig, ivs) = (header.ref_names[tid].as_str(), &trav[tid]);
                        let mut s = SiteSel::new(n_ref, min_af, max_af);
                        tl.for_each_line(contig, ivs, |line| {
                            let Some((c, start, end)) = peek_sites_line(line) else { return };
                            if c != contig {
                                return;
                            }
                            let j = ivs.partition_point(|iv| iv.1 < start);
                            if !(j < ivs.len() && ivs[j].0 <= start) {
                                s.note_end(tid, end);
                                return;
                            }
                            if let Some((_, r)) = parse_sites_line(line) {
                                s.feed(tid, c, &r, Some(ivs));
                            }
                        });
                        parts.lock().unwrap().push((tid, s));
                    });
                }
            });
            let mut parts = parts.into_inner().unwrap();
            parts.sort_by_key(|p| p.0);
            for (tid, mut s) in parts {
                sel.sites[tid] = std::mem::take(&mut s.sites[tid]);
                sel.saw_with |= s.saw_with;
                sel.saw_without |= s.saw_without;
                if sel.first_noaf.is_none() {
                    sel.first_noaf = s.first_noaf;
                }
            }
        } else {
            let tids: std::collections::HashMap<&str, usize> =
                header.ref_names.iter().enumerate().filter(|(t, _)| !trav[*t].is_empty()).map(|(t, n)| (n.as_str(), t)).collect();
            for_each_sites_rec(&vcf, |c, r| {
                if let Some(&t) = tids.get(c) {
                    sel.feed(t, c, &r, Some(&trav[t]));
                }
            });
        }
    }
    if let Some(at) = &sel.first_noaf {
        eprintln!("WARN Variant context at {at} lacks allele frequency (AF) field.");
    }
    let (saw_with, saw_without) = (sel.saw_with, sel.saw_without);
    let t_sel = t0.elapsed().as_secs_f64();
    let work: Vec<(usize, Vec<(i64, usize, usize, f64)>)> = sel.sites.into_iter().enumerate().filter(|(_, s)| !s.is_empty()).collect();
    // pileups per contig, in parallel (VARIANT_RS_THREADS); output in dictionary order
    let n_threads: usize = std::env::var("VARIANT_RS_THREADS").ok().and_then(|v| v.parse().ok()).unwrap_or_else(|| std::thread::available_parallelism().map_or(4, |n| n.get()));
    let results: std::sync::Mutex<Vec<Option<Vec<Site>>>> = std::sync::Mutex::new(vec![None; work.len()]);
    let next = std::sync::atomic::AtomicUsize::new(0);
    std::thread::scope(|scope| {
        for _ in 0..n_threads.clamp(1, work.len().max(1)) {
            scope.spawn(|| loop {
                let w = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                let Some((tid, sites)) = work.get(w) else { break };
                let tid = *tid;
                let contig = &header.ref_names[tid];
                let mut out_sites: Vec<Site> = Vec::new();
        let pos: Vec<i64> = sites.iter().map(|s| s.0).collect();
        let mut counts = vec![[0i32; 4]; sites.len()];
        let mut covered = vec![false; sites.len()];
        let query: Vec<(i64, i64)> = merge_intervals(pos.iter().map(|&p| (p, p)).collect());
        let mut stream = inp.stream(tid, &query, false, Arc::new(AtomicU64::new(0)));
        let mut buf = Vec::new();
        let mut rec = Rec::default();
        while stream.next_record(&mut buf) {
            let (rid, _, _) = bam::peek_core(&buf);
            if rid != tid as i32 {
                continue;
            }
            bam::parse_record(&buf, &mut rec);
            if !passes_gps_filters(&rec, &header.ref_lens) {
                continue;
            }
            // LocusIteratorByState leaves out bases inside the adaptor (ReadUtils.isBaseInsideAdaptor)
            let adaptor = adaptor_boundary(&rec).filter(|_| rec.tlen <= DEFAULT_ADAPTOR_SIZE);
            let rev = rec.flag & 0x10 != 0;
            let inside = |p: i64| adaptor.map_or(false, |b| if rev { p <= b } else { p >= b });
            let (mut rp, mut qp) = (rec.start, 0usize);
            for &(op, len) in &rec.cigar {
                let len64 = len as i64;
                match op {
                    OP_M | OP_EQ | OP_X | OP_D => {
                        let (a, b) = (rp, rp + len64 - 1);
                        let mut k = pos.partition_point(|&p| p < a);
                        while k < pos.len() && pos[k] <= b {
                            if inside(pos[k]) {
                                k += 1;
                                continue;
                            }
                            covered[k] = true;
                            if op != OP_D {
                                if let Some(bi) = base_index(rec.seq[qp + (pos[k] - rp) as usize]) {
                                    counts[k][bi] += 1;
                                }
                            }
                            k += 1;
                        }
                        rp += len64;
                        if op != OP_D {
                            qp += len as usize;
                        }
                    }
                    OP_N => rp += len64,
                    OP_I | OP_S => qp += len as usize,
                    _ => {}
                }
            }
        }
        for (k, s) in sites.iter().enumerate() {
            if !covered[k] {
                continue;
            }
            let c = counts[k];
            let total = c[0] + c[1] + c[2] + c[3];
            let altc = if s.2 < 4 { c[s.2] } else { panic!("bad alt base") };
            let refc = if s.1 < 4 { c[s.1] } else { panic!("bad ref base") };
            out_sites.push(Site { contig: contig.clone(), pos: s.0, refc, altc, other: total - altc - refc, total, af: s.3 });
        }
                    results.lock().unwrap()[w] = Some(out_sites);
            });
        }
    });
    let out_sites: Vec<Site> = results.into_inner().unwrap().into_iter().flatten().flatten().collect();
    if saw_without && !saw_with {
        panic!("No variants in population vcf had an allele frequency (AF) field.");
    }
    write_pileups(&output, &sample, &out_sites);
    eprintln!("getpileupsummaries: {} sites, {:.2}s (site selection {:.2}s)", out_sites.len(), t0.elapsed().as_secs_f64(), t_sel);
}

// ------------------------------------------------------------------------------------------------
// KernelSegmenter (ContaminationSegmenter's use: one window size, changepoints sorted by index)

const SEG_EPSILON: f64 = 1e-10;

fn seg_kernel(a: f64, b: f64) -> f64 {
    // a, b are alt fractions
    let m1 = jmin(a, 1.0 - a);
    let m2 = jmin(b, 1.0 - b);
    let d = m1 - m2;
    fexp(-(d * d) / (2.0 * 0.025))
}

struct Cost {
    d: f64,
    w: Vec<f64>,
    v: f64,
    c: f64,
}

fn segment_cost(start: usize, end: usize, z: &[Vec<f64>], diag: &[f64]) -> Cost {
    let n = z.len();
    let mut d = diag[start];
    let mut w = z[start].clone();
    let mut v = stream_sum(w.iter().map(|x| x * x));
    let indices: Vec<usize> = if start <= end { (start + 1..end + 1).collect() } else { (start + 1..n).chain(0..end + 1).collect() };
    for &t in &indices {
        d += diag[t];
        let mut zdotw = 0.0;
        for j in 0..w.len() {
            zdotw += z[t][j] * w[j];
            w[j] += z[t][j];
        }
        v += 2.0 * zdotw + diag[t];
    }
    let c = d - v / (indices.len() + 1) as f64;
    Cost { d, w, v, c }
}

fn window_costs(z: &[Vec<f64>], diag: &[f64], ws: usize) -> Vec<f64> {
    let n = z.len();
    let p = z[0].len();
    let mut center = 0usize;
    let mut start = (center + n + 1 - ws) % n; // (center - ws + 1 + N) % N
    let mut end = (center + ws) % n;
    let l = segment_cost(start, center, z, diag);
    let r = segment_cost(center + 1, end, z, diag);
    let t = segment_cost(start, end, z, diag);
    let (mut ld, mut lw, mut lv, mut lc) = (l.d, l.w, l.v, l.c);
    let (mut rd, mut rw, mut rv, mut rc) = (r.d, r.w, r.v, r.c);
    let (mut td, mut tw, mut tv, mut tc) = (t.d, t.w, t.v, t.c);
    let mut costs = vec![0.0f64; n];
    costs[center] = lc + rc - tc;
    let rec = 1.0 / ws as f64;
    while center < n {
        let cn = (center + 1) % n;
        let en = (end + 1) % n;
        let sub = |dd: &mut f64, ww: &mut Vec<f64>, vv: &mut f64, i: usize| {
            *dd -= diag[i];
            let mut zd = 0.0;
            for j in 0..p {
                zd += z[i][j] * ww[j];
                ww[j] -= z[i][j];
            }
            *vv += -2.0 * zd + diag[i];
        };
        let add = |dd: &mut f64, ww: &mut Vec<f64>, vv: &mut f64, i: usize| {
            *dd += diag[i];
            let mut zd = 0.0;
            for j in 0..p {
                zd += z[i][j] * ww[j];
                ww[j] += z[i][j];
            }
            *vv += 2.0 * zd + diag[i];
        };
        sub(&mut ld, &mut lw, &mut lv, start);
        add(&mut ld, &mut lw, &mut lv, cn);
        lc = ld - lv * rec;
        sub(&mut rd, &mut rw, &mut rv, cn);
        add(&mut rd, &mut rw, &mut rv, en);
        rc = rd - rv * rec;
        sub(&mut td, &mut tw, &mut tv, start);
        add(&mut td, &mut tw, &mut tv, en);
        tc = td - 0.5 * tv * rec;
        costs[cn] = lc + rc - tc;
        start = (start + 1) % n;
        end = en;
        center += 1;
    }
    costs
}

/// `List.indexOf(Collections.max(list))` with Double natural ordering.
fn index_of_max(v: &[f64]) -> usize {
    let mut best = 0;
    for i in 1..v.len() {
        if java_compare(v[i], v[best]) == Ordering::Greater {
            best = i;
        }
    }
    v.iter().position(|&x| x.to_bits() == v[best].to_bits() || (x.is_nan() && v[best].is_nan())).unwrap()
}

fn index_of_min(v: &[f64]) -> usize {
    let mut best = 0;
    for i in 1..v.len() {
        if java_compare(v[i], v[best]) == Ordering::Less {
            best = i;
        }
    }
    v.iter().position(|&x| x.to_bits() == v[best].to_bits() || (x.is_nan() && v[best].is_nan())).unwrap()
}

#[derive(Clone, Copy)]
struct Seg {
    start: usize,
    end: usize,
    cost: f64,
}

fn select_changepoints(cands: &[usize], max_cp: usize, lin: f64, loglin: f64, z: &[Vec<f64>], diag: &[f64]) -> Vec<usize> {
    let mut changepoints: Vec<usize> = Vec::new();
    let nd = z.len();
    let penalties: Vec<f64> =
        (0..=max_cp).map(|k| lin * k as f64 + loglin * k as f64 * jlog(nd as f64 / (k as f64 + SEG_EPSILON))).collect();
    let mut sorted: Vec<usize> = cands.to_vec();
    sorted.sort();
    sorted.dedup();
    let mut starts: Vec<usize> = sorted.iter().map(|&i| (i + 1).min(nd - 1)).collect();
    starts.insert(0, 0);
    let mut ends = sorted.clone();
    ends.push(nd - 1);
    let ns = starts.len();
    let mut segs: Vec<Seg> = (0..ns).map(|i| Seg { start: starts[i], end: ends[i], cost: segment_cost(starts[i], ends[i], z, diag).c }).collect();
    let mut totals: Vec<f64> = vec![stream_sum(segs.iter().map(|s| s.cost))];
    let mut pair: Vec<f64> = (0..ns.saturating_sub(1)).map(|i| segs[i].cost + segs[i + 1].cost).collect();
    let mut merged: Vec<f64> = (0..ns.saturating_sub(1)).map(|i| segment_cost(starts[i], ends[i + 1], z, diag).c).collect();
    let mut merging: Vec<f64> = (0..ns.saturating_sub(1)).map(|i| pair[i] - merged[i]).collect();
    for _ in 0..ns.saturating_sub(1) {
        let li = index_of_max(&merging);
        let new_cost = merged[li];
        let new_start = segs[li].start;
        let mergepoint = segs[li].end;
        let new_end = segs[li + 1].end;
        segs.remove(li);
        segs.remove(li);
        segs.insert(li, Seg { start: new_start, end: new_end, cost: new_cost });
        pair.remove(li);
        merged.remove(li);
        merging.remove(li);
        if li > 0 {
            pair[li - 1] = segs[li - 1].cost + segs[li].cost;
            merged[li - 1] = segment_cost(segs[li - 1].start, new_end, z, diag).c;
            merging[li - 1] = pair[li - 1] - merged[li - 1];
        }
        if li + 1 < segs.len() {
            pair[li] = segs[li].cost + segs[li + 1].cost;
            merged[li] = segment_cost(new_start, segs[li + 1].end, z, diag).c;
            merging[li] = pair[li] - merged[li];
        }
        totals.insert(0, stream_sum(segs.iter().map(|s| s.cost)));
        changepoints.insert(0, mergepoint);
    }
    let eff = max_cp.min(changepoints.len());
    let plus: Vec<f64> = (0..=eff).map(|i| totals[i] + penalties[i]).collect();
    let k = index_of_min(&plus);
    changepoints.truncate(k);
    changepoints
}

/// `KernelSegmenter.findChangepoints(10, SEGMENTATION_KERNEL, 100, [50], 1, 1, INDEX)` on alt fractions.
fn find_changepoints(data: &[f64]) -> Vec<usize> {
    const MAX_CP: usize = 10;
    const DIM: usize = 100;
    const WINDOW: usize = 50;
    if data.is_empty() {
        return Vec::new();
    }
    let n = data.len();
    let mut rng = hc_downsample::JavaRandom::new(1216);
    let ns = DIM.min(n);
    let sub: Vec<f64> = if ns == n { data.to_vec() } else { (0..ns).map(|_| data[rng.next_int(n as i32) as usize]).collect() };
    let mut k = vec![vec![0.0f64; ns]; ns];
    for i in 0..ns {
        for j in 0..i {
            let v = seg_kernel(sub[i], sub[j]);
            k[i][j] = v;
            k[j][i] = v;
        }
        k[i][i] = seg_kernel(sub[i], sub[i]);
    }
    let (sv, u) = svd(k);
    let inv: Vec<f64> = sv.iter().map(|&x| 1.0 / (x.sqrt() + SEG_EPSILON)).collect();
    let uu: Vec<Vec<f64>> = (0..ns).map(|i| (0..ns).map(|j| u[i][j] * inv[j]).collect()).collect();
    let red: Vec<Vec<f64>> = (0..n).map(|i| (0..ns).map(|j| seg_kernel(data[i], sub[j])).collect()).collect();
    // reducedKernelMatrix.multiply(subKernelUMatrix)
    let z: Vec<Vec<f64>> = red
        .iter()
        .map(|row| {
            (0..ns)
                .map(|c| {
                    let mut s = 0.0;
                    for i in 0..ns {
                        s += row[i] * uu[i][c];
                    }
                    s
                })
                .collect()
        })
        .collect();
    let diag: Vec<f64> = z
        .iter()
        .map(|row| {
            let mut s = 0.0;
            for &a in row {
                s += a * a;
            }
            let norm = s.sqrt();
            norm * norm
        })
        .collect();
    let mut cands: Vec<usize> = Vec::new();
    if 2 * WINDOW <= n {
        let costs = window_costs(&z, &diag, WINDOW);
        let mut minima = persistence_minima(&costs);
        if let Some(i) = minima.iter().position(|&x| x == 0) {
            minima.remove(i);
        }
        if let Some(i) = minima.iter().position(|&x| x == n - 1) {
            minima.remove(i);
        }
        cands.extend(minima.iter().take(MAX_CP));
    }
    let mut cps = select_changepoints(&cands, MAX_CP, 1.0, 1.0, &z, &diag);
    cps.sort();
    cps
}

/// Iteration order of a `java.util.HashMap<String, _>` filled by `computeIfAbsent` (as `Collectors.groupingBy`)
/// with the distinct `keys` in order: new nodes go to the head of their bin, the table doubles before an
/// insertion once `size > threshold`, and resizing splits bins keeping their order. Returns key indices.
fn java_hashmap_compute_if_absent_order(keys: &[&str]) -> Vec<usize> {
    let mut cap = 16usize;
    let mut table: Vec<Vec<(usize, u32)>> = vec![Vec::new(); cap];
    for (i, k) in keys.iter().enumerate() {
        if i > cap * 3 / 4 {
            let new_cap = cap * 2;
            let mut t2: Vec<Vec<(usize, u32)>> = vec![Vec::new(); new_cap];
            for bin in table {
                for (idx, h) in bin {
                    t2[(h as usize) & (new_cap - 1)].push((idx, h));
                }
            }
            table = t2;
            cap = new_cap;
        }
        let h = java_string_hash(k) as u32;
        let h = h ^ (h >> 16);
        table[(h as usize) & (cap - 1)].insert(0, (i, h));
    }
    table.into_iter().flatten().map(|(i, _)| i).collect()
}

/// `ContaminationSegmenter.findSegments`: segments (as site index lists into `sites`).
fn find_segments(sites: &[Site]) -> Vec<Vec<usize>> {
    // Collectors.groupingBy(contig) -> HashMap; values() in java.util.HashMap order
    let mut contigs: Vec<&str> = Vec::new();
    let mut members: Vec<Vec<usize>> = Vec::new();
    for (i, s) in sites.iter().enumerate() {
        match contigs.iter().position(|&c| c == s.contig) {
            Some(k) => members[k].push(i),
            None => {
                contigs.push(&s.contig);
                members.push(vec![i]);
            }
        }
    }
    let order = java_hashmap_compute_if_absent_order(&contigs);
    let mut out = Vec::new();
    for k in order {
        let cs = &members[k];
        let hets: Vec<usize> = cs.iter().copied().filter(|&i| (0.1..=0.9).contains(&sites[i].alt_fraction())).collect();
        if hets.is_empty() {
            continue;
        }
        let af: Vec<f64> = hets.iter().map(|&i| sites[i].alt_fraction()).collect();
        let mut cps: Vec<i64> = vec![-1];
        cps.extend(find_changepoints(&af).into_iter().map(|c| c as i64));
        cps.push(hets.len() as i64 - 1);
        for n in 0..cps.len() - 1 {
            let first = &sites[hets[(cps[n] + 1) as usize]];
            let last = &sites[hets[cps[n + 1] as usize]];
            let (c, a, b) = (&first.contig, first.pos, last.pos);
            let mut seg: Vec<usize> = (0..sites.len()).filter(|&i| sites[i].contig == *c && sites[i].pos >= a && sites[i].pos <= b).collect();
            seg.sort_by_key(|&i| sites[i].pos);
            out.push(seg);
        }
    }
    out
}

// ------------------------------------------------------------------------------------------------
// ContaminationModel

const INITIAL_MAF_THRESHOLD: f64 = 0.40;
const MAF_TO_SWITCH_TO_HOM_REF: f64 = 0.25;
const MAF_TO_SWITCH_TO_UNSCRUPULOUS_HOM_REF: f64 = 0.20;
const UNSCRUPULOUS_HOM_REF_ALLELE_FRACTION: f64 = 0.15;
const UNSCRUPULOUS_HOM_REF_PERCENTILE: f64 = 100.0 * (1.0 - 0.1);
const MINIMUM_UNSCRUPULOUS_HOM_REF_ALT_FRACTION_THRESHOLD: f64 = 0.1;
const MAF_STEP_SIZE: f64 = 0.04;
const PRECISION_FOR_STANDARD_ERROR: f64 = 0.000001;
const HOM_REF: usize = 0;
const HOM_ALT: usize = 3;
const NUM_ITERATIONS: usize = 3;
const MIN_FRACTION_OF_SITES_TO_USE: f64 = 0.25;
const MIN_RELATIVE_ERROR: f64 = 0.2;
const MIN_ABSOLUTE_ERROR: f64 = 0.001;
const CONTAMINATION_INITIAL_GUESSES: [f64; 4] = [0.02, 0.05, 0.1, 0.2];

pub struct Model<'a> {
    sites: &'a [Site],
    contamination: f64,
    error_rate: f64,
    mafs: Vec<f64>,
    segments: Vec<Vec<usize>>,
}

fn error_rate(sites: &[Site]) -> f64 {
    // mapToInt(...).sum() is an int sum
    let mut total: i32 = 0;
    let mut other: i32 = 0;
    for s in sites {
        total = total.wrapping_add(s.total);
        other = other.wrapping_add(s.other);
    }
    1.5 * (other as i64 as f64 / total as i64 as f64)
}

fn genotype_likelihoods(s: &Site, c: f64, err: f64, maf: f64) -> [f64; 4] {
    let f = s.af;
    let k = s.altc;
    let n = k + s.refc;
    let priors = [(1.0 - f) * (1.0 - f), f * (1.0 - f), f * (1.0 - f), f * f];
    let afs = [err / 3.0, maf, 1.0 - maf, 1.0 - err];
    let mut out = [0.0; 4];
    for g in 0..4 {
        out[g] = priors[g] * binomial_probability(n, k, (1.0 - c) * afs[g] + c * f);
    }
    out
}

fn sum4(v: &[f64; 4]) -> f64 {
    let mut s = 0.0;
    for &x in v {
        s += x;
    }
    s
}

fn segment_log_likelihood(sites: &[Site], seg: &[usize], c: f64, err: f64, maf: f64) -> f64 {
    stream_sum(seg.iter().map(|&i| jlog(sum4(&genotype_likelihoods(&sites[i], c, err, maf)))))
}

fn model_log_likelihood(sites: &[Site], segs: &[&Vec<usize>], c: f64, err: f64, mafs: &[f64]) -> f64 {
    let mut s = 0.0;
    for n in 0..segs.len() {
        s += segment_log_likelihood(sites, segs[n], c, err, mafs[n]);
    }
    s
}

fn calc_maf(sites: &[Site], seg: &[usize], c: f64, err: f64) -> f64 {
    brent_max_pv(|maf| segment_log_likelihood(sites, seg, c, err, maf), 0.1, 0.5, 0.4, 0.01, 0.01, 20).0
}

fn calc_contamination(sites: &[Site], err: f64, segs: &[&Vec<usize>], mafs: &[f64]) -> f64 {
    let optima: Vec<(f64, f64)> = CONTAMINATION_INITIAL_GUESSES
        .iter()
        .map(|&g| brent_max_pv(|c| model_log_likelihood(sites, segs, c, err, mafs), 0.0, 0.5, g, 1.0e-4, 1.0e-4, 30))
        .collect();
    // Collections.max(comparingDouble(value)): the first maximal element
    let mut best = optima[0];
    for &o in &optima[1..] {
        if java_compare(o.1, best.1) == Ordering::Greater {
            best = o;
        }
    }
    best.0
}

fn non_loh<'b>(segments: &'b [Vec<usize>], mafs: &[f64]) -> (Vec<&'b Vec<usize>>, Vec<f64>) {
    let num_sites: i32 = segments.iter().map(|s| s.len() as i32).sum();
    let mut min_maf = INITIAL_MAF_THRESHOLD;
    while min_maf > 0.0 {
        let idx: Vec<usize> = (0..segments.len()).filter(|&s| mafs[s] > min_maf).collect();
        let n: i32 = idx.iter().map(|&s| segments[s].len() as i32).sum();
        if n as f64 / num_sites as f64 > MIN_FRACTION_OF_SITES_TO_USE {
            return (idx.iter().map(|&s| &segments[s]).collect(), idx.iter().map(|&s| mafs[s]).collect());
        }
        min_maf -= MAF_STEP_SIZE;
    }
    (segments.iter().collect(), mafs.to_vec())
}

#[derive(Clone, Copy, PartialEq)]
enum Strategy {
    HomAlt,
    HomRef,
    UnscrupulousHomRef,
}

/// `MathUtils.binarySearchFindZero`
fn binary_search_find_zero(f: impl Fn(f64) -> f64, lower: f64, upper: f64, precision: f64) -> Option<f64> {
    let signum = |x: f64| if x < 0.0 { -1.0 } else if x > 0.0 { 1.0 } else { x };
    let (mut bottom, mut top) = (lower, upper);
    while top - bottom > precision {
        let mid = (bottom + top) / 2.0;
        let bv = f(bottom);
        let tv = f(top);
        let mv = f(mid);
        if signum(bv) == signum(tv) {
            return None;
        }
        if signum(bv) == signum(mv) {
            bottom = mid;
        } else {
            top = mid;
        }
    }
    Some((bottom + top) / 2.0)
}

/// `Math.round(double)` for the values seen here (NaN -> 0).
fn java_round(x: f64) -> i64 {
    if x.is_nan() {
        0
    } else {
        (x + 0.5).floor() as i64
    }
}

impl<'a> Model<'a> {
    pub fn new(sites: &'a [Site]) -> Model<'a> {
        let err = error_rate(sites);
        let segments = find_segments(sites);
        let mut mafs = vec![0.5; segments.len()];
        let mut c = 0.0;
        for _ in 0..NUM_ITERATIONS {
            for s in 0..segments.len() {
                mafs[s] = calc_maf(sites, &segments[s], c, err);
            }
            let (segs, m) = non_loh(&segments, &mafs);
            c = calc_contamination(sites, err, &segs, &m);
        }
        Model { sites, contamination: c, error_rate: err, mafs, segments }
    }

    fn probability(&self, s: &Site, maf: f64, genotype: usize) -> f64 {
        let l = genotype_likelihoods(s, self.contamination, self.error_rate, maf);
        l[genotype] / sum4(&l)
    }

    fn get_type(&self, genotype: usize, min_maf: f64) -> Vec<&'a Site> {
        let mut out = Vec::new();
        for n in 0..self.segments.len() {
            if self.mafs[n] > min_maf {
                for &i in &self.segments[n] {
                    let s = &self.sites[i];
                    if self.probability(s, self.mafs[n], genotype) > 0.5 {
                        out.push(s);
                    }
                }
            }
        }
        out
    }

    pub fn segmentation_records(&self) -> Vec<(String, i64, i64, f64)> {
        self.segments
            .iter()
            .zip(&self.mafs)
            .map(|(seg, &m)| {
                let (a, b) = (&self.sites[seg[0]], &self.sites[seg[seg.len() - 1]]);
                (a.contig.clone(), a.pos, b.pos, m)
            })
            .collect()
    }

    fn calculate(&self, strategy: Strategy, tumor: &[Site], min_maf: f64) -> (f64, f64) {
        let use_hom_alt = strategy == Strategy::HomAlt;
        let genotyping: Vec<&Site> = match strategy {
            Strategy::HomAlt => self.get_type(HOM_ALT, min_maf),
            Strategy::HomRef => self.get_type(HOM_REF, min_maf),
            Strategy::UnscrupulousHomRef => {
                let cands: Vec<&Site> = tumor.iter().filter(|s| s.alt_fraction() < UNSCRUPULOUS_HOM_REF_ALLELE_FRACTION).collect();
                let afs: Vec<f64> = cands.iter().map(|s| s.alt_fraction()).collect();
                let thr = jmax(MINIMUM_UNSCRUPULOUS_HOM_REF_ALT_FRACTION_THRESHOLD, percentile(&afs, UNSCRUPULOUS_HOM_REF_PERCENTILE));
                cands.into_iter().filter(|s| s.alt_fraction() <= thr).collect()
            }
        };
        let keys: std::collections::HashSet<(&str, i64)> = genotyping.iter().map(|s| (s.contig.as_str(), s.pos)).collect();
        let homs: Vec<&Site> = tumor.iter().filter(|s| keys.contains(&(s.contig.as_str(), s.pos))).collect();
        let tumor_err = error_rate(tumor);
        let opp_count = |s: &Site| if use_hom_alt { s.refc } else { s.altc };
        let opp_freq = |s: &Site| if use_hom_alt { s.ref_frequency() } else { s.af };
        let total_depth: i64 = homs.iter().map(|s| s.total as i64).sum();
        let opp_depth: i64 = homs.iter().map(|s| opp_count(s) as i64).sum();
        let error_depth = java_round(total_depth as f64 * tumor_err / 3.0);
        let contam_opp_depth = (opp_depth - error_depth).max(0);
        let weighted = stream_sum(homs.iter().map(|s| s.total as f64 * opp_freq(s)));
        let estimate = contam_opp_depth as f64 / weighted;
        let coeff1 = stream_sum(homs.iter().map(|s| opp_freq(s) * s.total as f64));
        let coeff2 = stream_sum(homs.iter().map(|s| {
            let t = s.total as f64;
            opp_freq(s) * (1.0 - opp_freq(s)) * (t * t)
        }));
        let empty = homs.is_empty();
        let error_func = |c: f64| if empty { 1.0 } else { (coeff1 * c * (1.0 - c) + coeff2 * c * c).sqrt() / weighted };
        let upper = binary_search_find_zero(|c| c - error_func(c) - estimate, estimate, 1.0, PRECISION_FOR_STANDARD_ERROR);
        let lower = binary_search_find_zero(|c| c + error_func(c) - estimate, 0.0, estimate, PRECISION_FOR_STANDARD_ERROR);
        let std_err = match (upper, lower) {
            (Some(u), Some(l)) => jmax(u - estimate, estimate - l),
            _ => error_func(estimate),
        };
        (jmin(estimate, 1.0), std_err)
    }

    pub fn contamination_from_homs(&self, tumor: &[Site]) -> (f64, f64) {
        let mut min_maf = INITIAL_MAF_THRESHOLD;
        while min_maf >= 0.0 {
            let strategy = if min_maf > MAF_TO_SWITCH_TO_HOM_REF {
                Strategy::HomAlt
            } else if min_maf > MAF_TO_SWITCH_TO_UNSCRUPULOUS_HOM_REF {
                Strategy::HomRef
            } else {
                Strategy::UnscrupulousHomRef
            };
            let r = self.calculate(strategy, tumor, min_maf);
            if !r.0.is_nan() && r.1 < r.0 * MIN_RELATIVE_ERROR + MIN_ABSOLUTE_ERROR {
                return r;
            }
            min_maf -= MAF_STEP_SIZE;
        }
        let r = self.calculate(Strategy::UnscrupulousHomRef, tumor, 0.0);
        if r.0.is_nan() { (0.0, 1.0) } else { r }
    }
}

// ------------------------------------------------------------------------------------------------
// CalculateContamination

const CC_USAGE: &str = "usage: variant-rs calculatecontamination -I tumor.pileups.table [-matched normal.pileups.table] -O contamination.table \
[--tumor-segmentation segments.table] [--low-coverage-ratio-threshold 0.5] [--high-coverage-ratio-threshold 3.0]";

const MIN_COVERAGE: i32 = 10;

fn filter_sites_by_coverage(all: Vec<Site>, low_ratio: f64, high_ratio: f64) -> Vec<Site> {
    let covered: Vec<Site> = all.into_iter().filter(|s| s.total > MIN_COVERAGE).collect();
    let cov: Vec<f64> = covered.iter().map(|s| s.total as f64).collect();
    let median = percentile(&cov, 50.0);
    let mean = mean(&cov);
    let low = median * low_ratio;
    let high = mean * high_ratio;
    covered.into_iter().filter(|s| (s.total as f64) > low && (s.total as f64) < high).collect()
}

pub fn run_cc(argv: &[String]) {
    let t0 = std::time::Instant::now();
    let (mut input, mut matched, mut output, mut segs_out) = (None::<String>, None::<String>, None::<String>, None::<String>);
    let (mut low, mut high) = (0.5f64, 3.0f64);
    let mut i = 0;
    while i < argv.len() {
        let v = || argv.get(i + 1).cloned().unwrap_or_else(|| panic!("{CC_USAGE}"));
        match argv[i].as_str() {
            "-I" | "--input" => input = Some(v()),
            "-matched" | "--matched-normal" => matched = Some(v()),
            "-O" | "--output" => output = Some(v()),
            "-segments" | "--tumor-segmentation" => segs_out = Some(v()),
            "--low-coverage-ratio-threshold" => low = v().parse().unwrap(),
            "--high-coverage-ratio-threshold" => high = v().parse().unwrap(),
            a => panic!("unknown argument {a}\n{CC_USAGE}"),
        }
        i += 2;
    }
    let (input, output) = (input.expect(CC_USAGE), output.expect(CC_USAGE));
    let (sample, all) = read_pileups(&input);
    let sites = filter_sites_by_coverage(all, low, high);
    let genotyping_sites = matched.as_ref().map(|m| filter_sites_by_coverage(read_pileups(m).1, low, high));
    let genotyping_model = Model::new(genotyping_sites.as_deref().unwrap_or(&sites));
    if let Some(path) = &segs_out {
        let tumor_model;
        let m = if genotyping_sites.is_none() {
            &genotyping_model
        } else {
            tumor_model = Model::new(&sites);
            &tumor_model
        };
        let mut w = std::io::BufWriter::new(std::fs::File::create(path).unwrap_or_else(|e| panic!("{path}: {e}")));
        writeln!(w, "#<METADATA>SAMPLE={sample}").unwrap();
        writeln!(w, "contig\tstart\tend\tminor_allele_fraction").unwrap();
        for (c, s, e, maf) in m.segmentation_records() {
            writeln!(w, "{c}\t{s}\t{e}\t{}", java_double_to_string(maf)).unwrap();
        }
    }
    let (c, e) = genotyping_model.contamination_from_homs(&sites);
    let mut w = std::io::BufWriter::new(std::fs::File::create(&output).unwrap_or_else(|e| panic!("{output}: {e}")));
    writeln!(w, "sample\tcontamination\terror").unwrap();
    writeln!(w, "{sample}\t{}\t{}", java_double_to_string(c), java_double_to_string(e)).unwrap();
    eprintln!("calculatecontamination: {sample} contamination {c} error {e}, {:.2}s", t0.elapsed().as_secs_f64());
}

