//! ApplyBQSR (`variant-rs applybqsr`): GATK 4.7 `ApplyBQSR` / `BQSRReadTransformer` (no dynamic quantization,
//! `--preserve-qscores-less-than 6`, OQ tags not written), with `--use-original-qualities`,
//! `--static-quantized-quals` / `--round-down-quantized` and `--create-output-bam-md5`. Records are copied
//! byte for byte except the quality string (and BI/BD tags, which GATK clears), as htsjdk does.

use super::covariates::{self, ReadGroups, Scratch, CONTEXT_DIM, CYCLE_DIM, PRESERVE_QSCORES_LESS_THAN, QUAL_DIM};
use super::datum::{Datum, MULTIPLIER};
use crate::bamwrite;
use crate::input::Input;
use crate::reference::Fasta;
use hc_preproc::Read;
use rayon::prelude::*;
use std::collections::{BTreeSet, HashMap};
use std::sync::atomic::AtomicU64;
use std::sync::Arc;
use std::time::Instant;

/// The recalibration model of a report, reduced to what `BQSRReadTransformer` needs: the
/// posterior quality per (read group, reported quality) and the context/cycle deltas.
pub struct Model {
    pub rg_names: Vec<String>, // sorted (RecalibrationReport keys read groups in TreeSet order)
    n_rg: usize,
    post: Vec<f64>,       // [rg][q]
    /// [rg]: the read group's estimate, the posterior of a quality with no cell in the tables (94-127)
    emp_rg: Vec<f64>,
    d_context: Vec<f64>,  // [rg][q][key], NaN = no datum
    d_cycle: Vec<f64>,    // [rg][q][key]
    quantized: [u8; QUAL_DIM],
    /// `--static-quantized-quals`: applied after `quantized` (BQSRReadTransformer quantizes twice)
    static_map: Option<Vec<u8>>,
}

/// `BQSRReadTransformer.constructStaticQuantizedMapping` (QualityUtils.MAX_QUAL = 254 entries, MIN_USABLE_Q_SCORE
/// = 6; `qualToProb` of an int is the double overload, 1 - 10^(q/-10)).
pub fn static_quantized_mapping(quals: &[i32], round_down: bool) -> Vec<u8> {
    const MAX_QUAL: i32 = 254;
    const MIN_USABLE: i32 = 6;
    let mut q: Vec<i32> = quals.to_vec();
    if q.is_empty() {
        return (0..MAX_QUAL).map(|i| i as u8).collect();
    }
    let prob = |x: i32| 1.0 - crate::jmath::pow(10.0, x as f64 / -10.0);
    let mut mapping = vec![0u8; MAX_QUAL as usize];
    q.sort();
    for i in 0..MIN_USABLE {
        mapping[i as usize] = i as u8;
    }
    if q.len() == 1 {
        for m in mapping.iter_mut().skip(MIN_USABLE as usize) {
            *m = q[0] as u8;
        }
        return mapping;
    }
    let mut prev = MIN_USABLE;
    let mut prev_p = prob(prev);
    for &next in &q {
        let next_p = prob(next);
        for i in prev..next {
            mapping[i as usize] = if round_down {
                prev as u8
            } else {
                let ip = prob(i);
                if ip - prev_p > next_p - ip { next as u8 } else { prev as u8 }
            };
        }
        prev = next;
        prev_p = next_p;
    }
    for i in prev..MAX_QUAL {
        mapping[i as usize] = prev as u8;
    }
    mapping
}

struct ParsedTable {
    cols: Vec<String>,
    rows: Vec<Vec<String>>,
}

fn parse_report(text: &str) -> HashMap<String, ParsedTable> {
    let mut tables = HashMap::new();
    let mut lines = text.lines().peekable();
    while let Some(l) = lines.next() {
        if !l.starts_with("#:GATKTable:") {
            continue;
        }
        let name_line = lines.next().expect("table name line");
        let name = name_line.trim_start_matches("#:GATKTable:").split(':').next().unwrap().to_string();
        let cols: Vec<String> = lines.next().expect("column line").split_whitespace().map(String::from).collect();
        let mut rows = Vec::new();
        while let Some(&r) = lines.peek() {
            if r.trim().is_empty() {
                break;
            }
            lines.next();
            let f: Vec<String> = r.split_whitespace().map(String::from).collect();
            assert_eq!(f.len(), cols.len(), "table {name}: malformed row {r:?}");
            rows.push(f);
        }
        tables.insert(name, ParsedTable { cols, rows });
    }
    tables
}

impl ParsedTable {
    fn col(&self, name: &str) -> usize {
        self.cols.iter().position(|c| c == name).unwrap_or_else(|| panic!("column {name} missing"))
    }
}

fn datum(obs: &str, err: &str, reported: f64) -> Datum {
    // RecalibrationReport.getRecalDatum: RecalDatum(nObs, nErrors, 1) then setReportedQuality
    Datum { obs: obs.parse().expect("Observations"), mism: err.parse::<f64>().expect("Errors") * MULTIPLIER, reported }
}

impl Model {
    pub fn load(path: &str) -> Model {
        let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
        let t = parse_report(&text);
        let get = |n: &str| t.get(n).unwrap_or_else(|| panic!("{path}: no {n} table"));
        let (t0, t1, t2, tq) = (get("RecalTable0"), get("RecalTable1"), get("RecalTable2"), get("Quantized"));
        let rg_set: BTreeSet<String> = t0.rows.iter().map(|r| r[t0.col("ReadGroup")].clone()).collect();
        let rg_names: Vec<String> = rg_set.into_iter().collect();
        let rg_key: HashMap<&str, usize> = rg_names.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
        let n_rg = rg_names.len();
        // quantization: noQuantization() (default --quantize-quals 0) keeps 0..92, 93 from the table
        let mut quantized = [0u8; QUAL_DIM];
        let (cq, cv) = (tq.col("QualityScore"), tq.col("QuantizedScore"));
        for r in &tq.rows {
            quantized[r[cq].parse::<usize>().unwrap()] = r[cv].parse().unwrap();
        }
        for (i, q) in quantized.iter_mut().enumerate().take(QUAL_DIM - 1) {
            *q = i as u8;
        }
        let only_m = |r: &Vec<String>, c: usize| r[c] == "M";
        // read group datums
        let mut rg_d: Vec<Option<Datum>> = vec![None; n_rg];
        let (c_rg, c_ev, c_eq, c_obs, c_err) = (t0.col("ReadGroup"), t0.col("EventType"), t0.col("EstimatedQReported"), t0.col("Observations"), t0.col("Errors"));
        for r in t0.rows.iter().filter(|r| only_m(r, c_ev)) {
            rg_d[rg_key[r[c_rg].as_str()]] = Some(datum(&r[c_obs], &r[c_err], r[c_eq].parse().unwrap()));
        }
        // empirical quality of each read group, with its reported quality as the prior
        let emp_rg: Vec<Option<f64>> = rg_d.iter().map(|d| d.map(|d| d.empirical_quality(d.reported) as f64)).collect();
        let mut post = vec![f64::NAN; n_rg * QUAL_DIM];
        for rg in 0..n_rg {
            if let Some(e) = emp_rg[rg] {
                for q in 0..QUAL_DIM {
                    post[rg * QUAL_DIM + q] = e;
                }
            }
        }
        let (c_rg, c_q, c_ev, c_obs, c_err) = (t1.col("ReadGroup"), t1.col("QualityScore"), t1.col("EventType"), t1.col("Observations"), t1.col("Errors"));
        for r in t1.rows.iter().filter(|r| only_m(r, c_ev)) {
            let rg = rg_key[r[c_rg].as_str()];
            let q: usize = r[c_q].parse().unwrap();
            let d = datum(&r[c_obs], &r[c_err], q as f64);
            if let Some(e) = emp_rg[rg] {
                post[rg * QUAL_DIM + q] = d.empirical_quality(e) as f64;
            }
        }
        let (c_rg, c_q, c_v, c_n, c_ev, c_obs, c_err) =
            (t2.col("ReadGroup"), t2.col("QualityScore"), t2.col("CovariateValue"), t2.col("CovariateName"), t2.col("EventType"), t2.col("Observations"), t2.col("Errors"));
        let deltas: Vec<(bool, usize, f64)> = t2
            .rows
            .par_iter()
            .filter(|r| only_m(r, c_ev))
            .map(|r| {
                let rg = rg_key[r[c_rg].as_str()];
                let q: usize = r[c_q].parse().unwrap();
                let d = datum(&r[c_obs], &r[c_err], q as f64);
                let p = post[rg * QUAL_DIM + q];
                let delta = d.empirical_quality(p) as f64 - p;
                let is_ctx = r[c_n] == "Context";
                let key = if is_ctx { context_key(&r[c_v]) } else { cycle_key(r[c_v].parse().unwrap()) };
                let dim = if is_ctx { CONTEXT_DIM } else { CYCLE_DIM };
                (is_ctx, (rg * QUAL_DIM + q) * dim + key, delta)
            })
            .collect();
        let mut d_context = vec![f64::NAN; n_rg * QUAL_DIM * CONTEXT_DIM];
        let mut d_cycle = vec![f64::NAN; n_rg * QUAL_DIM * CYCLE_DIM];
        for (is_ctx, i, d) in deltas {
            if is_ctx {
                d_context[i] = d;
            } else {
                d_cycle[i] = d;
            }
        }
        let emp_rg = emp_rg.iter().map(|e| e.unwrap_or(f64::NAN)).collect();
        Model { rg_names, n_rg, post, emp_rg, d_context, d_cycle, quantized, static_map: None }
    }

    /// `BQSRReadTransformer.apply` on the quality string of one read.
    fn recalibrate(&self, rg: usize, bases: &[u8], quals: &mut [u8], flag: u16, s: &mut Scratch) {
        // a byte >= 128 is a negative quality to Java: GATK sets the read's qualities and fails
        if let Some(&b) = quals.iter().find(|&&q| q >= 128) {
            panic!("Base quality score {} is invalid", b as i8);
        }
        covariates::compute(bases, quals, flag, s);
        for i in 0..quals.len() {
            let q = quals[i];
            if q < PRESERVE_QSCORES_LESS_THAN {
                continue;
            }
            if q as usize >= QUAL_DIM {
                // 94-127: the tables have no cell (their lookups return null): the read group's estimate alone
                let raw = self.emp_rg[rg];
                let r = if raw > 0.0 { (raw + 0.5) as i32 } else { (raw - 0.5) as i32 };
                let qq = self.quantized[r.clamp(1, 93) as usize];
                quals[i] = match &self.static_map {
                    Some(m) => m[qq as usize],
                    None => qq,
                };
                continue;
            }
            let base = rg * QUAL_DIM + q as usize;
            let post = self.post[base];
            let mut delta = 0.0;
            let ctx = s.context[i];
            if ctx >= 0 {
                let d = self.d_context[base * CONTEXT_DIM + ctx as usize];
                if !d.is_nan() {
                    delta += d;
                }
            }
            let d = self.d_cycle[base * CYCLE_DIM + s.cycle[i] as usize];
            if !d.is_nan() {
                delta += d;
            }
            let raw = post + delta;
            let r = if raw > 0.0 { (raw + 0.5) as i32 } else { (raw - 0.5) as i32 }; // MathUtils.fastRound
            let qq = self.quantized[r.clamp(1, 93) as usize];
            quals[i] = match &self.static_map {
                Some(m) => m[qq as usize],
                None => qq,
            };
        }
    }
}

fn context_key(s: &str) -> usize {
    let mut key = s.len();
    let mut off = 4;
    for b in s.bytes() {
        key |= (covariates::base_index(b) as usize) << off;
        off += 2;
    }
    key
}

fn cycle_key(c: i32) -> usize {
    ((c.unsigned_abs() as usize) << 1) + (c < 0) as usize
}

struct Args {
    bam: String,
    fasta: String,
    table: String,
    out: String,
    intervals: Vec<String>,
    threads: usize,
    level: i32,
    chunk: i64,
    use_oq: bool,
    static_quals: Vec<i32>,
    round_down: bool,
    md5: bool,
}

const USAGE: &str = "usage: variant-rs applybqsr -I in.bam -R ref.fa --bqsr-recal-file recal.table -O out.bam [-L interval ...] \
[--use-original-qualities] [--static-quantized-quals Q ...] [--round-down-quantized] [--create-output-bam-md5] \
[--threads N (16)] [--compression-level N (5)]";

fn parse_args(a: &[String]) -> Args {
    let mut r = Args { bam: String::new(), fasta: String::new(), table: String::new(), out: String::new(), intervals: vec![], threads: 16, level: 5, chunk: 500_000, use_oq: false, static_quals: vec![], round_down: false, md5: false };
    let mut i = 0;
    while i < a.len() {
        let v = || a.get(i + 1).cloned().unwrap_or_else(|| panic!("{} needs a value\n{USAGE}", a[i]));
        let flag = match a[i].as_str() {
            "--use-original-qualities" | "-OQ" => Some(&mut r.use_oq),
            "--round-down-quantized" => Some(&mut r.round_down),
            "--create-output-bam-md5" => Some(&mut r.md5),
            _ => None,
        };
        if let Some(f) = flag {
            *f = true;
            i += 1;
            continue;
        }
        match a[i].as_str() {
            "-I" | "--input" => r.bam = v(),
            "-R" | "--reference" => r.fasta = v(),
            "--static-quantized-quals" => r.static_quals.push(v().parse().expect("--static-quantized-quals")),
            "--bqsr-recal-file" | "-bqsr" => r.table = v(),
            "-O" | "--output" => r.out = v(),
            "-L" | "--intervals" => r.intervals.push(v()),
            "--threads" => r.threads = v().parse().expect("--threads"),
            "--compression-level" => r.level = v().parse().expect("--compression-level"),
            "--chunk" => r.chunk = v().parse().expect("--chunk"),
            x => panic!("unknown argument {x}\n{USAGE}"),
        }
        i += 2;
    }
    assert!(!r.bam.is_empty() && !r.fasta.is_empty() && !r.out.is_empty() && !r.table.is_empty(), "{USAGE}");
    r
}

/// GATK's default `WellformedReadFilter` (ApplyBQSR drops reads failing it).
fn wellformed(r: &Read, has_rg: bool, qual_missing: bool, ref_lens: &[i64]) -> bool {
    let unmapped = r.flag & 0x4 != 0 || r.ref_id < 0 || r.start < 1;
    if !unmapped {
        if r.ref_id as usize >= ref_lens.len() || r.start > ref_lens[r.ref_id as usize] {
            return false;
        }
        let cl: usize = r.cigar.iter().filter(|e| matches!(e.0, b'M' | b'I' | b'S' | b'=' | b'X')).map(|e| e.1 as usize).sum();
        if cl != r.seq.len() || r.end - r.start + 1 < 0 {
            return false;
        }
    }
    has_rg && !qual_missing && r.qual.len() == r.seq.len() && !r.seq.is_empty() && !r.cigar.iter().any(|e| e.0 == b'N')
}

pub fn run(argv: &[String]) {
    let args = parse_args(argv);
    rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();
    let t0 = Instant::now();
    let mut model = Model::load(&args.table);
    if !args.static_quals.is_empty() {
        model.static_map = Some(static_quantized_mapping(&args.static_quals, args.round_down));
    }
    let (input, header) = Input::open(&args.bam, &args.fasta);
    let rgs = ReadGroups::from_header(&header.text);
    // read group ID -> model key (via the PU/ID identifier)
    let mut rg_map: HashMap<Vec<u8>, usize> = HashMap::new();
    for (id, &k) in &rgs.key_by_id {
        let name = &rgs.names[k];
        if let Some(mk) = model.rg_names.iter().position(|n| n == name) {
            rg_map.insert(id.clone(), mk);
        }
    }
    let _ = Fasta::open(&args.fasta); // reference must exist (CRAM decoding)
    eprintln!("model: {} read group(s) loaded in {:.1}s", model.n_rg, t0.elapsed().as_secs_f64());

    let traversal: Vec<(usize, i64, i64)> = if args.intervals.is_empty() {
        header.ref_lens.iter().enumerate().map(|(t, &l)| (t, 1, l)).collect()
    } else {
        let dict = crate::intervals::Dict::from_reference(&args.fasta);
        let mut v = Vec::new();
        for l in &args.intervals {
            v.extend(crate::intervals::parse_l_arg(l, &dict).unwrap_or_else(|e| panic!("-L {l}: {e}")));
        }
        let mut t: Vec<(usize, i64, i64)> = crate::intervals::sort_merge(v)
            .into_iter()
            .map(|iv| {
                let name = &dict.names[iv.tid];
                (header.ref_names.iter().position(|n| n == name).unwrap_or_else(|| panic!("{name} not in the BAM header")), iv.start, iv.end)
            })
            .collect();
        t.sort_by_key(|x| (x.0, x.1));
        t
    };
    let chunks = bamwrite::plan_chunks(&traversal, header.ref_lens.len(), args.chunk, args.intervals.is_empty());
    let tail = bamwrite::tail_range(&input, &args.bam);
    let cmdline = std::env::args().collect::<Vec<_>>().join(" ");
    let pg = format!("@PG\tID:variant-rs ApplyBQSR\tVN:{}\tCL:{}\tPN:variant-rs ApplyBQSR", env!("CARGO_PKG_VERSION"), cmdline);
    let text = bamwrite::header_with_pg(&header.text, &[pg]);
    let mut bam = bamwrite::BamFile::create(&args.out, &text, &header.ref_names, &header.ref_lens, args.level, args.md5);
    let big = bam.too_large.clone();

    let io_ns = Arc::new(AtomicU64::new(0));
    let n_dropped = AtomicU64::new(0);
    let process = |ch: &bamwrite::Chunk| -> bamwrite::ChunkOut {
        let mut w = bamwrite::RecWriter::new(args.level, &big);
        let Some(mut stream) = bamwrite::open_chunk(&input, ch, tail, &io_ns) else { return w.finish() };
        let mut rec = Vec::new();
        let mut r = super::new_read();
        let mut rg = Vec::new();
        let mut oq = Vec::new();
        let mut s = Scratch::default();
        let mut ts = bamwrite::TagScratch::default();
        while stream.next_record(&mut rec) {
            if !bamwrite::in_chunk(&rec, ch) {
                continue;
            }
            let (has_rg, qual_missing) = super::decode(&rec, &mut r, &mut rg, &mut oq);
            if !wellformed(&r, has_rg, qual_missing, &header.ref_lens) {
                n_dropped.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                continue;
            }
            let rgk = *rg_map.get(&rg[..]).unwrap_or_else(|| panic!("Read group {} not found in the recalibration table", String::from_utf8_lossy(&rg)));
            let l_read_name = rec[8] as usize;
            let n_cigar = u16::from_le_bytes([rec[12], rec[13]]) as usize;
            let l_seq = r.seq.len();
            let qo = 32 + l_read_name + 4 * n_cigar + (l_seq + 1) / 2;
            if args.use_oq && !oq.is_empty() {
                // ReadUtils.resetOriginalBaseQualities before recalibrating (the OQ tag stays)
                let (head, quals) = rec.split_at_mut(qo);
                super::reset_original_qualities(&mut quals[..l_seq], &oq, head);
            }
            model.recalibrate(rgk, &r.seq, &mut rec[qo..qo + l_seq], r.flag, &mut s);
            bamwrite::normalize_tags(&mut rec, qo + l_seq, &[*b"BI", *b"BD"], None, &mut ts);
            w.push(&rec);
        }
        w.finish()
    };
    let t1 = Instant::now();
    let window = 4 * args.threads.max(1);
    for group in chunks.chunks(window) {
        let outs: Vec<bamwrite::ChunkOut> = group.par_iter().map(process).collect();
        for (ch, o) in group.iter().zip(outs.iter()) {
            bam.append(ch.tid, o);
        }
    }
    let n_records = bam.n_records;
    bam.finish(&args.out);
    eprintln!(
        "records written: {n_records} (dropped by WellformedReadFilter: {}) in {:.1}s; total {:.1}s",
        n_dropped.into_inner(),
        t1.elapsed().as_secs_f64(),
        t0.elapsed().as_secs_f64()
    );
}

/// ApplyBQSR fused into another tool's read stream (`variant-rs --bqsr-recal-file`): recalibrates the
/// quality bytes of BAM records in place, as `BQSRReadTransformer` would before the tool sees them.
pub struct Recalibrator {
    model: Model,
    rg_map: HashMap<Vec<u8>, usize>,
}

#[derive(Default)]
pub struct RecScratch {
    read: Option<Read>,
    rg: Vec<u8>,
    oq: Vec<u8>,
    cov: Scratch,
}

impl Recalibrator {
    pub fn new(table: &str, header_text: &str) -> Recalibrator {
        let model = Model::load(table);
        let rgs = ReadGroups::from_header(header_text);
        let mut rg_map = HashMap::new();
        for (id, &k) in &rgs.key_by_id {
            if let Some(mk) = model.rg_names.iter().position(|n| n == &rgs.names[k]) {
                rg_map.insert(id.clone(), mk);
            }
        }
        Recalibrator { model, rg_map }
    }

    /// `rec` = record bytes after block_size. Reads without a read group or qualities are left alone
    /// (GATK's WellformedReadFilter drops them before and after ApplyBQSR).
    pub fn recalibrate_record(&self, rec: &mut [u8], s: &mut RecScratch) {
        let r = s.read.get_or_insert_with(super::new_read);
        let (has_rg, qual_missing) = super::decode(rec, r, &mut s.rg, &mut s.oq);
        if !has_rg || qual_missing || r.seq.is_empty() || r.qual.len() != r.seq.len() {
            return;
        }
        let rgk = *self.rg_map.get(&s.rg[..]).unwrap_or_else(|| panic!("Read group {} not found in the recalibration table", String::from_utf8_lossy(&s.rg)));
        let l_read_name = rec[8] as usize;
        let n_cigar = u16::from_le_bytes([rec[12], rec[13]]) as usize;
        let l_seq = r.seq.len();
        let qo = 32 + l_read_name + 4 * n_cigar + (l_seq + 1) / 2;
        self.model.recalibrate(rgk, &r.seq, &mut rec[qo..qo + l_seq], r.flag, &mut s.cov);
    }
}

#[cfg(test)]
mod tests {
    /// Mappings from GATK's own `BQSRReadTransformer.constructStaticQuantizedMapping` (tests/data/*.jsh, run with
    /// `jshell --class-path gatk.jar`).
    #[test]
    fn static_quantized_mapping_matches_gatk() {
        let text = include_str!("../../tests/data/static_quantized_mapping.txt");
        let mut n = 0;
        for l in text.lines() {
            let f: Vec<&str> = l.split(' ').collect();
            let quals: Vec<i32> = f[0].trim_matches(|c| c == '[' || c == ']').split(',').map(|x| x.parse().unwrap()).collect();
            let want: Vec<u8> = f[2].trim_end_matches(',').split(',').map(|x| x.parse().unwrap()).collect();
            assert_eq!(super::static_quantized_mapping(&quals, f[1] == "true"), want, "{l}");
            n += 1;
        }
        assert_eq!(n, 14);
    }
}
