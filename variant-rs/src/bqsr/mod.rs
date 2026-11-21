//! Base quality score recalibration: GATK 4.7 `BaseRecalibrator` (`variant-rs bqsr`), default
//! arguments, standard covariates, no BAQ, no indel tables (GATK4 defaults).

pub mod apply;
pub mod covariates;
pub mod datum;
pub mod gather;
pub mod known;
pub mod recal;
pub mod report;

use crate::input::Input;
use crate::intervals::{self, Dict};
use crate::reference::Fasta;
use covariates::ReadGroups;
use hc_preproc::Read;
use rayon::prelude::*;
use recal::{Counts, Work};
use std::sync::atomic::AtomicU64;
use std::sync::Arc;
use std::time::Instant;

struct Args {
    bam: String,
    fasta: String,
    known: Vec<String>,
    out: String,
    intervals: Vec<String>,
    threads: usize,
    chunk: i64,
    use_oq: bool,
}

const USAGE: &str = "usage: variant-rs bqsr -I in.bam|cram -R ref.fa --known-sites a.vcf.gz [--known-sites ...] -O recal.table \
[-L interval ...] [--use-original-qualities] [--threads N (16)] [--chunk BP (2000000)]";

fn parse_args(a: &[String]) -> Args {
    let mut r = Args { bam: String::new(), fasta: String::new(), known: vec![], out: String::new(), intervals: vec![], threads: 16, chunk: 2_000_000, use_oq: false };
    let mut i = 0;
    while i < a.len() {
        let v = || a.get(i + 1).cloned().unwrap_or_else(|| panic!("{} needs a value\n{USAGE}", a[i]));
        match a[i].as_str() {
            "-I" | "--input" => r.bam = v(),
            "-R" | "--reference" => r.fasta = v(),
            "--known-sites" => r.known.push(v()),
            "--use-original-qualities" | "-OQ" => {
                r.use_oq = true;
                i += 1;
                continue;
            }
            "-O" | "--output" => r.out = v(),
            "-L" | "--intervals" => r.intervals.push(v()),
            "--threads" => r.threads = v().parse().expect("--threads"),
            "--chunk" => r.chunk = v().parse().expect("--chunk"),
            x => panic!("unknown argument {x}\n{USAGE}"),
        }
        i += 2;
    }
    assert!(!r.bam.is_empty() && !r.fasta.is_empty() && !r.out.is_empty() && !r.known.is_empty(), "{USAGE}");
    r
}


/// Decoded record fields BQSR needs (plus the read group ID, and the OQ tag's qualities in `oq`: empty when the
/// read has no OQ tag or an empty one, as `ReadUtils.getOriginalBaseQualities` returns null for both).
fn decode(b: &[u8], r: &mut Read, rg: &mut Vec<u8>, oq: &mut Vec<u8>) -> (bool, bool) {
    const SEQ_NT16: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
    const OPS: &[u8; 16] = b"MIDNSHP=X???????";
    r.ref_id = i32::from_le_bytes(b[0..4].try_into().unwrap());
    r.start = i32::from_le_bytes(b[4..8].try_into().unwrap()) as i64 + 1;
    let l_read_name = b[8] as usize;
    r.mapq = b[9];
    let n_cigar = u16::from_le_bytes([b[12], b[13]]) as usize;
    r.flag = u16::from_le_bytes([b[14], b[15]]);
    let l_seq = i32::from_le_bytes(b[16..20].try_into().unwrap()) as usize;
    r.mate_ref_id = i32::from_le_bytes(b[20..24].try_into().unwrap());
    r.mate_start = i32::from_le_bytes(b[24..28].try_into().unwrap()) as i64 + 1;
    r.tlen = i32::from_le_bytes(b[28..32].try_into().unwrap()) as i64;
    let mut o = 32 + l_read_name;
    r.cigar.clear();
    let mut rl = 0i64;
    for _ in 0..n_cigar {
        let v = u32::from_le_bytes(b[o..o + 4].try_into().unwrap());
        let op = OPS[(v & 0xf) as usize];
        if matches!(op, b'M' | b'D' | b'N' | b'=' | b'X') {
            rl += (v >> 4) as i64;
        }
        r.cigar.push((op, v >> 4));
        o += 4;
    }
    r.end = r.start + rl - 1;
    r.seq.clear();
    for i in 0..l_seq / 2 {
        let byte = b[o + i];
        r.seq.push(SEQ_NT16[(byte >> 4) as usize]);
        r.seq.push(SEQ_NT16[(byte & 0xf) as usize]);
    }
    if l_seq % 2 == 1 {
        r.seq.push(SEQ_NT16[(b[o + l_seq / 2] >> 4) as usize]);
    }
    o += (l_seq + 1) / 2;
    r.qual.clear();
    r.qual.extend_from_slice(&b[o..o + l_seq]);
    let qual_missing = l_seq > 0 && b[o] == 0xff;
    o += l_seq;
    rg.clear();
    oq.clear();
    let mut has_rg = false;
    while o + 3 <= b.len() {
        let tag = [b[o], b[o + 1]];
        let ty = b[o + 2];
        o += 3;
        let size = |t: u8| -> usize {
            match t {
                b'A' | b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                _ => 0,
            }
        };
        match ty {
            b'Z' | b'H' => {
                let s = o;
                while b[o] != 0 {
                    o += 1;
                }
                if tag == *b"RG" && ty == b'Z' {
                    has_rg = true;
                    rg.extend_from_slice(&b[s..o]);
                } else if tag == *b"OQ" && ty == b'Z' {
                    oq.extend(b[s..o].iter().map(|&c| c.wrapping_sub(33))); // SAMUtils.fastqToPhred
                }
                o += 1;
            }
            b'B' => {
                let sub = b[o];
                let n = u32::from_le_bytes(b[o + 1..o + 5].try_into().unwrap()) as usize;
                o += 5 + n * size(sub);
            }
            t => o += size(t),
        }
    }
    (has_rg, qual_missing)
}

/// `ReadUtils.resetOriginalBaseQualities`: the OQ tag's qualities replace the read's (`rec`: the record, for the
/// read name in the error).
pub(super) fn reset_original_qualities(qual: &mut [u8], oq: &[u8], rec: &[u8]) {
    if qual.len() != oq.len() {
        let name = String::from_utf8_lossy(&rec[32..32 + rec[8] as usize - 1]);
        panic!("read {name}: OQ tag has {} qualities for {} bases", oq.len(), qual.len());
    }
    qual.copy_from_slice(oq);
}

fn new_read() -> Read {
    Read {
        name: String::new(),
        flag: 0,
        start: 0,
        end: 0,
        cigar: vec![],
        seq: vec![],
        qual: vec![],
        mapq: 0,
        mate_start: 0,
        tlen: 0,
        ref_id: -1,
        mate_ref_id: -1,
        bi: Vec::new(),
        bd: Vec::new(),
        sample: 0,
        oa_other_contig: false,
        bam_rec: None,
    }
}

fn make_chunks(args: &Args, ref_names: &[String], ref_lens: &[i64], fasta: &str) -> Vec<crate::bamwrite::Chunk> {
    // traversal intervals per BAM contig (dictionary order is irrelevant: counts are summed)
    let ivs: Vec<(usize, i64, i64)> = if args.intervals.is_empty() {
        ref_lens.iter().enumerate().map(|(t, &l)| (t, 1, l)).collect()
    } else {
        let dict = Dict::from_reference(fasta);
        let mut v = Vec::new();
        for l in &args.intervals {
            v.extend(intervals::parse_l_arg(l, &dict).unwrap_or_else(|e| panic!("-L {l}: {e}")));
        }
        let v = intervals::sort_merge(v);
        let mut out = Vec::new();
        for iv in v {
            let name = &dict.names[iv.tid];
            let t = ref_names.iter().position(|n| n == name).unwrap_or_else(|| panic!("{name} not in the BAM header"));
            out.push((t, iv.start, iv.end));
        }
        out.sort_by_key(|x| (x.0, x.1));
        out
    };
    crate::bamwrite::plan_chunks(&ivs, ref_lens.len(), args.chunk, false)
}

pub fn run(argv: &[String]) {
    let args = parse_args(argv);
    rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();
    let t0 = Instant::now();
    let (input, header) = Input::open(&args.bam, &args.fasta);
    let rgs = ReadGroups::from_header(&header.text);
    if let Some(id) = rgs.missing_pl.first() {
        panic!("read group {id} has no platform (PL); GATK BaseRecalibrator requires one");
    }
    let known = known::KnownSites::load(&args.known, &header.ref_names, args.threads);
    eprintln!("known sites: {} records from {} file(s) in {:.1}s", known.n_records, args.known.len(), t0.elapsed().as_secs_f64());
    let fasta = Fasta::open(&args.fasta);
    let chunks = make_chunks(&args, &header.ref_names, &header.ref_lens, &args.fasta);
    let n_rg = rgs.names.len();
    let t1 = Instant::now();
    let io_ns = Arc::new(AtomicU64::new(0));
    let counts = chunks
        .par_iter()
        .fold(
            || (Counts::new(n_rg), Work::default(), new_read(), Vec::new(), Vec::new(), Vec::new()),
            |(mut c, mut w, mut r, mut rec, mut rg, mut oq), ch| {
                let tid = ch.tid as usize;
                let name = &header.ref_names[tid];
                let len = header.ref_lens[tid];
                // reference window: from the first possible read start to well past the chunk end
                let win_from = (ch.from - 1).max(1);
                let mut laps = crate::cpuprof::Laps::new();
                let (mut refseq, _) = fasta.fetch(name, win_from, ch.to + 20_000);
                let mut ref_from = win_from;
                let mut stream = input.stream(tid, &ch.query(), false, io_ns.clone());
                laps.lap(4);
                while {
                    let more = stream.next_record(&mut rec);
                    laps.lap(8);
                    more
                } {
                    let (has_rg, qual_missing) = decode(&rec, &mut r, &mut rg, &mut oq);
                    laps.lap(2);
                    if r.ref_id != ch.tid || r.start < ch.from || r.start > ch.to {
                        continue;
                    }
                    if !ch.overlaps(r.start, r.end) {
                        continue; // does not overlap the traversal intervals
                    }
                    if !recal::passes_filters(&r, has_rg, qual_missing, len) {
                        continue;
                    }
                    let rgk = *rgs.key_by_id.get(&rg[..]).unwrap_or_else(|| panic!("read group {} not in the header", String::from_utf8_lossy(&rg)));
                    // --use-original-qualities: BaseRecalibrationEngine.resetOriginalBaseQualities, after the
                    // read filters and before the clipping
                    if args.use_oq && !oq.is_empty() {
                        reset_original_qualities(&mut r.qual, &oq, &rec);
                    }
                    let orig_span = (r.start, r.end);
                    // a quality byte >= 128 is negative to Java: GATK fails if clipping sets it on the read
                    if r.qual.iter().any(|&q| q >= 128) {
                        if let Some(&b) = hc_preproc::bqsr_clip_first_quals(&r).and_then(|q| q.into_iter().find(|&q| q >= 128)).as_ref() {
                            panic!("Base quality score {} is invalid (read {})", b as i8, r.name);
                        }
                    }
                    if !hc_preproc::bqsr_clip(&mut r, &mut w.cigar_scratch) {
                        continue;
                    }
                    laps.lap(3);
                    if r.start < ref_from || r.end >= ref_from + refseq.len() as i64 {
                        let f = r.start.max(1);
                        let (s, _) = fasta.fetch(name, f, r.end.max(f) + 20_000);
                        refseq = s;
                        ref_from = f;
                        laps.lap(4);
                    }
                    recal::process_read(&r, rgk, tid, orig_span, &known, &refseq, ref_from, &mut c, &mut w);
                    laps.lap(5);
                }
                (c, w, r, rec, rg, oq)
            },
        )
        .map(|x| x.0)
        .reduce(
            || Counts::new(n_rg),
            |mut a, b| {
                a.add(&b);
                a
            },
        );
    eprintln!("reads processed: {} in {:.1}s ({} chunks)", counts.reads_processed, t1.elapsed().as_secs_f64(), chunks.len());
    let t2 = Instant::now();
    let fin = report::finalize(&counts, 16);
    let text = report::report_text(&fin, &rgs, 16);
    std::fs::write(&args.out, text).unwrap_or_else(|e| panic!("{}: {e}", args.out));
    eprintln!("report written in {:.1}s; total {:.1}s", t2.elapsed().as_secs_f64(), t0.elapsed().as_secs_f64());
    crate::cpuprof::report_laps(
        "BaseRecalibrator",
        &[
            (0, "input: file reads (reader threads)"),
            (1, "input: BGZF inflate (reader threads)"),
            (8, "records from the stream (copy)"),
            (2, "decode (CIGAR, bases, qualities, tags)"),
            (3, "read filters + clipping (reads skipped: their filter time goes to the stream)"),
            (4, "reference fetch"),
            (5, "covariates + known sites + counts (process_read)"),
        ],
    );
}
