//! SortSam (`variant-rs sortsam`): Picard `SortSam SORT_ORDER=coordinate`. htsjdk's sorting writer
//! (`SAMFileWriterImpl` + `SortingCollection`) orders records with `SAMRecordCoordinateComparator` stably: chunks
//! of `MAX_RECORDS_IN_RAM` sorted with TimSort, spilled, merged with ties to the earlier chunk. Any stable sort by
//! that comparator gives the same order, so the chunk size here (a memory budget) cannot change the output.
//! Records are written with their bytes as read (htsjdk re-encodes only records whose attributes were set);
//! the header is htsjdk's rewrite of the input's with `SO:coordinate`.
//!
//! ```text
//! input (BGZF, parallel inflate) -> chunks of records -> stable parallel sort -> spill (BGZF level 1, temp dir)
//!   -> k-way merge -> batches per reference -> parallel deflate -> BAM (+ .bai, .md5)
//! ```

use crate::bamwrite::{self, BamFile, ChunkOut, RecWriter};
use crate::bgzf::ChunkStream;
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::fs::File;
use std::sync::atomic::AtomicU64;
use std::sync::Arc;
use std::time::Instant;

const USAGE: &str = "usage: variant-rs sortsam -I in.bam -O out.bam [--SORT_ORDER coordinate] [--CREATE_INDEX true] \
[--CREATE_MD5_FILE true] [--tmp-dir dir] [--compression-level 2] [--sort-buffer-mb 4096] [--threads n]";

struct Args {
    input: String,
    out: String,
    index: bool,
    md5: bool,
    tmp: Option<String>,
    level: i32,
    buffer: usize,
    threads: usize,
}

fn bool_arg(v: &str) -> bool {
    match v.to_ascii_lowercase().as_str() {
        "true" | "t" | "yes" | "y" => true,
        "false" | "f" | "no" | "n" => false,
        _ => panic!("not a boolean: {v}"),
    }
}

fn parse_args(a: &[String]) -> Args {
    let mut r = Args {
        input: String::new(),
        out: String::new(),
        index: false,
        md5: false,
        tmp: None,
        level: 2,
        buffer: 4096 << 20,
        threads: std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1),
    };
    let mut i = 0;
    while i < a.len() {
        let v = || a.get(i + 1).cloned().unwrap_or_else(|| panic!("{} needs a value\n{USAGE}", a[i]));
        match a[i].as_str() {
            "-I" | "--INPUT" => r.input = v(),
            "-O" | "--OUTPUT" => r.out = v(),
            "-SO" | "--SORT_ORDER" => assert_eq!(v(), "coordinate", "variant-rs SortSam sorts by coordinate only"),
            "--CREATE_INDEX" => r.index = bool_arg(&v()),
            "--CREATE_MD5_FILE" => r.md5 = bool_arg(&v()),
            "--tmp-dir" | "--TMP_DIR" => r.tmp = Some(v()),
            "--compression-level" => r.level = v().parse().expect("--compression-level"),
            "--sort-buffer-mb" => r.buffer = v().parse::<usize>().expect("--sort-buffer-mb") << 20,
            "--threads" => r.threads = v().parse().expect("--threads"),
            x => panic!("unknown argument {x}\n{USAGE}"),
        }
        i += 2;
    }
    assert!(!r.input.is_empty() && !r.out.is_empty(), "{USAGE}");
    r
}

#[inline]
fn i32_at(b: &[u8], o: usize) -> i32 {
    i32::from_le_bytes(b[o..o + 4].try_into().unwrap())
}

/// htsjdk `SAMRecordCoordinateComparator` on two BAM records (bytes after `block_size`).
pub fn coord_cmp(a: &[u8], b: &[u8]) -> Ordering {
    // fileOrderCompare: unplaced reads (reference index -1) last and equal to each other
    let (ra, rb) = (i32_at(a, 0), i32_at(b, 0));
    let c = match (ra == -1, rb == -1) {
        (true, true) => Ordering::Equal,
        (true, false) => Ordering::Greater,
        (false, true) => Ordering::Less,
        _ => ra.cmp(&rb).then(i32_at(a, 4).cmp(&i32_at(b, 4))),
    };
    if c != Ordering::Equal {
        return c;
    }
    let (fa, fb) = (u16::from_le_bytes([a[14], a[15]]), u16::from_le_bytes([b[14], b[15]]));
    let (rev_a, rev_b) = (fa & 0x10 != 0, fb & 0x10 != 0);
    if rev_a != rev_b {
        return if rev_a { Ordering::Greater } else { Ordering::Less };
    }
    // read name (without its NUL), flags, MAPQ, mate reference, mate start, insert size
    let (na, nb) = (&a[32..32 + a[8] as usize - 1], &b[32..32 + b[8] as usize - 1]);
    na.cmp(nb)
        .then(fa.cmp(&fb))
        .then(a[9].cmp(&b[9]))
        .then(i32_at(a, 20).cmp(&i32_at(b, 20)))
        .then(i32_at(a, 24).cmp(&i32_at(b, 24)))
        .then(i32_at(a, 28).cmp(&i32_at(b, 28)))
}

/// The header htsjdk writes: `header_with_pg` (no new @PG) with `SO:coordinate` in place of the input's `SO`
/// (`SAMFileHeader.setSortOrder`: same position, or appended).
fn header_text(input: &str) -> String {
    let t = bamwrite::header_with_pg(input, &[]);
    let mut out = String::new();
    for (k, l) in t.lines().enumerate() {
        if k == 0 {
            let mut f: Vec<String> = l.split('\t').map(|s| s.to_string()).collect();
            match f.iter_mut().find(|x| x.starts_with("SO:")) {
                Some(x) => *x = "SO:coordinate".into(),
                None => f.push("SO:coordinate".into()),
            }
            out.push_str(&f.join("\t"));
        } else {
            out.push_str(l);
        }
        out.push('\n');
    }
    out
}

/// One in-memory chunk: record bytes back to back, and (offset, length) per record.
#[derive(Default)]
struct Chunk {
    data: Vec<u8>,
    recs: Vec<(usize, u32)>,
}

impl Chunk {
    fn rec(&self, k: usize) -> &[u8] {
        let (o, l) = self.recs[k];
        &self.data[o..o + l as usize]
    }
    fn sort(&mut self) {
        let data = &self.data;
        // rayon's par_sort_by is stable (TimSort's order on ties)
        self.recs.par_sort_by(|&(oa, la), &(ob, lb)| coord_cmp(&data[oa..oa + la as usize], &data[ob..ob + lb as usize]));
    }
}

/// Merged output side: records in batches of one reference, deflated in parallel, appended in order.
struct Out {
    bam: BamFile,
    level: i32,
    batches: Vec<(i32, Chunk)>,
    cur: Chunk,
    cur_tid: i32,
}

const BATCH: usize = 8 << 20;
const GROUP: usize = 64;

impl Out {
    fn push(&mut self, rec: &[u8]) {
        let tid = i32_at(rec, 0);
        if tid != self.cur_tid || self.cur.data.len() >= BATCH {
            self.cut();
            self.cur_tid = tid;
        }
        self.cur.recs.push((self.cur.data.len(), rec.len() as u32));
        self.cur.data.extend_from_slice(rec);
    }
    fn cut(&mut self) {
        if !self.cur.recs.is_empty() {
            self.batches.push((self.cur_tid, std::mem::take(&mut self.cur)));
        }
        if self.batches.len() >= GROUP {
            self.flush();
        }
    }
    fn flush(&mut self) {
        let level = self.level;
        let big = self.bam.too_large.clone();
        let outs: Vec<(i32, ChunkOut)> = std::mem::take(&mut self.batches)
            .into_par_iter()
            .map(|(tid, c)| {
                let mut w = RecWriter::new(level, &big);
                for k in 0..c.recs.len() {
                    w.push(c.rec(k));
                }
                (tid, w.finish())
            })
            .collect();
        for (tid, o) in &outs {
            self.bam.append(*tid, o);
        }
    }
    fn finish(mut self, path: &str) -> u64 {
        self.cut();
        self.flush();
        let n = self.bam.n_records;
        self.bam.finish(path);
        n
    }
}

/// A spilled chunk being merged: its next record.
struct Head {
    rec: Vec<u8>,
    src: usize,
}

impl PartialEq for Head {
    fn eq(&self, o: &Self) -> bool {
        self.cmp(o) == Ordering::Equal
    }
}
impl Eq for Head {}
impl PartialOrd for Head {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Head {
    // reversed: BinaryHeap pops the greatest; ties to the earlier chunk
    fn cmp(&self, o: &Self) -> Ordering {
        coord_cmp(&o.rec, &self.rec).then(o.src.cmp(&self.src))
    }
}

pub fn run(argv: &[String]) {
    let a = parse_args(argv);
    let t0 = Instant::now();
    rayon::ThreadPoolBuilder::new().num_threads(a.threads).build_global().ok();
    let inp = crate::markdup::qname::open_input(&a.input);
    let header = header_text(&inp.header.text);
    let tmp_dir = std::path::PathBuf::from(a.tmp.clone().unwrap_or_else(|| {
        std::path::Path::new(&a.out).parent().map(|p| p.to_string_lossy().into_owned()).filter(|p| !p.is_empty()).unwrap_or_else(|| ".".into())
    }));
    let io_ns = Arc::new(AtomicU64::new(0));
    let mut stream = ChunkStream::new(inp.file.clone(), vec![(inp.first, inp.len << 16)], true, io_ns.clone());

    // read and sort chunks; spill all but the last
    let mut spills: Vec<(std::path::PathBuf, u64)> = Vec::new();
    let mut chunk = Chunk::default();
    let mut body = Vec::new();
    let mut n_in = 0u64;
    let spill = |c: &Chunk, k: usize| -> (std::path::PathBuf, u64) {
        let p = tmp_dir.join(format!("variant-rs-sortsam.{}.{k}.tmp", std::process::id()));
        // compress in parallel pieces; a spill file is BGZF blocks of records without a header
        let pieces: Vec<ChunkOut> = c
            .recs
            .par_chunks(1 << 16)
            .map(|rs| {
                let mut w = RecWriter::new(1, &crate::bamwrite::too_large(&[]));
                for &(o, l) in rs {
                    w.push(&c.data[o..o + l as usize]);
                }
                w.finish()
            })
            .collect();
        let mut f = std::io::BufWriter::new(File::create(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display())));
        let mut len = 0u64;
        for o in &pieces {
            std::io::Write::write_all(&mut f, &o.data).unwrap();
            len += o.data.len() as u64;
        }
        std::io::Write::flush(&mut f).unwrap();
        (p, len)
    };
    while stream.next_record(&mut body) {
        n_in += 1;
        chunk.recs.push((chunk.data.len(), body.len() as u32));
        chunk.data.extend_from_slice(&body);
        if chunk.data.len() + chunk.recs.len() * 16 >= a.buffer {
            chunk.sort();
            spills.push(spill(&chunk, spills.len()));
            chunk = Chunk::default();
        }
    }
    chunk.sort();
    let t_read = t0.elapsed().as_secs_f64();

    let bam = BamFile::create(&a.out, &header, &inp.header.ref_names, &inp.header.ref_lens, a.level, a.md5);
    let mut out = Out { bam, level: a.level, batches: Vec::new(), cur: Chunk::default(), cur_tid: i32::MIN };
    out.bam.index = a.index;
    if spills.is_empty() {
        for k in 0..chunk.recs.len() {
            out.push(chunk.rec(k));
        }
    } else {
        // the in-memory chunk is the last source (the latest records: ties go to the spilled ones first)
        let mut streams: Vec<ChunkStream> = spills
            .iter()
            .map(|(p, len)| ChunkStream::new(Arc::new(File::open(p).unwrap()), vec![(0, len << 16)], true, io_ns.clone()))
            .collect();
        let mem_src = streams.len();
        let mut mem_next = 0usize;
        let mut heap = BinaryHeap::new();
        for (k, s) in streams.iter_mut().enumerate() {
            let mut rec = Vec::new();
            if s.next_record(&mut rec) {
                heap.push(Head { rec, src: k });
            }
        }
        if mem_next < chunk.recs.len() {
            heap.push(Head { rec: chunk.rec(0).to_vec(), src: mem_src });
            mem_next = 1;
        }
        while let Some(mut h) = heap.pop() {
            out.push(&h.rec);
            let more = if h.src == mem_src {
                if mem_next < chunk.recs.len() {
                    h.rec.clear();
                    h.rec.extend_from_slice(chunk.rec(mem_next));
                    mem_next += 1;
                    true
                } else {
                    false
                }
            } else {
                streams[h.src].next_record(&mut h.rec)
            };
            if more {
                heap.push(h);
            }
        }
        for (p, _) in &spills {
            let _ = std::fs::remove_file(p);
        }
    }
    let n_out = out.finish(&a.out);
    assert_eq!(n_in, n_out, "records in / out");
    eprintln!(
        "sortsam: {n_out} records, {} spill file(s), read+sort {:.1}s, total {:.1}s",
        spills.len(),
        t_read,
        t0.elapsed().as_secs_f64()
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn header_so() {
        let h = header_text("@HD\tVN:1.5\tSO:queryname\tGO:query\n@SQ\tSN:c\tLN:5\n");
        assert_eq!(h, "@HD\tVN:1.6\tSO:coordinate\tGO:query\n@SQ\tSN:c\tLN:5\n");
        assert_eq!(header_text("@SQ\tSN:c\tLN:5\n"), "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:5\n");
    }
}
