//! Read input: BAM + BAI (own BGZF reader) or CRAM + CRAI (htslib, `cram.rs`). Both yield BAM record bytes.

use crate::bam::{Bai, Header};
use crate::bgzf::ChunkStream;
use crate::cram::{self, CramPool, RecordStream};
use std::fs::File;
use std::io::Read;
use std::sync::atomic::AtomicU64;
use std::sync::Arc;

pub enum Input {
    Bam { file: Arc<File>, bai: Bai },
    Cram { pool: Arc<CramPool>, crai: String },
}

pub enum Stream {
    Bam(ChunkStream),
    Cram(RecordStream),
}

impl Stream {
    /// Next record (the bytes after `block_size`) into `out`. False at the end.
    pub fn next_record(&mut self, out: &mut Vec<u8>) -> bool {
        match self {
            Stream::Bam(s) => s.next_record(out),
            Stream::Cram(s) => s.next_record(out),
        }
    }
}

impl Input {
    /// Opens `path` as BAM or CRAM (by content); CRAM needs the reference `fasta`.
    pub fn open(path: &str, fasta: &str) -> (Input, Header) {
        let mut magic = [0u8; 4];
        File::open(path).and_then(|mut f| f.read_exact(&mut magic)).unwrap_or_else(|e| panic!("{path}: {e}"));
        if &magic == b"CRAM" {
            let pool = Arc::new(CramPool::open(path, fasta));
            let (text, ref_names, ref_lens) = pool.header();
            (Input::Cram { pool, crai: cram::crai_path(path) }, Header { text, ref_names, ref_lens })
        } else {
            let file = Arc::new(File::open(path).unwrap_or_else(|e| panic!("{path}: {e}")));
            let header = Header::read(&file);
            // no .bai (as nf-core/sarek hands MarkDuplicates each FASTQ chunk's BAM): one is built in memory
            let bai = match crate::plan::find_bai(path) {
                Some(p) => Bai::load(&p).unwrap(),
                None => {
                    let t = std::time::Instant::now();
                    let b = Bai::from_bytes(&crate::gatherbam::bai_of(path));
                    eprintln!("{path}: no .bai; indexed in memory in {:.1}s", t.elapsed().as_secs_f64());
                    b
                }
            };
            (Input::Bam { file, bai }, header)
        }
    }

    /// Records of `tid` that may overlap the (padded, merged) intervals, in file order. Records outside
    /// the intervals or on other contigs may be included; the caller filters as htsjdk does.
    pub fn stream(&self, tid: usize, intervals: &[(i64, i64)], parallel: bool, io_ns: Arc<AtomicU64>) -> Stream {
        match self {
            Input::Bam { file, bai } => Stream::Bam(ChunkStream::new(file.clone(), bai.query(tid, intervals), parallel, io_ns)),
            Input::Cram { pool, .. } => Stream::Cram(pool.stream(tid, intervals, io_ns)),
        }
    }

    /// Whether the index points to any data for `tid` over `[from, to]` (CRAM: assumed yes).
    pub fn has_data(&self, tid: usize, from: i64, to: i64) -> bool {
        match self {
            Input::Bam { bai, .. } => !bai.query(tid, &[(from, to)]).is_empty(),
            Input::Cram { .. } => true,
        }
    }

    /// Estimated read data per 16 kb window, per reference (index-only; no reads decoded).
    pub fn window_costs(&self, lens: &[i64]) -> Vec<Vec<f64>> {
        match self {
            Input::Bam { bai, .. } => lens.iter().enumerate().map(|(t, &l)| crate::plan::window_costs(bai, t, l)).collect(),
            Input::Cram { crai, .. } => cram::crai_window_costs(crai, lens),
        }
    }
}
