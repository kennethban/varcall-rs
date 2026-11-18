//! Coordinate-sorted BAM + .bai output (`-o x.bam`), in place of `| samtools sort --write-index`.
//!
//! Every SAM line is parsed by htslib's sam_parse1 (hts/alnhts.c), the parser samtools uses on bwa's output, so
//! the records are the ones samtools would store. Per batch, records are split into buckets of the linear genome
//! (unplaced reads last), and each bucket's share is deflated and appended to one spill file. At the end each
//! bucket is inflated and stably sorted by (tid, pos, strand), samtools sort's order (bam1_cmp_core; ties keep
//! input order), then written through htslib's BGZF writer with the index built on the fly (sam_idx_init), as
//! `samtools sort --write-index` does.

use rayon::prelude::*;
use std::ffi::c_void;
use std::io::{Read, Write};
use std::os::unix::fs::FileExt;

extern "C" {
    fn ah_hdr(text: *const u8, len: usize) -> *mut c_void;
    fn ah_parse(h: *mut c_void, line: *const u8, len: usize, out: *mut u8, cap: usize, tid: *mut i32, pos: *mut i64, flag: *mut u16) -> i64;
    fn ah_open(path: *const std::ffi::c_char, mode: *const std::ffi::c_char, text: *const u8, len: usize, threads: i32) -> *mut c_void;
    fn ah_write(w: *mut c_void, blob: *const u8, n: usize) -> i32;
    fn ah_close(w: *mut c_void) -> i32;
}

/// spill frame: u32 blob length, u32 tid (as samtools compares it: unsigned, -1 last), i64 pos << 1 | reverse
const FRAME: usize = 16;
const N_BUCKETS: i64 = 1024;
/// buckets sorted at a time in the final pass
const WINDOW: usize = 64;

#[derive(Clone, Copy)]
struct Hdr(*mut c_void);
unsafe impl Send for Hdr {}
unsafe impl Sync for Hdr {}
impl Hdr {
    fn get(&self) -> *mut c_void {
        self.0
    }
}

pub struct BamSink {
    hdr: Hdr,
    text: String,
    path: String,
    spill_path: String,
    spill: std::fs::File,
    spill_len: u64,
    /// per bucket: (offset, deflated length, raw length) of its pieces, in input order
    chunks: Vec<Vec<(u64, u64, u64)>>,
    /// linear genome offset of each reference
    cum: Vec<i64>,
    bsize: i64,
    level: u32,
    pub n_records: u64,
}

impl BamSink {
    /// `text`: the SAM header; `ref_lens`: reference lengths in header order
    pub fn new(path: &str, text: &str, ref_lens: &[i64], level: u32) -> BamSink {
        let hdr = unsafe { ah_hdr(text.as_ptr(), text.len()) };
        assert!(!hdr.is_null(), "htslib could not parse the SAM header");
        let mut cum = Vec::with_capacity(ref_lens.len() + 1);
        let mut acc = 0i64;
        for &l in ref_lens {
            cum.push(acc);
            acc += l;
        }
        let spill_path = format!("{path}.spill.tmp");
        let spill = std::fs::OpenOptions::new().read(true).write(true).create(true).truncate(true).open(&spill_path).unwrap_or_else(|e| panic!("{spill_path}: {e}"));
        BamSink {
            hdr: Hdr(hdr),
            text: text.to_string(),
            path: path.to_string(),
            spill_path,
            spill,
            spill_len: 0,
            chunks: vec![Vec::new(); N_BUCKETS as usize + 1],
            cum,
            bsize: acc / N_BUCKETS + 1,
            level,
            n_records: 0,
        }
    }

    pub fn path(&self) -> &str {
        &self.path
    }

    fn bucket(&self, tid: i32, pos: i64) -> usize {
        if tid < 0 {
            N_BUCKETS as usize
        } else {
            (((self.cum[tid as usize] + pos.max(0)) / self.bsize) as usize).min(N_BUCKETS as usize - 1)
        }
    }

    /// One batch of SAM text (whole lines), in output order
    pub fn add_batch(&mut self, text: &[Vec<u8>]) {
        let hdr = self.hdr;
        // per text: framed records and (bucket, start, length) of each
        let parsed: Vec<(Vec<u8>, Vec<(u32, u64, u32)>)> = text
            .par_iter()
            .map(|t| crate::prof::cpu(crate::prof::C_BPARSE, || {
                let mut buf: Vec<u8> = Vec::with_capacity(t.len() + 256);
                let mut recs = Vec::new();
                for line in t.split(|&c| c == b'\n').filter(|l| !l.is_empty()) {
                    let st = buf.len();
                    let (mut tid, mut pos, mut flag) = (0i32, 0i64, 0u16);
                    let mut cap = buf.capacity() - st;
                    if cap < FRAME + 64 {
                        buf.reserve(FRAME + 1024);
                        cap = buf.capacity() - st;
                    }
                    let mut n = unsafe { ah_parse(hdr.get(), line.as_ptr(), line.len(), buf.as_mut_ptr().add(st + FRAME), cap - FRAME, &mut tid, &mut pos, &mut flag) };
                    assert!(n >= 0, "htslib could not parse: {}", String::from_utf8_lossy(line));
                    if n as usize > cap - FRAME {
                        buf.reserve(FRAME + n as usize);
                        cap = buf.capacity() - st;
                        n = unsafe { ah_parse(hdr.get(), line.as_ptr(), line.len(), buf.as_mut_ptr().add(st + FRAME), cap - FRAME, &mut tid, &mut pos, &mut flag) };
                    }
                    let n = n as usize;
                    // SAFETY: ah_parse wrote n bytes after the frame, within the capacity
                    unsafe { buf.set_len(st + FRAME + n) };
                    buf[st..st + 4].copy_from_slice(&(n as u32).to_le_bytes());
                    buf[st + 4..st + 8].copy_from_slice(&(tid as u32).to_le_bytes());
                    let k = (pos << 1) | ((flag & 0x10) >> 4) as i64;
                    buf[st + 8..st + 16].copy_from_slice(&k.to_le_bytes());
                    recs.push((self.bucket(tid, pos) as u32, st as u64, (FRAME + n) as u32));
                }
                (buf, recs)
            }))
            .collect();
        // records by bucket, input order kept within each
        let mut all: Vec<(u32, u32, u64, u32)> = parsed.iter().enumerate().flat_map(|(i, (_, r))| r.iter().map(move |&(b, s, l)| (b, i as u32, s, l))).collect();
        all.par_sort_by_key(|x| x.0);
        self.n_records += all.len() as u64;
        let mut groups = Vec::new();
        let mut s = 0;
        while s < all.len() {
            let mut e = s;
            while e < all.len() && all[e].0 == all[s].0 {
                e += 1;
            }
            groups.push((s, e));
            s = e;
        }
        let level = flate2::Compression::new(1);
        let pieces: Vec<(usize, u64, Vec<u8>)> = groups
            .par_iter()
            .map(|&(s, e)| crate::prof::cpu(crate::prof::C_BDEFL, || {
                let raw: u64 = all[s..e].iter().map(|x| x.3 as u64).sum();
                let mut enc = flate2::write::DeflateEncoder::new(Vec::with_capacity(raw as usize / 3), level);
                for &(_, i, st, l) in &all[s..e] {
                    enc.write_all(&parsed[i as usize].0[st as usize..st as usize + l as usize]).unwrap();
                }
                (all[s].0 as usize, raw, enc.finish().unwrap())
            }))
            .collect();
        drop(all);
        drop(parsed);
        for (b, raw, z) in pieces {
            self.spill.write_all_at(&z, self.spill_len).unwrap_or_else(|e| panic!("{}: {e}", self.spill_path));
            self.chunks[b].push((self.spill_len, z.len() as u64, raw));
            self.spill_len += z.len() as u64;
        }
    }

    /// Sorts the buckets and writes the BAM and its .bai
    pub fn finish(self, threads: usize) {
        let path = std::ffi::CString::new(self.path.clone()).unwrap();
        let mode = std::ffi::CString::new(format!("wb{}", self.level)).unwrap();
        let w = unsafe { ah_open(path.as_ptr(), mode.as_ptr(), self.text.as_ptr(), self.text.len(), threads as i32) };
        assert!(!w.is_null(), "{}: could not open for writing", self.path);
        let me = &self;
        // one bucket: its records sorted, as (buffer, (start, length) of each blob). Offsets are 64-bit: a bucket can
        // pass 4 GiB (all unplaced reads share one)
        let sort_bucket = |b: usize| -> (Vec<u8>, Vec<(u64, u32)>) {
            let raw: usize = me.chunks[b].iter().map(|c| c.2 as usize).sum();
            let mut buf = Vec::with_capacity(raw);
            for &(off, clen, _) in &me.chunks[b] {
                let mut z = vec![0u8; clen as usize];
                me.spill.read_exact_at(&mut z, off).unwrap();
                flate2::read::DeflateDecoder::new(&z[..]).read_to_end(&mut buf).unwrap();
            }
            let mut recs: Vec<(u32, i64, u64, u32)> = Vec::new();
            let mut p = 0usize;
            while p < buf.len() {
                let n = u32::from_le_bytes(buf[p..p + 4].try_into().unwrap());
                let tid = u32::from_le_bytes(buf[p + 4..p + 8].try_into().unwrap());
                let k = i64::from_le_bytes(buf[p + 8..p + 16].try_into().unwrap());
                recs.push((tid, k, (p + FRAME) as u64, n));
                p += FRAME + n as usize;
            }
            recs.sort_by_key(|r| (r.0, r.1)); // stable: ties keep input order
            (buf, recs.into_iter().map(|r| (r.2, r.3)).collect())
        };
        std::thread::scope(|sc| {
            let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<(Vec<u8>, Vec<(u64, u32)>)>>(1);
            sc.spawn(move || {
                let nb = me.chunks.len();
                for w0 in (0..nb).step_by(WINDOW) {
                    let v: Vec<_> = (w0..(w0 + WINDOW).min(nb)).into_par_iter().map(sort_bucket).collect();
                    if tx.send(v).is_err() {
                        break;
                    }
                }
            });
            for v in rx {
                for (buf, recs) in v {
                    for (s, n) in recs {
                        let rc = unsafe { ah_write(w, buf.as_ptr().add(s as usize), n as usize) };
                        assert_eq!(rc, 0, "{}: write failed", me.path);
                    }
                }
            }
        });
        assert_eq!(unsafe { ah_close(w) }, 0, "{}: closing / writing the index failed", self.path);
        let _ = std::fs::remove_file(&self.spill_path);
    }
}
