//! Coordinate-sorted BAM + .bai output (`-o x.bam`), in place of `| samtools sort --write-index`.
//!
//! Every SAM line is parsed by htslib's sam_parse1 (hts/alnhts.c), the parser samtools uses on bwa's output, and
//! stored as bam_write1 writes it, so the records are the ones samtools would store. Per batch, records are split
//! into buckets of the linear genome (unplaced reads last), and each bucket's share is deflated and appended to one
//! spill file. At the end each bucket is inflated and stably sorted by (tid, pos, strand), samtools sort's order
//! (bam1_cmp_core; ties keep input order). The sorted records are cut into BGZF blocks as htslib's writer cuts them
//! (a record starts a new block when it does not fit; a full block is flushed), the blocks compressed in parallel
//! by htslib's bgzf_compress and written after the header, and the .bai built by htslib's hts_idx_push from the
//! records' virtual offsets as `samtools index` reads them (a record ending a block ends at the next block's start).
//! One thread passing every record through htslib's writer took ~25 s for a 30x lane.

use rayon::prelude::*;
use std::ffi::c_void;
use std::io::{Read, Write};
use std::os::unix::fs::FileExt;

extern "C" {
    fn ah_hdr(text: *const u8, len: usize) -> *mut c_void;
    fn ah_parse(h: *mut c_void, line: *const u8, len: usize, out: *mut u8, cap: usize, tid: *mut i32, pos: *mut i64, flag: *mut u16, end: *mut i64) -> i64;
    fn ah_open(path: *const std::ffi::c_char, mode: *const std::ffi::c_char, text: *const u8, len: usize) -> *mut c_void;
    fn ah_tell(w: *mut c_void) -> u64;
    fn ah_raw(w: *mut c_void, data: *const u8, n: usize) -> i32;
    fn ah_close(w: *mut c_void) -> i32;
    fn bgzf_compress(dst: *mut c_void, dlen: *mut usize, src: *const c_void, slen: usize, level: i32) -> i32;
    fn hts_idx_init(n: i32, fmt: i32, offset0: u64, min_shift: i32, n_lvls: i32) -> *mut c_void;
    fn hts_idx_push(idx: *mut c_void, tid: i32, beg: i64, end: i64, offset: u64, is_mapped: i32) -> i32;
    fn hts_idx_finish(idx: *mut c_void, final_offset: u64) -> i32;
    fn hts_idx_save_as(idx: *mut c_void, fn_: *const std::ffi::c_char, fnidx: *const std::ffi::c_char, fmt: i32) -> i32;
    fn hts_idx_destroy(idx: *mut c_void);
}

/// spill frame: u32 record length (top bit: unmapped), u32 tid (as samtools compares it: unsigned, -1 last),
/// i64 pos << 1 | reverse, i64 end (bam_endpos, for the index)
const FRAME: usize = 24;
const UNMAPPED: u32 = 1 << 31;
/// BGZF_BLOCK_SIZE (htslib bgzf.h): uncompressed bytes per block
const BLOCK: usize = 0xff00;
const HTS_FMT_BAI: i32 = 1;
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
                    let (mut tid, mut pos, mut flag, mut end) = (0i32, 0i64, 0u16, 0i64);
                    let mut cap = buf.capacity() - st;
                    if cap < FRAME + 64 {
                        buf.reserve(FRAME + 1024);
                        cap = buf.capacity() - st;
                    }
                    let mut n = unsafe { ah_parse(hdr.get(), line.as_ptr(), line.len(), buf.as_mut_ptr().add(st + FRAME), cap - FRAME, &mut tid, &mut pos, &mut flag, &mut end) };
                    assert!(n != -1, "htslib could not parse: {}", String::from_utf8_lossy(line));
                    assert!(n != -2, "BAM cannot store this record (bam_write1 refuses it): {}", String::from_utf8_lossy(line));
                    if n as usize > cap - FRAME {
                        buf.reserve(FRAME + n as usize);
                        cap = buf.capacity() - st;
                        n = unsafe { ah_parse(hdr.get(), line.as_ptr(), line.len(), buf.as_mut_ptr().add(st + FRAME), cap - FRAME, &mut tid, &mut pos, &mut flag, &mut end) };
                    }
                    let n = n as usize;
                    // SAFETY: ah_parse wrote n bytes after the frame, within the capacity
                    unsafe { buf.set_len(st + FRAME + n) };
                    let unmapped = if flag & 0x4 != 0 { UNMAPPED } else { 0 };
                    buf[st..st + 4].copy_from_slice(&(n as u32 | unmapped).to_le_bytes());
                    buf[st + 4..st + 8].copy_from_slice(&(tid as u32).to_le_bytes());
                    let k = (pos << 1) | ((flag & 0x10) >> 4) as i64;
                    buf[st + 8..st + 16].copy_from_slice(&k.to_le_bytes());
                    buf[st + 16..st + 24].copy_from_slice(&end.to_le_bytes());
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

    /// Sorts the buckets and writes the BAM and its .bai. Four steps overlap, window by window (WINDOW buckets):
    /// read back and sort; cut into blocks and compress them (in parallel); write them; index the records.
    pub fn finish(self) {
        let path = std::ffi::CString::new(self.path.clone()).unwrap();
        let mode = std::ffi::CString::new(format!("wb{}", self.level)).unwrap();
        let w = Ptr(unsafe { ah_open(path.as_ptr(), mode.as_ptr(), self.text.as_ptr(), self.text.len()) });
        assert!(!w.0.is_null(), "{}: could not open for writing", self.path);
        let off0 = unsafe { ah_tell(w.0) };
        let idx = Ptr(unsafe { hts_idx_init(self.cum.len() as i32, HTS_FMT_BAI, off0, 14, 5) });
        assert!(!idx.0.is_null(), "{}: hts_idx_init failed", self.path);
        // compressed address of the next block
        let mut addr = off0 >> 16;
        let level = self.level as i32;
        let me = &self;
        // one bucket: its records sorted (buffer, and per record where it is and what the index needs). Offsets are
        // 64-bit: a bucket can pass 4 GiB (all unplaced reads share one)
        let sort_bucket = |b: usize| -> (Vec<u8>, Vec<Rec>) {
            let raw: usize = me.chunks[b].iter().map(|c| c.2 as usize).sum();
            let mut buf = Vec::with_capacity(raw);
            for &(off, clen, _) in &me.chunks[b] {
                let mut z = vec![0u8; clen as usize];
                me.spill.read_exact_at(&mut z, off).unwrap();
                flate2::read::DeflateDecoder::new(&z[..]).read_to_end(&mut buf).unwrap();
            }
            let mut recs: Vec<Rec> = Vec::new();
            let mut p = 0usize;
            while p < buf.len() {
                let n = u32::from_le_bytes(buf[p..p + 4].try_into().unwrap());
                let r = Rec {
                    at: (p + FRAME) as u64,
                    n: n & !UNMAPPED,
                    mapped: n & UNMAPPED == 0,
                    tid: u32::from_le_bytes(buf[p + 4..p + 8].try_into().unwrap()),
                    k: i64::from_le_bytes(buf[p + 8..p + 16].try_into().unwrap()),
                    end: i64::from_le_bytes(buf[p + 16..p + 24].try_into().unwrap()),
                };
                p += FRAME + r.n as usize;
                recs.push(r);
            }
            recs.sort_by_key(|r| (r.tid, r.k)); // stable: ties keep input order
            (buf, recs)
        };
        std::thread::scope(|sc| {
            let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<(Vec<u8>, Vec<Rec>)>>(1);
            let (wtx, wrx) = std::sync::mpsc::sync_channel::<Vec<Vec<u8>>>(1);
            let (itx, irx) = std::sync::mpsc::sync_channel::<Vec<(Rec, u64)>>(1);
            sc.spawn(move || {
                let nb = me.chunks.len();
                for w0 in (0..nb).step_by(WINDOW) {
                    let v: Vec<_> = (w0..(w0 + WINDOW).min(nb)).into_par_iter().map(sort_bucket).collect();
                    if tx.send(v).is_err() {
                        break;
                    }
                }
            });
            // the blocks, in one write a window
            sc.spawn(move || {
                let w = w;
                for z in wrx {
                    let all = z.concat();
                    drop(z);
                    let rc = unsafe { ah_raw(w.0, all.as_ptr(), all.len()) };
                    assert_eq!(rc, 0, "{}: write failed", me.path);
                }
            });
            // the index, record by record in file order
            sc.spawn(move || {
                let idx = idx;
                for e in irx {
                    for (r, off) in e {
                        let rc = unsafe { hts_idx_push(idx.0, r.tid as i32, r.k >> 1, r.end, off, r.mapped as i32) };
                        assert!(rc >= 0, "{}: hts_idx_push failed", me.path);
                    }
                }
            });
            for v in rx {
                // the window's records cut into blocks as bgzf_write / bgzf_flush_try cut them: a record that does
                // not fit the block's room starts a new one, a full block is flushed (a record longer than a block
                // spans several). The window's last block is flushed with it.
                let mut pieces: Vec<(u32, u64, u32)> = Vec::new(); // (bucket in the window, start, length)
                let mut blocks: Vec<(usize, usize, usize)> = Vec::new(); // (first piece, end piece, length)
                let mut ends: Vec<(u32, u32)> = Vec::new(); // per record: (block, offset in it after the record)
                let (mut p0, mut len) = (0usize, 0usize);
                for (bi, (_, recs)) in v.iter().enumerate() {
                    for r in recs {
                        let n = r.n as usize;
                        if len > 0 && len + n > BLOCK {
                            blocks.push((p0, pieces.len(), len));
                            (p0, len) = (pieces.len(), 0);
                        }
                        let (mut src, mut rem) = (r.at, n);
                        while rem > 0 {
                            let take = (BLOCK - len).min(rem);
                            pieces.push((bi as u32, src, take as u32));
                            (src, rem, len) = (src + take as u64, rem - take, len + take);
                            if len == BLOCK {
                                blocks.push((p0, pieces.len(), len));
                                (p0, len) = (pieces.len(), 0);
                            }
                        }
                        // a record that filled its block ends at that block's end
                        ends.push(if len == 0 { (blocks.len() as u32 - 1, BLOCK as u32) } else { (blocks.len() as u32, len as u32) });
                    }
                }
                if len > 0 {
                    blocks.push((p0, pieces.len(), len));
                }
                let z: Vec<Vec<u8>> = blocks
                    .par_iter()
                    .map(|&(a, b, l)| {
                        let mut raw = Vec::with_capacity(l);
                        for &(bi, s, n) in &pieces[a..b] {
                            raw.extend_from_slice(&v[bi as usize].0[s as usize..s as usize + n as usize]);
                        }
                        let mut out = vec![0u8; 0x10000];
                        let mut dlen = out.len();
                        let rc = unsafe { bgzf_compress(out.as_mut_ptr() as *mut c_void, &mut dlen, raw.as_ptr() as *const c_void, raw.len(), level) };
                        assert_eq!(rc, 0, "bgzf_compress failed");
                        out.truncate(dlen);
                        out
                    })
                    .collect();
                // block addresses, then each record's virtual offset as a reader's bgzf_tell gives it after the record
                let mut baddr = Vec::with_capacity(z.len() + 1);
                for zb in &z {
                    baddr.push(addr);
                    addr += zb.len() as u64;
                }
                baddr.push(addr);
                let e: Vec<(Rec, u64)> = v
                    .iter()
                    .flat_map(|(_, recs)| recs.iter())
                    .zip(&ends)
                    .map(|(r, &(b, o))| (*r, if o as usize == blocks[b as usize].2 { baddr[b as usize + 1] << 16 } else { baddr[b as usize] << 16 | o as u64 }))
                    .collect();
                drop(v);
                wtx.send(z).unwrap();
                itx.send(e).unwrap();
            }
            drop(wtx);
            drop(itx);
        });
        assert_eq!(unsafe { hts_idx_finish(idx.0, addr << 16) }, 0, "{}: hts_idx_finish failed", self.path);
        assert_eq!(unsafe { ah_close(w.0) }, 0, "{}: closing failed", self.path);
        let bai = std::ffi::CString::new(format!("{}.bai", self.path)).unwrap();
        assert_eq!(unsafe { hts_idx_save_as(idx.0, std::ptr::null(), bai.as_ptr(), HTS_FMT_BAI) }, 0, "{}: writing the index failed", self.path);
        unsafe { hts_idx_destroy(idx.0) };
        let _ = std::fs::remove_file(&self.spill_path);
    }
}

/// A sorted record: where its bytes are in the bucket's buffer, and what the sort and the index need
#[derive(Clone, Copy)]
struct Rec {
    at: u64,
    n: u32,
    mapped: bool,
    /// as samtools compares it: unsigned, -1 last
    tid: u32,
    /// pos << 1 | reverse
    k: i64,
    end: i64,
}

/// the BAM writer or the index, on the thread that owns it
#[derive(Clone, Copy)]
struct Ptr(*mut c_void);
unsafe impl Send for Ptr {}
