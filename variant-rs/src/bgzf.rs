//! BGZF input (random access by virtual offset, parallel inflate) and output (parallel deflate,
//! block bookkeeping for the tabix index).

use libdeflater::{CompressionLvl, Compressor, Decompressor};
use rayon::prelude::*;
use std::cell::RefCell;
use std::fs::File;
use std::io::{self, Write};
use std::os::unix::fs::FileExt;
use std::sync::mpsc::{Receiver, SyncSender};
use std::sync::Arc;
use std::time::Instant;

/// Uncompressed bytes per output block (what htslib/htsjdk use: 0xff00).
pub const BLOCK_DATA: usize = 0xff00;
/// The empty BGZF block that terminates a file.
pub const EOF_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

thread_local! {
    static DECOMP: RefCell<Decompressor> = RefCell::new(Decompressor::new());
}

/// One inflated BGZF block.
pub struct Block {
    /// file offset of the block (upper 48 bits of a virtual offset)
    pub coff: u64,
    pub data: Vec<u8>,
}

/// Parses the BGZF block headers in `buf` (read from file offset `base`). Returns
/// (offset in buf, total block size) of every complete block.
fn split_blocks(buf: &[u8]) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let mut p = 0;
    while p + 18 <= buf.len() {
        let h = &buf[p..];
        assert!(h[0] == 31 && h[1] == 139 && h[3] & 4 != 0, "not a BGZF block");
        let xlen = u16::from_le_bytes([h[10], h[11]]) as usize;
        let mut bsize = None;
        let mut x = 12;
        while x + 4 <= 12 + xlen {
            let slen = u16::from_le_bytes([h[x + 2], h[x + 3]]) as usize;
            if h[x] == 66 && h[x + 1] == 67 {
                bsize = Some(u16::from_le_bytes([h[x + 4], h[x + 5]]) as usize + 1);
            }
            x += 4 + slen;
        }
        let bsize = bsize.expect("BGZF block without BSIZE");
        if p + bsize > buf.len() {
            break;
        }
        out.push((p, bsize));
        p += bsize;
    }
    out
}

/// Inflates one raw BGZF block.
pub fn inflate_raw_block(raw: &[u8]) -> Vec<u8> {
    inflate_block(raw)
}

fn inflate_block(raw: &[u8]) -> Vec<u8> {
    let xlen = u16::from_le_bytes([raw[10], raw[11]]) as usize;
    let n = raw.len();
    let isize = u32::from_le_bytes(raw[n - 4..].try_into().unwrap()) as usize;
    let mut out = vec![0u8; isize];
    if isize > 0 {
        let cdata = &raw[12 + xlen..n - 8];
        DECOMP.with(|d| {
            let got = d.borrow_mut().deflate_decompress(cdata, &mut out).expect("BGZF inflate failed");
            assert_eq!(got, isize, "BGZF block size mismatch");
        });
    }
    out
}

/// Reads and inflates the blocks of `ranges` (each a file-offset range [first block, last block])
/// on a background thread, sending them in order, one `None` after each range.
pub fn spawn_block_reader(file: Arc<File>, ranges: Vec<(u64, u64)>, parallel: bool, inflate_ns: Arc<std::sync::atomic::AtomicU64>) -> Receiver<Option<Vec<Block>>> {
    let (tx, rx) = std::sync::mpsc::sync_channel::<Option<Vec<Block>>>(4);
    std::thread::Builder::new()
        .name("bgzf-read".into())
        .spawn(move || read_ranges(&file, &ranges, parallel, &tx, &inflate_ns))
        .unwrap();
    rx
}

fn read_ranges(file: &File, ranges: &[(u64, u64)], parallel: bool, tx: &SyncSender<Option<Vec<Block>>>, inflate_ns: &std::sync::atomic::AtomicU64) {
    const READ: usize = 4 << 20;
    let mut buf = vec![0u8; READ + 70000];
    let mut laps = crate::cpuprof::Laps::new();
    for &(first, last) in ranges {
        let mut pos = first;
        // blocks with offset <= last are needed
        while pos <= last {
            // only up to the end of the last needed block (a BGZF block is at most 64 KiB)
            let want = ((last - pos) as usize + 65_600).min(READ + 70000);
            let n = match file.read_at(&mut buf[..want], pos) {
                Ok(n) => n,
                Err(e) => panic!("BAM read error: {e}"),
            };
            laps.lap(0);
            // every range ends at a block the index points to: ending early means a damaged file
            assert!(n > 0, "BAM file truncated: end of file at offset {pos}, index expects data up to block {last}");
            let mut blocks = split_blocks(&buf[..n]);
            blocks.retain(|&(o, _)| pos + o as u64 <= last);
            // BGZF blocks are <= 64 KiB, so an empty list means an incomplete or corrupt block here
            assert!(!blocks.is_empty(), "BAM file truncated or corrupt at offset {pos} (index expects data up to block {last})");
            let t = Instant::now();
            let base = pos;
            let b = &buf;
            let out: Vec<Block> = if parallel && blocks.len() > 4 {
                blocks.par_iter().map(|&(o, s)| Block { coff: base + o as u64, data: inflate_block(&b[o..o + s]) }).collect()
            } else {
                blocks.iter().map(|&(o, s)| Block { coff: base + o as u64, data: inflate_block(&b[o..o + s]) }).collect()
            };
            inflate_ns.fetch_add(t.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
            laps.lap(1);
            let &(lo, ls) = blocks.last().unwrap();
            pos = base + (lo + ls) as u64;
            if tx.send(Some(out)).is_err() {
                return;
            }
        }
        if tx.send(None).is_err() {
            return;
        }
    }
}

/// Sequential reader over the decompressed bytes of a list of BAM index chunks
/// (virtual-offset ranges), positioned exactly like htsjdk's chunk iteration.
pub struct ChunkStream {
    rx: Receiver<Option<Vec<Block>>>,
    chunks: Vec<(u64, u64)>,
    ci: usize,
    blocks: std::collections::VecDeque<Block>,
    off: usize,
    range_done: bool,
    started: bool,
}

impl ChunkStream {
    /// `chunks` must be sorted and non-overlapping.
    pub fn new(file: Arc<File>, chunks: Vec<(u64, u64)>, parallel: bool, inflate_ns: Arc<std::sync::atomic::AtomicU64>) -> Self {
        let ranges: Vec<(u64, u64)> = chunks
            .iter()
            .map(|&(b, e)| {
                let last = if e & 0xffff == 0 { (e >> 16).saturating_sub(1) } else { e >> 16 };
                (b >> 16, last.max(b >> 16))
            })
            .collect();
        let rx = spawn_block_reader(file, ranges, parallel, inflate_ns);
        ChunkStream { rx, chunks, ci: 0, blocks: Default::default(), off: 0, range_done: false, started: false }
    }

    /// Next block of the current range; false at the end of the range.
    fn pull(&mut self) -> bool {
        if self.range_done {
            return false;
        }
        match self.rx.recv() {
            Ok(Some(v)) => {
                self.blocks.extend(v);
                true
            }
            Ok(None) => {
                self.range_done = true;
                false
            }
            Err(_) => panic!("BAM reader thread failed"),
        }
    }

    /// Moves to the next chunk: drains the rest of the current range.
    fn next_chunk(&mut self) -> bool {
        if self.started {
            while self.pull() {}
            self.blocks.clear();
            self.ci += 1;
        }
        self.started = true;
        if self.ci >= self.chunks.len() {
            return false;
        }
        self.range_done = false;
        if self.blocks.is_empty() && !self.pull() {
            return self.next_chunk();
        }
        self.off = (self.chunks[self.ci].0 & 0xffff) as usize;
        true
    }

    /// Ensures the cursor is on a block with unread data; false if the current range is exhausted.
    fn normalize(&mut self) -> bool {
        loop {
            if let Some(b) = self.blocks.front() {
                if self.off < b.data.len() {
                    return true;
                }
                self.blocks.pop_front();
                self.off = 0;
            }
            if self.blocks.is_empty() && !self.pull() {
                return false;
            }
        }
    }

    fn voff(&self) -> u64 {
        let b = self.blocks.front().unwrap();
        (b.coff << 16) | self.off as u64
    }

    fn read_into(&mut self, mut n: usize, out: &mut Vec<u8>) -> bool {
        while n > 0 {
            if !self.normalize() {
                return false;
            }
            let b = self.blocks.front().unwrap();
            let k = n.min(b.data.len() - self.off);
            out.extend_from_slice(&b.data[self.off..self.off + k]);
            self.off += k;
            n -= k;
        }
        true
    }

    /// Reads the next BAM record (the bytes after `block_size`) into `out`. False at the end.
    pub fn next_record(&mut self, out: &mut Vec<u8>) -> bool {
        loop {
            if !self.started {
                if !self.next_chunk() {
                    return false;
                }
            }
            if self.ci >= self.chunks.len() {
                return false;
            }
            if self.normalize() && self.voff() < self.chunks[self.ci].1 {
                out.clear();
                assert!(self.read_into(4, out), "truncated BAM record");
                let bs = u32::from_le_bytes(out[..4].try_into().unwrap()) as usize;
                out.clear();
                assert!(self.read_into(bs, out), "truncated BAM record");
                return true;
            }
            if !self.next_chunk() {
                return false;
            }
        }
    }
}

/// Reads the whole decompressed stream from the file start up to `max` bytes (for the BAM header).
pub fn read_prefix(file: &File, need: impl Fn(&[u8]) -> bool) -> Vec<u8> {
    let mut out = Vec::new();
    let mut pos = 0u64;
    let mut buf = vec![0u8; 1 << 20];
    loop {
        let n = file.read_at(&mut buf, pos).expect("read BAM");
        if n == 0 {
            return out;
        }
        let blocks = split_blocks(&buf[..n]);
        if blocks.is_empty() {
            return out;
        }
        for &(o, s) in &blocks {
            out.extend(inflate_block(&buf[o..o + s]));
            if need(&out) {
                return out;
            }
        }
        let &(lo, ls) = blocks.last().unwrap();
        pos += (lo + ls) as u64;
    }
}

// ------------------------------------------------------------------------------------------------
// Output
// ------------------------------------------------------------------------------------------------

thread_local! {
    static COMP: RefCell<Option<(i32, Compressor)>> = RefCell::new(None);
}

/// Compresses one BGZF block (data <= 64 KiB).
pub fn deflate_block(data: &[u8], level: i32) -> Vec<u8> {
    COMP.with(|c| {
        let mut c = c.borrow_mut();
        if c.as_ref().map_or(true, |(l, _)| *l != level) {
            *c = Some((level, Compressor::new(CompressionLvl::new(level).expect("compression level"))));
        }
        let comp = &mut c.as_mut().unwrap().1;
        let mut out = vec![0u8; 18 + comp.deflate_compress_bound(data.len()) + 8];
        let mut n = comp.deflate_compress(data, &mut out[18..]).expect("deflate");
        if 18 + n + 8 > 65536 {
            // incompressible: store (cannot happen for <= 0xff00 bytes of text, kept for safety)
            out = vec![0u8; 18 + data.len() + 5 + 8];
            out[18] = 1;
            out[19..21].copy_from_slice(&(data.len() as u16).to_le_bytes());
            out[21..23].copy_from_slice(&(!(data.len() as u16)).to_le_bytes());
            out[23..23 + data.len()].copy_from_slice(data);
            n = data.len() + 5;
        }
        out.truncate(18 + n + 8);
        let total = out.len();
        out[..18].copy_from_slice(&[0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, 0x42, 0x43, 2, 0, 0, 0]);
        out[16..18].copy_from_slice(&((total - 1) as u16).to_le_bytes());
        let mut crc = libdeflater::Crc::new();
        crc.update(data);
        out[total - 8..total - 4].copy_from_slice(&crc.sum().to_le_bytes());
        out[total - 4..].copy_from_slice(&(data.len() as u32).to_le_bytes());
        out
    })
}

/// Called for every complete line written to a [`BgzfWriter`]: the line, and the pseudo virtual
/// offsets (block index << 16 | offset) of its start and of its end.
pub trait LineSink {
    fn line(&mut self, line: &[u8], start: u64, end: u64);
    /// false: binary data, no line tracking
    fn wants_lines(&self) -> bool {
        true
    }
}

pub struct NoLines;
impl LineSink for NoLines {
    fn line(&mut self, _: &[u8], _: u64, _: u64) {}
    fn wants_lines(&self) -> bool {
        false
    }
}

/// BGZF writer. Blocks are compressed in parallel in groups. Every complete line is reported to the
/// [`LineSink`] with *pseudo* virtual offsets (block index instead of file offset); `block_sizes` maps
/// block indices to compressed sizes once the writer is finished.
pub struct BgzfWriter<W: Write, L: LineSink> {
    out: W,
    level: i32,
    pending: Vec<Vec<u8>>,
    cur: Vec<u8>,
    /// compressed size of every block written so far
    pub block_sizes: Vec<u32>,
    pub sink: L,
    line: Vec<u8>,
    line_start: u64,
    parallel: bool,
    pub bytes_in: u64,
    pub deflate_secs: f64,
}

const GROUP: usize = 64;

impl<W: Write, L: LineSink> BgzfWriter<W, L> {
    pub fn new(out: W, level: i32, sink: L, parallel: bool) -> Self {
        BgzfWriter { out, level, pending: Vec::new(), cur: Vec::with_capacity(BLOCK_DATA), block_sizes: Vec::new(), sink, line: Vec::new(), line_start: 0, parallel, bytes_in: 0, deflate_secs: 0.0 }
    }

    fn block_index(&self) -> u64 {
        (self.block_sizes.len() + self.pending.len()) as u64
    }

    /// pseudo virtual offset of the next byte
    fn pos(&self) -> u64 {
        if self.cur.len() == BLOCK_DATA {
            (self.block_index() + 1) << 16
        } else {
            (self.block_index() << 16) | self.cur.len() as u64
        }
    }

    fn flush_group(&mut self) -> io::Result<()> {
        if self.pending.is_empty() {
            return Ok(());
        }
        let t = Instant::now();
        let level = self.level;
        let comp: Vec<Vec<u8>> = if self.parallel && self.pending.len() > 1 {
            self.pending.par_iter().map(|b| deflate_block(b, level)).collect()
        } else {
            self.pending.iter().map(|b| deflate_block(b, level)).collect()
        };
        self.deflate_secs += t.elapsed().as_secs_f64();
        for c in &comp {
            self.out.write_all(c)?;
            self.block_sizes.push(c.len() as u32);
        }
        self.pending.clear();
        Ok(())
    }

    fn push_bytes(&mut self, mut data: &[u8]) -> io::Result<()> {
        while !data.is_empty() {
            if self.cur.len() == BLOCK_DATA {
                let full = std::mem::replace(&mut self.cur, Vec::with_capacity(BLOCK_DATA));
                self.pending.push(full);
                if self.pending.len() >= GROUP {
                    self.flush_group()?;
                }
            }
            let k = data.len().min(BLOCK_DATA - self.cur.len());
            self.cur.extend_from_slice(&data[..k]);
            data = &data[k..];
        }
        Ok(())
    }

    /// Ends the current block (the next byte starts a new block). Used between file parts.
    pub fn end_block(&mut self) -> io::Result<()> {
        if !self.cur.is_empty() {
            let full = std::mem::replace(&mut self.cur, Vec::with_capacity(BLOCK_DATA));
            self.pending.push(full);
        }
        Ok(())
    }

    /// Flushes everything (without the EOF marker) and returns the underlying writer.
    pub fn finish(mut self) -> io::Result<(W, Vec<u32>, L, f64)> {
        assert!(self.line.is_empty(), "unterminated line");
        self.end_block()?;
        self.flush_group()?;
        self.out.flush()?;
        Ok((self.out, self.block_sizes, self.sink, self.deflate_secs))
    }
}

impl<W: Write, L: LineSink> Write for BgzfWriter<W, L> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.bytes_in += buf.len() as u64;
        if !self.sink.wants_lines() {
            self.push_bytes(buf)?;
            return Ok(buf.len());
        }
        let mut rest = buf;
        while let Some(i) = memchr::memchr(b'\n', rest) {
            let (head, tail) = rest.split_at(i + 1);
            if self.line.is_empty() {
                self.line_start = self.pos();
            }
            self.line.extend_from_slice(head);
            self.push_bytes(head)?;
            let end = self.pos();
            let line = std::mem::take(&mut self.line);
            self.sink.line(&line[..line.len() - 1], self.line_start, end);
            self.line = line;
            self.line.clear();
            rest = tail;
        }
        if !rest.is_empty() {
            if self.line.is_empty() {
                self.line_start = self.pos();
            }
            self.line.extend_from_slice(rest);
            self.push_bytes(rest)?;
        }
        Ok(buf.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// Converts a pseudo virtual offset (block index << 16 | offset) to a real one, given the cumulative
/// compressed offsets of the blocks (`cum[i]` = file offset of block i; `cum.len()` = #blocks + 1).
pub fn real_voff(pseudo: u64, cum: &[u64]) -> u64 {
    let bi = (pseudo >> 16) as usize;
    (cum[bi] << 16) | (pseudo & 0xffff)
}

pub fn cumulative(base: u64, sizes: &[u32]) -> Vec<u64> {
    let mut cum = Vec::with_capacity(sizes.len() + 1);
    let mut o = base;
    cum.push(o);
    for &s in sizes {
        o += s as u64;
        cum.push(o);
    }
    cum
}

/// Virtual offset of uncompressed offset `u` of a BGZF file (block sizes from the BSIZE fields, uncompressed
/// sizes from ISIZE; nothing is inflated). An offset at a block boundary maps to the start of the next block.
pub fn voff_at(file: &File, u: u64) -> u64 {
    let len = file.metadata().map(|m| m.len()).unwrap_or(u64::MAX);
    let mut pos = 0u64;
    let mut cum = 0u64;
    let mut h = [0u8; 18];
    loop {
        assert!(pos < len, "BGZF: offset {u} is past the end of the data");
        file.read_exact_at(&mut h, pos).expect("BGZF block header");
        let bsize = split_blocks_header(&h, file, pos);
        let mut isz = [0u8; 4];
        file.read_exact_at(&mut isz, pos + bsize as u64 - 4).expect("BGZF ISIZE");
        let isize = u32::from_le_bytes(isz) as u64;
        // the end of the data (a BAM with a header and no records): the last block, the EOF block
        if cum + isize > u || (cum + isize == u && pos + bsize as u64 >= len) {
            return (pos << 16) | (u - cum);
        }
        cum += isize;
        pos += bsize as u64;
    }
}

fn split_blocks_header(h: &[u8; 18], file: &File, pos: u64) -> usize {
    assert!(h[0] == 31 && h[1] == 139 && h[3] & 4 != 0, "not a BGZF block");
    let xlen = u16::from_le_bytes([h[10], h[11]]) as usize;
    let mut extra = vec![0u8; xlen];
    file.read_exact_at(&mut extra, pos + 12).expect("BGZF extra field");
    let mut x = 0;
    while x + 4 <= xlen {
        let slen = u16::from_le_bytes([extra[x + 2], extra[x + 3]]) as usize;
        if extra[x] == 66 && extra[x + 1] == 67 {
            return u16::from_le_bytes([extra[x + 4], extra[x + 5]]) as usize + 1;
        }
        x += 4 + slen;
    }
    panic!("BGZF block without BSIZE");
}
