//! Tabix (.tbi) index for bgzipped VCF, built incrementally while the records are written
//! (htslib's `hts_idx_push` scheme: min_shift 14, 5 levels, binning + linear index).

use crate::bgzf::{self, BgzfWriter, LineSink};
use std::collections::BTreeMap;
use std::io::{self, Write};

const MIN_SHIFT: u32 = 14;
const META_BIN: u32 = 37450;

/// `hts_reg2bin(beg, end, 14, 5)`, 0-based half-open.
fn reg2bin(beg: i64, end: i64) -> u32 {
    let end = end - 1;
    let mut s = MIN_SHIFT;
    let mut t: i64 = ((1i64 << 15) - 1) / 7;
    let mut l = 5;
    while l > 0 {
        if beg >> s == end >> s {
            return (t + (beg >> s)) as u32;
        }
        s += 3;
        t -= 1i64 << ((l - 1) * 3);
        l -= 1;
    }
    0
}

/// Index of one contig. Offsets are whatever the caller pushes (pseudo or real virtual offsets);
/// `map_offsets` converts them.
#[derive(Default, Clone)]
pub struct ContigIndex {
    pub name: String,
    bins: BTreeMap<u32, Vec<(u64, u64)>>,
    linear: Vec<u64>,
    n_rec: u64,
    off_beg: u64,
    off_end: u64,
    save_bin: Option<u32>,
    save_off: u64,
    last_off: u64,
}

impl ContigIndex {
    pub fn new(name: &str) -> Self {
        ContigIndex { name: name.to_string(), ..Default::default() }
    }

    /// One record covering [beg, end) (0-based) stored at virtual offsets [start, stop).
    pub fn push(&mut self, beg: i64, end: i64, start: u64, stop: u64) {
        let beg = beg.max(0);
        let end = end.max(beg + 1);
        if self.n_rec == 0 {
            self.off_beg = start;
        }
        let (w0, w1) = ((beg >> MIN_SHIFT) as usize, ((end - 1) >> MIN_SHIFT) as usize);
        if self.linear.len() <= w1 {
            self.linear.resize(w1 + 1, u64::MAX);
        }
        for w in w0..=w1 {
            if self.linear[w] == u64::MAX {
                self.linear[w] = start;
            }
        }
        let bin = reg2bin(beg, end);
        if self.save_bin != Some(bin) {
            self.close_chunk();
            self.save_bin = Some(bin);
            self.save_off = start;
        }
        self.last_off = stop;
        self.off_end = stop;
        self.n_rec += 1;
    }

    fn close_chunk(&mut self) {
        if let Some(b) = self.save_bin {
            let v = self.bins.entry(b).or_default();
            match v.last_mut() {
                Some(l) if l.1 >= self.save_off => l.1 = l.1.max(self.last_off),
                _ => v.push((self.save_off, self.last_off)),
            }
        }
        self.save_bin = None;
    }

    pub fn finish(&mut self) {
        self.close_chunk();
        // empty windows take the previous window's offset (the first set one for leading windows)
        let first = self.linear.iter().copied().find(|&x| x != u64::MAX).unwrap_or(0);
        let mut prev = first;
        for x in self.linear.iter_mut() {
            if *x == u64::MAX {
                *x = prev;
            } else {
                prev = *x;
            }
        }
    }

    pub fn is_empty(&self) -> bool {
        self.n_rec == 0
    }

    /// Appends the index of a later piece of the same contig (a contig split into several parts).
    /// Both must be finished and use real offsets, with all of `other` stored after `self`.
    pub fn merge(&mut self, other: ContigIndex) {
        debug_assert_eq!(self.name, other.name);
        if other.n_rec == 0 {
            return;
        }
        if self.n_rec == 0 {
            *self = other;
            return;
        }
        for (bin, chunks) in other.bins {
            let v = self.bins.entry(bin).or_default();
            for ch in chunks {
                match v.last_mut() {
                    Some(l) if l.1 >= ch.0 => l.1 = l.1.max(ch.1),
                    _ => v.push(ch),
                }
            }
        }
        // windows this index already covers keep their (smaller) offsets
        if other.linear.len() > self.linear.len() {
            let n = self.linear.len();
            self.linear.extend_from_slice(&other.linear[n..]);
        }
        self.n_rec += other.n_rec;
        self.off_end = other.off_end;
    }

    pub fn map_offsets(&mut self, f: impl Fn(u64) -> u64) {
        for v in self.bins.values_mut() {
            for c in v.iter_mut() {
                *c = (f(c.0), f(c.1));
            }
        }
        for x in self.linear.iter_mut() {
            *x = f(*x);
        }
        self.off_beg = f(self.off_beg);
        self.off_end = f(self.off_end);
    }
}

/// Builds the index of the VCF lines written through a [`BgzfWriter`] (header lines are skipped).
#[derive(Default)]
pub struct VcfIndexer {
    pub contigs: Vec<ContigIndex>,
    pub n_lines: u64,
}

impl LineSink for VcfIndexer {
    fn line(&mut self, line: &[u8], voff_start: u64, voff_end: u64) {
        if line.first() == Some(&b'#') {
            return;
        }
        self.n_lines += 1;
        let mut it = line.splitn(9, |&b| b == b'\t');
        let chrom = it.next().unwrap();
        let pos: i64 = std::str::from_utf8(it.next().unwrap()).unwrap().parse().unwrap();
        let _id = it.next();
        let reflen = it.next().unwrap().len() as i64;
        let (_alt, _qual, _filter) = (it.next(), it.next(), it.next());
        let info = it.next().unwrap_or(b".");
        let beg = pos - 1;
        let mut end = beg + reflen;
        // htslib uses INFO/END when present
        for kv in info.split(|&b| b == b';') {
            if let Some(v) = kv.strip_prefix(b"END=") {
                if let Ok(e) = std::str::from_utf8(v).unwrap_or("").parse::<i64>() {
                    if e > beg {
                        end = e;
                    }
                }
                break;
            }
        }
        if self.contigs.last().map_or(true, |c| c.name.as_bytes() != chrom) {
            if let Some(c) = self.contigs.last_mut() {
                c.finish();
            }
            self.contigs.push(ContigIndex::new(std::str::from_utf8(chrom).unwrap()));
        }
        self.contigs.last_mut().unwrap().push(beg, end, voff_start, voff_end);
    }
}

impl VcfIndexer {
    pub fn finish(&mut self) {
        if let Some(c) = self.contigs.last_mut() {
            c.finish();
        }
    }
}

/// Writes the .tbi for `contigs` (offsets already real).
pub fn write_tbi(path: &str, contigs: &[ContigIndex]) -> io::Result<()> {
    let contigs: Vec<&ContigIndex> = contigs.iter().filter(|c| !c.is_empty()).collect();
    let mut b: Vec<u8> = Vec::new();
    b.extend_from_slice(b"TBI\x01");
    let i32le = |b: &mut Vec<u8>, v: i32| b.extend_from_slice(&v.to_le_bytes());
    i32le(&mut b, contigs.len() as i32);
    for v in [2, 1, 2, 0, b'#' as i32, 0] {
        i32le(&mut b, v);
    }
    let names: Vec<u8> = contigs.iter().flat_map(|c| c.name.bytes().chain(std::iter::once(0))).collect();
    i32le(&mut b, names.len() as i32);
    b.extend_from_slice(&names);
    for c in &contigs {
        i32le(&mut b, (c.bins.len() + 1) as i32);
        for (&bin, ch) in &c.bins {
            b.extend_from_slice(&bin.to_le_bytes());
            i32le(&mut b, ch.len() as i32);
            for &(s, e) in ch {
                b.extend_from_slice(&s.to_le_bytes());
                b.extend_from_slice(&e.to_le_bytes());
            }
        }
        b.extend_from_slice(&META_BIN.to_le_bytes());
        i32le(&mut b, 2);
        for v in [c.off_beg, c.off_end, c.n_rec, 0] {
            b.extend_from_slice(&v.to_le_bytes());
        }
        i32le(&mut b, c.linear.len() as i32);
        for &x in &c.linear {
            b.extend_from_slice(&x.to_le_bytes());
        }
    }
    b.extend_from_slice(&0u64.to_le_bytes());
    let f = std::fs::File::create(path)?;
    let mut w = BgzfWriter::new(std::io::BufWriter::new(f), 6, bgzf::NoLines, false);
    w.write_all(&b)?;
    let (mut out, _, _, _) = w.finish()?;
    out.write_all(&bgzf::EOF_BLOCK)?;
    out.flush()
}

/// Writes a BAI (`.bai`) for `n_ref` references. `contigs` holds (tid, index, n_mapped, n_unmapped)
/// for the references with records (offsets already real); `n_no_coor` = unplaced unmapped reads.
pub fn write_bai(path: &str, n_ref: usize, contigs: &[(usize, ContigIndex, u64, u64)], n_no_coor: u64) -> io::Result<()> {
    std::fs::write(path, bai_bytes(n_ref, contigs, n_no_coor))
}

/// The `.bai` file's bytes.
pub fn bai_bytes(n_ref: usize, contigs: &[(usize, ContigIndex, u64, u64)], n_no_coor: u64) -> Vec<u8> {
    let mut b: Vec<u8> = Vec::new();
    b.extend_from_slice(b"BAI\x01");
    b.extend_from_slice(&(n_ref as i32).to_le_bytes());
    let mut by_tid: Vec<Option<&(usize, ContigIndex, u64, u64)>> = vec![None; n_ref];
    for c in contigs {
        by_tid[c.0] = Some(c);
    }
    for slot in by_tid {
        match slot {
            Some((_, c, n_mapped, n_unmapped)) if !c.is_empty() => {
                b.extend_from_slice(&((c.bins.len() + 1) as i32).to_le_bytes());
                for (&bin, ch) in &c.bins {
                    b.extend_from_slice(&bin.to_le_bytes());
                    b.extend_from_slice(&(ch.len() as i32).to_le_bytes());
                    for &(s, e) in ch {
                        b.extend_from_slice(&s.to_le_bytes());
                        b.extend_from_slice(&e.to_le_bytes());
                    }
                }
                b.extend_from_slice(&META_BIN.to_le_bytes());
                b.extend_from_slice(&2i32.to_le_bytes());
                for v in [c.off_beg, c.off_end, *n_mapped, *n_unmapped] {
                    b.extend_from_slice(&v.to_le_bytes());
                }
                b.extend_from_slice(&(c.linear.len() as i32).to_le_bytes());
                for &x in &c.linear {
                    b.extend_from_slice(&x.to_le_bytes());
                }
            }
            _ => {
                b.extend_from_slice(&0i32.to_le_bytes());
                b.extend_from_slice(&0i32.to_le_bytes());
            }
        }
    }
    b.extend_from_slice(&n_no_coor.to_le_bytes());
    b
}
