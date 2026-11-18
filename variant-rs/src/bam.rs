//! BAM header, record decoding and the BAI index (chunk lists for interval queries).

use crate::bgzf;
use std::fs::File;
use std::io::Read;

pub struct Header {
    pub text: String,
    pub ref_names: Vec<String>,
    pub ref_lens: Vec<i64>,
}

impl Header {
    pub fn read(file: &File) -> Header {
        // header = magic, l_text, text, n_ref, refs; read enough decompressed bytes to cover it
        let need = |b: &[u8]| -> bool {
            if b.len() < 12 {
                return false;
            }
            let l_text = i32::from_le_bytes(b[4..8].try_into().unwrap()) as usize;
            let mut o = 8 + l_text;
            if b.len() < o + 4 {
                return false;
            }
            let n_ref = i32::from_le_bytes(b[o..o + 4].try_into().unwrap()) as usize;
            o += 4;
            for _ in 0..n_ref {
                if b.len() < o + 4 {
                    return false;
                }
                let l = i32::from_le_bytes(b[o..o + 4].try_into().unwrap()) as usize;
                o += 4 + l + 4;
                if b.len() < o {
                    return false;
                }
            }
            true
        };
        let b = bgzf::read_prefix(file, need);
        assert_eq!(&b[..4], b"BAM\x01", "not a BAM file");
        let l_text = i32::from_le_bytes(b[4..8].try_into().unwrap()) as usize;
        let text = String::from_utf8_lossy(&b[8..8 + l_text]).trim_end_matches('\0').to_string();
        let mut o = 8 + l_text;
        let n_ref = i32::from_le_bytes(b[o..o + 4].try_into().unwrap()) as usize;
        o += 4;
        let (mut ref_names, mut ref_lens) = (Vec::with_capacity(n_ref), Vec::with_capacity(n_ref));
        for _ in 0..n_ref {
            let l = i32::from_le_bytes(b[o..o + 4].try_into().unwrap()) as usize;
            o += 4;
            ref_names.push(String::from_utf8_lossy(&b[o..o + l - 1]).to_string());
            o += l;
            ref_lens.push(i32::from_le_bytes(b[o..o + 4].try_into().unwrap()) as i64);
            o += 4;
        }
        Header { text, ref_names, ref_lens }
    }

    /// `SM` of every `@RG` line.
    /// (ID, SM) of every `@RG` line that has both
    pub fn read_group_samples(&self) -> Vec<(String, String)> {
        self.text
            .lines()
            .filter(|l| l.starts_with("@RG"))
            .filter_map(|l| {
                let id = l.split('\t').find_map(|f| f.strip_prefix("ID:"))?;
                let sm = l.split('\t').find_map(|f| f.strip_prefix("SM:"))?;
                Some((id.to_string(), sm.to_string()))
            })
            .collect()
    }

    pub fn sample_names(&self) -> Vec<String> {
        self.text.lines().filter(|l| l.starts_with("@RG")).filter_map(|l| l.split('\t').find_map(|f| f.strip_prefix("SM:")).map(|v| v.to_string())).collect()
    }

    /// `AS` tag of every `@SQ` line, by sequence name (htsjdk puts it in the VCF contig line).
    pub fn sq_assembly(&self) -> std::collections::HashMap<String, String> {
        let mut m = std::collections::HashMap::new();
        for l in self.text.lines().filter(|l| l.starts_with("@SQ")) {
            let (mut sn, mut asm) = (None, None);
            for f in l.split('\t') {
                if let Some(v) = f.strip_prefix("SN:") {
                    sn = Some(v.to_string());
                }
                if let Some(v) = f.strip_prefix("AS:") {
                    asm = Some(v.to_string());
                }
            }
            if let (Some(s), Some(a)) = (sn, asm) {
                m.insert(s, a);
            }
        }
        m
    }
}

/// One BAM record, decoded just enough for HaplotypeCaller.
pub struct Rec {
    pub name: Vec<u8>,
    pub ref_id: i32,
    pub start: i64, // 1-based
    pub mapq: u8,
    pub flag: u16,
    pub cigar: Vec<(u8, u32)>,
    pub seq: Vec<u8>, // ASCII bases
    pub qual: Vec<u8>,
    pub has_rg: bool,
    pub qual_missing: bool,
    pub next_ref_id: i32,
    pub mate_start: i64, // 1-based, 0 if unset
    pub tlen: i32,
}

impl Default for Rec {
    fn default() -> Self {
        Rec { name: vec![], ref_id: -1, start: 0, mapq: 0, flag: 0, cigar: vec![], seq: vec![], qual: vec![], has_rg: false, qual_missing: false, next_ref_id: -1, mate_start: 0, tlen: 0 }
    }
}

const SEQ_NT16: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// Fixed-size part of a record: (ref_id, 1-based start, flag, reference span end) without decoding the rest.
#[inline]
pub fn peek_core(b: &[u8]) -> (i32, i64, u16) {
    let ref_id = i32::from_le_bytes(b[0..4].try_into().unwrap());
    let start = i32::from_le_bytes(b[4..8].try_into().unwrap()) as i64 + 1;
    let flag = u16::from_le_bytes([b[14], b[15]]);
    (ref_id, start, flag)
}

/// The bin field htsjdk's BAMRecordCodec writes: 0 for an unplaced read (POS -1) and on a reference too long for
/// the BAI (`too_large`, > 2^29 bp); otherwise SAMRecord.computeIndexingBin, regionToBin(start - 1, alignment end)
/// with the end of an unmapped read (0) taken as start. htslib differs for unplaced reads (4680) and for a CIGAR
/// without reference bases at a 16 kb boundary, so records copied from an htslib file need it recomputed.
pub fn htsjdk_bin(rec: &[u8], too_large: bool) -> u16 {
    let (_, start, flag) = peek_core(rec);
    if start == 0 || too_large {
        return 0;
    }
    let beg = start - 1;
    let mut end = if flag & 0x4 != 0 { 0 } else { query_end(rec) };
    if end <= 0 {
        end = beg + 1;
    }
    // GenomicIndexUtil.regionToBin(beg, end), end exclusive
    let e = end - 1;
    for (sh, off) in [(14, 4681i64), (17, 585), (20, 73), (23, 9), (26, 1)] {
        if beg >> sh == e >> sh {
            return ((off + (beg >> sh)) & 0xffff) as u16;
        }
    }
    0
}

/// htsjdk `getAlignmentEnd()` for the query filter: start + reference length - 1 (start - 1 for a
/// cigar without reference bases); for an unmapped read placed at its mate: its start.
pub fn query_end(b: &[u8]) -> i64 {
    let (_, start, flag) = peek_core(b);
    if flag & 0x4 != 0 {
        return start;
    }
    let l_read_name = b[8] as usize;
    let n_cigar = u16::from_le_bytes([b[12], b[13]]) as usize;
    let mut o = 32 + l_read_name;
    let mut rl = 0i64;
    for _ in 0..n_cigar {
        let v = u32::from_le_bytes(b[o..o + 4].try_into().unwrap());
        if matches!(v & 0xf, 0 | 2 | 3 | 7 | 8) {
            rl += (v >> 4) as i64;
        }
        o += 4;
    }
    start + rl - 1
}

/// Parses one record from `b` (the bytes after block_size). Reuses `rec` buffers.
pub fn parse_record(b: &[u8], rec: &mut Rec) {
    rec.ref_id = i32::from_le_bytes(b[0..4].try_into().unwrap());
    rec.start = i32::from_le_bytes(b[4..8].try_into().unwrap()) as i64 + 1;
    let l_read_name = b[8] as usize;
    rec.mapq = b[9];
    let n_cigar = u16::from_le_bytes([b[12], b[13]]) as usize;
    rec.flag = u16::from_le_bytes([b[14], b[15]]);
    let l_seq = i32::from_le_bytes(b[16..20].try_into().unwrap()) as usize;
    rec.next_ref_id = i32::from_le_bytes(b[20..24].try_into().unwrap());
    rec.mate_start = i32::from_le_bytes(b[24..28].try_into().unwrap()) as i64 + 1;
    rec.tlen = i32::from_le_bytes(b[28..32].try_into().unwrap());
    rec.name.clear();
    rec.name.extend_from_slice(&b[32..32 + l_read_name - 1]);
    let mut o = 32 + l_read_name;
    rec.cigar.clear();
    for _ in 0..n_cigar {
        let v = u32::from_le_bytes(b[o..o + 4].try_into().unwrap());
        rec.cigar.push(((v & 0xf) as u8, v >> 4));
        o += 4;
    }
    rec.seq.clear();
    rec.seq.reserve(l_seq);
    for i in 0..l_seq / 2 {
        let byte = b[o + i];
        rec.seq.push(SEQ_NT16[(byte >> 4) as usize]);
        rec.seq.push(SEQ_NT16[(byte & 0xf) as usize]);
    }
    if l_seq % 2 == 1 {
        rec.seq.push(SEQ_NT16[(b[o + l_seq / 2] >> 4) as usize]);
    }
    o += (l_seq + 1) / 2;
    rec.qual.clear();
    rec.qual.extend_from_slice(&b[o..o + l_seq]);
    rec.qual_missing = l_seq > 0 && b[o] == 0xff;
    o += l_seq;
    // aux tags: only need to know whether RG is present
    rec.has_rg = false;
    while o + 3 <= b.len() {
        let tag = [b[o], b[o + 1]];
        let ty = b[o + 2];
        o += 3;
        if tag == *b"RG" {
            rec.has_rg = true;
        }
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
                while b[o] != 0 {
                    o += 1;
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
}

// ------------------------------------------------------------------------------------------------
// BAI
// ------------------------------------------------------------------------------------------------

pub struct RefIndex {
    /// bin -> chunks (virtual offsets)
    pub bins: std::collections::HashMap<u32, Vec<(u64, u64)>>,
    pub linear: Vec<u64>,
}

pub struct Bai {
    pub refs: Vec<RefIndex>,
}

const META_BIN: u32 = 37450;

impl Bai {
    pub fn load(path: &str) -> std::io::Result<Bai> {
        let mut b = Vec::new();
        File::open(path)?.read_to_end(&mut b)?;
        assert_eq!(&b[..4], b"BAI\x01", "{path} is not a BAI index");
        Ok(Bai::from_bytes(&b))
    }

    /// A `.bai` held in memory.
    pub fn from_bytes(b: &[u8]) -> Bai {
        assert_eq!(&b[..4], b"BAI\x01", "not a BAI index");
        let rd_i32 = |o: &mut usize| -> i32 {
            let v = i32::from_le_bytes(b[*o..*o + 4].try_into().unwrap());
            *o += 4;
            v
        };
        let rd_u64 = |o: &mut usize| -> u64 {
            let v = u64::from_le_bytes(b[*o..*o + 8].try_into().unwrap());
            *o += 8;
            v
        };
        let mut o = 4;
        let n_ref = rd_i32(&mut o) as usize;
        let mut refs = Vec::with_capacity(n_ref);
        for _ in 0..n_ref {
            let n_bin = rd_i32(&mut o) as usize;
            let mut bins = std::collections::HashMap::with_capacity(n_bin);
            for _ in 0..n_bin {
                let bin = rd_i32(&mut o) as u32;
                let n_chunk = rd_i32(&mut o) as usize;
                let mut ch = Vec::with_capacity(n_chunk);
                for _ in 0..n_chunk {
                    let s = rd_u64(&mut o);
                    let e = rd_u64(&mut o);
                    ch.push((s, e));
                }
                if bin != META_BIN {
                    bins.insert(bin, ch);
                }
            }
            let n_intv = rd_i32(&mut o) as usize;
            let mut linear = Vec::with_capacity(n_intv);
            for _ in 0..n_intv {
                linear.push(rd_u64(&mut o));
            }
            refs.push(RefIndex { bins, linear });
        }
        Bai { refs }
    }

    /// Chunks that may contain records of `tid` overlapping any of `intervals` (1-based inclusive),
    /// sorted and merged.
    pub fn query(&self, tid: usize, intervals: &[(i64, i64)]) -> Vec<(u64, u64)> {
        let Some(ri) = self.refs.get(tid) else { return vec![] };
        let mut chunks: Vec<(u64, u64)> = Vec::new();
        for &(s, e) in intervals {
            let (beg, end) = ((s - 1).max(0) as u64, e.max(s) as u64); // 0-based half-open
            let min_off = if ri.linear.is_empty() { 0 } else { ri.linear[((beg >> 14) as usize).min(ri.linear.len() - 1)] };
            for bin in reg2bins(beg, end) {
                if let Some(ch) = ri.bins.get(&bin) {
                    for &c in ch {
                        if c.1 > min_off {
                            chunks.push(c);
                        }
                    }
                }
            }
        }
        chunks.sort();
        let mut merged: Vec<(u64, u64)> = Vec::with_capacity(chunks.len());
        for c in chunks {
            match merged.last_mut() {
                Some(l) if c.0 <= l.1 => l.1 = l.1.max(c.1),
                _ => merged.push(c),
            }
        }
        merged
    }
}

/// Bins overlapping [beg, end) (0-based), as in the SAM specification.
fn reg2bins(beg: u64, end: u64) -> Vec<u32> {
    let end = end - 1;
    let mut v = vec![0u32];
    for (off, sh) in [(1u64, 26), (9, 23), (73, 20), (585, 17), (4681, 14)] {
        for k in (off + (beg >> sh))..=(off + (end >> sh)) {
            v.push(k as u32);
        }
    }
    v
}
