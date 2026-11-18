//! Known sites (`--known-sites` VCFs) for BaseRecalibrator: every record's span as htsjdk's VCF codec
//! gives it (`POS .. POS + len(REF) - 1`, or `INFO/END`), per contig, for overlap queries.
//! GATK passes all records overlapping the read's alignment span (`FeatureContext.getValues`),
//! whatever their FILTER; so do we.

use rayon::prelude::*;
use std::collections::HashMap;

pub struct KnownSites {
    /// per BAM reference index: sorted starts, matching ends, running max of ends
    contigs: Vec<Option<Contig>>,
    pub n_records: u64,
}

struct Contig {
    starts: Vec<i32>,
    ends: Vec<i32>,
    max_end: Vec<i32>,
}

/// BGZF block offsets of a whole file in memory.
fn block_offsets(data: &[u8]) -> Vec<usize> {
    let mut offs = Vec::new();
    let mut o = 0usize;
    while o + 18 <= data.len() {
        assert!(data[o] == 31 && data[o + 1] == 139, "not a BGZF file (bad block at offset {o})");
        let bsize = u16::from_le_bytes([data[o + 16], data[o + 17]]) as usize + 1;
        offs.push(o);
        o += bsize;
    }
    offs
}

fn inflate_block(data: &[u8], off: usize, d: &mut libdeflater::Decompressor, out: &mut Vec<u8>) {
    let bsize = u16::from_le_bytes([data[off + 16], data[off + 17]]) as usize + 1;
    let isize = u32::from_le_bytes(data[off + bsize - 4..off + bsize].try_into().unwrap()) as usize;
    let start = out.len();
    out.resize(start + isize, 0);
    if isize > 0 {
        let n = d.deflate_decompress(&data[off + 18..off + bsize - 8], &mut out[start..]).expect("BGZF inflate");
        assert_eq!(n, isize);
    }
}

fn parse_line(l: &[u8], out: &mut HashMap<Vec<u8>, Vec<(i32, i32)>>, n: &mut u64) {
    if l.is_empty() || l[0] == b'#' {
        return;
    }
    let mut f = l.splitn(9, |&b| b == b'\t');
    let chrom = f.next().unwrap();
    let pos: i32 = std::str::from_utf8(f.next().unwrap()).unwrap().parse().expect("VCF POS");
    let _id = f.next();
    let reflen = f.next().unwrap().len() as i32;
    let (_alt, _qual, _filter) = (f.next(), f.next(), f.next());
    let mut end = pos + reflen - 1;
    if let Some(info) = f.next() {
        let info = info.split(|&b| b == b'\t').next().unwrap();
        for kv in info.split(|&b| b == b';') {
            if let Some(v) = kv.strip_prefix(b"END=") {
                if let Ok(e) = std::str::from_utf8(v).unwrap_or("").parse::<i32>() {
                    end = e;
                }
            }
        }
    }
    *n += 1;
    match out.get_mut(chrom) {
        Some(v) => v.push((pos, end)),
        None => {
            out.insert(chrom.to_vec(), vec![(pos, end)]);
        }
    }
}

/// Parses one bgzipped VCF with `nthreads` workers. Each worker owns the lines that start in its
/// block range; it peeks at the previous block to know whether its first byte starts a line, and
/// reads past its range to finish its last line.
/// Plain text (an uncompressed VCF, as the hg38 bundle's dbsnp138.vcf that WARP passes, or one inflated from plain
/// gzip): `nthreads` workers, each parsing the lines that start in its byte range.
fn load_text(data: &[u8], nthreads: usize) -> (HashMap<Vec<u8>, Vec<(i32, i32)>>, u64) {
    let nw = nthreads.max(1);
    // range w starts after the first newline at or after w * len / nw (range 0 at 0)
    let cut = |w: usize| -> usize {
        if w == 0 {
            return 0;
        }
        if w >= nw {
            return data.len();
        }
        let at = data.len() * w / nw;
        memchr::memchr(b'\n', &data[at.saturating_sub(1)..]).map_or(data.len(), |p| at.saturating_sub(1) + p + 1)
    };
    let parts: Vec<(HashMap<Vec<u8>, Vec<(i32, i32)>>, u64)> = (0..nw)
        .into_par_iter()
        .map(|w| {
            let (a, b) = (cut(w), cut(w + 1));
            let mut out = HashMap::new();
            let mut n = 0u64;
            if a < b {
                for l in data[a..b].split(|&c| c == b'\n') {
                    parse_line(l.strip_suffix(b"\r").unwrap_or(l), &mut out, &mut n);
                }
            }
            (out, n)
        })
        .collect();
    merge_parts(parts)
}

fn merge_parts(parts: Vec<(HashMap<Vec<u8>, Vec<(i32, i32)>>, u64)>) -> (HashMap<Vec<u8>, Vec<(i32, i32)>>, u64) {
    let mut all: HashMap<Vec<u8>, Vec<(i32, i32)>> = HashMap::new();
    let mut n = 0;
    for (p, k) in parts {
        n += k;
        for (c, mut v) in p {
            all.entry(c).or_default().append(&mut v);
        }
    }
    (all, n)
}

fn load_vcf(path: &str, nthreads: usize) -> (HashMap<Vec<u8>, Vec<(i32, i32)>>, u64) {
    let data = std::fs::read(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let gzip = data.len() >= 2 && data[0] == 31 && data[1] == 139;
    // BGZF: gzip with the BC extra subfield
    let bgzf = gzip && data.len() >= 18 && data[3] & 4 != 0 && data[12] == b'B' && data[13] == b'C';
    if !bgzf {
        if gzip {
            let mut text = Vec::new();
            std::io::Read::read_to_end(&mut flate2::read::MultiGzDecoder::new(&data[..]), &mut text).unwrap_or_else(|e| panic!("{path}: {e}"));
            return load_text(&text, nthreads);
        }
        return load_text(&data, nthreads);
    }
    let offs = block_offsets(&data);
    let nb = offs.len();
    let nw = nthreads.max(1).min(nb.max(1));
    let parts: Vec<(HashMap<Vec<u8>, Vec<(i32, i32)>>, u64)> = (0..nw)
        .into_par_iter()
        .map(|w| {
            let (b0, b1) = (nb * w / nw, nb * (w + 1) / nw);
            let mut out = HashMap::new();
            let mut n = 0u64;
            if b0 >= b1 {
                return (out, n);
            }
            let mut dec = libdeflater::Decompressor::new();
            let mut buf = Vec::new();
            // does our first byte start a line?
            let mut skip_partial = false;
            if b0 > 0 {
                let mut prev = Vec::new();
                let mut k = b0;
                // the previous non-empty block decides
                while k > 0 {
                    k -= 1;
                    prev.clear();
                    inflate_block(&data, offs[k], &mut dec, &mut prev);
                    if !prev.is_empty() {
                        break;
                    }
                }
                skip_partial = prev.last().map_or(false, |&c| c != b'\n');
            }
            // our blocks, one at a time: parse the complete lines, carry the partial last one
            let mut skipping = skip_partial;
            let mut b = b0;
            while b < b1 {
                inflate_block(&data, offs[b], &mut dec, &mut buf);
                b += 1;
                if skipping {
                    // the line in progress started in an earlier range: not ours
                    match memchr::memchr(b'\n', &buf) {
                        Some(p) => {
                            buf.drain(..p + 1);
                            skipping = false;
                        }
                        None => {
                            buf.clear();
                            continue;
                        }
                    }
                }
                let mut consumed = 0;
                while let Some(p) = memchr::memchr(b'\n', &buf[consumed..]) {
                    parse_line(&buf[consumed..consumed + p], &mut out, &mut n);
                    consumed += p + 1;
                }
                buf.drain(..consumed);
            }
            // a line that started in our range: finish it from the following blocks (only that line)
            if !skipping && !buf.is_empty() {
                while b < nb {
                    let before = buf.len();
                    inflate_block(&data, offs[b], &mut dec, &mut buf);
                    b += 1;
                    if let Some(p) = memchr::memchr(b'\n', &buf[before..]) {
                        buf.truncate(before + p);
                        break;
                    }
                }
                parse_line(&buf, &mut out, &mut n);
            }
            (out, n)
        })
        .collect();
    merge_parts(parts)
}

impl KnownSites {
    /// Loads all `paths` and indexes them by the BAM's reference names.
    pub fn load(paths: &[String], ref_names: &[String], nthreads: usize) -> KnownSites {
        let mut merged: HashMap<Vec<u8>, Vec<(i32, i32)>> = HashMap::new();
        let mut n_records = 0;
        for p in paths {
            let (m, n) = load_vcf(p, nthreads);
            n_records += n;
            for (c, mut v) in m {
                merged.entry(c).or_default().append(&mut v);
            }
        }
        let contigs = ref_names
            .par_iter()
            .map(|name| {
                merged.get(name.as_bytes()).map(|v| {
                    let mut v = v.clone();
                    v.sort_unstable();
                    let starts: Vec<i32> = v.iter().map(|x| x.0).collect();
                    let ends: Vec<i32> = v.iter().map(|x| x.1).collect();
                    let mut max_end = Vec::with_capacity(ends.len());
                    let mut m = i32::MIN;
                    for &e in &ends {
                        m = m.max(e);
                        max_end.push(m);
                    }
                    Contig { starts, ends, max_end }
                })
            })
            .collect();
        KnownSites { contigs, n_records }
    }

    /// Calls `f(start, end)` for every site overlapping `[s, e]` (1-based, inclusive) on `tid`.
    #[inline]
    pub fn for_each_overlapping(&self, tid: usize, s: i64, e: i64, mut f: impl FnMut(i64, i64)) {
        let c = match self.contigs.get(tid) {
            Some(Some(c)) => c,
            _ => return,
        };
        let hi = c.starts.partition_point(|&x| (x as i64) <= e);
        let lo = c.max_end[..hi].partition_point(|&m| (m as i64) < s);
        for i in lo..hi {
            if c.ends[i] as i64 >= s {
                f(c.starts[i] as i64, c.ends[i] as i64);
            }
        }
    }
}
