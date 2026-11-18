//! `--germline-resource`: population allele frequencies (INFO/AF) for the POPAF annotation
//! (`SomaticGenotypingEngine.getGermlineAltAlleleFrequencies`). A bgzipped VCF with a tabix index (`.tbi`) is read
//! on demand through the index (16 kb windows, cached per thread), as GATK does; a plain VCF is read once for the
//! called contigs and kept in memory.

use flate2::read::MultiGzDecoder;
use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Read};
use std::os::unix::fs::FileExt;
use std::sync::Arc;

pub(crate) struct Rec {
    pub(crate) start: i64,
    pub(crate) alleles: Vec<Vec<u8>>,
    /// raw AF values (strings as in the file); None if absent
    af: Option<Vec<String>>,
    /// ID column
    pub(crate) id: String,
    /// FILTER is neither "." nor "PASS" (`VariantContext.isFiltered`)
    pub(crate) filtered: bool,
}

/// htsjdk SimpleAllele: bases are upper-cased unless the allele is symbolic (`<...>`, a breakend with `[` or `]`,
/// or a single breakend `.X` / `X.`); VCF allows lower case, and alleles are compared after this
fn allele_bases(a: &str) -> Vec<u8> {
    let b = a.as_bytes();
    let symbolic = b.len() > 1 && (b[0] == b'<' || b[b.len() - 1] == b'>' || b.iter().any(|&c| c == b'[' || c == b']') || b[0] == b'.' || b[b.len() - 1] == b'.');
    if symbolic {
        b.to_vec()
    } else {
        b.to_ascii_uppercase()
    }
}

/// One VCF data line -> (contig, record).
fn parse_line(line: &str) -> Option<(&str, Rec)> {
    if line.starts_with('#') || line.is_empty() {
        return None;
    }
    let mut it = line.split('\t');
    let chrom = it.next()?;
    let pos: i64 = it.next()?.parse().ok()?;
    let id = it.next()?.to_string();
    let r = allele_bases(it.next()?);
    let alt = it.next()?;
    let _qual = it.next();
    let filter = it.next().unwrap_or(".");
    let info = it.next().unwrap_or(".");
    let mut alleles = vec![r];
    if alt != "." {
        alleles.extend(alt.split(',').map(allele_bases));
    }
    let mut af = None;
    for kv in info.split(';') {
        if let Some(v) = kv.strip_prefix("AF=") {
            af = Some(v.split(',').map(|s| s.to_string()).collect());
        }
    }
    Some((chrom, Rec { start: pos, alleles, af, id, filtered: filter != "." && filter != "PASS" }))
}

/// Tabix-indexed bgzipped VCF.
struct Indexed {
    id: usize,
    file: std::fs::File,
    index: crate::bam::Bai,
    tid: HashMap<String, usize>,
    /// drop repeated identical records (the germline resource); off for dbSNP, where GATK sees every line
    dedup: bool,
}

const WINDOW: i64 = 1 << 14;

thread_local! {
    /// (resource id, contig, window start, records starting in [start, start + WINDOW))
    static CACHE: RefCell<Option<(usize, String, i64, Arc<Vec<Rec>>)>> = const { RefCell::new(None) };
    /// (resource id, contig, window start, records overlapping [start, start + WINDOW)) for overlap queries
    static OVL_CACHE: RefCell<Option<(usize, String, i64, Arc<Vec<Rec>>)>> = const { RefCell::new(None) };
}

impl Indexed {
    fn open(path: &str, tbi: &str, id: usize, dedup: bool) -> Indexed {
        let mut b = Vec::new();
        MultiGzDecoder::new(std::fs::File::open(tbi).unwrap_or_else(|e| panic!("{tbi}: {e}"))).read_to_end(&mut b).unwrap_or_else(|e| panic!("{tbi}: {e}"));
        assert_eq!(&b[..4], b"TBI\x01", "{tbi} is not a tabix index");
        let i32_at = |o: usize| i32::from_le_bytes(b[o..o + 4].try_into().unwrap());
        let n_ref = i32_at(4) as usize;
        let l_nm = i32_at(32) as usize;
        let names: Vec<String> = b[36..36 + l_nm].split(|&c| c == 0).filter(|n| !n.is_empty()).map(|n| String::from_utf8_lossy(n).into_owned()).collect();
        assert_eq!(names.len(), n_ref, "{tbi}: bad sequence name block");
        let mut o = 36 + l_nm;
        let mut refs = Vec::with_capacity(n_ref);
        for _ in 0..n_ref {
            let n_bin = i32_at(o) as usize;
            o += 4;
            let mut bins = HashMap::with_capacity(n_bin);
            for _ in 0..n_bin {
                let bin = i32_at(o) as u32;
                let n_chunk = i32_at(o + 4) as usize;
                o += 8;
                let mut ch = Vec::with_capacity(n_chunk);
                for _ in 0..n_chunk {
                    let s = u64::from_le_bytes(b[o..o + 8].try_into().unwrap());
                    let e = u64::from_le_bytes(b[o + 8..o + 16].try_into().unwrap());
                    ch.push((s, e));
                    o += 16;
                }
                if bin != 37450 {
                    bins.insert(bin, ch);
                }
            }
            let n_intv = i32_at(o) as usize;
            o += 4;
            let linear = (0..n_intv).map(|k| u64::from_le_bytes(b[o + 8 * k..o + 8 * k + 8].try_into().unwrap())).collect();
            o += 8 * n_intv;
            refs.push(crate::bam::RefIndex { bins, linear });
        }
        let tid = names.into_iter().enumerate().map(|(i, n)| (n, i)).collect();
        Indexed { id, file: std::fs::File::open(path).unwrap_or_else(|e| panic!("{path}: {e}")), index: crate::bam::Bai { refs }, tid, dedup }
    }

    /// Uncompressed bytes of a BGZF block at file offset `coff` and its compressed size.
    fn block(&self, coff: u64, dec: &mut libdeflater::Decompressor) -> (Vec<u8>, u64) {
        let mut h = [0u8; 18];
        self.file.read_exact_at(&mut h, coff).expect("germline resource: read BGZF header");
        let bsize = u16::from_le_bytes([h[16], h[17]]) as usize + 1;
        let mut blk = vec![0u8; bsize];
        self.file.read_exact_at(&mut blk, coff).expect("germline resource: read BGZF block");
        let isize = u32::from_le_bytes(blk[bsize - 4..].try_into().unwrap()) as usize;
        let mut out = vec![0u8; isize];
        if isize > 0 {
            dec.deflate_decompress(&blk[18..bsize - 8], &mut out).expect("germline resource: bad BGZF block");
        }
        (out, bsize as u64)
    }

    /// Records of `contig` starting in [beg, end), in file order.
    fn fetch(&self, contig: &str, beg: i64, end: i64) -> Vec<Rec> {
        self.fetch_with(contig, beg, end, false)
    }

    /// Records of `contig` starting in [beg, end) or, with `overlapping`, whose span (POS .. POS + |REF| - 1)
    /// overlaps [beg, end), in file order.
    fn fetch_with(&self, contig: &str, beg: i64, end: i64, overlapping: bool) -> Vec<Rec> {
        let Some(&t) = self.tid.get(contig) else { return Vec::new() };
        let chunks = self.index.query(t, &[(beg, end - 1)]);
        let mut dec = libdeflater::Decompressor::new();
        let mut text: Vec<u8> = Vec::new();
        for (vs, ve) in chunks {
            let (mut coff, uoff) = (vs >> 16, (vs & 0xffff) as usize);
            let (cend, uend) = (ve >> 16, (ve & 0xffff) as usize);
            let mut first = true;
            while coff <= cend {
                let (data, csize) = self.block(coff, &mut dec);
                let lo = if first { uoff } else { 0 };
                let hi = if coff == cend { uend.min(data.len()) } else { data.len() };
                if lo < hi {
                    text.extend_from_slice(&data[lo..hi]);
                }
                first = false;
                coff += csize;
                if data.is_empty() {
                    break;
                }
            }
            if !text.ends_with(b"\n") {
                text.push(b'\n');
            }
        }
        let mut out = Vec::new();
        for line in text.split(|&c| c == b'\n') {
            let line = std::str::from_utf8(line).expect("germline resource: not UTF-8");
            if let Some((c, r)) = parse_line(line) {
                let keep = if overlapping { r.start < end && r.start + r.alleles[0].len() as i64 > beg } else { r.start >= beg && r.start < end };
                if c == contig && keep {
                    out.push(r);
                }
            }
        }
        out.sort_by_key(|r| r.start);
        if self.dedup {
            out.dedup_by(|a, b| a.start == b.start && a.alleles == b.alleles && a.af == b.af);
        }
        out
    }
}

/// A bgzipped, tabix-indexed text file read line by line over regions (GetPileupSummaries' resource).
pub(crate) struct TabixLines(Indexed);

impl TabixLines {
    /// Some if `path` is bgzipped with a `.tbi` beside it.
    pub(crate) fn open(path: &str) -> Option<TabixLines> {
        let tbi = format!("{path}.tbi");
        ((path.ends_with(".gz") || path.ends_with(".bgz")) && std::path::Path::new(&tbi).exists())
            .then(|| TabixLines(Indexed::open(path, &tbi, usize::MAX, false)))
    }

    /// Calls `f` on every line of the index chunks covering `ivs` (1-based, inclusive) of `contig`, in file order,
    /// each line once. Lines outside the intervals can be included; every line overlapping them is.
    pub(crate) fn for_each_line(&self, contig: &str, ivs: &[(i64, i64)], mut f: impl FnMut(&str)) {
        let Some(&t) = self.0.tid.get(contig) else { return };
        let mut dec = libdeflater::Decompressor::new();
        let mut carry: Vec<u8> = Vec::new();
        // the last block read: with many small intervals (exome targets) consecutive chunks mostly start in the block
        // the previous one ended in, which would otherwise be read and inflated once per interval
        let mut last: Option<(u64, std::rc::Rc<Vec<u8>>, u64)> = None;
        let emit = |bytes: &[u8], f: &mut dyn FnMut(&str)| f(std::str::from_utf8(bytes).expect("VCF: not UTF-8"));
        for (vs, ve) in self.0.index.query(t, ivs) {
            let (mut coff, uoff) = (vs >> 16, (vs & 0xffff) as usize);
            let (cend, uend) = (ve >> 16, (ve & 0xffff) as usize);
            let mut first = true;
            while coff <= cend {
                let (data, csize) = match &last {
                    Some((c, d, z)) if *c == coff => (d.clone(), *z),
                    _ => {
                        let (d, z) = self.0.block(coff, &mut dec);
                        let d = std::rc::Rc::new(d);
                        last = Some((coff, d.clone(), z));
                        (d, z)
                    }
                };
                let lo = if first { uoff } else { 0 };
                let hi = if coff == cend { uend.min(data.len()) } else { data.len() };
                if lo < hi {
                    carry.extend_from_slice(&data[lo..hi]);
                    if let Some(k) = carry.iter().rposition(|&c| c == b'\n') {
                        for line in carry[..k].split(|&c| c == b'\n') {
                            emit(line, &mut f);
                        }
                        carry.drain(..=k);
                    }
                }
                first = false;
                coff += csize;
                if data.is_empty() {
                    break;
                }
            }
            if !carry.is_empty() {
                emit(&carry, &mut f);
                carry.clear();
            }
        }
    }
}

enum Store {
    Memory(HashMap<String, Vec<Rec>>),
    Indexed(Indexed),
}

pub struct Germline {
    store: Store,
}

static NEXT_ID: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

impl Germline {
    /// Opens the resource: through its tabix index if `path` is bgzipped with a `.tbi` next to it, else loaded into
    /// memory for `contigs`.
    pub fn load(path: &str, contigs: &HashSet<String>) -> Germline {
        Self::load_with(path, contigs, true, "germline resource")
    }

    /// `load` for any sites VCF; `dedup` as for the germline resource.
    pub(crate) fn load_with(path: &str, contigs: &HashSet<String>, dedup: bool, what: &str) -> Germline {
        let tbi = format!("{path}.tbi");
        if (path.ends_with(".gz") || path.ends_with(".bgz")) && std::path::Path::new(&tbi).exists() {
            let id = NEXT_ID.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            eprintln!("{what}: {path} (tabix index)");
            return Germline { store: Store::Indexed(Indexed::open(path, &tbi, id, dedup)) };
        }
        let f = std::fs::File::open(path).unwrap_or_else(|e| panic!("{path}: {e}"));
        let r: Box<dyn Read> = if path.ends_with(".gz") || path.ends_with(".bgz") { Box::new(MultiGzDecoder::new(f)) } else { Box::new(f) };
        let mut by_contig: HashMap<String, Vec<Rec>> = HashMap::new();
        for line in BufReader::with_capacity(1 << 20, r).lines() {
            let line = line.unwrap();
            if let Some((chrom, rec)) = parse_line(&line) {
                if contigs.contains(chrom) {
                    by_contig.entry(chrom.to_string()).or_default().push(rec);
                }
            }
        }
        for v in by_contig.values_mut() {
            v.sort_by_key(|r| r.start);
        }
        let n: usize = by_contig.values().map(|v| v.len()).sum();
        eprintln!("{what}: {n} records on the called contigs (in memory)");
        Germline { store: Store::Memory(by_contig) }
    }

    /// Calls `f` with the records starting at `loc` (`featureContext.getValues(germlineResource, loc)` keeps only
    /// features that start at that position), in file order.
    pub(crate) fn with_records_at<T>(&self, contig: &str, loc: i64, f: impl FnOnce(&[Rec]) -> T) -> T {
        let at = |v: &[Rec]| -> (usize, usize) { (v.partition_point(|r| r.start < loc), v.partition_point(|r| r.start <= loc)) };
        match &self.store {
            Store::Memory(m) => match m.get(contig) {
                Some(v) => {
                    let (lo, hi) = at(v);
                    f(&v[lo..hi])
                }
                None => f(&[]),
            },
            Store::Indexed(ix) => {
                let recs = CACHE.with(|c| {
                    let mut c = c.borrow_mut();
                    let hit = c.as_ref().filter(|(id, ctg, w, _)| *id == ix.id && ctg == contig && loc >= *w && loc < *w + WINDOW).map(|x| x.3.clone());
                    hit.unwrap_or_else(|| {
                        let recs = Arc::new(ix.fetch(contig, loc, loc + WINDOW));
                        *c = Some((ix.id, contig.to_string(), loc, recs.clone()));
                        recs
                    })
                });
                let (lo, hi) = at(&recs);
                f(&recs[lo..hi])
            }
        }
    }

    /// Calls `f` with the records whose span covers `pos` (`featureContext.getValues(germlineResource, locus)`), in
    /// file order.
    pub(crate) fn with_records_overlapping<T>(&self, contig: &str, pos: i64, f: impl FnOnce(&[&Rec]) -> T) -> T {
        let covers = |r: &&Rec| r.start <= pos && pos < r.start + r.alleles[0].len() as i64;
        match &self.store {
            Store::Memory(m) => {
                let v: Vec<&Rec> = m.get(contig).map_or(Vec::new(), |v| v.iter().filter(covers).collect());
                f(&v)
            }
            Store::Indexed(ix) => {
                let recs = OVL_CACHE.with(|c| {
                    let mut c = c.borrow_mut();
                    let hit = c.as_ref().filter(|(id, ctg, w, _)| *id == ix.id && ctg == contig && pos >= *w && pos < *w + WINDOW).map(|x| x.3.clone());
                    hit.unwrap_or_else(|| {
                        let recs = Arc::new(ix.fetch_with(contig, pos, pos + WINDOW, true));
                        *c = Some((ix.id, contig.to_string(), pos, recs.clone()));
                        recs
                    })
                });
                let v: Vec<&Rec> = recs.iter().filter(covers).collect();
                f(&v)
            }
        }
    }

    /// Tumor-only `Mutect2Engine.isActive`: whether a germline record covering `pos` with a population AF of at
    /// least `max_af` makes the locus inactive for the tumor's likeliest alt index (`PileupQualBuffer`: 0-3 a base,
    /// 4 another substitution, 5 an indel). As GATK, a substitution is compared with the germline *reference*
    /// allele's first base.
    pub(crate) fn blocks_tumor_only(&self, contig: &str, pos: i64, likeliest: usize, max_af: f64) -> bool {
        self.with_records_overlapping(contig, pos, |recs| {
            for g in recs {
                // getAttributeAsDoubleList(vc, AF, 0.0): '.' -> 0.0
                let afs: Vec<f64> = g.af.as_ref().map_or(Vec::new(), |v| v.iter().map(|s| if s == "." { 0.0 } else { s.parse::<f64>().unwrap() }).collect());
                let r = &g.alleles[0];
                for (k, &af) in afs.iter().enumerate() {
                    if af < max_af {
                        continue;
                    }
                    let alt = &g.alleles[k + 1];
                    if likeliest == 5 && alt.len() != r.len() {
                        return true;
                    }
                    if likeliest < 4 && alt.len() == r.len() && b"ACGT"[likeliest] == r[0].to_ascii_uppercase() {
                        return true;
                    }
                }
            }
            false
        })
    }

    /// `getGermlineAltAlleleFrequencies(allAlleles, germlineVCs, default)`: one AF per alt of `alleles`
    /// (reference first).
    pub fn alt_allele_frequencies(&self, contig: &str, loc: i64, alleles: &[Vec<u8>], default: f64) -> Vec<f64> {
        self.with_records_at(contig, loc, |recs| Self::frequencies(recs, alleles, default))
    }

    fn frequencies(recs: &[Rec], alleles: &[Vec<u8>], default: f64) -> Vec<f64> {
        let mut af = vec![default; alleles.len()];
        for g in recs {
            let vals: Vec<f64> = match &g.af {
                Some(v) => v.iter().map(|s| if s == "." { default } else { s.parse::<f64>().unwrap() }).collect(),
                None => Vec::new(), // missing AF: GATK warns; the size check below skips the record
            };
            if vals.len() != g.alleles.len() - 1 {
                continue;
            }
            let idx = allele_indices(alleles, &g.alleles);
            for k in 1..alleles.len() {
                if let Some(gi) = idx[k] {
                    if gi >= 1 {
                        af[k] = vals[gi - 1];
                    }
                }
            }
        }
        af[1..].to_vec()
    }
}

/// `GATKVariantContextUtils.alleleIndices(alleles1, alleles2)`: index in `a2` of each allele of `a1`
/// after extending both to a common reference allele.
fn allele_indices(a1: &[Vec<u8>], a2: &[Vec<u8>]) -> Vec<Option<usize>> {
    let (r1, r2) = (&a1[0], &a2[0]);
    let common: &Vec<u8> = if r1.len() < r2.len() {
        r2
    } else if r2.len() < r1.len() {
        r1
    } else {
        assert!(r1 == r2, "germline resource and call have different reference alleles");
        r1
    };
    let extend = |a: &Vec<u8>, r: &Vec<u8>| -> Option<Vec<u8>> {
        if a.first() == Some(&b'<') || a.as_slice() == b"*" || a.as_slice() == b"." {
            return None;
        }
        let mut b = a.clone();
        b.extend_from_slice(&common[r.len()..]);
        Some(b)
    };
    // commonToIndex2: the reference maps through `null` to index 0
    let mut common_to_2: HashMap<Option<Vec<u8>>, usize> = HashMap::new();
    common_to_2.insert(None, 0);
    for (n, a) in a2.iter().enumerate().skip(1) {
        if let Some(e) = extend(a, r2) {
            common_to_2.insert(Some(e), n);
        }
    }
    a1.iter()
        .enumerate()
        .map(|(n, a)| if n == 0 { Some(0) } else { extend(a, r1).and_then(|e| common_to_2.get(&Some(e)).copied()) })
        .collect()
}

/// `--panel-of-normals` (`-pon`): record spans per contig. Mutect2 uses it twice: isActive returns 0 at loci that
/// any PoN record overlaps, and a call gets the PON flag when a PoN record starts at its (merged, untrimmed) start.
pub struct Pon {
    /// merged spans per contig (overlap queries)
    merged: HashMap<String, Vec<(i64, i64)>>,
    starts: HashMap<String, HashSet<i64>>,
}

impl Pon {
    pub fn load(path: &str, contigs: &HashSet<String>) -> Pon {
        let f = std::fs::File::open(path).unwrap_or_else(|e| panic!("{path}: {e}"));
        let r: Box<dyn Read> = if path.ends_with(".gz") || path.ends_with(".bgz") { Box::new(MultiGzDecoder::new(f)) } else { Box::new(f) };
        let mut spans: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
        for line in BufReader::with_capacity(1 << 20, r).lines() {
            let line = line.unwrap();
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            let f: Vec<&str> = line.splitn(9, '\t').collect();
            if !contigs.contains(f[0]) {
                continue;
            }
            let pos: i64 = f[1].parse().unwrap();
            let end = f[7].split(';').find_map(|kv| kv.strip_prefix("END=")).and_then(|v| v.parse().ok()).unwrap_or(pos + f[3].len() as i64 - 1);
            spans.entry(f[0].to_string()).or_default().push((pos, end));
        }
        let mut merged = HashMap::new();
        let mut starts = HashMap::new();
        for (c, mut v) in spans {
            starts.insert(c.clone(), v.iter().map(|s| s.0).collect());
            v.sort();
            let mut m: Vec<(i64, i64)> = Vec::new();
            for (s, e) in v {
                match m.last_mut() {
                    Some(l) if s <= l.1 => l.1 = l.1.max(e),
                    _ => m.push((s, e)),
                }
            }
            merged.insert(c, m);
        }
        Pon { merged, starts }
    }

    pub fn overlaps(&self, contig: &str, pos: i64) -> bool {
        self.merged.get(contig).is_some_and(|m| {
            let k = m.partition_point(|iv| iv.1 < pos);
            k < m.len() && m[k].0 <= pos
        })
    }

    pub fn starts_at(&self, contig: &str, pos: i64) -> bool {
        self.starts.get(contig).is_some_and(|s| s.contains(&pos))
    }
}
