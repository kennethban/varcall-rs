//! Writing BAM files the way GATK/Picard (htsjdk) write modified records: coordinate-ordered chunks
//! processed in parallel, BGZF-compressed per chunk, concatenated in order, with a `.bai` built from
//! the records' virtual offsets. Also htsjdk's tag re-encoding and the shared chunk plan.

use crate::bgzf::{self, BLOCK_DATA};
use crate::input::{Input, Stream};
use crate::tabix::{self, ContigIndex};
use std::io::Write;
use std::sync::atomic::AtomicU64;
use std::sync::Arc;

/// A unit of output: records of contig `tid` whose start lies in `[from, to]` and that overlap one of `ivs` (the
/// traversal intervals the chunk owns: one long interval cut into pieces, or a run of short ones such as exome
/// targets, so that they share one index query); `tid == -1` is the unplaced unmapped tail.
#[derive(Clone, Debug)]
pub struct Chunk {
    pub tid: i32,
    pub from: i64,
    pub to: i64,
    pub ivs: Vec<(i64, i64)>,
    /// whole-file traversal, first chunk of its contig: also the contig's records at POS 0 (unmapped reads with an
    /// RNAME, which Picard / GATK copy through; samtools sort puts them first on the contig)
    pub pos0: bool,
}

impl Chunk {
    /// The regions to query: the intervals clipped to `[from, to]`.
    pub fn query(&self) -> Vec<(i64, i64)> {
        self.ivs.iter().map(|&(s, e)| (s.max(self.from), e.min(self.to))).filter(|&(s, e)| s <= e).collect()
    }

    /// Does `[start, end]` overlap one of the chunk's intervals?
    #[inline]
    pub fn overlaps(&self, start: i64, end: i64) -> bool {
        let k = self.ivs.partition_point(|iv| iv.1 < start);
        k < self.ivs.len() && self.ivs[k].0 <= end
    }
}

/// Chunks covering `traversal` (tid, start, end; sorted by tid then start, merged) in file order. An interval longer
/// than `chunk` is cut into pieces; consecutive shorter ones on a contig share a chunk while it spans at most
/// `chunk`. The first chunk of a run also takes reads starting before it (after the previous run on the contig,
/// which owns reads overlapping it). With `tail`, a final chunk for the unplaced unmapped reads.
pub fn plan_chunks(traversal: &[(usize, i64, i64)], n_ref: usize, chunk: i64, tail: bool) -> Vec<Chunk> {
    let mut chunks = Vec::new();
    let mut prev_end = vec![0i64; n_ref];
    let mut i = 0;
    while i < traversal.len() {
        let (t, s, e) = traversal[i];
        if e - s + 1 > chunk {
            let mut a = s;
            let mut first = true;
            while a <= e {
                let b = (a + chunk - 1).min(e);
                chunks.push(Chunk { tid: t as i32, from: if first { prev_end[t] + 1 } else { a }, to: b, ivs: vec![(s, e)], pos0: tail && first && prev_end[t] == 0 });
                first = false;
                a = b + 1;
            }
            prev_end[t] = e;
            i += 1;
            continue;
        }
        let mut ivs = Vec::new();
        while i < traversal.len() && traversal[i].0 == t && traversal[i].2 - s + 1 <= chunk {
            ivs.push((traversal[i].1, traversal[i].2));
            i += 1;
        }
        let to = ivs.last().unwrap().1;
        chunks.push(Chunk { tid: t as i32, from: prev_end[t] + 1, to, ivs, pos0: tail && prev_end[t] == 0 });
        prev_end[t] = to;
    }
    if tail {
        chunks.push(Chunk { tid: -1, from: 0, to: 0, ivs: Vec::new(), pos0: false });
    }
    chunks
}

/// Where the unplaced unmapped reads start (end of the last indexed chunk) and the file length.
pub fn tail_range(input: &Input, path: &str) -> (u64, u64) {
    let start = match input {
        // after the last placed record; with none (an empty BAM, or only unplaced reads) the first record, right
        // after the header (offset 0 would read the header as records)
        Input::Bam { bai, .. } => bai
            .refs
            .iter()
            .flat_map(|r| r.bins.values().flatten().map(|c| c.1))
            .max()
            .unwrap_or_else(|| crate::markdup::qname::open_input(path).first),
        Input::Cram { .. } => 0,
    };
    (start, std::fs::metadata(path).map(|m| m.len()).unwrap_or(0))
}

/// Record stream of a chunk; the tail chunk (tid < 0) is the unplaced unmapped records (always Some: before
/// 2026-11-17 a CRAM's tail was None and its unplaced reads were dropped from the output without a word).
pub fn open_chunk(input: &Input, ch: &Chunk, tail: (u64, u64), io_ns: &Arc<AtomicU64>) -> Option<Stream> {
    if ch.tid >= 0 {
        Some(input.stream(ch.tid as usize, &ch.query(), false, io_ns.clone()))
    } else {
        match input {
            Input::Bam { file, .. } => Some(Stream::Bam(bgzf::ChunkStream::new(file.clone(), vec![(tail.0, tail.1 << 16)], false, io_ns.clone()))),
            Input::Cram { pool, .. } => Some(Stream::Cram(pool.stream_unplaced(io_ns.clone()))),
        }
    }
}

/// Does record `b` (bytes after block_size) belong to chunk `ch`?
#[inline]
pub fn in_chunk(b: &[u8], ch: &Chunk) -> bool {
    let (ref_id, start, flag) = crate::bam::peek_core(b);
    if ch.tid < 0 {
        return ref_id == -1;
    }
    if ref_id != ch.tid {
        return false;
    }
    if start < 1 {
        return ch.pos0;
    }
    if start < ch.from || start > ch.to {
        return false;
    }
    let qend = if flag & 0x4 != 0 { start } else { crate::bam::query_end(b) };
    ch.overlaps(start, qend)
}

#[inline]
fn pseudo(u: usize) -> u64 {
    (((u / BLOCK_DATA) as u64) << 16) | (u % BLOCK_DATA) as u64
}

/// A finished chunk: compressed BGZF blocks plus index entries with pseudo virtual offsets.
#[derive(Default)]
pub struct ChunkOut {
    pub data: Vec<u8>,
    pub sizes: Vec<u32>,
    /// (beg 0-based, end exclusive, pseudo start, pseudo end, mapped)
    pub index: Vec<(i64, i64, u64, u64, bool)>,
    pub n_records: u64,
}

/// Per reference: too long for the BAI (> 2^29 bp), where htsjdk writes bin 0
pub type TooLarge = std::sync::Arc<Vec<bool>>;

pub fn too_large(ref_lens: &[i64]) -> TooLarge {
    std::sync::Arc::new(ref_lens.iter().map(|&l| l > 1 << 29).collect())
}

/// Accumulates the records of one chunk, compressing full blocks as it goes. Every record gets the bin field
/// htsjdk would write (`bam::htsjdk_bin`), as GATK / Picard re-encode the records they write.
pub struct RecWriter {
    level: i32,
    too_large: TooLarge,
    ubuf: Vec<u8>,
    ubase: usize,
    pub out: ChunkOut,
}

impl RecWriter {
    pub fn new(level: i32, too_large: &TooLarge) -> RecWriter {
        RecWriter { level, too_large: too_large.clone(), ubuf: Vec::with_capacity(4 << 20), ubase: 0, out: ChunkOut::default() }
    }

    fn flush(&mut self, all: bool) {
        let mut k = 0;
        while self.ubuf.len() - k >= BLOCK_DATA || (all && k < self.ubuf.len()) {
            let n = (self.ubuf.len() - k).min(BLOCK_DATA);
            let c = bgzf::deflate_block(&self.ubuf[k..k + n], self.level);
            self.out.data.extend_from_slice(&c);
            self.out.sizes.push(c.len() as u32);
            k += n;
        }
        self.ubuf.drain(..k);
        self.ubase += k;
    }

    /// Appends record `rec` (bytes after block_size); indexed unless it is an unplaced read.
    pub fn push(&mut self, rec: &[u8]) {
        let start = self.ubase + self.ubuf.len();
        self.ubuf.extend_from_slice(&(rec.len() as u32).to_le_bytes());
        let at = self.ubuf.len();
        self.ubuf.extend_from_slice(rec);
        let tid = i32::from_le_bytes(rec[0..4].try_into().unwrap());
        let big = tid >= 0 && self.too_large.get(tid as usize).copied().unwrap_or(false);
        let bin = crate::bam::htsjdk_bin(rec, big);
        self.ubuf[at + 10..at + 12].copy_from_slice(&bin.to_le_bytes());
        let end = self.ubase + self.ubuf.len();
        let (ref_id, pos1, flag) = crate::bam::peek_core(rec);
        if ref_id >= 0 && pos1 >= 1 {
            let beg = pos1 - 1;
            let mapped = flag & 0x4 == 0;
            let e = if mapped { crate::bam::query_end(rec).max(beg + 1) } else { beg + 1 };
            self.out.index.push((beg, e, pseudo(start), pseudo(end), mapped));
        }
        self.out.n_records += 1;
        if self.ubuf.len() >= 16 * BLOCK_DATA {
            self.flush(false);
        }
    }

    pub fn finish(mut self) -> ChunkOut {
        self.flush(true);
        self.out
    }
}

/// The output file: header, then chunks appended in order, then EOF block and `.bai`.
pub struct BamFile {
    file: std::io::BufWriter<std::fs::File>,
    offset: u64,
    idx: Vec<(usize, ContigIndex, u64, u64)>,
    n_no_coor: u64,
    names: Vec<String>,
    /// per reference: too long for a BAI bin (htsjdk writes bin 0 there)
    pub too_large: TooLarge,
    pub n_records: u64,
    /// MD5 of the file as written (`--create-output-bam-md5` / `CREATE_MD5_FILE`)
    md5: Option<crate::md5::Md5>,
    /// write the `.bai` (htsjdk indexes coordinate-sorted output only)
    pub index: bool,
}

impl BamFile {
    /// `md5`: also write `<path>.md5`, as htsjdk's `Md5CalculatingOutputStream` does.
    pub fn create(path: &str, text: &str, ref_names: &[String], ref_lens: &[i64], level: i32, md5: bool) -> BamFile {
        let mut hdr: Vec<u8> = Vec::new();
        hdr.extend_from_slice(b"BAM\x01");
        hdr.extend_from_slice(&(text.len() as i32).to_le_bytes());
        hdr.extend_from_slice(text.as_bytes());
        hdr.extend_from_slice(&(ref_names.len() as i32).to_le_bytes());
        for (n, &l) in ref_names.iter().zip(ref_lens) {
            hdr.extend_from_slice(&((n.len() + 1) as i32).to_le_bytes());
            hdr.extend_from_slice(n.as_bytes());
            hdr.push(0);
            hdr.extend_from_slice(&(l as i32).to_le_bytes());
        }
        let mut file = std::io::BufWriter::with_capacity(8 << 20, std::fs::File::create(path).unwrap_or_else(|e| panic!("{path}: {e}")));
        let mut offset = 0u64;
        let mut md5 = md5.then(crate::md5::Md5::default);
        for blk in hdr.chunks(BLOCK_DATA) {
            let c = bgzf::deflate_block(blk, level);
            file.write_all(&c).unwrap();
            if let Some(m) = &mut md5 {
                m.update(&c);
            }
            offset += c.len() as u64;
        }
        BamFile { file, offset, idx: Vec::new(), n_no_coor: 0, names: ref_names.to_vec(), too_large: too_large(ref_lens), n_records: 0, md5, index: true }
    }

    pub fn append(&mut self, tid: i32, o: &ChunkOut) {
        self.file.write_all(&o.data).unwrap();
        if let Some(m) = &mut self.md5 {
            m.update(&o.data);
        }
        let cum = bgzf::cumulative(self.offset, &o.sizes);
        if tid >= 0 {
            if !o.index.is_empty() {
                let t = tid as usize;
                if self.idx.last().map_or(true, |x| x.0 != t) {
                    self.idx.push((t, ContigIndex::new(&self.names[t]), 0, 0));
                }
                let entry = self.idx.last_mut().unwrap();
                let mut ci = ContigIndex::new(&self.names[t]);
                for &(b, e, s, en, mapped) in &o.index {
                    ci.push(b, e, bgzf::real_voff(s, &cum), bgzf::real_voff(en, &cum));
                    if mapped {
                        entry.2 += 1;
                    } else {
                        entry.3 += 1;
                    }
                }
                ci.finish();
                entry.1.merge(ci);
            }
        } else {
            self.n_no_coor += o.n_records;
        }
        self.offset = *cum.last().unwrap();
        self.n_records += o.n_records;
    }

    /// EOF block, flush, and the index (`<out minus .bam>.bai`, as GATK/Picard name it).
    pub fn finish(mut self, path: &str) {
        self.file.write_all(&bgzf::EOF_BLOCK).unwrap();
        self.file.flush().unwrap();
        // a stream (`-O /dev/stdout` piped to samtools, as nf-core's ApplyBQSR does): no `.bai` / `.md5` beside it
        if !std::fs::metadata(path).is_ok_and(|m| m.is_file()) {
            return;
        }
        if let Some(mut m) = self.md5.take() {
            m.update(&bgzf::EOF_BLOCK);
            std::fs::write(format!("{path}.md5"), m.hex()).unwrap_or_else(|e| panic!("{path}.md5: {e}"));
        }
        if self.index {
            tabix::write_bai(&format!("{}.bai", path.trim_end_matches(".bam")), self.names.len(), &self.idx, self.n_no_coor).unwrap();
        }
    }
}

/// The `@HD` line htsjdk writes (SAMTextHeaderCodec.writeHDLine without keepExistingVersionNumber): always first,
/// `VN:1.6` (SAMFileHeader.CURRENT_VERSION) whatever the input's version, then the input's other `@HD` tags in
/// their order; present even when the input has no `@HD` line.
fn htsjdk_hd_line(input_hd: Option<&str>) -> String {
    let mut l = String::from("@HD\tVN:1.6");
    if let Some(hd) = input_hd {
        for f in hd.split('\t').skip(1).filter(|f| !f.starts_with("VN:")) {
            l.push('\t');
            l.push_str(f);
        }
    }
    l
}

/// The header text htsjdk writes for `text` with `pg_lines` added (SAMTextHeaderCodec.encode): `@HD` as htsjdk
/// writes it, then the @SQ, @RG, @PG (the new ones last) and @CO lines, each kind in input order, whatever order
/// the input had them in.
pub fn header_with_pg(text: &str, pg_lines: &[String]) -> String {
    let text = text.trim_end_matches('\0');
    let mut out = htsjdk_hd_line(text.lines().find(|l| l.starts_with("@HD")));
    out.push('\n');
    for kind in ["@SQ", "@RG", "@PG", "@CO"] {
        for l in text.lines().filter(|l| l.starts_with(kind)) {
            out.push_str(&if kind == "@RG" { htsjdk_rg_line(l) } else { l.to_string() });
            out.push('\n');
        }
        if kind == "@PG" {
            for p in pg_lines {
                out.push_str(p);
                out.push('\n');
            }
        }
    }
    out
}

/// An `@RG` line as htsjdk writes it after reading it: `DT` is parsed (`DateParser`, ISO 8601) and written back as
/// `Iso8601Date.toString()`, i.e. `yyyy-MM-dd'T'HH:mm:ssZ` in the JVM's default time zone (that of the machine);
/// a value that does not parse is kept as it is.
pub fn htsjdk_rg_line(l: &str) -> String {
    if !l.starts_with("@RG\t") {
        return l.to_string();
    }
    l.split('\t')
        .map(|f| match f.strip_prefix("DT:").and_then(iso8601_epoch) {
            Some(t) => format!("DT:{}", local_iso8601(t)),
            None => f.to_string(),
        })
        .collect::<Vec<_>>()
        .join("\t")
}

/// `DateParser.parse` (UTC unless a zone is given), seconds since the epoch, milliseconds dropped.
fn iso8601_epoch(s: &str) -> Option<i64> {
    // StringTokenizer(isodate, "-T:.+Z", true)
    let mut toks: Vec<&str> = Vec::new();
    let mut st = 0;
    for (i, c) in s.char_indices() {
        if "-T:.+Z".contains(c) {
            if st < i {
                toks.push(&s[st..i]);
            }
            toks.push(&s[i..i + 1]);
            st = i + 1;
        }
    }
    if st < s.len() {
        toks.push(&s[st..]);
    }
    let mut it = toks.into_iter().peekable();
    let num = |t: Option<&str>| -> Option<i64> { t?.parse::<i64>().ok() };
    // check(st, token): false at the end, an error for another token
    fn check<'a>(it: &mut std::iter::Peekable<impl Iterator<Item = &'a str>>, t: &str) -> Option<bool> {
        match it.next() {
            None => Some(false),
            Some(x) if x == t => Some(true),
            Some(_) => None,
        }
    }
    let (mut y, mut mo, mut d, mut h, mut mi, mut se) = (1970, 0, 1, 0, 0, 0);
    let mut zone: i64 = 0;
    'p: {
        match it.next() {
            Some(t) => y = t.parse().ok()?,
            None => break 'p,
        }
        if check(&mut it, "-")? && it.peek().is_some() {
            mo = num(it.next())? - 1;
        } else {
            break 'p;
        }
        if check(&mut it, "-")? && it.peek().is_some() {
            d = num(it.next())?;
        } else {
            break 'p;
        }
        if check(&mut it, "T")? && it.peek().is_some() {
            h = num(it.next())?;
        } else {
            break 'p;
        }
        if check(&mut it, ":")? && it.peek().is_some() {
            mi = num(it.next())?;
        } else {
            break 'p;
        }
        let Some(mut tok) = it.next() else { break 'p };
        if tok == ":" {
            se = num(it.next())?;
            let Some(t) = it.next() else { break 'p };
            tok = t;
            if tok == "." {
                it.next()?; // milliseconds (truncated by Iso8601Date)
                let Some(t) = it.next() else { break 'p };
                tok = t;
            }
        }
        if tok != "Z" {
            if tok != "+" && tok != "-" {
                return None;
            }
            let mut tzh = num(it.next())?;
            let mut tzm = 0;
            if check(&mut it, ":")? && it.peek().is_some() {
                tzm = num(it.next())?;
            } else if tzh >= 100 {
                tzm = tzh % 100;
                tzh /= 100;
            }
            zone = (tzh * 60 + tzm) * 60 * if tok == "+" { 1 } else { -1 };
        }
    }
    // lenient GregorianCalendar fields (proleptic Gregorian here)
    let (yy, mm) = (y + mo.div_euclid(12), mo.rem_euclid(12) + 1);
    let days = days_from_civil(yy, mm, 1) + d - 1;
    Some(days * 86400 + h * 3600 + mi * 60 + se - zone)
}

fn days_from_civil(y: i64, m: i64, d: i64) -> i64 {
    let y = if m <= 2 { y - 1 } else { y };
    let era = y.div_euclid(400);
    let yoe = y - era * 400;
    let doy = (153 * (if m > 2 { m - 3 } else { m + 9 }) + 2) / 5 + d - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146097 + doe - 719468
}

#[repr(C)]
struct Tm {
    sec: i32,
    min: i32,
    hour: i32,
    mday: i32,
    mon: i32,
    year: i32,
    wday: i32,
    yday: i32,
    isdst: i32,
    gmtoff: i64,
    zone: *const u8,
}

extern "C" {
    fn localtime_r(t: *const i64, out: *mut Tm) -> *mut Tm;
}

/// `new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ssZ").format(date)` in the local time zone (glibc's, which reads
/// `TZ` / `/etc/localtime` as the JVM does).
fn local_iso8601(t: i64) -> String {
    let mut tm = Tm { sec: 0, min: 0, hour: 0, mday: 0, mon: 0, year: 0, wday: 0, yday: 0, isdst: 0, gmtoff: 0, zone: std::ptr::null() };
    // SAFETY: localtime_r only writes the struct we pass
    let ok = unsafe { !localtime_r(&t, &mut tm).is_null() };
    assert!(ok, "localtime_r failed");
    let off = tm.gmtoff / 60;
    format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:02}{}{:02}{:02}",
        tm.year + 1900,
        tm.mon + 1,
        tm.mday,
        tm.hour,
        tm.min,
        tm.sec,
        if off < 0 { '-' } else { '+' },
        off.abs() / 60,
        off.abs() % 60
    )
}

/// Reusable buffers for [`normalize_tags`].
#[derive(Default)]
pub struct TagScratch {
    list: Vec<(u16, usize, usize, i64)>,
    out: Vec<u8>,
}

#[inline]
fn tag_size(t: u8) -> usize {
    match t {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        _ => 0,
    }
}

/// Offset of the aux data in record `rec` (bytes after block_size).
#[inline]
pub fn aux_start(rec: &[u8]) -> usize {
    let l_read_name = rec[8] as usize;
    let n_cigar = u16::from_le_bytes([rec[12], rec[13]]) as usize;
    let l_seq = i32::from_le_bytes(rec[16..20].try_into().unwrap()) as usize;
    32 + l_read_name + 4 * n_cigar + (l_seq + 1) / 2 + l_seq
}

/// htsjdk's re-encoding of a modified record's tags: decoded into a list kept sorted by binary tag
/// (`(c1 << 8) | c0`, a repeated tag replacing the earlier value), `remove` tags cleared, `set_z`
/// set as a string tag, and integers retyped (`BinaryTagCodec.getIntegerType`: the smallest of c, C,
/// s, S, i, I, signed first).
pub fn normalize_tags(rec: &mut Vec<u8>, aux: usize, remove: &[[u8; 2]], set_z: Option<([u8; 2], &[u8])>, s: &mut TagScratch) {
    s.list.clear();
    let mut o = aux;
    while o + 3 <= rec.len() {
        let st = o;
        let tag = [rec[o], rec[o + 1]];
        let ty = rec[o + 2];
        o += 3;
        let mut int = i64::MIN;
        match ty {
            b'Z' | b'H' => {
                while rec[o] != 0 {
                    o += 1;
                }
                o += 1;
            }
            b'B' => {
                let sub = rec[o];
                let n = u32::from_le_bytes(rec[o + 1..o + 5].try_into().unwrap()) as usize;
                o += 5 + n * tag_size(sub);
            }
            b'c' => int = rec[o] as i8 as i64,
            b'C' => int = rec[o] as i64,
            b's' => int = i16::from_le_bytes([rec[o], rec[o + 1]]) as i64,
            b'S' => int = u16::from_le_bytes([rec[o], rec[o + 1]]) as i64,
            b'i' => int = i32::from_le_bytes(rec[o..o + 4].try_into().unwrap()) as i64,
            b'I' => int = u32::from_le_bytes(rec[o..o + 4].try_into().unwrap()) as i64,
            _ => {}
        }
        if !matches!(ty, b'Z' | b'H' | b'B') {
            o += tag_size(ty);
        }
        if remove.contains(&tag) || set_z.map_or(false, |(t, _)| t == tag) {
            continue;
        }
        let key = ((tag[1] as u16) << 8) | tag[0] as u16;
        if let Some(p) = s.list.iter().position(|x| x.0 == key) {
            s.list.remove(p);
        }
        s.list.push((key, st, o, int));
    }
    // the string tag to set goes in as a pseudo entry (start = usize::MAX)
    if let Some((t, _)) = set_z {
        s.list.push((((t[1] as u16) << 8) | t[0] as u16, usize::MAX, 0, i64::MIN));
    }
    s.list.sort_by_key(|x| x.0);
    s.out.clear();
    for &(_, st, e, int) in s.list.iter() {
        if st == usize::MAX {
            let (t, v) = set_z.unwrap();
            s.out.extend_from_slice(&t);
            s.out.push(b'Z');
            s.out.extend_from_slice(v);
            s.out.push(0);
        } else if int == i64::MIN && rec[st + 2] == b'H' {
            // htsjdk decodes H to a byte[] (hexStringToBytes: an odd length is an error), which it writes back as a
            // signed byte array, B:c
            let v = &rec[st + 3..e - 1];
            assert!(v.len() % 2 == 0, "Hex representation of byte string does not have even number of hex chars: {}", String::from_utf8_lossy(v));
            let hex = |c: u8| (c as char).to_digit(16).unwrap_or_else(|| panic!("H tag: not a hex digit: {}", c as char)) as u8;
            s.out.extend_from_slice(&rec[st..st + 2]);
            s.out.extend_from_slice(b"Bc");
            s.out.extend_from_slice(&((v.len() / 2) as u32).to_le_bytes());
            s.out.extend(v.chunks(2).map(|p| hex(p[0]) << 4 | hex(p[1])));
        } else if int == i64::MIN {
            s.out.extend_from_slice(&rec[st..e]);
        } else {
            s.out.extend_from_slice(&rec[st..st + 2]);
            let v = int;
            if v > i32::MAX as i64 {
                s.out.push(b'I');
                s.out.extend_from_slice(&(v as u32).to_le_bytes());
            } else if v > 65535 {
                s.out.push(b'i');
                s.out.extend_from_slice(&(v as i32).to_le_bytes());
            } else if v > 32767 {
                s.out.push(b'S');
                s.out.extend_from_slice(&(v as u16).to_le_bytes());
            } else if v > 255 {
                s.out.push(b's');
                s.out.extend_from_slice(&(v as i16).to_le_bytes());
            } else if v > 127 {
                s.out.push(b'C');
                s.out.push(v as u8);
            } else if v >= -128 {
                s.out.push(b'c');
                s.out.push(v as i8 as u8);
            } else if v >= -32768 {
                s.out.push(b's');
                s.out.extend_from_slice(&(v as i16).to_le_bytes());
            } else {
                s.out.push(b'i');
                s.out.extend_from_slice(&(v as i32).to_le_bytes());
            }
        }
    }
    rec.truncate(aux);
    rec.extend_from_slice(&s.out);
}

#[cfg(test)]
mod date_tests {
    use super::*;
    #[test]
    fn iso8601() {
        assert_eq!(iso8601_epoch("2012-02-24T05:00:00+0000"), Some(1330059600));
        assert_eq!(iso8601_epoch("2012-02-24T13:00:00+08:00"), Some(1330059600));
        assert_eq!(iso8601_epoch("2012-02-24"), Some(1330041600));
        assert_eq!(iso8601_epoch("2012-02-24T05:00:00.123Z"), Some(1330059600));
        assert_eq!(iso8601_epoch("2012-02-24T00:00:00-0500"), Some(1330059600));
        assert_eq!(iso8601_epoch("yesterday"), None);
    }
}
