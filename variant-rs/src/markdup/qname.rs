//! MarkDuplicates on query-grouped input (`ASSUME_SORT_ORDER=queryname`, or a `SO:queryname` header), as WARP
//! runs it on MergeBamAlignment output; several inputs are merged first.
//!
//! - Inputs: htsjdk's `SamFileHeaderMerger` (identical sequence dictionaries; @PG trees and @RG merged, colliding
//!   IDs renamed `ID.<base36>`, both sorted by ID; comments concatenated) and `MergingSamRecordIterator` (a
//!   `java.util.PriorityQueue` of the inputs ordered by `SAMRecordQueryNameComparator`, RG / PG tags translated
//!   when IDs collided). One input is read as is.
//! - Pass 1 (`buildSortedReadEndLists`): the read ends of every primary mapped record carry the index of the
//!   first record of its query-name group; ends go to the same buckets and duplicate sets as in coordinate mode.
//! - Pass 2: Picard's walk over the sorted duplicate indexes (`nextIndexIfNeeded`): every record of a duplicate
//!   group is flagged, secondary / supplementary / unmapped ones included. Output in input order, `@HD
//!   SO:unknown GO:query` when the order was assumed, no index.

use super::*;
use crate::bam::Header;
use crate::bgzf::ChunkStream;
use std::cmp::Ordering as Cmp;
use std::collections::HashSet;
use std::fs::File;

// ------------------------------------------------------------------------------------------------
// inputs and htsjdk's header merger
// ------------------------------------------------------------------------------------------------

pub(crate) struct InFile {
    pub(crate) file: Arc<File>,
    pub(crate) header: Header,
    /// virtual offset of the first record, file length
    pub(crate) first: u64,
    pub(crate) len: u64,
}

pub(crate) fn open_input(path: &str) -> InFile {
    let file = File::open(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut magic = [0u8; 4];
    std::os::unix::fs::FileExt::read_exact_at(&file, &mut magic, 0).unwrap_or_else(|e| panic!("{path}: {e}"));
    assert!(&magic != b"CRAM", "{path}: variant-rs MarkDuplicates reads query-grouped input from BAM only");
    let header = Header::read(&file);
    // uncompressed header length: magic, l_text, text, n_ref, (l_name, name, l_ref) per reference
    let hb = crate::bgzf::read_prefix(&file, |b| b.len() >= 8 && b.len() >= 12 + i32::from_le_bytes(b[4..8].try_into().unwrap()) as usize);
    let l_text = i32::from_le_bytes(hb[4..8].try_into().unwrap()) as u64;
    let hlen = 8 + l_text + 4 + header.ref_names.iter().map(|n| 4 + n.len() as u64 + 1 + 4).sum::<u64>();
    let first = crate::bgzf::voff_at(&file, hlen);
    let len = file.metadata().map(|m| m.len()).unwrap_or(0);
    InFile { file: Arc::new(file), header, first, len }
}

/// A header record (@RG / @PG): ID and the other attributes in order.
#[derive(Clone, Debug)]
struct HRec {
    id: String,
    attrs: Vec<(String, String)>,
}

impl HRec {
    fn parse(line: &str) -> HRec {
        let mut id = String::new();
        let mut attrs = Vec::new();
        for f in line.split('\t').skip(1) {
            let (k, v) = f.split_at(2.min(f.len()));
            let v = v.strip_prefix(':').unwrap_or(v);
            if k == "ID" {
                id = v.to_string();
            } else {
                attrs.push((k.to_string(), v.to_string()));
            }
        }
        HRec { id, attrs }
    }
    fn get(&self, k: &str) -> Option<&str> {
        self.attrs.iter().find(|a| a.0 == k).map(|a| a.1.as_str())
    }
    fn set(&mut self, k: &str, v: &str) {
        match self.attrs.iter_mut().find(|a| a.0 == k) {
            Some(a) => a.1 = v.to_string(),
            None => self.attrs.push((k.to_string(), v.to_string())),
        }
    }
    /// `equals` of SAMReadGroupRecord / SAMProgramRecord: same ID and the same attribute map
    fn same(&self, o: &HRec) -> bool {
        let mut a = self.attrs.clone();
        let mut b = o.attrs.clone();
        a.sort();
        b.sort();
        self.id == o.id && a == b
    }
    fn line(&self, kind: &str) -> String {
        let mut l = format!("@{kind}\tID:{}", self.id);
        for (k, v) in &self.attrs {
            l.push_str(&format!("\t{k}:{v}"));
        }
        l
    }
}

/// `SamFileHeaderMerger.positiveFourDigitBase36Str`
fn base36(mut v: i32) -> String {
    if v == 0 {
        return "0".into();
    }
    let mut s = Vec::new();
    while v > 0 {
        s.push(b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"[(v % 36) as usize]);
        v /= 36;
    }
    s.reverse();
    String::from_utf8(s).unwrap()
}

/// `mergeHeaderRecords`: `recs` = (input, record); returns whether IDs collided.
fn merge_records(recs: &[(usize, HRec)], taken: &mut HashSet<String>, counter: &mut i32, tr: &mut [HashMap<String, String>], result: &mut Vec<HRec>) -> bool {
    // id -> unique records (first-seen order) -> inputs
    let mut by_id: Vec<(String, Vec<(HRec, Vec<usize>)>)> = Vec::new();
    for (f, r) in recs {
        let slot = match by_id.iter().position(|x| x.0 == r.id) {
            Some(p) => p,
            None => {
                by_id.push((r.id.clone(), Vec::new()));
                by_id.len() - 1
            }
        };
        let list = &mut by_id[slot].1;
        match list.iter_mut().find(|x| x.0.same(r)) {
            Some(x) => x.1.push(*f),
            None => list.push((r.clone(), vec![*f])),
        }
    }
    let mut coll = false;
    for (id, uniq) in by_id {
        for (r, files) in uniq {
            let new_id = if !taken.contains(&id) {
                taken.insert(id.clone());
                *counter += 1;
                id.clone()
            } else {
                coll = true;
                let n = loop {
                    let c = format!("{id}.{}", base36(*counter));
                    *counter += 1;
                    if !taken.contains(&c) {
                        break c;
                    }
                };
                taken.insert(n.clone());
                n
            };
            for f in files {
                tr[f].insert(id.clone(), new_id.clone());
            }
            let mut nr = r.clone();
            nr.id = new_id;
            result.push(nr);
        }
    }
    coll
}

/// `translateIds`
fn translate_ids(list: Vec<(usize, HRec)>, tr: &[HashMap<String, String>], pp: bool) -> Vec<(usize, HRec)> {
    list.into_iter()
        .map(|(f, mut r)| {
            if let Some(t) = tr[f].get(&r.id) {
                r.id = t.clone();
            }
            if pp {
                if let Some(t) = r.get("PP").and_then(|p| tr[f].get(p)).cloned() {
                    r.set("PP", &t);
                }
            }
            (f, r)
        })
        .collect()
}

/// The merged input: header text (before MarkDuplicates' own changes) and per-input RG / PG translations.
pub(crate) struct Merged {
    pub(crate) text: String,
    ref_names: Vec<String>,
    ref_lens: Vec<i64>,
    rg_tr: Option<Vec<HashMap<Vec<u8>, Vec<u8>>>>,
    pg_tr: Option<Vec<HashMap<Vec<u8>, Vec<u8>>>>,
}

fn merge_headers(ins: &[InFile], sort_order: &str) -> Merged {
    let hs: Vec<&Header> = ins.iter().map(|i| &i.header).collect();
    merge_header_list(&hs, sort_order)
}

/// htsjdk's `SamFileHeaderMerger(sortOrder, headers, mergeDictionaries)` for identical dictionaries.
pub(crate) fn merge_header_list(hs: &[&Header], sort_order: &str) -> Merged {
    struct In<'a> {
        header: &'a Header,
    }
    let ins: Vec<In> = hs.iter().map(|&h| In { header: h }).collect();
    let n = ins.len();
    let lines = |i: usize, kind: &str| -> Vec<String> { ins[i].header.text.lines().filter(|l| l.starts_with(kind)).map(String::from).collect() };
    // sequence dictionaries must agree (SequenceUtil.assertSequenceDictionariesEqual)
    for i in 1..n {
        assert!(
            ins[i].header.ref_names == ins[0].header.ref_names && ins[i].header.ref_lens == ins[0].header.ref_lens,
            "the inputs' sequence dictionaries differ (Picard requires them to be identical)"
        );
    }
    // program groups: tree roots first, then their children, level by level
    let mut taken = HashSet::new();
    let mut counter = 0i32;
    let mut pg_tr_s: Vec<HashMap<String, String>> = vec![HashMap::new(); n];
    let mut left: Vec<(usize, HRec)> = Vec::new();
    for i in 0..n {
        let mut ids = HashSet::new();
        for l in lines(i, "@PG\t") {
            let r = HRec::parse(&l);
            assert!(ids.insert(r.id.clone()), "input {i} contains more than one PG with the same id ({})", r.id);
            left.push((i, r));
        }
    }
    let (mut current, rest): (Vec<_>, Vec<_>) = left.into_iter().partition(|(_, r)| r.get("PP").is_none());
    let mut left = rest;
    let mut pgs: Vec<HRec> = Vec::new();
    let mut pg_coll = false;
    while !current.is_empty() {
        let mut cur_result = Vec::new();
        pg_coll |= merge_records(&current, &mut taken, &mut counter, &mut pg_tr_s, &mut cur_result);
        pgs.extend(cur_result);
        current = translate_ids(current, &pg_tr_s, false);
        left = translate_ids(left, &pg_tr_s, true);
        let mut next = Vec::new();
        let mut keep = Vec::new();
        for (f, r) in left {
            let pp = r.get("PP").unwrap().to_string();
            if current.iter().any(|(cf, cr)| *cf == f && cr.id == pp) {
                next.push((f, r));
            } else {
                keep.push((f, r));
            }
        }
        left = keep;
        current = next;
    }
    assert!(left.is_empty(), "{} program groups weren't processed. Do their PP ids point to existing PGs?", left.len());
    pgs.sort_by(|a, b| a.id.cmp(&b.id));
    // read groups
    let mut rgs_in: Vec<(usize, HRec)> = Vec::new();
    for i in 0..n {
        let mut ids = HashSet::new();
        for l in lines(i, "@RG\t") {
            let r = HRec::parse(&l);
            assert!(ids.insert(r.id.clone()), "input {i} contains more than one RG with the same id ({})", r.id);
            rgs_in.push((i, r));
        }
    }
    let mut rg_tr_s: Vec<HashMap<String, String>> = vec![HashMap::new(); n];
    let mut rgs = Vec::new();
    let mut taken = HashSet::new();
    let mut counter = 0i32;
    let rg_coll = merge_records(&rgs_in, &mut taken, &mut counter, &mut rg_tr_s, &mut rgs);
    rgs.sort_by(|a, b| a.id.cmp(&b.id));
    // new SAMFileHeader(): VN; setGroupOrder(none); setSortOrder(sortOrder)
    let mut text = format!("@HD\tVN:1.6\tGO:none\tSO:{sort_order}\n");
    for l in lines(0, "@SQ\t") {
        text.push_str(&l);
        text.push('\n');
    }
    for r in &rgs {
        text.push_str(&r.line("RG"));
        text.push('\n');
    }
    for r in &pgs {
        text.push_str(&r.line("PG"));
        text.push('\n');
    }
    for i in 0..n {
        for l in lines(i, "@CO\t") {
            text.push_str(&l);
            text.push('\n');
        }
    }
    let bytes = |t: Vec<HashMap<String, String>>| -> Vec<HashMap<Vec<u8>, Vec<u8>>> {
        t.into_iter().map(|m| m.into_iter().map(|(k, v)| (k.into_bytes(), v.into_bytes())).collect()).collect()
    };
    Merged {
        text,
        ref_names: ins[0].header.ref_names.clone(),
        ref_lens: ins[0].header.ref_lens.clone(),
        rg_tr: rg_coll.then(|| bytes(rg_tr_s)),
        pg_tr: pg_coll.then(|| bytes(pg_tr_s)),
    }
}

/// htsjdk's in-place updates of the `@HD` attributes (a LinkedHashMap): existing keys keep their position.
fn hd_put(text: &str, kv: &[(&str, &str)]) -> String {
    let hd = text.lines().find(|l| l.starts_with("@HD"));
    let mut attrs: Vec<(String, String)> = match hd {
        Some(l) => l.split('\t').skip(1).map(|f| (f[..2.min(f.len())].to_string(), f.get(3..).unwrap_or("").to_string())).collect(),
        None => vec![("VN".into(), "1.6".into())],
    };
    for (k, v) in kv {
        match attrs.iter_mut().find(|a| a.0 == *k) {
            Some(a) => a.1 = v.to_string(),
            None => attrs.push((k.to_string(), v.to_string())),
        }
    }
    let mut line = String::from("@HD");
    for (k, v) in attrs {
        line.push_str(&format!("\t{k}:{v}"));
    }
    let mut out = format!("{line}\n");
    for l in text.trim_end_matches('\0').lines().filter(|l| !l.starts_with("@HD")) {
        out.push_str(l);
        out.push('\n');
    }
    out
}

// ------------------------------------------------------------------------------------------------
// the merged record stream
// ------------------------------------------------------------------------------------------------

fn rec_name(r: &[u8]) -> &[u8] {
    &r[32..32 + r[8] as usize - 1]
}

fn rec_flag(r: &[u8]) -> u16 {
    u16::from_le_bytes([r[14], r[15]])
}

/// Integer value of tag `tag` (any integer type), as `getIntegerAttribute`.
fn find_int(rec: &[u8], tag: &[u8; 2]) -> Option<i64> {
    let mut o = bamwrite::aux_start(rec);
    while o + 3 <= rec.len() {
        let t = [rec[o], rec[o + 1]];
        let ty = rec[o + 2];
        o += 3;
        let v: Option<i64> = match ty {
            b'c' => Some(rec[o] as i8 as i64),
            b'C' => Some(rec[o] as i64),
            b's' => Some(i16::from_le_bytes([rec[o], rec[o + 1]]) as i64),
            b'S' => Some(u16::from_le_bytes([rec[o], rec[o + 1]]) as i64),
            b'i' => Some(i32::from_le_bytes(rec[o..o + 4].try_into().unwrap()) as i64),
            b'I' => Some(u32::from_le_bytes(rec[o..o + 4].try_into().unwrap()) as i64),
            _ => None,
        };
        if t == *tag && v.is_some() {
            return v;
        }
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
                let sz = match sub {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    _ => 4,
                };
                o += 5 + n * sz;
            }
            b'A' | b'c' | b'C' => o += 1,
            b's' | b'S' => o += 2,
            _ => o += 4,
        }
    }
    None
}

/// `SAMRecordQueryNameComparator.compare`
fn qname_cmp(a: &[u8], b: &[u8]) -> Cmp {
    let c = rec_name(a).cmp(rec_name(b));
    if c != Cmp::Equal {
        return c;
    }
    let (fa, fb) = (rec_flag(a), rec_flag(b));
    let (pa, pb) = (fa & 0x1 != 0, fb & 0x1 != 0);
    if pa || pb {
        if !pa {
            return Cmp::Greater;
        } else if !pb {
            return Cmp::Less;
        } else if fa & 0x40 != 0 && fb & 0x80 != 0 {
            return Cmp::Less;
        } else if fa & 0x80 != 0 && fb & 0x40 != 0 {
            return Cmp::Greater;
        }
    }
    if (fa & 0x10 != 0) != (fb & 0x10 != 0) {
        return if fa & 0x10 != 0 { Cmp::Greater } else { Cmp::Less };
    }
    if (fa & 0x100 != 0) != (fb & 0x100 != 0) {
        return if fb & 0x100 != 0 { Cmp::Less } else { Cmp::Greater };
    }
    if (fa & 0x800 != 0) != (fb & 0x800 != 0) {
        return if fb & 0x800 != 0 { Cmp::Less } else { Cmp::Greater };
    }
    match (find_int(a, b"HI"), find_int(b, b"HI")) {
        (Some(x), Some(y)) => x.cmp(&y),
        (Some(_), None) => Cmp::Greater,
        (None, Some(_)) => Cmp::Less,
        (None, None) => Cmp::Equal,
    }
}

/// Replaces (or adds) string tag `tag` of a record; htsjdk re-encodes all tags later anyway.
fn set_z(rec: &mut Vec<u8>, tag: &[u8; 2], val: &[u8]) {
    let aux = bamwrite::aux_start(rec);
    let mut ts = bamwrite::TagScratch::default();
    bamwrite::normalize_tags(rec, aux, &[], Some((*tag, val)), &mut ts);
}

/// The records of all inputs in htsjdk's merge order (one input: file order).
struct Source<'a> {
    streams: Vec<ChunkStream>,
    heads: Vec<Vec<u8>>,
    /// java.util.PriorityQueue of input indexes (binary heap, same sift rules)
    pq: Vec<usize>,
    m: &'a Merged,
}

impl<'a> Source<'a> {
    fn new(ins: &[InFile], m: &'a Merged, io_ns: &Arc<AtomicU64>) -> Source<'a> {
        let mut s = Source { streams: Vec::new(), heads: Vec::new(), pq: Vec::new(), m };
        for f in ins {
            s.streams.push(ChunkStream::new(f.file.clone(), vec![(f.first, f.len << 16)], true, io_ns.clone()));
            s.heads.push(Vec::new());
        }
        for i in 0..ins.len() {
            if s.advance(i) {
                s.offer(i);
            }
        }
        s
    }
    fn advance(&mut self, i: usize) -> bool {
        let mut buf = std::mem::take(&mut self.heads[i]);
        let ok = self.streams[i].next_record(&mut buf);
        self.heads[i] = buf;
        ok
    }
    fn cmp(&self, a: usize, b: usize) -> Cmp {
        qname_cmp(&self.heads[a], &self.heads[b])
    }
    fn offer(&mut self, x: usize) {
        let mut k = self.pq.len();
        self.pq.push(x);
        while k > 0 {
            let parent = (k - 1) >> 1;
            let e = self.pq[parent];
            if self.cmp(x, e) != Cmp::Less {
                break;
            }
            self.pq[k] = e;
            k = parent;
        }
        self.pq[k] = x;
    }
    fn poll(&mut self) -> Option<usize> {
        if self.pq.is_empty() {
            return None;
        }
        let result = self.pq[0];
        let x = self.pq.pop().unwrap();
        let n = self.pq.len();
        if n > 0 {
            let mut k = 0;
            let half = n >> 1;
            while k < half {
                let mut child = 2 * k + 1;
                let mut c = self.pq[child];
                let right = child + 1;
                if right < n && self.cmp(c, self.pq[right]) == Cmp::Greater {
                    child = right;
                    c = self.pq[child];
                }
                if self.cmp(x, c) != Cmp::Greater {
                    break;
                }
                self.pq[k] = c;
                k = child;
            }
            self.pq[k] = x;
        }
        Some(result)
    }
    /// The next record (bytes after block_size) into `out`.
    fn next(&mut self, out: &mut Vec<u8>) -> bool {
        let Some(i) = self.poll() else { return false };
        std::mem::swap(out, &mut self.heads[i]);
        if self.advance(i) {
            self.offer(i);
        }
        translate_record(self.m, i, out);
        true
    }
}

/// `MergingSamRecordIterator`: the RG / PG tags of a record from input `i` translated to the merged header's IDs
/// (only when some ID collided).
pub(crate) fn translate_record(m: &Merged, i: usize, out: &mut Vec<u8>) {
    let tr = |t: &Option<Vec<HashMap<Vec<u8>, Vec<u8>>>>, tag: &[u8; 2], out: &mut Vec<u8>| {
        if let Some(t) = t {
            let aux = bamwrite::aux_start(out);
            if let Some(old) = find_z(out, aux, tag) {
                // getReadGroupId / getProgramGroupId: null for an ID not in the input's header, which
                // setAttribute(tag, null) removes
                match t[i].get(old).cloned() {
                    Some(new) => set_z(out, tag, &new),
                    None => {
                        let mut ts = bamwrite::TagScratch::default();
                        bamwrite::normalize_tags(out, aux, &[*tag], None, &mut ts);
                    }
                }
            }
        }
    };
    tr(&m.rg_tr, b"RG", out);
    tr(&m.pg_tr, b"PG", out);
}

// ------------------------------------------------------------------------------------------------
// driver
// ------------------------------------------------------------------------------------------------

pub(super) fn run(args: &Args) {
    let t0 = Instant::now();
    let ins: Vec<InFile> = args.inputs.iter().map(|p| open_input(p)).collect();
    let assumed = args.assume.as_deref() == Some("queryname");
    // openInputs: an assumed order is written into the first header
    let first_text = if assumed { hd_put(&ins[0].header.text, &[("SO", "queryname")]) } else { ins[0].header.text.clone() };
    let merged = if ins.len() == 1 {
        Merged { text: first_text, ref_names: ins[0].header.ref_names.clone(), ref_lens: ins[0].header.ref_lens.clone(), rg_tr: None, pg_tr: None }
    } else {
        if !assumed {
            for f in &ins {
                assert!(sort_order(&f.header.text) == "queryname", "Files are not compatible with sort order");
            }
        }
        merge_headers(&ins, "queryname")
    };
    let header_text = merged.text.clone();
    let hi = HeaderInfo::new(&header_text);
    let io_ns = Arc::new(AtomicU64::new(0));
    let bmap = BucketMap::new(&merged.ref_lens, args.bucket_bp);
    let pid = std::process::id();
    let frag_b: Buckets<FragDisk> = Buckets::new(&args.tmp, &format!("ardmdq{pid}.frag"), bmap.n);
    let pair_b: Buckets<PairDisk> = Buckets::new(&args.tmp, &format!("ardmdq{pid}.pair"), bmap.n);

    // ---- pass 1: read ends indexed by the first record of their query-name group
    let mut pg_seen: Vec<Option<Vec<u8>>> = Vec::new();
    let mut pg_set: HashSet<Option<Vec<u8>>> = HashSet::new();
    let mut total = 0u64;
    {
        let mut src = Source::new(&ins, &merged, &io_ns);
        let mut rec = Vec::new();
        let mut lf = frag_b.local();
        let mut lp = pair_b.local();
        let mut waiting: HashMap<(Vec<u8>, i32), Half> = HashMap::new();
        let mut group_name: Vec<u8> = Vec::new();
        let mut group_idx = 0u64;
        let mut first = true;
        while src.next(&mut rec) {
            let idx = total;
            total += 1;
            let aux = bamwrite::aux_start(&rec);
            let pg = find_z(&rec, aux, b"PG").map(|v| v.to_vec());
            if !pg_set.contains(&pg) {
                pg_set.insert(pg.clone());
                pg_seen.push(pg);
            }
            let name = rec_name(&rec);
            if first || name != group_name.as_slice() {
                group_name.clear();
                group_name.extend_from_slice(name);
                group_idx = idx;
                first = false;
            }
            let v = view(&rec);
            if v.flag & 0x4 != 0 || v.flag & 0x900 != 0 {
                continue;
            }
            let rg = find_z(&rec, aux, b"RG");
            let (lib, rg_ord) = hi.lib_and_rg(rg);
            let mut e = End { lib, orientation: if v.flag & 0x10 != 0 { R } else { F }, ofo: u8::MAX, rg: -1, tile: -1, ref1: v.ref_id, coord1: v.unclipped_5p as i32, ref2: -1, coord2: -1, x: -1, y: -1, score: v.score, _pad: [0; 6], idx1: group_idx, idx2: u64::MAX };
            let paired = v.flag & 0x1 != 0 && v.flag & 0x8 == 0;
            if paired {
                e.ref2 = v.mate_ref;
            }
            if let Some((tile, x, y)) = parse_location(name) {
                e.tile = tile;
                e.x = x;
                e.y = y;
                e.rg = rg_ord;
            }
            Buckets::push(&mut lf, bmap.bucket(e.ref1, e.coord1), &e);
            if paired {
                let mut key = rg.expect("paired read without a read group").to_vec();
                key.extend_from_slice(name);
                let h = Half { end: e, neg: v.flag & 0x10 != 0, first_of_pair: v.flag & 0x40 != 0, target_ref: v.mate_ref };
                match waiting.remove(&(key.clone(), v.ref_id)) {
                    Some(a) => {
                        let pe = make_pair(&a, &h);
                        Buckets::push(&mut lp, bmap.bucket(pe.ref1, pe.coord1), &pe);
                    }
                    None => {
                        waiting.insert((key, v.mate_ref), h);
                    }
                }
            }
            if total % 4096 == 0 {
                frag_b.flush_local(&mut lf, false);
                pair_b.flush_local(&mut lp, false);
            }
        }
        frag_b.flush_local(&mut lf, true);
        pair_b.flush_local(&mut lp, true);
        eprintln!("pass 1: {total} records ({} inputs), {} pairs never matched, in {:.1}s", ins.len(), waiting.len(), t0.elapsed().as_secs_f64());
    }
    frag_b.finish();
    pair_b.finish();

    // ---- duplicate sets (as in coordinate mode)
    let t1 = Instant::now();
    let res = (0..bmap.n)
        .into_par_iter()
        .map(|b| group_from_disk(&frag_b, &pair_b, b, args.optical_distance))
        .reduce(GroupResult::default, |mut a, b| {
            a.merge(b);
            a
        });
    drop(frag_b);
    drop(pair_b);
    let mut dups = res.dups.clone();
    dups.par_sort_unstable();
    eprintln!("duplicates: {} indexes in {:.1}s", dups.len(), t1.elapsed().as_secs_f64());

    // ---- pass 2: flags in input order (decided sequentially), BGZF compression in parallel
    let t2 = Instant::now();
    let cmdline = std::env::args().collect::<Vec<_>>().join(" ");
    let (pg_map, pg_lines) = chain_pg_ids(&pg_seen, &header_text, &cmdline);
    // MarkDuplicates.doWork: an assumed queryname order is written as SO:unknown GO:query
    let out_text = if assumed { hd_put(&header_text, &[("GO", "query"), ("SO", "unknown")]) } else { header_text.clone() };
    let text = bamwrite::header_with_pg(&out_text, &pg_lines);
    let mut bam = bamwrite::BamFile::create(&args.out, &text, &merged.ref_names, &merged.ref_lens, args.level, args.md5);
    let big = bam.too_large.clone();
    bam.index = false;
    let n_threads = args.threads.max(1);
    let lm: Vec<LibMetrics> = std::thread::scope(|sc| {
        let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<u8>>(2 * n_threads);
        let producer = sc.spawn(|| {
            let tx = tx;
            let mut src = Source::new(&ins, &merged, &io_ns);
            let mut lm = vec![LibMetrics::default(); hi.libs.len()];
            let mut ts = bamwrite::TagScratch::default();
            let mut rec = Vec::new();
            let mut batch: Vec<u8> = Vec::with_capacity(32 << 20);
            let mut it = dups.iter().copied();
            let mut next_dup = it.next().unwrap_or(u64::MAX); // NO_SUCH_INDEX = Long.MAX_VALUE
            let mut dup_name: Option<Vec<u8>> = None;
            let mut i = 0u64;
            while src.next(&mut rec) {
                let same_as_dup = dup_name.as_deref() == Some(rec_name(&rec));
                // nextIndexIfNeeded (queryname)
                if i > next_dup && !same_as_dup {
                    next_dup = it.next().unwrap_or(u64::MAX);
                }
                let is_dup = i == next_dup || (i > next_dup && same_as_dup);
                let mut flag = rec_flag(&rec);
                flag = if is_dup { flag | 0x400 } else { flag & !0x400 };
                rec[14..16].copy_from_slice(&flag.to_le_bytes());
                let aux = bamwrite::aux_start(&rec);
                let (lib, _) = hi.lib_and_rg(find_z(&rec, aux, b"RG"));
                lm[lib as usize].add_read(flag, is_dup);
                if is_dup {
                    dup_name = Some(rec_name(&rec).to_vec());
                }
                retag(&mut rec, &pg_map, args, &mut ts);
                batch.extend_from_slice(&(rec.len() as u32).to_le_bytes());
                batch.extend_from_slice(&rec);
                if batch.len() >= 32 << 20 {
                    tx.send(std::mem::replace(&mut batch, Vec::with_capacity(32 << 20))).unwrap();
                }
                i += 1;
            }
            if !batch.is_empty() {
                tx.send(batch).unwrap();
            }
            lm
        });
        let compress = |b: &Vec<u8>| -> bamwrite::ChunkOut {
            let mut w = bamwrite::RecWriter::new(args.level, &big);
            let mut o = 0;
            while o < b.len() {
                let l = u32::from_le_bytes(b[o..o + 4].try_into().unwrap()) as usize;
                w.push(&b[o + 4..o + 4 + l]);
                o += 4 + l;
            }
            w.finish()
        };
        let mut pending: Vec<Vec<u8>> = Vec::new();
        loop {
            let got = rx.recv();
            let done = got.is_err();
            if let Ok(b) = got {
                pending.push(b);
            }
            if pending.len() >= n_threads || (done && !pending.is_empty()) {
                let outs: Vec<bamwrite::ChunkOut> = pending.par_iter().map(compress).collect();
                for o in &outs {
                    bam.append(-1, o);
                }
                pending.clear();
            }
            if done {
                break;
            }
        }
        producer.join().unwrap()
    });
    let n_written = bam.n_records;
    bam.finish(&args.out);
    eprintln!("pass 2: {n_written} records written in {:.1}s", t2.elapsed().as_secs_f64());
    write_metrics(args, &cmdline, &header_text, &hi, &lm, &res);
    eprintln!("total {:.1}s", t0.elapsed().as_secs_f64());
}
