//! MarkDuplicates (`variant-rs markdup`): Picard 3.5.0 `MarkDuplicates` on a coordinate-sorted BAM
//! (SUM_OF_BASE_QUALITIES scoring, optical duplicates from the default read name parsing, no tagging), with
//! OPTICAL_DUPLICATE_PIXEL_DISTANCE, CLEAR_DT, ADD_PG_TAG_TO_READS and CREATE_MD5_FILE; several coordinate-sorted
//! inputs merged per chunk ([`ChunkRecords`]); query-grouped input (one or several) in `qname.rs`.
//!
//! Pass 1 (parallel over 2 Mb chunks in file order): one fragment end per primary mapped read and,
//! for pairs with both mates mapped, a pair end once both mates are seen (joined by read group +
//! name within a chunk, then across chunks). Ends go to on-disk buckets by (contig, 5' position).
//! Grouping (parallel over buckets): each bucket is sorted with Picard's comparator and duplicate
//! sets are marked exactly as Picard does; the duplicate record indexes go into a bitset.
//! Pass 2 (parallel over chunks): flags set/cleared, PG chained, tags re-encoded the htsjdk way,
//! BGZF + BAI written in order.

pub mod metrics;
pub(crate) mod qname;

use crate::bamwrite::{self, Chunk};
use crate::input::Input;
use metrics::LibMetrics;
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap};
use std::io::{BufWriter, Read as _, Write};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;

struct Args {
    /// the first input
    bam: String,
    /// all inputs (`INPUT` repeated: merged as htsjdk's `SamFileHeaderMerger` + `MergingSamRecordIterator`)
    inputs: Vec<String>,
    /// `ASSUME_SORT_ORDER` (queryname: query-grouped input, output in input order, `SO:unknown GO:query`)
    assume: Option<String>,
    optical_distance: i32,
    clear_dt: bool,
    add_pg: bool,
    /// `CREATE_INDEX`: write `<out>.bai` (coordinate-sorted output only, as htsjdk)
    index: bool,
    md5: bool,
    out: String,
    metrics: String,
    fasta: String,
    threads: usize,
    chunk: i64,
    level: i32,
    tmp: String,
    bucket_bp: i64,
}

const USAGE: &str = "usage: variant-rs markdup -I in.bam -O out.bam -M metrics.txt [--threads N (16)] [--tmp-dir DIR (.)] \
[--compression-level N (5)] [-R ref.fa (CRAM input)] [-I in2.bam ...] [--ASSUME_SORT_ORDER queryname] \
[--OPTICAL_DUPLICATE_PIXEL_DISTANCE N (100)] [--CLEAR_DT true|false] [--ADD_PG_TAG_TO_READS true|false] [--CREATE_MD5_FILE true|false]";

fn parse_args(a: &[String]) -> Args {
    let mut r = Args { bam: String::new(), inputs: vec![], assume: None, optical_distance: 100, clear_dt: true, add_pg: true, md5: false, index: false, out: String::new(), metrics: String::new(), fasta: String::new(), threads: 16, chunk: 2_000_000, level: 5, tmp: ".".into(), bucket_bp: 20_000_000 };
    let mut i = 0;
    while i < a.len() {
        let v = || a.get(i + 1).cloned().unwrap_or_else(|| panic!("{} needs a value\n{USAGE}", a[i]));
        match a[i].as_str() {
            "-I" | "--INPUT" | "--input" => r.inputs.push(v()),
            "--ASSUME_SORT_ORDER" | "-ASO" => r.assume = Some(v()),
            "--OPTICAL_DUPLICATE_PIXEL_DISTANCE" => r.optical_distance = v().parse().expect("--OPTICAL_DUPLICATE_PIXEL_DISTANCE"),
            "--CLEAR_DT" => r.clear_dt = bool_arg(&v()),
            "--ADD_PG_TAG_TO_READS" => r.add_pg = bool_arg(&v()),
            "--CREATE_MD5_FILE" => r.md5 = bool_arg(&v()),
            "--CREATE_INDEX" => r.index = bool_arg(&v()),
            "-O" | "--OUTPUT" | "--output" => r.out = v(),
            "-M" | "--METRICS_FILE" | "--metrics" => r.metrics = v(),
            "-R" | "--reference" => r.fasta = v(),
            "--threads" => r.threads = v().parse().expect("--threads"),
            "--chunk" => r.chunk = v().parse().expect("--chunk"),
            "--compression-level" => r.level = v().parse().expect("--compression-level"),
            "--tmp-dir" | "--TMP_DIR" => r.tmp = v(),
            "--bucket-bp" => r.bucket_bp = v().parse().expect("--bucket-bp"),
            x => panic!("unknown argument {x}\n{USAGE}"),
        }
        i += 2;
    }
    assert!(!r.inputs.is_empty() && !r.out.is_empty() && !r.metrics.is_empty(), "{USAGE}");
    r.bam = r.inputs[0].clone();
    r
}

fn bool_arg(v: &str) -> bool {
    match v.to_ascii_lowercase().as_str() {
        "true" | "t" | "yes" | "y" => true,
        "false" | "f" | "no" | "n" => false,
        _ => panic!("not a boolean: {v}"),
    }
}

// ------------------------------------------------------------------------------------------------
// read ends
// ------------------------------------------------------------------------------------------------

const F: u8 = 0;
const R: u8 = 1;
const FF: u8 = 2;
const FR: u8 = 3;
const RR: u8 = 4;
const RF: u8 = 5;

#[inline]
fn orientation_byte(r1_neg: bool, r2_neg: bool) -> u8 {
    match (r1_neg, r2_neg) {
        (true, true) => RR,
        (true, false) => RF,
        (false, true) => FR,
        (false, false) => FF,
    }
}

/// `ReadEndsForMarkDuplicates` (fragment or pair). `idx1/idx2` are packed (chunk << 32 | local),
/// which orders like the index in the file.
#[derive(Clone, Copy, Debug, Default)]
#[repr(C)]
struct End {
    lib: u16,
    orientation: u8,
    ofo: u8, // orientationForOpticalDuplicates (u8::MAX = -1)
    rg: i16,
    tile: i16,
    ref1: i32,
    coord1: i32,
    ref2: i32,
    coord2: i32,
    x: i32,
    y: i32,
    score: i16,
    _pad: [u8; 6],
    idx1: u64,
    idx2: u64, // u64::MAX = -1
}

const END_SIZE: usize = std::mem::size_of::<End>();

impl End {
    fn is_paired(&self) -> bool {
        self.ref2 != -1
    }
    fn has_location(&self) -> bool {
        self.tile != -1
    }
    /// `ReadEndsMDComparator`
    fn key(&self) -> (u16, i32, i32, u8, i32, i32, i16, i32, i32, u64, i64) {
        let i2 = if self.idx2 == u64::MAX { -1 } else { self.idx2 as i64 };
        (self.lib, self.ref1, self.coord1, self.orientation, self.ref2, self.coord2, self.tile, self.x, self.y, self.idx1, i2)
    }
}

fn ends_to_bytes(v: &[End]) -> &[u8] {
    unsafe { std::slice::from_raw_parts(v.as_ptr() as *const u8, v.len() * END_SIZE) }
}

fn bytes_to_ends(b: &[u8]) -> Vec<End> {
    assert_eq!(b.len() % END_SIZE, 0);
    let n = b.len() / END_SIZE;
    let mut v: Vec<End> = Vec::with_capacity(n);
    unsafe {
        std::ptr::copy_nonoverlapping(b.as_ptr(), v.as_mut_ptr() as *mut u8, b.len());
        v.set_len(n);
    }
    v
}

/// `ReadNameParser.getLastThreeFields` + `rapidParseInt` (default READ_NAME_REGEX): (tile, x, y).
fn parse_location(name: &[u8]) -> Option<(i16, i32, i32)> {
    fn rapid(s: &[u8]) -> Option<i32> {
        let mut i = 0;
        let mut neg = false;
        if !s.is_empty() && s[0] == b'-' {
            i = 1;
            neg = true;
        }
        let mut val: i32 = 0;
        let mut has = false;
        while i < s.len() && s[i].is_ascii_digit() {
            val = val.wrapping_mul(10).wrapping_add((s[i] - b'0') as i32);
            has = true;
            i += 1;
        }
        if !has {
            return None;
        }
        Some(if neg { val.wrapping_neg() } else { val })
    }
    let mut tokens = [0i32; 3];
    let mut ti: i32 = 2;
    let mut num_fields = 0;
    let mut end = name.len();
    let mut i: i64 = name.len() as i64 - 1;
    while i >= 0 && ti >= 0 {
        let iu = i as usize;
        if name[iu] == b':' || iu == 0 {
            num_fields += 1;
            let from = if iu == 0 { 0 } else { iu + 1 };
            tokens[ti as usize] = rapid(&name[from..end])?;
            ti -= 1;
            end = iu;
        }
        i -= 1;
    }
    while i >= 0 {
        let iu = i as usize;
        if name[iu] == b':' || iu == 0 {
            num_fields += 1;
        }
        i -= 1;
    }
    if num_fields < 3 || !(num_fields == 5 || num_fields == 7) {
        return None;
    }
    Some((tokens[0] as i16, tokens[1], tokens[2]))
}

// ------------------------------------------------------------------------------------------------
// header information
// ------------------------------------------------------------------------------------------------

struct HeaderInfo {
    /// RG ID -> (header ordinal, library index)
    rg: HashMap<Vec<u8>, (i16, u16)>,
    libs: Vec<String>, // library index -> name ("Unknown Library" last)
    unknown_lib: u16,
}

impl HeaderInfo {
    fn new(text: &str) -> HeaderInfo {
        let mut rgs: Vec<(String, String)> = Vec::new();
        for l in text.lines().filter(|l| l.starts_with("@RG\t")) {
            let mut id = String::new();
            let mut lb = None;
            for f in l.split('\t').skip(1) {
                if let Some(v) = f.strip_prefix("ID:") {
                    id = v.to_string();
                } else if let Some(v) = f.strip_prefix("LB:") {
                    lb = Some(v.to_string());
                }
            }
            rgs.push((id, lb.unwrap_or_else(|| "Unknown Library".into())));
        }
        let mut libs: Vec<String> = rgs.iter().map(|x| x.1.clone()).collect::<std::collections::BTreeSet<_>>().into_iter().collect();
        let unknown_lib = match libs.iter().position(|l| l == "Unknown Library") {
            Some(p) => p as u16,
            None => {
                libs.push("Unknown Library".into());
                (libs.len() - 1) as u16
            }
        };
        let mut rg = HashMap::new();
        for (i, (id, lb)) in rgs.iter().enumerate() {
            let li = libs.iter().position(|l| l == lb).unwrap() as u16;
            rg.entry(id.clone().into_bytes()).or_insert((i as i16, li));
        }
        HeaderInfo { rg, libs, unknown_lib }
    }
    /// (library index, read group ordinal as `buildReadEnds` computes it)
    fn lib_and_rg(&self, rg: Option<&[u8]>) -> (u16, i16) {
        match rg {
            Some(id) => match self.rg.get(id) {
                Some(&(o, l)) => (l, o),
                None => (self.unknown_lib, self.rg.len() as i16),
            },
            None => (self.unknown_lib, 0),
        }
    }
}

// ------------------------------------------------------------------------------------------------
// record access
// ------------------------------------------------------------------------------------------------

/// Value of string tag `tag` in record `rec` (aux data from `aux`).
fn find_z<'a>(rec: &'a [u8], aux: usize, tag: &[u8; 2]) -> Option<&'a [u8]> {
    let mut o = aux;
    while o + 3 <= rec.len() {
        let t = [rec[o], rec[o + 1]];
        let ty = rec[o + 2];
        o += 3;
        match ty {
            b'Z' | b'H' => {
                let s = o;
                while rec[o] != 0 {
                    o += 1;
                }
                if t == *tag && ty == b'Z' {
                    return Some(&rec[s..o]);
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

struct RecView {
    ref_id: i32,
    flag: u16,
    mate_ref: i32,
    unclipped_5p: i64,
    score: i16,
}

fn view(rec: &[u8]) -> RecView {
    let ref_id = i32::from_le_bytes(rec[0..4].try_into().unwrap());
    let start = i32::from_le_bytes(rec[4..8].try_into().unwrap()) as i64 + 1;
    let l_read_name = rec[8] as usize;
    let n_cigar = u16::from_le_bytes([rec[12], rec[13]]) as usize;
    let flag = u16::from_le_bytes([rec[14], rec[15]]);
    let l_seq = i32::from_le_bytes(rec[16..20].try_into().unwrap()) as usize;
    let mate_ref = i32::from_le_bytes(rec[20..24].try_into().unwrap());
    let co = 32 + l_read_name;
    let op = |k: usize| -> (u32, u32) {
        let v = u32::from_le_bytes(rec[co + 4 * k..co + 4 * k + 4].try_into().unwrap());
        (v & 0xf, v >> 4)
    };
    let neg = flag & 0x10 != 0;
    let unclipped_5p = if !neg {
        let mut s = start;
        for k in 0..n_cigar {
            let (o, l) = op(k);
            if o == 4 || o == 5 {
                s -= l as i64;
            } else {
                break;
            }
        }
        s
    } else {
        let mut rl = 0i64;
        for k in 0..n_cigar {
            let (o, l) = op(k);
            if matches!(o, 0 | 2 | 3 | 7 | 8) {
                rl += l as i64;
            }
        }
        let mut e = start + rl - 1;
        for k in (0..n_cigar).rev() {
            let (o, l) = op(k);
            if o == 4 || o == 5 {
                e += l as i64;
            } else {
                break;
            }
        }
        e
    };
    // DuplicateScoringStrategy SUM_OF_BASE_QUALITIES
    let qo = co + 4 * n_cigar + (l_seq + 1) / 2;
    let mut sum: i32 = 0;
    if l_seq > 0 && rec[qo] != 0xff {
        // getBaseQualities() is a Java byte[]: 128..254 are negative and not counted
        for &q in &rec[qo..qo + l_seq] {
            if q as i8 >= 15 {
                sum += q as i32;
            }
        }
    }
    let mut score = sum.min(i16::MAX as i32 / 2) as i16;
    if flag & 0x200 != 0 {
        score = score.wrapping_add(i16::MIN / 2);
    }
    RecView { ref_id, flag, mate_ref, unclipped_5p, score }
}

/// One mate waiting for the other: its fragment end plus what the pairing needs.
#[derive(Clone)]
struct Half {
    end: End,
    neg: bool,
    first_of_pair: bool,
    target_ref: i32, // the reference index it is filed under (its mate's)
}

/// The pairing in `buildSortedReadEndLists`: `a` was seen first in the file, `b` completes it.
fn make_pair(a: &Half, b: &Half) -> End {
    let mut p = a.end;
    let a_neg = a.end.orientation == R;
    p.ofo = if b.first_of_pair { orientation_byte(b.neg, a_neg) } else { orientation_byte(a_neg, b.neg) };
    let (mref, mcoord) = (b.end.ref1, b.end.coord1);
    if mref > p.ref1 || (mref == p.ref1 && mcoord >= p.coord1) {
        p.ref2 = mref;
        p.coord2 = mcoord;
        p.idx2 = b.end.idx1;
        p.orientation = orientation_byte(a_neg, b.neg);
        if p.ref2 == p.ref1 && p.coord2 == p.coord1 && p.orientation == RF {
            p.orientation = FR;
        }
    } else {
        p.ref2 = p.ref1;
        p.coord2 = p.coord1;
        p.idx2 = p.idx1;
        p.ref1 = mref;
        p.coord1 = mcoord;
        p.idx1 = b.end.idx1;
        p.orientation = orientation_byte(b.neg, a_neg);
    }
    p.score = p.score.wrapping_add(b.end.score);
    p
}

// ------------------------------------------------------------------------------------------------
// buckets on disk
// ------------------------------------------------------------------------------------------------

struct Buckets {
    files: Vec<Mutex<BufWriter<std::fs::File>>>,
    paths: Vec<String>,
}

impl Buckets {
    fn new(dir: &str, prefix: &str, n: usize) -> Buckets {
        let mut files = Vec::with_capacity(n);
        let mut paths = Vec::with_capacity(n);
        for i in 0..n {
            let p = format!("{dir}/{prefix}.{i:05}");
            files.push(Mutex::new(BufWriter::with_capacity(1 << 20, std::fs::File::create(&p).unwrap_or_else(|e| panic!("{p}: {e}")))));
            paths.push(p);
        }
        Buckets { files, paths }
    }
    fn flush_local(&self, local: &mut HashMap<usize, Vec<End>>, force: bool) {
        for (b, v) in local.iter_mut() {
            if !v.is_empty() && (force || v.len() >= 8192) {
                self.files[*b].lock().unwrap().write_all(ends_to_bytes(v)).unwrap();
                v.clear();
            }
        }
    }
    fn finish(&self) {
        for f in &self.files {
            f.lock().unwrap().flush().unwrap();
        }
    }
    fn size(&self, b: usize) -> u64 {
        std::fs::metadata(&self.paths[b]).map_or(0, |m| m.len())
    }
    fn read(&self, b: usize) -> Vec<End> {
        let mut data = Vec::new();
        std::fs::File::open(&self.paths[b]).unwrap().read_to_end(&mut data).unwrap();
        bytes_to_ends(&data)
    }
}

impl Drop for Buckets {
    fn drop(&mut self) {
        for p in &self.paths {
            let _ = std::fs::remove_file(p);
        }
    }
}

struct BucketMap {
    cum: Vec<i64>,
    lens: Vec<i64>,
    bp: i64,
    n: usize,
}

impl BucketMap {
    fn new(lens: &[i64], bp: i64) -> BucketMap {
        let mut cum = Vec::with_capacity(lens.len());
        let mut c = 0i64;
        for &l in lens {
            cum.push(c);
            c += l.max(1);
        }
        BucketMap { cum, lens: lens.to_vec(), bp, n: (c / bp) as usize + 1 }
    }
    fn bucket(&self, r: i32, coord: i32) -> usize {
        let r = r as usize;
        let c = (coord as i64).clamp(1, self.lens[r].max(1));
        ((self.cum[r] + c - 1) / self.bp) as usize
    }
}

// ------------------------------------------------------------------------------------------------
// optical duplicates (OpticalDuplicateFinder)
// ------------------------------------------------------------------------------------------------

const MAX_OPTICAL_SET: usize = 300000;

fn close_enough(a: &End, b: &End, dist: i32) -> bool {
    a.has_location() && b.has_location() && a.rg == b.rg && a.tile == b.tile && (a.x - b.x).abs() <= dist && (a.y - b.y).abs() <= dist
}

/// Iteration order of a Java `HashMap<Integer, _>` holding the keys 0..n-1 (inserted ascending).
fn java_int_hashmap_order(n: usize) -> Vec<usize> {
    let mut cap = 16usize;
    while n as f64 > cap as f64 * 0.75 {
        cap *= 2;
    }
    let mut v: Vec<(usize, usize)> = (0..n)
        .map(|k| {
            let h = k as u32;
            let h = h ^ (h >> 16);
            ((h as usize) & (cap - 1), k)
        })
        .collect();
    v.sort();
    v.into_iter().map(|x| x.1).collect()
}

/// `findOpticalDuplicates(list, keeper)`; `keeper` = index of the keeper in `list` if it is there
/// (`keeper_given`: the caller passed a non-null keeper).
/// `dist`: OPTICAL_DUPLICATE_PIXEL_DISTANCE.
fn find_optical(list: &[End], keeper: Option<usize>, keeper_given: bool, dist: i32) -> Vec<bool> {
    let n = list.len();
    let mut flags = vec![false; n];
    if n < 2 || n > MAX_OPTICAL_SET {
        return flags;
    }
    let actual = keeper.filter(|&k| list[k].has_location());
    if n >= if keeper_given { 4 } else { 3 } {
        // graph: nodes 0..n-1, edges within (read group, tile) groups
        let mut parent: Vec<usize> = (0..n).collect();
        fn find(p: &mut [usize], mut x: usize) -> usize {
            let mut r = x;
            while p[r] != r {
                r = p[r];
            }
            while x != r {
                let nx = p[x];
                p[x] = r;
                x = nx;
            }
            r
        }
        let mut groups: HashMap<i32, Vec<usize>> = HashMap::new();
        for (i, e) in list.iter().enumerate() {
            if e.has_location() {
                groups.entry(((e.rg as i32) << 16) + e.tile as i32).or_default().push(i);
            }
        }
        for g in groups.values() {
            for a in 0..g.len() {
                for b in a + 1..g.len() {
                    let (i, j) = (g[a], g[b]);
                    if (list[i].x - list[j].x).abs() <= dist && (list[i].y - list[j].y).abs() <= dist {
                        let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                        if ri != rj {
                            parent[rj] = ri;
                        }
                    }
                }
            }
        }
        let cluster: Vec<usize> = (0..n).map(|i| find(&mut parent, i)).collect();
        let mut rep: HashMap<usize, usize> = HashMap::new();
        let keeper_cluster = actual.map(|k| {
            rep.insert(cluster[k], k);
            cluster[k]
        });
        for i in java_int_hashmap_order(n) {
            let c = cluster[i];
            if rep.contains_key(&c) && Some(i) != actual {
                let r = rep[&c];
                let (cur, rl) = (&list[i], &list[r]);
                if !(actual.is_some() && Some(c) == keeper_cluster) && (cur.x < rl.x || (cur.x == rl.x && cur.y < rl.y)) {
                    flags[r] = true;
                    rep.insert(c, i);
                } else {
                    flags[i] = true;
                }
            } else {
                rep.insert(c, i);
            }
        }
    } else {
        if let Some(k) = actual {
            for i in 0..n {
                flags[i] = i != k && close_enough(&list[k], &list[i], dist);
            }
        }
        for i in 0..n {
            if Some(i) == actual {
                continue;
            }
            for j in i + 1..n {
                if Some(j) == actual {
                    continue;
                }
                if flags[i] && flags[j] {
                    continue;
                }
                if close_enough(&list[i], &list[j], dist) {
                    let idx = if flags[j] { i } else { j };
                    flags[idx] = true;
                }
            }
        }
    }
    flags
}

#[derive(Default)]
struct GroupResult {
    dups: Vec<u64>,
    n_dups: u64,
    optical_by_lib: HashMap<u16, i64>,
    all_sets: HashMap<i64, i64>,
    optical_sets: HashMap<i64, i64>,
    non_optical_sets: HashMap<i64, i64>,
}

impl GroupResult {
    fn merge(&mut self, o: GroupResult) {
        self.dups.extend(o.dups);
        self.n_dups += o.n_dups;
        for (k, v) in o.optical_by_lib {
            *self.optical_by_lib.entry(k).or_default() += v;
        }
        for (a, b) in [(&mut self.all_sets, o.all_sets), (&mut self.optical_sets, o.optical_sets), (&mut self.non_optical_sets, o.non_optical_sets)] {
            for (k, v) in b {
                *a.entry(k).or_default() += v;
            }
        }
    }
}

fn best_index(list: &[End]) -> usize {
    let mut best = 0;
    let mut max = 0i16;
    for (i, e) in list.iter().enumerate() {
        if e.score > max || i == 0 {
            max = e.score;
            best = i;
        }
    }
    best
}

/// `markDuplicatePairs` (+ optical tracking and set-size histograms) for one duplicate set.
fn mark_pairs(list: &[End], res: &mut GroupResult, dist: i32) {
    let best = best_index(list);
    // trackOpticalDuplicates
    let has_fr = list.iter().any(|e| e.ofo == FR);
    let has_rf = list.iter().any(|e| e.ofo == RF);
    let mut n_optical = 0i64;
    let count = |sub: &[usize], res: &mut GroupResult| -> i64 {
        let l: Vec<End> = sub.iter().map(|&i| list[i]).collect();
        let k = sub.iter().position(|&i| i == best);
        let flags = find_optical(&l, k, true, dist);
        let n = flags.iter().filter(|&&f| f).count() as i64;
        if n > 0 {
            *res.optical_by_lib.entry(l[0].lib).or_default() += n;
        }
        n
    };
    if has_fr && has_rf {
        let fr: Vec<usize> = (0..list.len()).filter(|&i| list[i].ofo == FR).collect();
        let rf: Vec<usize> = (0..list.len()).filter(|&i| list[i].ofo == RF).collect();
        assert_eq!(fr.len() + rf.len(), list.len(), "Found an unexpected orientation");
        n_optical += count(&fr, res);
        n_optical += count(&rf, res);
    } else {
        let all: Vec<usize> = (0..list.len()).collect();
        n_optical += count(&all, res);
    }
    let size = list.len() as i64;
    *res.all_sets.entry(size).or_default() += 1;
    if size - n_optical > 0 {
        *res.non_optical_sets.entry(size - n_optical).or_default() += 1;
    }
    if n_optical > 0 {
        *res.optical_sets.entry(n_optical + 1).or_default() += 1;
    }
    for (i, e) in list.iter().enumerate() {
        if i != best {
            res.dups.push(e.idx1);
            if e.idx2 != e.idx1 {
                res.dups.push(e.idx2);
            }
        }
    }
}

/// `markDuplicateFragments`.
fn mark_frags(list: &[End], res: &mut GroupResult) {
    if list.iter().any(|e| e.is_paired()) {
        for e in list {
            if !e.is_paired() {
                res.dups.push(e.idx1);
            }
        }
    } else {
        let best = best_index(list);
        for (i, e) in list.iter().enumerate() {
            if i != best {
                res.dups.push(e.idx1);
            }
        }
    }
}

fn group_bucket(frags: Vec<End>, pairs: Vec<End>, dist: i32) -> GroupResult {
    let mut res = GroupResult::default();
    let mut pairs = pairs;
    pairs.sort_unstable_by_key(|e| e.key());
    let mut i = 0;
    while i < pairs.len() {
        let a = &pairs[i];
        let mut j = i + 1;
        while j < pairs.len() {
            let b = &pairs[j];
            if !(a.lib == b.lib && a.ref1 == b.ref1 && a.coord1 == b.coord1 && a.orientation == b.orientation && a.ref2 == b.ref2 && a.coord2 == b.coord2) {
                break;
            }
            j += 1;
        }
        if j - i > 1 {
            mark_pairs(&pairs[i..j], &mut res, dist);
        } else {
            *res.all_sets.entry(1).or_default() += 1;
            *res.non_optical_sets.entry(1).or_default() += 1;
        }
        i = j;
    }
    let mut frags = frags;
    frags.sort_unstable_by_key(|e| e.key());
    let mut i = 0;
    while i < frags.len() {
        let a = &frags[i];
        let mut j = i + 1;
        while j < frags.len() {
            let b = &frags[j];
            if !(a.lib == b.lib && a.ref1 == b.ref1 && a.coord1 == b.coord1 && a.orientation == b.orientation) {
                break;
            }
            j += 1;
        }
        if j - i > 1 {
            mark_frags(&frags[i..j], &mut res);
        }
        i = j;
    }
    res
}

// ------------------------------------------------------------------------------------------------
// PG chaining
// ------------------------------------------------------------------------------------------------

fn java_string_hash(s: &[u8]) -> i32 {
    let mut h: i32 = 0;
    for &c in s {
        h = h.wrapping_mul(31).wrapping_add(c as i32);
    }
    h
}

/// `getChainedPgIds`: PG values seen (first-occurrence order) -> new program record IDs, in the
/// iteration order of Picard's `HashSet<String>`, plus the new @PG lines.
fn chain_pg_ids(seen: &[Option<Vec<u8>>], header_text: &str, cmdline: &str) -> (HashMap<Option<Vec<u8>>, String>, Vec<String>) {
    let mut cap = 16usize;
    while seen.len() as f64 > cap as f64 * 0.75 {
        cap *= 2;
    }
    let mut order: Vec<(usize, usize)> = seen
        .iter()
        .enumerate()
        .map(|(i, s)| {
            let h = match s {
                None => 0u32,
                Some(v) => java_string_hash(v) as u32,
            };
            (((h ^ (h >> 16)) as usize) & (cap - 1), i)
        })
        .collect();
    order.sort();
    let mut taken: std::collections::HashSet<String> = header_text
        .lines()
        .filter(|l| l.starts_with("@PG\t"))
        .filter_map(|l| l.split('\t').find_map(|f| f.strip_prefix("ID:")).map(String::from))
        .collect();
    let mut counter = taken.len() as i32;
    let base36 = |mut v: i32| -> String {
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
    };
    let mut map = HashMap::new();
    let mut lines = Vec::new();
    for (_, i) in order {
        let id = if !taken.contains("MarkDuplicates") {
            counter += 1;
            "MarkDuplicates".to_string()
        } else {
            loop {
                let c = format!("MarkDuplicates.{}", base36(counter));
                counter += 1;
                if !taken.contains(&c) {
                    break c;
                }
            }
        };
        taken.insert(id.clone());
        let mut line = format!("@PG\tID:{id}\tVN:variant-rs {}\tCL:{cmdline}\tPN:MarkDuplicates", env!("CARGO_PKG_VERSION"));
        if let Some(pp) = &seen[i] {
            line.push_str(&format!("\tPP:{}", String::from_utf8_lossy(pp)));
        }
        lines.push(line);
        map.insert(seen[i].clone(), id);
    }
    (map, lines)
}

// ------------------------------------------------------------------------------------------------
// driver
// ------------------------------------------------------------------------------------------------

struct Pass1 {
    n_records: u64,
    pg_seen: Vec<Option<Vec<u8>>>,
    cross: Vec<(Vec<u8>, Half)>,
}

/// The `SO` value of a header's `@HD` line (htsjdk: `unsorted` when missing).
fn sort_order(text: &str) -> String {
    text.lines()
        .find(|l| l.starts_with("@HD"))
        .and_then(|l| l.split('\t').find_map(|f| f.strip_prefix("SO:")))
        .unwrap_or("unsorted")
        .to_string()
}

/// Pass 2's tag handling for one record: `CLEAR_DT` removes DT; `ADD_PG_TAG_TO_READS` sets PG to the chained
/// MarkDuplicates program ID (and removes it for a PG value pass 1 never saw, as `chainedPgIds.get` returns null).
/// Either makes htsjdk re-encode all tags (`normalize_tags`); a record whose attributes were not set keeps its
/// bytes (BAMRecord.getVariableBinaryRepresentation).
fn retag(rec: &mut Vec<u8>, pg_map: &HashMap<Option<Vec<u8>>, String>, args: &Args, ts: &mut bamwrite::TagScratch) {
    // with neither, Picard sets no attribute and htsjdk writes the record's tag bytes as read
    if !args.clear_dt && !args.add_pg {
        return;
    }
    let aux = bamwrite::aux_start(rec);
    let dt: &[[u8; 2]] = if args.clear_dt { &[*b"DT"] } else { &[] };
    if !args.add_pg {
        bamwrite::normalize_tags(rec, aux, dt, None, ts);
        return;
    }
    let pg = find_z(rec, aux, b"PG").map(|v| v.to_vec());
    match pg_map.get(&pg) {
        Some(id) => bamwrite::normalize_tags(rec, aux, dt, Some((*b"PG", id.as_bytes())), ts),
        None => {
            let mut rm = dt.to_vec();
            rm.push(*b"PG");
            bamwrite::normalize_tags(rec, aux, &rm, None, ts)
        }
    }
}

pub fn run(argv: &[String]) {
    let args = parse_args(argv);
    rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();
    // Picard's sort order: ASSUME_SORT_ORDER, else the (first) input header's
    let first_so = {
        let f = std::fs::File::open(&args.bam).unwrap_or_else(|e| panic!("{}: {e}", args.bam));
        let mut magic = [0u8; 4];
        std::os::unix::fs::FileExt::read_exact_at(&f, &mut magic, 0).unwrap_or_else(|e| panic!("{}: {e}", args.bam));
        if &magic == b"CRAM" { "coordinate".to_string() } else { sort_order(&crate::bam::Header::read(&f).text) }
    };
    let so = args.assume.clone().unwrap_or(first_so.clone());
    match so.as_str() {
        "queryname" => return qname::run(&args),
        "coordinate" => {
            if first_so != "coordinate" {
                crate::gatk_cli::refuse("ASSUME_SORT_ORDER=coordinate on an input that is not coordinate-sorted is not supported by variant-rs");
            }
        }
        _ => panic!(
            "This program requires input that are either coordinate or query sorted (according to the header, or at least ASSUME_SORT_ORDER and the content.) Found ASSUME_SORT_ORDER={} and header sortorder={first_so}",
            args.assume.as_deref().unwrap_or("null")
        ),
    }
    let t0 = Instant::now();
    // several inputs: SamFileHeaderMerger + MergingSamRecordIterator (SAMRecordCoordinateComparator), per chunk
    let (inputs, headers): (Vec<Input>, Vec<crate::bam::Header>) = args.inputs.iter().map(|p| Input::open(p, &args.fasta)).unzip();
    for (p, h) in args.inputs.iter().zip(&headers) {
        if !h.text.lines().next().is_some_and(|l| l.contains("SO:coordinate")) {
            crate::gatk_cli::refuse(&format!("{p}: not coordinate-sorted (@HD SO:coordinate), while the first input is"));
        }
    }
    let merged: Option<qname::Merged> = (inputs.len() > 1).then(|| qname::merge_header_list(&headers.iter().collect::<Vec<_>>(), "coordinate"));
    let header = crate::bam::Header {
        text: merged.as_ref().map_or_else(|| headers[0].text.clone(), |m| m.text.clone()),
        ref_names: headers[0].ref_names.clone(),
        ref_lens: headers[0].ref_lens.clone(),
    };
    let hi = HeaderInfo::new(&header.text);
    let traversal: Vec<(usize, i64, i64)> = header.ref_lens.iter().enumerate().map(|(t, &l)| (t, 1, l)).collect();
    let chunks = bamwrite::plan_chunks(&traversal, header.ref_lens.len(), args.chunk, true);
    let chunks = split_dense(chunks, &inputs, &header.ref_lens, args.chunk);
    eprintln!("{} chunks", chunks.len());
    let tails: Vec<(u64, u64)> = inputs.iter().zip(&args.inputs).map(|(i, p)| bamwrite::tail_range(i, p)).collect();
    let io_ns = Arc::new(AtomicU64::new(0));
    let open = |ch: &Chunk| ChunkRecords::open(&inputs, &tails, merged.as_ref(), ch, &io_ns);
    let bmap = BucketMap::new(&header.ref_lens, args.bucket_bp);
    let pid = std::process::id();
    let frag_b = Buckets::new(&args.tmp, &format!("ardmd{pid}.frag"), bmap.n);
    let pair_b = Buckets::new(&args.tmp, &format!("ardmd{pid}.pair"), bmap.n);

    // ---- pass 1
    let mut pass1: Vec<Pass1> = chunks
        .par_iter()
        .enumerate()
        .map(|(ci, ch)| {
            let mut p = Pass1 { n_records: 0, pg_seen: Vec::new(), cross: Vec::new() };
            if ch.tid < 0 {
                // Picard stops at the first unplaced unmapped read, after adding its PG to pgIdsSeen (which matters
                // when no placed read came before it)
                if let Some(mut stream) = open(ch) {
                    let mut rec = Vec::new();
                    if stream.next(&mut rec) {
                        let aux = bamwrite::aux_start(&rec);
                        p.pg_seen.push(find_z(&rec, aux, b"PG").map(|v| v.to_vec()));
                    }
                }
                return p;
            }
            let Some(mut stream) = open(ch) else { return p };
            let mut rec = Vec::new();
            let mut lf: HashMap<usize, Vec<End>> = HashMap::new();
            let mut lp: HashMap<usize, Vec<End>> = HashMap::new();
            let mut waiting: HashMap<(Vec<u8>, i32), Half> = HashMap::new();
            let mut seen_pg: std::collections::HashSet<Option<Vec<u8>>> = std::collections::HashSet::new();
            let mut local: u64 = 0;
            while stream.next(&mut rec) {
                let idx = ((ci as u64) << 32) | local;
                local += 1;
                let aux = bamwrite::aux_start(&rec);
                let pg = find_z(&rec, aux, b"PG").map(|v| v.to_vec());
                if !seen_pg.contains(&pg) {
                    seen_pg.insert(pg.clone());
                    p.pg_seen.push(pg);
                }
                let v = view(&rec);
                if v.flag & 0x4 != 0 || v.flag & 0x900 != 0 {
                    continue;
                }
                let rg = find_z(&rec, aux, b"RG");
                let (lib, rg_ord) = hi.lib_and_rg(rg);
                let l_read_name = rec[8] as usize;
                let name = &rec[32..32 + l_read_name - 1];
                let mut e = End { lib, orientation: if v.flag & 0x10 != 0 { R } else { F }, ofo: u8::MAX, rg: -1, tile: -1, ref1: v.ref_id, coord1: v.unclipped_5p as i32, ref2: -1, coord2: -1, x: -1, y: -1, score: v.score, _pad: [0; 6], idx1: idx, idx2: u64::MAX };
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
                lf.entry(bmap.bucket(e.ref1, e.coord1)).or_default().push(e);
                if paired {
                    let mut key = rg.expect("paired read without a read group").to_vec();
                    key.extend_from_slice(name);
                    // tmp.remove(rec.getReferenceIndex(), key)
                    match waiting.remove(&(key.clone(), v.ref_id)) {
                        Some(first) => {
                            let second = Half { end: e, neg: v.flag & 0x10 != 0, first_of_pair: v.flag & 0x40 != 0, target_ref: v.mate_ref };
                            let pe = make_pair(&first, &second);
                            lp.entry(bmap.bucket(pe.ref1, pe.coord1)).or_default().push(pe);
                        }
                        None => {
                            let h = Half { end: e, neg: v.flag & 0x10 != 0, first_of_pair: v.flag & 0x40 != 0, target_ref: v.mate_ref };
                            waiting.insert((key, v.mate_ref), h);
                        }
                    }
                }
                if local % 4096 == 0 {
                    frag_b.flush_local(&mut lf, false);
                    pair_b.flush_local(&mut lp, false);
                }
            }
            frag_b.flush_local(&mut lf, true);
            pair_b.flush_local(&mut lp, true);
            p.n_records = local;
            p.cross = waiting.into_iter().map(|((k, _), h)| (k, h)).collect();
            p
        })
        .collect();
    let t_cross = Instant::now();
    // pairs whose mates are in different chunks, matched in file order. A name's halves only ever meet each other,
    // so the matching runs in parallel over partitions by name (each partition in file order): the same pairs.
    const PARTS: usize = 256;
    let part_of = |k: &[u8]| -> usize {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        for &b in k {
            h = (h ^ b as u64).wrapping_mul(0x100_0000_01b3);
        }
        (h % PARTS as u64) as usize
    };
    let n_cross: usize = pass1.iter().map(|p| p.cross.len()).sum();
    let mut parts: Vec<Vec<(Vec<u8>, Half)>> = (0..PARTS).map(|_| Vec::new()).collect();
    for p in pass1.iter_mut() {
        for (k, h) in std::mem::take(&mut p.cross) {
            let i = part_of(&k);
            parts[i].push((k, h));
        }
    }
    let n_unmatched: u64 = parts
        .into_par_iter()
        .map(|mut cross| {
            cross.sort_unstable_by_key(|x| x.1.end.idx1);
            let mut waiting: HashMap<(Vec<u8>, i32), Half> = HashMap::new();
            let mut lp: HashMap<usize, Vec<End>> = HashMap::new();
            for (key, h) in cross {
                match waiting.remove(&(key.clone(), h.end.ref1)) {
                    Some(first) => {
                        let pe = make_pair(&first, &h);
                        lp.entry(bmap.bucket(pe.ref1, pe.coord1)).or_default().push(pe);
                    }
                    None => {
                        let t = h.target_ref;
                        waiting.insert((key, t), h);
                    }
                }
            }
            pair_b.flush_local(&mut lp, true);
            waiting.len() as u64
        })
        .sum();
    frag_b.finish();
    pair_b.finish();
    let cross_secs = t_cross.elapsed().as_secs_f64();
    let offsets: Vec<u64> = {
        let mut o = Vec::with_capacity(chunks.len());
        let mut c = 0u64;
        for p in &pass1 {
            o.push(c);
            c += p.n_records;
        }
        o
    };
    let total: u64 = pass1.iter().map(|p| p.n_records).sum();
    eprintln!(
        "pass 1: {total} records, {} cross-chunk mates ({} without a mate) in {:.1}s (cross-chunk mates {:.1}s)",
        n_cross,
        n_unmatched,
        t0.elapsed().as_secs_f64(),
        cross_secs
    );

    // ---- grouping
    let t1 = Instant::now();
    let bits: Vec<AtomicU64> = (0..total / 64 + 1).map(|_| AtomicU64::new(0)).collect();
    // largest buckets first (their sorts are the longest tasks); the results are sums and a bitset, in any order
    let mut order: Vec<(u64, usize)> = (0..bmap.n).map(|b| (frag_b.size(b) + pair_b.size(b), b)).collect();
    order.sort_unstable_by(|a, b| b.cmp(a));
    let res = order
        .into_par_iter()
        .with_max_len(1)
        .map(|(_, b)| {
            let mut r = group_bucket(frag_b.read(b), pair_b.read(b), args.optical_distance);
            for &p in &r.dups {
                let g = offsets[(p >> 32) as usize] + (p & 0xffff_ffff);
                bits[(g / 64) as usize].fetch_or(1 << (g % 64), Ordering::Relaxed);
            }
            r.n_dups = r.dups.len() as u64;
            r.dups = Vec::new();
            r
        })
        .reduce(GroupResult::default, |mut a, b| {
            a.merge(b);
            a
        });
    drop(frag_b);
    drop(pair_b);
    eprintln!("duplicates: {} records in {:.1}s", res.n_dups, t1.elapsed().as_secs_f64());

    // ---- pass 2
    let t2 = Instant::now();
    let mut pg_seen: Vec<Option<Vec<u8>>> = Vec::new();
    {
        let mut s = std::collections::HashSet::new();
        for p in &pass1 {
            for v in &p.pg_seen {
                if s.insert(v.clone()) {
                    pg_seen.push(v.clone());
                }
            }
        }
    }
    let cmdline = std::env::args().collect::<Vec<_>>().join(" ");
    let (pg_map, pg_lines) = chain_pg_ids(&pg_seen, &header.text, &cmdline);
    let text = bamwrite::header_with_pg(&header.text, &pg_lines);
    let mut bam = bamwrite::BamFile::create(&args.out, &text, &header.ref_names, &header.ref_lens, args.level, args.md5);
    let big = bam.too_large.clone();
    bam.index = args.index; // CREATE_INDEX (Picard's default: false)
    let lib_metrics: Mutex<Vec<LibMetrics>> = Mutex::new(vec![LibMetrics::default(); hi.libs.len()]);
    let process = |(ci, ch): (usize, &Chunk)| -> bamwrite::ChunkOut {
        let mut w = bamwrite::RecWriter::new(args.level, &big);
        let Some(mut stream) = open(ch) else { return w.finish() };
        let mut rec = Vec::new();
        let mut ts = bamwrite::TagScratch::default();
        let mut lm = vec![LibMetrics::default(); hi.libs.len()];
        let mut local = 0u64;
        while stream.next(&mut rec) {
            let dup = if ch.tid >= 0 {
                let g = offsets[ci] + local;
                bits.get((g / 64) as usize).map_or(false, |w| w.load(Ordering::Relaxed) & (1 << (g % 64)) != 0)
            } else {
                false
            };
            local += 1;
            let mut flag = u16::from_le_bytes([rec[14], rec[15]]);
            flag = if dup { flag | 0x400 } else { flag & !0x400 };
            rec[14..16].copy_from_slice(&flag.to_le_bytes());
            let aux = bamwrite::aux_start(&rec);
            let rg = find_z(&rec, aux, b"RG");
            let (lib, _) = hi.lib_and_rg(rg);
            lm[lib as usize].add_read(flag, dup);
            retag(&mut rec, &pg_map, &args, &mut ts);
            w.push(&rec);
        }
        let mut g = lib_metrics.lock().unwrap();
        for (a, b) in g.iter_mut().zip(lm.iter()) {
            a.add(b);
        }
        w.finish()
    };
    // chunks run on the pool, at most `window` ahead of the writer; this thread appends them in order as they finish
    let window = 4 * args.threads.max(1);
    // a worker's panic travels to this thread (else the writer would wait for its chunk forever)
    type Res = std::thread::Result<bamwrite::ChunkOut>;
    let (tx, rx) = std::sync::mpsc::channel::<(usize, Res)>();
    let mut write_ns = 0u64;
    rayon::in_place_scope_fifo(|s| {
        let process = &process;
        let mut submitted = 0usize;
        let submit = |k: usize| {
            let tx = tx.clone();
            let ch = &chunks[k];
            s.spawn_fifo(move |_| {
                let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| process((k, ch))));
                let _ = tx.send((k, r));
            });
        };
        while submitted < chunks.len().min(window) {
            submit(submitted);
            submitted += 1;
        }
        let mut done: BTreeMap<usize, bamwrite::ChunkOut> = BTreeMap::new();
        let mut next = 0usize;
        while next < chunks.len() {
            let (k, r) = rx.recv().expect("pass 2 worker failed");
            match r {
                Ok(o) => {
                    done.insert(k, o);
                }
                Err(e) => std::panic::resume_unwind(e),
            }
            while let Some(o) = done.remove(&next) {
                let tw = Instant::now();
                bam.append(chunks[next].tid, &o);
                write_ns += tw.elapsed().as_nanos() as u64;
                drop(o);
                next += 1;
                if submitted < chunks.len() {
                    submit(submitted);
                    submitted += 1;
                }
            }
        }
    });
    let n_written = bam.n_records;
    bam.finish(&args.out);
    eprintln!("pass 2: {n_written} records written in {:.1}s (writer busy {:.1}s)", t2.elapsed().as_secs_f64(), write_ns as f64 * 1e-9);

    write_metrics(&args, &cmdline, &header.text, &hi, &lib_metrics.into_inner().unwrap(), &res);
    eprintln!("total {:.1}s", t0.elapsed().as_secs_f64());
}

/// Cuts chunks holding much more read data than a `chunk`-sized piece of the genome at mean density (by the index's
/// per-16 kb estimate, summed over the inputs) into that many equal pieces. Small contigs are one chunk each
/// whatever their depth (on GRCh38 + decoys, a 33 kb decoy and chrM hold 6x a 2 Mb chunk's reads), and one thread
/// each would finish both passes alone. Records are numbered in file order whatever the chunks, so the duplicates
/// are the same; only where the output's BGZF blocks end changes.
fn split_dense(chunks: Vec<Chunk>, inputs: &[Input], lens: &[i64], chunk: i64) -> Vec<Chunk> {
    let mut costs = inputs[0].window_costs(lens);
    for i in &inputs[1..] {
        for (a, b) in costs.iter_mut().zip(i.window_costs(lens)) {
            a.iter_mut().zip(b).for_each(|(x, y)| *x += y);
        }
    }
    let total: f64 = costs.iter().flatten().sum();
    let genome: i64 = lens.iter().sum();
    if total <= 0.0 || genome <= 0 {
        return chunks;
    }
    let target = total * chunk as f64 / genome as f64;
    let mut out = Vec::with_capacity(chunks.len());
    for ch in chunks {
        let from = ch.from.max(1);
        let k = if ch.tid >= 0 && ch.ivs.len() == 1 {
            (crate::plan::interval_cost(&costs[ch.tid as usize], from, ch.to) / target).ceil() as i64
        } else {
            1
        };
        let k = k.min((ch.to - from + 1) / 1000);
        if k < 2 {
            out.push(ch);
            continue;
        }
        let step = (ch.to - from + 1 + k - 1) / k;
        let mut a = from;
        let mut first = true;
        while a <= ch.to {
            let b = (a + step - 1).min(ch.to);
            out.push(Chunk { tid: ch.tid, from: if first { ch.from } else { a }, to: b, ivs: ch.ivs.clone(), pos0: first && ch.pos0 });
            first = false;
            a = b + 1;
        }
    }
    out
}

/// The records of one chunk in Picard's order: one input's in file order; several inputs merged as htsjdk's
/// `MergingSamRecordIterator` does (`SAMRecordCoordinateComparator`, ties to the earlier input), RG / PG IDs
/// translated to the merged header.
struct ChunkRecords<'a> {
    streams: Vec<crate::input::Stream>,
    /// the input each stream reads (inputs without data in the chunk get no stream: opening one starts a thread)
    which: Vec<usize>,
    /// next in-chunk record of each stream (None: exhausted)
    heads: Vec<Option<Vec<u8>>>,
    ch: Chunk,
    merged: Option<&'a qname::Merged>,
}

impl<'a> ChunkRecords<'a> {
    fn open(inputs: &[Input], tails: &[(u64, u64)], merged: Option<&'a qname::Merged>, ch: &Chunk, io_ns: &Arc<AtomicU64>) -> Option<ChunkRecords<'a>> {
        let mut streams = Vec::with_capacity(inputs.len());
        let mut which = Vec::with_capacity(inputs.len());
        for (k, (i, t)) in inputs.iter().zip(tails).enumerate() {
            if inputs.len() > 1 && ch.tid >= 0 && ch.query().iter().all(|&(a, b)| !i.has_data(ch.tid as usize, a, b)) {
                continue;
            }
            streams.push(bamwrite::open_chunk(i, ch, *t, io_ns)?);
            which.push(k);
        }
        let ch = ch.clone();
        let mut cr = ChunkRecords { heads: vec![None; streams.len()], streams, which, ch, merged };
        if cr.merged.is_some() {
            for k in 0..cr.streams.len() {
                cr.heads[k] = cr.pull(k);
            }
        }
        Some(cr)
    }

    /// next record of stream `k` that belongs to the chunk
    fn pull(&mut self, k: usize) -> Option<Vec<u8>> {
        let mut rec = Vec::new();
        while self.streams[k].next_record(&mut rec) {
            if bamwrite::in_chunk(&rec, &self.ch) {
                return Some(rec);
            }
        }
        None
    }

    fn next(&mut self, out: &mut Vec<u8>) -> bool {
        if self.merged.is_none() {
            while self.streams[0].next_record(out) {
                if bamwrite::in_chunk(out, &self.ch) {
                    return true;
                }
            }
            return false;
        }
        let mut best: Option<usize> = None;
        for k in 0..self.heads.len() {
            if let Some(r) = &self.heads[k] {
                if best.is_none_or(|b| crate::sortsam::coord_cmp(r, self.heads[b].as_ref().unwrap()) == std::cmp::Ordering::Less) {
                    best = Some(k);
                }
            }
        }
        let Some(k) = best else { return false };
        let next = self.pull(k);
        *out = std::mem::replace(&mut self.heads[k], next).unwrap();
        if let Some(m) = self.merged {
            qname::translate_record(m, self.which[k], out);
        }
        true
    }
}

/// The metrics file: one row per library of the header (Picard creates them up front) or with reads.
fn write_metrics(args: &Args, cmdline: &str, header_text: &str, hi: &HeaderInfo, lm: &[LibMetrics], res: &GroupResult) {
    let mut libs: BTreeMap<String, LibMetrics> = BTreeMap::new();
    // libraries of the header (Picard creates metrics for them up front), plus any library with reads
    let header_libs: std::collections::BTreeSet<&str> = header_text
        .lines()
        .filter(|l| l.starts_with("@RG\t"))
        .map(|l| l.split('\t').find_map(|f| f.strip_prefix("LB:")).unwrap_or("Unknown Library"))
        .collect();
    for (i, name) in hi.libs.iter().enumerate() {
        let m = &lm[i];
        let any = m.unpaired_reads_examined + m.read_pairs_examined + m.secondary_or_supplementary_rds + m.unmapped_reads > 0;
        if header_libs.contains(name.as_str()) || any {
            let mut m = m.clone();
            m.read_pair_optical_duplicates = *res.optical_by_lib.get(&(i as u16)).unwrap_or(&0);
            libs.insert(name.clone(), m);
        }
    }
    let to_b = |h: &HashMap<i64, i64>| -> BTreeMap<i64, i64> { h.iter().map(|(&k, &v)| (k, v)).collect() };
    metrics::write(&args.metrics, cmdline, &libs, &to_b(&res.all_sets), &to_b(&res.optical_sets), &to_b(&res.non_optical_sets));
}
