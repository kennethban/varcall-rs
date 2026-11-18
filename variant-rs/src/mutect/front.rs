//! Mutect2 front end of one contig shard (GATK `MultiIntervalLocalReadShard` + `AssemblyRegionIterator` with
//! `Mutect2Engine.isActive`):
//!
//! ```text
//! per BAM: BAI chunks -> inflate -> htsjdk interval filter -> Mutect2 read filters -> PalindromeArtifactClip
//!   -> merge (htsjdk MergingSamRecordIterator, SAMRecordCoordinateComparator) -> MutectDownsampler
//!   -> per-locus tumor/normal PileupQualBuffers -> isActive (0/1) -> band-pass activity profile -> regions
//!   -> region filling
//! ```
//! Filtering and clipping per file before the merge is equivalent to GATK's order (both are per read and
//! leave the merge keys unchanged).

use super::math::M2Math;
use crate::active::*;
use crate::bam::{self, Rec};
use crate::frontend::{FilledRegion, M2Quals, ReadRec, RefSeq, RegionFiller, REGION_PADDING};
use crate::input::Input;
use hc_downsample::{reservoir_downsample, JavaRandom};
use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::atomic::AtomicU64;
use std::sync::Arc;

pub const MIN_READ_LENGTH: usize = 30;
pub const READ_QUALITY_FILTER_THRESHOLD: u8 = 20;
pub const MIN_PALINDROME_SIZE: i64 = 5;
pub const MINIMUM_BASE_QUALITY: u8 = 6;
pub const INDEL_START_QUAL: i32 = 30;
pub const INDEL_CONTINUATION_QUAL: i32 = 10;
pub const MAX_ALT_FRACTION_IN_NORMAL: f64 = 0.3;
pub const MAX_NORMAL_QUAL_SUM: f64 = 100.0;
pub const PCR_SNV_QUAL: u8 = 40;
pub const MULTIPLE_SUBSTITUTION_QUAL_CORRECTION: f64 = 5.0; // round(10 * log10(3))
pub const SUSPICIOUS_MAPPING_QUALITY: u8 = 50;
pub const CALLABLE_DEPTH: u32 = 10;

/// One input file: its reader and read-group -> sample index map.
pub struct M2Input {
    pub input: Input,
    pub rg_sample: std::collections::HashMap<Vec<u8>, u16>,
    pub ref_lens: Vec<i64>,
    pub ref_names: Vec<String>,
    /// `--read-filter MateOnSameContigOrNoMappedMateReadFilter`
    pub mate_same_contig: bool,
    /// `--read-filter MateUnmappedAndUnmappedReadFilter`
    pub mate_mapped: bool,
    /// keep each read's original record (`--bam-output`)
    pub keep_records: bool,
    /// several inputs: their merged header (RG / PG tags of kept records are translated to it) and this input's
    /// index
    pub merged: Option<(Arc<crate::markdup::qname::Merged>, usize)>,
}

#[derive(Default, Debug, Clone)]
pub struct M2FilterCounts {
    pub total: u64,
    pub passed: u64,
    pub mapq: u64,
    pub mapq_unavailable: u64,
    pub unmapped: u64,
    pub secondary: u64,
    pub duplicate: u64,
    pub qcfail: u64,
    pub chimeric_oa: u64,
    pub zero_ref_len: u64,
    pub short: u64,
    pub bad_cigar: u64,
    pub malformed: u64,
    pub palindrome_clipped: u64,
    /// `--read-filter` filters
    pub user: u64,
}

impl M2FilterCounts {
    pub fn add(&mut self, o: &M2FilterCounts) {
        self.total += o.total;
        self.passed += o.passed;
        self.user += o.user;
        self.mapq += o.mapq;
        self.mapq_unavailable += o.mapq_unavailable;
        self.unmapped += o.unmapped;
        self.secondary += o.secondary;
        self.duplicate += o.duplicate;
        self.qcfail += o.qcfail;
        self.chimeric_oa += o.chimeric_oa;
        self.zero_ref_len += o.zero_ref_len;
        self.short += o.short;
        self.bad_cigar += o.bad_cigar;
        self.malformed += o.malformed;
        self.palindrome_clipped += o.palindrome_clipped;
    }
}

/// The aux tags the front end needs: RG, OA, XM (`AddOriginalAlignmentTags.MATE_CONTIG_TAG_NAME`).
pub(crate) struct Aux<'a> {
    pub(crate) rg: Option<&'a [u8]>,
    oa: Option<&'a [u8]>,
    xm: Option<&'a [u8]>,
}

pub(crate) fn aux_tags(b: &[u8]) -> Aux<'_> {
    let mut o = crate::bamwrite::aux_start(b);
    let mut a = Aux { rg: None, oa: None, xm: None };
    let size = |t: u8| -> usize {
        match t {
            b'A' | b'c' | b'C' => 1,
            b's' | b'S' => 2,
            b'i' | b'I' | b'f' => 4,
            _ => 0,
        }
    };
    while o + 3 <= b.len() {
        let tag = [b[o], b[o + 1]];
        let ty = b[o + 2];
        o += 3;
        match ty {
            b'Z' | b'H' => {
                let s = o;
                while b[o] != 0 {
                    o += 1;
                }
                if ty == b'Z' {
                    match &tag {
                        b"RG" => a.rg = Some(&b[s..o]),
                        b"OA" => a.oa = Some(&b[s..o]),
                        b"XM" => a.xm = Some(&b[s..o]),
                        _ => {}
                    }
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
    a
}

/// `Mutect2Engine.makeStandardMutect2ReadFilters()` (counts attribute a read to the first failing filter).
fn passes_m2_filters(r: &Rec, aux: &Aux, inp: &M2Input, fc: &mut M2FilterCounts) -> bool {
    let ref_lens = &inp.ref_lens;
    fc.total += 1;
    let unmapped = r.flag & 0x4 != 0 || r.ref_id < 0 || r.start < 1;
    if r.mapq < READ_QUALITY_FILTER_THRESHOLD {
        if unmapped { fc.unmapped += 1 } else { fc.mapq += 1 }
        return false;
    }
    if r.mapq == 255 {
        fc.mapq_unavailable += 1;
        return false;
    }
    if unmapped {
        fc.unmapped += 1;
        return false;
    }
    if r.flag & 0x100 != 0 {
        fc.secondary += 1;
        return false;
    }
    if r.flag & 0x400 != 0 {
        fc.duplicate += 1;
        return false;
    }
    if r.flag & 0x200 != 0 {
        fc.qcfail += 1;
        return false;
    }
    // NON_CHIMERIC_ORIGINAL_ALIGNMENT_READ_FILTER: OA contig must equal the XM (original mate contig) tag
    if let (Some(oa), Some(xm)) = (aux.oa, aux.xm) {
        let oa_contig = oa.split(|&c| c == b',').next().unwrap_or(b"");
        if oa_contig != xm {
            fc.chimeric_oa += 1;
            return false;
        }
    }
    if ref_length(&r.cigar) == 0 {
        fc.zero_ref_len += 1;
        return false;
    }
    if r.seq.len() < MIN_READ_LENGTH {
        fc.short += 1;
        return false;
    }
    if !good_cigar(&r.cigar) {
        fc.bad_cigar += 1;
        return false;
    }
    let wellformed = r.has_rg
        && !r.seq.is_empty()
        && !r.qual_missing
        && r.qual.len() == r.seq.len()
        && read_length(&r.cigar) == r.seq.len()
        && !r.cigar.iter().any(|e| e.0 == OP_N)
        && (r.ref_id as usize) < ref_lens.len()
        && r.start <= ref_lens[r.ref_id as usize];
    if !wellformed {
        fc.malformed += 1;
        return false;
    }
    // user read filters (after the tool's defaults)
    let paired = r.flag & 0x1 != 0;
    let mate_unmapped = r.flag & 0x8 != 0;
    if inp.mate_same_contig && !(!paired || mate_unmapped || r.ref_id == r.next_ref_id) {
        fc.user += 1;
        return false;
    }
    if inp.mate_mapped && paired && mate_unmapped {
        fc.user += 1;
        return false;
    }
    fc.passed += 1;
    true
}

/// `ReadUtils.getAdaptorBoundary` (None = CANNOT_COMPUTE_ADAPTOR_BOUNDARY).
fn adaptor_boundary(r: &ReadRec) -> Option<i64> {
    let paired = r.flag & 0x1 != 0;
    if r.tlen == 0 || !paired || r.flag & 0x4 != 0 || r.flag & 0x8 != 0 || r.mate_start == 0 {
        return None;
    }
    let rev = r.flag & 0x10 != 0;
    if rev == (r.flag & 0x20 != 0) {
        return None;
    }
    let well = if rev { r.end > r.mate_start } else { r.start <= r.mate_start + r.tlen as i64 };
    if !well {
        return None;
    }
    Some(if rev { r.mate_start - 1 } else { r.start + (r.tlen as i64).abs() })
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        x => x,
    }
}

/// `PalindromeArtifactClipReadTransformer(ref, 5)`. Returns true if the read was clipped.
fn palindrome_clip(r: &mut ReadRec, refseq: &RefSeq, scratch: &mut Vec<hc_preproc::CigarElem>) -> bool {
    let Some(boundary) = adaptor_boundary(r) else { return false };
    if r.flag & 0x3 != 0x3 {
        return false;
    }
    let (first, last) = (r.cigar[0], *r.cigar.last().unwrap());
    let upstream = r.tlen > 0;
    if (upstream && first.0 != OP_S && first.0 != OP_I) || (!upstream && last.0 != OP_S && last.0 != OP_I) {
        return false;
    }
    let count = if upstream { first.1 } else { last.1 } as i64;
    let len = r.seq.len() as i64;
    let n_cmp = (count + MIN_PALINDROME_SIZE).min(len);
    let (ref_start, ref_end) = if upstream { (boundary - n_cmp, boundary - 1) } else { (boundary + 1, boundary + n_cmp) };
    if ref_start < 1 || ref_end > refseq.contig_len {
        return false;
    }
    if (upstream && ref_start < r.start) || (!upstream && r.end < ref_end) {
        return false;
    }
    let mut idx = if upstream { n_cmp - 1 } else { len - 1 };
    let mut matches = 0i64;
    for pos in ref_start..=ref_end {
        if complement(refseq.base(pos)) == r.seq[idx as usize] {
            matches += 1;
        }
        idx -= 1;
    }
    if (matches as f64) / (n_cmp as f64) < 0.9 {
        return false;
    }
    let op = if upstream { (0, count - 1) } else { (len - count, len) };
    let mut p = r.to_preproc();
    let ok = hc_preproc::clip_read_hard(&mut p, &[op], scratch);
    assert!(ok, "palindrome clipping emptied read {}", p.name);
    r.cigar = p.cigar.iter().map(|&(c, l)| (b"MIDNSHP=X".iter().position(|&x| x == c).unwrap() as u8, l)).collect();
    r.seq = p.seq;
    r.qual = p.qual;
    r.start = p.start;
    r.end = p.end;
    true
}

/// `MutectDownsampler(maxReadsPerAlignmentStart, maxSuspicious = 0 (off), stride = 1)`.
struct MutectDownsampler {
    max: usize,
    pending: Vec<Arc<ReadRec>>,
    stride_start: Option<i64>,
    pub draws_groups: u64,
    /// dry run: only record, per drawing group, the number of reads the reservoir sees (no reads are emitted)
    sizes: Option<Vec<usize>>,
}

impl MutectDownsampler {
    fn submit(&mut self, r: Arc<ReadRec>, rng: &mut JavaRandom, out: &mut Vec<Arc<ReadRec>>) {
        match self.stride_start {
            Some(s) if r.start >= s + 1 => {
                self.finalize(rng, out);
                self.stride_start = Some(r.start);
            }
            None => self.stride_start = Some(r.start),
            _ => {}
        }
        self.pending.push(r);
    }

    fn finalize(&mut self, rng: &mut JavaRandom, out: &mut Vec<Arc<ReadRec>>) {
        if let Some(sizes) = self.sizes.as_mut() {
            if self.pending.len() > self.max {
                sizes.push(self.pending.iter().filter(|r| r.mapq > SUSPICIOUS_MAPPING_QUALITY).count());
            }
            self.pending.clear();
            return;
        }
        if self.pending.len() <= self.max {
            out.append(&mut self.pending);
        } else {
            self.draws_groups += 1;
            if std::env::var_os("VARIANT_RS_DEBUG_DS").is_some() {
                eprintln!("drawing group at {} ({} reads)", self.pending[0].start, self.pending.len());
            }
            let dbg = std::env::var("VARIANT_RS_DS_AT").ok().is_some_and(|v| v.split(',').any(|p| p.parse::<i64>().ok() == Some(self.pending[0].start)));
            if dbg {
                eprintln!("dsgroup {} pending={}", self.pending[0].start, self.pending.len());
                for r in &self.pending {
                    eprintln!("dsin {} {} {} mq={}", String::from_utf8_lossy(&r.name), r.flag, r.start, r.mapq);
                }
            }
            let well: Vec<Arc<ReadRec>> = self.pending.drain(..).filter(|r| r.mapq > SUSPICIOUS_MAPPING_QUALITY).collect();
            let kept = reservoir_downsample(well, self.max, rng);
            if dbg {
                for r in &kept {
                    eprintln!("dskeep {} {} {}", String::from_utf8_lossy(&r.name), r.flag, r.start);
                }
            }
            out.extend(kept);
        }
    }
}

/// Dry run of [`MutectDownsampler`]: per stride (alignment start) the number of reads and of reads with MQ > 50;
/// records the latter for every stride that needs a draw.
struct DryGroups {
    max: usize,
    stride_start: Option<i64>,
    n: usize,
    well: usize,
    sizes: Vec<usize>,
}

impl DryGroups {
    fn submit(&mut self, start: i64, mapq: u8) {
        match self.stride_start {
            Some(s) if start >= s + 1 => {
                self.finalize();
                self.stride_start = Some(start);
            }
            None => self.stride_start = Some(start),
            _ => {}
        }
        self.n += 1;
        if mapq > SUSPICIOUS_MAPPING_QUALITY {
            self.well += 1;
        }
    }

    fn finalize(&mut self) {
        if self.n > self.max {
            self.sizes.push(self.well);
        }
        self.n = 0;
        self.well = 0;
    }
}

/// Per-locus pileup accumulator: tumor / normal pileup sizes, `PileupQualBuffer` sums, and the tumor's
/// alt quals in pileup order (buffer index, Java byte qual).
#[derive(Clone, Default)]
struct MAcc {
    n_t: u32,
    n_n: u32,
    t_sum: [i64; 6],
    n_sum: [i64; 6],
    n_cnt: [u32; 6],
    t_alts: Vec<(u8, i8)>,
    /// substitution entries whose qual may since have been lowered by overlap cleaning:
    /// (tumor, index into t_alts or buffer index, read, read offset)
    lazy: Vec<(bool, u32, Arc<ReadRec>, u32)>,
    /// F1R2 counts of the tumor elements with qual > 20 (as of now)
    f1r2c: super::orientation::LocusCounts,
}

/// The per-locus accumulators, shared by the read stage (which scatters reads into the slots of loci not yet
/// evaluated) and the evaluation stage (which takes, evaluates and resets the slots of loci up to its frontier).
/// The protocol in [`run`] guarantees that the two stages never use the same slot at the same time.
/// The ring holds `window` loci (a power of two, initially WINDOW), more than twice the longest read span: a longer
/// read grows it first (`Front::grow`).
struct Ring(std::cell::UnsafeCell<(Box<[std::cell::UnsafeCell<MAcc>]>, usize)>);
// SAFETY: slots are only accessed under the read / evaluation protocol of `run` (disjoint slots, ordered by atomics);
// the buffer is only replaced by `grow`, while neither stage uses it.
unsafe impl Sync for Ring {}
impl Ring {
    fn new() -> Ring {
        Ring(std::cell::UnsafeCell::new((Self::cells(WINDOW), WINDOW)))
    }
    fn cells(n: usize) -> Box<[std::cell::UnsafeCell<MAcc>]> {
        (0..n).map(|_| std::cell::UnsafeCell::new(MAcc::default())).collect()
    }
    fn window(&self) -> usize {
        // SAFETY: `grow` does not run concurrently with any other use of the ring
        unsafe { (*self.0.get()).1 }
    }
    /// Whether a read spanning `span` bp needs a larger ring first.
    fn too_short_for(&self, span: i64) -> bool {
        span >= self.window() as i64 / 2
    }
    /// Grows the ring to hold a read spanning `span` bp, keeping the slots of the live loci `lo..lo + window` (every
    /// slot of a locus < `lo` has been taken, nothing beyond that range has been scattered). Returns the old window.
    /// SAFETY: neither stage may be using the ring (a quiescent point of the protocol of `run`).
    unsafe fn grow(&self, lo: i64, span: i64) -> usize {
        let (old, w0) = &mut *self.0.get();
        let w = (2 * (span as usize + 1)).next_power_of_two().max(*w0);
        let mut new = Self::cells(w);
        for pos in lo..lo + *w0 as i64 {
            *new[(pos as usize) % w].get_mut() = std::mem::take(old[(pos as usize) % *w0].get_mut());
        }
        let w_old = *w0;
        *self.0.get() = (new, w);
        w_old
    }
    /// SAFETY: the caller must own the slot of `pos` under the protocol of `run`.
    #[allow(clippy::mut_from_ref)]
    unsafe fn slot(&self, pos: i64) -> &mut MAcc {
        let (b, w) = &*self.0.get();
        &mut *b[(pos as usize) % *w].get()
    }
}

/// `PileupQualBuffer` qual of a substitution with (possibly cleaned) base qual `q` (not overlapping its mate).
fn substitution_qual(idx: usize, q: u8) -> i8 {
    if idx < OTHER_SUBSTITUTION {
        ((q as i8 as f64 + MULTIPLE_SUBSTITUTION_QUAL_CORRECTION).min(254.0) as i32) as i8
    } else {
        q as i8
    }
}

const WINDOW: usize = 1 << 18;

/// Read stage -> evaluation stage (two-stage front end).
enum Batch {
    Reads(Vec<Arc<ReadRec>>),
    /// the next read spans this many bp: grow the ring, then reply
    Grow(i64),
}
const OTHER_SUBSTITUTION: usize = 4;
const INDEL: usize = 5;

fn indel_qual(len: u32) -> i8 {
    (INDEL_START_QUAL as i64 + (len as i64 - 1) * INDEL_CONTINUATION_QUAL as i64).min(127) as i8
}

fn best_index(sums: &[i64; 6]) -> usize {
    let (mut best, mut best_sum) = (0, 0i64);
    for (n, &s) in sums.iter().enumerate() {
        if s > best_sum {
            best_sum = s;
            best = n;
        }
    }
    best
}

/// Where filled regions go.
pub trait M2Sink {
    /// A region at `callRegion` time, with its reads' qualities as of then (None = unchanged).
    fn region(&mut self, fr: FilledRegion, quals: Vec<Option<Arc<M2Quals>>>);
}

#[derive(Default)]
pub struct M2FrontResult {
    pub regions: Vec<Region>,
    pub filters: M2FilterCounts,
    pub n_loci: u64,
    pub n_active_loci: u64,
    pub callable: u64,
    pub assigned: u64,
    pub reads_in: u64,
    pub reads_kept: u64,
    pub drawing_groups: u64,
    /// dry run: reservoir input sizes of the drawing groups, in order
    pub group_sizes: Vec<usize>,
}

struct Front<'a, S: M2Sink> {
    t: &'a Tables,
    math: &'a M2Math,
    is_tumor: &'a [bool],
    has_normal: bool,
    contig: &'a str,
    intervals: &'a [(i64, i64)],
    refseq: &'a RefSeq,
    sink: S,
    dump: Option<&'a mut BufWriter<File>>,
    res: M2FrontResult,
    ring: &'a Ring,
    /// F1R2 elements taken out by overlap cleaning before their locus was evaluated (evaluation stage only):
    /// (base, mapq, F1R2, indel-adjacent) by locus slot
    f1r2_delta: Vec<Vec<(u8, u8, bool, bool)>>,
    /// highest locus whose slot the read stage may not reuse yet (published after each read)
    frontier: &'a std::sync::atomic::AtomicI64,
    fin_iv: usize,
    next_final: i64,
    profile: Profile,
    regions: Vec<Region>,
    filler: RegionFiller,
    alt_scratch: Vec<i8>,
    filled: VecDeque<FilledRegion>,
    polled: usize,
    applied: usize,
    f1r2: Option<&'a mut super::orientation::F1R2Counts>,
    /// last locus whose isActive (and F1R2 counts) has been evaluated
    evaluated_upto: i64,
    /// VARIANT_RS_TRACE_FRONT: inside `finalize`
    in_fin: bool,
    pon: Option<&'a super::germline::Pon>,
    germline: Option<&'a super::germline::Germline>,
    opts: super::M2Opts,
}

impl<'a, S: M2Sink> Front<'a, S> {
    fn emit(&mut self, upto_stream: i64, eof: bool) {
        let t = ft_now();
        self.emit_inner(upto_stream, eof);
        ft_add(8, t);
        if self.in_fin {
            ft_add(9, t);
        }
    }

    fn emit_inner(&mut self, upto_stream: i64, eof: bool) {
        let Front { filler, regions, filled, res, .. } = self;
        filler.drain(regions, upto_stream, eof, &mut |fr: FilledRegion| {
            res.assigned += fr.reads.len() as u64;
            filled.push_back(fr)
        });
    }

    /// `Mutect2Engine.callRegion` up to its first branch: `cleanOverlappingReadPairs` (in place, on reads shared
    /// with other regions and with the pileups of later loci), then hand the region to the sink with a snapshot of
    /// its reads' qualities.
    fn apply(&mut self, k: usize) {
        let t = ft_now();
        self.apply_inner(k);
        ft_add(10, t);
        if self.in_fin {
            ft_add(11, t);
        }
    }

    fn apply_inner(&mut self, k: usize) {
        let fr = self.filled.pop_front().expect("region applied before it was filled");
        assert_eq!(fr.id, k, "regions applied out of order");
        assert_eq!(self.applied, k);
        self.applied += 1;
        let t_clean = ft_now();
        let capped = clean_overlapping_read_pairs(&fr.reads, self.is_tumor.len(), self.f1r2.is_some().then_some(self.evaluated_upto));
        let t_u = ft_now();
        self.f1r2_uncount(&capped);
        ft_add(16, t_u);
        ft_add(13, t_clean);
        let snap: Vec<Option<Arc<M2Quals>>> =
            fr.reads.iter().map(|r| if r.cleaned.load(std::sync::atomic::Ordering::Relaxed) { r.m2q.lock().unwrap().clone() } else { None }).collect();
        let t_send = ft_now();
        self.sink.region(fr, snap);
        ft_add(14, t_send);
    }

    /// `Mutect2Engine.isActive` for one locus (tumor-normal, no PoN, no force-calling).
    fn is_active(&mut self, acc: &mut MAcc, pos: i64) -> f64 {
        // qualities lowered in place by the overlap cleaning of regions already processed
        for (tumor, i, r, off) in acc.lazy.drain(..) {
            let g = r.m2q.lock().unwrap();
            let Some(st) = g.as_ref() else { continue };
            let q = st.qual[off as usize];
            if q == r.qual[off as usize] {
                continue;
            }
            if tumor {
                let (idx, old) = acc.t_alts[i as usize];
                let new = substitution_qual(idx as usize, q);
                acc.t_alts[i as usize].1 = new;
                acc.t_sum[idx as usize] += new as i64 - old as i64;
            } else {
                let idx = (i >> 8) as usize;
                let old = substitution_qual(idx, r.qual[off as usize]);
                acc.n_sum[idx] += substitution_qual(idx, q) as i64 - old as i64;
            }
        }
        let depth = acc.n_t + acc.n_n;
        if depth == 0 {
            return 0.0;
        }
        if depth >= CALLABLE_DEPTH {
            self.res.callable += 1;
        }
        // F1R2CountsCollector.process(tumorPileup, ref), with the qualities as of now
        if let Some(c) = self.f1r2.as_deref_mut() {
            let rs = self.refseq;
            let ctx = if pos > 1 && pos < rs.contig_len { Some([rs.base(pos - 1), rs.base(pos), rs.base(pos + 1)]) } else { None };
            c.process(ctx, &acc.f1r2c);
        }
        let tb = best_index(&acc.t_sum);
        self.alt_scratch.clear();
        self.alt_scratch.extend(acc.t_alts.iter().filter(|a| a.0 as usize == tb).map(|a| a.1));
        let n_alt = self.alt_scratch.len() as i64;
        let tumor_log_odds = self.math.log_likelihood_ratio(acc.n_t as i64 - n_alt, &self.alt_scratch);
        if tumor_log_odds < self.opts.initial_log_odds {
            return 0.0;
        }
        if self.opts.genotype_germline {
        } else if self.has_normal {
            let nb = best_index(&acc.n_sum);
            if nb == tb {
                let normal_alt_count = acc.n_cnt[nb];
                if (normal_alt_count as f64) > acc.n_n as f64 * MAX_ALT_FRACTION_IN_NORMAL && (acc.n_sum[nb] as f64) > MAX_NORMAL_QUAL_SUM {
                    return 0.0;
                }
            }
        } else if self.germline.is_some_and(|g| g.blocks_tumor_only(self.contig, pos, tb, self.opts.max_population_af)) {
            return 0.0;
        }
        if !self.opts.genotype_pon && self.pon.is_some_and(|p| p.overlaps(self.contig, pos)) {
            return 0.0;
        }
        if depth < CALLABLE_DEPTH {
            self.res.callable += 1;
        }
        1.0
    }

    /// Finalize the interval loci up to `upto`: isActive, then the activity profile.
    fn finalize(&mut self, upto: i64) {
        let t = ft_now();
        self.in_fin = true;
        self.finalize_inner(upto);
        self.in_fin = false;
        ft_add(7, t);
    }

    fn finalize_inner(&mut self, upto: i64) {
        while self.fin_iv < self.intervals.len() {
            let (s, e) = self.intervals[self.fin_iv];
            if self.next_final < s {
                self.next_final = s;
            }
            if self.next_final > e {
                self.fin_iv += 1;
                continue;
            }
            if self.next_final > upto {
                return;
            }
            let pos = self.next_final;
            // SAFETY: the read stage has scattered every read starting before `pos` and writes no slot <= the frontier
            let mut acc = std::mem::take(unsafe { self.ring.slot(pos) });
            let w = self.f1r2_delta.len();
            for (b, mq, f, adj) in self.f1r2_delta[(pos as usize) % w].drain(..) {
                acc.f1r2c.update(b, mq, f, adj, false);
            }
            let prob = self.is_active(&mut acc, pos);
            self.evaluated_upto = pos;
            self.res.n_loci += 1;
            if prob > 0.0 {
                self.res.n_active_loci += 1;
            }
            if let Some(d) = self.dump.as_mut() {
                writeln!(d, "{}\t{}\t{}\t{}\t{}", self.contig, pos, prob, acc.n_t, acc.n_n).unwrap();
            }
            if !self.profile.is_empty() {
                let force = pos != self.profile.stop + 1;
                self.profile.pop_ready(force, &mut self.regions);
            }
            self.profile.add(pos, prob, None);
            // AssemblyRegionIterator: after this locus, poll the next pending region if the locus is past its
            // padded span; the walker then processes the region polled before it (callRegion)
            if self.polled < self.regions.len() {
                let r = self.regions[self.polled];
                if pos > (r.end + REGION_PADDING).min(self.refseq.contig_len) {
                    self.polled += 1;
                    if self.polled >= 2 {
                        self.emit(pos, false);
                        self.apply(self.polled - 2);
                    }
                }
            }
            self.next_final += 1;
            // reset the slot, keeping the allocation of its alt list
            let mut alts = acc.t_alts;
            alts.clear();
            let mut lazy = acc.lazy;
            lazy.clear();
            let mut f1r2c = acc.f1r2c;
            f1r2c.clear();
            *unsafe { self.ring.slot(pos) } = MAcc { t_alts: alts, lazy, f1r2c, ..MAcc::default() };
        }
    }

    /// One read leaving the downsampler: region read cache, loci before it, pileup scatter.
    fn process(&mut self, r: Arc<ReadRec>) {
        let t = ft_now();
        self.process_inner(r);
        ft_add(5, t);
    }

    fn process_inner(&mut self, r: Arc<ReadRec>) {
        let (rstart, rend) = (r.start, r.end);
        assert!(!self.ring.too_short_for(rend - rstart + 1), "read {}: ring not grown", String::from_utf8_lossy(&r.name));
        self.emit(rstart, false);
        self.filler.cache.push_back(r);
        self.finalize(rstart - 1);
        self.frontier.store(rstart - 1, std::sync::atomic::Ordering::Release);
    }
}

/// The read stage's scatter of one read into the locus accumulators (`MAcc`) of the loci it covers.
struct Scatter<'a> {
    ring: &'a Ring,
    intervals: &'a [(i64, i64)],
    is_tumor: &'a [bool],
    refseq: &'a RefSeq,
    want_f1r2: bool,
    scatter_iv: usize,
    local_ivs: Vec<(i64, i64)>,
}

impl Scatter<'_> {
    /// SAFETY (ring): the caller has waited until the evaluation frontier allows writing the slots of `r`'s loci.
    fn scatter(&mut self, r: &Arc<ReadRec>) {
        let (rstart, rend) = (r.start, r.end);
        while self.scatter_iv < self.intervals.len() && self.intervals[self.scatter_iv].1 < rstart {
            self.scatter_iv += 1;
        }
        self.local_ivs.clear();
        for &(s, e) in &self.intervals[self.scatter_iv..] {
            if s > rend {
                break;
            }
            self.local_ivs.push((s.max(rstart), e.min(rend)));
        }
        if self.local_ivs.is_empty() {
            return;
        }
        let local = &self.local_ivs;
        let (lo, hi) = (local[0].0, local.last().unwrap().1);
        let gaps = local.len() > 1;
        let in_iv = |pos: i64| -> bool { pos >= lo && pos <= hi && (!gaps || local.iter().any(|&(s, e)| pos >= s && pos <= e)) };

        let tumor = self.is_tumor[r.sample as usize];
        let adaptor = r.adaptor_exclusion();
        let want_f1r2 = tumor && self.want_f1r2;
        // ReadUtils.isF1R2: isReverseStrand != isFirstOfPair (first of pair = paired and flag 0x40)
        let is_f1r2 = is_f1r2(r);
        // PileupQualBuffer: overlapsMate = mateStart <= pos < mateStart + read length (properly paired, mate mapped)
        let mate_start = if r.flag & 0x3 != 0x3 || r.flag & 0x8 != 0 || r.mate_ref_id < 0 || r.mate_start == 0 { i64::MAX } else { r.mate_start };
        let read_len = r.seq.len() as i64;
        let (refseq, ring) = (self.refseq, self.ring);
        let n = r.cigar.len();
        let mut rp = rstart;
        let mut ro = 0usize;
        for (ci, &(op, len)) in r.cigar.iter().enumerate() {
            let len_u = len as usize;
            match op {
                OP_M | OP_EQ | OP_X | OP_D => {
                    let is_del = op == OP_D;
                    let adj_prev = if ci > 0 { Some(r.cigar[ci - 1].0) } else { None };
                    let adj_next = if ci + 1 < n { Some(r.cigar[ci + 1]) } else { None };
                    // length of the indel immediately after the element's last base
                    let following = if is_del {
                        0
                    } else if let Some(d) = r.cigar[ci + 1..].iter().find(|e| on_genome(e.0)).filter(|e| e.0 == OP_D) {
                        d.1
                    } else if let Some((OP_I, l)) = adj_next {
                        l
                    } else {
                        0
                    };
                    for k in 0..len_u {
                        let pos = rp + k as i64;
                        if !in_iv(pos) {
                            continue;
                        }
                        if let Some((rev, b)) = adaptor {
                            if (rev && pos <= b) || (!rev && pos >= b) {
                                continue;
                            }
                        }
                        // SAFETY: see `scatter`
                        let acc = unsafe { ring.slot(pos) };
                        if tumor {
                            acc.n_t += 1;
                        } else {
                            acc.n_n += 1;
                        }
                        let (at_s, at_e) = (k == 0, k == len_u - 1);
                        if want_f1r2 && !is_del {
                            let i = ro + k;
                            // original qualities: overlap cleaning takes capped bases out again (`f1r2_uncount`)
                            if r.qual[i] > F1R2_MIN_BASE_QUALITY {
                                let after_ins = at_s && adj_prev == Some(OP_I);
                                let before_del = at_e && r.cigar[ci + 1..].iter().find(|e| on_genome(e.0)).is_some_and(|e| e.0 == OP_D);
                                acc.f1r2c.update(r.seq[i], r.mapq, is_f1r2, after_ins || before_del, true);
                            }
                        }
                        let entry: Option<(usize, i8)> = if is_del {
                            Some((INDEL, indel_qual(len)))
                        } else if at_e && following > 0 {
                            Some((INDEL, indel_qual(following)))
                        } else {
                            let i = ro + k;
                            let q = r.qual[i];
                            let useful_clip = q > MINIMUM_BASE_QUALITY
                                && ((at_e && adj_next.map(|e| e.0) == Some(OP_S) && r.qual[i + 1] > MINIMUM_BASE_QUALITY)
                                    || (at_s && adj_prev == Some(OP_S) && r.qual[i - 1] > MINIMUM_BASE_QUALITY));
                            if useful_clip {
                                Some((INDEL, indel_qual(1)))
                            } else {
                                let base = r.seq[i];
                                if base != refseq.base(pos) && q > MINIMUM_BASE_QUALITY {
                                    let overlaps = mate_start <= pos && pos < mate_start.saturating_add(read_len);
                                    let q2 = if overlaps { q.min(PCR_SNV_QUAL / 2) } else { q };
                                    let idx = match base {
                                        b'A' | b'a' => 0,
                                        b'C' | b'c' => 1,
                                        b'G' | b'g' => 2,
                                        b'T' | b't' => 3,
                                        _ => OTHER_SUBSTITUTION,
                                    };
                                    if !overlaps && q > PCR_SNV_QUAL / 2 {
                                        // cleaning may cap this base at 20 before the locus is evaluated
                                        let li = if tumor { acc.t_alts.len() as u32 } else { (idx as u32) << 8 };
                                        acc.lazy.push((tumor, li, r.clone(), i as u32));
                                    }
                                    Some((idx, substitution_qual(idx, q2)))
                                } else {
                                    None
                                }
                            }
                        };
                        if let Some((idx, q)) = entry {
                            if tumor {
                                acc.t_sum[idx] += q as i64;
                                acc.t_alts.push((idx as u8, q));
                            } else {
                                acc.n_sum[idx] += q as i64;
                                acc.n_cnt[idx] += 1;
                            }
                        }
                    }
                    rp += len as i64;
                    if !is_del {
                        ro += len_u;
                    }
                }
                OP_I | OP_S => ro += len_u,
                OP_N => rp += len as i64,
                _ => {}
            }
        }
    }
}

impl<'a, S: M2Sink> Front<'a, S> {
    /// Grows the ring and the F1R2 deltas (same slots) for a read spanning `span` bp.
    /// SAFETY: as [`Ring::grow`]; every read before it has been processed here and scattered.
    unsafe fn grow(&mut self, span: i64) {
        let lo = self.next_final;
        let w0 = self.ring.grow(lo, span);
        let w = self.ring.window();
        let mut d = vec![Vec::new(); w];
        for pos in lo..lo + w0 as i64 {
            d[(pos as usize) % w] = std::mem::take(&mut self.f1r2_delta[(pos as usize) % w0]);
        }
        self.f1r2_delta = d;
    }

    /// Overlap cleaning capped the qualities of `capped` (read, read index) at 20: take those bases out of the F1R2
    /// counts of loci that have not been evaluated yet (the read stage counts every base with its original quality,
    /// before or after this cleaning).
    fn f1r2_uncount(&mut self, capped: &[(Arc<ReadRec>, Vec<usize>)]) {
        for (r, idxs) in capped {
            if !self.is_tumor[r.sample as usize] {
                continue;
            }
            let adaptor = r.adaptor_exclusion();
            let f1r2 = is_f1r2(r);
            for &idx in idxs {
                let Some((pos, indel_adj)) = element_at(r, idx) else { continue };
                if pos <= self.evaluated_upto {
                    continue;
                }
                let k = self.intervals.partition_point(|iv| iv.1 < pos);
                if !(k < self.intervals.len() && self.intervals[k].0 <= pos) {
                    continue;
                }
                if let Some((rev, b)) = adaptor {
                    if (rev && pos <= b) || (!rev && pos >= b) {
                        continue;
                    }
                }
                let w = self.f1r2_delta.len();
                self.f1r2_delta[(pos as usize) % w].push((r.seq[idx], r.mapq, f1r2, indel_adj));
            }
        }
    }
}

/// `F1R2CountsCollector` keeps bases with qual > `CollectF1R2CountsArgumentCollection.minBaseQuality` (20).
const F1R2_MIN_BASE_QUALITY: u8 = 20;

/// `ReadUtils.isF1R2`: isReverseStrand != isFirstOfPair (first of pair = paired and flag 0x40)
fn is_f1r2(r: &ReadRec) -> bool {
    (r.flag & 0x10 != 0) != (r.flag & 0x1 != 0 && r.flag & 0x40 != 0)
}

/// The pileup element of read index `idx`: (reference position, isAfterInsertion || isBeforeDeletionStart), or
/// None if the base is not aligned (insertion / soft clip).
fn element_at(r: &ReadRec, idx: usize) -> Option<(i64, bool)> {
    let (mut rp, mut ro) = (r.start, 0usize);
    for (ci, &(op, len)) in r.cigar.iter().enumerate() {
        let l = len as usize;
        match op {
            OP_M | OP_EQ | OP_X => {
                if idx < ro + l {
                    let k = idx - ro;
                    let after_ins = k == 0 && ci > 0 && r.cigar[ci - 1].0 == OP_I;
                    let before_del = k == l - 1 && r.cigar[ci + 1..].iter().find(|e| on_genome(e.0)).is_some_and(|e| e.0 == OP_D);
                    return Some((rp + k as i64, after_ins || before_del));
                }
                rp += len as i64;
                ro += l;
            }
            OP_D | OP_N => rp += len as i64,
            OP_I | OP_S => {
                if idx < ro + l {
                    return None;
                }
                ro += l;
            }
            _ => {}
        }
    }
    None
}

/// `VARIANT_RS_TRACE_FRONT=1`: front-end time split (nanoseconds; see `FT_NAMES`), summed over contigs and threads.
pub static FT: [AtomicU64; 17] = [const { AtomicU64::new(0) }; 17];
pub const FT_NAMES: [&str; 17] = [
    "inflate (bgzf threads)",
    "decode+filter+clip (decoder threads)",
    "front: wait for decoders",
    "front: merge loop total",
    "front: downsample",
    "front: process total",
    "front: scatter",
    "front: finalize total",
    "front: emit (region filling) total",
    "front: emit inside finalize",
    "front: apply (cleaning, snapshot, send) total",
    "front: apply inside finalize",
    "front: end-of-contig flush",
    "front: apply: overlap cleaning",
    "front: apply: hand-off to region workers",
    "front: apply: cleaning: pair adjustments",
    "front: apply: F1R2 subtraction",
];
/// VARIANT_RS_TRACE_FRONT counters: pairs skipped (already adjusted), pairs adjusted
pub static FTC: [AtomicU64; 2] = [const { AtomicU64::new(0) }; 2];
pub fn ft_on() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("VARIANT_RS_TRACE_FRONT").is_some())
}
#[inline]
fn ft_add(i: usize, t: Option<std::time::Instant>) {
    if let Some(t) = t {
        FT[i].fetch_add(t.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
    }
}
#[inline]
fn ft_now() -> Option<std::time::Instant> {
    if ft_on() { Some(std::time::Instant::now()) } else { None }
}

/// One raw BAM record -> merge item: Mutect2 read filters, read record, palindrome clipping (dry run: only the
/// merge keys and the start after clipping).
#[allow(clippy::too_many_arguments)]
fn decode_record(body: &[u8], rec: &mut Rec, scratch: &mut Vec<hc_preproc::CigarElem>, fc: &mut M2FilterCounts, inp: &M2Input, refseq: &RefSeq, dry: bool) -> Item {
        bam::parse_record(body, rec);
        let aux = aux_tags(body);
        let key = MergeKey::of(rec);
        if !passes_m2_filters(rec, &aux, inp, fc) {
            // (filtered reads never reach the downsampler, so they may span chunks)
            return (None, key);
        }
        let sample = *aux.rg.and_then(|g| inp.rg_sample.get(g)).unwrap_or_else(|| panic!("read {} has no known read group", String::from_utf8_lossy(&rec.name)));
        if dry {
            // only the start after palindrome clipping and the MQ matter; clipping needs a properly
            // paired read with a soft clip / insertion at its fragment-start end
            let mut key = key;
            key.passed = true;
            let up = rec.tlen > 0;
            let end_op = if up { rec.cigar.first() } else { rec.cigar.last() };
            if rec.flag & 0x3 == 0x3 && rec.tlen != 0 && end_op.is_some_and(|e| e.0 == OP_S || e.0 == OP_I) {
                let rend = rec.start + ref_length(&rec.cigar) - 1;
                let mut rr = ReadRec::from_rec(rec, rend);
                rr.sample = sample;
                palindrome_clip(&mut rr, refseq, scratch);
                key.eff_start = rr.start;
            }
            return (None, key);
        }
        let rend = rec.start + ref_length(&rec.cigar) - 1;
        let mut rr = ReadRec::from_rec(rec, rend);
        rr.sample = sample;
        // OriginalAlignment: AddOriginalAlignmentTags.getOAContig (OA up to the first ',') vs the current contig
        if inp.keep_records {
            rr.bam_rec = Some(match &inp.merged {
                Some((m, i)) => {
                    let mut b = body.to_vec();
                    crate::markdup::qname::translate_record(m, *i, &mut b);
                    Arc::from(b)
                }
                None => Arc::from(body),
            });
        }
        rr.oa_other_contig = aux.oa.is_some_and(|oa| oa.split(|&c| c == b',').next().unwrap_or(b"") != inp.ref_names[rec.ref_id as usize].as_bytes());
        if palindrome_clip(&mut rr, refseq, scratch) {
            fc.palindrome_clipped += 1;
        }
        (Some(Arc::new(rr)), key)
}

/// Raw BAM records (the bytes after `block_size`), `offs` = record boundaries.
struct RawBatch {
    data: Vec<u8>,
    offs: Vec<usize>,
}

/// The decoders' result batches of one input, read back in round-robin (= record) order.
struct OrderedRx {
    rxs: Vec<std::sync::mpsc::Receiver<Vec<Item>>>,
    next: usize,
}

impl OrderedRx {
    fn recv(&mut self) -> Result<Vec<Item>, std::sync::mpsc::RecvError> {
        let k = self.next % self.rxs.len();
        let r = self.rxs[k].recv();
        if r.is_ok() {
            self.next += 1;
        }
        r
    }
}


/// htsjdk `SAMRecordCoordinateComparator` keys of a raw record (before read filters and palindrome clipping:
/// GATK merges the BAMs first, so reads that are filtered out still decide the interleaving of the inputs).
pub(crate) struct MergeKey {
    start: i64,
    rev: bool,
    name: Box<[u8]>,
    flag: u16,
    mapq: u8,
    mate_ref_id: i32,
    mate_start: i64,
    tlen: i32,
    /// dry run: passes the read filters, and the alignment start after palindrome clipping
    passed: bool,
    eff_start: i64,
}

impl MergeKey {
    pub(crate) fn of(rec: &Rec) -> MergeKey {
        MergeKey {
            start: rec.start,
            rev: rec.flag & 0x10 != 0,
            name: rec.name.clone().into_boxed_slice(),
            flag: rec.flag,
            mapq: rec.mapq,
            mate_ref_id: rec.next_ref_id,
            mate_start: rec.mate_start,
            tlen: rec.tlen,
            passed: false,
            eff_start: rec.start,
        }
    }
}

/// One merge candidate: the read (None if a read filter rejects it) and its comparator keys.
type Item = (Option<Arc<ReadRec>>, MergeKey);

fn coord_cmp(a: &Item, b: &Item) -> std::cmp::Ordering {
    key_cmp(&a.1, &b.1)
}

/// htsjdk `SAMRecordCoordinateComparator` on records of one contig.
pub(crate) fn key_cmp(x: &MergeKey, y: &MergeKey) -> std::cmp::Ordering {
    x.start.cmp(&y.start).then_with(|| {
        if x.rev != y.rev {
            return if x.rev { std::cmp::Ordering::Greater } else { std::cmp::Ordering::Less };
        }
        x.name
            .cmp(&y.name)
            .then(x.flag.cmp(&y.flag))
            .then(x.mapq.cmp(&y.mapq))
            .then(x.mate_ref_id.cmp(&y.mate_ref_id))
            .then(x.mate_start.cmp(&y.mate_start))
            .then(x.tlen.cmp(&y.tlen))
    })
}

/// Runs the Mutect2 front end of one contig (`tids[f]`: the contig's index in input `f`).
#[allow(clippy::too_many_arguments)]
pub fn run<S: M2Sink>(
    inputs: &[M2Input],
    tids: &[usize],
    is_tumor: &[bool],
    contig: &str,
    intervals: &[(i64, i64)],
    refseq: &RefSeq,
    tables: &Tables,
    math: &M2Math,
    max_reads_per_start: usize,
    rng: &mut JavaRandom,
    dump: Option<&mut BufWriter<File>>,
    sink: S,
    f1r2: Option<&mut super::orientation::F1R2Counts>,
    pon: Option<&super::germline::Pon>,
    germline: Option<&super::germline::Germline>,
    opts: super::M2Opts,
    dry: bool,
    decode_threads: usize,
    two_stage: bool,
) -> (M2FrontResult, S) {
    let contig_len = refseq.contig_len;
    let padded = crate::intervals::pad_merge(intervals, REGION_PADDING, contig_len);
    let has_normal = is_tumor.iter().any(|&t| !t);
    // read stage <-> evaluation stage (see the protocol below)
    let ring = Ring::new();
    let frontier = std::sync::atomic::AtomicI64::new(intervals.first().map_or(0, |i| i.0 - 1));
    let scattered = AtomicU64::new(0);
    let want_f1r2 = f1r2.is_some();
    let mut f = Front {
        t: tables,
        math,
        is_tumor,
        has_normal,
        contig,
        intervals,
        refseq,
        sink,
        dump,
        res: M2FrontResult::default(),
        ring: &ring,
        f1r2_delta: vec![Vec::new(); WINDOW],
        frontier: &frontier,
        fin_iv: 0,
        next_final: intervals.first().map_or(1, |i| i.0),
        profile: Profile { probs: VecDeque::new(), start: 0, stop: 0, contig_len, kernel: tables.kernel.clone(), fsize: tables.filter_size },
        regions: Vec::new(),
        filler: RegionFiller::new(contig_len),
        alt_scratch: Vec::new(),
        filled: VecDeque::new(),
        polled: 0,
        applied: 0,
        f1r2,
        evaluated_upto: 0,
        in_fin: false,
        pon,
        germline,
        opts,
    };
    let _ = f.t;
    let mut ds = if max_reads_per_start > 0 {
        Some(MutectDownsampler { max: max_reads_per_start, pending: Vec::new(), stride_start: None, draws_groups: 0, sizes: None })
    } else {
        None
    };
    let mut dry_groups = (dry && max_reads_per_start > 0).then(|| DryGroups { max: max_reads_per_start, stride_start: None, n: 0, well: 0, sizes: Vec::new() });
    let mut fcs = vec![M2FilterCounts::default(); inputs.len()];
    std::thread::scope(|scope| {
        // per input: a splitter thread (records -> htsjdk interval filter) hands batches of raw records round-robin to
        // `decode_threads` workers (filters -> read records -> palindrome clip); results are read back in the same
        // round-robin order, so the record stream is unchanged
        const RAW_BATCH: usize = 2048;
        let k_dec = decode_threads.max(1);
        let mut rxs: Vec<OrderedRx> = Vec::new();
        let mut handles = Vec::new();
        for (fi, inp) in inputs.iter().enumerate() {
            let tid = tids[fi];
            let padded = &padded;
            // records of `tid` overlapping the padded intervals, in file order (htsjdk query semantics)
            let read_records = move |mut each: Box<dyn FnMut(&[u8]) -> bool + '_>| {
                let inflate_ns = Arc::new(AtomicU64::new(0));
                let mut cs = inp.input.stream(tid, padded, false, inflate_ns.clone());
                let mut body = Vec::new();
                let mut q_idx = 0usize;
                'reads: while cs.next_record(&mut body) {
                    let (rid, rstart, _) = bam::peek_core(&body);
                    if rid != tid as i32 {
                        if rid > tid as i32 {
                            break;
                        }
                        continue;
                    }
                    loop {
                        if q_idx >= padded.len() {
                            break 'reads;
                        }
                        if padded[q_idx].1 < rstart {
                            q_idx += 1;
                            continue;
                        }
                        break;
                    }
                    if bam::query_end(&body) < padded[q_idx].0 {
                        continue;
                    }
                    if !each(&body) {
                        break;
                    }
                }
                FT[0].fetch_add(inflate_ns.load(std::sync::atomic::Ordering::Relaxed), std::sync::atomic::Ordering::Relaxed);
            };
            if k_dec == 1 {
                // one thread reads and decodes
                let (rtx, rrx) = std::sync::mpsc::sync_channel::<Vec<Item>>(8);
                rxs.push(OrderedRx { rxs: vec![rrx], next: 0 });
                handles.push((fi, scope.spawn(move || {
                    let (mut rec, mut scratch, mut fc, mut dec_ns) = (Rec::default(), Vec::new(), M2FilterCounts::default(), 0u64);
                    let mut batch: Vec<Item> = Vec::with_capacity(4096);
                    read_records(Box::new(|body: &[u8]| {
                        let t_dec = ft_now();
                        batch.push(decode_record(body, &mut rec, &mut scratch, &mut fc, inp, refseq, dry));
                        if let Some(t) = t_dec {
                            dec_ns += t.elapsed().as_nanos() as u64;
                        }
                        batch.len() < 4096 || rtx.send(std::mem::replace(&mut batch, Vec::with_capacity(4096))).is_ok()
                    }));
                    if !batch.is_empty() {
                        let _ = rtx.send(batch);
                    }
                    FT[1].fetch_add(dec_ns, std::sync::atomic::Ordering::Relaxed);
                    fc
                })));
                continue;
            }
            let mut work_txs = Vec::new();
            let mut res_rxs = Vec::new();
            for _ in 0..k_dec {
                let (wtx, wrx) = std::sync::mpsc::sync_channel::<RawBatch>(2);
                let (rtx, rrx) = std::sync::mpsc::sync_channel::<Vec<Item>>(2);
                work_txs.push(wtx);
                res_rxs.push(rrx);
                handles.push((fi, scope.spawn(move || {
                    let (mut rec, mut scratch, mut fc, mut dec_ns) = (Rec::default(), Vec::new(), M2FilterCounts::default(), 0u64);
                    while let Ok(raw) = wrx.recv() {
                        let t_dec = ft_now();
                        let batch: Vec<Item> =
                            raw.offs.windows(2).map(|w| decode_record(&raw.data[w[0]..w[1]], &mut rec, &mut scratch, &mut fc, inp, refseq, dry)).collect();
                        if let Some(t) = t_dec {
                            dec_ns += t.elapsed().as_nanos() as u64;
                        }
                        if rtx.send(batch).is_err() {
                            break;
                        }
                    }
                    FT[1].fetch_add(dec_ns, std::sync::atomic::Ordering::Relaxed);
                    fc
                })));
            }
            rxs.push(OrderedRx { rxs: res_rxs, next: 0 });
            scope.spawn(move || {
                let mut seq = 0usize;
                let mut raw = RawBatch { data: Vec::new(), offs: vec![0] };
                let mut ok = true;
                read_records(Box::new(|body: &[u8]| {
                    raw.data.extend_from_slice(body);
                    raw.offs.push(raw.data.len());
                    if raw.offs.len() > RAW_BATCH {
                        let full = std::mem::replace(&mut raw, RawBatch { data: Vec::new(), offs: vec![0] });
                        ok = work_txs[seq % k_dec].send(full).is_ok();
                        seq += 1;
                    }
                    ok
                }));
                if ok && raw.offs.len() > 1 {
                    let _ = work_txs[seq % k_dec].send(raw);
                }
            });
        }
        // k-way merge (htsjdk MergingSamRecordIterator); the closure gets each record in merged order
        let merge = |rxs: &mut Vec<OrderedRx>, each: &mut dyn FnMut(Option<Arc<ReadRec>>, MergeKey)| {
            let mut heads: Vec<VecDeque<Item>> = (0..rxs.len()).map(|_| VecDeque::new()).collect();
            let mut done = vec![false; rxs.len()];
            loop {
                for k in 0..rxs.len() {
                    if heads[k].is_empty() && !done[k] {
                        let t_w = ft_now();
                        let got = rxs[k].recv();
                        ft_add(2, t_w);
                        match got {
                            Ok(b) => heads[k].extend(b),
                            Err(_) => done[k] = true,
                        }
                    }
                }
                let mut best: Option<usize> = None;
                for k in 0..rxs.len() {
                    if let Some(h) = heads[k].front() {
                        if best.is_none_or(|b| coord_cmp(h, heads[b].front().unwrap()) == std::cmp::Ordering::Less) {
                            best = Some(k);
                        }
                    }
                }
                let Some(b) = best else { break };
                let (r, key) = heads[b].pop_front().unwrap();
                each(r, key);
            }
        };
        if let Some(g) = dry_groups.as_mut() {
            merge(&mut rxs, &mut |_, key| {
                if key.passed {
                    g.submit(key.eff_start, key.mapq);
                }
            });
            g.finalize();
            f.res.group_sizes = std::mem::take(&mut g.sizes);
        } else if !two_stage {
            // single thread: per read, the evaluation stage's work, then its scatter (the order of the two-stage
            // protocol below, run serially)
            let t_loop = ft_now();
            let mut sc = Scatter { ring: &ring, intervals, is_tumor, refseq, want_f1r2, scatter_iv: 0, local_ivs: Vec::new() };
            let mut out: Vec<Arc<ReadRec>> = Vec::new();
            let mut ds = ds.take();
            let (mut reads_in, mut kept) = (0u64, 0u64);
            let handle = |r: Arc<ReadRec>, f: &mut Front<S>, sc: &mut Scatter, kept: &mut u64| {
                if f.ring.too_short_for(r.end - r.start + 1) {
                    // SAFETY: single thread, between reads
                    unsafe { f.grow(r.end - r.start + 1) };
                }
                *kept += 1;
                f.process(r.clone());
                let t = ft_now();
                sc.scatter(&r);
                ft_add(6, t);
            };
            merge(&mut rxs, &mut |r, _| {
                let Some(r) = r else { return };
                reads_in += 1;
                let t_ds = ft_now();
                match ds.as_mut() {
                    None => out.push(r),
                    Some(d) => d.submit(r, rng, &mut out),
                }
                ft_add(4, t_ds);
                for r in out.drain(..) {
                    handle(r, &mut f, &mut sc, &mut kept);
                }
            });
            if let Some(d) = ds.as_mut() {
                d.finalize(rng, &mut out);
                f.res.drawing_groups = d.draws_groups;
            }
            for r in out.drain(..) {
                handle(r, &mut f, &mut sc, &mut kept);
            }
            f.res.reads_in = reads_in;
            f.res.reads_kept = kept;
            ft_add(3, t_loop);
        } else {
            // Read stage (own thread): merge, downsample, scatter into the locus accumulators. Evaluation stage (this
            // thread): per read, in order, region filling, evaluation of the loci before its start, regions (overlap
            // cleaning, hand-off), as the single-threaded front end did. Protocol: the read stage sends a batch of
            // reads, then scatters them one by one, counting them in `scattered`; before scattering a read it waits
            // until its last locus is < frontier + 1 + window (a read spanning half the window or more first grows it:
            // `Batch::Grow`, at a point where neither stage uses the ring). The evaluation stage processes read n only when
            // `scattered` >= n (every earlier read is in the accumulators, so all loci before its start are complete)
            // and then publishes frontier = its start - 1 (slots up to there are evaluated and reset). The read stage
            // only writes slots of loci >= the start of the read it scatters, i.e. above the frontier of every read the
            // evaluation stage can be processing, and within one window of it: the two never share a slot.
            let (btx, brx) = std::sync::mpsc::sync_channel::<Batch>(64);
            let (gtx, grx) = std::sync::mpsc::sync_channel::<()>(1);
            let (ring, frontier, scattered) = (&ring, &frontier, &scattered);
            let mut rxs_s = std::mem::take(&mut rxs);
            let rng = &mut *rng;
            let mut ds = ds.take();
            let reader = scope.spawn(move || {
                let mut sc = Scatter { ring, intervals, is_tumor, refseq, want_f1r2, scatter_iv: 0, local_ivs: Vec::new() };
                let (mut reads_in, mut groups) = (0u64, 0u64);
                let mut out: Vec<Arc<ReadRec>> = Vec::new();
                let mut batch: Vec<Arc<ReadRec>> = Vec::with_capacity(256);
                let mut flush = |batch: &mut Vec<Arc<ReadRec>>| -> bool {
                    let b = std::mem::take(batch);
                    let keep: Vec<Arc<ReadRec>> = b.iter().cloned().collect();
                    if btx.send(Batch::Reads(b)).is_err() {
                        return false;
                    }
                    let w = ring.window() as i64;
                    for r in &keep {
                        let mut spins = 0u32;
                        while r.end >= frontier.load(std::sync::atomic::Ordering::Acquire) + 1 + w {
                            spins += 1;
                            if spins < 64 { std::hint::spin_loop() } else { std::thread::yield_now() }
                        }
                        let t = ft_now();
                        sc.scatter(r);
                        ft_add(6, t);
                        scattered.fetch_add(1, std::sync::atomic::Ordering::Release);
                    }
                    true
                };
                // kept reads into the batch; a read too long for the ring: everything before it is sent and
                // scattered, then the evaluation stage grows the ring while this stage waits
                let deliver = |batch: &mut Vec<Arc<ReadRec>>, out: &mut Vec<Arc<ReadRec>>, flush: &mut dyn FnMut(&mut Vec<Arc<ReadRec>>) -> bool| -> bool {
                    for r in out.drain(..) {
                        let span = r.end - r.start + 1;
                        if ring.too_short_for(span) {
                            if !batch.is_empty() && !flush(batch) {
                                return false;
                            }
                            if btx.send(Batch::Grow(span)).is_err() || grx.recv().is_err() {
                                return false;
                            }
                        }
                        batch.push(r);
                    }
                    batch.len() < 256 || flush(batch)
                };
                let mut ok = true;
                merge(&mut rxs_s, &mut |r, _| {
                    let Some(r) = r else { return };
                    if !ok {
                        return;
                    }
                    reads_in += 1;
                    let t_ds = ft_now();
                    match ds.as_mut() {
                        None => out.push(r),
                        Some(d) => d.submit(r, rng, &mut out),
                    }
                    ft_add(4, t_ds);
                    ok = deliver(&mut batch, &mut out, &mut flush);
                });
                if let Some(d) = ds.as_mut() {
                    d.finalize(rng, &mut out);
                    groups = d.draws_groups;
                }
                if ok && deliver(&mut batch, &mut out, &mut flush) && !batch.is_empty() {
                    flush(&mut batch);
                }
                drop(rxs_s);
                (reads_in, groups)
            });
            let t_loop = ft_now();
            let mut n = 0u64;
            for b in brx.iter() {
                let b = match b {
                    Batch::Reads(b) => b,
                    Batch::Grow(span) => {
                        // every earlier read is processed here and scattered by the read stage, which now waits
                        // SAFETY: neither stage uses the ring until the reply
                        unsafe { f.grow(span) };
                        if gtx.send(()).is_err() {
                            break;
                        }
                        continue;
                    }
                };
                for r in b {
                    let mut spins = 0u32;
                    while scattered.load(std::sync::atomic::Ordering::Acquire) < n {
                        spins += 1;
                        if spins < 64 { std::hint::spin_loop() } else { std::thread::yield_now() }
                    }
                    f.res.reads_kept += 1;
                    f.process(r);
                    n += 1;
                }
            }
            let (reads_in, groups) = reader.join().expect("read stage failed");
            f.res.reads_in = reads_in;
            f.res.drawing_groups = groups;
            ft_add(3, t_loop);
        }
        drop(rxs);
        for (k, h) in handles {
            fcs[k].add(&h.join().expect("BAM decoding failed"));
        }
    });
    for fc in &fcs {
        f.res.filters.add(fc);
    }
    if dry {
        return (f.res, f.sink);
    }
    let t_end = ft_now();
    f.finalize(i64::MAX);
    if !f.profile.is_empty() {
        f.profile.pop_ready(true, &mut f.regions);
    }
    f.emit(i64::MAX, true);
    while f.applied < f.regions.len() {
        f.apply(f.applied);
    }
    ft_add(12, t_end);
    f.res.regions = std::mem::take(&mut f.regions);
    (f.res, f.sink)
}

/// `ReadUtils.getReadIndexForReferenceCoordinate(read, refCoord)` (soft clips aligned).
fn read_index(r: &ReadRec, cig: &[hc_preproc::CigarElem], ref_coord: i64) -> (i64, Option<u8>) {
    let mut soft_start = r.start;
    for &(op, len) in cig {
        if op == b'S' {
            soft_start -= len as i64;
        } else if op != b'H' {
            break;
        }
    }
    hc_preproc::read_index_for_ref(soft_start, cig, ref_coord)
}

fn char_cigar(r: &ReadRec) -> Vec<hc_preproc::CigarElem> {
    r.cigar.iter().map(|&(op, l)| (b"MIDNSHP=X"[op as usize], l)).collect()
}

/// `FragmentUtils.adjustQualsOfOverlappingPairedFragments(pair, setConflictingToZero = false, 20, 20)`
/// applied to the reads' shared quality state.
/// `capped`: (reads ending after this position, list) receiving, per read, the read indices whose quality was capped
/// from above 20 to 20.
fn adjust_pair(left: &Arc<ReadRec>, right: &Arc<ReadRec>, capped: Option<(i64, &mut Vec<(Arc<ReadRec>, Vec<usize>)>)>) {
    const HALF: u8 = PCR_SNV_QUAL / 2; // pcrSnvQual / 2 = pcrIndelQual / 2 = 20
    let soft = |r: &ReadRec| {
        let mut s = r.start;
        for &(op, len) in &r.cigar {
            if op == OP_S {
                s -= len as i64;
            } else if op != OP_H {
                break;
            }
        }
        s
    };
    let in_order = soft(left) < soft(right);
    let (first, second) = if in_order { (left, right) } else { (right, left) };
    if first.end < second.start {
        return;
    }
    // character cigars on the stack (heap only for unusually long cigars)
    let mut bf = [(0u8, 0u32); 32];
    let mut bs = [(0u8, 0u32); 32];
    let (vf, vs);
    let cf: &[hc_preproc::CigarElem] = if first.cigar.len() <= 32 {
        for (k, &(op, l)) in first.cigar.iter().enumerate() {
            bf[k] = (b"MIDNSHP=X"[op as usize], l);
        }
        &bf[..first.cigar.len()]
    } else {
        vf = char_cigar(first);
        &vf
    };
    let cs: &[hc_preproc::CigarElem] = if second.cigar.len() <= 32 {
        for (k, &(op, l)) in second.cigar.iter().enumerate() {
            bs[k] = (b"MIDNSHP=X"[op as usize], l);
        }
        &bs[..second.cigar.len()]
    } else {
        vs = char_cigar(second);
        &vs
    };
    let (offset, op) = read_index(first, cf, second.start);
    if offset == hc_preproc::READ_INDEX_NOT_FOUND || matches!(op, Some(b'S') | Some(b'H')) {
        return;
    }
    let first_end_base = read_index(first, cf, first.end).0;
    let second_end_base = read_index(second, cs, second.end).0;
    let second_offset = read_index(second, cs, second.start).0;
    let n = (first_end_base - offset).min(second_end_base - second_offset) + 1;
    if n <= 0 {
        return;
    }
    // re-cleaning a pair (it is in every region overlapping it) usually changes nothing: check before copying
    {
        let (gf, gs) = (first.m2q.lock().unwrap(), second.m2q.lock().unwrap());
        let (Some(sf), Some(ss)) = (gf.as_ref(), gs.as_ref()) else {
            // no state yet: BI / BD start at 45 and are capped
            drop((gf, gs));
            return adjust_pair_apply(first, second, offset, second_offset, n, capped);
        };
        let mut any = false;
        for i in 0..n {
            let (a, b) = ((offset + i) as usize, (second_offset + i) as usize);
            if (first.seq[a] == second.seq[b] && (sf.qual[a] > HALF || ss.qual[b] > HALF)) || sf.bi[a] > HALF || sf.bd[a] > HALF || ss.bi[b] > HALF || ss.bd[b] > HALF {
                any = true;
                break;
            }
        }
        if !any {
            return;
        }
    }
    adjust_pair_apply(first, second, offset, second_offset, n, capped)
}

/// The capping part of [`adjust_pair`] (overlap of `n` bases from read index `offset` of `first` and
/// `second_offset` of `second`).
fn adjust_pair_apply(first: &Arc<ReadRec>, second: &Arc<ReadRec>, offset: i64, second_offset: i64, n: i64, capped: Option<(i64, &mut Vec<(Arc<ReadRec>, Vec<usize>)>)>) {
    const HALF: u8 = PCR_SNV_QUAL / 2;
    let state = |r: &ReadRec| -> M2Quals {
        match r.m2q.lock().unwrap().as_ref() {
            Some(s) => (**s).clone(),
            None => M2Quals { qual: r.qual.clone(), bi: vec![45; r.seq.len()], bd: vec![45; r.seq.len()] },
        }
    };
    let (mut fs, mut ss) = (state(first), state(second));
    let n = n as usize;
    let (o1, o2) = (offset as usize, second_offset as usize);
    let (fb, sb) = (&first.seq[o1..o1 + n], &second.seq[o2..o2 + n]);
    // read indices whose base quality gets capped (for the F1R2 counts), before capping
    let want = capped.as_ref().is_some_and(|(min_end, _)| first.end > *min_end || second.end > *min_end);
    let (mut caps_f, mut caps_s) = (Vec::new(), Vec::new());
    if want {
        for i in 0..n {
            if fb[i] == sb[i] {
                if fs.qual[o1 + i] > HALF {
                    caps_f.push(o1 + i);
                }
                if ss.qual[o2 + i] > HALF {
                    caps_s.push(o2 + i);
                }
            }
        }
    }
    // branch-free capping (vectorizes); `changed` = some value was above HALF
    let mut over = 0u8;
    {
        let (fq, sq) = (&mut fs.qual[o1..o1 + n], &mut ss.qual[o2..o2 + n]);
        for i in 0..n {
            let eq = fb[i] == sb[i];
            let (a, b) = (fq[i], sq[i]);
            over |= (eq & (a > HALF)) as u8 | (eq & (b > HALF)) as u8;
            fq[i] = if eq { a.min(HALF) } else { a };
            sq[i] = if eq { b.min(HALF) } else { b };
        }
    }
    for v in [&mut fs.bi[o1..o1 + n], &mut fs.bd[o1..o1 + n], &mut ss.bi[o2..o2 + n], &mut ss.bd[o2..o2 + n]] {
        for x in v.iter_mut() {
            over |= (*x > HALF) as u8;
            *x = (*x).min(HALF);
        }
    }
    // setBaseQualities / BI / BD tags are written even when nothing changed (no observable difference)
    if over != 0 {
        *first.m2q.lock().unwrap() = Some(Arc::new(fs));
        *second.m2q.lock().unwrap() = Some(Arc::new(ss));
        first.cleaned.store(true, std::sync::atomic::Ordering::Relaxed);
        second.cleaned.store(true, std::sync::atomic::Ordering::Relaxed);
        if let Some((min_end, c)) = capped {
            for (r, idx) in [(first, caps_f), (second, caps_s)] {
                if r.end > min_end && !idx.is_empty() {
                    c.push((r.clone(), idx));
                }
            }
        }
    }
}

/// `AssemblyBasedCallerUtils.cleanOverlappingReadPairs(reads, samples, header, false, 20, 20)`: per sample,
/// `FragmentCollection.create` pairing (reads in region order, which is start order).
/// FNV-1a (read names in the pairing map).
#[derive(Default)]
struct Fnv(u64);
impl std::hash::Hasher for Fnv {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        let mut h = if self.0 == 0 { 0xcbf29ce484222325 } else { self.0 };
        for &b in bytes {
            h ^= b as u64;
            h = h.wrapping_mul(0x100000001b3);
        }
        self.0 = h;
    }
}

fn clean_overlapping_read_pairs(reads: &[Arc<ReadRec>], n_samples: usize, want_capped: Option<i64>) -> Vec<(Arc<ReadRec>, Vec<usize>)> {
    let mut capped = Vec::new();
    for s in 0..n_samples {
        let mut name_map: std::collections::HashMap<&[u8], usize, std::hash::BuildHasherDefault<Fnv>> = Default::default();
        let mut pairs: Vec<(usize, usize)> = Vec::new();
        for (i, r) in reads.iter().enumerate() {
            if r.sample as usize != s {
                continue;
            }
            // singletons: unpaired, mate unmapped (mate_start 0), or mate starting after this read's end
            if r.flag & 0x1 == 0 || r.mate_start == 0 || r.mate_start > r.end {
                continue;
            }
            match name_map.remove(&r.name[..]) {
                Some(j) => pairs.push((j, i)),
                None => {
                    name_map.insert(&r.name[..], i);
                }
            }
        }
        use std::sync::atomic::Ordering::Relaxed;
        for (j, i) in pairs {
            let (a, b) = (&reads[j], &reads[i]);
            let (pa, pb) = (Arc::as_ptr(a) as usize, Arc::as_ptr(b) as usize);
            if a.cleaned_with.load(Relaxed) == pb && b.cleaned_with.load(Relaxed) == pa {
                if ft_on() {
                    FTC[0].fetch_add(1, Relaxed);
                }
                continue;
            }
            if ft_on() {
                FTC[1].fetch_add(1, Relaxed);
            }
            let t = ft_now();
            adjust_pair(a, b, want_capped.map(|m| (m, &mut capped)));
            ft_add(15, t);
            a.cleaned_with.store(pb, Relaxed);
            b.cleaned_with.store(pa, Relaxed);
        }
    }
    capped
}
