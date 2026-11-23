//! Per-contig front end: one GATK read shard (`MultiIntervalLocalReadShard`: all intervals of one contig,
//! reads queried over the intervals padded by 100) through `AssemblyRegionIterator`:
//!
//! ```text
//! BAI chunks -> BGZF inflate (reader thread) -> htsjdk interval-overlap filter -> HC read filters (decoder thread)
//!   -> read stage: PositionalDownsampler -> pileup scatter into the locus accumulators
//!   -> evaluation stage: region read cache / isActive / activity profile / region cutting -> region filling
//!   -> batches of FilledRegion to the back end
//! ```
//!
//! The two stages run in one thread, or on two threads under the protocol described in [`run`] (`two_stage`).

use crate::active::*;
use crate::bam::{self, Rec};
use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::atomic::{AtomicBool, AtomicI64, AtomicU64};
use std::sync::Arc;
use hc_downsample::{lookahead_locus, pileup_span, Horizon, JavaRandom, PositionalDownsampler, ReadyAt, ReadyTracker};
use std::cmp::Reverse;
use std::collections::BinaryHeap;

/// A filtered read kept for the rest of the pipeline.
pub struct ReadRec {
    pub name: Vec<u8>,
    pub flag: u16,
    pub start: i64,
    pub end: i64,
    pub mapq: u8,
    pub cigar: Vec<(u8, u32)>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    /// mate start, 0 if unpaired or the mate is unmapped (GATK `getMateStart()`)
    pub mate_start: i64,
    pub tlen: i32,
    pub ref_id: i32,
    /// BAM reference index of the mate (RNEXT); -1 if none
    pub mate_ref_id: i32,
    /// sample index (Mutect2; 0 for HaplotypeCaller)
    pub sample: u16,
    /// Mutect2 `OriginalAlignment`: the read's OA tag names another contig
    pub oa_other_contig: bool,
    /// the original BAM record (bytes after block_size; Mutect2 `--bam-output` only)
    pub bam_rec: Option<Arc<[u8]>>,
    /// Mutect2: base / BI / BD qualities as changed in place by `cleanOverlappingReadPairs` (None = unchanged)
    pub m2q: std::sync::Mutex<Option<Arc<M2Quals>>>,
    /// Mutect2: `m2q` has been set (lets readers skip the lock)
    pub cleaned: std::sync::atomic::AtomicBool,
    /// Mutect2 overlap cleaning: address of the read this one was last adjusted with (0 = none); adjusting the same
    /// pair again cannot change anything
    pub cleaned_with: std::sync::atomic::AtomicUsize,
}

/// A read's qualities after Mutect2's in-place overlap cleaning (BI/BD start at the default Q45).
#[derive(Clone, Debug)]
pub struct M2Quals {
    pub qual: Vec<u8>,
    pub bi: Vec<u8>,
    pub bd: Vec<u8>,
}

impl ReadRec {
    pub(crate) fn from_rec(rec: &Rec, end: i64) -> ReadRec {
        let paired_mapped_mate = rec.flag & 0x1 != 0 && rec.flag & 0x8 == 0 && rec.next_ref_id >= 0;
        ReadRec {
            name: rec.name.clone(),
            flag: rec.flag,
            start: rec.start,
            end,
            mapq: rec.mapq,
            cigar: rec.cigar.clone(),
            seq: rec.seq.clone(),
            qual: rec.qual.clone(),
            mate_start: if paired_mapped_mate { rec.mate_start } else { 0 },
            tlen: rec.tlen,
            ref_id: rec.ref_id,
            mate_ref_id: rec.next_ref_id,
            sample: 0,
            oa_other_contig: false,
            bam_rec: None,
            m2q: std::sync::Mutex::new(None),
            cleaned: std::sync::atomic::AtomicBool::new(false),
            cleaned_with: std::sync::atomic::AtomicUsize::new(0),
        }
    }

    /// The read as `hc_preproc` takes it.
    pub fn to_preproc(&self) -> hc_preproc::Read {
        hc_preproc::Read {
            name: String::from_utf8_lossy(&self.name).into_owned(),
            flag: self.flag,
            start: self.start,
            end: self.end,
            cigar: self.cigar.iter().map(|&(op, l)| (b"MIDNSHP=X"[op as usize], l)).collect(),
            seq: self.seq.clone(),
            qual: self.qual.clone(),
            mapq: self.mapq,
            mate_start: self.mate_start,
            tlen: self.tlen as i64,
            ref_id: self.ref_id,
            mate_ref_id: self.mate_ref_id,
            bi: Vec::new(),
            bd: Vec::new(),
            sample: self.sample,
            oa_other_contig: self.oa_other_contig,
            bam_rec: self.bam_rec.clone(),
        }
    }

    /// `ReadUtils.getAdaptorBoundary` + `isBaseInsideAdaptor` (applied by `LocusIteratorByState`): Some((reverse,
    /// boundary)) if bases beyond the boundary are excluded from pileups. GATK compares the *signed* fragment
    /// length with 100, so the size cutoff never applies to reverse-strand reads. Mates on another contig:
    /// `hasWellDefinedFragmentSize` has no contig test, only `TLEN == 0` (what aligners write for such pairs);
    /// with a non-zero TLEN GATK uses the mate position on the other contig, and so does this port.
    pub(crate) fn adaptor_exclusion(&self) -> Option<(bool, i64)> {
        const DEFAULT_ADAPTOR_SIZE: i32 = 100;
        let paired = self.flag & 0x1 != 0;
        let rev = self.flag & 0x10 != 0;
        let mate_rev = self.flag & 0x20 != 0;
        if self.tlen == 0 || !paired {
            return None;
        }
        // mateIsUnmapped(): flag, RNEXT '*' or mate position 0 (all folded into mate_start == 0)
        if self.flag & 0x4 != 0 || self.flag & 0x8 != 0 || self.mate_start == 0 {
            return None;
        }
        if rev == mate_rev {
            return None;
        }
        let well = if rev { self.end > self.mate_start } else { self.start <= self.mate_start + self.tlen as i64 };
        if !well {
            return None;
        }
        if self.tlen > DEFAULT_ADAPTOR_SIZE {
            return None;
        }
        Some((rev, if rev { self.mate_start - 1 } else { self.start + (self.tlen as i64).abs() }))
    }
}

pub struct FilledRegion {
    /// index of the region in the contig (apply order)
    pub id: usize,
    pub start: i64,
    pub end: i64,
    pub padded_start: i64,
    pub padded_end: i64,
    pub active: bool,
    pub reads: Vec<Arc<ReadRec>>,
}

pub const REGION_PADDING: i64 = 100;

/// GATK semantics: a region gets (1) the previous region's reads that overlap its padded span, in that
/// order, then (2) reads from the cache, in file order, until a read starts after the padded span.
pub(crate) struct RegionFiller {
    pub(crate) cache: VecDeque<Arc<ReadRec>>,
    prev: Vec<Arc<ReadRec>>,
    next_region: usize,
    contig_len: i64,
}

impl RegionFiller {
    pub(crate) fn new(contig_len: i64) -> Self {
        RegionFiller { cache: VecDeque::new(), prev: Vec::new(), next_region: 0, contig_len }
    }

    pub(crate) fn drain(&mut self, regions: &[Region], stream_pos: i64, eof: bool, sink: &mut impl FnMut(FilledRegion)) {
        while self.next_region < regions.len() {
            let r = regions[self.next_region];
            let (ps, pe) = ((r.start - REGION_PADDING).max(1), (r.end + REGION_PADDING).min(self.contig_len));
            if !eof && stream_pos <= pe {
                return;
            }
            let mut reads: Vec<Arc<ReadRec>> = Vec::new();
            for rr in &self.prev {
                if rr.start <= pe && rr.end >= ps {
                    reads.push(rr.clone());
                }
            }
            while let Some(front) = self.cache.front() {
                if front.start > pe {
                    break;
                }
                let rr = self.cache.pop_front().unwrap();
                if rr.end >= ps {
                    reads.push(rr);
                }
            }
            self.prev = reads.clone();
            sink(FilledRegion { id: self.next_region, start: r.start, end: r.end, padded_start: ps, padded_end: pe, active: r.active, reads });
            self.next_region += 1;
        }
    }
}

/// Reference bases of a window of one contig.
pub struct RefSeq {
    pub from: i64,
    pub seq: Vec<u8>,
    pub contig_len: i64,
}

impl RefSeq {
    #[inline]
    pub fn base(&self, pos: i64) -> u8 {
        self.seq[(pos - self.from) as usize]
    }
    pub fn window(&self, s: i64, e: i64) -> &[u8] {
        &self.seq[(s - self.from) as usize..=(e - self.from) as usize]
    }
}

/// What the front end needs from the run.
pub struct FrontCtx<'a> {
    /// the inputs (`-I`, one sample), with the same sequence dictionary
    pub inputs: &'a [crate::input::Input],
    pub ref_lens: &'a [i64],
    pub tables: &'a Tables,
    pub parallel_inflate: bool,
    pub inflate_ns: Arc<AtomicU64>,
    /// ApplyBQSR on the fly (`--bqsr-recal-file`)
    pub bqsr: Option<Arc<crate::bqsr::apply::Recalibrator>>,
    /// several samples: per input, read group -> sample index (into the sorted sample names); None with one
    /// sample (every read is sample 0)
    pub rg_sample: Option<&'a [std::collections::HashMap<Vec<u8>, u16>]>,
    /// the sample names, sorted (`ReadUtils.getSamplesFromHeader`)
    pub samples: &'a [String],
}

/// `hc_downsample::LibsCoverage` with the same answers in O(log n): which loci have a non-empty
/// LocusIteratorByState pileup. Queries come with non-decreasing `from`, and a span added after a query
/// starts after that query's `from` (its read starts after the query's `limit`), so spans that start at or
/// before `from` can be folded into a running maximum end.
#[derive(Default)]
struct Coverage {
    heap: BinaryHeap<Reverse<(i64, i64)>>,
    max_end: i64,
    next_prune: usize,
}

impl Coverage {
    /// `lower_bound`: no future query has `from` below it.
    fn add(&mut self, span: Option<(i64, i64)>, lower_bound: i64) {
        let Some((s, e)) = span else { return };
        if e < lower_bound {
            return;
        }
        self.heap.push(Reverse((s, e)));
        if self.heap.len() > self.next_prune {
            self.heap.retain(|Reverse((_, e))| *e >= lower_bound);
            self.next_prune = 2 * self.heap.len() + (1 << 16);
        }
    }
    /// Smallest covered locus >= `from` among loci <= `limit`.
    fn first_covered(&mut self, from: i64, limit: i64) -> Option<i64> {
        while let Some(&Reverse((s, e))) = self.heap.peek() {
            if s > from {
                break;
            }
            self.max_end = self.max_end.max(e);
            self.heap.pop();
        }
        let best = if self.max_end >= from { from } else { self.heap.peek()?.0 .0 };
        if best <= limit { Some(best) } else { None }
    }
}

/// Where the front end sends its output, and where it gets GATK's global random generator from.
pub trait FrontSink {
    /// A filled region, in apply order.
    fn region(&mut self, fr: FilledRegion);
    /// The generator for the draws of a drawing group (a start position with more than `target` reads).
    /// Before it returns, every QD jitter draw of the regions with id < `resolved_before` (the regions that
    /// GATK applies before this group's draws) must have been made.
    /// `None`: this segment was absorbed by the previous one (its results are void; the front end stops).
    fn rng_for_draw(&mut self, resolved_before: usize) -> Option<&mut JavaRandom>;
    /// Segments: this segment takes over the range of the next one (`rng::DomainRng::absorb`).
    fn absorb_next(&mut self) {}
    /// Segments: has the previous segment taken over this one's range?
    fn is_absorbed(&self) -> bool {
        false
    }
}

/// A segment of a contig inside one RNG stream (single-run mode, VCF): the front end reads from `own_from` -
/// `SEG_WARMUP` and owns the active regions starting in [b(own_from), b(cuts[0])), where b(c) is the first locus
/// x >= c whose isActive value is 0 at every locus of [x - f, x + f] (f: the band-pass filter size). Splitting a
/// contig there changes nothing:
///
/// * The isActive value of a locus depends only on the kept reads covering it. The smoothed value at x only receives
///   the spread of loci within f of x (the HQ soft-clip states too, as GATK centres them on the added locus), and
///   none of them adds anything: it is 0, x is inactive in every run, and no active region contains it. Every
///   popped value is final (a region pops with >= 350 + f loci of lookahead).
/// * After b, the activity status of every locus is the same in both runs (same kept reads), so the first active
///   stretch after b starts at the same locus: an inactive region always ends at a status change, whatever the
///   phase of the inactive chunks before it, and regions are then cut from the probabilities alone. Before b, the
///   previous segment is the continuous run. So both segments together give the whole-contig active regions.
///   Inactive regions give no VCF line and no draw (VCF mode only).
/// * Kept reads: no drawing group may lie in [cut - `SEG_ZONE`, stop of the previous segment], nor have a read
///   ending at or after cut - `SEG_ZONE`. The previous segment, which sees those groups (it is the continuous run
///   there), takes the next range over when one does (`absorb`), so this segment's results are then void. Earlier
///   drawing groups of the warm-up (reads ending before the zone) are drawn from a private generator: they change
///   nothing at or after the zone.
/// * Draw order: the previous segment stops only once the read-pull horizon X of its last owned region (GATK draws
///   a group before that region's QD jitter iff the group starts <= X) is below the current read start s, and every
///   group up to s has been seen (no absorption). So all its QD draws precede this segment's first group draw, as in
///   GATK, and the generator is handed on in segment order.
#[derive(Clone)]
pub struct SegSpec {
    /// None: from the start of the contig (no warm-up)
    pub own_from: Option<i64>,
    /// the cuts after this segment's range, in order (the first ends it; later ones in case it absorbs)
    pub cuts: Vec<i64>,
}

/// Reads before a segment's cut that its front end evaluates (warm-up), so that loci >= cut - SEG_ZONE are complete.
pub const SEG_WARMUP: i64 = 2000;
/// A drawing group starting or with a read ending at or after cut - SEG_ZONE (and up to the previous segment's stop)
/// makes the previous segment absorb the next one.
pub const SEG_ZONE: i64 = 1000;
/// The previous segment reads at least this far past the first region it does not own.
pub const SEG_STOP_MIN: i64 = 1000;

struct SegState {
    own_from: Option<i64>,
    cuts: VecDeque<i64>,
    /// boundary loci: the segment owns the active regions starting in [from_b, end_b)
    from_b: Option<i64>,
    end_b: Option<i64>,
    /// first region starting at or after `end_b`
    end_idx: Option<usize>,
    scanned: usize,
    /// consecutive loci (ending at `last_locus`) with isActive 0
    zero_run: i64,
    last_locus: i64,
    /// filled regions starting at or after `end_b`: dropped at the stop, released if the next range is absorbed
    held: Vec<FilledRegion>,
    /// stop rule: lookahead locus N of the first region not owned, an upper bound on its first covered locus C,
    /// distinct kept-read starts seen after it
    n_look: Option<i64>,
    c_bound: Option<i64>,
    starts_after: usize,
    last_start: i64,
    stopped: bool,
    n_proc: u64,
    dummy: JavaRandom,
    absorbed: usize,
}

impl SegState {
    fn reset_stop(&mut self) {
        self.n_look = None;
        self.c_bound = None;
        self.starts_after = 0;
        self.last_start = i64::MIN;
    }
}

/// `VARIANT_RS_TRACE_FRONT=1`: HaplotypeCaller front-end time split (nanoseconds, summed over contigs and threads).
pub static HFT: [std::sync::atomic::AtomicU64; 14] = [const { std::sync::atomic::AtomicU64::new(0) }; 14];
pub const HFT_NAMES: [&str; 14] = [
    "decoder thread: decode + filters (excl. hand-off waits)",
    "front: wait for decoder",
    "front: process total",
    "scatter (read stage)",
    "front: finalize (isActive + profile + region cutting)",
    "front: emit (region filling + hand-off)",
    "front: downsampler advance (incl. nested process)",
    "front: end-of-contig flush",
    "front: whole read loop",
    "front: wait for back end (QD draws before a drawing group)",
    "front: hand-off of region batches to the back end (incl. waiting when it is busy)",
    "two-stage read stage: wait for the evaluation frontier",
    "two-stage evaluation: wait for the read stage's scatter",
    "two-stage read stage: drawing groups (wait for evaluation catch-up, QD order, draws)",
];
pub fn hft_now() -> Option<std::time::Instant> {
    if crate::mutect::front::ft_on() { Some(std::time::Instant::now()) } else { None }
}
pub fn hft_add(i: usize, t: Option<std::time::Instant>) {
    if let Some(t) = t {
        HFT[i].fetch_add(t.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
    }
}

#[derive(Default)]
pub struct FrontResult {
    pub regions: Vec<Region>,
    pub filters: FilterCounts,
    pub n_loci: u64,
    pub n_active_loci: u64,
    pub assigned: u64,
    pub ds: hc_downsample::DownsampleStats,
    /// segments: (start of the first owned region, start of the first region not owned, segments absorbed, stopped)
    pub seg_range: Option<(i64, i64, usize, bool)>,
}

/// Initial size of the ring buffer of per-locus accumulators. It must exceed twice the longest read reference span;
/// a longer read grows it first (`Ring::grow`).
const WINDOW: usize = 1 << 18;

type Ds = PositionalDownsampler<Arc<ReadRec>>;

/// The per-locus accumulators, shared by the read stage (which scatters reads into the slots of loci not yet
/// evaluated) and the evaluation stage (which takes and evaluates the slots of loci up to its frontier). The
/// protocol in [`run`] guarantees that the two stages never use the same slot at the same time.
struct Ring {
    buf: std::cell::UnsafeCell<RingBuf>,
    ns: usize,
}
struct RingBuf {
    /// `ns` accumulators per locus (one per sample)
    accs: Box<[std::cell::UnsafeCell<Acc>]>,
    /// several samples: per locus and sample, the HQ soft-clip counts of the elements that feed GATK's running
    /// average, in pileup order (the average runs over the samples in `splitContextBySampleName` order)
    hq: Box<[std::cell::UnsafeCell<Vec<u32>>]>,
    /// loci per lap (a power of two)
    window: usize,
}
impl RingBuf {
    fn new(window: usize, ns: usize) -> RingBuf {
        let hq_n = if ns > 1 { window * ns } else { 0 };
        RingBuf {
            accs: (0..window * ns).map(|_| std::cell::UnsafeCell::new(Acc::default())).collect(),
            hq: (0..hq_n).map(|_| std::cell::UnsafeCell::new(Vec::new())).collect(),
            window,
        }
    }
}
// SAFETY: slots are only accessed under the read / evaluation protocol of `run` (disjoint slots, ordered by atomics
// and channel messages); the buffer is only replaced by `grow`, while neither stage uses it.
unsafe impl Sync for Ring {}
impl Ring {
    fn new(ns: usize) -> Ring {
        Ring { buf: std::cell::UnsafeCell::new(RingBuf::new(WINDOW, ns)), ns }
    }
    fn window(&self) -> usize {
        // SAFETY: `grow` does not run concurrently with any other use of the ring
        unsafe { (*self.buf.get()).window }
    }
    /// Whether a read spanning `span` bp needs a larger ring first.
    fn too_short_for(&self, span: i64) -> bool {
        span >= self.window() as i64 / 2
    }
    /// Grows the ring to hold a read spanning `span` bp, keeping the accumulators of the live loci `lo..lo + window`
    /// (every slot of a locus < `lo` has been taken, and nothing beyond that range has been scattered).
    /// SAFETY: neither stage may be using the ring (a quiescent point of the protocol of `run`).
    unsafe fn grow(&self, lo: i64, span: i64) {
        let ns = self.ns;
        let old = &mut *self.buf.get();
        let window = (2 * (span as usize + 1)).next_power_of_two().max(old.window);
        let mut new = RingBuf::new(window, ns);
        for pos in lo..lo + old.window as i64 {
            let (i, j) = ((pos as usize) % old.window * ns, (pos as usize) % window * ns);
            for s in 0..ns {
                *new.accs[j + s].get_mut() = std::mem::take(old.accs[i + s].get_mut());
                if ns > 1 {
                    *new.hq[j + s].get_mut() = std::mem::take(old.hq[i + s].get_mut());
                }
            }
        }
        *old = new;
    }
    /// SAFETY: the caller must own the slot of `pos` under the protocol of `run`.
    #[allow(clippy::mut_from_ref)]
    unsafe fn slot(&self, pos: i64, s: usize) -> &mut Acc {
        let b = &*self.buf.get();
        &mut *b.accs[((pos as usize) % b.window) * self.ns + s].get()
    }
    /// SAFETY: as [`Ring::slot`]; several samples only.
    #[allow(clippy::mut_from_ref)]
    unsafe fn hq(&self, pos: i64, s: usize) -> &mut Vec<u32> {
        let b = &*self.buf.get();
        &mut *b.hq[((pos as usize) % b.window) * self.ns + s].get()
    }
}

/// The evaluation stage of one contig shard: region read cache, isActive, activity profile, region cutting,
/// readiness / horizons of the downsampler's RNG schedule, region filling, hand-off to the sink.
struct Front<'a, S: FrontSink> {
    t: &'a Tables,
    contig: &'a str,
    intervals: &'a [(i64, i64)],
    refseq: &'a RefSeq,
    sink: S,
    dump: Option<&'a mut BufWriter<File>>,
    res: FrontResult,
    ring: &'a Ring,
    /// loci <= frontier have been evaluated (their slots are free for the read stage's next lap)
    frontier: &'a AtomicI64,
    fin_iv: usize,
    next_final: i64,
    profile: Profile,
    regions: Vec<Region>,
    filler: RegionFiller,
    // RNG schedule of the downsampling (see variant-rs/crates/downsample/README.md)
    downsampling: bool,
    cov: Coverage,
    tracker: ReadyTracker,
    ready_at: Vec<Option<ReadyAt>>,
    horizons: Vec<Horizon>,
    all_popped: bool,
    resolved: usize,
    /// several samples: `HashMap.hash` of each sample name (Java `String.hashCode` spread)
    sample_hash: Vec<i32>,
    af: hc_genotype::genotyper::AfCalculator,
    /// a segment of the contig (single-run mode), see `SegSpec`
    seg: Option<SegState>,
    /// set when a segment stops (or is void): the read stage stops reading
    stop: &'a AtomicBool,
}

/// `ReferenceConfidenceModel.calcGenotypeLikelihoodsOfRefVsAny` of one sample's pileup, as the integer PLs of
/// `GenotypeBuilder.PL(double[])`.
fn ref_vs_any_pls(acc: &Acc, t: &Tables) -> [i32; 3] {
    let denom = acc.read_count as f64 * t.log10_ploidy;
    let raw = [acc.gl[0] - denom, acc.gl[1] - denom, acc.gl[2] - denom];
    let m = raw.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    raw.map(|g| (-10.0 * (g - m)).min(2147483647.0).round() as i64 as i32)
}

/// `AlleleFrequencyCalculator.calculateSingleSampleBiallelicNonRefPosterior(gl, true)`
fn single_sample_active_prob(gl: &[f64; 3], t: &Tables) -> f64 {
    if max_element_index(gl) == 0 {
        return 0.0;
    }
    let post = [gl[0] + t.prior[0], gl[1] + t.prior[1], gl[2] + t.prior[2]];
    if max_element_index(&post) == 0 {
        0.0
    } else {
        let l = log10_sum_log10(&post);
        1.0 - crate::jmath::pow(10.0, post[0] - l)
    }
}

/// Table size of a `java.util.HashMap` filled from empty with `n` entries (default capacity 16, load factor 0.75).
fn java_hashmap_capacity(n: usize) -> usize {
    let mut cap = 16;
    while n > cap * 3 / 4 {
        cap *= 2;
    }
    cap
}

/// `HashMap.hash(key)` of a Java String: `h ^ (h >>> 16)` with `h = String.hashCode()`.
pub(crate) fn java_string_hash(s: &str) -> i32 {
    let mut h: i32 = 0;
    for c in s.encode_utf16() {
        h = h.wrapping_mul(31).wrapping_add(c as i32);
    }
    h ^ ((h as u32) >> 16) as i32
}

impl<'a, S: FrontSink> Front<'a, S> {
    fn emit(&mut self, upto_stream: i64, eof: bool) {
        let t = hft_now();
        self.emit_inner(upto_stream, eof);
        hft_add(5, t);
    }

    fn emit_inner(&mut self, upto_stream: i64, eof: bool) {
        let Front { filler, regions, sink, res, seg, .. } = self;
        filler.drain(regions, upto_stream, eof, &mut |fr: FilledRegion| {
            if let Some(sg) = seg.as_mut() {
                // regions are scanned for the segment bounds when they pop, before they are filled
                // a region pops >= 350 loci after its start, so a boundary at or before its start is known by then
                if sg.from_b.is_none_or(|b| fr.start < b) || sg.stopped {
                    return;
                }
                if sg.end_b.is_some_and(|b| fr.start >= b) {
                    sg.held.push(fr);
                    return;
                }
            }
            res.assigned += fr.reads.len() as u64;
            sink.region(fr)
        });
    }

    /// The first region starting at or after the end boundary, among the regions popped so far.
    fn seg_scan(&mut self) {
        let Some(sg) = self.seg.as_mut() else { return };
        let Some(b) = sg.end_b else { return };
        while sg.end_idx.is_none() && sg.scanned < self.regions.len() {
            if self.regions[sg.scanned].start >= b {
                sg.end_idx = Some(sg.scanned);
            } else {
                sg.scanned += 1;
            }
        }
    }

    /// Boundary loci (see `SegSpec`), after the isActive value of `pos`.
    fn seg_locus(&mut self, pos: i64, prob: f64) {
        let f = self.profile.fsize;
        let Some(sg) = self.seg.as_mut() else { return };
        if prob != 0.0 {
            sg.zero_run = 0;
        } else if pos == sg.last_locus + 1 {
            sg.zero_run += 1;
        } else {
            sg.zero_run = 1;
        }
        sg.last_locus = pos;
        if sg.zero_run < 2 * f + 1 {
            return;
        }
        let x = pos - f;
        if sg.from_b.is_none() && x >= sg.own_from.unwrap() {
            sg.from_b = Some(x);
        }
        // (possibly the same locus: then the segment owns nothing)
        if sg.from_b.is_some() && sg.end_b.is_none() && sg.cuts.front().is_some_and(|&c| x >= c) {
            sg.end_b = Some(x);
        }
    }

    /// Takes the next segment's range over (a drawing group near the cut).
    fn seg_absorb(&mut self) {
        let sg = self.seg.as_mut().unwrap();
        sg.cuts.pop_front();
        sg.absorbed += 1;
        sg.end_b = None;
        if let Some(e) = sg.end_idx.take() {
            sg.scanned = e;
        }
        sg.reset_stop();
        let held = std::mem::take(&mut sg.held);
        self.sink.absorb_next();
        for fr in held {
            self.res.assigned += fr.reads.len() as u64;
            self.sink.region(fr);
        }
        self.seg_scan();
    }

    /// The segment stops reading (its range is done, or it was absorbed).
    fn seg_stop(&mut self) {
        let sg = self.seg.as_mut().unwrap();
        sg.stopped = true;
        sg.held.clear();
        self.stop.store(true, std::sync::atomic::Ordering::Release);
        // the read stage may wait for loci to be evaluated: nothing is evaluated any more
        self.frontier.store(i64::MAX / 4, std::sync::atomic::Ordering::Release);
    }

    fn seg_stopped(&self) -> bool {
        self.seg.as_ref().is_some_and(|s| s.stopped)
    }

    /// Stop rule of a segment, after processing a kept read starting at `rstart` (see `SegSpec`).
    fn seg_check_stop(&mut self, r: &ReadRec) {
        let Some(sg) = self.seg.as_mut() else { return };
        if sg.stopped {
            return;
        }
        sg.n_proc += 1;
        if sg.n_proc % 1024 == 0 && sg.own_from.is_some() && self.sink.is_absorbed() {
            self.seg_stop();
            return;
        }
        let sg = self.seg.as_mut().unwrap();
        let (Some(e), Some(z)) = (sg.end_idx, sg.end_b) else { return };
        let rstart = r.start;
        if self.downsampling {
            if sg.n_look.is_none() {
                match self.ready_at.get(e).copied().flatten() {
                    Some(a @ ReadyAt::Locus(_)) => sg.n_look = lookahead_locus(Some(a), self.intervals),
                    _ => return,
                }
                if sg.n_look.is_none() {
                    return;
                }
            }
            let n = sg.n_look.unwrap();
            if let Some((ps, pe)) = pileup_span(r.start, r.end, r.flag, r.mate_start, r.tlen) {
                if pe >= n {
                    let c = ps.max(n);
                    sg.c_bound = Some(sg.c_bound.map_or(c, |b| b.min(c)));
                }
            }
            let Some(c) = sg.c_bound else { return };
            if rstart > c && rstart != sg.last_start {
                sg.starts_after += 1;
                sg.last_start = rstart;
            }
            if sg.starts_after < 3 {
                return;
            }
        }
        if rstart >= z + SEG_STOP_MIN {
            self.seg_stop();
        }
    }

    /// Registers the regions popped at `locus` with the readiness tracker.
    fn track(&mut self, locus: Option<i64>, n_before: usize) {
        if !self.downsampling {
            return;
        }
        let ends: Vec<i64> = self.regions[n_before..].iter().map(|r| (r.end + REGION_PADDING).min(self.refseq.contig_len)).collect();
        self.ready_at.resize(self.regions.len(), None);
        match locus {
            Some(l) => {
                if let Some(id) = self.tracker.on_locus(l, &ends) {
                    self.ready_at[id] = Some(ReadyAt::Locus(l));
                }
            }
            None => {
                for id in self.tracker.on_end(&ends) {
                    self.ready_at[id] = Some(ReadyAt::End);
                }
                self.all_popped = true;
            }
        }
    }

    /// isActive with several samples: `splitContextBySampleName` (samples with pileup elements, in Java HashMap
    /// order), the HQ soft-clip running average over them in that order, one sample: the exact single-sample
    /// posterior, else `MinimalGenotypingEngine.calculateGenotypes` on the biallelic N / `<FAKE_ALT>` site.
    /// SAFETY: as the slot reads of `finalize_inner`.
    unsafe fn multi_sample_active(&mut self, pos: i64) -> (f64, Option<f64>) {
        let t = self.t;
        let ns = self.ring.ns;
        let mut present: Vec<(usize, [i32; 3])> = Vec::new();
        let mut hq_vals: Vec<Vec<u32>> = vec![Vec::new(); ns];
        for (s, hv) in hq_vals.iter_mut().enumerate() {
            let acc = std::mem::take(self.ring.slot(pos, s));
            std::mem::swap(hv, self.ring.hq(pos, s));
            if acc.n_elems > 0 {
                present.push((s, ref_vs_any_pls(&acc, t)));
            }
        }
        if present.is_empty() {
            return (0.0, None);
        }
        let cap = java_hashmap_capacity(present.len());
        present.sort_by_key(|&(s, _)| ((self.sample_hash[s] as u32 as usize) & (cap - 1), s));
        let (mut n, mut mean) = (0u64, 0f64);
        for &(s, _) in &present {
            for &v in &hq_vals[s] {
                n += 1;
                mean += (v as f64 - mean) / n as f64;
            }
        }
        let hq = if mean > AVG_HQ_SOFTCLIPS_THRESHOLD { Some(mean) } else { None };
        if present.len() == 1 {
            let gl = present[0].1.map(|p| p as f64 / -10.0);
            return (single_sample_active_prob(&gl, t), hq);
        }
        // GenotypingEngine.calculateGenotypes: alleles N (ref) and <FAKE_ALT> (symbolic: length 0, an indel prior),
        // stand-call-conf min(4.0, 30)
        let pls: Vec<&[i32]> = present.iter().map(|(_, p)| &p[..]).collect();
        let af = self.af.calculate_multi(&pls, &[1, 0], None);
        let conf = 4.0f64;
        let plausible = af.log10_p_ref_by_allele[1] + 1.0e-10 < conf * -0.1;
        let log10_conf = af.log10_p_no_variant + 0.0;
        let phred = (-10.0 * log10_conf) + 0.0;
        if !plausible || !(phred >= conf) {
            return (0.0, hq);
        }
        // QualityUtils.qualToProb(vc.getPhredScaledQual())
        let qual = (log10_conf * -10.0) + 0.0;
        (1.0 - crate::jmath::pow(10.0, qual / -10.0), hq)
    }

    /// Finalize the interval loci up to `upto`: isActive, then the activity profile. Every kept read starting at or
    /// before `upto` must have been scattered.
    fn finalize(&mut self, upto: i64) {
        let tm = hft_now();
        self.finalize_inner(upto);
        hft_add(4, tm);
    }

    fn finalize_inner(&mut self, upto: i64) {
        let t = self.t;
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
            // SAFETY: the read stage has scattered every kept read starting at or before `upto` and only writes slots
            // of loci > the frontier of the read being evaluated (see `run`)
            // HaplotypeCallerEngine.isActive
            let (prob, hq) = if self.ring.ns == 1 {
                let acc = std::mem::take(unsafe { self.ring.slot(pos, 0) });
                if acc.n_elems == 0 {
                    (0.0, None)
                } else {
                    let gl = ref_vs_any_pls(&acc, t).map(|p| p as f64 / -10.0);
                    let hq = if acc.hq_mean > AVG_HQ_SOFTCLIPS_THRESHOLD { Some(acc.hq_mean) } else { None };
                    (single_sample_active_prob(&gl, t), hq)
                }
            } else {
                // SAFETY: as above
                unsafe { self.multi_sample_active(pos) }
            };
            self.res.n_loci += 1;
            if prob > 0.0 {
                self.res.n_active_loci += 1;
            }
            self.seg_locus(pos, prob);
            if let Some(d) = self.dump.as_mut() {
                match hq {
                    Some(v) => writeln!(d, "{}\t{}\t{:.16e}\tHIGH_QUALITY_SOFT_CLIPS\t{:.16e}", self.contig, pos, prob, v).unwrap(),
                    None => writeln!(d, "{}\t{}\t{:.16e}\tNONE", self.contig, pos, prob).unwrap(),
                }
            }
            // AssemblyRegionIterator loop body: pops (forced at a non-contiguous locus) come first
            let n_before = self.regions.len();
            if !self.profile.is_empty() {
                let force = pos != self.profile.stop + 1;
                self.profile.pop_ready(force, &mut self.regions);
            }
            self.seg_scan();
            self.track(Some(pos), n_before);
            self.profile.add(pos, prob, hq);
            self.next_final += 1;
        }
    }

    /// One read leaving the downsampler (kept order), evaluation side: region read cache, loci before it. Every
    /// earlier kept read must have been scattered.
    fn process(&mut self, r: &Arc<ReadRec>) {
        if self.seg_stopped() {
            return;
        }
        let t = hft_now();
        self.process_inner(r);
        self.seg_check_stop(r);
        hft_add(2, t);
    }

    fn process_inner(&mut self, r: &Arc<ReadRec>) {
        let (rstart, rend) = (r.start, r.end);
        if self.downsampling {
            // future horizon queries start at or after the lookahead locus of the first region without a horizon
            let i = self.horizons.len();
            let lb = match self.ready_at.get(i + 1).copied().flatten() {
                Some(a) => lookahead_locus(Some(a), self.intervals).unwrap_or(i64::MAX),
                None => self.next_final,
            };
            self.cov.add(pileup_span(rstart, rend, r.flag, r.mate_start, r.tlen), lb);
        }
        self.emit(rstart, false);
        self.filler.cache.push_back(r.clone());
        // all loci before this read's start are complete (input is coordinate sorted)
        self.finalize(rstart - 1);
        self.frontier.store(rstart - 1, std::sync::atomic::Ordering::Release);
    }

    /// Number of leading regions whose read-pull horizon X is known and lies before position `p`
    /// (GATK applies them, and makes their QD draws, before the draws of the group at `p`).
    fn regions_before(&mut self, p: i64, ds: &mut Ds) -> usize {
        while self.horizons.len() < self.regions.len() {
            let i = self.horizons.len();
            let next: Option<ReadyAt> = if i + 1 < self.regions.len() {
                match self.ready_at.get(i + 1).copied().flatten() {
                    Some(a) => Some(a),
                    None => break,
                }
            } else if self.all_popped {
                None
            } else {
                break;
            };
            let h = match lookahead_locus(next, self.intervals) {
                None => Horizon::All,
                Some(n) => {
                    let limit = ds.first_pending_start().map_or(i64::MAX, |s| s - 1);
                    match self.cov.first_covered(n, limit) {
                        Some(c) => ds.horizon(Some(c)),
                        None if ds.first_pending_start().is_none() && ds.is_finished() => Horizon::All,
                        None => break,
                    }
                }
            };
            if h == Horizon::NeedInput {
                break;
            }
            self.horizons.push(h);
        }
        while self.resolved < self.horizons.len() {
            match self.horizons[self.resolved] {
                Horizon::Through(x) if x < p => self.resolved += 1,
                _ => break,
            }
        }
        self.resolved
    }

    /// The draws of drawing group `p` (the head group of `ds`), placed in GATK's order relative to the QD draws of
    /// the regions; its kept reads are appended to `out`. Every kept read before the group must have been
    /// processed and scattered.
    fn draw(&mut self, p: i64, ds: &mut Ds, out: &mut Vec<Arc<ReadRec>>) {
        if let Some(sg) = self.seg.as_mut() {
            if sg.stopped {
                ds.finalize_head(&mut sg.dummy, out);
                return;
            }
            // a group near the cut after this range: the next segment's warm-up would not see GATK's kept reads
            let max_end = ds.head_group().map_or(p, |(_, rs)| rs.iter().map(|r| r.end).max().unwrap_or(p));
            while let Some(&c) = self.seg.as_ref().unwrap().cuts.front() {
                if p >= c - SEG_ZONE || max_end >= c - SEG_ZONE {
                    self.seg_absorb();
                } else {
                    break;
                }
            }
            let sg = self.seg.as_mut().unwrap();
            // a group of the warm-up, ending before the zone: it changes nothing this segment owns
            if sg.own_from.is_some_and(|c0| p < c0 - SEG_ZONE) {
                ds.finalize_head(&mut sg.dummy, out);
                return;
            }
        }
        // loci < p are complete; fill the regions whose padded span ends before p
        self.finalize(p - 1);
        self.emit(p, false);
        let k = self.regions_before(p, ds);
        let absorbed = match self.sink.rng_for_draw(k) {
            Some(rng) => {
                ds.finalize_head(rng, out);
                false
            }
            None => true,
        };
        if absorbed {
            // the previous segment took this range over: void
            self.seg_stop();
            ds.finalize_head(&mut self.seg.as_mut().unwrap().dummy, out);
        } else if let Some(sg) = self.seg.as_ref() {
            // not absorbed, so the previous segment stopped before p: the first owned region is known
            assert!(sg.from_b.is_some(), "segment of {} at {p}: drawing group before its first region", self.contig);
        }
    }
}

/// The read stage's scatter of kept reads into the locus accumulators (pileup elements for isActive).
struct Scatter<'a> {
    t: &'a Tables,
    intervals: &'a [(i64, i64)],
    refseq: &'a RefSeq,
    ring: &'a Ring,
    scatter_iv: usize,
    local_ivs: Vec<(i64, i64)>,
}

impl Scatter<'_> {
    fn scatter(&mut self, r: &ReadRec) {
        let t_sc = hft_now();
        self.scatter_inner(r);
        hft_add(3, t_sc);
    }

    /// SAFETY (ring): the caller owns the slots of `r`'s loci under the protocol of `run`.
    fn scatter_inner(&mut self, r: &ReadRec) {
        let (rstart, rend) = (r.start, r.end);
        assert!(!self.ring.too_short_for(rend - rstart + 1), "read {}: ring not grown", String::from_utf8_lossy(&r.name));

        // interval parts this read touches
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
        let hq_softclips: u32 = {
            let mut n = 0;
            let mut ro = 0usize;
            for &(op, len) in &r.cigar {
                if op == OP_S {
                    n += r.qual[ro..ro + len as usize].iter().filter(|&&q| q > HQ_SOFTCLIP_QUAL).count() as u32;
                }
                if consumes_read(op) {
                    ro += len as usize;
                }
            }
            n
        };
        // loci excluded as inside the adaptor: <= b (reverse) or >= b (forward)
        let (adapt_lo, adapt_hi) = match r.adaptor_exclusion() {
            Some((true, b)) => (b + 1, i64::MAX),
            Some((false, b)) => (i64::MIN, b - 1),
            None => (i64::MIN, i64::MAX),
        };
        let (t, refseq, ring, local) = (self.t, self.refseq, self.ring, &self.local_ivs);
        let (smp, multi) = (r.sample as usize, ring.ns > 1);
        // neighbour flags (alt, next to soft clip) of the element at offset k of a cigar element of length len
        let edge = |fl: &ElemFlags, k: usize, len: usize| -> (bool, bool) {
            match (k == 0, k == len - 1) {
                (true, true) => (fl.alt_both, fl.sc_both),
                (true, false) => (fl.alt_start, fl.sc_start),
                (false, true) => (fl.alt_end, fl.sc_end),
                (false, false) => (false, false),
            }
        };
        // walk the alignment, scattering each pileup element (loci within the intervals, outside the adaptor) into
        // its locus accumulator; every locus gets one element per read, so the order within the read is free
        let mut rp = rstart;
        let mut ro = 0usize;
        for (ci, &(op, len)) in r.cigar.iter().enumerate() {
            let len = len as usize;
            match op {
                OP_M | OP_EQ | OP_X | OP_D => {
                    let is_del = op == OP_D;
                    let (e_lo, e_hi) = (rp.max(adapt_lo), (rp + len as i64 - 1).min(adapt_hi));
                    if e_lo <= e_hi {
                        let fl = element_flags(&r.cigar, ci);
                        for &(s, e) in local.iter() {
                            let (a, b) = (s.max(e_lo), e.min(e_hi));
                            if a > b {
                                continue;
                            }
                            let k0 = (a - rp) as usize;
                            if is_del {
                                for (i, pos) in (a..=b).enumerate() {
                                    // SAFETY: see the function comment
                                    let acc = unsafe { ring.slot(pos, smp) };
                                    acc.n_elems += 1;
                                    acc.read_count += 1;
                                    let (_, nb_sc) = edge(&fl, k0 + i, len);
                                    add_element(acc, t, true, REF_MODEL_DELETION_QUAL, nb_sc, hq_softclips);
                                    if multi && nb_sc {
                                        // SAFETY: see the function comment
                                        unsafe { ring.hq(pos, smp) }.push(hq_softclips);
                                    }
                                }
                            } else {
                                let n = (b - a + 1) as usize;
                                let (qs, bs) = (&r.qual[ro + k0..ro + k0 + n], &r.seq[ro + k0..ro + k0 + n]);
                                let rf = refseq.window(a, b);
                                for (i, ((&q, &base), &rb)) in qs.iter().zip(bs).zip(rf).enumerate() {
                                    // SAFETY: see the function comment
                                    let acc = unsafe { ring.slot(a + i as i64, smp) };
                                    acc.n_elems += 1;
                                    if q <= MIN_BASE_QUAL {
                                        continue;
                                    }
                                    acc.read_count += 1;
                                    let (nb_alt, nb_sc) = edge(&fl, k0 + i, len);
                                    let is_alt = base != rb || nb_alt;
                                    add_element(acc, t, is_alt, q, nb_sc, hq_softclips);
                                    if multi && is_alt && nb_sc {
                                        // SAFETY: see the function comment
                                        unsafe { ring.hq(a + i as i64, smp) }.push(hq_softclips);
                                    }
                                }
                            }
                        }
                    }
                    rp += len as i64;
                    if !is_del {
                        ro += len;
                    }
                }
                OP_I | OP_S => ro += len,
                OP_N => rp += len as i64,
                _ => {}
            }
        }
    }
}

/// How the read stage hands kept reads and drawing groups to the evaluation stage.
trait Link {
    /// Kept reads, in kept order (drains `out`).
    fn deliver(&mut self, out: &mut Vec<Arc<ReadRec>>);
    /// Drawing group `p` (the head group of `ds`): its kept reads are appended to `out`.
    fn draw(&mut self, p: i64, ds: &mut Ds, out: &mut Vec<Arc<ReadRec>>);
    /// End of the kept reads.
    fn finish(&mut self) {}
}

/// Single thread: per read, the evaluation stage's work, then its scatter.
struct Serial<'f, 'a, S: FrontSink> {
    f: &'f mut Front<'a, S>,
    sc: Scatter<'a>,
}

impl<S: FrontSink> Link for Serial<'_, '_, S> {
    fn deliver(&mut self, out: &mut Vec<Arc<ReadRec>>) {
        for r in out.drain(..) {
            if self.f.ring.too_short_for(r.end - r.start + 1) {
                // SAFETY: single thread, between reads
                unsafe { self.f.ring.grow(self.f.next_final, r.end - r.start + 1) };
            }
            self.f.process(&r);
            self.sc.scatter(&r);
        }
    }
    fn draw(&mut self, p: i64, ds: &mut Ds, out: &mut Vec<Arc<ReadRec>>) {
        self.f.draw(p, ds, out);
    }
}

/// Read stage -> evaluation stage.
enum Msg {
    Reads(Vec<Arc<ReadRec>>),
    /// a drawing group at this start: the evaluation stage makes its draws and sends the downsampler back
    Draw(i64, Ds),
    /// the next read spans this many bp: every earlier read is scattered, and the read stage waits (on `grown`)
    /// while the evaluation stage, after processing them, grows the ring
    Grow(i64),
}

/// Sets the flag when dropped (also on panic), so that the other stage stops waiting.
struct AbortOnDrop<'a>(&'a AtomicBool);
impl Drop for AbortOnDrop<'_> {
    fn drop(&mut self) {
        self.0.store(true, std::sync::atomic::Ordering::Release);
    }
}

/// Spins (then yields) until `ready()`; false if the other stage has stopped first. The stop flag is checked after a
/// failed `ready()`, so a stage that publishes its last progress and then stops (the read stage at the end of a
/// contig) could be seen as stopped with the progress unseen: `ready()` is asked once more after the flag
/// (the flag's Release store follows the progress store, so it is visible then).
fn wait_until(abort: &AtomicBool, mut ready: impl FnMut() -> bool) -> bool {
    let mut spins = 0u32;
    while !ready() {
        if abort.load(std::sync::atomic::Ordering::Acquire) {
            return ready();
        }
        spins += 1;
        if spins < 64 { std::hint::spin_loop() } else { std::thread::yield_now() }
    }
    true
}

/// The read stage of the two-stage front end (see `run`).
struct Staged<'a> {
    sc: Scatter<'a>,
    batch: Vec<Arc<ReadRec>>,
    tx: std::sync::mpsc::SyncSender<Msg>,
    reply: std::sync::mpsc::Receiver<(Ds, Vec<Arc<ReadRec>>)>,
    grown: std::sync::mpsc::Receiver<()>,
    frontier: &'a AtomicI64,
    scattered: &'a AtomicU64,
    abort: &'a AtomicBool,
}

impl Staged<'_> {
    fn flush(&mut self) {
        if self.batch.is_empty() {
            return;
        }
        let b = std::mem::take(&mut self.batch);
        let keep = b.clone();
        self.tx.send(Msg::Reads(b)).expect("evaluation stage stopped");
        for r in &keep {
            // the slots of r's loci must be free: its last locus within one window of the evaluation frontier
            let t = hft_now();
            let w = self.sc.ring.window() as i64;
            let ok = wait_until(self.abort, || r.end < self.frontier.load(std::sync::atomic::Ordering::Acquire) + 1 + w);
            hft_add(11, t);
            assert!(ok, "evaluation stage stopped");
            self.sc.scatter(r);
            self.scattered.fetch_add(1, std::sync::atomic::Ordering::Release);
        }
        self.batch = keep;
        self.batch.clear();
    }
}

impl Link for Staged<'_> {
    fn deliver(&mut self, out: &mut Vec<Arc<ReadRec>>) {
        if out.iter().any(|r| self.sc.ring.too_short_for(r.end - r.start + 1)) {
            for r in out.drain(..) {
                let span = r.end - r.start + 1;
                if self.sc.ring.too_short_for(span) {
                    self.flush();
                    self.tx.send(Msg::Grow(span)).expect("evaluation stage stopped");
                    self.grown.recv().expect("evaluation stage stopped");
                }
                self.batch.push(r);
            }
        }
        self.batch.append(out);
        if self.batch.len() >= 256 {
            self.flush();
        }
    }
    fn draw(&mut self, p: i64, ds: &mut Ds, out: &mut Vec<Arc<ReadRec>>) {
        // everything before the group is scattered once `flush` returns (and sent ahead of the Draw message)
        self.flush();
        let t = hft_now();
        let owned = std::mem::replace(ds, PositionalDownsampler::new(1));
        self.tx.send(Msg::Draw(p, owned)).expect("evaluation stage stopped");
        let (back, kept) = self.reply.recv().expect("evaluation stage stopped");
        hft_add(13, t);
        *ds = back;
        out.extend(kept);
    }
    fn finish(&mut self) {
        self.flush();
    }
}

/// Releases the downsampled groups with start <= `upto` to `link`, drawing groups via `link.draw`.
fn release(ds: &mut Ds, upto: i64, out: &mut Vec<Arc<ReadRec>>, link: &mut impl Link) {
    let t = hft_now();
    loop {
        let stop = ds.finalize_free_upto(upto, out);
        link.deliver(out);
        let Some(p) = stop else { break };
        link.draw(p, ds, out);
        link.deliver(out);
    }
    hft_add(6, t);
}

/// The read stage: decoded reads -> PositionalDownsampler -> `link`, in BAM order. Returns the downsampling stats.
fn read_loop(rrx: std::sync::mpsc::Receiver<Vec<Arc<ReadRec>>>, mut ds: Option<Ds>, link: &mut impl Link, stop: &AtomicBool) -> hc_downsample::DownsampleStats {
    let t_loop = hft_now();
    let mut out: Vec<Arc<ReadRec>> = Vec::new();
    let mut rrx = rrx.into_iter();
    let mut stopped = false;
    'reads: loop {
        let t_w = hft_now();
        let Some(batch) = rrx.next() else { break };
        hft_add(1, t_w);
        for r in batch {
            // a segment that has done its range (or is void) reads no further: no more draws
            if stop.load(std::sync::atomic::Ordering::Acquire) {
                stopped = true;
                break 'reads;
            }
            match ds.as_mut() {
                None => {
                    out.push(r);
                    link.deliver(&mut out);
                }
                Some(d) => {
                    // PositionalDownsampler: groups before this start are complete now
                    let start = r.start;
                    d.push(start, r);
                    release(d, start - 1, &mut out, link);
                }
            }
        }
    }
    let stats = match ds.as_mut() {
        Some(d) if stopped => d.stats,
        Some(d) => {
            d.finish();
            release(d, i64::MAX, &mut out, link);
            d.stats
        }
        None => Default::default(),
    };
    link.finish();
    hft_add(8, t_loop);
    stats
}

/// The decoders' own thread pool (a quarter of the cores, at least 2): in the global pool, a decoder's batch would
/// wait behind the back end's long tasks, and the decoder paces the whole contig.
fn decode_pool() -> &'static rayon::ThreadPool {
    static POOL: std::sync::OnceLock<rayon::ThreadPool> = std::sync::OnceLock::new();
    POOL.get_or_init(|| {
        let n = std::thread::available_parallelism().map_or(4, |n| n.get()) / 4;
        rayon::ThreadPoolBuilder::new().num_threads(n.max(2)).thread_name(|i| format!("decode-{i}")).build().unwrap()
    })
}

/// The interval filter's side of `decode_stream`: reads records in file order and keeps those htsjdk's
/// BAMQueryMultipleIntervalsIteratorFilter returns, a batch at a time.
struct IntervalReader<'a> {
    cs: crate::input::Stream,
    tid: usize,
    padded: &'a [(i64, i64)],
    q_idx: usize,
    body: Vec<u8>,
    spare: Vec<Vec<u8>>,
    end: bool,
}

impl IntervalReader<'_> {
    const BATCH: usize = 1024;
    /// The next record the interval filter keeps.
    fn next(&mut self) -> Option<Vec<u8>> {
        let mut one = self.take(1);
        one.pop()
    }
    fn batch(&mut self) -> Vec<Vec<u8>> {
        self.take(Self::BATCH)
    }
    fn take(&mut self, n: usize) -> Vec<Vec<u8>> {
        let mut out = Vec::with_capacity(n);
        while !self.end && out.len() < n {
            if !self.cs.next_record(&mut self.body) {
                self.end = true;
                break;
            }
            let (rid, rstart, _flag) = bam::peek_core(&self.body);
            if rid != self.tid as i32 {
                if rid > self.tid as i32 {
                    self.end = true;
                }
                continue;
            }
            while self.q_idx < self.padded.len() && self.padded[self.q_idx].1 < rstart {
                self.q_idx += 1;
            }
            if self.q_idx >= self.padded.len() {
                self.end = true;
                break;
            }
            if bam::query_end(&self.body) < self.padded[self.q_idx].0 {
                continue;
            }
            out.push(std::mem::replace(&mut self.body, self.spare.pop().unwrap_or_default()));
        }
        out
    }
}

/// One input's records of `tid` over `padded`: htsjdk's BAMQueryMultipleIntervalsIteratorFilter, then BQSR on the
/// fly and the HC read filters. `each` gets every record the interval filter returns, with the read when it passes
/// the read filters (and, with `merge_keys`, its key for merging several inputs); it returns false to stop.
///
/// With BQSR, the records the interval filter keeps go in batches to the decode pool to be recalibrated (each
/// record's recalibration is independent), up to `DEPTH` batches at a time, and come back in file order; parsing,
/// the read filters and `each` stay on this thread, in file order. (Batching the parsing and filters too, or waiting
/// for each batch before reading on, was slower.)
fn decode_stream(
    cs: crate::input::Stream,
    tid: usize,
    padded: &[(i64, i64)],
    ref_lens: &[i64],
    bqsr: Option<Arc<crate::bqsr::apply::Recalibrator>>,
    rg_sample: Option<&std::collections::HashMap<Vec<u8>, u16>>,
    merge_keys: bool,
    mut each: impl FnMut(Option<Arc<ReadRec>>, Option<crate::mutect::front::MergeKey>) -> bool,
) -> FilterCounts {
    const DEPTH: usize = 8;
    let mut fc = FilterCounts::default();
    let mut rec = Rec::default();
    // parse, read filters, `each` for one (recalibrated) record; false: stop
    let mut handle = |body: &[u8], rec: &mut Rec, fc: &mut FilterCounts| -> bool {
        let t_dec = hft_now();
        bam::parse_record(body, rec);
        let read = if passes_filters(rec, ref_lens, fc) {
            let rend = rec.start + ref_length(&rec.cigar) - 1;
            let mut rr = ReadRec::from_rec(rec, rend);
            if let Some(m) = rg_sample {
                let rg = crate::mutect::front::aux_tags(body).rg;
                rr.sample = *rg.and_then(|g| m.get(g)).unwrap_or_else(|| {
                    crate::gatk_cli::refuse(&format!("read {} has no read group of the header (several samples need one per read)", String::from_utf8_lossy(&rec.name)))
                });
            }
            Some(Arc::new(rr))
        } else {
            None
        };
        hft_add(0, t_dec);
        each(read, merge_keys.then(|| crate::mutect::front::MergeKey::of(rec)))
    };
    let mut reader = IntervalReader { cs, tid, padded, q_idx: 0, body: Vec::new(), spare: Vec::new(), end: false };
    let Some(q) = bqsr else {
        while let Some(body) = reader.next() {
            let go = handle(&body, &mut rec, &mut fc);
            reader.spare.push(body);
            if !go {
                break;
            }
        }
        return fc;
    };
    let (tx, rx) = std::sync::mpsc::channel::<(usize, Vec<Vec<u8>>)>();
    let mut back: std::collections::HashMap<usize, Vec<Vec<u8>>> = std::collections::HashMap::new();
    let (mut sent, mut next) = (0usize, 0usize);
    loop {
        while sent - next < DEPTH && !reader.end {
            let mut batch = reader.batch();
            if batch.is_empty() {
                break;
            }
            let (q, tx, k) = (q.clone(), tx.clone(), sent);
            decode_pool().spawn(move || {
                let mut scratch = crate::bqsr::apply::RecScratch::default();
                for b in batch.iter_mut() {
                    q.recalibrate_record(b, &mut scratch);
                }
                let _ = tx.send((k, batch));
            });
            sent += 1;
        }
        if next == sent {
            break;
        }
        let batch = loop {
            if let Some(b) = back.remove(&next) {
                break b;
            }
            let (k, b) = rx.recv().expect("BQSR worker failed");
            back.insert(k, b);
        };
        next += 1;
        for b in &batch {
            if !handle(b, &mut rec, &mut fc) {
                return fc;
            }
        }
        reader.spare.extend(batch);
    }
    fc
}

/// The downsampler's drawing groups among the reads a contig job over `padded` (its padded, merged intervals) pushes
/// to it, for the reads starting in [lo, hi]: (start, number of reads) of every start with more than `target` reads,
/// in order. The same records as `decode_stream` (htsjdk's interval filter, then the HC read filters, summed over
/// the inputs), without BQSR (the filters do not read base qualities) and without decoding the reads further.
/// Chunks [lo, hi] that partition the starts give the job's groups (a record's fate under the interval filter only
/// depends on the record and `padded`).
pub fn count_drawing_groups(inputs: &[crate::input::Input], tid: usize, padded: &[(i64, i64)], ref_lens: &[i64], lo: i64, hi: i64, target: usize, io_ns: Arc<AtomicU64>) -> Vec<(i64, u32)> {
    let (Some(p0), Some(p1)) = (padded.first().map(|p| p.0), padded.last().map(|p| p.1)) else { return Vec::new() };
    if hi < p0 || lo > p1 {
        return Vec::new();
    }
    // the stream: every padded interval from the chunk's start on (a read starting in the chunk may reach a later
    // interval only); the loop stops past `hi`
    let from = lo.max(p0);
    let query: Vec<(i64, i64)> = padded.iter().filter(|p| p.1 >= from).map(|&(a, b)| (a.max(from), b)).collect();
    let mut fc = FilterCounts::default();
    let mut rec = Rec::default();
    let mut per_input: Vec<Vec<(i64, u32)>> = Vec::with_capacity(inputs.len());
    for inp in inputs {
        let cs = inp.stream(tid, &query, false, io_ns.clone());
        let mut reader = IntervalReader { cs, tid, padded, q_idx: 0, body: Vec::new(), spare: Vec::new(), end: false };
        let mut counts: Vec<(i64, u32)> = Vec::new();
        'stream: loop {
            let batch = reader.batch();
            if batch.is_empty() {
                break;
            }
            for body in &batch {
                let (_, start, _) = bam::peek_core(body);
                if start > hi {
                    break 'stream;
                }
                if start >= lo {
                    bam::parse_record(body, &mut rec);
                    if passes_filters(&rec, ref_lens, &mut fc) {
                        match counts.last_mut() {
                            Some(c) if c.0 == rec.start => c.1 += 1,
                            _ => counts.push((rec.start, 1)),
                        }
                    }
                }
            }
            reader.spare.extend(batch);
        }
        per_input.push(counts);
    }
    let mut all = if per_input.len() == 1 {
        per_input.pop().unwrap()
    } else {
        let mut v: Vec<(i64, u32)> = per_input.into_iter().flatten().collect();
        v.sort_unstable_by_key(|c| c.0);
        let mut m: Vec<(i64, u32)> = Vec::with_capacity(v.len());
        for (s, n) in v {
            match m.last_mut() {
                Some(c) if c.0 == s => c.1 += n,
                _ => m.push((s, n)),
            }
        }
        m
    };
    all.retain(|c| c.1 as usize > target);
    all
}

/// Runs the front end of one contig (`tid` in the BAM, unpadded merged `intervals`). `max_reads_per_start`
/// is `--max-reads-per-alignment-start` (0 = no downsampling). `two_stage`: read stage (downsampling, scatter) on
/// its own thread.
#[allow(clippy::too_many_arguments)]
pub fn run<S: FrontSink>(
    ctx: &FrontCtx,
    tid: usize,
    contig: &str,
    intervals: &[(i64, i64)],
    refseq: &RefSeq,
    max_reads_per_start: usize,
    dump: Option<&mut BufWriter<File>>,
    sink: S,
    two_stage: bool,
    seg: Option<SegSpec>,
) -> (FrontResult, S) {
    let t = ctx.tables;
    let stop = AtomicBool::new(false);
    let contig_len = refseq.contig_len;
    let padded = crate::intervals::pad_merge(intervals, REGION_PADDING, contig_len);
    let ring = Ring::new(ctx.samples.len().max(1));
    let frontier = AtomicI64::new(intervals.first().map_or(0, |i| i.0 - 1));
    let scattered = AtomicU64::new(0);
    let abort = AtomicBool::new(false);
    let mut f = Front {
        t,
        contig,
        intervals,
        refseq,
        sink,
        dump,
        res: FrontResult::default(),
        ring: &ring,
        frontier: &frontier,
        fin_iv: 0,
        next_final: intervals.first().map_or(1, |i| i.0),
        profile: Profile { probs: VecDeque::new(), start: 0, stop: 0, contig_len, kernel: t.kernel.clone(), fsize: t.filter_size },
        regions: Vec::new(),
        filler: RegionFiller { cache: VecDeque::new(), prev: Vec::new(), next_region: 0, contig_len },
        downsampling: max_reads_per_start > 0,
        cov: Coverage::default(),
        tracker: ReadyTracker::new(),
        ready_at: Vec::new(),
        horizons: Vec::new(),
        all_popped: false,
        resolved: 0,
        sample_hash: ctx.samples.iter().map(|s| java_string_hash(s)).collect(),
        af: Default::default(),
        seg: seg.map(|sp| SegState {
            own_from: sp.own_from,
            cuts: sp.cuts.into(),
            from_b: if sp.own_from.is_none() { Some(1) } else { None },
            end_b: None,
            end_idx: None,
            scanned: 0,
            zero_run: 0,
            last_locus: i64::MIN,
            held: Vec::new(),
            n_look: None,
            c_bound: None,
            starts_after: 0,
            last_start: i64::MIN,
            stopped: false,
            n_proc: 0,
            dummy: JavaRandom::gatk(),
            absorbed: 0,
        }),
        stop: &stop,
    };
    let ds: Option<Ds> = if max_reads_per_start > 0 { Some(PositionalDownsampler::new(max_reads_per_start)) } else { None };
    let new_scatter = || Scatter { t, intervals, refseq, ring: &ring, scatter_iv: 0, local_ivs: Vec::new() };
    // decoding thread: BAM records -> htsjdk interval filter -> HC read filters, in BAM order; several inputs:
    // one decoder each, merged as htsjdk's MergingSamRecordIterator does (SAMRecordCoordinateComparator over every
    // record the interval filter returns: GATK merges before the read filters)
    let (rtx, rrx) = std::sync::mpsc::sync_channel::<Vec<Arc<ReadRec>>>(8);
    let ref_lens = ctx.ref_lens.to_vec();
    let bqsr = ctx.bqsr.clone();
    let mut streams: Vec<crate::input::Stream> = ctx.inputs.iter().map(|inp| inp.stream(tid, &padded, ctx.parallel_inflate, ctx.inflate_ns.clone())).collect();
    let rg_maps: Option<Vec<std::collections::HashMap<Vec<u8>, u16>>> = ctx.rg_sample.map(|m| m.to_vec());
    let decoder = if streams.len() == 1 {
        let cs = streams.pop().unwrap();
        let rg0 = rg_maps.map(|mut m| m.swap_remove(0));
        let padded = padded.clone();
        std::thread::Builder::new()
            .name("bam-decode".into())
            .spawn(move || {
                let mut batch: Vec<Arc<ReadRec>> = Vec::with_capacity(4096);
                let mut ok = true;
                let fc = decode_stream(cs, tid, &padded, &ref_lens, bqsr, rg0.as_ref(), false, |r, _| {
                    if let Some(r) = r {
                        batch.push(r);
                        if batch.len() == 4096 {
                            ok = rtx.send(std::mem::replace(&mut batch, Vec::with_capacity(4096))).is_ok();
                        }
                    }
                    ok
                });
                if ok && !batch.is_empty() {
                    let _ = rtx.send(batch);
                }
                fc
            })
            .unwrap()
    } else {
        let padded = padded.clone();
        std::thread::Builder::new()
            .name("bam-merge".into())
            .spawn(move || {
                type Item = (Option<Arc<ReadRec>>, crate::mutect::front::MergeKey);
                std::thread::scope(|scope| {
                    let mut rxs = Vec::new();
                    let mut handles = Vec::new();
                    for (k, cs) in streams.into_iter().enumerate() {
                        let (itx, irx) = std::sync::mpsc::sync_channel::<Vec<Item>>(8);
                        rxs.push(irx);
                        let (padded, ref_lens, bqsr) = (&padded, &ref_lens, bqsr.clone());
                        let rg = rg_maps.as_ref().map(|m| &m[k]);
                        handles.push(std::thread::Builder::new().name("bam-decode".into()).spawn_scoped(scope, move || {
                            let mut batch: Vec<Item> = Vec::with_capacity(4096);
                            let mut ok = true;
                            let fc = decode_stream(cs, tid, padded, ref_lens, bqsr, rg, true, |r, key| {
                                batch.push((r, key.unwrap()));
                                if batch.len() == 4096 {
                                    ok = itx.send(std::mem::replace(&mut batch, Vec::with_capacity(4096))).is_ok();
                                }
                                ok
                            });
                            if ok && !batch.is_empty() {
                                let _ = itx.send(batch);
                            }
                            fc
                        }).unwrap());
                    }
                    {
                        // k-way merge; ties go to the earlier input
                        let mut heads: Vec<VecDeque<Item>> = (0..rxs.len()).map(|_| VecDeque::new()).collect();
                        let mut done = vec![false; rxs.len()];
                        let mut out: Vec<Arc<ReadRec>> = Vec::with_capacity(4096);
                        let mut last: Option<(crate::mutect::front::MergeKey, usize)> = None;
                        'merge: loop {
                            for k in 0..rxs.len() {
                                if heads[k].is_empty() && !done[k] {
                                    match rxs[k].recv() {
                                        Ok(b) => heads[k].extend(b),
                                        Err(_) => done[k] = true,
                                    }
                                }
                            }
                            let mut best: Option<usize> = None;
                            for k in 0..heads.len() {
                                if let Some(h) = heads[k].front() {
                                    if best.is_none_or(|b| crate::mutect::front::key_cmp(&h.1, &heads[b].front().unwrap().1) == std::cmp::Ordering::Less) {
                                        best = Some(k);
                                    }
                                }
                            }
                            let Some(b) = best else { break };
                            if let (Some(r), key) = heads[b].pop_front().unwrap() {
                                // the same read in two inputs: GATK's reads are equal() by their fields and
                                // collapse in its hash-keyed read maps (not reproduced)
                                if let Some((lk, li)) = &last {
                                    if *li != b && crate::mutect::front::key_cmp(lk, &key) == std::cmp::Ordering::Equal {
                                        crate::gatk_cli::refuse(&format!(
                                            "read {} is in two inputs: GATK merges identical reads in some steps only; not supported",
                                            String::from_utf8_lossy(&r.name)
                                        ));
                                    }
                                }
                                last = Some((key, b));
                                out.push(r);
                                if out.len() == 4096 && rtx.send(std::mem::replace(&mut out, Vec::with_capacity(4096))).is_err() {
                                    break 'merge;
                                }
                            }
                        }
                        if !out.is_empty() {
                            let _ = rtx.send(out);
                        }
                        drop(rxs);
                    }
                    let mut fc = FilterCounts::default();
                    for h in handles {
                        fc.add(&h.join().expect("BAM decoding failed"));
                    }
                    fc
                })
            })
            .unwrap()
    };
    if !two_stage {
        let stats = read_loop(rrx, ds, &mut Serial { f: &mut f, sc: new_scatter() }, &stop);
        f.res.ds = stats;
    } else {
        // Read stage (own thread): downsampling and scatter into the locus accumulators. Evaluation stage (this
        // thread): per kept read, in order, region read cache, isActive of the loci before its start, activity
        // profile, regions, hand-off, as the single-threaded front end does. Protocol (as in the Mutect2 front end):
        // the read stage sends a batch of kept reads, then scatters them one by one, counting them in `scattered`;
        // before scattering a read it waits until its last locus is < frontier + 1 + window (a read spanning half the
        // window or more first grows it: `Msg::Grow`, at a point where neither stage uses the ring). The evaluation stage
        // processes read n only when `scattered` >= n (every earlier read is in the accumulators, so all loci before
        // its start are complete) and then publishes frontier = its start - 1. The read stage only writes slots of
        // loci >= the start of the read it scatters, above the frontier and within one window of it.
        // Drawing groups: the read stage scatters everything before the group, then moves the downsampler to the
        // evaluation stage (`Msg::Draw`, after the reads), which evaluates the loci before the group, fills regions,
        // waits for the QD draws GATK makes first (`rng_for_draw`), makes the group's draws and sends the downsampler
        // and the group's kept reads back. The draws therefore happen in the single-threaded order.
        let (mtx, mrx) = std::sync::mpsc::sync_channel::<Msg>(64);
        let (ptx, prx) = std::sync::mpsc::sync_channel::<(Ds, Vec<Arc<ReadRec>>)>(1);
        let (gtx, grx) = std::sync::mpsc::sync_channel::<()>(1);
        let stats = std::thread::scope(|scope| {
            // owned here, so that a failing stage unblocks the other one
            let (mrx, ptx, gtx) = (mrx, ptx, gtx);
            let (frontier, scattered, abort, stop) = (&frontier, &scattered, &abort, &stop);
            let sc = new_scatter();
            let reader = std::thread::Builder::new()
                .name("hc-read".into())
                .spawn_scoped(scope, move || {
                    let _guard = AbortOnDrop(abort);
                    let mut link = Staged { sc, batch: Vec::with_capacity(256), tx: mtx, reply: prx, grown: grx, frontier, scattered, abort };
                    read_loop(rrx, ds, &mut link, stop)
                })
                .unwrap();
            {
                let _guard = AbortOnDrop(abort);
                let mut n = 0u64;
                for m in mrx.iter() {
                    match m {
                        Msg::Reads(b) => {
                            for r in b {
                                if f.seg_stopped() {
                                    n += 1;
                                    continue;
                                }
                                let t = hft_now();
                                let ok = wait_until(abort, || scattered.load(std::sync::atomic::Ordering::Acquire) >= n);
                                hft_add(12, t);
                                assert!(ok, "read stage stopped");
                                f.process(&r);
                                n += 1;
                            }
                        }
                        Msg::Draw(p, mut d) => {
                            let mut out = Vec::new();
                            f.draw(p, &mut d, &mut out);
                            if ptx.send((d, out)).is_err() {
                                break;
                            }
                        }
                        Msg::Grow(span) => {
                            // every earlier read is processed here and scattered by the read stage, which now waits
                            // SAFETY: neither stage uses the ring until the reply
                            unsafe { f.ring.grow(f.next_final, span) };
                            if gtx.send(()).is_err() {
                                break;
                            }
                        }
                    }
                }
            }
            reader.join().expect("read stage failed")
        });
        f.res.ds = stats;
    }
    let t_end = hft_now();
    f.res.filters = decoder.join().expect("BAM decoding failed");
    if !f.seg_stopped() {
        f.finalize(i64::MAX);
        let n_before = f.regions.len();
        if !f.profile.is_empty() {
            f.profile.pop_ready(true, &mut f.regions);
        }
        f.seg_scan();
        f.track(None, n_before);
        f.emit(i64::MAX, true);
    }
    hft_add(7, t_end);
    if let Some(sg) = &f.seg {
        f.res.seg_range = Some((sg.from_b.unwrap_or(-1), sg.end_b.unwrap_or(contig_len + 1), sg.absorbed, sg.stopped));
    }
    f.res.regions = std::mem::take(&mut f.regions);
    (f.res, f.sink)
}
