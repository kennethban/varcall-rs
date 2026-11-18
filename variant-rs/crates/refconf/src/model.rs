//! `ReferenceConfidenceModel.calculateRefConfidence` (GATK 4.7, default arguments, diploid, one sample).
//!
//! For every position of a region span this builds the pileup GATK's `LocusIteratorByState` would build
//! over the given reads (`AssemblyBasedCallerUtils.getPileupsOverReference`), computes the ref-vs-any
//! genotype likelihoods (`calcGenotypeLikelihoodsOfRefVsAny`, `readsWereRealigned = true`), the
//! indel-informative read count (`calcNReadsWithNoPlausibleIndelsReads`, with GATK's per-read cache) and
//! picks the least confident of the SNP and indel likelihoods (`getGLwithWorstGQ`).

use std::cmp::Ordering;
use std::sync::OnceLock;

use hc_assembly::Read;
use hc_genotype::javamath::java_round;
use hc_genotype::mathutils::{approximate_log10_sum_log10, jacobian};

/// `ReferenceConfidenceModel.REF_MODEL_DELETION_QUAL` (`--reference-model-deletion-quality`, 30)
pub const REF_MODEL_DELETION_QUAL: u8 = 30;
/// `ReferenceConfidenceModel.BASE_QUAL_THRESHOLD`
pub const BASE_QUAL_THRESHOLD: u8 = 6;
/// `--indel-size-to-eliminate-in-ref-model` (10)
pub const INDEL_SIZE_TO_ELIMINATE: usize = 10;
/// `ReferenceConfidenceModel.MAX_N_INDEL_INFORMATIVE_READS`
pub const MAX_N_INDEL_INFORMATIVE_READS: usize = 40;
/// `ReadUtils.DEFAULT_ADAPTOR_SIZE`
const DEFAULT_ADAPTOR_SIZE: i64 = 100;
/// `AlignmentUtils.GAP_CHARACTER`
const GAP: u8 = b'-';

/// GATK numbers dumped by `tables/RefconfTables.java` from the same gatk.jar.
pub struct Tables {
    /// `QualityUtils.qualToProbLog10(byte)`
    pub qual_to_prob_log10: [f64; 256],
    /// `MathUtils.LOG10_ONE_THIRD`
    pub log10_one_third: f64,
    /// `Math.log10(2)`
    pub log10_2: f64,
    /// `ReferenceConfidenceModel.getIndelPLs(2, n).getAsVector()`, n = 0..=40
    pub indel_gls: [[f64; 3]; MAX_N_INDEL_INFORMATIVE_READS + 1],
    /// the same as PLs
    pub indel_pls: [[i32; 3]; MAX_N_INDEL_INFORMATIVE_READS + 1],
    /// per (qual, is_alt): the three genotype-likelihood addends of one pileup element
    /// (`applyPileupElementRefVsNonRefLikelihoodAndCount`, ploidy 2, read weight 1), computed with the
    /// exact same operations as GATK, once per quality value instead of once per element
    pub contrib: Vec<[[f64; 3]; 2]>,
}

pub fn tables() -> &'static Tables {
    static T: OnceLock<Tables> = OnceLock::new();
    T.get_or_init(|| {
        let b = include_bytes!("../tables/refconf_tables.bin");
        let mut it = b.chunks_exact(8).map(|c| f64::from_le_bytes(c.try_into().unwrap()));
        let mut q = [0f64; 256];
        for v in q.iter_mut() {
            *v = it.next().unwrap();
        }
        let log10_one_third = it.next().unwrap();
        let log10_2 = it.next().unwrap();
        let mut indel_gls = [[0f64; 3]; MAX_N_INDEL_INFORMATIVE_READS + 1];
        let mut indel_pls = [[0i32; 3]; MAX_N_INDEL_INFORMATIVE_READS + 1];
        for n in 0..=MAX_N_INDEL_INFORMATIVE_READS {
            for k in 0..3 {
                indel_gls[n][k] = it.next().unwrap();
            }
            indel_pls[n] = gls_to_pls3(&indel_gls[n]);
        }
        let jac = jacobian();
        let contrib = (0..256usize)
            .map(|qual| {
                let mut c = [[0f64; 3]; 2];
                for (alt, cc) in c.iter_mut().enumerate() {
                    let (ref_lk, nonref_lk) = if alt == 1 {
                        ((qual as f64) * -0.1 + log10_one_third, q[qual])
                    } else {
                        (q[qual], (qual as f64) * -0.1 + log10_one_third)
                    };
                    *cc = [
                        1.0 * (ref_lk + log10_2),
                        1.0 * approximate_log10_sum_log10(ref_lk + 0.0, nonref_lk + 0.0, jac),
                        1.0 * (nonref_lk + log10_2),
                    ];
                }
                c
            })
            .collect();
        Tables { qual_to_prob_log10: q, log10_one_third, log10_2, indel_gls, indel_pls, contrib }
    })
}

/// htsjdk `GenotypeLikelihoods.GLsToPLs`
#[inline]
pub fn gls_to_pls3(gls: &[f64; 3]) -> [i32; 3] {
    let mut max = f64::NEG_INFINITY;
    for &g in gls {
        max = max.max(g);
    }
    let f = |g: f64| java_round((-10.0 * (g - max)).min(2147483647.0)) as i32;
    [f(gls[0]), f(gls[1]), f(gls[2])]
}

/// htsjdk `GenotypeLikelihoods.getGQLog10FromLikelihoods(0, gls)`
#[inline]
fn gq_log10_hom_ref(gls: &[f64; 3]) -> f64 {
    let mut qual = f64::NEG_INFINITY;
    for &l in &gls[1..] {
        if l >= qual {
            qual = l;
        }
    }
    qual = gls[0] - qual;
    if qual < 0.0 {
        hc_genotype::genotyper::gq_log10_from_likelihoods(0, gls)
    } else {
        -1.0 * qual
    }
}

/// `GATKVariantContextUtils.calculateGQFromPLs`
pub fn gq_from_pls(pls: &[i32]) -> i32 {
    let (mut first, mut second) = (pls[0], pls[1]);
    if first > second {
        second = first;
        first = pls[1];
    }
    for &c in &pls[2..] {
        if c >= second {
            continue;
        }
        if c <= first {
            second = first;
            first = c;
        } else {
            second = c;
        }
    }
    second - first
}

/// One reference-confidence record (`makeReferenceConfidenceVariantContext`): REF/<NON_REF>, GT 0/0.
#[derive(Clone, Debug)]
pub struct SiteRecord {
    pub pos: i64,
    pub ref_base: u8,
    pub ad: [i32; 2],
    pub gq: i32,
    pub pl: [i32; 3],
}

impl SiteRecord {
    #[inline]
    pub fn dp(&self) -> i32 {
        self.ad[0] + self.ad[1]
    }
}

/// Output of `calculateRefConfidence` for one position: a reference-confidence record or the variant
/// call (index into the caller's call list) that starts there.
#[derive(Clone, Debug)]
pub enum Entry {
    Site(SiteRecord),
    Call(usize),
}

// ------------------------------------------------------------------------------------------
// reads
// ------------------------------------------------------------------------------------------
#[inline]
fn consumes_ref(op: u8) -> bool {
    matches!(op, b'M' | b'D' | b'N' | b'=' | b'X')
}
#[inline]
fn consumes_read(op: u8) -> bool {
    matches!(op, b'M' | b'I' | b'S' | b'=' | b'X')
}
#[inline]
fn on_genome(op: u8) -> bool {
    matches!(op, b'M' | b'D' | b'=' | b'X')
}

/// `ReadCoordinateComparator` on one contig (same as `hc_likelihood::clip::read_coordinate_cmp`).
pub fn read_coordinate_cmp(a: &Read, b: &Read) -> Ordering {
    let rev = |r: &Read| r.flag & 0x10 != 0;
    let paired = |r: &Read| r.flag & 0x1 != 0;
    let mate_start = |r: &Read| if r.mate_is_unmapped() { 0 } else { r.mate_start };
    a.start
        .cmp(&b.start)
        .then_with(|| rev(a).cmp(&rev(b)))
        .then_with(|| a.name.as_bytes().cmp(b.name.as_bytes()))
        .then_with(|| a.gatk_sam_flags().cmp(&b.gatk_sam_flags()))
        .then_with(|| a.mapq.cmp(&b.mapq))
        .then_with(|| {
            if paired(a) && paired(b) {
                a.mate_reference_index().cmp(&b.mate_reference_index()).then_with(|| mate_start(a).cmp(&mate_start(b)))
            } else {
                Ordering::Equal
            }
        })
        .then_with(|| a.tlen.cmp(&b.tlen))
}

/// `ReadUtils.getAdaptorBoundary` (None = CANNOT_COMPUTE_ADAPTOR_BOUNDARY)
fn adaptor_boundary(r: &Read) -> Option<i64> {
    if r.tlen == 0 || r.flag & 0x1 == 0 {
        return None;
    }
    if r.flag & 0x4 != 0 || r.flag & 0x8 != 0 || r.mate_start == 0 {
        return None;
    }
    let rev = r.flag & 0x10 != 0;
    if rev == (r.flag & 0x20 != 0) {
        return None;
    }
    let well = if rev { r.end > r.mate_start } else { r.start <= r.mate_start + r.tlen };
    if !well {
        return None;
    }
    Some(if rev { r.mate_start - 1 } else { r.start + r.tlen.abs() })
}

/// One on-genome cigar element (M/=/X/D) of a read, with what the pileup walk needs.
#[derive(Clone, Copy)]
struct OnGenome {
    /// reference start (1-based) and length
    ref_start: i64,
    len: i64,
    /// read offset of the element's first base (for D: of the base after the deletion)
    read_start: i64,
    /// `getCigarModifiedOffset` of the element's first position: lengths of the preceding elements
    /// that consume reference or are soft clips
    mod_off: i64,
    is_del: bool,
    /// the next on-genome element is a deletion (`isBeforeDeletionStart` at the element's last base)
    next_on_genome_del: bool,
    /// the next cigar element is an insertion (`isBeforeInsertion` at the element's last base)
    next_is_ins: bool,
}

/// Pre-computed per-read data for the pileup walk.
struct PRead<'a> {
    r: &'a Read,
    elems: Vec<OnGenome>,
    /// `isBaseInsideAdaptor`: Some((reverse, boundary))
    adaptor: Option<(bool, i64)>,
    cursor: usize,
    /// indel-informativeness cache (`IDL` transient attribute)
    idl: Option<Vec<bool>>,
}

impl<'a> PRead<'a> {
    fn new(r: &'a Read) -> Self {
        let c = &r.cigar;
        let n = c.len();
        let mut elems = Vec::with_capacity(n);
        let (mut rp, mut ro, mut mo) = (r.start, 0i64, 0i64);
        for (i, &(op, len)) in c.iter().enumerate() {
            let len = len as i64;
            if consumes_ref(op) && op != b'N' {
                let next_on_genome_del = c[i + 1..].iter().find(|e| on_genome(e.0)).map_or(false, |e| e.0 == b'D');
                let next_is_ins = i + 1 < n && c[i + 1].0 == b'I';
                elems.push(OnGenome { ref_start: rp, len, read_start: ro, mod_off: mo, is_del: op == b'D', next_on_genome_del, next_is_ins });
            }
            if consumes_ref(op) || op == b'S' {
                mo += len;
            }
            if consumes_ref(op) {
                rp += len;
            }
            if consumes_read(op) {
                ro += len;
            }
        }
        let adaptor = match adaptor_boundary(r) {
            Some(b) if r.tlen <= DEFAULT_ADAPTOR_SIZE => Some((r.flag & 0x10 != 0, b)),
            _ => None,
        };
        PRead { r, elems, adaptor, cursor: 0, idl: None }
    }

    /// Aligned one-to-one bases and quals (`AlignmentUtils.getBasesAndBaseQualitiesAlignedOneToOne`).
    fn aligned_one_to_one(&self) -> (std::borrow::Cow<'a, [u8]>, std::borrow::Cow<'a, [u8]>) {
        use std::borrow::Cow;
        let r = self.r;
        if !r.cigar.iter().any(|e| e.0 == b'I' || e.0 == b'D') {
            return (Cow::Borrowed(&r.seq), Cow::Borrowed(&r.qual));
        }
        let n: usize = r.cigar.iter().filter(|e| consumes_ref(e.0) || e.0 == b'S').map(|e| e.1 as usize).sum();
        let mut b = vec![0u8; n];
        let mut q = vec![0u8; n];
        let (mut lit, mut pad) = (0usize, 0usize);
        for &(op, len) in &r.cigar {
            let len = len as usize;
            if consumes_read(op) {
                if !consumes_ref(op) {
                    lit += len; // inserted (and soft-clipped) bases are skipped
                } else {
                    b[pad..pad + len].copy_from_slice(&r.seq[lit..lit + len]);
                    q[pad..pad + len].copy_from_slice(&r.qual[lit..lit + len]);
                    lit += len;
                    pad += len;
                }
            } else if consumes_ref(op) {
                for _ in 0..len {
                    b[pad] = GAP;
                    q[pad] = 0;
                    pad += 1;
                }
            }
        }
        (Cow::Owned(b), Cow::Owned(q))
    }
}

// ------------------------------------------------------------------------------------------
// Nucleotide.intersect
// ------------------------------------------------------------------------------------------
fn nucleotide_masks() -> &'static [u8; 256] {
    static M: OnceLock<[u8; 256]> = OnceLock::new();
    M.get_or_init(|| {
        let mut m = [0u8; 256];
        let codes: [(u8, u8); 15] = [
            (b'A', 0b0001),
            (b'C', 0b0010),
            (b'G', 0b0100),
            (b'T', 0b1000),
            (b'R', 0b0101),
            (b'Y', 0b1010),
            (b'S', 0b0110),
            (b'W', 0b1001),
            (b'K', 0b1100),
            (b'M', 0b0011),
            (b'B', 0b1110),
            (b'D', 0b1101),
            (b'H', 0b1011),
            (b'V', 0b0111),
            (b'N', 0b1111),
        ];
        for (c, v) in codes {
            m[c as usize] = v;
            m[c.to_ascii_lowercase() as usize] = v;
        }
        m[b'U' as usize] = 0b1000;
        m[b'u' as usize] = 0b1000;
        m
    })
}

#[inline]
fn is_mismatch_and_not_gap(read_base: u8, ref_base: u8) -> bool {
    let m = nucleotide_masks();
    (m[read_base as usize] & m[ref_base as usize]) == 0 && read_base != GAP
}

// ------------------------------------------------------------------------------------------
// indel informativeness (readHasNoPlausibleIdealsOfSize), ported line by line
// ------------------------------------------------------------------------------------------
fn calculate_baseline_mm_qualities(rb: &[u8], rq: &[u8], read_start: usize, refb: &[u8], ref_start: usize) -> Vec<i32> {
    let n = (rb.len() as i64 - read_start as i64).min(refb.len() as i64 - ref_start as i64).max(0) as usize;
    let mut out = vec![0i32; n];
    let mut sum = 0i32;
    for i in (0..n).rev() {
        if is_mismatch_and_not_gap(rb[read_start + i], refb[ref_start + i]) {
            sum += rq[read_start + i] as i8 as i32;
        }
        out[i] = sum;
    }
    out
}

#[allow(clippy::too_many_arguments)]
fn traverse_end_of_read_for_indel_mismatches(
    inf: &mut Vec<bool>,
    read_start: i64,
    rb: &[u8],
    rq: &[u8],
    last_relevant: i64,
    secondary_break: i64,
    ref_start: i64,
    refb: &[u8],
    baseline: &[i32],
    indel_size: i64,
    insertion: bool,
) {
    let global = baseline[0];
    let mut sum = 0i32;
    let ins_len = if insertion { indel_size } else { 0 };
    let del_len = if insertion { 0 } else { indel_size };
    let n_cmp = (rb.len() as i64 - read_start - ins_len).min(refb.len() as i64 - ref_start - del_len);
    let mut ro = n_cmp + ins_len - 1;
    let mut fo = n_cmp + del_len - 1;
    while ro >= 0 && fo >= 0 {
        let read_base = rb[(read_start + ro) as usize];
        let ref_base = refb[(ref_start + fo) as usize];
        if is_mismatch_and_not_gap(read_base, ref_base) {
            sum += rq[(read_start + ro) as usize] as i8 as i32;
            if sum > global {
                break;
            }
        }
        let site = ro.min(fo);
        if rb[(read_start + site) as usize] != GAP
            && read_start + site < last_relevant
            && read_start + site <= secondary_break
            && baseline[site as usize] >= sum
        {
            bitset_set(inf, (read_start + site) as usize, true);
        }
        ro -= 1;
        fo -= 1;
    }
}

#[inline]
fn bitset_set(b: &mut Vec<bool>, i: usize, v: bool) {
    if i >= b.len() {
        if !v {
            return;
        }
        b.resize(i + 1, false);
    }
    b[i] = v;
}

#[inline]
fn bitset_flip(b: &mut Vec<bool>, from: usize, to: usize) {
    if to > b.len() {
        b.resize(to, false);
    }
    for x in &mut b[from..to] {
        *x = !*x;
    }
}

/// Builds the `IDL` cache of a read the first time it is queried (at `read_start` / `ref_start`).
fn compute_informative_bases(pr: &PRead, read_start: i64, refb: &[u8], ref_start: i64, max_indel: i64) -> Vec<bool> {
    let mut inf: Vec<bool> = Vec::new();
    let read_len = pr.r.seq.len() as i64;
    if !(read_len - read_start < max_indel) && !(refb.len() as i64 - ref_start < max_indel) {
        let secondary = read_len - max_indel;
        let (rb, rq) = pr.aligned_one_to_one();
        if rb.len() as i64 - read_start > max_indel {
            let (last, ref_shorter) = if (rb.len() as i64) < refb.len() as i64 - ref_start + read_start + 1 {
                (rb.len() as i64 - max_indel, false)
            } else {
                (refb.len() as i64 - ref_start + read_start - max_indel + 1, true)
            };
            let baseline = calculate_baseline_mm_qualities(&rb, &rq, read_start as usize, refb, ref_start as usize);
            for size in 1..=max_indel {
                traverse_end_of_read_for_indel_mismatches(&mut inf, read_start, &rb, &rq, last, secondary, ref_start, refb, &baseline, size, false);
                traverse_end_of_read_for_indel_mismatches(&mut inf, read_start, &rb, &rq, last, secondary, ref_start, refb, &baseline, size, true);
            }
            if last <= secondary {
                bitset_flip(&mut inf, 0, last as usize);
                if ref_shorter {
                    crate::stats::inc(&crate::stats::IDL_REF_SHORTER);
                    bitset_set(&mut inf, (last - 1) as usize, false);
                }
            } else {
                bitset_flip(&mut inf, 0, (secondary + 1) as usize);
            }
        }
    }
    inf
}

// ------------------------------------------------------------------------------------------
// calculateRefConfidence
// ------------------------------------------------------------------------------------------
/// A call as `getOverlappingVariantContext` sees it: [start, end] (VC end = start + ref length - 1).
#[derive(Clone, Copy, Debug)]
pub struct CallSpan {
    pub start: i64,
    pub end: i64,
}

/// `ReferenceConfidenceModel.calculateRefConfidence` for one region.
///
/// * `ref_hap`: the reference haplotype = reference bases of the region's padded span, starting at `padded_start`.
/// * `span_start..=span_end`: the region span (one output entry per position).
/// * `reads`: `readLikelihoods.sampleEvidence(0)` (any order; sorted here with `ReadCoordinateComparator`).
/// * `calls`: the variant calls of the region (positions where a call starts emit `Entry::Call`).
pub fn calculate_ref_confidence(
    ref_hap: &[u8],
    padded_start: i64,
    span_start: i64,
    span_end: i64,
    reads: &[&Read],
    calls: &[CallSpan],
) -> Vec<Entry> {
    let t = tables();
    let mut order: Vec<usize> = (0..reads.len()).collect();
    order.sort_by(|&a, &b| read_coordinate_cmp(reads[a], reads[b])); // stable
    let mut pr: Vec<PRead> = order.iter().map(|&i| PRead::new(reads[i])).collect();
    let mut out = Vec::with_capacity((span_end - span_start + 1).max(0) as usize);
    let mut active: Vec<usize> = Vec::new();
    let mut next = 0usize;
    let log10_2 = t.log10_2;
    let mut c = [0u64; 12]; // local path counters, flushed into crate::stats at the end
    for pos in span_start..=span_end {
        // getOverlappingVariantContext: the overlapping call with the largest start (first on ties)
        let mut ov: Option<usize> = None;
        for (k, c) in calls.iter().enumerate() {
            if c.start <= pos && pos <= c.end && ov.map_or(true, |o| c.start > calls[o].start) {
                ov = Some(k);
            }
        }
        // update the active reads (LocusIteratorByState keeps them in insertion = sorted order)
        while next < pr.len() && pr[next].r.start <= pos {
            active.push(next);
            next += 1;
        }
        active.retain(|&i| pr[i].r.end >= pos);
        if let Some(k) = ov {
            if calls[k].start == pos {
                c[0] += 1;
                out.push(Entry::Call(k));
                continue;
            }
            c[1] += 1; // a reference-confidence record inside a call (e.g. under a deletion)
        }
        // One pass over the pileup at pos (LocusIteratorByState order): the ref-vs-any likelihoods
        // (calcGenotypeLikelihoodsOfRefVsAny, ploidy 2, minBaseQual 6, readsWereRealigned = true) and the
        // indel-informative read count (calcNReadsWithNoPlausibleIndelsReads), which looks at every element
        // in the same order and stops after 41 informative reads -- reads after that point do not get their
        // IDL cache built at this position.
        let ref_offset = (pos - padded_start) as usize;
        let ref_base = ref_hap[ref_offset];
        let mut gl = [0f64; 3];
        let (mut ref_depth, mut nonref_depth, mut read_count) = (0i32, 0i32, 0i32);
        let mut n_inf = 0usize;
        let mut inf_done = false;
        for &i in &active {
            let p = &mut pr[i];
            while p.cursor < p.elems.len() && p.elems[p.cursor].ref_start + p.elems[p.cursor].len <= pos {
                p.cursor += 1;
            }
            if p.cursor >= p.elems.len() {
                continue;
            }
            let e = p.elems[p.cursor];
            if e.ref_start > pos {
                continue; // inside an N
            }
            if let Some((rev, b)) = p.adaptor {
                if (rev && pos <= b) || (!rev && pos >= b) {
                    c[2] += 1;
                    continue;
                }
            }
            c[9] += 1;
            let off = pos - e.ref_start;
            let r = p.r;
            let read_off = (e.read_start + off) as usize;
            let qual = if e.is_del { REF_MODEL_DELETION_QUAL } else { r.qual[read_off] };
            if qual <= BASE_QUAL_THRESHOLD && !e.is_del {
                c[3] += 1;
            } else {
                read_count += 1;
                let is_alt = e.is_del || r.seq[read_off] != ref_base;
                c[4] += e.is_del as u64;
                if is_alt {
                    nonref_depth += 1;
                } else {
                    ref_depth += 1;
                }
                let k = &t.contrib[qual as usize][is_alt as usize];
                gl[0] += k[0];
                gl[2] += k[2];
                gl[1] += k[1];
            }
            if inf_done {
                continue;
            }
            let at_end = off == e.len - 1;
            let before_del = !e.is_del && at_end && e.next_on_genome_del;
            let before_ins = at_end && e.next_is_ins;
            if before_del || before_ins || e.is_del {
                continue;
            }
            // getCigarModifiedOffset (the element consumes reference here)
            let offset = e.mod_off + off;
            if p.idl.is_none() {
                let v = compute_informative_bases(p, offset, ref_hap, ref_offset as i64, INDEL_SIZE_TO_ELIMINATE as i64);
                c[5] += 1;
                c[6] += v.is_empty() as u64;
                p.idl = Some(v);
            }
            let v = p.idl.as_ref().unwrap();
            if (offset as usize) < v.len() && v[offset as usize] {
                n_inf += 1;
                if n_inf > MAX_N_INDEL_INFORMATIVE_READS {
                    n_inf = MAX_N_INDEL_INFORMATIVE_READS;
                    c[7] += 1;
                    inf_done = true;
                }
            }
        }
        let denom = read_count as f64 * log10_2;
        for g in gl.iter_mut() {
            *g -= denom;
        }
        // doIndelRefConfCalc
        let snp = [gl[0].min(gl[0]), gl[1].min(gl[0]), gl[2].min(gl[0])];
        let indel = &t.indel_gls[n_inf.min(MAX_N_INDEL_INFORMATIVE_READS)];
        let use_indel = gq_log10_hom_ref(indel) > gq_log10_hom_ref(&snp);
        c[8] += use_indel as u64;
        let pl = if use_indel { t.indel_pls[n_inf] } else { gls_to_pls3(&snp) };
        c[10] += 1;
        out.push(Entry::Site(SiteRecord { pos, ref_base, ad: [ref_depth, nonref_depth], gq: gq_from_pls(&pl), pl }));
    }
    use crate::stats as st;
    let add = |k: &st::Counter, v: u64| { if v > 0 { k.1.fetch_add(v, std::sync::atomic::Ordering::Relaxed); } };
    add(&st::CALL_STARTS, c[0]);
    add(&st::SITES, c[10]);
    add(&st::SITES_UNDER_CALL, c[1]);
    add(&st::ADAPTOR_EXCLUDED, c[2]);
    add(&st::LOW_QUAL_BASES_SKIPPED, c[3]);
    add(&st::DELETION_ELEMENTS, c[4]);
    add(&st::IDL_CACHE_BUILT, c[5]);
    add(&st::IDL_CACHE_EMPTY, c[6]);
    add(&st::INDEL_INFORMATIVE_CAPPED, c[7]);
    add(&st::INDEL_GLS_CHOSEN, c[8]);
    add(&st::PILEUP_ELEMENTS, c[9]);
    out
}

/// `HaplotypeCallerEngine.filterNonPassingReads` (READ_LENGTH_FILTER_THRESHOLD 10, mapping quality 20,
/// `MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE` from the read's `ref_id` / `mate_ref_id`).
pub fn passes_hc_read_filter(r: &Read) -> bool {
    let unclipped = r.seq.len() - r.cigar.iter().filter(|e| e.0 == b'S').map(|e| e.1 as usize).sum::<usize>();
    !(unclipped < 10 || r.mapq < 20 || r.mate_on_other_contig())
}
