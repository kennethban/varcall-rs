//! Port of GATK 4.7.0.0 HaplotypeCaller per-region read preprocessing:
//! `AssemblyBasedCallerUtils.finalizeRegion(...)` as called from `assembleReads(...)` with
//! HaplotypeCaller's default arguments.
//!
//! Pipeline per raw read (in input order):
//!  1. `ReadUtils.hasWellDefinedFragmentSize(original)` ? `revertSoftClippedBases` : `ReadClipper.hardClipSoftClippedBases`
//!  2. `ReadClipper.hardClipLowQualEnds(read, minBaseQualityScore - 1 = 9)`
//!  3. `HardClipAndPossiblyAddToCollection`: `hardClipAdaptorSequence`, then `hardClipToRegion(paddedSpan)`
//! then a stable sort with `ReadCoordinateComparator` and `cleanOverlappingReadPairs`
//! (overlapping mates: agreeing bases capped at Q20, disagreeing bases set to Q0).
//!
//! All coordinates are 1-based inclusive, as in GATK.  Reads are mutated in place; "empty reads"
//! (GATK's `ReadUtils.emptyRead`) are represented by a `false` return and are always dropped.

use std::cmp::Ordering;
use std::collections::HashMap;

/// One CIGAR element: (operator char such as `b'M'`, length).
pub type CigarElem = (u8, u32);

pub const FLAG_PAIRED: u16 = 0x1;
pub const FLAG_UNMAPPED: u16 = 0x4;
pub const FLAG_MATE_UNMAPPED: u16 = 0x8;
pub const FLAG_REVERSE: u16 = 0x10;
pub const FLAG_MATE_REVERSE: u16 = 0x20;

/// A self-contained aligned read.
///
/// * `start`/`end`: 1-based inclusive alignment start/end (`end = start + refLen(cigar) - 1`).
///   `end` is recomputed from `start` + `cigar` by [`finalize_region`], so callers may leave it stale.
/// * `qual`: raw phred values (NOT +33).
/// * `mate_start`: mate alignment start, 0 if unpaired or mate unmapped.  For a mate on another
///   contig this is the mate's position on *that* contig (GATK uses it as is).
/// * `tlen`: SAM TLEN (signed).
/// * `ref_id` / `mate_ref_id`: BAM reference indices of the read and of its mate (`RNEXT`).  They
///   only matter for `ReadCoordinateComparator` ties (mate reference index) and for the
///   mate-on-another-contig test; `Default` (0 / 0) means "mate on the read's contig".
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Read {
    pub name: String,
    pub flag: u16,
    pub start: i64,
    pub end: i64,
    pub cigar: Vec<CigarElem>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub mapq: u8,
    pub mate_start: i64,
    pub tlen: i64,
    pub ref_id: i32,
    pub mate_ref_id: i32,
    /// BI / BD base insertion / deletion qualities (raw phred); empty = absent (GATK default 45 per base)
    pub bi: Vec<u8>,
    pub bd: Vec<u8>,
    /// sample index (Mutect2; 0 for single-sample callers)
    pub sample: u16,
    /// Mutect2 `OriginalAlignment`: the read's OA tag names another contig
    pub oa_other_contig: bool,
    /// the original BAM record (bytes after block_size; Mutect2 `--bam-output` only)
    pub bam_rec: Option<std::sync::Arc<[u8]>>,
}

impl Read {
    #[inline]
    pub fn is_paired(&self) -> bool {
        self.flag & FLAG_PAIRED != 0
    }
    #[inline]
    pub fn is_reverse(&self) -> bool {
        self.flag & FLAG_REVERSE != 0
    }
    #[inline]
    pub fn mate_is_reverse(&self) -> bool {
        self.flag & FLAG_MATE_REVERSE != 0
    }
    #[inline]
    pub fn is_unmapped(&self) -> bool {
        self.flag & FLAG_UNMAPPED != 0
    }
    /// GATK `mateIsUnmapped()` (only meaningful for paired reads).
    #[inline]
    pub fn mate_is_unmapped(&self) -> bool {
        self.flag & FLAG_MATE_UNMAPPED != 0 || self.mate_start == 0
    }
    /// `ReadUtils.getMateReferenceIndex`: -1 when the mate is unmapped.
    #[inline]
    pub fn mate_reference_index(&self) -> i32 {
        if self.mate_is_unmapped() {
            -1
        } else {
            self.mate_ref_id
        }
    }
    /// Fails `ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE`: paired, mate mapped, and the
    /// mate on another contig (the read itself is mapped here).
    #[inline]
    pub fn mate_on_other_contig(&self) -> bool {
        self.is_paired() && !self.mate_is_unmapped() && self.mate_ref_id != self.ref_id
    }
    /// `ReadUtils.getSAMFlagsForRead`: the flags GATK rebuilds from the read's state (used by
    /// `ReadCoordinateComparator`).  Differs from the raw flag when the mate is "unmapped" by position
    /// only (bit 0x8 added, 0x20 dropped) or for an unpaired read with pair bits set.
    #[inline]
    pub fn gatk_sam_flags(&self) -> u16 {
        sam_flags_for_read(self.flag, self.mate_is_unmapped())
    }
    /// GATK `getMateStart()`: 0 (UNSET_POSITION) when the mate is unmapped.
    #[inline]
    fn mate_start_gatk(&self) -> i64 {
        if self.mate_is_unmapped() {
            0
        } else {
            self.mate_start
        }
    }
    #[inline]
    fn recompute_end(&mut self) {
        self.end = self.start + ref_len(&self.cigar) - 1;
    }
    /// `ReadUtils.getSoftStart`
    #[inline]
    pub fn soft_start(&self) -> i64 {
        let mut s = self.start;
        for &(op, len) in &self.cigar {
            if op == b'S' {
                s -= len as i64;
            } else if op != b'H' {
                break;
            }
        }
        s
    }
}

/// Parameters of finalizeRegion.  `Default` gives HaplotypeCaller's defaults.
#[derive(Clone, Copy, Debug)]
pub struct FinalizeParams {
    /// `minBaseQualityScore - 1` (default 10 - 1 = 9).  Bases with qual <= this are clipped from the ends.
    pub min_tail_quality: u8,
    /// `--dont-use-soft-clipped-bases` (default false).
    pub dont_use_soft_clipped_bases: bool,
    /// `--override-fragment-softclip-check` (default false).
    pub override_softclip_fragment_check: bool,
    /// `!--do-not-correct-overlapping-base-qualities` (default true).
    pub correct_overlapping_base_qualities: bool,
}

impl Default for FinalizeParams {
    fn default() -> Self {
        FinalizeParams {
            min_tail_quality: 9,
            dont_use_soft_clipped_bases: false,
            override_softclip_fragment_check: false,
            correct_overlapping_base_qualities: true,
        }
    }
}

/// `FragmentUtils.HALF_OF_DEFAULT_PCR_SNV_ERROR_QUAL` = round(-10*log10(1e-4)) / 2 = 20.
pub const HALF_OF_DEFAULT_PCR_SNV_ERROR_QUAL: u8 = 20;

// ---------------------------------------------------------------------------------------------
// CIGAR helpers (htsjdk CigarOperator semantics)
// ---------------------------------------------------------------------------------------------

#[inline]
fn consumes_read(op: u8) -> bool {
    matches!(op, b'M' | b'I' | b'S' | b'=' | b'X')
}
#[inline]
fn consumes_ref(op: u8) -> bool {
    matches!(op, b'M' | b'D' | b'N' | b'=' | b'X')
}
#[inline]
fn is_clipping(op: u8) -> bool {
    op == b'S' || op == b'H'
}

/// htsjdk `Cigar.getReferenceLength`.
#[inline]
pub fn ref_len(cigar: &[CigarElem]) -> i64 {
    cigar.iter().filter(|e| consumes_ref(e.0)).map(|e| e.1 as i64).sum()
}

/// htsjdk `Cigar.getReadLength`.
#[inline]
fn cigar_read_len(cigar: &[CigarElem]) -> i64 {
    cigar.iter().filter(|e| consumes_read(e.0)).map(|e| e.1 as i64).sum()
}

/// Port of GATK `CigarBuilder` (removeDeletionsAtEnds = true), writing into a caller-owned Vec.
/// Section-order validation (which only throws) is omitted.
struct CigarBuilder<'a> {
    els: &'a mut Vec<CigarElem>,
    /// GATK's `lastOperator`.  NOTE: can be stale (differ from the op of the last element) after the
    /// "deletion+insertion followed by clip" branch -- replicated on purpose.
    last: Option<u8>,
}

impl<'a> CigarBuilder<'a> {
    fn new(els: &'a mut Vec<CigarElem>) -> Self {
        els.clear();
        CigarBuilder { els, last: None }
    }

    fn last_two_were_deletion_and_insertion(&self) -> bool {
        let n = self.els.len();
        self.last == Some(b'I') && n > 1 && self.els[n - 2].0 == b'D'
    }

    fn add(&mut self, op: u8, len: u32) {
        if len == 0 {
            return;
        }
        if op == b'D' {
            match self.last {
                None => return,
                Some(l) if is_clipping(l) => return,
                Some(b'I') => {
                    let n = self.els.len();
                    if n == 1 || is_clipping(self.els[n - 2].0) {
                        return;
                    }
                }
                _ => {}
            }
        }
        let n = self.els.len();
        match self.last {
            Some(l) if l == op => {
                // merge into the last element (takes `op` as its operator, as GATK does)
                let e = self.els[n - 1];
                self.els[n - 1] = (op, e.1 + len);
            }
            None => {
                self.els.push((op, len));
                self.last = Some(op);
            }
            Some(l) => {
                if is_clipping(op) {
                    if !consumes_read(l) && !is_clipping(l) {
                        // trailing deletion (or N/P) before right clip: drop it
                        self.els[n - 1] = (op, len);
                        self.last = Some(op);
                    } else if self.last_two_were_deletion_and_insertion() {
                        // [.., D, I] + clip -> [.., I, clip]; lastOperator is NOT updated (GATK quirk)
                        self.els[n - 2] = self.els[n - 1];
                        self.els[n - 1] = (op, len);
                    } else {
                        self.els.push((op, len));
                        self.last = Some(op);
                    }
                } else if op == b'D' && l == b'I' {
                    // shift deletions left of insertions; lastOperator stays I
                    if n > 1 && self.els[n - 2].0 == b'D' {
                        self.els[n - 2].1 += len;
                    } else {
                        self.els.insert(n - 1, (op, len));
                    }
                } else {
                    self.els.push((op, len));
                    self.last = Some(op);
                }
            }
        }
    }

    fn make(self) {
        if self.last == Some(b'D') {
            self.els.pop();
        } else if self.last_two_were_deletion_and_insertion() {
            let n = self.els.len();
            self.els.remove(n - 2);
        }
    }
}

/// `CigarUtils.clipCigar(cigar, start, stop, clippingOperator)` (stop exclusive).
fn clip_cigar(cigar: &[CigarElem], start: i64, stop: i64, clip_op: u8, out: &mut Vec<CigarElem>) {
    let clip_left = start == 0;
    let mut b = CigarBuilder::new(out);
    let mut element_start: i64 = 0;
    for &(op, len) in cigar {
        if op == b'H' {
            b.add(op, len);
            continue;
        }
        let element_end = element_start + if consumes_read(op) { len as i64 } else { 0 };
        if element_end <= start || element_start >= stop {
            // deletions at the edge of the clip are meaningless and skipped
            if consumes_read(op) || (element_start != start && element_start != stop) {
                b.add(op, len);
            }
        } else {
            let unclipped = if clip_left { element_end - stop } else { start - element_start };
            let clipped = len as i64 - unclipped;
            if unclipped <= 0 {
                if consumes_read(op) {
                    b.add(clip_op, len);
                }
            } else if clip_left {
                b.add(clip_op, clipped as u32);
                b.add(op, unclipped as u32);
            } else {
                b.add(op, unclipped as u32);
                b.add(clip_op, clipped as u32);
            }
        }
        element_start = element_end;
    }
    b.make();
}

/// `CigarUtils.alignmentStartShift(cigar, numClipped)`.
fn alignment_start_shift(cigar: &[CigarElem], num_clipped: i64) -> i64 {
    let mut ref_clipped = 0i64;
    let mut element_start = 0i64;
    for &(op, len) in cigar {
        if op == b'H' {
            continue;
        }
        let element_end = element_start + if consumes_read(op) { len as i64 } else { 0 };
        if element_end <= num_clipped {
            if consumes_ref(op) {
                ref_clipped += len as i64;
            }
        } else if element_start < num_clipped {
            if consumes_ref(op) {
                ref_clipped += num_clipped - element_start;
            }
            break;
        }
        element_start = element_end;
    }
    ref_clipped
}

/// `CigarUtils.revertSoftClips`.
fn revert_soft_clips(cigar: &[CigarElem], out: &mut Vec<CigarElem>) {
    let mut b = CigarBuilder::new(out);
    for &(op, len) in cigar {
        b.add(if op == b'S' { b'M' } else { op }, len);
    }
    b.make();
}

pub const READ_INDEX_NOT_FOUND: i64 = -1;

/// `ReadUtils.getReadIndexForReferenceCoordinate(alignmentStart=softStart, cigar, refCoord)`.
/// Soft clips are treated as if aligned to the reference.  Returns (index, operator).
pub fn read_index_for_ref(soft_start: i64, cigar: &[CigarElem], ref_coord: i64) -> (i64, Option<u8>) {
    if ref_coord < soft_start {
        return (READ_INDEX_NOT_FOUND, None);
    }
    let mut last_read = 0i64;
    let mut last_ref = soft_start;
    for &(op, len) in cigar {
        let first_read = last_read;
        let first_ref = last_ref;
        if consumes_read(op) {
            last_read += len as i64;
        }
        if consumes_ref(op) || op == b'S' {
            last_ref += len as i64;
        }
        if first_ref <= ref_coord && ref_coord < last_ref {
            let pos = first_read + if consumes_read(op) { ref_coord - first_ref } else { 0 };
            return (pos, Some(op));
        }
    }
    (READ_INDEX_NOT_FOUND, None)
}

#[inline]
fn read_index_of(r: &Read, ref_coord: i64) -> (i64, Option<u8>) {
    read_index_for_ref(r.soft_start(), &r.cigar, ref_coord)
}

/// htsjdk `CoordMath.overlaps(start, end, start2, end2)` (as used by `Locatable.overlaps`).
#[inline]
fn coord_overlaps(start: i64, end: i64, start2: i64, end2: i64) -> bool {
    (start2 >= start && start2 <= end)
        || (end2 >= start && end2 <= end)
        || (start >= start2 && end <= end2)
}

// ---------------------------------------------------------------------------------------------
// Clipping (ReadClipper / ClippingOp, HARDCLIP_BASES and REVERT_SOFTCLIPPED_BASES only)
// Every function returns `false` when GATK would produce an empty (unmapped) read.
// ---------------------------------------------------------------------------------------------

/// `ClippingOp.applyHardClipBases(read, start, stop)` (inclusive read coordinates).
fn apply_hard_clip(r: &mut Read, start: i64, stop: i64, scratch: &mut Vec<CigarElem>) -> bool {
    let len = r.seq.len() as i64;
    let new_len = len - (stop - start + 1);
    if new_len == 0 {
        return false;
    }
    clip_cigar(&r.cigar, start, stop + 1, b'H', scratch);
    if start == 0 {
        r.start += alignment_start_shift(&r.cigar, stop + 1);
        let n = (stop + 1) as usize;
        r.seq.drain(..n);
        r.qual.drain(..n);
        if !r.bi.is_empty() {
            r.bi.drain(..n);
        }
        if !r.bd.is_empty() {
            r.bd.drain(..n);
        }
    } else {
        // GATK copies newLength bases from index 0 (correct for right-tail clips)
        r.seq.truncate(new_len as usize);
        r.qual.truncate(new_len as usize);
        r.bi.truncate(new_len as usize);
        r.bd.truncate(new_len as usize);
    }
    std::mem::swap(&mut r.cigar, scratch);
    r.recompute_end();
    true
}

/// `ReadClipper.clipRead(HARDCLIP_BASES)` for a list of ops applied in order.
pub fn clip_read_hard(r: &mut Read, ops: &[(i64, i64)], scratch: &mut Vec<CigarElem>) -> bool {
    for &(start, stop) in ops {
        let rl = r.seq.len() as i64;
        if start < rl {
            let stop = if stop >= rl { rl - 1 } else { stop };
            if !apply_hard_clip(r, start, stop, scratch) {
                return false;
            }
        }
    }
    !r.seq.is_empty()
}

/// `ReadClipper.hardClipSoftClippedBases(read)`.
fn hard_clip_soft_clipped_bases(r: &mut Read, scratch: &mut Vec<CigarElem>) -> bool {
    if r.seq.is_empty() {
        return true;
    }
    let mut read_index = 0i64;
    let mut cut_left = -1i64;
    let mut cut_right = -1i64;
    let mut right_tail = false;
    for &(op, len) in &r.cigar {
        if op == b'S' {
            if right_tail {
                cut_right = read_index;
            } else {
                cut_left = read_index + len as i64 - 1;
            }
        } else if op != b'H' {
            right_tail = true;
        }
        if consumes_read(op) {
            read_index += len as i64;
        }
    }
    let mut ops = [(0i64, 0i64); 2];
    let mut n = 0;
    if cut_right >= 0 {
        ops[n] = (cut_right, r.seq.len() as i64 - 1);
        n += 1;
    }
    if cut_left >= 0 {
        ops[n] = (0, cut_left);
        n += 1;
    }
    if n == 0 {
        return true;
    }
    clip_read_hard(r, &ops[..n], scratch)
}

/// `AssemblyBasedCallerUtils.revertSoftClippedBases` -> `ReadClipper.revertSoftClippedBases`
/// (ClippingOp(0,0) with REVERT_SOFTCLIPPED_BASES).
fn revert_soft_clipped_bases(r: &mut Read, scratch: &mut Vec<CigarElem>) -> bool {
    if r.seq.is_empty() {
        return true;
    }
    let first_clip = r.cigar.first().is_some_and(|e| is_clipping(e.0));
    let last_clip = r.cigar.last().is_some_and(|e| is_clipping(e.0));
    if r.cigar.is_empty() || !(first_clip || last_clip) {
        return true;
    }
    let new_start = r.soft_start();
    revert_soft_clips(&r.cigar, scratch);
    std::mem::swap(&mut r.cigar, scratch);
    if new_start <= 0 {
        r.start = 1;
        r.recompute_end();
        if !apply_hard_clip(r, 0, -new_start, scratch) {
            return false;
        }
        r.start = 1;
        r.recompute_end();
    } else {
        r.start = new_start;
        r.recompute_end();
    }
    !r.seq.is_empty()
}

/// `ReadClipper.hardClipLowQualEnds(read, lowQual)`.
fn hard_clip_low_qual_ends(r: &mut Read, low_qual: u8, scratch: &mut Vec<CigarElem>) -> bool {
    if r.seq.is_empty() {
        return true;
    }
    let len = r.qual.len() as i64;
    let mut right = len - 1;
    while right >= 0 && r.qual[right as usize] <= low_qual {
        right -= 1;
    }
    let mut left = 0i64;
    while left < len && r.qual[left as usize] <= low_qual {
        left += 1;
    }
    if left > right {
        return false;
    }
    let mut ops = [(0i64, 0i64); 2];
    let mut n = 0;
    if right < len - 1 {
        ops[n] = (right + 1, len - 1);
        n += 1;
    }
    if left > 0 {
        ops[n] = (0, left - 1);
        n += 1;
    }
    if n == 0 {
        return true;
    }
    clip_read_hard(r, &ops[..n], scratch)
}

/// `ReadClipper.clipByReferenceCoordinates(-1, refStop, HARDCLIP_BASES)` (left tail).
fn hard_clip_left_tail(r: &mut Read, ref_stop: i64, scratch: &mut Vec<CigarElem>) -> bool {
    if r.seq.is_empty() {
        return true;
    }
    let (idx, op) = read_index_of(r, ref_stop);
    let op = match op {
        Some(op) => op,
        None => return true, // GATK would NPE here; unreachable for our callers
    };
    // if refStop falls in a deletion, the index after the deletion is returned: step back one base
    let stop = idx - if consumes_read(op) { 0 } else { 1 };
    if stop == READ_INDEX_NOT_FOUND {
        return true;
    }
    if stop < 0 || stop > r.seq.len() as i64 - 1 {
        return true; // GATK throws; unreachable
    }
    clip_read_hard(r, &[(0, stop)], scratch)
}

/// `ReadClipper.clipByReferenceCoordinates(refStart, -1, HARDCLIP_BASES)` (right tail).
fn hard_clip_right_tail(r: &mut Read, ref_start: i64, scratch: &mut Vec<CigarElem>) -> bool {
    if r.seq.is_empty() {
        return true;
    }
    let start = read_index_of(r, ref_start).0;
    let stop = r.seq.len() as i64 - 1;
    if start == READ_INDEX_NOT_FOUND {
        return true;
    }
    if start < 0 || start > stop {
        return true; // GATK throws; unreachable
    }
    clip_read_hard(r, &[(start, stop)], scratch)
}

/// `ReadClipper.hardClipBothEndsByReferenceCoordinates(read, left, right)`.
fn hard_clip_both_ends(r: &mut Read, left: i64, right: i64, scratch: &mut Vec<CigarElem>) -> bool {
    if r.seq.is_empty() || left == right {
        return false;
    }
    if !hard_clip_right_tail(r, right, scratch) {
        return false; // empty read: getEnd() = 0 < left
    }
    if left > r.end {
        return false;
    }
    hard_clip_left_tail(r, left, scratch)
}

/// `ReadClipper.hardClipToRegion(read, refStart, refStop)`.
fn hard_clip_to_region(r: &mut Read, ref_start: i64, ref_stop: i64, scratch: &mut Vec<CigarElem>) -> bool {
    let a_start = r.start;
    let a_stop = r.end;
    if a_start <= ref_stop && a_stop >= ref_start {
        if a_start < ref_start && a_stop > ref_stop {
            hard_clip_both_ends(r, ref_start - 1, ref_stop + 1, scratch)
        } else if a_start < ref_start {
            hard_clip_left_tail(r, ref_start - 1, scratch)
        } else if a_stop > ref_stop {
            hard_clip_right_tail(r, ref_stop + 1, scratch)
        } else {
            true
        }
    } else {
        false
    }
}

/// `ReadUtils.hasWellDefinedFragmentSize(read)`.
pub fn has_well_defined_fragment_size(r: &Read) -> bool {
    if r.tlen == 0 {
        return false;
    }
    if !r.is_paired() {
        return false;
    }
    if r.is_unmapped() || r.mate_is_unmapped() {
        return false;
    }
    if r.is_reverse() == r.mate_is_reverse() {
        return false;
    }
    if r.is_reverse() {
        r.end > r.mate_start_gatk()
    } else {
        r.start <= r.mate_start_gatk() + r.tlen
    }
}

/// `ReadClipper.hardClipAdaptorSequence(read)`.
fn hard_clip_adaptor_sequence(r: &mut Read, scratch: &mut Vec<CigarElem>) -> bool {
    // ReadUtils.getAdaptorBoundary -- evaluated on the *current* (already clipped/reverted) read
    if !has_well_defined_fragment_size(r) {
        return true;
    }
    let boundary = if r.is_reverse() {
        r.mate_start_gatk() - 1
    } else {
        r.start + r.tlen.abs()
    };
    if !(boundary >= r.start && boundary <= r.end) {
        return true;
    }
    if r.is_reverse() {
        hard_clip_left_tail(r, boundary, scratch)
    } else {
        hard_clip_right_tail(r, boundary, scratch)
    }
}

/// `BaseRecalibrationEngine.makeReadTransform()` (defaults: no default base qualities, no OQ reset):
/// consolidate the cigar (`CigarBuilder`), `ReadClipper.hardClipAdaptorSequence`,
/// `ReadClipper.hardClipSoftClippedBases`. False if the read became empty (BQSR then skips it).
pub fn bqsr_clip(r: &mut Read, scratch: &mut Vec<CigarElem>) -> bool {
    {
        let mut b = CigarBuilder::new(scratch);
        for &(op, len) in &r.cigar {
            b.add(op, len);
        }
        b.make();
    }
    std::mem::swap(&mut r.cigar, scratch);
    r.recompute_end();
    if !hard_clip_adaptor_sequence(r, scratch) || r.seq.is_empty() {
        return false;
    }
    if !hard_clip_soft_clipped_bases(r, scratch) {
        return false;
    }
    !r.seq.is_empty()
}

/// The base qualities GATK sets on the read after the first clipping step of `bqsr_clip` that applies (adaptor,
/// then left soft clip, then right soft clip), or None when none applies. GATK's `setBaseQualities` refuses a
/// negative (Java byte >= 128) quality, and later steps only remove bases: so BaseRecalibrator fails on a read
/// exactly when this contains a byte >= 128.
pub fn bqsr_clip_first_quals(r: &Read) -> Option<Vec<u8>> {
    let mut c = r.clone();
    let mut scratch = Vec::new();
    {
        let mut b = CigarBuilder::new(&mut scratch);
        for &(op, len) in &c.cigar {
            b.add(op, len);
        }
        b.make();
    }
    std::mem::swap(&mut c.cigar, &mut scratch);
    c.recompute_end();
    let len0 = c.seq.len();
    if !hard_clip_adaptor_sequence(&mut c, &mut scratch) {
        return None;
    }
    if c.seq.len() != len0 {
        return Some(c.qual);
    }
    let soft = |e: Option<&(u8, u32)>| e.filter(|e| e.0 == b'S').map_or(0, |e| e.1 as usize);
    let (left, right) = (soft(c.cigar.iter().find(|e| e.0 != b'H')), soft(c.cigar.iter().rev().find(|e| e.0 != b'H')));
    if left > 0 {
        Some(c.qual[left..].to_vec())
    } else if right > 0 {
        Some(c.qual[..c.qual.len() - right].to_vec())
    } else {
        None
    }
}

// ---------------------------------------------------------------------------------------------
// Sorting and overlapping-pair correction
// ---------------------------------------------------------------------------------------------

/// `ReadUtils.getSAMFlagsForRead(read)` from the raw SAM flag and GATK's `mateIsUnmapped()`.
pub fn sam_flags_for_read(flag: u16, mate_unmapped: bool) -> u16 {
    let paired = flag & FLAG_PAIRED != 0;
    let unmapped = flag & FLAG_UNMAPPED != 0;
    let mut f = 0u16;
    if paired {
        f |= 0x1;
        if flag & 0x2 != 0 {
            f |= 0x2;
        }
        if mate_unmapped {
            f |= 0x8;
        } else if flag & FLAG_MATE_REVERSE != 0 {
            f |= 0x20;
        }
        f |= flag & 0xc0; // first / second of pair
    }
    if unmapped {
        f |= 0x4;
    } else if flag & FLAG_REVERSE != 0 {
        f |= 0x10;
    }
    f | (flag & 0xf00) // secondary, QC fail, duplicate, supplementary
}

/// `ReadCoordinateComparator` for reads of one contig: start, strand, name, GATK SAM flags, MAPQ,
/// (both paired) mate reference index and mate start, TLEN.
fn read_coordinate_cmp(a: &Read, b: &Read) -> Ordering {
    a.start
        .cmp(&b.start)
        .then_with(|| a.is_reverse().cmp(&b.is_reverse()))
        .then_with(|| a.name.as_bytes().cmp(b.name.as_bytes()))
        .then_with(|| a.gatk_sam_flags().cmp(&b.gatk_sam_flags()))
        .then_with(|| a.mapq.cmp(&b.mapq))
        .then_with(|| {
            if a.is_paired() && b.is_paired() {
                a.mate_reference_index()
                    .cmp(&b.mate_reference_index())
                    .then_with(|| a.mate_start_gatk().cmp(&b.mate_start_gatk()))
            } else {
                Ordering::Equal
            }
        })
        .then_with(|| a.tlen.cmp(&b.tlen))
}

/// `FragmentUtils.adjustQualsOfOverlappingPairedFragments(pair, setConflictingToZero=true,
/// halfOfPcrSnvQual=empty, halfOfPcrIndelQual=empty)`.
fn adjust_overlapping_pair(left: &mut Read, right: &mut Read) {
    let in_order = left.soft_start() < right.soft_start();
    let (first, second) = if in_order { (left, right) } else { (right, left) };
    if first.end < second.start {
        return;
    }
    let (offset, op) = read_index_of(first, second.start);
    match op {
        None => return,
        Some(op) if is_clipping(op) => return,
        _ => {}
    }
    if offset == READ_INDEX_NOT_FOUND {
        return;
    }
    let first_end_base = read_index_of(first, first.end).0;
    let second_end_base = read_index_of(second, second.end).0;
    let second_offset = read_index_of(second, second.start).0;
    let n = (first_end_base - offset).min(second_end_base - second_offset) + 1;
    if n <= 0 || second_offset < 0 {
        return;
    }
    let half = HALF_OF_DEFAULT_PCR_SNV_ERROR_QUAL;
    let (o1, o2, n) = (offset as usize, second_offset as usize, n as usize);
    let fb = &first.seq[o1..o1 + n];
    let sb = &second.seq[o2..o2 + n];
    let fq = &mut first.qual[o1..o1 + n];
    let sq = &mut second.qual[o2..o2 + n];
    for i in 0..n {
        if fb[i] == sb[i] {
            fq[i] = fq[i].min(half);
            sq[i] = sq[i].min(half);
        } else {
            fq[i] = 0;
            sq[i] = 0;
        }
    }
}

/// `AssemblyBasedCallerUtils.cleanOverlappingReadPairs` for a single sample, on reads already
/// sorted by `ReadCoordinateComparator` (FragmentCollection.create pairing rules).
fn clean_overlapping_read_pairs(reads: &mut [Read]) {
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    {
        // one FragmentCollection per sample (splitReadsBySample): mates are paired within their sample
        let mut name_map: HashMap<(u16, &str), usize> = HashMap::new();
        for (i, r) in reads.iter().enumerate() {
            if !r.is_paired() || r.mate_is_unmapped() || r.mate_start_gatk() == 0 || r.mate_start_gatk() > r.end {
                continue; // singleton
            }
            if let Some(j) = name_map.remove(&(r.sample, r.name.as_str())) {
                pairs.push((j, i));
            } else {
                name_map.insert((r.sample, r.name.as_str()), i);
            }
        }
    }
    for (j, i) in pairs {
        // j < i always
        let (a, b) = reads.split_at_mut(i);
        adjust_overlapping_pair(&mut a[j], &mut b[0]);
    }
}

// ---------------------------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------------------------

/// Per-read part of finalizeRegion (everything before sorting).  Returns false if the read is
/// dropped.  `scratch` is a reusable CIGAR buffer.
fn process_read(r: &mut Read, padded_start: i64, padded_end: i64, p: &FinalizeParams, scratch: &mut Vec<CigarElem>) -> bool {
    r.recompute_end();
    // hasWellDefinedFragmentSize is evaluated on the ORIGINAL read
    let ok = if p.dont_use_soft_clipped_bases
        || !(p.override_softclip_fragment_check || has_well_defined_fragment_size(r))
    {
        hard_clip_soft_clipped_bases(r, scratch)
    } else {
        revert_soft_clipped_bases(r, scratch)
    };
    if !ok || !hard_clip_low_qual_ends(r, p.min_tail_quality, scratch) {
        return false;
    }
    // HardClipAndPossiblyAddToCollection
    if !(r.start <= r.end && !r.is_unmapped()) {
        return false;
    }
    if !hard_clip_adaptor_sequence(r, scratch) {
        return false;
    }
    if r.seq.is_empty() || cigar_read_len(&r.cigar) <= 0 {
        return false;
    }
    let (as_, ae) = (r.start, r.end);
    if !hard_clip_to_region(r, padded_start, padded_end, scratch) {
        return false;
    }
    r.start <= r.end && !r.seq.is_empty() && coord_overlaps(as_, ae, padded_start, padded_end)
}

/// finalizeRegion on owned reads (no cloning).  `reads` must be in GATK's region order
/// (`region.getReads()`); the result is in GATK's final order.
pub fn finalize_region_owned_with(mut reads: Vec<Read>, padded_start: i64, padded_end: i64, p: &FinalizeParams) -> Vec<Read> {
    let mut scratch: Vec<CigarElem> = Vec::with_capacity(16);
    reads.retain_mut(|r| process_read(r, padded_start, padded_end, p, &mut scratch));
    // List.sort is a stable merge sort (TimSort); slice::sort_by is stable too.
    reads.sort_by(read_coordinate_cmp);
    if p.correct_overlapping_base_qualities {
        clean_overlapping_read_pairs(&mut reads);
    }
    reads
}

/// finalizeRegion on owned reads with HaplotypeCaller defaults.
pub fn finalize_region_owned(reads: Vec<Read>, padded_start: i64, padded_end: i64) -> Vec<Read> {
    finalize_region_owned_with(reads, padded_start, padded_end, &FinalizeParams::default())
}

/// finalizeRegion with explicit parameters (clones the input reads).
pub fn finalize_region_with(raw_reads: &[Read], padded_start: i64, padded_end: i64, p: &FinalizeParams) -> Vec<Read> {
    let mut scratch: Vec<CigarElem> = Vec::with_capacity(16);
    let mut out: Vec<Read> = Vec::with_capacity(raw_reads.len());
    for raw in raw_reads {
        let mut r = raw.clone();
        if process_read(&mut r, padded_start, padded_end, p, &mut scratch) {
            out.push(r);
        }
    }
    out.sort_by(read_coordinate_cmp);
    if p.correct_overlapping_base_qualities {
        clean_overlapping_read_pairs(&mut out);
    }
    out
}

/// GATK `AssemblyBasedCallerUtils.finalizeRegion` with HaplotypeCaller defaults.
/// `padded_start`/`padded_end`: the region's padded span (1-based inclusive).
pub fn finalize_region(raw_reads: &[Read], padded_start: i64, padded_end: i64) -> Vec<Read> {
    finalize_region_with(raw_reads, padded_start, padded_end, &FinalizeParams::default())
}

/// Parse a CIGAR string ("" or "*" -> empty).
pub fn parse_cigar(s: &str) -> Vec<CigarElem> {
    let mut v = Vec::new();
    let mut n: u32 = 0;
    for &c in s.as_bytes() {
        if c.is_ascii_digit() {
            n = n * 10 + (c - b'0') as u32;
        } else if c != b'*' {
            v.push((c, n));
            n = 0;
        }
    }
    v
}

/// Format a CIGAR as htsjdk does ("" for an empty cigar).
pub fn cigar_to_string(c: &[CigarElem]) -> String {
    use std::fmt::Write;
    let mut s = String::new();
    for &(op, len) in c {
        let _ = write!(s, "{}{}", len, op as char);
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cig(s: &str) -> Vec<CigarElem> {
        parse_cigar(s)
    }

    #[test]
    fn builder_quirk_stale_last_operator() {
        // 50M2D3I45M, hard clip right from read index 53 -> 50M3I45H (deletion dropped)
        let mut out = Vec::new();
        clip_cigar(&cig("50M2D3I45M"), 53, 98, b'H', &mut out);
        assert_eq!(cigar_to_string(&out), "50M3I45H");
    }

    #[test]
    fn builder_quirk_unmerged_hard_clips() {
        let mut out = Vec::new();
        clip_cigar(&cig("50M2D3I45M5H"), 53, 98, b'H', &mut out);
        assert_eq!(cigar_to_string(&out), "50M3I45H5H");
    }

    #[test]
    fn left_clip_shift() {
        let c = cig("10M5D20M");
        assert_eq!(alignment_start_shift(&c, 10), 15);
        let mut out = Vec::new();
        clip_cigar(&c, 0, 10, b'H', &mut out);
        assert_eq!(cigar_to_string(&out), "10H20M");
    }
}
