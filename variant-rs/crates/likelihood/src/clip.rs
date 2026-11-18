//! Read clipping needed by `AssemblyRegion.trim` (`ReadClipper.hardClipToRegion`) and
//! `ReadCoordinateComparator`.  Copied from `hc_preproc` (validated there against GATK by
//! differential fuzzing) and adapted to `hc_assembly::Read`; the code is unchanged.

use crate::Read;
use std::cmp::Ordering;

pub type CigarElem = (u8, u32);

const FLAG_PAIRED: u16 = 0x1;
const FLAG_MATE_UNMAPPED: u16 = 0x8;
const FLAG_REVERSE: u16 = 0x10;

#[inline]
fn is_paired(r: &Read) -> bool { r.flag & FLAG_PAIRED != 0 }
#[inline]
fn is_reverse(r: &Read) -> bool { r.flag & FLAG_REVERSE != 0 }
/// GATK `mateIsUnmapped()`.
#[inline]
pub fn mate_is_unmapped(r: &Read) -> bool { r.flag & FLAG_MATE_UNMAPPED != 0 || r.mate_start == 0 }
#[inline]
fn mate_start_gatk(r: &Read) -> i64 { if mate_is_unmapped(r) { 0 } else { r.mate_start } }
#[inline]
fn recompute_end(r: &mut Read) { r.end = r.start + ref_len(&r.cigar) - 1; }
/// `ReadUtils.getSoftStart`
#[inline]
fn soft_start(r: &Read) -> i64 {
    let mut s = r.start;
    for &(op, len) in &r.cigar {
        if op == b'S' { s -= len as i64; } else if op != b'H' { break; }
    }
    s
}

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
/// Port of GATK `CigarBuilder` (removeDeletionsAtEnds = true), writing into a caller-owned Vec.
/// Section-order validation (which only throws) is omitted.
pub(crate) struct CigarBuilder<'a> {
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

const READ_INDEX_NOT_FOUND: i64 = -1;

/// `ReadUtils.getReadIndexForReferenceCoordinate(alignmentStart=softStart, cigar, refCoord)`.
/// Soft clips are treated as if aligned to the reference.  Returns (index, operator).
fn read_index_for_ref(soft_start: i64, cigar: &[CigarElem], ref_coord: i64) -> (i64, Option<u8>) {
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
    read_index_for_ref(soft_start(r), &r.cigar, ref_coord)
}

/// htsjdk `CoordMath.overlaps(start, end, start2, end2)` (as used by `Locatable.overlaps`).
#[inline]
pub fn coord_overlaps(start: i64, end: i64, start2: i64, end2: i64) -> bool {
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
    recompute_end(r);
    true
}

/// `ReadClipper.clipRead(HARDCLIP_BASES)` for a list of ops applied in order.
fn clip_read_hard(r: &mut Read, ops: &[(i64, i64)], scratch: &mut Vec<CigarElem>) -> bool {
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
pub fn hard_clip_to_region(r: &mut Read, ref_start: i64, ref_stop: i64, scratch: &mut Vec<CigarElem>) -> bool {
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


/// `ReadCoordinateComparator` for reads of one contig: start, strand, name, GATK SAM flags
/// (`getSAMFlagsForRead`), MAPQ, (both paired) mate reference index and mate start, TLEN.
pub fn read_coordinate_cmp(a: &Read, b: &Read) -> Ordering {
    a.start
        .cmp(&b.start)
        .then_with(|| is_reverse(a).cmp(&is_reverse(b)))
        .then_with(|| a.name.as_bytes().cmp(b.name.as_bytes()))
        .then_with(|| a.gatk_sam_flags().cmp(&b.gatk_sam_flags()))
        .then_with(|| a.mapq.cmp(&b.mapq))
        .then_with(|| {
            if is_paired(a) && is_paired(b) {
                a.mate_reference_index().cmp(&b.mate_reference_index()).then_with(|| mate_start_gatk(a).cmp(&mate_start_gatk(b)))
            } else {
                Ordering::Equal
            }
        })
        .then_with(|| a.tlen.cmp(&b.tlen))
}
