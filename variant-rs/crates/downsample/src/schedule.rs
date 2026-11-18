//! When GATK's downsampler draws, relative to `apply()` (QD jitter), and the helpers a driver needs
//! to reproduce it. Full derivation and evidence: README.md.
//!
//! GATK's lazy pipeline per shard (contig):
//! ```text
//! BAM -> filters -> ReadsDownsamplingIterator (PushToPullIterator over PositionalDownsampler)
//!     -> ReadCachingIterator -> LocusIteratorByState (PeekableIterator) -> IntervalAlignmentContextIterator
//!     -> AssemblyRegionIterator.loadNextAssemblyRegion -> AssemblyRegionWalker: apply(region)
//! ```
//! `AssemblyRegionIterator.next()` returns region R_i only after it has already loaded R_{i+1},
//! so when `apply(R_i)` runs, the locus iterator stands at the locus A_{i+1} where R_{i+1} became
//! ready, and the reads pulled through the downsampler are fixed by that.
//!
//! **Rule.** Let the regions of a shard be R_0, R_1, ... in apply order, and A_j the interval locus at
//! which R_j became ready (see [`ReadyTracker`]). For apply(R_i):
//! 1. if R_{i+1} does not exist, or became ready at the shard's last interval locus or after the
//!    loci ran out (`ReadyAt::End`): every group of the shard has drawn (`Horizon::All`);
//! 2. else N = the interval locus following A_{i+1} ([`lookahead_locus`]); C = the first locus >= N
//!    whose LocusIteratorByState pileup is non-empty (kept reads only; deletions count, bases inside
//!    the adaptor do not; see [`LibsCoverage`]). No such locus: `Horizon::All`;
//! 3. else the horizon is `PositionalDownsampler::horizon(Some(C))`: P1 = first start > C if P1 keeps
//!    >= 2 reads, else P2 = the start after P1.
//!
//! The draws of every group with start <= horizon happen BEFORE apply(R_i)'s QD draws, the draws of
//! every later group AFTER them. Horizons are non-decreasing in i. Groups that do not draw
//! (n <= target) can be finalized at any time.

use std::collections::VecDeque;

/// When a region became ready in `AssemblyRegionIterator.loadNextAssemblyRegion`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReadyAt {
    /// inside the locus loop, at this interval locus
    Locus(i64),
    /// in the "out of loci" branch (after the final forced popReadyAssemblyRegions)
    End,
}

/// Emulates which locus each region becomes ready at. GATK checks once per locus (one region per
/// locus at most): `pendingRegions` non-empty and `locus > paddedSpan.end` of its head.
/// Feed it, for every interval locus in order, the padded ends of the regions that the activity
/// profile popped at the top of that locus's loop iteration (the pop happens BEFORE the locus's own
/// isActive state is added, as in GATK and variant-rs's `finalize`).
#[derive(Default)]
pub struct ReadyTracker {
    pending: VecDeque<(usize, i64)>,
    next_id: usize,
}

impl ReadyTracker {
    pub fn new() -> Self {
        Self::default()
    }
    /// Registers regions popped at this locus (ids are assigned in pop order = apply order,
    /// starting at 0) and returns the id of the region that becomes ready at `locus`, if any.
    pub fn on_locus(&mut self, locus: i64, popped_padded_ends: &[i64]) -> Option<usize> {
        for &pe in popped_padded_ends {
            self.pending.push_back((self.next_id, pe));
            self.next_id += 1;
        }
        match self.pending.front() {
            Some(&(id, pe)) if locus > pe => {
                self.pending.pop_front();
                Some(id)
            }
            _ => None,
        }
    }
    /// After the last interval locus: registers the final forced pop and returns the ids of all
    /// remaining regions, which become ready one per `next()` call with `ReadyAt::End`.
    pub fn on_end(&mut self, popped_padded_ends: &[i64]) -> Vec<usize> {
        for &pe in popped_padded_ends {
            self.pending.push_back((self.next_id, pe));
            self.next_id += 1;
        }
        self.pending.drain(..).map(|(id, _)| id).collect()
    }
}

/// N: the locus IntervalAlignmentContextIterator advances to after returning the pileup at which
/// R_{i+1} became ready. `intervals` = the shard's (merged, unpadded, sorted) intervals, 1-based
/// inclusive. `None` = the read stream has been drained (`Horizon::All`): R_{i+1} absent, ready at
/// `End`, or ready at the last interval locus.
pub fn lookahead_locus(ready_next: Option<ReadyAt>, intervals: &[(i64, i64)]) -> Option<i64> {
    let a = match ready_next {
        None | Some(ReadyAt::End) => return None,
        Some(ReadyAt::Locus(a)) => a,
    };
    for &(s, e) in intervals {
        if a < s {
            return Some(s);
        }
        if a < e {
            return Some(a + 1);
        }
    }
    None
}

/// `ReadUtils.isBaseInsideAdaptor` as applied by LocusIteratorByState: the part of a read's
/// reference span [start, end] that is NOT excluded from pileups. `None` if all of it is excluded.
/// Arguments are raw BAM fields: flag, 1-based mate start (pnext+1, 0 if unavailable), signed TLEN,
/// and whether the mate is on the same contig (needed only to treat a missing mate as unmapped).
pub fn pileup_span(start: i64, end: i64, flag: u16, mate_start: i64, tlen: i32) -> Option<(i64, i64)> {
    const DEFAULT_ADAPTOR_SIZE: i32 = 100;
    let paired = flag & 0x1 != 0;
    let rev = flag & 0x10 != 0;
    let mate_rev = flag & 0x20 != 0;
    let full = Some((start, end));
    // ReadUtils.hasWellDefinedFragmentSize
    if tlen == 0 || !paired || flag & 0x4 != 0 || flag & 0x8 != 0 || mate_start == 0 || rev == mate_rev {
        return full;
    }
    let well = if rev { end > mate_start } else { start <= mate_start + tlen as i64 };
    if !well {
        return full;
    }
    // signed TLEN comparison: never true for reverse-strand reads (negative TLEN)
    if tlen > DEFAULT_ADAPTOR_SIZE {
        return full;
    }
    let (s, e) = if rev {
        // positions <= mateStart - 1 are inside the adaptor
        (start.max(mate_start), end)
    } else {
        // positions >= start + |tlen| are inside the adaptor
        (start, end.min(start + (tlen as i64).abs() - 1))
    };
    if s > e { None } else { Some((s, e)) }
}

/// Which loci have a non-empty LocusIteratorByState pileup, from the kept reads (each read
/// contributes the single contiguous span returned by [`pileup_span`]; N operators cannot occur
/// because WellformedReadFilter removes them). Add reads in kept order, then query loci in
/// non-decreasing order; a locus' answer is final once every kept read with start <= locus is added.
#[derive(Default)]
pub struct LibsCoverage {
    spans: Vec<(i64, i64)>,
}

impl LibsCoverage {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn add(&mut self, span: Option<(i64, i64)>) {
        if let Some(s) = span {
            self.spans.push(s);
        }
    }
    /// Smallest locus >= `from` covered by an added span, among loci <= `limit`
    /// (`limit` = the largest locus for which all reads starting at or before it have been added).
    /// Also forgets spans ending before `from`.
    pub fn first_covered(&mut self, from: i64, limit: i64) -> Option<i64> {
        self.spans.retain(|&(_, e)| e >= from);
        let best = self.spans.iter().map(|&(s, _)| s.max(from)).min()?;
        if best <= limit { Some(best) } else { None }
    }
}
