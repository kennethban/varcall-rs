//! `PositionalDownsampler` wrapping a `ReservoirDownsampler` (GATK 4.7, random mode), as
//! attached by `AssemblyRegionWalker.createDownsampler()` to each per-contig read shard.
//!
//! Semantics reproduced (see README for the GATK evidence):
//! * input = the shard's reads AFTER the HaplotypeCaller read filters (the downsampler is the
//!   last wrapper of `MultiIntervalLocalReadShard.iterator()`), in BAM order;
//! * one downsampler for all reads of the shard: NOT per sample / read group / strand;
//! * a "position" is `(contig, alignment start)` = `getAssignedStart()` = BAM POS (1-based,
//!   soft clips not included). Reads with no assigned position never reach it on HC's path
//!   (MappedReadFilter runs first), so they are not modelled;
//! * per position: Vitter's Algorithm R with `Utils.getRandomGenerator().nextInt(totalReadsSeen)`
//!   for every read after the first `target`; reads are kept in RESERVOIR SLOT order, i.e. a
//!   replacing read takes the slot (list index) of the read it evicts. So the output is not in
//!   BAM order within a start position once the position overflowed;
//! * the draws for a position all happen inside one `PushToPullIterator.fillCache` call, so they
//!   form one contiguous block of `nextInt` calls ("drawing group"). A group with
//!   `n <= target` reads makes no draw at all.
//!
//! This type separates *buffering* (`push`) from *finalizing* (running the reservoir, which is
//! when the random draws happen), so a driver can place each drawing group's draws exactly where
//! GATK places them relative to the QD-jitter draws of `apply()` (see [`crate::schedule`]).

use crate::javarandom::JavaRandom;
use std::collections::VecDeque;

/// Algorithm R exactly as `ReservoirDownsampler.submit` (random mode) runs it over the reads of one
/// alignment start, in submission order. Returns the kept reads in GATK output order (slot order).
/// Draws `nextInt(i)` for i = target+1 ..= n (none if n <= target).
pub fn reservoir_downsample<T>(reads: impl IntoIterator<Item = T>, target: usize, rng: &mut JavaRandom) -> Vec<T> {
    assert!(target > 0, "targetCoverage must be > 0");
    let mut reservoir: Vec<T> = Vec::new();
    let mut seen: usize = 0;
    for r in reads {
        seen += 1;
        if seen <= target {
            reservoir.push(r);
        } else {
            assert!(seen <= i32::MAX as usize);
            let slot = rng.next_int(seen as i32) as usize;
            if slot < target {
                reservoir[slot] = r;
            }
        }
    }
    reservoir
}

/// The draws `reservoir_downsample` makes over `n` reads, without the reads: `nextInt(i)` for i = target+1 ..= n.
/// (How many reads a group has does not depend on the generator, so a replay of the group sizes gives the
/// generator's state after them.)
pub fn reservoir_draws(n: usize, target: usize, rng: &mut JavaRandom) {
    assert!(target > 0, "targetCoverage must be > 0");
    for seen in target + 1..=n {
        assert!(seen <= i32::MAX as usize);
        rng.next_int(seen as i32);
    }
}

struct Group<T> {
    start: i64,
    reads: Vec<T>,
}

/// Where GATK's read pull stands when `apply(R_i)` is called (see [`crate::schedule`]):
/// every group with start <= the given position has been finalized (drawn), no later one has.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Horizon {
    /// groups with start <= this position
    Through(i64),
    /// the whole shard (input exhausted)
    All,
    /// not decidable yet: push more reads (or call `finish`) first
    NeedInput,
}

impl Horizon {
    /// does a group starting at `start` lie within the horizon?
    pub fn covers(&self, start: i64) -> bool {
        match *self {
            Horizon::Through(p) => start <= p,
            Horizon::All => true,
            Horizon::NeedInput => panic!("horizon not known yet"),
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct DownsampleStats {
    pub submitted: u64,
    pub kept: u64,
    pub groups: u64,
    /// groups with more than `target` reads (the only ones that draw)
    pub drawing_groups: u64,
    /// `nextInt` calls (= sum over drawing groups of n - target)
    pub draws: u64,
}

/// One per read shard (GATK: one per contig, recreated by `AssemblyRegionWalker.traverse` for each
/// shard). Generic over the driver's read type.
pub struct PositionalDownsampler<T> {
    target: usize,
    pending: VecDeque<Group<T>>,
    /// (start, kept count) of recently finalized groups, kept so that `horizon` still sees them
    history: VecDeque<(i64, usize)>,
    last_start: i64,
    input_done: bool,
    pub stats: DownsampleStats,
}

impl<T> PositionalDownsampler<T> {
    /// `target` = `--max-reads-per-alignment-start` (HaplotypeCaller default 50; 0 disables
    /// downsampling in GATK, in which case the driver should not create a downsampler).
    pub fn new(target: usize) -> Self {
        assert!(target > 0, "targetCoverage must be > 0");
        PositionalDownsampler {
            target,
            pending: VecDeque::new(),
            history: VecDeque::new(),
            last_start: i64::MIN,
            input_done: false,
            stats: DownsampleStats::default(),
        }
    }

    pub fn target(&self) -> usize {
        self.target
    }

    /// Submit the next read of the shard that passed the read filters, in BAM order.
    /// `start` = 1-based alignment start (POS). Never draws.
    pub fn push(&mut self, start: i64, read: T) {
        assert!(!self.input_done, "push after finish");
        // GATK: IllegalStateException("Reads must be coordinate sorted")
        assert!(start >= self.last_start, "reads must be coordinate sorted ({} after {})", start, self.last_start);
        self.stats.submitted += 1;
        if start == self.last_start {
            self.pending.back_mut().unwrap().reads.push(read);
        } else {
            self.last_start = start;
            self.pending.push_back(Group { start, reads: vec![read] });
        }
    }

    /// End of the shard's read stream.
    pub fn finish(&mut self) {
        self.input_done = true;
    }

    pub fn is_finished(&self) -> bool {
        self.input_done
    }

    /// Start of the most recently pushed read.
    pub fn last_pushed_start(&self) -> Option<i64> {
        if self.last_start == i64::MIN { None } else { Some(self.last_start) }
    }

    /// Start of the first group not yet finalized, if any.
    pub fn first_pending_start(&self) -> Option<i64> {
        self.pending.front().map(|g| g.start)
    }

    /// Start of the first not-yet-finalized group that is known to draw (has > target reads),
    /// or that may still draw. Returns `(start, complete)`.
    pub fn first_drawing_group(&self) -> Option<(i64, bool)> {
        let n = self.pending.len();
        for (i, g) in self.pending.iter().enumerate() {
            let complete = i + 1 < n || self.input_done;
            if g.reads.len() > self.target {
                return Some((g.start, complete));
            }
        }
        None
    }

    fn group_complete(&self, idx: usize) -> bool {
        idx + 1 < self.pending.len() || self.input_done
    }

    /// Start and reads of the head group (the next one to be finalized), if any.
    pub fn head_group(&self) -> Option<(i64, &[T])> {
        self.pending.front().map(|g| (g.start, &g.reads[..]))
    }

    /// Does finalizing the head group draw from the RNG? (`None` if there is no complete head group.)
    pub fn head_draws(&self) -> Option<bool> {
        if self.pending.is_empty() || !self.group_complete(0) {
            return None;
        }
        Some(self.pending[0].reads.len() > self.target)
    }

    /// Finalize the head group (it must be complete) and append its kept reads to `out`.
    /// Returns the group's start.
    pub fn finalize_head(&mut self, rng: &mut JavaRandom, out: &mut Vec<T>) -> i64 {
        assert!(self.group_complete(0), "head group is not complete");
        let g = self.pending.pop_front().unwrap();
        let n = g.reads.len();
        self.stats.groups += 1;
        if n > self.target {
            self.stats.drawing_groups += 1;
            self.stats.draws += (n - self.target) as u64;
        }
        let before = out.len();
        if n <= self.target {
            out.extend(g.reads);
        } else {
            out.extend(reservoir_downsample(g.reads, self.target, rng));
        }
        let kept = out.len() - before;
        self.stats.kept += kept as u64;
        self.history.push_back((g.start, kept));
        g.start
    }

    /// Finalize every complete group with start <= `upto` (drawing ones included), in order.
    pub fn finalize_upto(&mut self, upto: i64, rng: &mut JavaRandom, out: &mut Vec<T>) {
        while let Some(g) = self.pending.front() {
            if g.start > upto || !self.group_complete(0) {
                break;
            }
            self.finalize_head(rng, out);
        }
    }

    /// Finalize complete groups with start <= `upto` as long as they do not draw.
    /// Returns `Some(start)` of the drawing group it stopped at (start <= upto), else `None`.
    pub fn finalize_free_upto(&mut self, upto: i64, out: &mut Vec<T>) -> Option<i64> {
        let mut dummy = JavaRandom::gatk();
        while let Some(g) = self.pending.front() {
            if g.start > upto || !self.group_complete(0) {
                break;
            }
            if g.reads.len() > self.target {
                return Some(g.start);
            }
            self.finalize_head(&mut dummy, out);
        }
        None
    }

    /// Finalize everything (requires `finish()`): GATK's `signalEndOfInput`.
    pub fn finalize_all(&mut self, rng: &mut JavaRandom, out: &mut Vec<T>) {
        assert!(self.input_done, "finalize_all before finish");
        self.finalize_upto(i64::MAX, rng, out);
    }

    /// GATK's read-pull horizon when LocusIteratorByState has just produced its alignment context
    /// at locus `c` (`None` = LIBS exhausted): LIBS has consumed every kept read with start <= c
    /// (= k_1..k_m), its PeekableIterator holds k_{m+1}, and PushToPullIterator holds k_{m+2}, so the
    /// groups up to and including the start of k_{m+2} have been finalized (their draws done).
    /// That start is P1 (the first start > c) if P1 keeps >= 2 reads, else P2 (the next start).
    ///
    /// Group sizes after `c` are taken from pushed input; the kept count of a group is
    /// min(n, target), which does not depend on the RNG.
    pub fn horizon(&mut self, c: Option<i64>) -> Horizon {
        let c = match c {
            None => return Horizon::All,
            Some(c) => c,
        };
        // drop history that can never matter again (callers pass non-decreasing c)
        while let Some(&(s, _)) = self.history.front() {
            if s <= c { self.history.pop_front(); } else { break; }
        }
        // (start, kept-count-known, kept) of groups after c, in order
        let np = self.pending.len();
        let mut it = self.history.iter().map(|&(s, k)| (s, Some(k))).chain(self.pending.iter().enumerate().map(|(i, g)| {
            let complete = i + 1 < np || self.input_done;
            let n = g.reads.len();
            // kept = min(n, target); only "1 or >= 2" matters, which is known once n >= min(target, 2)
            let k = if complete || n >= self.target.min(2) { Some(n.min(self.target)) } else { None };
            (g.start, k)
        })).filter(|&(s, _)| s > c);
        match it.next() {
            None => if self.input_done { Horizon::All } else { Horizon::NeedInput },
            Some((p1, k1)) => match k1 {
                None => Horizon::NeedInput,
                Some(k) if k >= 2 => Horizon::Through(p1),
                Some(_) => match it.next() {
                    Some((p2, _)) => Horizon::Through(p2),
                    None => if self.input_done { Horizon::All } else { Horizon::NeedInput },
                },
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn draws_without_reads_match_the_reservoir() {
        let (mut a, mut b) = (JavaRandom::gatk(), JavaRandom::gatk());
        for n in [0usize, 1, 50, 51, 77, 1000, 64, 65, 129] {
            reservoir_downsample(0..n, 50, &mut a);
            reservoir_draws(n, 50, &mut b);
            assert_eq!(a.state(), b.state());
        }
    }

    #[test]
    fn no_draw_when_under_target() {
        let mut rng = JavaRandom::gatk();
        let v = reservoir_downsample(0..50, 50, &mut rng);
        assert_eq!(v, (0..50).collect::<Vec<_>>());
        assert_eq!(rng.steps, 0);
    }

    #[test]
    fn slot_order_and_draw_count() {
        let mut rng = JavaRandom::gatk();
        let v = reservoir_downsample(0..10, 3, &mut rng);
        assert_eq!(v.len(), 3);
        assert_eq!(rng.steps, 7); // nextInt(4..=10), no rejection for such bounds here
        // replay by hand
        let mut r2 = JavaRandom::gatk();
        let mut res = vec![0, 1, 2];
        for i in 3..10 {
            let s = r2.next_int(i + 1) as usize;
            if s < 3 { res[s] = i; }
        }
        assert_eq!(v, res);
    }

    #[test]
    fn horizon_rules() {
        let mut ds: PositionalDownsampler<u32> = PositionalDownsampler::new(2);
        for (s, r) in [(10, 0), (12, 1), (13, 2), (13, 3), (15, 4)] { ds.push(s, r); }
        // P1 = 12 keeps 1 read -> P2 = 13
        assert_eq!(ds.horizon(Some(10)), Horizon::Through(13));
        // P1 = 13 keeps 2 -> 13
        assert_eq!(ds.horizon(Some(12)), Horizon::Through(13));
        // P1 = 15 is incomplete with 1 read -> unknown
        assert_eq!(ds.horizon(Some(13)), Horizon::NeedInput);
        ds.finish();
        assert_eq!(ds.horizon(Some(13)), Horizon::All);
        assert_eq!(ds.horizon(None), Horizon::All);
        // target 1: always P2
        let mut d1: PositionalDownsampler<u32> = PositionalDownsampler::new(1);
        for (s, r) in [(10, 0), (13, 2), (13, 3), (15, 4), (16, 5)] { d1.push(s, r); }
        assert_eq!(d1.horizon(Some(10)), Horizon::Through(15));
        // finalized groups are still seen by horizon
        let mut rng = JavaRandom::gatk();
        let mut out = Vec::new();
        d1.finalize_upto(15, &mut rng, &mut out);
        assert_eq!(d1.horizon(Some(10)), Horizon::Through(15));
    }
}
