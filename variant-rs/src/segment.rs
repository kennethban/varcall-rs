//! Segments of contigs inside one RNG stream (single-run mode, VCF): a contig is cut into pieces that run in
//! parallel, with the cuts placed away from drawing groups so that most pieces never wait for GATK's generator.
//!
//! In single-run mode the generator is one stream over the whole genome: a drawing group (a start position with
//! more than `--max-reads-per-alignment-start` reads) needs the generator's state after every earlier QD jitter
//! and group draw. With one job per contig, every chromosome reaches its first drawing group within seconds and
//! then waits until the previous chromosome is done: the run is the chromosomes one after another (HG002 30x,
//! 128 threads: 1,064 s against 165 s for 48 shard files). Cutting a contig at points where the cut provably
//! changes nothing (`frontend::SegSpec`) and isolating the clusters of drawing groups in short segments leaves
//! only those short segments waiting for their turn.
//!
//! The drawing groups are located by a pre-scan of the reads' starts (flags and mapping quality only, so a superset
//! of the true groups); the front ends check the cuts themselves and a segment that meets a drawing group near its
//! cut takes over the next segment's range, so the pre-scan only guides where the cuts go.

use crate::frontend::{SEG_WARMUP, SEG_ZONE};
use std::sync::atomic::AtomicU64;
use std::sync::Arc;

/// Margin between a cut and the nearest drawing group found by the pre-scan: room for the previous segment's stop
/// (about 1.2 kb past the cut) and for the next segment's zone (`SEG_ZONE` plus a read).
pub const GROUP_MARGIN: i64 = 5_000;
/// Drawing groups closer than this share a segment (a segment of its own costs more than reading the gap).
pub const CLUSTER_GAP: i64 = 100_000;

/// Start positions of `tid` where more than `cap / n_inputs` reads (of one input) start, among reads that are mapped,
/// primary, not duplicates, not QC-failed, with mapping quality >= 20 and != 255: every drawing group of the
/// downsampler is among them (the HC read filters only remove more reads; with several inputs one of them holds
/// more than its share of a group's reads).
/// Only starts in [from, to] (a chunk of the contig; reads starting before it are skipped).
pub fn scan_groups(inputs: &[crate::input::Input], tid: usize, from: i64, to: i64, cap: usize, io_ns: Arc<AtomicU64>) -> Vec<i64> {
    let thr = (cap / inputs.len().max(1)).max(1);
    let mut out: Vec<i64> = Vec::new();
    for inp in inputs {
        let mut cs = inp.stream(tid, &[(from, to)], true, io_ns.clone());
        let mut body = Vec::new();
        let (mut cur, mut n) = (i64::MIN, 0usize);
        while cs.next_record(&mut body) {
            let (rid, start, flag) = crate::bam::peek_core(&body);
            if rid != tid as i32 {
                if rid > tid as i32 || rid < 0 {
                    break;
                }
                continue;
            }
            if start < from {
                continue;
            }
            if start > to {
                break;
            }
            let mapq = body[9];
            if flag & (0x4 | 0x100 | 0x200 | 0x400) != 0 || mapq < 20 || mapq == 255 {
                continue;
            }
            if start != cur {
                cur = start;
                n = 0;
            }
            n += 1;
            if n == thr + 1 {
                out.push(start);
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

/// One segment of a contig: it owns [from, next cut) (see `frontend::SegSpec`); `groups`: the pre-scan found
/// drawing groups in it (it will wait for its generator turn).
#[derive(Clone, Debug)]
pub struct Seg {
    pub from: Option<i64>,
    pub cuts: Vec<i64>,
    pub groups: bool,
}

/// Cuts of a contig of length `len` with (pre-scan) drawing groups `groups`: each cluster of groups (closer than
/// CLUSTER_GAP) gets a segment of its own from GROUP_MARGIN before its first group to GROUP_MARGIN after its last;
/// the stretches between are cut into pieces of about `piece` bp. Returns the segments in order (one segment: the
/// contig is not cut).
pub fn plan(len: i64, groups: &[i64], piece: i64) -> Vec<Seg> {
    let min_piece = 2 * GROUP_MARGIN;
    let lo = SEG_WARMUP + SEG_ZONE + 2; // a cut leaves room for the warm-up
    let hi = len - GROUP_MARGIN; // and for the previous segment's stop
    if hi - lo < min_piece {
        return vec![Seg { from: None, cuts: vec![], groups: !groups.is_empty() }];
    }
    // clusters, extended by the margin and merged when closer than a minimal piece
    let mut clusters: Vec<(i64, i64)> = Vec::new();
    for &g in groups {
        let (a, b) = (g - GROUP_MARGIN, g + GROUP_MARGIN);
        match clusters.last_mut() {
            Some(c) if a - c.1 < CLUSTER_GAP.max(min_piece) => c.1 = c.1.max(b),
            _ => clusters.push((a, b)),
        }
    }
    // boundaries between segments
    let mut bounds: Vec<i64> = Vec::new();
    let stretch = |bounds: &mut Vec<i64>, from: i64, to: i64| {
        // group-free stretch [from, to): pieces of about `piece`
        let n = ((to - from) as f64 / piece as f64).round().max(1.0) as i64;
        for i in 1..n {
            bounds.push(from + (to - from) * i / n);
        }
    };
    let mut cur = 1i64;
    for &(a, b) in &clusters {
        // the cluster's segment [a, b)
        let (a, b) = (a.max(lo), (b + 1).min(hi));
        if a - cur >= min_piece {
            stretch(&mut bounds, cur, a);
            bounds.push(a);
            cur = a;
        } // else the cluster joins the current segment
        if hi - b < min_piece {
            cur = hi; // the rest of the contig belongs to it
            break;
        }
        if b > cur {
            bounds.push(b);
            cur = b;
        }
    }
    if cur < hi {
        stretch(&mut bounds, cur, len + 1);
    }
    bounds.retain(|&c| c >= lo && c <= hi);
    bounds.dedup();
    let mut segs = Vec::with_capacity(bounds.len() + 1);
    for i in 0..=bounds.len() {
        let from = if i == 0 { 1 } else { bounds[i - 1] };
        let to = bounds.get(i).copied().unwrap_or(len + 1);
        let groups = groups.iter().any(|&g| g >= from - SEG_ZONE && g < to);
        segs.push(Seg { from: if i == 0 { None } else { Some(from) }, cuts: bounds[i..].to_vec(), groups });
    }
    segs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plan_isolates_clusters() {
        let len = 100_000_000;
        let g = vec![10_000_000, 10_005_000, 50_000_000];
        let s = plan(len, &g, 8_000_000);
        // contiguous, ordered
        for w in s.windows(2) {
            assert_eq!(w[0].cuts[0], w[1].from.unwrap());
        }
        // every group inside a segment flagged as holding groups, at least GROUP_MARGIN from its cuts
        for &x in &g {
            let i = s.iter().rposition(|sg| sg.from.is_none_or(|f| f <= x)).unwrap();
            assert!(s[i].groups);
            assert!(s[i].from.is_none_or(|f| x - f >= GROUP_MARGIN));
            assert!(s[i].cuts.first().is_none_or(|&c| c - x >= GROUP_MARGIN));
            assert!(s[i].cuts.first().unwrap() - s[i].from.unwrap_or(1) < 1_000_000);
        }
        assert!(s.iter().filter(|x| !x.groups).count() >= 10);
    }

    #[test]
    fn plan_small_contig_not_cut() {
        assert_eq!(plan(30_000, &[], 1_000_000).len(), 1);
        let s = plan(5_000_000, &[], 8_000_000);
        assert_eq!(s.len(), 1);
    }
}
