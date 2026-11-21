//! Seed extension for one read (bwamem.cpp `mem_chain2aln_across_reads_V2`, `mem_kernel2_core`,
//! `mem_sort_dedup_patch`, `mem_patch_reg`; bandedSWA.cpp `scalarBandedSWA` = ksw_extend2; ksw.cpp `ksw_global2`).
//!
//! bwa-mem2 extends every seed of every chain (left, then right with the left score as h0, each with up to two band
//! widths), then replays bwa's skip test in seed order and purges the regions of seeds bwa would not have extended.

use crate::chain::Chain;
use crate::index::Bns;

pub const H0: i64 = -99;

#[derive(Clone)]
pub struct ExtOpt {
    pub a: i32,
    pub o_del: i32,
    pub e_del: i32,
    pub o_ins: i32,
    pub e_ins: i32,
    pub w: i32,
    pub zdrop: i32,
    pub pen_clip5: i32,
    pub pen_clip3: i32,
    pub mat: [i8; 25],
    pub max_chain_gap: i64,
    pub mask_level_redun: f32,
}

impl ExtOpt {
    pub fn bwa_default() -> ExtOpt {
        let (a, b) = (1i8, 4i8);
        // bwa_fill_scmat
        let mut mat = [0i8; 25];
        let mut k = 0;
        for i in 0..4 {
            for j in 0..4 {
                mat[k] = if i == j { a } else { -b };
                k += 1;
            }
            mat[k] = -1;
            k += 1;
        }
        for _ in 0..5 {
            mat[k] = -1;
            k += 1;
        }
        ExtOpt { a: 1, o_del: 6, e_del: 1, o_ins: 6, e_ins: 1, w: 100, zdrop: 100, pen_clip5: 5, pen_clip3: 5, mat, max_chain_gap: 10000, mask_level_redun: 0.95 }
    }
}

/// mem_alnreg_t
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Reg {
    pub rb: i64,
    pub re: i64,
    pub qb: i32,
    pub qe: i32,
    pub rid: i32,
    pub chain: usize,
    pub score: i32,
    pub truesc: i32,
    pub sub: i32,
    pub alt_sc: i32,
    pub csub: i32,
    pub sub_n: i32,
    pub w: i32,
    pub seedcov: i32,
    pub secondary: i32,
    pub secondary_all: i32,
    pub seedlen0: i32,
    pub n_comp: i32,
    pub is_alt: i32,
    pub frac_rep: f32,
    pub hash: u64,
}

/// cal_max_gap
fn cal_max_gap(o: &ExtOpt, qlen: i32) -> i32 {
    let l_del = ((qlen * o.a - o.o_del) as f64 / o.e_del as f64 + 1.0) as i32;
    let l_ins = ((qlen * o.a - o.o_ins) as f64 / o.e_ins as f64 + 1.0) as i32;
    let l = l_del.max(l_ins).max(1);
    l.min(o.w << 1)
}

pub struct ExtOut {
    pub score: i32,
    pub qle: i32,
    pub tle: i32,
    pub gtle: i32,
    pub gscore: i32,
    pub max_off: i32,
}

/// scalarBandedSWA (= ksw_extend2 with end_bonus)
pub fn ksw_extend2(o: &ExtOpt, query: &[u8], target: &[u8], mut w: i32, end_bonus: i32, h0: i32) -> ExtOut {
    let qlen = query.len() as i32;
    let tlen = target.len() as i32;
    let m = 5usize;
    let (o_del, e_del, o_ins, e_ins) = (o.o_del, o.e_del, o.o_ins, o.e_ins);
    let (oe_del, oe_ins) = (o_del + e_del, o_ins + e_ins);
    // query profile
    let mut qp = vec![0i8; query.len() * m];
    let mut i = 0;
    for k in 0..m {
        let p = &o.mat[k * m..k * m + m];
        for &qj in query {
            qp[i] = p[qj as usize];
            i += 1;
        }
    }
    let mut eh_h = vec![0i32; query.len() + 1];
    let mut eh_e = vec![0i32; query.len() + 1];
    eh_h[0] = h0;
    eh_h[1] = if h0 > oe_ins { h0 - oe_ins } else { 0 };
    let mut j = 2;
    while j <= qlen && eh_h[(j - 1) as usize] > e_ins {
        eh_h[j as usize] = eh_h[(j - 1) as usize] - e_ins;
        j += 1;
    }
    let mut max = o.mat.iter().map(|&x| x as i32).fold(0, |a, b| a.max(b));
    let max_ins = (((qlen * max + end_bonus - o_ins) as f64) / e_ins as f64 + 1.0) as i32;
    let max_ins = max_ins.max(1);
    w = w.min(max_ins);
    let max_del = (((qlen * max + end_bonus - o_del) as f64) / e_del as f64 + 1.0) as i32;
    let max_del = max_del.max(1);
    w = w.min(max_del);
    max = h0;
    let (mut max_i, mut max_j, mut max_ie, mut gscore) = (-1i32, -1i32, -1i32, -1i32);
    let mut max_off = 0i32;
    let (mut beg, mut end) = (0i32, qlen);
    let mut i = 0i32;
    while i < tlen {
        let (mut f, mut h1, mut mm, mut mj) = (0i32, 0i32, 0i32, -1i32);
        let q = &qp[target[i as usize] as usize * query.len()..];
        if beg < i - w {
            beg = i - w;
        }
        if end > i + w + 1 {
            end = i + w + 1;
        }
        if end > qlen {
            end = qlen;
        }
        if beg == 0 {
            h1 = h0 - (o_del + e_del * (i + 1));
            if h1 < 0 {
                h1 = 0;
            }
        }
        let mut j = beg;
        while j < end {
            let ju = j as usize;
            let mut mv = eh_h[ju];
            let mut e = eh_e[ju];
            eh_h[ju] = h1;
            mv = if mv != 0 { mv + q[ju] as i32 } else { 0 };
            let mut h = if mv > e { mv } else { e };
            h = if h > f { h } else { f };
            h1 = h;
            mj = if mm > h { mj } else { j };
            mm = if mm > h { mm } else { h };
            let mut t = mv - oe_del;
            t = if t > 0 { t } else { 0 };
            e -= e_del;
            e = if e > t { e } else { t };
            eh_e[ju] = e;
            let mut t = mv - oe_ins;
            t = if t > 0 { t } else { 0 };
            f -= e_ins;
            f = if f > t { f } else { t };
            j += 1;
        }
        eh_h[end as usize] = h1;
        eh_e[end as usize] = 0;
        if j == qlen {
            max_ie = if gscore > h1 { max_ie } else { i };
            gscore = if gscore > h1 { gscore } else { h1 };
        }
        if mm == 0 {
            break;
        }
        if mm > max {
            max = mm;
            max_i = i;
            max_j = mj;
            max_off = max_off.max((mj - i).abs());
        } else if o.zdrop > 0 {
            if i - max_i > mj - max_j {
                if max - mm - ((i - max_i) - (mj - max_j)) * e_del > o.zdrop {
                    break;
                }
            } else if max - mm - ((mj - max_j) - (i - max_i)) * e_ins > o.zdrop {
                break;
            }
        }
        let mut j = beg;
        while j < end && eh_h[j as usize] == 0 && eh_e[j as usize] == 0 {
            j += 1;
        }
        beg = j;
        let mut j = end;
        while j >= beg && eh_h[j as usize] == 0 && eh_e[j as usize] == 0 {
            j -= 1;
        }
        end = if j + 2 < qlen { j + 2 } else { qlen };
        i += 1;
    }
    ExtOut { score: max, qle: max_j + 1, tle: max_i + 1, gtle: max_ie + 1, gscore, max_off }
}

const MINUS_INF: i32 = -0x4000_0000;

/// ksw_global2 without a CIGAR (score only)
pub fn ksw_global2_score(o: &ExtOpt, query: &[u8], target: &[u8], w: i32) -> i32 {
    let qlen = query.len() as i32;
    let tlen = target.len() as i32;
    let m = 5usize;
    let (oe_del, oe_ins) = (o.o_del + o.e_del, o.o_ins + o.e_ins);
    let mut qp = vec![0i8; query.len() * m];
    let mut i = 0;
    for k in 0..m {
        let p = &o.mat[k * m..k * m + m];
        for &qj in query {
            qp[i] = p[qj as usize];
            i += 1;
        }
    }
    let mut eh_h = vec![0i32; query.len() + 1];
    let mut eh_e = vec![0i32; query.len() + 1];
    eh_h[0] = 0;
    eh_e[0] = MINUS_INF;
    let mut j = 1;
    while j <= qlen && j <= w {
        eh_h[j as usize] = -(o.o_ins + o.e_ins * j);
        eh_e[j as usize] = MINUS_INF;
        j += 1;
    }
    while j <= qlen {
        eh_h[j as usize] = MINUS_INF;
        eh_e[j as usize] = MINUS_INF;
        j += 1;
    }
    for i in 0..tlen {
        let mut f = MINUS_INF;
        let q = &qp[target[i as usize] as usize * query.len()..];
        let beg = if i > w { i - w } else { 0 };
        let end = if i + w + 1 < qlen { i + w + 1 } else { qlen };
        let mut h1 = if beg == 0 { -(o.o_del + o.e_del * (i + 1)) } else { MINUS_INF };
        for j in beg..end {
            let ju = j as usize;
            let mut mv = eh_h[ju];
            let mut e = eh_e[ju];
            eh_h[ju] = h1;
            mv += q[ju] as i32;
            let mut h = if mv >= e { mv } else { e };
            h = if h >= f { h } else { f };
            h1 = h;
            let t = mv - oe_del;
            e -= o.e_del;
            e = if e > t { e } else { t };
            eh_e[ju] = e;
            let t = mv - oe_ins;
            f -= o.e_ins;
            f = if f > t { f } else { t };
        }
        eh_h[end as usize] = h1;
        eh_e[end as usize] = MINUS_INF;
    }
    eh_h[qlen as usize]
}

/// bns_get_seq from the 2-bit pac (forward bases; reverse strand complemented and reversed)
pub fn bns_get_seq(bns: &Bns, mut beg: i64, mut end: i64) -> Vec<u8> {
    let l_pac = bns.l_pac;
    if end < beg {
        std::mem::swap(&mut beg, &mut end);
    }
    if end > l_pac << 1 {
        end = l_pac << 1;
    }
    if beg < 0 {
        beg = 0;
    }
    let mut seq = Vec::new();
    if beg >= l_pac || end <= l_pac {
        seq.reserve((end - beg).max(0) as usize);
        let pac = bns.pac_bytes();
        let base = |i: i64| (pac[(i >> 2) as usize] >> ((!(i as u64) & 3) << 1)) & 3;
        if beg >= l_pac {
            let beg_f = (l_pac << 1) - 1 - end;
            let end_f = (l_pac << 1) - 1 - beg;
            seq.extend(((beg_f + 1)..=end_f).rev().map(|k| 3 - base(k)));
        } else {
            seq.extend((beg..end).map(base));
        }
    }
    seq
}

/// bwa_gen_cigar2's score (no CIGAR): global alignment of query against [rb, re)
pub fn gen_cigar2_score(o: &ExtOpt, bns: &Bns, w_: i32, query: &[u8], rb: i64, re: i64) -> Option<i32> {
    let l_pac = bns.l_pac;
    let l_query = query.len() as i32;
    if l_query <= 0 || rb >= re || (rb < l_pac && re > l_pac) {
        return None;
    }
    let mut rseq = bns_get_seq(bns, rb, re);
    if re - rb != rseq.len() as i64 {
        return None;
    }
    let mut q = query.to_vec();
    if rb >= l_pac {
        q.reverse();
        rseq.reverse();
    }
    let rlen = rseq.len() as i32;
    if l_query == (re - rb) as i32 && w_ == 0 {
        return Some((0..query.len()).map(|i| o.mat[rseq[i] as usize * 5 + q[i] as usize] as i32).sum());
    }
    let max_ins = ((((l_query + 1) >> 1) * o.mat[0] as i32 - o.o_ins) as f64 / o.e_ins as f64 + 1.0) as i32;
    let max_del = ((((l_query + 1) >> 1) * o.mat[0] as i32 - o.o_del) as f64 / o.e_del as f64 + 1.0) as i32;
    let max_gap = max_ins.max(max_del).max(1);
    let mut w = (max_gap + (rlen - l_query).abs() + 1) >> 1;
    w = w.min(w_);
    let min_w = (rlen - l_query).abs() + 3;
    w = w.max(min_w);
    Some(ksw_global2_score(o, &q, &rseq, w))
}

/// bns_fetch_seq_v2: clamp [beg, end) to the reference sequence holding `mid`; returns (beg, end, rid)
/// `rid`: the reference of `mid` (the chain's: bns_pos2rid of its first seed, as mem_chain2aln recomputes it)
pub(crate) fn fetch_window(bns: &Bns, mut beg: i64, mid: i64, mut end: i64, rid: i32) -> (i64, i64, i32) {
    if end < beg {
        std::mem::swap(&mut beg, &mut end);
    }
    assert!(beg <= mid && mid < end);
    let (pf, is_rev) = bns.depos(mid);
    debug_assert_eq!(rid, bns.pos2rid(pf));
    let mut far_beg = bns.anns[rid as usize].offset;
    let mut far_end = far_beg + bns.anns[rid as usize].len as i64;
    if is_rev {
        let tmp = far_beg;
        far_beg = (bns.l_pac << 1) - far_end;
        far_end = (bns.l_pac << 1) - tmp;
    }
    beg = beg.max(far_beg);
    end = end.min(far_end);
    // bns_get_seq_v2: clipped to [0, 2 l_pac); a window bridging the strands would have length 0
    let (b2, e2) = (beg.max(0), end.min(bns.l_pac << 1));
    assert!((b2 >= bns.l_pac || e2 <= bns.l_pac) && e2 - b2 == end - beg, "reference window bridges the strands");
    (beg, end, rid)
}

fn seedcov(r: &Reg, c: &Chain) -> i32 {
    let mut cov = 0;
    for t in &c.seeds {
        if t.qbeg >= r.qb && t.qbeg + t.len <= r.qe && t.rbeg >= r.rb && t.rbeg + t.len as i64 <= r.re {
            cov += t.len;
        }
    }
    cov
}

const MAX_BAND_TRY: i32 = 2;

/// ALIGN_RS_PROFILE: wall time of the extension phases (ns): setup, job build, run, apply, replay+dedup; jobs run;
/// jobs the CPU extender took by the ungapped shortcut
pub static EPROF: [std::sync::atomic::AtomicU64; 7] = [const { std::sync::atomic::AtomicU64::new(0) }; 7];
fn eadd(i: usize, t: std::time::Instant) {
    EPROF[i].fetch_add(t.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
}

/// One banded extension (ksw_extend2) job: the query and the reference are read from `qstart` / `rstart` in
/// the round's direction for `qlen` / `rlen` bases. `rstart` is a BWT coordinate (both strands, as `.0123`).
/// 24 bytes (it is uploaded per round): what is the same for a whole round (band, end bonus, direction) is in
/// `Round`.
#[derive(Clone, Copy, Debug, Default)]
#[repr(C)]
pub struct Job {
    pub rstart: i64,
    pub read: u32,
    pub rlen: u32,
    pub qstart: u16,
    pub qlen: u16,
    pub h0: u16,
    pub pad: u16,
}

/// What a round's jobs share: band width, end bonus (pen_clip5 left, pen_clip3 right), direction (-1 left)
#[derive(Clone, Copy, Debug)]
#[repr(C)]
pub struct Round {
    pub w: i32,
    pub end_bonus: i32,
    pub dir: i32,
}

impl Job {
    fn new(read: u32, qstart: i32, qlen: i32, rlen: i64, rstart: i64, h0: i32) -> Job {
        let u16_ = |v: i32, what: &str| u16::try_from(v).unwrap_or_else(|_| panic!("extension job: {what} {v} out of range"));
        Job { rstart, read, rlen: u32::try_from(rlen).expect("extension job: reference length"), qstart: u16_(qstart, "query start"), qlen: u16_(qlen, "query length"), h0: u16_(h0, "h0"), pad: 0 }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[repr(C)]
pub struct JobOut {
    pub score: i32,
    pub qle: i32,
    pub tle: i32,
    pub gtle: i32,
    pub gscore: i32,
    pub max_off: i32,
}

/// ksw_extend2's result read off the main diagonal, when no path with a gap can change it; None otherwise (the job
/// runs the DP). The shortcut of bwa-mem3 (PRs #540, #544), with its conditions derived here for ksw_extend2:
/// * `loss` = sum over the diagonal's cells of (a - score), a the matrix maximum (a mismatch a + b, an N a + 1).
///   A path with a gap reaching row i (or column j) pairs at most i + 1 (j + 1) bases and pays at least
///   min(o_del + e_del, o_ins + e_ins) = `oe`, so it scores at most h0 + (i + 1) a - oe, while the diagonal cell
///   scores h0 + (i + 1) a - loss(0..i). With loss < oe over the whole query, every row's maximum and every
///   column's best cell are the diagonal's, strictly: the row maximum's last position (mj) is i, a new maximum is
///   a diagonal cell (max_off 0), and gscore is H(qlen-1, qlen-1) at row qlen-1 (later rows' column qlen-1 cells are
///   gapped, lower).
/// * the band always holds the diagonal (w >= 1, and the band's trimming keeps columns up to a non-zero cell + 1),
///   and cutting paths (band, zero cells) only lowers the gapped values.
/// * the walk must reach row qlen - 1 as the DP does: tlen >= qlen, every diagonal value > 0 (else `mm == 0` or a
///   dead cell), and z-drop not firing on a diagonal row (max - D(i) > zdrop; the rows' offsets cancel on the
///   diagonal). Rows past the query's end are gapped: they change nothing.
/// `q(k)` / `t(k)`: the query's / target's k-th base in the job's direction. The GPU runs the same walk in
/// `k_extkey` (cuda/aln_ext.cu); here it serves the CPU extender. ALIGN_RS_UNGAPPED=0 turns both off.
#[inline]
pub fn ungapped_extend(o: &ExtOpt, qlen: i32, tlen: i64, h0: i32, q: impl Fn(i32) -> u8, t: impl Fn(i32) -> u8) -> Option<JobOut> {
    if qlen < 1 || tlen < qlen as i64 {
        return None;
    }
    let a = o.mat.iter().copied().max().unwrap() as i32;
    let oe = (o.o_del + o.e_del).min(o.o_ins + o.e_ins);
    let (mut d, mut max, mut max_i, mut loss) = (h0, h0, -1i32, 0i32);
    for i in 0..qlen {
        let sc = o.mat[q(i) as usize * 5 + t(i) as usize] as i32;
        loss += a - sc;
        if loss >= oe {
            return None;
        }
        d += sc;
        if d <= 0 {
            return None;
        }
        if d > max {
            max = d;
            max_i = i;
        } else if o.zdrop > 0 && max - d > o.zdrop {
            return None;
        }
    }
    Some(JobOut { score: max, qle: max_i + 1, tle: max_i + 1, gtle: qlen, gscore: d, max_off: 0 })
}

fn ungapped_on() -> bool {
    static M: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *M.get_or_init(|| std::env::var("ALIGN_RS_UNGAPPED").as_deref() != Ok("0"))
}

/// Runs extension jobs (CPU: `CpuExtender`; GPU: crate::gpu)
pub trait Extender: Sync {
    /// `key`: the batch (gpu::batch_key), for the device's copy of its reads
    fn run(&self, o: &ExtOpt, rd: Round, jobs: &[Job], reads: &[Vec<u8>], key: u64) -> Vec<JobOut>;
}

pub struct CpuExtender<'a> {
    pub ref_string: &'a [u8],
}

impl Extender for CpuExtender<'_> {
    fn run(&self, o: &ExtOpt, rd: Round, jobs: &[Job], reads: &[Vec<u8>], _key: u64) -> Vec<JobOut> {
        use rayon::prelude::*;
        let dir = rd.dir as i64;
        jobs.par_iter()
            .with_min_len(256)
            .map(|j| {
                let q = &reads[j.read as usize];
                if ungapped_on() {
                    let r = ungapped_extend(o, j.qlen as i32, j.rlen as i64, j.h0 as i32, |k| q[(j.qstart as i64 + k as i64 * dir) as usize], |k| self.ref_string[(j.rstart + k as i64 * dir) as usize]);
                    if let Some(r) = r {
                        EPROF[6].fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        return r;
                    }
                }
                let qs: Vec<u8> = (0..j.qlen as i64).map(|i| q[(j.qstart as i64 + i * dir) as usize]).collect();
                let rs: Vec<u8> = (0..j.rlen as i64).map(|i| self.ref_string[(j.rstart + i * dir) as usize]).collect();
                let r = ksw_extend2(o, &qs, &rs, rd.w, rd.end_bonus, j.h0 as i32);
                JobOut { score: r.score, qle: r.qle, tle: r.tle, gtle: r.gtle, gscore: r.gscore, max_off: r.max_off }
            })
            .collect()
    }
}

/// Per read between the setup and the replay
struct ReadExt {
    av: Vec<Reg>,
    srts: Vec<Vec<u32>>,
    seed_aln: Vec<Vec<usize>>,
    /// per region: the left / right extension template (w and h0 filled per round)
    left: Vec<Option<Job>>,
    right: Vec<Option<Job>>,
    /// ALIGN_RS_PROFILE: extension jobs run per region
    nj: Vec<u8>,
}

fn setup_read(o: &ExtOpt, bns: &Bns, chains: &[Chain], read: u32, l_query: i32) -> ReadExt {
    let l_pac = bns.l_pac;
    let n_seeds: usize = chains.iter().map(|c| c.seeds.len()).sum();
    let mut st = ReadExt { av: Vec::with_capacity(n_seeds), srts: Vec::with_capacity(chains.len()), seed_aln: Vec::with_capacity(chains.len()), left: Vec::with_capacity(n_seeds), right: Vec::with_capacity(n_seeds), nj: Vec::new() };
    for (j, c) in chains.iter().enumerate() {
        let n = c.seeds.len();
        st.seed_aln.push(vec![usize::MAX; n]);
        if n == 0 {
            st.srts.push(Vec::new());
            continue;
        }
        let (mut rmax0, mut rmax1) = (l_pac << 1, 0i64);
        for t in &c.seeds {
            let b = t.rbeg - (t.qbeg as i64 + cal_max_gap(o, t.qbeg) as i64);
            let rest = l_query - t.qbeg - t.len;
            let e = t.rbeg + t.len as i64 + (rest as i64 + cal_max_gap(o, rest) as i64);
            rmax0 = rmax0.min(b);
            rmax1 = rmax1.max(e);
        }
        rmax0 = rmax0.max(0);
        rmax1 = rmax1.min(l_pac << 1);
        if rmax0 < l_pac && l_pac < rmax1 {
            if c.seeds[0].rbeg < l_pac {
                rmax1 = l_pac;
            } else {
                rmax0 = l_pac;
            }
        }
        let (rmax0, rmax1, rid) = crate::prof::cpu(crate::prof::C_FETCHW, || fetch_window(bns, rmax0, c.seeds[0].rbeg, rmax1, c.rid));
        assert_eq!(rid, c.rid);
        let mut srt: Vec<u64> = (0..n).map(|i| ((c.seeds[i].score as u64) << 32) | i as u64).collect();
        srt.sort_unstable();
        let srt: Vec<u32> = srt.iter().map(|&x| x as u32).collect();
        for k in (0..n).rev() {
            let si = srt[k] as usize;
            let s = c.seeds[si];
            st.seed_aln[j][si] = st.av.len();
            let mut a = Reg { w: o.w, score: -1, truesc: -1, rid: c.rid, frac_rep: c.frac_rep, seedlen0: s.len, chain: j, rb: H0, re: H0, qb: H0 as i32, qe: H0 as i32, ..Default::default() };
            let left = if s.qbeg != 0 {
                a.qb = s.qbeg;
                a.rb = s.rbeg;
                // query[qbeg-1 .. 0], reference rseq[tmp-1 .. 0] = ref[rbeg-1 ..] backwards
                Some(Job::new(read, s.qbeg - 1, s.qbeg, s.rbeg - rmax0, s.rbeg - 1, s.len * o.a))
            } else {
                a.score = s.len * o.a;
                a.truesc = a.score;
                a.qb = 0;
                a.rb = s.rbeg;
                None
            };
            let right = if s.qbeg + s.len != l_query {
                let qe = s.qbeg + s.len;
                let re = s.rbeg + s.len as i64 - rmax0;
                a.qe = qe;
                a.re = rmax0 + re;
                Some(Job::new(read, qe, l_query - qe, rmax1 - rmax0 - re, rmax0 + re, 0))
            } else {
                a.qe = l_query;
                a.re = s.rbeg + s.len as i64;
                None
            };
            st.av.push(a);
            st.left.push(left);
            st.right.push(right);
        }
        st.srts.push(srt);
    }
    if crate::prof::on() {
        st.nj = vec![0; st.av.len()];
    }
    st
}

/// One part of a batch (reads base..base + sts.len()) going through the extension rounds: left at w, left retried
/// at 2w, right at w with the left score as h0, right retried, as bwa-mem2 batches them.
struct Part {
    base: usize,
    sts: Vec<ReadExt>,
    side: usize,
    tr: i32,
    /// (read in the part, region) of the jobs of the current round
    pending: Vec<(u32, u32)>,
    /// the current round's jobs (sent to the device)
    jobs: Vec<Job>,
    done: bool,
}

impl Part {
    fn new(o: &ExtOpt, bns: &Bns, chains: &[Vec<Chain>], reads: &[Vec<u8>], base: usize, end: usize) -> Part {
        use rayon::prelude::*;
        let t = std::time::Instant::now();
        let sts: Vec<ReadExt> = (base..end).into_par_iter().map(|i| crate::prof::cpu(crate::prof::C_ESETUP, || setup_read(o, bns, &chains[i], i as u32, reads[i].len() as i32))).collect();
        eadd(0, t);
        let mut p = Part { base, sts, side: 0, tr: 0, pending: Vec::new(), jobs: Vec::new(), done: false };
        p.start_side();
        p
    }

    /// the regions with a template for the current side
    fn start_side(&mut self) {
        loop {
            if self.side == 2 {
                self.done = true;
                return;
            }
            self.tr = 0;
            self.pending.clear();
            let side = self.side;
            for (ri, st) in self.sts.iter_mut().enumerate() {
                for k in 0..st.av.len() {
                    if side == 1 {
                        // the right extension's h0: the score after the left extension (or the seed's), set once
                        // for all its rounds, as bwa-mem2 sets sp->h0 before the right rounds; a retry must not
                        // start from the previous try's score
                        let sc = st.av[k].score;
                        if let Some(j) = st.right[k].as_mut() {
                            j.h0 = u16::try_from(sc).expect("extension job: h0 out of range");
                        }
                    }
                    if (if side == 0 { &st.left[k] } else { &st.right[k] }).is_some() {
                        self.pending.push((ri as u32, k as u32));
                    }
                }
            }
            if !self.pending.is_empty() {
                return;
            }
            self.side += 1;
        }
    }

    /// what this round's jobs share
    fn round(&self, o: &ExtOpt) -> Round {
        Round { w: o.w << self.tr, end_bonus: if self.side == 0 { o.pen_clip5 } else { o.pen_clip3 }, dir: if self.side == 0 { -1 } else { 1 } }
    }

    /// this round's jobs (empty when the part is done)
    fn build(&mut self, _o: &ExtOpt) -> Vec<Job> {
        use rayon::prelude::*;
        if self.done {
            return Vec::new();
        }
        let t = std::time::Instant::now();
        let (side, sts) = (self.side, &self.sts);
        let jobs: Vec<Job> = self
            .pending
            .par_iter()
            .with_min_len(1024)
            .map(|&(ri, k)| {
                let st = &sts[ri as usize];
                // side 1: h0 was fixed when the side started (start_side)
                if side == 0 { st.left[k as usize] } else { st.right[k as usize] }.unwrap()
            })
            .collect();
        eadd(1, t);
        jobs
    }

    /// applies this round's results; the next round's regions (or the next side's, or done)
    fn apply(&mut self, o: &ExtOpt, reads: &[Vec<u8>], outs: &[JobOut]) {
        use rayon::prelude::*;
        let t = std::time::Instant::now();
        let (w, side, tr, base) = (o.w << self.tr, self.side, self.tr, self.base);
        let (pending, jobs) = (&self.pending, &self.jobs);
        // results applied per read in parallel: `pending` is ordered by read, each read's jobs contiguous
        let mut ranges: Vec<(usize, usize)> = Vec::new();
        let mut st0 = 0usize;
        for idx in 1..=pending.len() {
            if idx == pending.len() || pending[idx].0 != pending[st0].0 {
                ranges.push((st0, idx));
                st0 = idx;
            }
        }
        let next_parts: Vec<Vec<(u32, u32)>> = {
            #[derive(Clone, Copy)]
            struct P(*mut ReadExt);
            unsafe impl Send for P {}
            unsafe impl Sync for P {}
            let bp = P(self.sts.as_mut_ptr());
            ranges
                .par_iter()
                .map(|&(a0, a1)| crate::prof::cpu(crate::prof::C_EAPPLY, || {
                    let b = bp;
                    let ri = pending[a0].0 as usize;
                    // SAFETY: one range per read: no two tasks touch the same ReadExt
                    let st = unsafe { &mut *b.0.add(ri) };
                    let mut nx = Vec::new();
                    for idx in a0..a1 {
                        let k = pending[idx].1;
                        if !st.nj.is_empty() {
                            st.nj[k as usize] += 1;
                        }
                        let r = outs[idx];
                        let h0 = jobs[idx].h0 as i32;
                        let a = &mut st.av[k as usize];
                        let prev = a.score;
                        a.score = r.score;
                        if a.score == prev || r.max_off < (w >> 1) + (w >> 2) || tr + 1 == MAX_BAND_TRY {
                            if side == 0 {
                                if r.gscore <= 0 || r.gscore <= a.score - o.pen_clip5 {
                                    a.qb -= r.qle;
                                    a.rb -= r.tle as i64;
                                    a.truesc = a.score;
                                } else {
                                    a.qb = 0;
                                    a.rb -= r.gtle as i64;
                                    a.truesc = r.gscore;
                                }
                            } else if r.gscore <= 0 || r.gscore <= a.score - o.pen_clip3 {
                                a.qe += r.qle;
                                a.re += r.tle as i64;
                                a.truesc += a.score - h0;
                            } else {
                                a.qe = reads[base + ri].len() as i32;
                                a.re += r.gtle as i64;
                                a.truesc += r.gscore - h0;
                            }
                            a.w = a.w.max(w);
                        } else {
                            nx.push((ri as u32, k));
                        }
                    }
                    nx
                }))
                .collect()
        };
        self.pending = next_parts.into_iter().flatten().collect();
        self.tr += 1;
        if self.pending.is_empty() || self.tr == MAX_BAND_TRY {
            self.side += 1;
            self.start_side();
        }
        eadd(3, t);
    }

    /// the regions after the skip replay, the purge and mem_sort_dedup_patch
    fn finish(self, o: &ExtOpt, bns: &Bns, chains: &[Vec<Chain>], reads: &[Vec<u8>]) -> Vec<Vec<Reg>> {
        use rayon::prelude::*;
        let t = std::time::Instant::now();
        let base = self.base;
        let r = self
            .sts
            .into_par_iter()
            .enumerate()
            .map(|(i, st)| crate::prof::cpu(crate::prof::C_EREPLAY, || {
                let ri = base + i;
                let ReadExt { mut av, mut srts, seed_aln, nj, .. } = st;
                let chains = &chains[ri];
                for a in av.iter_mut() {
                    if a.rb != H0 && a.qb as i64 != H0 && a.qe as i64 != H0 && a.re != H0 {
                        a.seedcov = seedcov(a, &chains[a.chain]);
                    }
                }
                let mut v = replay_and_dedup(o, bns, chains, &reads[ri], av, &mut srts, &seed_aln, &nj);
                // the regions live on through mate rescue and SAM: drop the pre-dedup capacity
                v.shrink_to_fit();
                v
            }))
            .collect();
        eadd(4, t);
        r
    }
}

/// Extension of a batch: per read, the regions after the skip replay, the purge and mem_sort_dedup_patch.
/// `chains[i]` are read i's filtered chains. The batch goes through in two halves, alternating: while the
/// device runs one half's round, the CPU sets up, applies or replays the other (reads are independent, so the
/// regions are those of one pass over the whole batch).
/// `pool`: the thread pool for the CPU steps' parallel loops (None: the global one). The device calls stay outside it
/// (on their own thread): injected into a busy pool, a call would wait for its workers to finish their CPU work first.
pub fn extend_batch(o: &ExtOpt, bns: &Bns, ex: &dyn Extender, chains: &[Vec<Chain>], reads: &[Vec<u8>], key: u64, pool: Option<&rayon::ThreadPool>) -> Vec<Vec<Reg>> {
    match pool {
        Some(p) => p.install(|| extend_halves(o, bns, ex, chains, reads, key)),
        None => extend_halves(o, bns, ex, chains, reads, key),
    }
}

fn extend_halves(o: &ExtOpt, bns: &Bns, ex: &dyn Extender, chains: &[Vec<Chain>], reads: &[Vec<u8>], key: u64) -> Vec<Vec<Reg>> {
    let n = chains.len();
    let mid = n / 2;
    let run = |rd: Round, jobs: &[Job]| -> Vec<JobOut> {
        if jobs.is_empty() {
            return Vec::new();
        }
        let t = std::time::Instant::now();
        let r = ex.run(o, rd, jobs, reads, key);
        eadd(2, t);
        EPROF[5].fetch_add(jobs.len() as u64, std::sync::atomic::Ordering::Relaxed);
        r
    };
    let mut a = Part::new(o, bns, chains, reads, 0, mid);
    a.jobs = a.build(o);
    let mut b: Option<Part> = None;
    let mut b_outs: Vec<JobOut> = Vec::new();
    let mut a_res: Option<Vec<Vec<Reg>>> = None;
    let mut b_res: Option<Vec<Vec<Reg>>> = None;
    loop {
        // device: a's round; CPU: b's setup, or the results of b's last round and its next jobs (or its replay)
        let (a_outs, ()) = std::thread::scope(|sc| {
            let (a_round, a_jobs) = (a.round(o), &a.jobs);
            let h = sc.spawn(move || run(a_round, a_jobs));
            if b_res.is_none() {
                match b.as_mut() {
                    None => {
                        let mut p = Part::new(o, bns, chains, reads, mid, n);
                        p.jobs = p.build(o);
                        b = Some(p);
                    }
                    Some(p) => {
                        if !p.jobs.is_empty() {
                            p.apply(o, reads, &b_outs);
                            p.jobs = p.build(o);
                        }
                        if p.done {
                            b_res = Some(b.take().unwrap().finish(o, bns, chains, reads));
                        }
                    }
                }
            }
            (h.join().unwrap(), ())
        });
        // device: b's round; CPU: a's results and next jobs (or its replay)
        let b_jobs: Vec<Job> = b.as_mut().map(|p| std::mem::take(&mut p.jobs)).unwrap_or_default();
        let b_round = b.as_ref().map_or(Round { w: 0, end_bonus: 0, dir: 1 }, |p| p.round(o));
        let (outs, ()) = std::thread::scope(|sc| {
            let h = sc.spawn(|| run(b_round, &b_jobs));
            if a_res.is_none() {
                if !a.jobs.is_empty() {
                    a.apply(o, reads, &a_outs);
                    a.jobs = a.build(o);
                }
                if a.done {
                    let p = std::mem::replace(&mut a, Part { base: 0, sts: Vec::new(), side: 2, tr: 0, pending: Vec::new(), jobs: Vec::new(), done: true });
                    a_res = Some(p.finish(o, bns, chains, reads));
                }
            }
            (h.join().unwrap(), ())
        });
        b_outs = outs;
        if let Some(p) = b.as_mut() {
            p.jobs = b_jobs;
        }
        if a_res.is_some() && b_res.is_some() {
            break;
        }
    }
    let mut r = a_res.unwrap();
    r.extend(b_res.unwrap());
    r
}

/// ALIGN_RS_REPLAY_CHECK=1: the replay's position lookup also runs bwa's linear scan and must agree (a check)
fn replay_check() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("ALIGN_RS_REPLAY_CHECK").is_ok_and(|v| v == "1"))
}

fn replay_and_dedup(o: &ExtOpt, bns: &Bns, chains: &[Chain], query: &[u8], mut av: Vec<Reg>, srts: &mut [Vec<u32>], seed_aln: &[Vec<usize>], nj: &[u8]) -> Vec<Reg> {
    let l_query = query.len() as i32;
    let t_skip = if crate::prof::fine() { crate::prof::tcpu() } else { 0 };
    // replay bwa's skip test: seeds in the same order; a seed (almost) inside an earlier kept region, with no
    // overlapping longer seed leading elsewhere, loses its region. bwa scans the kept regions in order and stops at
    // the first that contains the seed, so what counts is whether any does: the kept regions are kept sorted by
    // start, and only those starting within the longest kept region's length before the seed are tested.
    let mut lim = 0i32;
    let mut kept: Vec<(i64, u32)> = Vec::new(); // (rb, index in av) of the kept regions, by rb
    let mut max_len = 0i64;
    for (j, c) in chains.iter().enumerate() {
        let n = c.seeds.len();
        let srt2 = &mut srts[j];
        for k in (0..n).rev() {
            let s = c.seeds[srt2[k] as usize];
            // a region containing the seed starts in [seed end - longest kept region, seed start]; a seed longer
            // than every kept region has none (the window is empty: lo can pass hi)
            let lo = kept.partition_point(|x| x.0 < s.rbeg + s.len as i64 - max_len);
            let hi = kept.partition_point(|x| x.0 <= s.rbeg);
            let test = |p: &Reg| {
                if s.rbeg < p.rb || s.rbeg + s.len as i64 > p.re || s.qbeg < p.qb || s.qbeg + s.len > p.qe {
                    return false;
                }
                if (s.len - p.seedlen0) as f64 > 0.1 * l_query as f64 {
                    return false;
                }
                let qd = (s.qbeg - p.qb) as i64;
                let rd = s.rbeg - p.rb;
                let max_gap = cal_max_gap(o, qd.min(rd) as i32);
                let w = max_gap.min(p.w) as i64;
                if qd - rd < w && rd - qd < w {
                    return true;
                }
                let qd = (p.qe - (s.qbeg + s.len)) as i64;
                let rd = p.re - (s.rbeg + s.len as i64);
                let max_gap = cal_max_gap(o, qd.min(rd) as i32);
                let w = max_gap.min(p.w) as i64;
                qd - rd < w && rd - qd < w
            };
            let contained = lo < hi && kept[lo..hi].iter().any(|&(_, pi)| test(&av[pi as usize]));
            if replay_check() {
                // bwa's scan: the first `lim` regions not purged, in order
                let (mut v, mut i, mut scan) = (0i32, 0usize, false);
                while i < av.len() && v < lim {
                    let p = &av[i];
                    i += 1;
                    if p.qb == -1 && p.qe == -1 {
                        continue;
                    }
                    if test(p) {
                        scan = true;
                        break;
                    }
                    v += 1;
                }
                assert_eq!(contained, scan, "replay lookup differs from bwa's scan");
            }
            let v = if contained { 0 } else { lim };
            if v < lim {
                let mut vv = k + 1;
                while vv < n {
                    if srt2[vv] == u32::MAX {
                        vv += 1;
                        continue;
                    }
                    let t = c.seeds[srt2[vv] as usize];
                    if (t.len as f64) < s.len as f64 * 0.95 {
                        vv += 1;
                        continue;
                    }
                    if s.qbeg <= t.qbeg && s.qbeg + s.len - t.qbeg >= s.len >> 2 && (t.qbeg - s.qbeg) as i64 != t.rbeg - s.rbeg {
                        break;
                    }
                    if t.qbeg <= s.qbeg && t.qbeg + t.len - s.qbeg >= s.len >> 2 && (s.qbeg - t.qbeg) as i64 != s.rbeg - t.rbeg {
                        break;
                    }
                    vv += 1;
                }
                if vv == n {
                    let ai = seed_aln[j][srt2[k] as usize];
                    av[ai].qb = -1;
                    av[ai].qe = -1;
                    srt2[k] = u32::MAX;
                    continue;
                }
            }
            lim += 1;
            let ai = seed_aln[j][srt2[k] as usize];
            let p = &av[ai];
            let at = kept.partition_point(|x| x.0 <= p.rb);
            kept.insert(at, (p.rb, ai as u32));
            max_len = max_len.max(p.re - p.rb);
        }
    }
    if !nj.is_empty() {
        use crate::prof::*;
        let kept: Vec<usize> = (0..av.len()).filter(|&i| !(av[i].qb == -1 && av[i].qe == -1)).collect();
        count(N_REG, av.len() as u64);
        count(N_KEPT, kept.len() as u64);
        count(N_JOBS, nj.iter().map(|&x| x as u64).sum());
        count(N_KJOBS, kept.iter().map(|&i| nj[i] as u64).sum());
        let h = match kept.len() {
            0..=4 => kept.len(),
            5..=8 => 5,
            9..=16 => 6,
            _ => 7,
        };
        count(N_HIST + h, 1);
    }
    if crate::prof::fine() {
        crate::prof::add_cpu(crate::prof::C_RSKIP, crate::prof::tcpu() - t_skip);
    }
    av.retain(|a| a.qe > a.qb);
    let mut av = crate::prof::cpu(crate::prof::C_RDEDUP, || sort_dedup_patch(o, bns, query, av));
    for r in av.iter_mut() {
        if r.rid >= 0 && bns.anns[r.rid as usize].is_alt {
            r.is_alt = 1;
        }
    }
    av
}

const PATCH_MAX_R_BW: f32 = 0.05;
const PATCH_MIN_SC_RATIO: f32 = 0.90;

/// mem_patch_reg: the merged score of a (left) and b (right), 0 if they should not be merged
fn patch_reg(o: &ExtOpt, bns: &Bns, query: &[u8], a: &Reg, b: &Reg) -> (i32, i32) {
    assert!(a.rid == b.rid && a.rb <= b.rb);
    if a.rb < bns.l_pac && b.rb >= bns.l_pac {
        return (0, 0);
    }
    if a.qb >= b.qb || a.qe >= b.qe || a.re >= b.re {
        return (0, 0);
    }
    let mut w = ((a.re - b.rb) - (a.qe - b.qb) as i64) as i32;
    w = w.abs();
    let mut r = (a.re - b.rb) as f64 / (b.re - a.rb) as f64 - (a.qe - b.qb) as f64 / (b.qe - a.qb) as f64;
    r = if r > 0.0 { r } else { -r };
    if a.re < b.rb || a.qe < b.qb {
        if w > o.w << 1 || r >= PATCH_MAX_R_BW as f64 {
            return (0, 0);
        }
    } else if w > o.w << 2 || r >= (PATCH_MAX_R_BW * 2.0) as f64 {
        return (0, 0);
    }
    w += a.w + b.w;
    w = w.min(o.w << 2);
    let score = match gen_cigar2_score(o, bns, w, &query[a.qb as usize..b.qe as usize], a.rb, b.re) {
        Some(s) => s,
        None => 0, // bwa_gen_cigar2 leaves *score unset when it rejects; mem_patch_reg's `score` is then
                   // uninitialised in bwa-mem2 (never seen: both regions are on one strand here)
    };
    let q_s = ((b.qe - a.qb) as f64 / ((b.qe - b.qb) + (a.qe - a.qb)) as f64 * (b.score + a.score) as f64 + 0.499) as i32;
    let r_s = ((b.re - a.rb) as f64 / ((b.re - b.rb) + (a.re - a.rb)) as f64 * (b.score + a.score) as f64 + 0.499) as i32;
    if (score as f64) / (q_s.max(r_s) as f64) < PATCH_MIN_SC_RATIO as f64 {
        return (0, 0);
    }
    (score, w)
}

/// mem_sort_dedup_patch(opt, 0, 0, 0, ...) as mate rescue calls it: no patching (mem_patch_reg returns 0 without bns)
pub fn sort_dedup_patch_nopatch(o: &ExtOpt, bns: &Bns, a: Vec<Reg>) -> Vec<Reg> {
    sort_dedup_patch_impl(o, bns, None, a)
}

/// mem_sort_dedup_patch
pub fn sort_dedup_patch(o: &ExtOpt, bns: &Bns, query: &[u8], a: Vec<Reg>) -> Vec<Reg> {
    sort_dedup_patch_impl(o, bns, Some(query), a)
}

fn sort_dedup_patch_impl(o: &ExtOpt, bns: &Bns, query: Option<&[u8]>, mut a: Vec<Reg>) -> Vec<Reg> {
    let n = a.len();
    if n <= 1 {
        return a;
    }
    a_idx_sort(&mut a, |r| r.re, |x, y| x < y); // ks_introsort(mem_ars2): by end position
    for r in a.iter_mut() {
        r.n_comp = 1;
    }
    for i in 1..n {
        if a[i].rid != a[i - 1].rid || a[i].rb >= a[i - 1].re + o.max_chain_gap {
            continue;
        }
        let mut j = i as isize - 1;
        while j >= 0 && a[i].rid == a[j as usize].rid && a[i].rb < a[j as usize].re + o.max_chain_gap {
            let ju = j as usize;
            if a[ju].qe == a[ju].qb {
                j -= 1;
                continue;
            }
            let (p, q) = (&a[i], &a[ju]);
            let or_ = q.re - p.rb;
            let oq = if q.qb < p.qb { (q.qe - p.qb) as i64 } else { (p.qe - q.qb) as i64 };
            let mr = (q.re - q.rb).min(p.re - p.rb);
            let mq = ((q.qe - q.qb) as i64).min((p.qe - p.qb) as i64);
            if (or_ as f32) > o.mask_level_redun * mr as f32 && (oq as f32) > o.mask_level_redun * mq as f32 {
                if p.score < q.score {
                    a[i].qe = a[i].qb;
                    break;
                } else {
                    a[ju].qe = a[ju].qb;
                }
            } else if q.rb < p.rb && query.is_some() {
                let (score, w) = patch_reg(o, bns, query.unwrap(), q, p);
                if score > 0 {
                    let q = a[ju].clone();
                    let p = &mut a[i];
                    p.n_comp += q.n_comp + 1;
                    p.seedcov = p.seedcov.max(q.seedcov);
                    p.sub = p.sub.max(q.sub);
                    p.csub = p.csub.max(q.csub);
                    p.qb = q.qb;
                    p.rb = q.rb;
                    p.score = score;
                    p.truesc = score;
                    p.w = w;
                    a[ju].qb = a[ju].qe;
                }
            }
            j -= 1;
        }
    }
    a.retain(|r| r.qe > r.qb);
    let n = a.len();
    // ks_introsort(mem_ars): score descending, then rb, then qb
    a_idx_sort(&mut a, |r| (r.score, r.rb, r.qb), |x, y| x.0 > y.0 || (x.0 == y.0 && (x.1 < y.1 || (x.1 == y.1 && x.2 < y.2))));
    for i in 1..n {
        if a[i].score == a[i - 1].score && a[i].rb == a[i - 1].rb && a[i].qb == a[i - 1].qb {
            a[i].qe = a[i].qb;
        }
    }
    let mut m = 1;
    for i in 1..n {
        if a[i].qe > a[i].qb {
            if m != i {
                a.swap(m, i);
            }
            m += 1;
        }
    }
    a.truncate(m.min(n));
    a
}

/// ks_introsort over regions (moving whole regions as bwa does). Only how it orders equal keys is particular to
/// it, so a faster sort is used first: when no two regions with equal keys differ, the order is the same whatever
/// sorted it. Otherwise ks_introsort runs on (key, index) pairs (the same comparisons, so the same order).
pub fn a_idx_sort<K: Copy>(a: &mut [Reg], key: impl Fn(&Reg) -> K, lt: impl Fn(&K, &K) -> bool) {
    if a.len() < 2 {
        return;
    }
    let mut ki: Vec<(K, u32)> = a.iter().enumerate().map(|(i, r)| (key(r), i as u32)).collect();
    let ord = |x: &K, y: &K| if lt(x, y) { std::cmp::Ordering::Less } else if lt(y, x) { std::cmp::Ordering::Greater } else { std::cmp::Ordering::Equal };
    let mut fast = ki.clone();
    fast.sort_unstable_by(|x, y| ord(&x.0, &y.0));
    // regions of equal keys must all be the same region
    let ok = fast.windows(2).all(|w| lt(&w[0].0, &w[1].0) || a[w[0].1 as usize] == a[w[1].1 as usize]);
    if ok {
        ki = fast;
    } else {
        crate::prof::count(crate::prof::N_SORTSLOW, 1);
        crate::ksort::introsort(&mut ki, |x, y| lt(&x.0, &y.0));
    }
    crate::prof::count(crate::prof::N_SORTS, 1);
    let old = a.to_vec();
    for (d, &(_, i)) in a.iter_mut().zip(&ki) {
        *d = old[i as usize];
    }
}
