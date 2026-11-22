//! Pairing and SAM output (bwamem_pair.cpp `mem_pestat`, `mem_matesw`, `mem_pair`, `mem_sam_pe`; bwamem.cpp
//! `mem_mark_primary_se`, `mem_approx_mapq_se`, `mem_reg2aln`, `mem_reg2sam`, `mem_aln2sam`; bwamem_extra.cpp
//! `mem_gen_alt`; bwa.cpp `bwa_gen_cigar2`), with bwa-mem2 2.2.1's settings: MATE_SORT 0, V17 1, `-Y`.
use std::fmt::Write as _;
use std::io::Write as _;

/// SAM text is bytes: read names and qualities are copied verbatim, as bwa-mem2's kputsn does.
trait PushStr {
    fn push_str(&mut self, s: &str);
}
impl PushStr for Vec<u8> {
    #[inline]
    fn push_str(&mut self, s: &str) {
        self.extend_from_slice(s.as_bytes());
    }
}

use crate::ext::{self, ExtOpt, Reg};
use crate::index::Bns;
use crate::ksw;

extern "C" {
    fn erfc(x: f64) -> f64;
    fn log(x: f64) -> f64;
}
fn c_log(x: f64) -> f64 {
    unsafe { log(x) }
}

/// ALIGN_RS_DUMP_PAIR=path: regions after mate rescue and primary marking (as bwa-mem2.dump's BWA_DUMP_PAIR)
pub static PDUMP: std::sync::OnceLock<std::sync::Mutex<std::fs::File>> = std::sync::OnceLock::new();

/// ALIGN_RS_PROFILE: thread time (ns) in mate SW, CIGAR (bwa_gen_cigar2), the rest of mem_reg2aln/XA, SAM text
pub static KSW_CALLS: [std::sync::atomic::AtomicU64; 3] = [const { std::sync::atomic::AtomicU64::new(0) }; 3]; // calls, u8, sum of qlen*tlen
pub static SPROF: [std::sync::atomic::AtomicU64; 4] = [const { std::sync::atomic::AtomicU64::new(0) }; 4];
fn sadd(i: usize, t: std::time::Instant) {
    SPROF[i].fetch_add(t.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
}

pub struct SamOpt {
    pub e: ExtOpt,
    pub b: i32,
    pub t: i32,
    pub min_seed_len: i32,
    pub mask_level: f32,
    pub drop_ratio: f32,
    pub xa_drop_ratio: f32,
    pub max_xa_hits: i32,
    pub max_xa_hits_alt: i32,
    pub pen_unpaired: i32,
    pub max_matesw: i32,
    pub max_ins: i32,
    pub mapq_coef_len: i32,
    pub mapq_coef_fac: f64,
    pub softclip: bool,
    pub rg_id: Option<String>,
}

impl SamOpt {
    pub fn bwa_default(softclip: bool, rg_id: Option<String>) -> SamOpt {
        SamOpt {
            e: ExtOpt::bwa_default(),
            b: 4,
            t: 30,
            min_seed_len: 19,
            mask_level: 0.5,
            drop_ratio: 0.5,
            xa_drop_ratio: 0.8,
            max_xa_hits: 5,
            max_xa_hits_alt: 200,
            pen_unpaired: 17,
            max_matesw: 50,
            max_ins: 10000,
            mapq_coef_len: 50,
            // `int mapQ_coef_fac` in mem_opt_t: log(50) = 3.91 is stored as 3
            mapq_coef_fac: c_log(50.0) as i32 as f64,
            softclip,
            rg_id,
        }
    }
}

#[derive(Clone, Copy, Default, Debug)]
pub struct Pestat {
    pub low: i32,
    pub high: i32,
    pub failed: bool,
    pub avg: f64,
    pub std: f64,
}

/// utils.h hash_64
pub fn hash_64(mut key: u64) -> u64 {
    key = key.wrapping_add(!(key << 32));
    key ^= key >> 22;
    key = key.wrapping_add(!(key << 13));
    key ^= key >> 8;
    key = key.wrapping_add(key << 3);
    key ^= key >> 15;
    key = key.wrapping_add(!(key << 27));
    key ^= key >> 31;
    key
}

fn infer_dir(l_pac: i64, b1: i64, b2: i64) -> (usize, i64) {
    let r1 = b1 >= l_pac;
    let r2 = b2 >= l_pac;
    let p2 = if r1 == r2 { b2 } else { (l_pac << 1) - 1 - b2 };
    let dist = if p2 > b1 { p2 - b1 } else { b1 - p2 };
    ((if r1 == r2 { 0 } else { 1 }) ^ (if p2 > b1 { 0 } else { 3 }), dist)
}

fn cal_sub(o: &SamOpt, r: &[Reg]) -> i32 {
    let mut j = 1;
    while j < r.len() {
        let b_max = r[j].qb.max(r[0].qb);
        let e_min = r[j].qe.min(r[0].qe);
        if e_min > b_max {
            let min_l = (r[j].qe - r[j].qb).min(r[0].qe - r[0].qb);
            if (e_min - b_max) as f32 >= min_l as f32 * o.mask_level {
                break;
            }
        }
        j += 1;
    }
    if j < r.len() {
        r[j].score
    } else {
        o.min_seed_len * o.e.a
    }
}

/// mem_pestat over a batch's regions (pairs: regs[2i], regs[2i+1])
pub fn pestat(o: &SamOpt, l_pac: i64, regs: &[Vec<Reg>]) -> [Pestat; 4] {
    let mut pes = [Pestat::default(); 4];
    let mut isize: [Vec<u64>; 4] = Default::default();
    for i in 0..regs.len() / 2 {
        let (r0, r1) = (&regs[2 * i], &regs[2 * i + 1]);
        if r0.is_empty() || r1.is_empty() {
            continue;
        }
        if cal_sub(o, r0) as f64 > 0.8 * r0[0].score as f64 {
            continue;
        }
        if cal_sub(o, r1) as f64 > 0.8 * r1[0].score as f64 {
            continue;
        }
        if r0[0].rid != r1[0].rid {
            continue;
        }
        let (dir, is) = infer_dir(l_pac, r0[0].rb, r1[0].rb);
        if is != 0 && is <= o.max_ins as i64 {
            isize[dir].push(is as u64);
        }
    }
    for d in 0..4 {
        let r = &mut pes[d];
        let q = &mut isize[d];
        if q.len() < 10 {
            r.failed = true;
            continue;
        }
        q.sort_unstable(); // ks_introsort_64 on plain values: any sort gives the same array
        let n = q.len() as f64;
        let p25 = q[(0.25 * n + 0.499) as usize] as i32;
        let p50 = q[(0.50 * n + 0.499) as usize] as i32;
        let p75 = q[(0.75 * n + 0.499) as usize] as i32;
        let _ = p50;
        r.low = ((p25 as f64) - 2.0 * (p75 - p25) as f64 + 0.499) as i32;
        if r.low < 1 {
            r.low = 1;
        }
        r.high = ((p75 as f64) + 2.0 * (p75 - p25) as f64 + 0.499) as i32;
        let (mut x, mut avg) = (0i32, 0f64);
        for &v in q.iter() {
            if v as i64 >= r.low as i64 && v as i64 <= r.high as i64 {
                avg += v as f64;
                x += 1;
            }
        }
        avg /= x as f64;
        r.avg = avg;
        let mut std = 0f64;
        for &v in q.iter() {
            if v as i64 >= r.low as i64 && v as i64 <= r.high as i64 {
                std += (v as f64 - avg) * (v as f64 - avg);
            }
        }
        r.std = (std / x as f64).sqrt();
        r.low = ((p25 as f64) - 3.0 * (p75 - p25) as f64 + 0.499) as i32;
        r.high = ((p75 as f64) + 3.0 * (p75 - p25) as f64 + 0.499) as i32;
        if (r.low as f64) > r.avg - 4.0 * r.std {
            r.low = (r.avg - 4.0 * r.std + 0.499) as i32;
        }
        if (r.high as f64) < r.avg + 4.0 * r.std {
            r.high = (r.avg + 4.0 * r.std + 0.499) as i32;
        }
        if r.low < 1 {
            r.low = 1;
        }
    }
    let max = isize.iter().map(|v| v.len()).max().unwrap_or(0);
    for d in 0..4 {
        if !pes[d].failed && (isize[d].len() as f64) < max as f64 * 0.05 {
            pes[d].failed = true;
        }
    }
    pes
}

// ------------------------------------------------------------------------------------------------
// primary marking, mapq
// ------------------------------------------------------------------------------------------------

fn mark_primary_core(o: &SamOpt, a: &mut [Reg], z: &mut Vec<usize>) {
    let mut tmp = o.e.a + o.b;
    tmp = tmp.max(o.e.o_del + o.e.e_del);
    tmp = tmp.max(o.e.o_ins + o.e.e_ins);
    z.clear();
    z.push(0);
    for i in 1..a.len() {
        let mut k = 0;
        while k < z.len() {
            let j = z[k];
            let b_max = a[j].qb.max(a[i].qb);
            let e_min = a[j].qe.min(a[i].qe);
            if e_min > b_max {
                let min_l = (a[i].qe - a[i].qb).min(a[j].qe - a[j].qb);
                if (e_min - b_max) as f32 >= min_l as f32 * o.mask_level {
                    if a[j].sub == 0 {
                        a[j].sub = a[i].score;
                    }
                    if a[j].score - a[i].score <= tmp && (a[j].is_alt != 0 || a[i].is_alt == 0) {
                        a[j].sub_n += 1;
                    }
                    break;
                }
            }
            k += 1;
        }
        if k == z.len() {
            z.push(i);
        } else {
            a[i].secondary = z[k] as i32;
        }
    }
}

fn sort_regs<K: Copy>(a: &mut [Reg], key: impl Fn(&Reg) -> K, lt: impl Fn(&K, &K) -> bool) {
    ext::a_idx_sort(a, key, lt)
}

/// mem_mark_primary_se: returns n_pri
pub fn mark_primary_se(o: &SamOpt, a: &mut Vec<Reg>, id: i64) -> usize {
    let n = a.len();
    if n == 0 {
        return 0;
    }
    let mut n_pri = 0;
    for (i, r) in a.iter_mut().enumerate() {
        r.sub = 0;
        r.alt_sc = 0;
        r.secondary = -1;
        r.secondary_all = -1;
        r.hash = hash_64((id + i as i64) as u64);
        if r.is_alt == 0 {
            n_pri += 1;
        }
    }
    sort_regs(a, |r| (r.score, r.is_alt, r.hash), |x, y| x.0 > y.0 || (x.0 == y.0 && (x.1 < y.1 || (x.1 == y.1 && x.2 < y.2))));
    let mut z = Vec::new();
    mark_primary_core(o, a, &mut z);
    for i in 0..n {
        a[i].secondary_all = i as i32;
        let s = a[i].secondary;
        if a[i].is_alt == 0 && s >= 0 && a[s as usize].is_alt != 0 {
            a[i].alt_sc = a[s as usize].score;
        }
    }
    if n_pri < n {
        let mut zz = vec![0usize; n];
        if n_pri > 0 {
            sort_regs(a, |r| (r.is_alt, r.score, r.hash), |x, y| x.0 < y.0 || (x.0 == y.0 && (x.1 > y.1 || (x.1 == y.1 && x.2 < y.2))));
        }
        for i in 0..n {
            zz[a[i].secondary_all as usize] = i;
        }
        for i in 0..n {
            if a[i].secondary >= 0 {
                a[i].secondary_all = zz[a[i].secondary as usize] as i32;
                if a[i].is_alt != 0 {
                    a[i].secondary = i32::MAX;
                }
            } else {
                a[i].secondary_all = -1;
            }
        }
        if n_pri > 0 {
            for r in a.iter_mut().take(n_pri) {
                r.sub = 0;
                r.secondary = -1;
            }
            mark_primary_core(o, &mut a[..n_pri], &mut z);
        }
    } else {
        for r in a.iter_mut() {
            r.secondary_all = r.secondary;
        }
    }
    n_pri
}

/// mem_approx_mapq_se
pub fn approx_mapq_se(o: &SamOpt, a: &Reg) -> i32 {

    let mut sub = if a.sub != 0 { a.sub } else { o.min_seed_len * o.e.a };
    sub = sub.max(a.csub);
    if sub >= a.score {
        return 0;
    }
    let l = ((a.qe - a.qb) as i64).max(a.re - a.rb) as i32;
    let identity = 1.0 - ((l * o.e.a - a.score) as f64) / (o.e.a + o.b) as f64 / l as f64;
    let mut mapq;
    if a.score == 0 {
        mapq = 0;
    } else {
        let mut tmp = if l < o.mapq_coef_len { 1.0 } else { o.mapq_coef_fac / c_log(l as f64) };
        tmp *= identity * identity;
        mapq = (6.02 * (a.score - sub) as f64 / o.e.a as f64 * tmp * tmp + 0.499) as i32;
    }
    if a.sub_n > 0 {
        mapq -= (4.343 * c_log((a.sub_n + 1) as f64) + 0.499) as i32;
    }
    if mapq > 60 {
        mapq = 60;
    }
    if mapq < 0 {
        mapq = 0;
    }
    ((mapq as f64) * (1.0 - a.frac_rep as f64) + 0.499) as i32
}

// ------------------------------------------------------------------------------------------------
// CIGAR (bwa_gen_cigar2 with ksw_global2's backtrack), mem_reg2aln
// ------------------------------------------------------------------------------------------------

const MINUS_INF: i32 = -0x4000_0000;

fn push_cigar(c: &mut Vec<u32>, op: u32, len: u32) {
    if c.is_empty() || op != (c[c.len() - 1] & 0xf) {
        c.push(len << 4 | op);
    } else {
        let n = c.len();
        c[n - 1] += len << 4;
    }
}

/// ksw_global2 with backtrack
fn ksw_global2(o: &ExtOpt, query: &[u8], target: &[u8], w: i32) -> (i32, Vec<u32>) {
    // per-thread scratch (the direction matrix, query profile and the H/E row), reused across calls
    thread_local! {
        static SCRATCH: std::cell::RefCell<(Vec<u8>, Vec<i8>, Vec<(i32, i32)>)> = const { std::cell::RefCell::new((Vec::new(), Vec::new(), Vec::new())) };
    }
    SCRATCH.with(|sc| {
        let mut sc = sc.borrow_mut();
        let (z, qp, eh) = &mut *sc;
        ksw_global2_in(o, query, target, w, z, qp, eh)
    })
}

fn ksw_global2_in(o: &ExtOpt, query: &[u8], target: &[u8], w: i32, z: &mut Vec<u8>, qp: &mut Vec<i8>, eh: &mut Vec<(i32, i32)>) -> (i32, Vec<u32>) {
    let qlen = query.len() as i32;
    let tlen = target.len() as i32;
    let m = 5usize;
    let (oe_del, oe_ins) = (o.o_del + o.e_del, o.o_ins + o.e_ins);
    let (e_del, e_ins) = (o.e_del, o.e_ins);
    let n_col = qlen.min(2 * w + 1) as usize;
    z.clear();
    z.resize(n_col * target.len(), 0);
    qp.clear();
    for k in 0..m {
        let p = &o.mat[k * m..k * m + m];
        qp.extend(query.iter().map(|&qj| p[qj as usize]));
    }
    // (H, E) per query position
    eh.clear();
    eh.resize(query.len() + 1, (MINUS_INF, MINUS_INF));
    eh[0].0 = 0;
    let mut j = 1;
    while j <= qlen && j <= w {
        eh[j as usize].0 = -(o.o_ins + o.e_ins * j);
        j += 1;
    }
    for i in 0..tlen {
        let mut f = MINUS_INF;
        let q = &qp[target[i as usize] as usize * query.len()..];
        let beg = if i > w { i - w } else { 0 };
        let end = if i + w + 1 < qlen { i + w + 1 } else { qlen };
        let mut h1 = if beg == 0 { -(o.o_del + o.e_del * (i + 1)) } else { MINUS_INF };
        let (b, e_) = (beg as usize, end as usize);
        let zi = &mut z[i as usize * n_col..i as usize * n_col + (e_ - b)];
        for ((cell, &qs), zc) in eh[b..e_].iter_mut().zip(&q[b..e_]).zip(zi.iter_mut()) {
            let mut mv = cell.0;
            let mut e = cell.1;
            cell.0 = h1;
            mv += qs as i32;
            let mut d: u8 = if mv >= e { 0 } else { 1 };
            let mut h = if mv >= e { mv } else { e };
            d = if h >= f { d } else { 2 };
            h = if h >= f { h } else { f };
            h1 = h;
            let t = mv - oe_del;
            e -= e_del;
            d |= if e > t { 1 << 2 } else { 0 };
            e = if e > t { e } else { t };
            cell.1 = e;
            let t = mv - oe_ins;
            f -= e_ins;
            d |= if f > t { 2 << 4 } else { 0 };
            f = if f > t { f } else { t };
            *zc = d;
        }
        eh[e_] = (h1, MINUS_INF);
    }
    let score = eh[qlen as usize].0;
    let mut cigar: Vec<u32> = Vec::new();
    let mut i = tlen - 1;
    let mut k = (if i + w + 1 < qlen { i + w + 1 } else { qlen }) - 1;
    let mut which = 0u32;
    while i >= 0 && k >= 0 {
        let beg = if i > w { i - w } else { 0 };
        which = (z[i as usize * n_col + (k - beg) as usize] as u32 >> (which << 1)) & 3;
        if which == 0 {
            push_cigar(&mut cigar, 0, 1);
            i -= 1;
            k -= 1;
        } else if which == 1 {
            push_cigar(&mut cigar, 2, 1);
            i -= 1;
        } else {
            push_cigar(&mut cigar, 1, 1);
            k -= 1;
        }
    }
    if i >= 0 {
        push_cigar(&mut cigar, 2, (i + 1) as u32);
    }
    if k >= 0 {
        push_cigar(&mut cigar, 1, (k + 1) as u32);
    }
    cigar.reverse();
    (score, cigar)
}

/// bwa_gen_cigar2 with CIGAR, NM and MD: (score, cigar, NM, MD); None if rejected
fn gen_cigar2(o: &ExtOpt, bns: &Bns, w_: i32, query: &[u8], rb: i64, re: i64) -> Option<(i32, Vec<u32>, i32, String)> {
    let l_pac = bns.l_pac;
    let l_query = query.len() as i32;
    if l_query <= 0 || rb >= re || (rb < l_pac && re > l_pac) {
        return None;
    }
    let mut rseq = ext::bns_get_seq(bns, rb, re);
    if re - rb != rseq.len() as i64 {
        return None; // bwa returns cigar 0 (n_cigar 0) here
    }
    let mut q = query.to_vec();
    if rb >= l_pac {
        q.reverse();
        rseq.reverse();
    }
    let rlen = rseq.len() as i32;
    crate::prof::count(crate::prof::N_CIGAR, 1);
    let (score, cigar) = if l_query == (re - rb) as i32 && w_ == 0 {
        (((0..q.len()).map(|i| o.mat[rseq[i] as usize * 5 + q[i] as usize] as i32).sum()), vec![(l_query as u32) << 4])
    } else {
        let max_ins = ((((l_query + 1) >> 1) * o.mat[0] as i32 - o.o_ins) as f64 / o.e_ins as f64 + 1.0) as i32;
        let max_del = ((((l_query + 1) >> 1) * o.mat[0] as i32 - o.o_del) as f64 / o.e_del as f64 + 1.0) as i32;
        let max_gap = max_ins.max(max_del).max(1);
        let mut w = (max_gap + (rlen - l_query).abs() + 1) >> 1;
        w = w.min(w_);
        let min_w = (rlen - l_query).abs() + 3;
        w = w.max(min_w);
        crate::prof::count(crate::prof::N_CIGDP, 1);
        let mode = ext::cert_mode();
        let kc = if mode != 0 { ext::cert_band(o, &q, &rseq) } else { -1 };
        if rlen != l_query {
            crate::prof::count(crate::prof::N_CIGUNEQ, 1);
        }
        let ungapped = || (q.iter().zip(&rseq).map(|(&a, &b)| o.mat[b as usize * 5 + a as usize] as i32).sum::<i32>(), vec![(l_query as u32) << 4]);
        let res = if kc == 0 {
            crate::prof::count(crate::prof::N_CERT0, 1);
            ungapped()
        } else {
            let wc = if kc > 0 && kc < w as i64 { crate::prof::count(crate::prof::N_CERTW, 1); kc as i32 } else { w };
            crate::prof::count(crate::prof::N_CIGCELLS, q.len() as u64 * (2 * wc as u64 + 1).min(rlen as u64));
            ksw_global2(o, &q, &rseq, wc)
        };
        if mode == 2 && kc >= 0 {
            let full = ksw_global2(o, &q, &rseq, w);
            assert!(res == full, "certified band {kc} (band {w}): {res:?} differs from {full:?}");
        }
        res
    };
    // NM and MD
    let int2base: &[u8; 5] = if rb < l_pac { b"ACGTN" } else { b"TGCAN" };
    let mut md = String::new();
    let (mut x, mut y, mut u, mut n_mm, mut n_gap) = (0usize, 0usize, 0i32, 0i32, 0i32);
    let nc = cigar.len();
    for (k, &c) in cigar.iter().enumerate() {
        let (op, len) = (c & 0xf, (c >> 4) as usize);
        if op == 0 {
            for i in 0..len {
                if q[x + i] != rseq[y + i] {
                    let _ = write!(md, "{u}");
                    md.push(int2base[rseq[y + i] as usize] as char);
                    n_mm += 1;
                    u = 0;
                } else {
                    u += 1;
                }
            }
            x += len;
            y += len;
        } else if op == 2 {
            if k > 0 && k < nc - 1 {
                let _ = write!(md, "{u}");
                md.push('^');
                for i in 0..len {
                    md.push(int2base[rseq[y + i] as usize] as char);
                }
                u = 0;
                n_gap += len as i32;
            }
            y += len;
        } else if op == 1 {
            x += len;
            n_gap += len as i32;
        }
    }
    let _ = write!(md, "{u}");
    Some((score, cigar, n_mm + n_gap, md))
}

/// mem_aln_t
#[derive(Clone, Debug, Default)]
pub struct Aln {
    pub pos: i64,
    pub rid: i32,
    pub flag: i32,
    pub is_rev: bool,
    pub is_alt: bool,
    pub mapq: i32,
    pub nm: i32,
    pub cigar: Vec<u32>,
    pub md: String,
    pub score: i32,
    pub sub: i32,
    pub alt_sc: i32,
    pub xa: Option<String>,
}

fn infer_bw(l1: i32, l2: i32, score: i32, a: i32, q: i32, r: i32) -> i32 {
    if l1 == l2 && l1 * a - score < (q + r - a) << 1 {
        return 0;
    }
    let mut w = ((((l1.min(l2)) * a - score - q) as f64) / r as f64 + 2.0) as i32;
    if w < (l1 - l2).abs() {
        w = (l1 - l2).abs();
    }
    w
}

/// mem_reg2aln (None region: unmapped record)
pub fn reg2aln(o: &SamOpt, bns: &Bns, query: &[u8], ar: Option<&Reg>) -> Aln {
    let mut a = Aln::default();
    let ar = match ar {
        Some(r) if r.rb >= 0 && r.re >= 0 => r,
        _ => {
            a.rid = -1;
            a.pos = -1;
            a.flag |= 0x4;
            return a;
        }
    };
    let l_query = query.len() as i32;
    let (qb, qe, rb, re) = (ar.qb, ar.qe, ar.rb, ar.re);
    a.mapq = if ar.secondary < 0 { approx_mapq_se(o, ar) } else { 0 };
    if ar.secondary >= 0 {
        a.flag |= 0x100;
    }
    let tmp = infer_bw(qe - qb, (re - rb) as i32, ar.truesc, o.e.a, o.e.o_del, o.e.e_del);
    let mut w2 = infer_bw(qe - qb, (re - rb) as i32, ar.truesc, o.e.a, o.e.o_ins, o.e.e_ins);
    w2 = w2.max(tmp);
    if w2 > o.e.w {
        w2 = w2.min(ar.w);
    }
    let mut last_sc = -(1 << 30);
    let mut i = 0;
    let mut res;
    loop {
        w2 = w2.min(o.e.w << 2);
        let tg = std::time::Instant::now();
        res = crate::prof::cpu(crate::prof::C_CIGAR, || gen_cigar2(&o.e, bns, w2, &query[qb as usize..qe as usize], rb, re));
        sadd(1, tg);
        let score = res.as_ref().map_or(0, |r| r.0);
        if score == last_sc || w2 == o.e.w << 2 {
            break;
        }
        last_sc = score;
        w2 <<= 1;
        i += 1;
        if !(i < 3 && score < ar.truesc - o.e.a) {
            break;
        }
    }
    let (_, mut cigar, nm, md) = res.expect("bwa_gen_cigar2 returned no CIGAR");
    a.nm = nm;
    let (mut pos, is_rev) = bns.depos(if rb < bns.l_pac { rb } else { re - 1 });
    a.is_rev = is_rev;
    if !cigar.is_empty() {
        if cigar[0] & 0xf == 2 {
            pos += (cigar[0] >> 4) as i64;
            cigar.remove(0);
        } else if cigar[cigar.len() - 1] & 0xf == 2 {
            cigar.pop();
        }
    }
    if qb != 0 || qe != l_query {
        let clip5 = if is_rev { l_query - qe } else { qb };
        let clip3 = if is_rev { qb } else { l_query - qe };
        if clip5 != 0 {
            cigar.insert(0, (clip5 as u32) << 4 | 3);
        }
        if clip3 != 0 {
            cigar.push((clip3 as u32) << 4 | 3);
        }
    }
    a.rid = bns.pos2rid(pos);
    assert_eq!(a.rid, ar.rid);
    a.pos = pos - bns.anns[a.rid as usize].offset;
    a.score = ar.score;
    a.sub = ar.sub.max(ar.csub);
    a.is_alt = ar.is_alt != 0;
    a.alt_sc = ar.alt_sc;
    a.cigar = cigar;
    a.md = md;
    a
}

fn get_pri_idx(xa_drop_ratio: f32, a: &[Reg], i: usize) -> i32 {
    let k = a[i].secondary_all;
    if k >= 0 && a[i].score as f64 >= a[k as usize].score as f64 * xa_drop_ratio as f64 {
        return k;
    }
    -1
}

/// mem_gen_alt: XA strings per region index
pub fn gen_alt(o: &SamOpt, bns: &Bns, a: &[Reg], query: &[u8]) -> Option<Vec<Option<String>>> {
    let n = a.len();
    let mut cnt = vec![0i32; n];
    let mut has_alt = vec![false; n];
    let mut tot = 0;
    for i in 0..n {
        let r = get_pri_idx(o.xa_drop_ratio, a, i);
        if r >= 0 {
            cnt[r as usize] += 1;
            tot += 1;
            if a[i].is_alt != 0 {
                has_alt[r as usize] = true;
            }
        }
    }
    if tot == 0 {
        return None;
    }
    let mut aln: Vec<Option<String>> = vec![None; n];
    for i in 0..n {
        let r = get_pri_idx(o.xa_drop_ratio, a, i);
        if r < 0 {
            continue;
        }
        let r = r as usize;
        if cnt[r] > o.max_xa_hits_alt || (!has_alt[r] && cnt[r] > o.max_xa_hits) {
            continue;
        }
        let t = reg2aln(o, bns, query, Some(&a[i]));
        let mut s = String::new();
        s.push_str(&bns.anns[t.rid as usize].name);
        s.push(',');
        s.push(if t.is_rev { '-' } else { '+' });
        let _ = write!(s, "{}", (t.pos + 1));
        s.push(',');
        for &c in &t.cigar {
            let _ = write!(s, "{}", (c >> 4));
            s.push(b"MIDSHN"[(c & 0xf) as usize] as char);
        }
        s.push(',');
        let _ = write!(s, "{}", t.nm);
        s.push(';');
        crate::prof::count(crate::prof::N_XA, 1);
        aln[r].get_or_insert_with(String::new).push_str(&s);
    }
    Some(aln)
}

fn get_rlen(cigar: &[u32]) -> i64 {
    cigar.iter().filter(|&&c| c & 0xf == 0 || c & 0xf == 2).map(|&c| (c >> 4) as i64).sum()
}

fn add_cigar(o: &SamOpt, p: &Aln, out: &mut Vec<u8>, which: usize) {
    if !p.cigar.is_empty() {
        for &c in &p.cigar {
            let mut op = c & 0xf;
            if !o.softclip && !p.is_alt && (op == 3 || op == 4) {
                op = if which != 0 { 4 } else { 3 };
            }
            let _ = write!(out, "{}", (c >> 4));
            out.push(b"MIDSH"[op as usize]);
        }
    } else {
        out.push(b'*');
    }
}

pub struct Read<'a> {
    pub name: &'a [u8],
    pub seq: &'a [u8], // 0..4
    pub qual: &'a [u8],
}

/// mem_aln2sam
pub fn aln2sam(o: &SamOpt, bns: &Bns, out: &mut Vec<u8>, s: &Read, list: &[Aln], which: usize, m_: Option<&Aln>) {
    let t = std::time::Instant::now();
    crate::prof::cpu(crate::prof::C_SAMTEXT, || aln2sam_(o, bns, out, s, list, which, m_));
    sadd(3, t);
}

fn aln2sam_(o: &SamOpt, bns: &Bns, out: &mut Vec<u8>, s: &Read, list: &[Aln], which: usize, m_: Option<&Aln>) {
    let mut p = list[which].clone();
    let mut m = m_.cloned();
    if m.is_some() {
        p.flag |= 0x1;
    }
    if p.rid < 0 {
        p.flag |= 0x4;
    }
    if let Some(mm) = &m {
        if mm.rid < 0 {
            p.flag |= 0x8;
        }
    }
    if let Some(mm) = &m {
        if p.rid < 0 && mm.rid >= 0 {
            p.rid = mm.rid;
            p.pos = mm.pos;
            p.is_rev = mm.is_rev;
            p.cigar.clear();
        }
    }
    if let Some(mm) = m.as_mut() {
        if mm.rid < 0 && p.rid >= 0 {
            mm.rid = p.rid;
            mm.pos = p.pos;
            mm.is_rev = p.is_rev;
            mm.cigar.clear();
        }
    }
    if p.is_rev {
        p.flag |= 0x10;
    }
    if let Some(mm) = &m {
        if mm.is_rev {
            p.flag |= 0x20;
        }
    }
    out.extend_from_slice(s.name);
    out.push(b'\t');
    let _ = write!(out, "{}", ((p.flag & 0xffff) | if p.flag & 0x10000 != 0 { 0x100 } else { 0 }));
    out.push(b'\t');
    if p.rid >= 0 {
        out.push_str(&bns.anns[p.rid as usize].name);
        out.push(b'\t');
        let _ = write!(out, "{}", (p.pos + 1));
        out.push(b'\t');
        let _ = write!(out, "{}", p.mapq);
        out.push(b'\t');
        add_cigar(o, &p, out, which);
    } else {
        out.push_str("*\t0\t0\t*");
    }
    out.push(b'\t');
    match &m {
        Some(mm) if mm.rid >= 0 => {
            if p.rid == mm.rid {
                out.push(b'=');
            } else {
                out.push_str(&bns.anns[mm.rid as usize].name);
            }
            out.push(b'\t');
            let _ = write!(out, "{}", (mm.pos + 1));
            out.push(b'\t');
            if p.rid == mm.rid {
                let p0 = p.pos + if p.is_rev { get_rlen(&p.cigar) - 1 } else { 0 };
                let p1 = mm.pos + if mm.is_rev { get_rlen(&mm.cigar) - 1 } else { 0 };
                if mm.cigar.is_empty() || p.cigar.is_empty() {
                    out.push(b'0');
                } else {
                    let v = -(p0 - p1 + if p0 > p1 { 1 } else if p0 < p1 { -1 } else { 0 });
                    let _ = write!(out, "{}", v);
                }
            } else {
                out.push(b'0');
            }
        }
        _ => out.push_str("*\t0\t0"),
    }
    out.push(b'\t');
    let l_seq = s.seq.len();
    if p.flag & 0x100 != 0 {
        out.push_str("*\t*");
    } else if !p.is_rev {
        let (mut qb, mut qe) = (0usize, l_seq);
        if !p.cigar.is_empty() && which != 0 && !o.softclip && !p.is_alt {
            let (c0, cn) = (p.cigar[0], p.cigar[p.cigar.len() - 1]);
            if c0 & 0xf == 4 || c0 & 0xf == 3 {
                qb += (c0 >> 4) as usize;
            }
            if cn & 0xf == 4 || cn & 0xf == 3 {
                qe -= (cn >> 4) as usize;
            }
        }
        for i in qb..qe {
            out.push(b"ACGTN"[s.seq[i] as usize]);
        }
        out.push(b'\t');
        if !s.qual.is_empty() {
            out.extend_from_slice(&s.qual[qb..qe]);
        } else {
            out.push(b'*');
        }
    } else {
        let (mut qb, mut qe) = (0usize, l_seq);
        if !p.cigar.is_empty() && which != 0 && !o.softclip && !p.is_alt {
            let (c0, cn) = (p.cigar[0], p.cigar[p.cigar.len() - 1]);
            if c0 & 0xf == 4 || c0 & 0xf == 3 {
                qe -= (c0 >> 4) as usize;
            }
            if cn & 0xf == 4 || cn & 0xf == 3 {
                qb += (cn >> 4) as usize;
            }
        }
        for i in (qb..qe).rev() {
            out.push(b"TGCAN"[s.seq[i] as usize]);
        }
        out.push(b'\t');
        if !s.qual.is_empty() {
            for i in (qb..qe).rev() {
                out.push(s.qual[i]);
            }
        } else {
            out.push(b'*');
        }
    }
    if !p.cigar.is_empty() {
        out.push_str("\tNM:i:");
        let _ = write!(out, "{}", p.nm);
        out.push_str("\tMD:Z:");
        out.push_str(&p.md);
    }
    if let Some(mm) = &m {
        if !mm.cigar.is_empty() {
            out.push_str("\tMC:Z:");
            add_cigar(o, mm, out, which);
        }
    }
    if p.score >= 0 {
        out.push_str("\tAS:i:");
        let _ = write!(out, "{}", p.score);
    }
    if p.sub >= 0 {
        out.push_str("\tXS:i:");
        let _ = write!(out, "{}", p.sub);
    }
    if let Some(id) = &o.rg_id {
        out.push_str("\tRG:Z:");
        out.push_str(id);
    }
    if p.flag & 0x100 == 0 {
        let others = (0..list.len()).any(|i| i != which && list[i].flag & 0x100 == 0);
        if others {
            out.push_str("\tSA:Z:");
            for (i, r) in list.iter().enumerate() {
                if i == which || r.flag & 0x100 != 0 {
                    continue;
                }
                out.push_str(&bns.anns[r.rid as usize].name);
                out.push(b',');
                let _ = write!(out, "{}", (r.pos + 1));
                out.push(b',');
                out.push(if r.is_rev { b'-' } else { b'+' });
                out.push(b',');
                for &c in &r.cigar {
                    let _ = write!(out, "{}", (c >> 4));
                    out.push(b"MIDSH"[(c & 0xf) as usize]);
                }
                out.push(b',');
                let _ = write!(out, "{}", r.mapq);
                out.push(b',');
                let _ = write!(out, "{}", r.nm);
                out.push(b';');
            }
        }
        if p.alt_sc > 0 {
            let _ = write!(out, "\tpa:f:{:.3}", p.score as f64 / p.alt_sc as f64);
        }
    }
    if let Some(xa) = &p.xa {
        out.push_str("\tXA:Z:");
        out.push_str(xa);
    }
    out.push(b'\n');
}

/// mem_reg2sam
pub fn reg2sam(o: &SamOpt, bns: &Bns, out: &mut Vec<u8>, s: &Read, a: &[Reg], extra_flag: i32, m: Option<&Aln>) {
    let xa = gen_alt(o, bns, a, s.seq);
    let mut aa: Vec<Aln> = Vec::new();
    let mut l = 0;
    for (k, p) in a.iter().enumerate() {
        if p.score < o.t {
            continue;
        }
        if p.secondary >= 0 && (p.is_alt != 0 || true) {
            continue; // MEM_F_ALL not set
        }
        if p.secondary >= 0 && p.secondary < i32::MAX && (p.score as f32) < a[p.secondary as usize].score as f32 * o.drop_ratio {
            continue;
        }
        let mut q = reg2aln(o, bns, s.seq, Some(p));
        q.xa = xa.as_ref().and_then(|x| x[k].clone());
        q.flag |= extra_flag;
        if p.secondary >= 0 {
            q.sub = -1;
        }
        if l != 0 && p.secondary < 0 {
            q.flag |= 0x800;
        }
        if l != 0 && p.is_alt == 0 && q.mapq > aa[0].mapq {
            q.mapq = aa[0].mapq;
        }
        aa.push(q);
        l += 1;
    }
    if aa.is_empty() {
        let mut t = reg2aln(o, bns, s.seq, None);
        t.flag |= extra_flag;
        aln2sam(o, bns, out, s, &[t], 0, m);
    } else {
        for k in 0..aa.len() {
            aln2sam(o, bns, out, s, &aa, k, m);
        }
    }
}

// ------------------------------------------------------------------------------------------------
// pairing
// ------------------------------------------------------------------------------------------------

/// A precomputed mate-rescue alignment: (side i, mate-hit index j, orientation r) -> result
pub type PreSw = (u8, u16, u8, ksw::Kswr);

/// A mate-rescue job for the GPU: which read is the mate, reverse-complemented or not, the window
#[cfg_attr(not(gpu), allow(dead_code))] // GPU build only
#[derive(Clone, Copy, Debug)]
pub struct SwSpec {
    pub read: u32,
    pub rev: bool,
    pub rb: i64,
    pub tlen: i32,
    pub xtra: i32,
}

fn skip_dirs(bns: &Bns, pes: &[Pestat; 4], a: &Reg, ma: &[Reg]) -> [bool; 4] {
    let mut skip = [false; 4];
    for r in 0..4 {
        skip[r] = pes[r].failed;
    }
    for x in ma.iter() {
        let (r, dist) = infer_dir(bns.l_pac, a.rb, x.rb);
        if dist >= pes[r].low as i64 && dist <= pes[r].high as i64 {
            skip[r] = true;
        }
    }
    skip
}

/// The reference window of orientation r for a mate of l_ms bases near `a` (bns_fetch_seq's clamping); None if
/// mem_matesw would not align there
fn mate_window(o: &SamOpt, bns: &Bns, pes: &[Pestat; 4], a: &Reg, l_ms: i64, r: usize) -> Option<(i64, i64)> {
    let l_pac = bns.l_pac;
    let is_rev = (r >> 1) != (r & 1);
    let is_larger = (r >> 1) == 0;
    let (low, high) = (pes[r].low as i64, pes[r].high as i64);
    let (mut rb, mut re);
    if !is_rev {
        rb = if is_larger { a.rb + low } else { a.rb - high };
        re = (if is_larger { a.rb + high } else { a.rb - low }) + l_ms;
    } else {
        rb = (if is_larger { a.rb + low } else { a.rb - high }) - l_ms;
        re = if is_larger { a.rb + high } else { a.rb - low };
    }
    if rb < 0 {
        rb = 0;
    }
    if re > l_pac << 1 {
        re = l_pac << 1;
    }
    let mut rid = -1;
    if rb < re {
        let mid = (rb + re) >> 1;
        let (pf, rev) = bns.depos(mid);
        rid = bns.pos2rid(pf);
        let mut far_beg = bns.anns[rid as usize].offset;
        let mut far_end = far_beg + bns.anns[rid as usize].len as i64;
        if rev {
            let t = far_beg;
            far_beg = (l_pac << 1) - far_end;
            far_end = (l_pac << 1) - t;
        }
        rb = rb.max(far_beg);
        re = re.min(far_end);
    }
    (a.rid == rid && re - rb >= o.min_seed_len as i64).then_some((rb, re))
}

fn sw_xtra(o: &SamOpt, l_ms: i64) -> i32 {
    ksw::KSW_XSUBO | ksw::KSW_XSTART | (if l_ms * (o.e.a as i64) < 250 { ksw::KSW_XBYTE } else { 0 }) | (o.min_seed_len * o.e.a)
}

/// mem_matesw_batch_pre: every rescue alignment a pair may need, from its regions before any rescue
#[cfg_attr(not(gpu), allow(dead_code))] // GPU build only
pub fn matesw_jobs(o: &SamOpt, bns: &Bns, pes: &[Pestat; 4], a: &[Vec<Reg>], l_seq: [usize; 2], read: [u32; 2], out: &mut Vec<(u8, u16, u8, SwSpec)>) {
    for i in 0..2 {
        if a[i].is_empty() {
            continue;
        }
        let top = a[i][0].score;
        let b: Vec<&Reg> = a[i].iter().filter(|r| r.score >= top - o.pen_unpaired).collect();
        for (j, bj) in b.iter().enumerate().take(o.max_matesw as usize) {
            let skip = skip_dirs(bns, pes, bj, &a[1 - i]);
            if skip.iter().all(|&x| x) {
                continue;
            }
            for r in 0..4 {
                if skip[r] {
                    continue;
                }
                let l_ms = l_seq[1 - i] as i64;
                if let Some((rb, re)) = mate_window(o, bns, pes, bj, l_ms, r) {
                    let rev = (r >> 1) != (r & 1);
                    out.push((i as u8, j as u16, r as u8, SwSpec { read: read[1 - i], rev, rb, tlen: (re - rb) as i32, xtra: sw_xtra(o, l_ms) }));
                }
            }
        }
    }
}

/// mem_matesw: rescue the mate (`ms`, 0..4) near `a`; adds hits to `ma`; returns the number of SW runs.
/// `pre`: precomputed alignments of this pair (looked up by side i, hit index j, orientation r).
/// ALIGN_RS_DEDUP_CHECK=1: each of mate rescue's region dedups also runs the full `sort_dedup_patch_nopatch` on a
/// copy and must agree (a check of `ext::dedup_add` and of skipping a dedup that changes nothing)
fn dedup_check() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("ALIGN_RS_DEDUP_CHECK").is_ok_and(|v| v == "1"))
}

/// bwa's insertion of a rescued region into ma: before the first region of a lower score
fn insert_by_score(ma: &mut Vec<Reg>, b: Reg) {
    let mut i = 0;
    while i < ma.len() {
        if ma[i].score < b.score {
            break;
        }
        i += 1;
    }
    ma.insert(i, b);
}

/// mem_matesw's mem_sort_dedup_patch(opt, 0, 0, 0, ...) after a direction, with `added` the region it rescued (if
/// any). `clean`: ma is clean for `ext::dedup_add` (as the last full dedup or dedup_add left it): then a direction
/// that rescued nothing leaves ma as it is (sorting and deduplicating it again changes nothing), and a rescued
/// region is added without sorting; otherwise the full function runs.
fn rescue_dedup(o: &SamOpt, bns: &Bns, ma: &mut Vec<Reg>, added: Option<Reg>, clean: &mut bool) {
    let want = dedup_check().then(|| {
        let mut v = ma.clone();
        if let Some(b) = added {
            insert_by_score(&mut v, b);
        }
        ext::sort_dedup_patch_nopatch(&o.e, bns, v)
    });
    let fast = if !*clean {
        None
    } else {
        match added {
            None => Some(true),
            Some(b) => ext::dedup_add(&o.e, ma, b),
        }
    };
    crate::prof::count(
        match (&fast, added) {
            (Some(_), Some(_)) => crate::prof::N_DDADD,
            (Some(_), None) => crate::prof::N_DDSKIP,
            _ => crate::prof::N_DDFULL,
        },
        1,
    );
    match fast {
        Some(c) => *clean = c,
        None => {
            if let Some(b) = added {
                insert_by_score(ma, b);
            }
            let v = std::mem::take(ma);
            (*ma, *clean) = ext::sort_dedup_nopatch_clean(&o.e, bns, v);
        }
    }
    if let Some(w) = want {
        assert!(w == *ma, "mate rescue: incremental region dedup differs from mem_sort_dedup_patch");
    }
}

#[allow(clippy::too_many_arguments)]
fn matesw(o: &SamOpt, bns: &Bns, pes: &[Pestat; 4], a: &Reg, ms: &[u8], ma: &mut Vec<Reg>, clean: &mut bool, pre: Option<&[PreSw]>, i_side: usize, j_hit: usize) -> i32 {
    let l_pac = bns.l_pac;
    let l_ms = ms.len() as i64;
    let skip = skip_dirs(bns, pes, a, ma);
    if skip.iter().all(|&s| s) {
        return 0;
    }
    let mut n = 0;
    for r in 0..4 {
        if skip[r] {
            continue;
        }
        let is_rev = (r >> 1) != (r & 1);
        let mut added = None;
        if let Some((rb, re)) = mate_window(o, bns, pes, a, l_ms, r) {
            let xtra = sw_xtra(o, l_ms);
            let found = pre.and_then(|p| p.iter().find(|x| x.0 as usize == i_side && x.1 as usize == j_hit && x.2 as usize == r).map(|x| x.3));
            let aln = match found {
                Some(k) => k,
                None => {
                    let seq: Vec<u8> = if is_rev { ms.iter().rev().map(|&c| if c < 4 { 3 - c } else { 4 }).collect() } else { ms.to_vec() };
                    let reff = ext::bns_get_seq(bns, rb, re);
                    assert_eq!(reff.len() as i64, re - rb);
                    let tk = std::time::Instant::now();
                    KSW_CALLS[0].fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if xtra & ksw::KSW_XBYTE != 0 {
                        KSW_CALLS[1].fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    }
                    KSW_CALLS[2].fetch_add((seq.len() * reff.len()) as u64, std::sync::atomic::Ordering::Relaxed);
                    let k = ksw::ksw_align2(&seq, &reff, &o.e.mat, o.e.o_del, o.e.e_del, o.e.o_ins, o.e.e_ins, xtra);
                    sadd(0, tk);
                    k
                }
            };
            if aln.score >= o.min_seed_len && aln.qb >= 0 {
                let mut b = Reg { rid: a.rid, is_alt: a.is_alt, score: aln.score, csub: aln.score2, secondary: -1, ..Default::default() };
                let l_ms = l_ms as i32;
                b.qb = if is_rev { l_ms - (aln.qe + 1) } else { aln.qb };
                b.qe = if is_rev { l_ms - aln.qb } else { aln.qe + 1 };
                b.rb = if is_rev { (l_pac << 1) - (rb + aln.te as i64 + 1) } else { rb + aln.tb as i64 };
                b.re = if is_rev { (l_pac << 1) - (rb + aln.tb as i64) } else { rb + aln.te as i64 + 1 };
                b.seedcov = (((b.re - b.rb) as i32).min(b.qe - b.qb)) >> 1;
                added = Some(b);
            }
            n += 1;
        }
        if n > 0 {
            rescue_dedup(o, bns, ma, added, clean);
        }
    }
    n
}

/// mem_pair: (o, subo, n_sub, z)
fn pair(o: &SamOpt, bns: &Bns, pes: &[Pestat; 4], a: &[Vec<Reg>; 2], id: i64, n_pri: [usize; 2]) -> (i32, i32, i32, [usize; 2]) {
    let l_pac = bns.l_pac;
    let mut v: Vec<(u64, u64)> = Vec::new();
    for r in 0..2 {
        for i in 0..n_pri[r] {
            let e = &a[r][i];
            let mut x: u64 = if e.rb < l_pac { e.rb as u64 } else { ((l_pac << 1) - 1 - e.rb) as u64 };
            x = ((e.rid as u64) << 32) | (x.wrapping_sub(bns.anns[e.rid as usize].offset as u64));
            let y = ((e.score as u64) << 32) | ((i as u64) << 2) | (((e.rb >= l_pac) as u64) << 1) | r as u64;
            v.push((x, y));
        }
    }
    v.sort_unstable(); // ks_introsort_128 on (x, y): unique keys
    let mut u: Vec<(u64, u64)> = Vec::new();
    let mut y = [-1i64; 4];
    for i in 0..v.len() {
        for r in 0..2 {
            let dir = (r << 1) | ((v[i].1 >> 1) & 1) as usize;
            if pes[dir].failed {
                continue;
            }
            let which = (r << 1) | ((v[i].1 & 1) ^ 1) as usize;
            if y[which] < 0 {
                continue;
            }
            let mut k = y[which];
            while k >= 0 {
                let ku = k as usize;
                if (v[ku].1 & 3) as usize != which {
                    k -= 1;
                    continue;
                }
                let dist = (v[i].0 as i64).wrapping_sub(v[ku].0 as i64);
                if dist > pes[dir].high as i64 {
                    break;
                }
                if dist < pes[dir].low as i64 {
                    k -= 1;
                    continue;
                }
                let ns = (dist as f64 - pes[dir].avg) / pes[dir].std;
                let mut q = (((v[i].1 >> 32) + (v[ku].1 >> 32)) as f64 + 0.721 * c_log(2.0 * unsafe { erfc(ns.abs() * std::f64::consts::FRAC_1_SQRT_2) }) * o.e.a as f64 + 0.499) as i32;
                if q < 0 {
                    q = 0;
                }
                let py = ((k as u64) << 32) | i as u64;
                // mem_pair's `int id`: the pair id truncated to 32 bits, `id<<8` computed in 32 bits (wrapping),
                // then sign-extended for the XOR with the 64-bit key
                let id8 = ((id as i32 as u32).wrapping_shl(8)) as i32 as i64 as u64;
                let px = ((q as u64) << 32) | (hash_64(py ^ id8) & 0xffff_ffff);
                u.push((px, py));
                k -= 1;
            }
        }
        y[(v[i].1 & 3) as usize] = i as i64;
    }
    if u.is_empty() {
        return (0, 0, 0, [0, 0]);
    }
    let mut tmp = o.e.a + o.b;
    tmp = tmp.max(o.e.o_del + o.e.e_del);
    tmp = tmp.max(o.e.o_ins + o.e.e_ins);
    u.sort_unstable();
    let last = u[u.len() - 1];
    let i = (last.1 >> 32) as usize;
    let k = (last.1 & 0xffff_ffff) as usize;
    let mut z = [0usize; 2];
    z[(v[i].1 & 1) as usize] = ((v[i].1 << 32) >> 34) as usize;
    z[(v[k].1 & 1) as usize] = ((v[k].1 << 32) >> 34) as usize;
    let ret = (last.0 >> 32) as i32;
    let sub = if u.len() > 1 { (u[u.len() - 2].0 >> 32) as i32 } else { 0 };
    let mut n_sub = 0;
    for j in (0..u.len() - 1).rev() {
        if sub - ((u[j].0 >> 32) as i32) <= tmp {
            n_sub += 1;
        }
    }
    (ret, sub, n_sub, z)
}

fn raw_mapq(diff: i32, a: i32) -> i32 {
    (6.02 * diff as f64 / a as f64 + 0.499) as i32
}

/// mem_sam_pe: SAM text of both reads of a pair
pub fn sam_pe(o: &SamOpt, bns: &Bns, pes: &[Pestat; 4], id: i64, s: [&Read; 2], mut a: [Vec<Reg>; 2], pre: Option<&[PreSw]>) -> (Vec<u8>, Vec<u8>) {
    // mate rescue
    let b: [Vec<Reg>; 2] = [
        a[0].iter().filter(|r| r.score >= a[0][0].score - o.pen_unpaired).cloned().collect(),
        a[1].iter().filter(|r| r.score >= a[1][0].score - o.pen_unpaired).cloned().collect(),
    ];
    crate::prof::cpu(crate::prof::C_PRESCUE, || {
        // per mate: its regions are clean for ext::dedup_add (none are before the first full dedup)
        let mut clean = [false; 2];
        for i in 0..2 {
            for j in 0..b[i].len().min(o.max_matesw as usize) {
                let (left, right) = a.split_at_mut(1);
                let ma = if i == 0 { &mut right[0] } else { &mut left[0] };
                matesw(o, bns, pes, &b[i][j], s[1 - i].seq, ma, &mut clean[1 - i], pre, i, j);
            }
        }
    });
    let n_pri = crate::prof::cpu(crate::prof::C_PMARK, || [mark_primary_se(o, &mut a[0], id << 1), mark_primary_se(o, &mut a[1], (id << 1) | 1)]);
    if let Some(d) = PDUMP.get() {
        let mut t = String::new();
        for i in 0..2 {
            for (j, p) in a[i].iter().enumerate() {
                t.push_str(&format!("PREG\t{}\t{i}\t{j}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", String::from_utf8_lossy(s[0].name), p.rid, p.rb, p.re, p.qb, p.qe, p.score, p.truesc, p.sub, p.csub, p.sub_n, p.secondary, p.secondary_all, p.is_alt, p.seedcov, n_pri[i]));
            }
        }
        use std::io::Write;
        d.lock().unwrap().write_all(t.as_bytes()).unwrap();
    }
    let mut extra_flag = 1;
    if n_pri[0] > 0 && n_pri[1] > 0 {
        let (op, mut subo, n_sub, mut z) = crate::prof::cpu(crate::prof::C_PPAIR, || pair(o, bns, pes, &a, id, n_pri));
        if op > 0 {
            let mut multi = false;
            for i in 0..2 {
                let mut j = 1;
                while j < n_pri[i] {
                    if a[i][j].secondary < 0 && a[i][j].score >= o.t {
                        break;
                    }
                    j += 1;
                }
                if j < n_pri[i] {
                    multi = true;
                }
            }
            if !multi {
                let score_un = a[0][0].score + a[1][0].score - o.pen_unpaired;
                subo = subo.max(score_un);
                let mut q_pe = raw_mapq(op - subo, o.e.a);
                if n_sub > 0 {
                    q_pe -= (4.343 * c_log((n_sub + 1) as f64) + 0.499) as i32;
                }
                q_pe = q_pe.clamp(0, 60);
                q_pe = ((q_pe as f64) * (1.0 - 0.5 * (a[0][0].frac_rep as f64 + a[1][0].frac_rep as f64)) + 0.499) as i32;
                let mut q_se = [0i32; 2];
                if op > score_un {
                    for i in 0..2 {
                        let zi = z[i];
                        if a[i][zi].secondary >= 0 {
                            let sec = a[i][zi].secondary as usize;
                            a[i][zi].sub = a[i][sec].score;
                            a[i][zi].secondary = -2;
                        }
                        q_se[i] = approx_mapq_se(o, &a[i][zi]);
                    }
                    for i in 0..2 {
                        q_se[i] = if q_se[i] > q_pe { q_se[i] } else if q_pe < q_se[i] + 40 { q_pe } else { q_se[i] + 40 };
                    }
                    extra_flag |= 2;
                    for i in 0..2 {
                        let c = &a[i][z[i]];
                        q_se[i] = q_se[i].min(raw_mapq(c.score - c.csub, o.e.a));
                    }
                } else {
                    z = [0, 0];
                    q_se[0] = approx_mapq_se(o, &a[0][0]);
                    q_se[1] = approx_mapq_se(o, &a[1][0]);
                }
                for i in 0..2 {
                    let k = a[i][z[i]].secondary_all;
                    if k >= 0 && (k as usize) < n_pri[i] {
                        assert!(a[i][k as usize].secondary_all < 0);
                        for j in 0..a[i].len() {
                            if a[i][j].secondary_all == k || j == k as usize {
                                a[i][j].secondary_all = z[i] as i32;
                            }
                        }
                        a[i][z[i]].secondary_all = -1;
                    }
                }
                let xa = crate::prof::cpu(crate::prof::C_PXA, || [gen_alt(o, bns, &a[0], s[0].seq), gen_alt(o, bns, &a[1], s[1].seq)]);
                let mut h: [Aln; 2] = Default::default();
                let mut aa: [Vec<Aln>; 2] = Default::default();
                for i in 0..2 {
                    h[i] = crate::prof::cpu(crate::prof::C_PALN, || reg2aln(o, bns, s[i].seq, Some(&a[i][z[i]])));
                    h[i].mapq = q_se[i];
                    h[i].flag |= (0x40 << i) | extra_flag;
                    h[i].xa = xa[i].as_ref().and_then(|x| x[z[i]].clone());
                    aa[i].push(h[i].clone());
                    if n_pri[i] < a[i].len() {
                        let p = &a[i][n_pri[i]];
                        if p.score < o.t || p.secondary >= 0 || p.is_alt == 0 {
                            continue;
                        }
                        let mut g = reg2aln(o, bns, s[i].seq, Some(p));
                        g.flag |= 0x800 | (0x40 << i) | extra_flag;
                        g.xa = xa[i].as_ref().and_then(|x| x[n_pri[i]].clone());
                        aa[i].push(g);
                    }
                }
                let mut s0 = Vec::new();
                for k in 0..aa[0].len() {
                    aln2sam(o, bns, &mut s0, s[0], &aa[0], k, Some(&h[1]));
                }
                let mut s1 = Vec::new();
                for k in 0..aa[1].len() {
                    aln2sam(o, bns, &mut s1, s[1], &aa[1], k, Some(&h[0]));
                }
                assert_eq!(s[0].name, s[1].name, "paired reads have different names");
                return (s0, s1);
            }
        }
    }
    // no_pairing
    let mut h: [Aln; 2] = Default::default();
    for i in 0..2 {
        let mut which: i64 = -1;
        if !a[i].is_empty() {
            if a[i][0].score >= o.t {
                which = 0;
            } else if n_pri[i] < a[i].len() && a[i][n_pri[i]].score >= o.t {
                which = n_pri[i] as i64;
            }
        }
        h[i] = if which >= 0 { reg2aln(o, bns, s[i].seq, Some(&a[i][which as usize])) } else { reg2aln(o, bns, s[i].seq, None) };
    }
    if h[0].rid == h[1].rid && h[0].rid >= 0 {
        let (d, dist) = infer_dir(bns.l_pac, a[0][0].rb, a[1][0].rb);
        if !pes[d].failed && dist >= pes[d].low as i64 && dist <= pes[d].high as i64 {
            extra_flag |= 2;
        }
    }
    let mut s0 = Vec::new();
    let mut s1 = Vec::new();
    reg2sam(o, bns, &mut s0, s[0], &a[0], 0x41 | extra_flag, Some(&h[1]));
    reg2sam(o, bns, &mut s1, s[1], &a[1], 0x81 | extra_flag, Some(&h[0]));
    assert_eq!(s[0].name, s[1].name, "paired reads have different names");
    (s0, s1)
}

/// single-end: mem_mark_primary_se + mem_reg2sam
pub fn sam_se(o: &SamOpt, bns: &Bns, id: i64, s: &Read, mut a: Vec<Reg>) -> Vec<u8> {
    mark_primary_se(o, &mut a, id);
    let mut out = Vec::new();
    reg2sam(o, bns, &mut out, s, &a, 0, None);
    out
}
