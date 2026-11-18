//! CPU PairHMM: a port of Intel GKL's `compute_full_prob` (the PairHMM GATK uses by default,
//! `VectorLoglessPairHMM` -> `IntelPairHmm`), used when no GPU is usable. `compute` is the scalar reference;
//! `compute_batched` runs 16 (read, haplotype) pairs per AVX-512 float vector (8 per double vector for the
//! double-precision redo), each lane doing exactly the scalar operation sequence, so results are bit-identical
//! (~11 Gcells/s per core vs 0.7 scalar on Zen 5; `hmmbench` measures and checks this).
//!
//! It reproduces GKL's numerics operation by operation: the same probability tables (built with
//! the same libm calls), float pass with the 2^120 scale and a double pass (2^1020) when the float
//! result is < 1e-28, flush-to-zero, the per-cell operation order, the separate M and X sums of the
//! last row, and `log10f(sum) - log10f(2^120)`.  GKL's AVX-512 build (what GATK uses on AVX-512
//! CPUs) is compiled with FMA contraction: each `a*b + c*d` becomes `fma(c, d, a*b)`.  The AVX2
//! build has no FMA.  [`Flavor`] selects which one to reproduce.

use crate::{HmmGroup, HmmRead};

/// Which GKL build to reproduce.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Flavor {
    /// `compute_full_prob_avx512{s,d}` (FMA-contracted)
    Avx512,
    /// `compute_full_prob_avx{s,d}` (separate multiply and add)
    Avx2,
}

const MAX_QUAL: usize = 254;
const NQ: usize = 128;
const MM_SIZE: usize = (NQ * (NQ + 1)) / 2;
const MIN_ACCEPTED: f32 = 1e-28;

trait Num: Copy + PartialOrd + std::ops::Add<Output = Self> + std::ops::Sub<Output = Self> + std::ops::Mul<Output = Self> + std::ops::Div<Output = Self> {
    const ZERO: Self;
    const ONE: Self;
    fn f(x: f64) -> Self;
    fn d(self) -> f64;
    fn fi(i: i32) -> Self;
    fn fma(a: Self, b: Self, c: Self) -> Self;
    fn is_inf(self) -> bool;
    fn round_half(self) -> i32;
}
impl Num for f32 {
    const ZERO: f32 = 0.0;
    const ONE: f32 = 1.0;
    #[inline(always)]
    fn f(x: f64) -> f32 {
        x as f32
    }
    #[inline(always)]
    fn d(self) -> f64 {
        self as f64
    }
    #[inline(always)]
    fn fi(i: i32) -> f32 {
        i as f32
    }
    #[inline(always)]
    fn fma(a: f32, b: f32, c: f32) -> f32 {
        a.mul_add(b, c)
    }
    fn is_inf(self) -> bool {
        self.is_infinite()
    }
    fn round_half(self) -> i32 {
        if self > 0.0 { (self + 0.5) as i32 } else { (self - 0.5) as i32 }
    }
}
impl Num for f64 {
    const ZERO: f64 = 0.0;
    const ONE: f64 = 1.0;
    #[inline(always)]
    fn f(x: f64) -> f64 {
        x
    }
    #[inline(always)]
    fn d(self) -> f64 {
        self
    }
    #[inline(always)]
    fn fi(i: i32) -> f64 {
        i as f64
    }
    #[inline(always)]
    fn fma(a: f64, b: f64, c: f64) -> f64 {
        a.mul_add(b, c)
    }
    fn is_inf(self) -> bool {
        self.is_infinite()
    }
    fn round_half(self) -> i32 {
        if self > 0.0 { (self + 0.5) as i32 } else { (self - 0.5) as i32 }
    }
}

/// GKL `Context<NUMBER>` tables.
struct Tables<T> {
    ph2pr: [T; NQ],
    mm: Vec<T>,
    init: T,
}

fn bb<X>(x: X) -> X {
    std::hint::black_box(x) // keep libm calls at run time (no compile-time folding)
}

impl<T: Num> Tables<T> {
    fn new(ph2pr: [T; NQ], init: T) -> Self {
        // ContextBase::initializeJacobianLogTable
        let size = (8.0f64 / 0.0001) as usize + 1;
        let jac: Vec<T> = (0..size).map(|k| T::f((1.0 + 10f64.powf(bb(-(k as f64) * 0.0001))).log10())).collect();
        let inv_step = T::f(1.0f64 / 0.0001);
        let approx = |mut small: T, mut big: T| -> T {
            if small > big {
                std::mem::swap(&mut small, &mut big);
            }
            if small.is_inf() || big.is_inf() {
                return big;
            }
            let diff = big - small;
            if diff >= T::f(8.0) {
                return big;
            }
            let ind = (diff * inv_step).round_half();
            big + jac[ind as usize]
        };
        // ContextBase::initializeMatchToMatchProb (only qualities < 128 are ever used)
        let inv_ln10 = 0.434294f64;
        let mut mm = vec![T::ZERO; MM_SIZE];
        let mut offset = 0usize;
        for i in 0..=MAX_QUAL {
            for j in 0..=i {
                let log10sum = approx(T::f(-0.1) * T::fi(i as i32), T::f(-0.1) * T::fi(j as i32)).d();
                let mm10 = (-(1.0f64.min(10f64.powf(bb(log10sum))))).ln_1p() * inv_ln10;
                if offset + j < MM_SIZE {
                    mm[offset + j] = T::f(10f64.powf(bb(mm10)));
                }
            }
            offset += i + 1;
        }
        Tables { ph2pr, mm, init }
    }
}

struct Ctx {
    f: Tables<f32>,
    d: Tables<f64>,
    log10_init_f: f32,
    log10_init_d: f64,
}

fn ctx() -> &'static Ctx {
    static C: std::sync::OnceLock<Ctx> = std::sync::OnceLock::new();
    C.get_or_init(|| {
        let mut pf = [0f32; NQ];
        let mut pd = [0f64; NQ];
        for x in 0..NQ {
            pf[x] = 10f32.powf(bb(-(x as f32) / 10.0));
            pd[x] = 10f64.powf(bb(-(x as f64) / 10.0));
        }
        let init_f = 2f32.powi(120);
        let init_d = 2f64.powi(1020);
        Ctx { f: Tables::new(pf, init_f), d: Tables::new(pd, init_d), log10_init_f: bb(init_f).log10(), log10_init_d: bb(init_d).log10() }
    })
}

/// GKL `ConvertChar`: A0 C1 T2 G3 N4, anything else 0.
#[inline(always)]
fn conv(b: u8) -> u8 {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'T' => 2,
        b'G' => 3,
        b'N' => 4,
        _ => 0,
    }
}

/// Reusable DP buffers.
#[derive(Default)]
pub struct Workspace {
    f: Rows<f32>,
    d: Rows<f64>,
    hap: Vec<u8>,
}
struct Rows<T> {
    m: Vec<T>,
    x: Vec<T>,
    y: Vec<T>,
}
impl<T> Default for Rows<T> {
    fn default() -> Self {
        Rows { m: Vec::new(), x: Vec::new(), y: Vec::new() }
    }
}

#[inline(always)]
fn full_prob<T: Num, const FMA: bool>(t: &Tables<T>, read: &HmmRead, hap: &[u8], w: &mut Rows<T>) -> T {
    let r_len = read.bases.len();
    let h = hap.len();
    let init_y = t.init / T::fi(h as i32);
    // row 0: M = X = 0, Y = initY (also at column 0, the diagonal of cell (1,1))
    w.m.clear();
    w.m.resize(h + 1, T::ZERO);
    w.x.clear();
    w.x.resize(h + 1, T::ZERO);
    w.y.clear();
    w.y.resize(h + 1, init_y);
    let (mut sum_m, mut sum_x) = (T::ZERO, T::ZERO);
    for i in 1..=r_len {
        let qi = (read.ins[i - 1] & 127) as usize;
        let qd = (read.del[i - 1] & 127) as usize;
        let qc = (read.gcp[i - 1] & 127) as usize;
        let qq = (read.quals[i - 1] & 127) as usize;
        let (mx, mn) = if qi <= qd { (qd, qi) } else { (qi, qd) };
        let p_mm = t.mm[((mx * (mx + 1)) >> 1) + mn];
        let p_gapm = T::ONE - t.ph2pr[qc];
        let p_mx = t.ph2pr[qi];
        let p_xx = t.ph2pr[qc];
        let p_my = t.ph2pr[qd];
        let p_yy = t.ph2pr[qc];
        let distm = t.ph2pr[qq];
        let one_m = T::ONE - distm;
        let distm3 = distm / T::f(3.0);
        let rc = conv(read.bases[i - 1]);
        // walk the row in place: (m2, x2, y2) = previous row at j-1 (diagonal)
        let (mut m2, mut x2, mut y2) = (w.m[0], w.x[0], w.y[0]);
        let (mut m_left, mut y_left) = (T::ZERO, T::ZERO); // column 0
        w.m[0] = T::ZERO;
        w.x[0] = T::ZERO;
        w.y[0] = T::ZERO;
        for j in 1..=h {
            let hc = hap[j - 1];
            let prior = if rc == hc || rc == 4 || hc == 4 { one_m } else { distm3 };
            let (m_up, x_up, y_up) = (w.m[j], w.x[j], w.y[j]);
            let (m, x, y);
            if FMA {
                m = T::fma(y2, p_gapm, T::fma(x2, p_gapm, m2 * p_mm)) * prior;
                x = T::fma(x_up, p_xx, m_up * p_mx);
                y = T::fma(y_left, p_yy, m_left * p_my);
            } else {
                m = ((m2 * p_mm + x2 * p_gapm) + y2 * p_gapm) * prior;
                x = m_up * p_mx + x_up * p_xx;
                y = m_left * p_my + y_left * p_yy;
            }
            w.m[j] = m;
            w.x[j] = x;
            w.y[j] = y;
            m2 = m_up;
            x2 = x_up;
            y2 = y_up;
            m_left = m;
            y_left = y;
        }
        if i == r_len {
            for j in 1..=h {
                sum_m = sum_m + w.m[j];
                sum_x = sum_x + w.x[j];
            }
        }
    }
    sum_m + sum_x
}

struct FtzGuard(u32);
impl FtzGuard {
    #[allow(deprecated)]
    fn on() -> Self {
        #[cfg(target_arch = "x86_64")]
        unsafe {
            use std::arch::x86_64::{_mm_getcsr, _mm_setcsr};
            let old = _mm_getcsr();
            _mm_setcsr(old | 0x8000); // flush-to-zero (as GKL's initNative)
            FtzGuard(old)
        }
        #[cfg(not(target_arch = "x86_64"))]
        FtzGuard(0)
    }
}
impl Drop for FtzGuard {
    #[allow(deprecated)]
    fn drop(&mut self) {
        #[cfg(target_arch = "x86_64")]
        unsafe {
            std::arch::x86_64::_mm_setcsr(self.0)
        }
    }
}

/// log10 likelihood of one (read, haplotype) pair, as `IntelPairHmm.computeLikelihoodsNative`.
pub fn pair_log10(read: &HmmRead, hap: &[u8], flavor: Flavor, w: &mut Workspace) -> f64 {
    let c = ctx();
    let _g = FtzGuard::on();
    if read.bases.is_empty() {
        return 0.0;
    }
    w.hap.clear();
    w.hap.extend(hap.iter().map(|&b| conv(b)));
    let hap = std::mem::take(&mut w.hap);
    let fl = match flavor {
        Flavor::Avx512 => full_prob::<f32, true>(&c.f, read, &hap, &mut w.f),
        Flavor::Avx2 => full_prob::<f32, false>(&c.f, read, &hap, &mut w.f),
    };
    let res = if fl < MIN_ACCEPTED {
        let dv = match flavor {
            Flavor::Avx512 => full_prob::<f64, true>(&c.d, read, &hap, &mut w.d),
            Flavor::Avx2 => full_prob::<f64, false>(&c.d, read, &hap, &mut w.d),
        };
        dv.log10() - c.log10_init_d
    } else {
        (fl.log10() - c.log10_init_f) as f64
    };
    w.hap = hap;
    res
}

/// Same contract as `hc_cudahmm::compute`: per group, a row-major (read x haplotype) matrix of
/// log10 likelihoods.  Single-threaded; parallelise across groups in the caller.
pub fn compute(groups: &[HmmGroup], flavor: Flavor) -> Vec<Vec<f64>> {
    let mut w = Workspace::default();
    groups
        .iter()
        .map(|g| {
            let mut out = Vec::with_capacity(g.reads.len() * g.haps.len());
            for r in &g.reads {
                for h in &g.haps {
                    out.push(pair_log10(r, h, flavor, &mut w));
                }
            }
            out
        })
        .collect()
}

// ------------------------------------------------------------------------------------------------
// Lane-batched version: 16 (read, haplotype) pairs advance through their DP matrices together
// (one SIMD lane each).  Every lane performs exactly the scalar operation sequence, so results are
// bit-identical to `pair_log10`; lanes only differ in which cells / rows they read out.
// ------------------------------------------------------------------------------------------------

const L: usize = 16;
type V = [f32; L];

#[derive(Clone, Copy)]
struct PairRef<'a> {
    read: &'a HmmRead<'a>,
    hap: &'a [u8],
    out: usize,
}

#[inline(always)]
fn tile_float<const FMA: bool>(t: &Tables<f32>, lanes: &[PairRef], ws: &mut TileWs) -> [f32; L] {
    let n = lanes.len();
    let r_max = lanes.iter().map(|p| p.read.bases.len()).max().unwrap();
    let h_max = lanes.iter().map(|p| p.hap.len()).max().unwrap();
    // haplotype codes, column-major by lane; padding never matters (cells right of a lane's H
    // are never read by that lane's own cells)
    ws.hap.clear();
    ws.hap.resize(h_max, [0u8; L]);
    let mut init = [0f32; L];
    for (l, p) in lanes.iter().enumerate() {
        for (j, &b) in p.hap.iter().enumerate() {
            ws.hap[j][l] = conv(b);
        }
        init[l] = t.init / (p.hap.len() as f32);
    }
    ws.m.clear();
    ws.m.resize(h_max + 1, [0f32; L]);
    ws.x.clear();
    ws.x.resize(h_max + 1, [0f32; L]);
    ws.y.clear();
    ws.y.resize(h_max + 1, init);
    let mut result = [0f32; L];
    for i in 1..=r_max {
        let (mut pmm, mut pgapm, mut pmx, mut pxx, mut pmy, mut pyy, mut one_m, mut d3) = ([0f32; L], [0f32; L], [0f32; L], [0f32; L], [0f32; L], [0f32; L], [0f32; L], [0f32; L]);
        let mut rc = [255u8; L];
        for (l, p) in lanes.iter().enumerate() {
            let r = p.read;
            if i > r.bases.len() {
                continue;
            }
            let qi = (r.ins[i - 1] & 127) as usize;
            let qd = (r.del[i - 1] & 127) as usize;
            let qc = (r.gcp[i - 1] & 127) as usize;
            let qq = (r.quals[i - 1] & 127) as usize;
            let (mx, mn) = if qi <= qd { (qd, qi) } else { (qi, qd) };
            pmm[l] = t.mm[((mx * (mx + 1)) >> 1) + mn];
            pgapm[l] = 1.0 - t.ph2pr[qc];
            pmx[l] = t.ph2pr[qi];
            pxx[l] = t.ph2pr[qc];
            pmy[l] = t.ph2pr[qd];
            pyy[l] = t.ph2pr[qc];
            let distm = t.ph2pr[qq];
            one_m[l] = 1.0 - distm;
            d3[l] = distm / 3.0;
            rc[l] = conv(r.bases[i - 1]);
        }
        let rc_n: [bool; L] = std::array::from_fn(|l| rc[l] == 4);
        let mut m2 = ws.m[0];
        let mut x2 = ws.x[0];
        let mut y2 = ws.y[0];
        let mut m_left = [0f32; L];
        let mut y_left = [0f32; L];
        ws.m[0] = [0f32; L];
        ws.x[0] = [0f32; L];
        ws.y[0] = [0f32; L];
        for j in 1..=h_max {
            let hc = ws.hap[j - 1];
            let m_up = ws.m[j];
            let x_up = ws.x[j];
            let y_up = ws.y[j];
            let mut m = [0f32; L];
            let mut x = [0f32; L];
            let mut y = [0f32; L];
            for l in 0..L {
                let prior = if rc[l] == hc[l] || rc_n[l] || hc[l] == 4 { one_m[l] } else { d3[l] };
                if FMA {
                    m[l] = y2[l].mul_add(pgapm[l], x2[l].mul_add(pgapm[l], m2[l] * pmm[l])) * prior;
                    x[l] = x_up[l].mul_add(pxx[l], m_up[l] * pmx[l]);
                    y[l] = y_left[l].mul_add(pyy[l], m_left[l] * pmy[l]);
                } else {
                    m[l] = ((m2[l] * pmm[l] + x2[l] * pgapm[l]) + y2[l] * pgapm[l]) * prior;
                    x[l] = m_up[l] * pmx[l] + x_up[l] * pxx[l];
                    y[l] = m_left[l] * pmy[l] + y_left[l] * pyy[l];
                }
            }
            ws.m[j] = m;
            ws.x[j] = x;
            ws.y[j] = y;
            m2 = m_up;
            x2 = x_up;
            y2 = y_up;
            m_left = m;
            y_left = y;
        }
        for (l, p) in lanes.iter().enumerate().take(n) {
            if p.read.bases.len() == i {
                let (mut sm, mut sx) = (0f32, 0f32);
                for j in 1..=p.hap.len() {
                    sm += ws.m[j][l];
                    sx += ws.x[j][l];
                }
                result[l] = sm + sx;
            }
        }
    }
    result
}

/// [`tile_float`] (FMA flavor) with the column sweep in AVX-512 registers. Same per-lane operation sequence, so
/// bit-identical results.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn tile_float_avx512(t: &Tables<f32>, lanes: &[PairRef], ws: &mut TileWs) -> [f32; L] {
    use std::arch::x86_64::*;
    let n = lanes.len();
    let r_max = lanes.iter().map(|p| p.read.bases.len()).max().unwrap();
    let h_max = lanes.iter().map(|p| p.hap.len()).max().unwrap();
    ws.hap.clear();
    ws.hap.resize(h_max, [0u8; L]);
    let mut init = [0f32; L];
    for (l, p) in lanes.iter().enumerate() {
        for (j, &b) in p.hap.iter().enumerate() {
            ws.hap[j][l] = conv(b);
        }
        init[l] = t.init / (p.hap.len() as f32);
    }
    ws.m.clear();
    ws.m.resize(h_max + 1, [0f32; L]);
    ws.x.clear();
    ws.x.resize(h_max + 1, [0f32; L]);
    ws.y.clear();
    ws.y.resize(h_max + 1, init);
    let mut result = [0f32; L];
    let four = _mm512_set1_epi32(4);
    let zero = _mm512_setzero_ps();
    for i in 1..=r_max {
        let (mut pmm, mut pgapm, mut pmx, mut pxx, mut pmy, mut pyy, mut one_m, mut d3) = ([0f32; L], [0f32; L], [0f32; L], [0f32; L], [0f32; L], [0f32; L], [0f32; L], [0f32; L]);
        let mut rc = [255i32; L];
        for (l, p) in lanes.iter().enumerate() {
            let r = p.read;
            if i > r.bases.len() {
                continue;
            }
            let qi = (r.ins[i - 1] & 127) as usize;
            let qd = (r.del[i - 1] & 127) as usize;
            let qc = (r.gcp[i - 1] & 127) as usize;
            let qq = (r.quals[i - 1] & 127) as usize;
            let (mx, mn) = if qi <= qd { (qd, qi) } else { (qi, qd) };
            pmm[l] = t.mm[((mx * (mx + 1)) >> 1) + mn];
            pgapm[l] = 1.0 - t.ph2pr[qc];
            pmx[l] = t.ph2pr[qi];
            pxx[l] = t.ph2pr[qc];
            pmy[l] = t.ph2pr[qd];
            pyy[l] = t.ph2pr[qc];
            let distm = t.ph2pr[qq];
            one_m[l] = 1.0 - distm;
            d3[l] = distm / 3.0;
            rc[l] = conv(r.bases[i - 1]) as i32;
        }
        let (pmm, pgapm, pmx, pxx, pmy, pyy, one_m, d3) = (
            _mm512_loadu_ps(pmm.as_ptr()),
            _mm512_loadu_ps(pgapm.as_ptr()),
            _mm512_loadu_ps(pmx.as_ptr()),
            _mm512_loadu_ps(pxx.as_ptr()),
            _mm512_loadu_ps(pmy.as_ptr()),
            _mm512_loadu_ps(pyy.as_ptr()),
            _mm512_loadu_ps(one_m.as_ptr()),
            _mm512_loadu_ps(d3.as_ptr()),
        );
        let rcv = _mm512_loadu_si512(rc.as_ptr() as *const _);
        let rc_n = _mm512_cmpeq_epi32_mask(rcv, four);
        let (mp, xp, yp, hp) = (ws.m.as_mut_ptr() as *mut f32, ws.x.as_mut_ptr() as *mut f32, ws.y.as_mut_ptr() as *mut f32, ws.hap.as_ptr() as *const u8);
        let mut m2 = _mm512_loadu_ps(mp);
        let mut x2 = _mm512_loadu_ps(xp);
        let mut y2 = _mm512_loadu_ps(yp);
        let mut m_left = zero;
        let mut y_left = zero;
        _mm512_storeu_ps(mp, zero);
        _mm512_storeu_ps(xp, zero);
        _mm512_storeu_ps(yp, zero);
        for j in 1..=h_max {
            let hc = _mm512_cvtepu8_epi32(_mm_loadu_si128(hp.add((j - 1) * L) as *const __m128i));
            let mask = _mm512_cmpeq_epi32_mask(rcv, hc) | rc_n | _mm512_cmpeq_epi32_mask(hc, four);
            let prior = _mm512_mask_blend_ps(mask, d3, one_m);
            let (pm, px, py) = (mp.add(j * L), xp.add(j * L), yp.add(j * L));
            let m_up = _mm512_loadu_ps(pm);
            let x_up = _mm512_loadu_ps(px);
            let y_up = _mm512_loadu_ps(py);
            // m = fma(y2, pgapm, fma(x2, pgapm, m2 * pmm)) * prior
            let m = _mm512_mul_ps(_mm512_fmadd_ps(y2, pgapm, _mm512_fmadd_ps(x2, pgapm, _mm512_mul_ps(m2, pmm))), prior);
            // x = fma(x_up, pxx, m_up * pmx); y = fma(y_left, pyy, m_left * pmy)
            let x = _mm512_fmadd_ps(x_up, pxx, _mm512_mul_ps(m_up, pmx));
            let y = _mm512_fmadd_ps(y_left, pyy, _mm512_mul_ps(m_left, pmy));
            _mm512_storeu_ps(pm, m);
            _mm512_storeu_ps(px, x);
            _mm512_storeu_ps(py, y);
            m2 = m_up;
            x2 = x_up;
            y2 = y_up;
            m_left = m;
            y_left = y;
        }
        for (l, p) in lanes.iter().enumerate().take(n) {
            if p.read.bases.len() == i {
                let (mut sm, mut sx) = (0f32, 0f32);
                for j in 1..=p.hap.len() {
                    sm += ws.m[j][l];
                    sx += ws.x[j][l];
                }
                result[l] = sm + sx;
            }
        }
    }
    result
}

const LD: usize = 8;

#[derive(Default)]
struct TileWsD {
    m: Vec<[f64; LD]>,
    x: Vec<[f64; LD]>,
    y: Vec<[f64; LD]>,
    hap: Vec<[u8; LD]>,
}

/// The double-precision redo (`full_prob::<f64, true>`) for up to 8 pairs at once, in AVX-512 registers. Same
/// per-lane operation sequence, so bit-identical results.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn tile_double_avx512(t: &Tables<f64>, lanes: &[PairRef], ws: &mut TileWsD) -> [f64; LD] {
    use std::arch::x86_64::*;
    let n = lanes.len();
    let r_max = lanes.iter().map(|p| p.read.bases.len()).max().unwrap();
    let h_max = lanes.iter().map(|p| p.hap.len()).max().unwrap();
    ws.hap.clear();
    ws.hap.resize(h_max, [0u8; LD]);
    let mut init = [0f64; LD];
    for (l, p) in lanes.iter().enumerate() {
        for (j, &b) in p.hap.iter().enumerate() {
            ws.hap[j][l] = conv(b);
        }
        init[l] = t.init / (p.hap.len() as i32 as f64);
    }
    ws.m.clear();
    ws.m.resize(h_max + 1, [0f64; LD]);
    ws.x.clear();
    ws.x.resize(h_max + 1, [0f64; LD]);
    ws.y.clear();
    ws.y.resize(h_max + 1, init);
    let mut result = [0f64; LD];
    let four = _mm512_set1_epi64(4);
    let zero = _mm512_setzero_pd();
    for i in 1..=r_max {
        let (mut pmm, mut pgapm, mut pmx, mut pxx, mut pmy, mut pyy, mut one_m, mut d3) = ([0f64; LD], [0f64; LD], [0f64; LD], [0f64; LD], [0f64; LD], [0f64; LD], [0f64; LD], [0f64; LD]);
        let mut rc = [255i64; LD];
        for (l, p) in lanes.iter().enumerate() {
            let r = p.read;
            if i > r.bases.len() {
                continue;
            }
            let qi = (r.ins[i - 1] & 127) as usize;
            let qd = (r.del[i - 1] & 127) as usize;
            let qc = (r.gcp[i - 1] & 127) as usize;
            let qq = (r.quals[i - 1] & 127) as usize;
            let (mx, mn) = if qi <= qd { (qd, qi) } else { (qi, qd) };
            pmm[l] = t.mm[((mx * (mx + 1)) >> 1) + mn];
            pgapm[l] = 1.0 - t.ph2pr[qc];
            pmx[l] = t.ph2pr[qi];
            pxx[l] = t.ph2pr[qc];
            pmy[l] = t.ph2pr[qd];
            pyy[l] = t.ph2pr[qc];
            let distm = t.ph2pr[qq];
            one_m[l] = 1.0 - distm;
            d3[l] = distm / 3.0;
            rc[l] = conv(r.bases[i - 1]) as i64;
        }
        let (pmm, pgapm, pmx, pxx, pmy, pyy, one_m, d3) = (
            _mm512_loadu_pd(pmm.as_ptr()),
            _mm512_loadu_pd(pgapm.as_ptr()),
            _mm512_loadu_pd(pmx.as_ptr()),
            _mm512_loadu_pd(pxx.as_ptr()),
            _mm512_loadu_pd(pmy.as_ptr()),
            _mm512_loadu_pd(pyy.as_ptr()),
            _mm512_loadu_pd(one_m.as_ptr()),
            _mm512_loadu_pd(d3.as_ptr()),
        );
        let rcv = _mm512_loadu_si512(rc.as_ptr() as *const _);
        let rc_n = _mm512_cmpeq_epi64_mask(rcv, four);
        let (mp, xp, yp, hp) = (ws.m.as_mut_ptr() as *mut f64, ws.x.as_mut_ptr() as *mut f64, ws.y.as_mut_ptr() as *mut f64, ws.hap.as_ptr() as *const u8);
        let mut m2 = _mm512_loadu_pd(mp);
        let mut x2 = _mm512_loadu_pd(xp);
        let mut y2 = _mm512_loadu_pd(yp);
        let mut m_left = zero;
        let mut y_left = zero;
        _mm512_storeu_pd(mp, zero);
        _mm512_storeu_pd(xp, zero);
        _mm512_storeu_pd(yp, zero);
        for j in 1..=h_max {
            let hc = _mm512_cvtepu8_epi64(_mm_loadl_epi64(hp.add((j - 1) * LD) as *const __m128i));
            let mask = _mm512_cmpeq_epi64_mask(rcv, hc) | rc_n | _mm512_cmpeq_epi64_mask(hc, four);
            let prior = _mm512_mask_blend_pd(mask, d3, one_m);
            let (pm, px, py) = (mp.add(j * LD), xp.add(j * LD), yp.add(j * LD));
            let m_up = _mm512_loadu_pd(pm);
            let x_up = _mm512_loadu_pd(px);
            let y_up = _mm512_loadu_pd(py);
            let m = _mm512_mul_pd(_mm512_fmadd_pd(y2, pgapm, _mm512_fmadd_pd(x2, pgapm, _mm512_mul_pd(m2, pmm))), prior);
            let x = _mm512_fmadd_pd(x_up, pxx, _mm512_mul_pd(m_up, pmx));
            let y = _mm512_fmadd_pd(y_left, pyy, _mm512_mul_pd(m_left, pmy));
            _mm512_storeu_pd(pm, m);
            _mm512_storeu_pd(px, x);
            _mm512_storeu_pd(py, y);
            m2 = m_up;
            x2 = x_up;
            y2 = y_up;
            m_left = m;
            y_left = y;
        }
        for (l, p) in lanes.iter().enumerate().take(n) {
            if p.read.bases.len() == i {
                let (mut sm, mut sx) = (0f64, 0f64);
                for j in 1..=p.hap.len() {
                    sm += ws.m[j][l];
                    sx += ws.x[j][l];
                }
                result[l] = sm + sx;
            }
        }
    }
    result
}

#[derive(Default)]
struct TileWs {
    m: Vec<V>,
    x: Vec<V>,
    y: Vec<V>,
    hap: Vec<[u8; L]>,
}

/// Same contract and bit-identical results as [`compute`], but 16 pairs at a time in SIMD lanes
/// (pairs are sorted by size so that a tile wastes few cells).  Single-threaded.
pub fn compute_batched(groups: &[HmmGroup], flavor: Flavor) -> Vec<Vec<f64>> {
    let c = ctx();
    let _g = FtzGuard::on();
    let mut out: Vec<Vec<f64>> = groups.iter().map(|g| vec![0f64; g.reads.len() * g.haps.len()]).collect();
    let mut pairs: Vec<(usize, PairRef)> = Vec::new();
    for (gi, g) in groups.iter().enumerate() {
        let nh = g.haps.len();
        for (ri, r) in g.reads.iter().enumerate() {
            if r.bases.is_empty() {
                continue; // result 0.0
            }
            for (hi, h) in g.haps.iter().enumerate() {
                pairs.push((gi, PairRef { read: r, hap: h, out: ri * nh + hi }));
            }
        }
    }
    pairs.sort_by_key(|(_, p)| (std::cmp::Reverse(p.read.bases.len()), std::cmp::Reverse(p.hap.len())));
    let mut ws = TileWs::default();
    let mut w = Workspace::default();
    let mut lanes: Vec<PairRef> = Vec::with_capacity(L);
    #[cfg(target_arch = "x86_64")]
    let avx512 = std::is_x86_feature_detected!("avx512f") && std::env::var_os("VARIANT_RS_NO_AVX512").is_none();
    #[cfg(target_arch = "x86_64")]
    let mut redo: Vec<(usize, PairRef)> = Vec::new();
    for tile in pairs.chunks(L) {
        lanes.clear();
        lanes.extend(tile.iter().map(|(_, p)| *p));
        let res = match flavor {
            #[cfg(target_arch = "x86_64")]
            Flavor::Avx512 if avx512 => unsafe { tile_float_avx512(&c.f, &lanes, &mut ws) },
            Flavor::Avx512 => tile_float::<true>(&c.f, &lanes, &mut ws),
            Flavor::Avx2 => tile_float::<false>(&c.f, &lanes, &mut ws),
        };
        for (l, (gi, p)) in tile.iter().enumerate() {
            let fl = res[l];
            #[cfg(target_arch = "x86_64")]
            if fl < MIN_ACCEPTED && avx512 && flavor == Flavor::Avx512 {
                redo.push((*gi, *p));
                continue;
            }
            let v = if fl < MIN_ACCEPTED {
                w.hap.clear();
                w.hap.extend(p.hap.iter().map(|&b| conv(b)));
                let hap = std::mem::take(&mut w.hap);
                let dv = match flavor {
                    Flavor::Avx512 => full_prob::<f64, true>(&c.d, p.read, &hap, &mut w.d),
                    Flavor::Avx2 => full_prob::<f64, false>(&c.d, p.read, &hap, &mut w.d),
                };
                w.hap = hap;
                dv.log10() - c.log10_init_d
            } else {
                (fl.log10() - c.log10_init_f) as f64
            };
            out[*gi][p.out] = v;
        }
    }
    // double-precision redo of the pairs whose float result was < 1e-28, 8 at a time
    #[cfg(target_arch = "x86_64")]
    if !redo.is_empty() {
        let mut wd = TileWsD::default();
        for tile in redo.chunks(LD) {
            lanes.clear();
            lanes.extend(tile.iter().map(|(_, p)| *p));
            let res = unsafe { tile_double_avx512(&c.d, &lanes, &mut wd) };
            for (l, (gi, p)) in tile.iter().enumerate() {
                out[*gi][p.out] = res[l].log10() - c.log10_init_d;
            }
        }
    }
    out
}

/// Debug: the float (AVX-512 flavor) pass for one pair, returning every DP cell as
/// (M, X, Y) interleaved over (R+1) x (H+1), and the raw float sum.  Same code as `full_prob`.
pub fn debug_cells_f32(read: &HmmRead, hap: &[u8]) -> (Vec<f32>, f32) {
    let c = ctx();
    let _g = FtzGuard::on();
    let hap: Vec<u8> = hap.iter().map(|&b| conv(b)).collect();
    let mut w = Rows::<f32>::default();
    let (r, h) = (read.bases.len(), hap.len());
    let mut cells = vec![0f32; 3 * (r + 1) * (h + 1)];
    let raw = full_prob_cells::<f32>(&c.f, read, &hap, &mut w, &mut cells);
    (cells, raw)
}

fn full_prob_cells<T: Num>(t: &Tables<T>, read: &HmmRead, hap: &[u8], w: &mut Rows<T>, cells: &mut [T]) -> T {
    let r_len = read.bases.len();
    let h = hap.len();
    let init_y = t.init / T::fi(h as i32);
    w.m = vec![T::ZERO; h + 1];
    w.x = vec![T::ZERO; h + 1];
    w.y = vec![init_y; h + 1];
    let (mut sum_m, mut sum_x) = (T::ZERO, T::ZERO);
    for i in 1..=r_len {
        let qi = (read.ins[i - 1] & 127) as usize;
        let qd = (read.del[i - 1] & 127) as usize;
        let qc = (read.gcp[i - 1] & 127) as usize;
        let qq = (read.quals[i - 1] & 127) as usize;
        let (mx, mn) = if qi <= qd { (qd, qi) } else { (qi, qd) };
        let p_mm = t.mm[((mx * (mx + 1)) >> 1) + mn];
        let p_gapm = T::ONE - t.ph2pr[qc];
        let (p_mx, p_xx, p_my, p_yy) = (t.ph2pr[qi], t.ph2pr[qc], t.ph2pr[qd], t.ph2pr[qc]);
        let distm = t.ph2pr[qq];
        let one_m = T::ONE - distm;
        let distm3 = distm / T::f(3.0);
        let rc = conv(read.bases[i - 1]);
        let (mut m2, mut x2, mut y2) = (w.m[0], w.x[0], w.y[0]);
        let (mut m_left, mut y_left) = (T::ZERO, T::ZERO);
        w.m[0] = T::ZERO;
        w.x[0] = T::ZERO;
        w.y[0] = T::ZERO;
        for j in 1..=h {
            let hc = hap[j - 1];
            let prior = if rc == hc || rc == 4 || hc == 4 { one_m } else { distm3 };
            let (m_up, x_up, y_up) = (w.m[j], w.x[j], w.y[j]);
            let m = T::fma(y2, p_gapm, T::fma(x2, p_gapm, m2 * p_mm)) * prior;
            let x = T::fma(x_up, p_xx, m_up * p_mx);
            let y = T::fma(y_left, p_yy, m_left * p_my);
            w.m[j] = m;
            w.x[j] = x;
            w.y[j] = y;
            let o = 3 * (i * (h + 1) + j);
            cells[o] = m;
            cells[o + 1] = x;
            cells[o + 2] = y;
            m2 = m_up;
            x2 = x_up;
            y2 = y_up;
            m_left = m;
            y_left = y;
        }
        if i == r_len {
            for j in 1..=h {
                sum_m = sum_m + w.m[j];
                sum_x = sum_x + w.x[j];
            }
        }
    }
    sum_m + sum_x
}
