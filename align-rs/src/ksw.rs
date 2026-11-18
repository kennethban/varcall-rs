//! ksw.cpp `ksw_align2` (mate rescue): striped Smith-Waterman (Farrar) with 16 x u8 or 8 x i16 lanes, emulated
//! lane by lane. The quirks that shape the result are kept: the row maximum is taken before the lazy-F pass, E is
//! not recomputed after it, the query end is the first maximal cell in lane order of the saved row, and the second
//! best score comes from the per-row maxima list `b`.

pub const KSW_XBYTE: i32 = 0x10000;
pub const KSW_XSTOP: i32 = 0x20000;
pub const KSW_XSUBO: i32 = 0x40000;
pub const KSW_XSTART: i32 = 0x80000;

#[derive(Clone, Copy, Debug)]
pub struct Kswr {
    pub score: i32,
    pub te: i32,
    pub qe: i32,
    pub score2: i32,
    pub te2: i32,
    pub tb: i32,
    pub qb: i32,
}

const DEFR: Kswr = Kswr { score: 0, te: -1, qe: -1, score2: -1, te2: -1, tb: -1, qb: -1 };

struct Query {
    size: usize, // 1 (u8) or 2 (i16)
    slen: usize,
    shift: u8,
    max: i32,
    /// profile per alphabet letter: slen vectors of `lanes` values (u8 values stored in i16)
    qp: Vec<i16>,
}

fn qinit(size: usize, query: &[u8], m: usize, mat: &[i8; 25]) -> Query {
    let size = if size > 1 { 2 } else { 1 };
    let p = 8 * (3 - size);
    let qlen = query.len();
    let slen = (qlen + p - 1) / p;
    let mut shift: i8 = 127;
    let mut mdiff: i8 = 0;
    for &v in mat.iter().take(m * m) {
        if v < shift {
            shift = v;
        }
        if v > mdiff {
            mdiff = v;
        }
    }
    let max = mdiff as i32;
    let shift_u8 = (256i32 - (shift as u8) as i32) as u8;
    let mut qp = Vec::with_capacity(m * slen * p);
    for a in 0..m {
        let ma = &mat[a * m..a * m + m];
        let nlen = slen * p;
        for i in 0..slen {
            let mut k = i;
            while k < nlen {
                let v = if k >= qlen { 0i32 } else { ma[query[k] as usize] as i32 };
                if size == 1 {
                    qp.push(((v + shift_u8 as i32) as i8 as u8) as i16);
                } else {
                    qp.push(v as i16);
                }
                k += slen;
            }
        }
    }
    Query { size, slen, shift: shift_u8, max, qp }
}

/// ksw_u8: 16 lanes of saturating u8
fn ksw_u8(q: &Query, target: &[u8], o_del: i32, e_del: i32, o_ins: i32, e_ins: i32, xtra: i32) -> Kswr {
    const L: usize = 16;
    let slen = q.slen;
    let mut r = DEFR;
    let minsc = if xtra & KSW_XSUBO != 0 { xtra & 0xffff } else { 0x10000 };
    let endsc = if xtra & KSW_XSTOP != 0 { xtra & 0xffff } else { 0x10000 };
    let (oe_del, e_del, oe_ins, e_ins, shift) = ((o_del + e_del) as u8, e_del as u8, (o_ins + e_ins) as u8, e_ins as u8, q.shift);
    let mut h0 = vec![[0u8; L]; slen];
    let mut h1 = vec![[0u8; L]; slen];
    let mut ee = vec![[0u8; L]; slen];
    let mut hmax = vec![[0u8; L]; slen];
    let mut b: Vec<u64> = Vec::new();
    let (mut te, mut gmax) = (-1i32, 0i32);
    for (i, &tc) in target.iter().enumerate() {
        let s = &q.qp[tc as usize * slen * L..];
        let mut f = [0u8; L];
        let mut mx = [0u8; L];
        // h = H0[slen-1] << 1 byte
        let mut h = [0u8; L];
        for l in 1..L {
            h[l] = h0[slen - 1][l - 1];
        }
        for j in 0..slen {
            let sj = &s[j * L..j * L + L];
            let mut hv = [0u8; L];
            for l in 0..L {
                let mut x = h[l].saturating_add(sj[l] as u8);
                x = x.saturating_sub(shift);
                let e = ee[j][l];
                x = x.max(e).max(f[l]);
                hv[l] = x;
                mx[l] = mx[l].max(x);
                let e2 = e.saturating_sub(e_del).max(x.saturating_sub(oe_del));
                ee[j][l] = e2;
                f[l] = f[l].saturating_sub(e_ins).max(x.saturating_sub(oe_ins));
            }
            h1[j] = hv;
            h = h0[j];
        }
        // lazy F
        'lazy: for _k in 0..16 {
            let mut fs = [0u8; L];
            for l in 1..L {
                fs[l] = f[l - 1];
            }
            f = fs;
            for j in 0..slen {
                let mut all = true;
                for l in 0..L {
                    let hh = h1[j][l].max(f[l]);
                    h1[j][l] = hh;
                    let hs = hh.saturating_sub(oe_ins);
                    f[l] = f[l].saturating_sub(e_ins);
                    if f[l].saturating_sub(hs) != 0 {
                        all = false;
                    }
                }
                if all {
                    break 'lazy;
                }
            }
        }
        let imax = *mx.iter().max().unwrap() as i32;
        let i32i = i as i32;
        if imax >= minsc {
            if b.is_empty() || (b[b.len() - 1] as u32 as i32) + 1 != i32i {
                b.push(((imax as u64) << 32) | i as u64);
            } else if ((b[b.len() - 1] >> 32) as i32) < imax {
                let n = b.len();
                b[n - 1] = ((imax as u64) << 32) | i as u64;
            }
        }
        if imax > gmax {
            gmax = imax;
            te = i32i;
            hmax.clone_from(&h1);
            if gmax + shift as i32 >= 255 || gmax >= endsc {
                break;
            }
        }
        std::mem::swap(&mut h0, &mut h1);
    }
    r.score = if gmax + (shift as i32) < 255 { gmax } else { 255 };
    r.te = te;
    if r.score != 255 {
        let mut max = -1i32;
        let qlen = slen * L;
        for i in 0..qlen {
            let v = hmax[i / L][i % L] as i32;
            let pos = (i / L + (i % L) * slen) as i32;
            if v > max {
                max = v;
                r.qe = pos;
            } else if v == max && pos < r.qe {
                r.qe = pos;
            }
        }
        if !b.is_empty() {
            let w = (r.score + q.max - 1) / q.max;
            let (low, high) = (te - w, te + w);
            for &x in &b {
                let e = x as u32 as i32;
                if (e < low || e > high) && ((x >> 32) as i32) > r.score2 {
                    r.score2 = (x >> 32) as i32;
                    r.te2 = e;
                }
            }
        }
    }
    r
}

/// ksw_i16: 8 lanes of i16 (signed saturating add and max, unsigned saturating subtraction)
fn ksw_i16(q: &Query, target: &[u8], o_del: i32, e_del: i32, o_ins: i32, e_ins: i32, xtra: i32) -> Kswr {
    const L: usize = 8;
    let slen = q.slen;
    let mut r = DEFR;
    let minsc = if xtra & KSW_XSUBO != 0 { xtra & 0xffff } else { 0x10000 };
    let endsc = if xtra & KSW_XSTOP != 0 { xtra & 0xffff } else { 0x10000 };
    let subs_epu16 = |a: i16, b: i16| -> i16 { (a as u16).saturating_sub(b as u16) as i16 };
    let (oe_del, e_del, oe_ins, e_ins) = ((o_del + e_del) as i16, e_del as i16, (o_ins + e_ins) as i16, e_ins as i16);
    let mut h0 = vec![[0i16; L]; slen];
    let mut h1 = vec![[0i16; L]; slen];
    let mut ee = vec![[0i16; L]; slen];
    let mut hmax = vec![[0i16; L]; slen];
    let mut b: Vec<u64> = Vec::new();
    let (mut te, mut gmax) = (-1i32, 0i32);
    for (i, &tc) in target.iter().enumerate() {
        let s = &q.qp[tc as usize * slen * L..];
        let mut f = [0i16; L];
        let mut mx = [0i16; L];
        let mut h = [0i16; L];
        for l in 1..L {
            h[l] = h0[slen - 1][l - 1];
        }
        for j in 0..slen {
            let sj = &s[j * L..j * L + L];
            let mut hv = [0i16; L];
            for l in 0..L {
                let mut x = h[l].saturating_add(sj[l]);
                let e = ee[j][l];
                x = x.max(e).max(f[l]);
                hv[l] = x;
                mx[l] = mx[l].max(x);
                let e2 = subs_epu16(e, e_del).max(subs_epu16(x, oe_del));
                ee[j][l] = e2;
                f[l] = subs_epu16(f[l], e_ins).max(subs_epu16(x, oe_ins));
            }
            h1[j] = hv;
            h = h0[j];
        }
        'lazy: for _k in 0..16 {
            let mut fs = [0i16; L];
            for l in 1..L {
                fs[l] = f[l - 1];
            }
            f = fs;
            for j in 0..slen {
                let mut any_gt = false;
                for l in 0..L {
                    let hh = h1[j][l].max(f[l]);
                    h1[j][l] = hh;
                    let hs = subs_epu16(hh, oe_ins);
                    f[l] = subs_epu16(f[l], e_ins);
                    if f[l] > hs {
                        any_gt = true;
                    }
                }
                if !any_gt {
                    break 'lazy;
                }
            }
        }
        let imax = *mx.iter().max().unwrap() as i32;
        let i32i = i as i32;
        if imax >= minsc {
            if b.is_empty() || (b[b.len() - 1] as u32 as i32) + 1 != i32i {
                b.push(((imax as u64) << 32) | i as u64);
            } else if ((b[b.len() - 1] >> 32) as i32) < imax {
                let n = b.len();
                b[n - 1] = ((imax as u64) << 32) | i as u64;
            }
        }
        if imax > gmax {
            gmax = imax;
            te = i32i;
            hmax.clone_from(&h1);
            if gmax >= endsc {
                break;
            }
        }
        std::mem::swap(&mut h0, &mut h1);
    }
    r.score = gmax;
    r.te = te;
    let mut max = -1i32;
    r.qe = -1;
    let qlen = slen * L;
    for i in 0..qlen {
        let v = hmax[i / L][i % L] as u16 as i32;
        let pos = (i / L + (i % L) * slen) as i32;
        if v > max {
            max = v;
            r.qe = pos;
        } else if v == max && pos < r.qe {
            r.qe = pos;
        }
    }
    if !b.is_empty() {
        let w = (r.score + q.max - 1) / q.max;
        let (low, high) = (te - w, te + w);
        for &x in &b {
            let e = x as u32 as i32;
            if (e < low || e > high) && ((x >> 32) as i32) > r.score2 {
                r.score2 = (x >> 32) as i32;
                r.te2 = e;
            }
        }
    }
    r
}

// ------------------------------------------------------------------------------------------------
// SSE2 versions, instruction for instruction as ksw.cpp (the lane-by-lane versions above are the reference)
// ------------------------------------------------------------------------------------------------
#[cfg(target_arch = "x86_64")]
mod sse {
    use super::*;
    use std::arch::x86_64::*;

    unsafe fn profile_vecs(q: &Query) -> Vec<__m128i> {
        // the profile as 16-byte vectors (u8: bytes; i16: words), straight from the i16 staging values
        let per = if q.size == 1 { 16 } else { 8 };
        let n = q.qp.len() / per;
        let mut v: Vec<__m128i> = Vec::with_capacity(n);
        if q.size == 1 {
            let bytes: Vec<u8> = q.qp.iter().map(|&x| x as u8).collect();
            for k in 0..n {
                v.push(_mm_loadu_si128(bytes.as_ptr().add(k * 16) as *const __m128i));
            }
        } else {
            for k in 0..n {
                v.push(_mm_loadu_si128(q.qp.as_ptr().add(k * 8) as *const __m128i));
            }
        }
        v
    }

    #[target_feature(enable = "sse2,sse4.1")]
    pub unsafe fn ksw_u8(q: &Query, target: &[u8], o_del: i32, e_del: i32, o_ins: i32, e_ins: i32, xtra: i32) -> Kswr {
        let slen = q.slen;
        let qp = profile_vecs(q);
        let mut r = DEFR;
        let minsc = if xtra & KSW_XSUBO != 0 { xtra & 0xffff } else { 0x10000 };
        let endsc = if xtra & KSW_XSTOP != 0 { xtra & 0xffff } else { 0x10000 };
        let zero = _mm_set1_epi32(0);
        let oe_del_v = _mm_set1_epi8((o_del + e_del) as i8);
        let e_del_v = _mm_set1_epi8(e_del as i8);
        let oe_ins_v = _mm_set1_epi8((o_ins + e_ins) as i8);
        let e_ins_v = _mm_set1_epi8(e_ins as i8);
        let shift = _mm_set1_epi8(q.shift as i8);
        let mut h0 = vec![zero; slen];
        let mut h1 = vec![zero; slen];
        let mut ee = vec![zero; slen];
        let mut hmax = vec![zero; slen];
        let mut b: Vec<u64> = Vec::new();
        let (mut te, mut gmax) = (-1i32, 0i32);
        for (i, &tc) in target.iter().enumerate() {
            let sp = &qp[tc as usize * slen..];
            let mut f = zero;
            let mut max = zero;
            let mut h = _mm_slli_si128(h0[slen - 1], 1);
            for j in 0..slen {
                h = _mm_adds_epu8(h, sp[j]);
                h = _mm_subs_epu8(h, shift);
                let mut e = ee[j];
                h = _mm_max_epu8(h, e);
                h = _mm_max_epu8(h, f);
                max = _mm_max_epu8(max, h);
                h1[j] = h;
                e = _mm_subs_epu8(e, e_del_v);
                let t = _mm_subs_epu8(h, oe_del_v);
                e = _mm_max_epu8(e, t);
                ee[j] = e;
                f = _mm_subs_epu8(f, e_ins_v);
                let t = _mm_subs_epu8(h, oe_ins_v);
                f = _mm_max_epu8(f, t);
                h = h0[j];
            }
            'lazy: for _k in 0..16 {
                f = _mm_slli_si128(f, 1);
                for j in 0..slen {
                    let mut h = h1[j];
                    h = _mm_max_epu8(h, f);
                    h1[j] = h;
                    h = _mm_subs_epu8(h, oe_ins_v);
                    f = _mm_subs_epu8(f, e_ins_v);
                    let cmp = _mm_movemask_epi8(_mm_cmpeq_epi8(_mm_subs_epu8(f, h), zero));
                    if cmp == 0xffff {
                        break 'lazy;
                    }
                }
            }
            // __max_16
            let mut xx = max;
            xx = _mm_max_epu8(xx, _mm_srli_si128(xx, 8));
            xx = _mm_max_epu8(xx, _mm_srli_si128(xx, 4));
            xx = _mm_max_epu8(xx, _mm_srli_si128(xx, 2));
            xx = _mm_max_epu8(xx, _mm_srli_si128(xx, 1));
            let imax = _mm_extract_epi16(xx, 0) & 0x00ff;
            let i32i = i as i32;
            if imax >= minsc {
                if b.is_empty() || (b[b.len() - 1] as u32 as i32) + 1 != i32i {
                    b.push(((imax as u64) << 32) | i as u64);
                } else if ((b[b.len() - 1] >> 32) as i32) < imax {
                    let n = b.len();
                    b[n - 1] = ((imax as u64) << 32) | i as u64;
                }
            }
            if imax > gmax {
                gmax = imax;
                te = i32i;
                hmax.copy_from_slice(&h1);
                if gmax + q.shift as i32 >= 255 || gmax >= endsc {
                    break;
                }
            }
            std::mem::swap(&mut h0, &mut h1);
        }
        r.score = if gmax + (q.shift as i32) < 255 { gmax } else { 255 };
        r.te = te;
        if r.score != 255 {
            let mut max = -1i32;
            let mut bytes = vec![0u8; slen * 16];
            for (k, v) in hmax.iter().enumerate() {
                _mm_storeu_si128(bytes[k * 16..].as_mut_ptr() as *mut __m128i, *v);
            }
            for i in 0..slen * 16 {
                let v = bytes[i] as i32;
                let pos = (i / 16 + (i % 16) * slen) as i32;
                if v > max {
                    max = v;
                    r.qe = pos;
                } else if v == max && pos < r.qe {
                    r.qe = pos;
                }
            }
            if !b.is_empty() {
                let w = (r.score + q.max - 1) / q.max;
                let (low, high) = (te - w, te + w);
                for &x in &b {
                    let e = x as u32 as i32;
                    if (e < low || e > high) && ((x >> 32) as i32) > r.score2 {
                        r.score2 = (x >> 32) as i32;
                        r.te2 = e;
                    }
                }
            }
        }
        r
    }

    #[target_feature(enable = "sse2,sse4.1")]
    pub unsafe fn ksw_i16(q: &Query, target: &[u8], o_del: i32, e_del: i32, o_ins: i32, e_ins: i32, xtra: i32) -> Kswr {
        let slen = q.slen;
        let qp = profile_vecs(q);
        let mut r = DEFR;
        let minsc = if xtra & KSW_XSUBO != 0 { xtra & 0xffff } else { 0x10000 };
        let endsc = if xtra & KSW_XSTOP != 0 { xtra & 0xffff } else { 0x10000 };
        let zero = _mm_set1_epi32(0);
        let oe_del_v = _mm_set1_epi16((o_del + e_del) as i16);
        let e_del_v = _mm_set1_epi16(e_del as i16);
        let oe_ins_v = _mm_set1_epi16((o_ins + e_ins) as i16);
        let e_ins_v = _mm_set1_epi16(e_ins as i16);
        let mut h0 = vec![zero; slen];
        let mut h1 = vec![zero; slen];
        let mut ee = vec![zero; slen];
        let mut hmax = vec![zero; slen];
        let mut b: Vec<u64> = Vec::new();
        let (mut te, mut gmax) = (-1i32, 0i32);
        for (i, &tc) in target.iter().enumerate() {
            let sp = &qp[tc as usize * slen..];
            let mut f = zero;
            let mut max = zero;
            let mut h = _mm_slli_si128(h0[slen - 1], 2);
            for j in 0..slen {
                h = _mm_adds_epi16(h, sp[j]);
                let mut e = ee[j];
                h = _mm_max_epi16(h, e);
                h = _mm_max_epi16(h, f);
                max = _mm_max_epi16(max, h);
                h1[j] = h;
                e = _mm_subs_epu16(e, e_del_v);
                let t = _mm_subs_epu16(h, oe_del_v);
                e = _mm_max_epi16(e, t);
                ee[j] = e;
                f = _mm_subs_epu16(f, e_ins_v);
                let t = _mm_subs_epu16(h, oe_ins_v);
                f = _mm_max_epi16(f, t);
                h = h0[j];
            }
            'lazy: for _k in 0..16 {
                f = _mm_slli_si128(f, 2);
                for j in 0..slen {
                    let mut h = h1[j];
                    h = _mm_max_epi16(h, f);
                    h1[j] = h;
                    h = _mm_subs_epu16(h, oe_ins_v);
                    f = _mm_subs_epu16(f, e_ins_v);
                    if _mm_movemask_epi8(_mm_cmpgt_epi16(f, h)) == 0 {
                        break 'lazy;
                    }
                }
            }
            let mut xx = max;
            xx = _mm_max_epi16(xx, _mm_srli_si128(xx, 8));
            xx = _mm_max_epi16(xx, _mm_srli_si128(xx, 4));
            xx = _mm_max_epi16(xx, _mm_srli_si128(xx, 2));
            let imax = _mm_extract_epi16(xx, 0) as i16 as i32 & 0xffff;
            let imax = imax as u16 as i32;
            let i32i = i as i32;
            if imax >= minsc {
                if b.is_empty() || (b[b.len() - 1] as u32 as i32) + 1 != i32i {
                    b.push(((imax as u64) << 32) | i as u64);
                } else if ((b[b.len() - 1] >> 32) as i32) < imax {
                    let n = b.len();
                    b[n - 1] = ((imax as u64) << 32) | i as u64;
                }
            }
            if imax > gmax {
                gmax = imax;
                te = i32i;
                hmax.copy_from_slice(&h1);
                if gmax >= endsc {
                    break;
                }
            }
            std::mem::swap(&mut h0, &mut h1);
        }
        r.score = gmax;
        r.te = te;
        let mut max = -1i32;
        r.qe = -1;
        let mut words = vec![0u16; slen * 8];
        for (k, v) in hmax.iter().enumerate() {
            _mm_storeu_si128(words[k * 8..].as_mut_ptr() as *mut __m128i, *v);
        }
        for i in 0..slen * 8 {
            let v = words[i] as i32;
            let pos = (i / 8 + (i % 8) * slen) as i32;
            if v > max {
                max = v;
                r.qe = pos;
            } else if v == max && pos < r.qe {
                r.qe = pos;
            }
        }
        if !b.is_empty() {
            let w = (r.score + q.max - 1) / q.max;
            let (low, high) = (te - w, te + w);
            for &x in &b {
                let e = x as u32 as i32;
                if (e < low || e > high) && ((x >> 32) as i32) > r.score2 {
                    r.score2 = (x >> 32) as i32;
                    r.te2 = e;
                }
            }
        }
        r
    }
}

/// ksw_align2 with m = 5
pub fn ksw_align2(query: &[u8], target: &[u8], mat: &[i8; 25], o_del: i32, e_del: i32, o_ins: i32, e_ins: i32, xtra: i32) -> Kswr {
    // an empty query (a zero-length read rescued as a mate): ksw.cpp runs its loops with no query segment and
    // returns g_defr (score 0, no end), below any rescue threshold; the striped code here would index H0[-1]
    if query.is_empty() {
        return DEFR;
    }
    let size = if xtra & KSW_XBYTE != 0 { 1 } else { 2 };
    let q = qinit(size, query, 5, mat);
    static EMULATE: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    let emulate = *EMULATE.get_or_init(|| std::env::var("ALIGN_RS_KSW_EMULATE").is_ok());
    let func = |q: &Query, t: &[u8], x: i32| {
        #[cfg(target_arch = "x86_64")]
        if !emulate {
            return unsafe { if q.size == 2 { sse::ksw_i16(q, t, o_del, e_del, o_ins, e_ins, x) } else { sse::ksw_u8(q, t, o_del, e_del, o_ins, e_ins, x) } };
        }
        let _ = emulate;
        if q.size == 2 { ksw_i16(q, t, o_del, e_del, o_ins, e_ins, x) } else { ksw_u8(q, t, o_del, e_del, o_ins, e_ins, x) }
    };
    let mut r = func(&q, target, xtra);
    if (xtra & KSW_XSTART) == 0 || ((xtra & KSW_XSUBO) != 0 && r.score < (xtra & 0xffff)) {
        return r;
    }
    let ql = (r.qe + 1) as usize;
    let tl = (r.te + 1) as usize;
    let mut qr: Vec<u8> = query[..ql].to_vec();
    qr.reverse();
    let mut tr: Vec<u8> = target.to_vec();
    tr[..tl].reverse();
    let q2 = qinit(size, &qr, 5, mat);
    let rr = func(&q2, &tr, KSW_XSTOP | r.score);
    if r.score == rr.score {
        r.tb = r.te - rr.te;
        r.qb = r.qe - rr.qe;
    }
    r
}
