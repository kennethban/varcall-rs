//! `java.lang.Math.exp` / `Math.log` as computed by HotSpot (JDK 17) on x86_64: the intrinsic stubs
//! `MacroAssembler::fast_exp` / `fast_log` (Intel LIBM), ported instruction by instruction (SSE2 lanes as
//! scalars, integer lane operations on the raw bits). They differ from glibc's exp/log in the last bit
//! for a small fraction of inputs. `log` uses `rcpps`, whose approximation is CPU-specific, through the
//! same instruction. Checked against the JVM by `somatic/javacheck/HmCheck.java`.

use crate::hotspot_tables::*;

#[inline]
fn qw(t: &[u32], byte_off: usize) -> u64 {
    let i = byte_off / 4;
    (t[i] as u64) | ((t[i + 1] as u64) << 32)
}

#[inline]
fn fd(t: &[u32], byte_off: usize) -> f64 {
    f64::from_bits(qw(t, byte_off))
}

#[inline]
fn f(b: u64) -> f64 {
    f64::from_bits(b)
}

/// packed 32-bit add / subtract on the two dwords of a qword
#[inline]
fn paddd(a: u64, b: u64) -> u64 {
    let lo = (a as u32).wrapping_add(b as u32) as u64;
    let hi = ((a >> 32) as u32).wrapping_add((b >> 32) as u32) as u64;
    lo | (hi << 32)
}

#[inline]
fn psubd(a: u64, b: u64) -> u64 {
    let lo = (a as u32).wrapping_sub(b as u32) as u64;
    let hi = ((a >> 32) as u32).wrapping_sub((b >> 32) as u32) as u64;
    lo | (hi << 32)
}

/// `Math.exp(x)` (HotSpot x86_64 intrinsic)
pub fn exp(x: f64) -> f64 {
    let xb = x.to_bits();
    let hi = (xb >> 32) as u32;
    let mut eax: i32 = ((xb >> 48) & 0xffff) as i32 & 32767;
    let mut edx: i32 = 16527 - eax;
    eax -= 15504;
    edx |= eax;
    if edx as u32 >= 0x8000_0000 {
        // L_2TAG_PACKET_0_0_2
        let a = hi & 0x7fff_ffff;
        if a >= 1083179008 {
            return exp_large(x);
        }
        return x + fd(&EXP_ONE_VAL, 0);
    }
    let c0 = fd(&EXP_CV, 0);
    let shifter = fd(&EXP_SHIFTER, 0);
    let c1 = fd(&EXP_CV, 16);
    let c2 = fd(&EXP_CV, 32);
    let t1 = x * c0 + shifter; // mulpd, addpd
    let x7_bits = t1.to_bits();
    let k = t1 - shifter;
    let x2 = c1 * k;
    let p4 = [fd(&EXP_CV, 64), fd(&EXP_CV, 72)];
    let x3 = c2 * k;
    let p5 = [fd(&EXP_CV, 80), fd(&EXP_CV, 88)];
    let mut r = x - x2;
    let mut eax = x7_bits as u32 as i32;
    let ecx = ((eax & 63) << 4) as usize;
    eax >>= 6;
    let mut edx = eax;
    let x7 = ((x7_bits & qw(&EXP_MMASK, 0)).wrapping_add(qw(&EXP_BIAS, 0))) << 46;
    r -= x3;
    let tlo = fd(&EXP_TBL_ADDR, ecx);
    let thi = qw(&EXP_TBL_ADDR, ecx + 8);
    let x4 = [p4[0] * r, p4[1] * r];
    let r2 = r * r;
    let mut x0 = [r * r2, r * r2];
    let x5 = [p5[0] + x4[0], p5[1] + x4[1]];
    x0[0] *= r2;
    let x6 = r2 * fd(&EXP_CV, 48);
    let mut x1 = r + tlo;
    let t = thi | x7;
    let x0 = [x0[0] * x5[0], x0[1] * x5[1]];
    x1 += x0[0];
    let mut res = x0[1] + x1;
    res += x6;
    edx = edx.wrapping_add(894);
    if (edx as u32) <= 1916 {
        return res * f(t) + f(t);
    }
    // L_2TAG_PACKET_1_0_2: result near overflow / underflow
    let n = eax;
    let shift = (-1022i32).wrapping_sub(n);
    let cnt = shift as u32 as u64;
    let x4 = if cnt > 63 { 0 } else { qw(&EXP_ALLONES, 0) << cnt };
    let half = n >> 1;
    let x3 = (((half as u32 & 0xffff) as u64) << 48) << 4;
    let ebias = qw(&EXP_EBIAS, 0);
    let t2 = psubd(t, x3);
    let mut x0 = res * f(t2);
    if shift > 52 {
        // L_2TAG_PACKET_2_0_2
        let s = f(paddd(x3, ebias));
        return (x0 + f(t2)) * s;
    }
    let x4 = x4 & t2;
    let s = f(paddd(x3, ebias));
    let d = f(t2) - f(x4);
    x0 += d;
    if n >= 1023 {
        // L_2TAG_PACKET_3_0_2 (overflow check only sets the error code)
        return (x0 + f(x4)) * s;
    }
    let sign = ((x0.to_bits() >> 48) as u32) & 32768;
    if (shift as u32 | sign) == 0 {
        // L_2TAG_PACKET_4_0_2
        return (x0 + f(x4)) * s;
    }
    let x6 = x0;
    let r = (x0 + f(x4)) * s;
    if ((r.to_bits() >> 48) as u32 & 32752) != 0 {
        return r;
    }
    // L_2TAG_PACKET_5_0_2: subnormal result, summed on the bit patterns
    let x6b = (x6 * s).to_bits();
    let x4b = (f(x4) * s).to_bits();
    let m: u64 = if ((x6b ^ x4b) >> 63) & 1 == 1 { u64::MAX } else { 0 };
    let mut x0b = (x6b << 1) >> 1;
    x0b ^= m;
    x0b = x0b.wrapping_add(m >> 63);
    x0b = x0b.wrapping_add(x4b);
    f(x0b)
}

/// |x| >= 1024, infinities and NaN
fn exp_large(x: f64) -> f64 {
    let xb = x.to_bits();
    let hi = (xb >> 32) as u32;
    let a = hi & 0x7fff_ffff;
    if a >= 2146435072 {
        let lo = xb as u32;
        if a > 2146435072 || lo != 0 {
            return x + x;
        }
        return if hi == 2146435072 { fd(&EXP_INF, 0) } else { fd(&EXP_ZERO, 0) };
    }
    if hi >= 0x8000_0000 {
        let m = fd(&EXP_XMIN, 0);
        m * m
    } else {
        let m = fd(&EXP_XMAX, 0);
        m * m
    }
}

/// `rcpps` on the low qword (two floats)
#[inline]
fn rcpps_lo(q: u64) -> u64 {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        use std::arch::x86_64::*;
        let v = _mm_castsi128_ps(_mm_set_epi64x(0, q as i64));
        let r = _mm_rcp_ps(v);
        _mm_cvtsi128_si64(_mm_castps_si128(r)) as u64
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let lo = 1.0f32 / f32::from_bits(q as u32);
        let hi = 1.0f32 / f32::from_bits((q >> 32) as u32);
        lo.to_bits() as u64 | ((hi.to_bits() as u64) << 32)
    }
}

/// psrlq 27 then psrld 2 then rcpps, of (x | 1.0)
#[inline]
fn log_rcp(xb: u64) -> u64 {
    let q = (xb | 0x3ff0_0000_0000_0000) >> 27;
    let q = (((q as u32) >> 2) as u64) | ((((q >> 32) as u32) >> 2) as u64) << 32;
    rcpps_lo(q)
}

/// `Math.log(x)` (HotSpot x86_64 intrinsic)
pub fn log(x: f64) -> f64 {
    let xb = x.to_bits();
    let w3 = ((xb >> 48) & 0xffff) as i32;
    let eax = w3 - 16;
    if (eax as u32) < 32736 {
        return log_core(log_rcp(xb), (xb << 12) >> 12, eax, 16352);
    }
    // L_2TAG_PACKET_0_0_2
    let eax = eax + 16;
    if eax >= 32768 {
        // negative (L_2TAG_PACKET_2_0_2)
        let edx = xb as u32;
        let ecx = ((xb >> 32) as u32).wrapping_add((xb >> 32) as u32);
        if ecx >= 0xffe0_0000 {
            if ecx > 0xffe0_0000 || edx > 0 {
                return x + x;
            }
            return 0.0 * f64::INFINITY; // L6: NaN
        }
        if (edx | ecx) == 0 {
            return -1.0 / 0.0; // L7
        }
        return 0.0 * f64::INFINITY; // L6
    }
    if eax < 16 {
        // zero or subnormal (L_2TAG_PACKET_3_0_2)
        let x1 = 0.0 + x;
        if x1.to_bits() == 0 {
            return -1.0 / 0.0;
        }
        let x0 = x * f(0x47f0_0000_0000_0000);
        let b = x0.to_bits();
        let w = ((b >> 48) & 0xffff) as i32;
        return log_core(log_rcp(b), (b << 12) >> 12, w, 18416);
    }
    x + x // +inf / NaN (L_2TAG_PACKET_4_0_2)
}

/// L_2TAG_PACKET_1_0_2 onwards: `rcp` = rcpps result (low qword), `mant` = mantissa bits of x
fn log_core(rcp: u64, mant: u64, eax: i32, ecx: i32) -> f64 {
    let mask = 0xffff_e000_0000_0000u64;
    let mut x0 = paddd(rcp, 32768);
    let x1 = mant | 0x77f0_0000_0000_0000;
    let edx = x0 as u32;
    x0 <<= 29;
    let x5 = mask & x1;
    x0 &= mask;
    let mut x1f = f(x1) - f(x5);
    let mut x5f = f(x5) * f(x0);
    let eax = (eax & 32752) - ecx;
    let mut x7 = eax as f64;
    x1f *= f(x0);
    let log2_hi = fd(&LOG_LOG2, 0);
    let c0 = fd(&LOG_COEFF, 0);
    let c1 = fd(&LOG_COEFF, 8);
    x5f -= 1.0;
    let idx = ((edx & 16711680) >> 12) as usize;
    let t_hi = fd(&LOG_L_TBL, idx);
    let t_lo = fd(&LOG_L_TBL, idx + 8);
    let c4 = [fd(&LOG_COEFF, 16), fd(&LOG_COEFF, 24)];
    let r = x1f + x5f;
    let c2 = [fd(&LOG_COEFF, 32), fd(&LOG_COEFF, 40)];
    let x6 = log2_hi * x7;
    x7 *= fd(&LOG_LOG2, 8);
    let x3 = [c0 * r, c1];
    let a = t_hi + x6;
    let x4 = [c4[0] * r, c4[1] * r];
    let x5 = [r * r, r * r];
    let x0lo = a + r;
    let x4 = [x4[0] + c2[0], x4[1] + c2[1]];
    let x3 = [x3[0] * x5[0], x3[1] * x5[1]];
    let x6 = a - x0lo;
    let x4 = [x4[0] * r, x4[1]];
    let mut x1 = r + x6;
    let x5 = [x5[0] * x5[0], x5[1]];
    x7 += t_lo;
    let x4 = [x4[0] + x3[0], x4[1] + x3[1]];
    x1 += x7;
    let x4 = [x4[0] * x5[0], x4[1] * x5[1]];
    x1 += x4[0];
    x1 += x4[1];
    x0lo + x1
}

#[cfg(test)]
mod tests {
    /// Against the JVM (`somatic/javacheck/HmCheck.java` -> hm.txt: "exp|log x_bits result_bits")
    #[test]
    fn matches_hotspot() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../somatic/javacheck/hm.txt");
        let text = match std::fs::read_to_string(path) {
            Ok(t) => t,
            Err(_) => return,
        };
        let (mut n, mut bad, mut libm) = (0, 0, 0);
        for l in text.lines() {
            let t: Vec<&str> = l.split(' ').collect();
            let x = f64::from_bits(u64::from_str_radix(t[1], 16).unwrap());
            let want = u64::from_str_radix(t[2], 16).unwrap();
            let got = if t[0] == "exp" { super::exp(x) } else { super::log(x) };
            n += 1;
            let lm = if t[0] == "exp" { x.exp() } else { x.ln() };
            if lm.to_bits() != want && !(lm.is_nan() && f64::from_bits(want).is_nan()) {
                libm += 1;
            }
            if got.to_bits() != want && !(got.is_nan() && f64::from_bits(want).is_nan()) {
                bad += 1;
                if bad <= 10 {
                    eprintln!("{} {:e}: got {:e} want {:e}", t[0], x, got, f64::from_bits(want));
                }
            }
        }
        eprintln!("hotspot exp/log: {bad} / {n} differ (libm: {libm})");
        assert_eq!(bad, 0);
    }
}
