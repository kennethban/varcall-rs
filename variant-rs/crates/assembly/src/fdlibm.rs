//! fdlibm 5.3 `log1p` / `expm1` as JDK 17 runs them: `Math.log1p` / `Math.expm1` are not HotSpot
//! intrinsics, they call `StrictMath`, which in JDK 17 is native fdlibm
//! (openjdk/jdk17u `src/java.base/share/native/libfdlibm/s_log1p.c`, `s_expm1.c`). Line-by-line
//! ports: same operations in the same order, `__HI`/`__LO` word access on the raw bits, C int
//! arithmetic as i32. Re-exported as `jmath::log1p` / `jmath::expm1`; verified against `java` by
//! `variant-rs/tests/jmath/run.sh`.

#[inline(always)]
fn hi(x: f64) -> i32 {
    (x.to_bits() >> 32) as u32 as i32
}

#[inline(always)]
fn lo(x: f64) -> i32 {
    x.to_bits() as u32 as i32
}

/// `__HI(x) = h` (low word kept)
#[inline(always)]
fn with_hi(x: f64, h: i32) -> f64 {
    f64::from_bits((x.to_bits() & 0xffff_ffff) | ((h as u32 as u64) << 32))
}

/// `StrictMath.log1p` (fdlibm `s_log1p.c`)
pub fn log1p(x: f64) -> f64 {
    const LN2_HI: f64 = 6.93147180369123816490e-01; // 3fe62e42 fee00000
    const LN2_LO: f64 = 1.90821492927058770002e-10; // 3dea39ef 35793c76
    const TWO54: f64 = 1.80143985094819840000e+16; // 43500000 00000000
    const LP1: f64 = 6.666666666666735130e-01; // 3FE55555 55555593
    const LP2: f64 = 3.999999999940941908e-01; // 3FD99999 9997FA04
    const LP3: f64 = 2.857142874366239149e-01; // 3FD24924 94229359
    const LP4: f64 = 2.222219843214978396e-01; // 3FCC71C5 1D8E78AF
    const LP5: f64 = 1.818357216161805012e-01; // 3FC74664 96CB03DE
    const LP6: f64 = 1.531383769920937332e-01; // 3FC39A09 D078C69F
    const LP7: f64 = 1.479819860511658591e-01; // 3FC2F112 DF3E5244
    let zero = 0.0f64;

    let (mut f, mut c) = (0.0f64, 0.0f64);
    let mut hu: i32 = 0;
    let hx = hi(x); // high word of x
    let ax = hx & 0x7fffffff;

    let mut k: i32 = 1;
    if hx < 0x3FDA827A {
        // x < 0.41422
        if ax >= 0x3ff00000 {
            // x <= -1.0
            if x == -1.0 && hx == 0xbff00000u32 as i32 {
                return -TWO54 / zero; // log1p(-1) = -inf
            } else {
                return (x - x) / (x - x); // log1p(x<-1) = NaN
            }
        }
        if ax < 0x3e200000 {
            // |x| < 2**-29
            if TWO54 + x > zero && ax < 0x3c900000 {
                // |x| < 2**-54
                return x;
            } else {
                return x - x * x * 0.5;
            }
        }
        if hx > 0 || hx <= 0xbfd2bec3u32 as i32 {
            // -0.2929 < x < 0.41422
            k = 0;
            f = x;
            hu = 1;
        }
    }
    if hx >= 0x7ff00000 {
        return x + x;
    }
    if k != 0 {
        let mut u;
        if hx < 0x43400000 {
            u = 1.0 + x;
            hu = hi(u);
            k = (hu >> 20) - 1023;
            c = if k > 0 { 1.0 - (u - x) } else { x - (u - 1.0) }; // correction term
            c /= u;
        } else {
            u = x;
            hu = hi(u);
            k = (hu >> 20) - 1023;
            c = 0.0;
        }
        hu &= 0x000fffff;
        if hu < 0x6a09e {
            u = with_hi(u, hu | 0x3ff00000); // normalize u
        } else {
            k += 1;
            u = with_hi(u, hu | 0x3fe00000); // normalize u/2
            hu = (0x00100000 - hu) >> 2;
        }
        f = u - 1.0;
    }
    let hfsq = 0.5 * f * f;
    if hu == 0 {
        // |f| < 2**-20
        if f == zero {
            if k == 0 {
                return zero;
            } else {
                c += k as f64 * LN2_LO;
                return k as f64 * LN2_HI + c;
            }
        }
        let r = hfsq * (1.0 - 0.66666666666666666 * f);
        if k == 0 {
            return f - r;
        } else {
            return k as f64 * LN2_HI - ((r - (k as f64 * LN2_LO + c)) - f);
        }
    }
    let s = f / (2.0 + f);
    let z = s * s;
    let r = z * (LP1 + z * (LP2 + z * (LP3 + z * (LP4 + z * (LP5 + z * (LP6 + z * LP7))))));
    if k == 0 {
        f - (hfsq - s * (hfsq + r))
    } else {
        k as f64 * LN2_HI - ((hfsq - (s * (hfsq + r) + (k as f64 * LN2_LO + c))) - f)
    }
}

/// `StrictMath.expm1` (fdlibm `s_expm1.c`)
pub fn expm1(x: f64) -> f64 {
    const ONE: f64 = 1.0;
    const HUGE: f64 = 1.0e+300;
    const TINY: f64 = 1.0e-300;
    const O_THRESHOLD: f64 = 7.09782712893383973096e+02; // 0x40862E42, 0xFEFA39EF
    const LN2_HI: f64 = 6.93147180369123816490e-01; // 0x3fe62e42, 0xfee00000
    const LN2_LO: f64 = 1.90821492927058770002e-10; // 0x3dea39ef, 0x35793c76
    const INVLN2: f64 = 1.44269504088896338700e+00; // 0x3ff71547, 0x652b82fe
    const Q1: f64 = -3.33333333333331316428e-02; // BFA11111 111110F4
    const Q2: f64 = 1.58730158725481460165e-03; // 3F5A01A0 19FE5585
    const Q3: f64 = -7.93650757867487942473e-05; // BF14CE19 9EAADBB7
    const Q4: f64 = 4.00821782732936239552e-06; // 3ED0CFCA 86E65239
    const Q5: f64 = -2.01099218183624371326e-07; // BE8AFDB7 6E09C32D

    let mut x = x;
    let mut c = 0.0f64;
    let k: i32;
    let mut hx = hi(x) as u32; // high word of x
    let xsb = hx & 0x80000000; // sign bit of x
    hx &= 0x7fffffff; // high word of |x|

    // filter out huge and non-finite argument
    if hx >= 0x4043687A {
        // |x| >= 56*ln2
        if hx >= 0x40862E42 {
            // |x| >= 709.78...
            if hx >= 0x7ff00000 {
                if ((hx & 0xfffff) as i32 | lo(x)) != 0 {
                    return x + x; // NaN
                } else {
                    return if xsb == 0 { x } else { -1.0 }; // exp(+-inf) = {inf, -1}
                }
            }
            if x > O_THRESHOLD {
                return HUGE * HUGE; // overflow
            }
        }
        if xsb != 0 {
            // x < -56*ln2, return -1.0 with inexact
            if x + TINY < 0.0 {
                return TINY - ONE; // return -1
            }
        }
    }

    // argument reduction
    if hx > 0x3fd62e42 {
        // |x| > 0.5 ln2
        let (hi_, lo_);
        if hx < 0x3FF0A2B2 {
            // and |x| < 1.5 ln2
            if xsb == 0 {
                hi_ = x - LN2_HI;
                lo_ = LN2_LO;
                k = 1;
            } else {
                hi_ = x + LN2_HI;
                lo_ = -LN2_LO;
                k = -1;
            }
        } else {
            k = (INVLN2 * x + if xsb == 0 { 0.5 } else { -0.5 }) as i32;
            let t = k as f64;
            hi_ = x - t * LN2_HI; // t*ln2_hi is exact here
            lo_ = t * LN2_LO;
        }
        x = hi_ - lo_;
        c = (hi_ - x) - lo_;
    } else if hx < 0x3c900000 {
        // when |x| < 2**-54, return x
        let t = HUGE + x; // return x with inexact flags when x != 0
        return x - (t - (HUGE + x));
    } else {
        k = 0;
    }

    // x is now in primary range
    let hfx = 0.5 * x;
    let hxs = x * hfx;
    let r1 = ONE + hxs * (Q1 + hxs * (Q2 + hxs * (Q3 + hxs * (Q4 + hxs * Q5))));
    let mut t = 3.0 - r1 * hfx;
    let mut e = hxs * ((r1 - t) / (6.0 - x * t));
    if k == 0 {
        return x - (x * e - hxs); // c is 0
    }
    e = x * (e - c) - c;
    e -= hxs;
    if k == -1 {
        return 0.5 * (x - e) - 0.5;
    }
    if k == 1 {
        if x < -0.25 {
            return -2.0 * (e - (x + 0.5));
        } else {
            return ONE + 2.0 * (x - e);
        }
    }
    let mut y;
    if k <= -2 || k > 56 {
        // suffice to return exp(x)-1
        y = ONE - (e - x);
        y = with_hi(y, hi(y).wrapping_add(k.wrapping_shl(20))); // add k to y's exponent
        return y - ONE;
    }
    t = ONE;
    if k < 20 {
        t = with_hi(t, 0x3ff00000 - (0x200000 >> k)); // t = 1-2^-k
        y = t - (e - x);
        y = with_hi(y, hi(y).wrapping_add(k.wrapping_shl(20))); // add k to y's exponent
    } else {
        t = with_hi(t, (0x3ff - k) << 20); // 2^-k
        y = x - (e + t);
        y += ONE;
        y = with_hi(y, hi(y).wrapping_add(k.wrapping_shl(20))); // add k to y's exponent
    }
    y
}
