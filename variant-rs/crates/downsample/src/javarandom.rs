//! `java.util.Random` (JDK 17), bit-exact: `next(bits)`, `nextInt(bound)`, `nextDouble()`,
//! `nextGaussian()` (polar method with the cached second value and fdlibm `StrictMath.log`).
//!
//! GATK has one process-wide generator, `Utils.getRandomGenerator()` = `new Random(47382911)`,
//! never reseeded on HaplotypeCaller's path. Both `ReservoirDownsampler.submit` (`nextInt`) and
//! `QualByDepth.fixTooHighQD` (`nextGaussian`) draw from it, so the Rust driver must use ONE
//! instance of [`JavaRandom`] for both (see the crate docs and README).
//!
//! `strict_log`, `format_fixed` and the LCG/`nextGaussian` code are copied verbatim from
//! `hc_genotype::javamath` (that crate is not a dependency, so it is left untouched).

/// `Utils.GATK_RANDOM_SEED`
pub const GATK_RANDOM_SEED: i64 = 47382911;

const MULT: u64 = 0x5_DEEC_E66D;
const MASK: u64 = (1 << 48) - 1;

/// Full state of a `java.util.Random`: the 48-bit LCG seed plus `nextGaussian`'s cached value.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JavaRandomState {
    pub seed: u64,
    pub next_next_gaussian: f64,
    pub have_next_next_gaussian: bool,
}

impl JavaRandomState {
    /// `(seed, nextNextGaussian, haveNextNextGaussian)`, the form used by the proposed
    /// `hc_genotype::GenotypingEngine::{rng_state, set_rng_state}` hook (integration/*.patch).
    pub fn to_tuple(self) -> (u64, f64, bool) {
        (self.seed, self.next_next_gaussian, self.have_next_next_gaussian)
    }
    pub fn from_tuple(t: (u64, f64, bool)) -> Self {
        JavaRandomState { seed: t.0, next_next_gaussian: t.1, have_next_next_gaussian: t.2 }
    }
}

#[derive(Clone, Debug)]
pub struct JavaRandom {
    seed: u64,
    next_next_gaussian: f64,
    have_next_next_gaussian: bool,
    /// number of `next(bits)` LCG steps taken since construction (diagnostics only)
    pub steps: u64,
}

impl JavaRandom {
    /// `new Random(seed)`
    pub fn new(seed: i64) -> Self {
        JavaRandom { seed: (seed as u64 ^ MULT) & MASK, next_next_gaussian: 0.0, have_next_next_gaussian: false, steps: 0 }
    }
    /// GATK's global generator in its initial state.
    pub fn gatk() -> Self {
        Self::new(GATK_RANDOM_SEED)
    }
    pub fn state(&self) -> JavaRandomState {
        JavaRandomState { seed: self.seed, next_next_gaussian: self.next_next_gaussian, have_next_next_gaussian: self.have_next_next_gaussian }
    }
    pub fn from_state(s: JavaRandomState) -> Self {
        JavaRandom { seed: s.seed, next_next_gaussian: s.next_next_gaussian, have_next_next_gaussian: s.have_next_next_gaussian, steps: 0 }
    }
    pub fn set_state(&mut self, s: JavaRandomState) {
        self.seed = s.seed;
        self.next_next_gaussian = s.next_next_gaussian;
        self.have_next_next_gaussian = s.have_next_next_gaussian;
    }
    /// `Random.next(bits)`
    #[inline]
    pub fn next(&mut self, bits: u32) -> i32 {
        self.seed = (self.seed.wrapping_mul(MULT).wrapping_add(0xB)) & MASK;
        self.steps += 1;
        (self.seed >> (48 - bits)) as i64 as i32
    }
    /// `Random.nextInt(int bound)` (legacy algorithm, kept by JDK 17's `java.util.Random`).
    #[inline]
    pub fn next_int(&mut self, bound: i32) -> i32 {
        assert!(bound > 0, "bound must be positive");
        let mut r = self.next(31);
        let m = bound - 1;
        if bound & m == 0 {
            ((bound as i64 * r as i64) >> 31) as i32
        } else {
            let mut u = r;
            loop {
                r = u % bound;
                if u.wrapping_sub(r).wrapping_add(m) < 0 {
                    u = self.next(31);
                } else {
                    return r;
                }
            }
        }
    }
    /// `Random.nextDouble()`
    pub fn next_double(&mut self) -> f64 {
        let a = (self.next(26) as i64) << 27;
        let b = self.next(27) as i64;
        (a + b) as f64 * (1.0f64 / (1u64 << 53) as f64)
    }
    /// `Random.nextGaussian()`
    pub fn next_gaussian(&mut self) -> f64 {
        if self.have_next_next_gaussian {
            self.have_next_next_gaussian = false;
            return self.next_next_gaussian;
        }
        loop {
            let v1 = 2.0 * self.next_double() - 1.0;
            let v2 = 2.0 * self.next_double() - 1.0;
            let s = v1 * v1 + v2 * v2;
            if s < 1.0 && s != 0.0 {
                let multiplier = (-2.0 * strict_log(s) / s).sqrt();
                self.next_next_gaussian = v2 * multiplier;
                self.have_next_next_gaussian = true;
                return v1 * multiplier;
            }
        }
    }
    /// `QualByDepth.fixTooHighQD` for a raw QD that is >= 35: `30 + nextGaussian() * 3`.
    pub fn qd_jitter(&mut self) -> f64 {
        30.0 + self.next_gaussian() * 3.0
    }
}

// ---------------------------------------------------------------------------------------------
// copied from hc_genotype::javamath

/// `String.format(Locale.US, "%.<prec>f", x)` as produced by JDK 17 (used for the QD text).
/// `String.format(Locale.US, "%.<prec>f", x)` as produced by JDK 17: the shortest decimal digits of
/// `x` (FloatingDecimal) rounded HALF_UP at the requested precision.
pub fn format_fixed(x: f64, prec: usize) -> String {
    if x.is_nan() {
        return "NaN".to_string();
    }
    if x.is_infinite() {
        return if x > 0.0 { "Infinity".to_string() } else { "-Infinity".to_string() };
    }
    let neg = x.is_sign_negative();
    let a = x.abs();
    let mut out = String::new();
    if neg {
        out.push('-');
    }
    if a == 0.0 {
        out.push('0');
        if prec > 0 {
            out.push('.');
            for _ in 0..prec {
                out.push('0');
            }
        }
        return out;
    }
    // shortest digits: a = 0.d1 d2 ... dn * 10^dec_exp
    let s = format!("{:e}", a);
    let (mant, exp) = s.split_once('e').unwrap();
    let exp: i32 = exp.parse().unwrap();
    let mut digits: Vec<u8> = mant.bytes().filter(|b| b.is_ascii_digit()).collect();
    let mut dec_exp = exp + 1;
    // FormattedFloatingDecimal.applyPrecision
    let keep = dec_exp + prec as i32;
    let n = digits.len() as i32;
    if keep < 0 {
        digits = vec![b'0'];
        dec_exp = 0;
    } else if keep == 0 {
        if digits[0] >= b'5' {
            digits = vec![b'1'];
            dec_exp += 1;
        } else {
            digits = vec![b'0'];
            dec_exp = 0;
        }
    } else if keep < n {
        let k = keep as usize;
        if digits[k] >= b'5' {
            let mut i = k;
            loop {
                if i == 0 {
                    // carry out
                    digits = vec![b'1'];
                    dec_exp += 1;
                    break;
                }
                i -= 1;
                if digits[i] == b'9' {
                    continue;
                }
                digits[i] += 1;
                digits.truncate(i + 1);
                break;
            }
        } else {
            digits.truncate(k);
        }
    }
    // render digits * 10^(dec_exp - len) in fixed notation with `prec` decimals
    let int_len = dec_exp;
    if int_len <= 0 {
        out.push('0');
    } else {
        for i in 0..int_len as usize {
            out.push(if i < digits.len() { digits[i] as char } else { '0' });
        }
    }
    if prec > 0 {
        out.push('.');
        for p in 0..prec as i32 {
            let idx = int_len + p; // digit index
            if idx < 0 || idx as usize >= digits.len() {
                out.push('0');
            } else {
                out.push(digits[idx as usize] as char);
            }
        }
    }
    out
}


// fdlibm StrictMath.log

#[inline]
fn hi(x: f64) -> i32 {
    (x.to_bits() >> 32) as i32
}
#[inline]
fn lo(x: f64) -> u32 {
    x.to_bits() as u32
}
#[inline]
fn with_hi(x: f64, h: i32) -> f64 {
    f64::from_bits(((h as u32 as u64) << 32) | (x.to_bits() & 0xffff_ffff))
}

pub fn strict_log(mut x: f64) -> f64 {
    const LN2_HI: f64 = 6.93147180369123816490e-01;
    const LN2_LO: f64 = 1.90821492927058770002e-10;
    const TWO54: f64 = 1.80143985094819840000e+16;
    const LG1: f64 = 6.666666666666735130e-01;
    const LG2: f64 = 3.999999999940941908e-01;
    const LG3: f64 = 2.857142874366239149e-01;
    const LG4: f64 = 2.222219843214978396e-01;
    const LG5: f64 = 1.818357216161805012e-01;
    const LG6: f64 = 1.531383769920937332e-01;
    const LG7: f64 = 1.479819860511658591e-01;
    let mut hx = hi(x);
    let lx = lo(x);
    let mut k: i32 = 0;
    if hx < 0x0010_0000 {
        if ((hx & 0x7fff_ffff) as u32 | lx) == 0 {
            return f64::NEG_INFINITY;
        }
        if hx < 0 {
            return f64::NAN;
        }
        k -= 54;
        x *= TWO54;
        hx = hi(x);
    }
    if hx >= 0x7ff0_0000 {
        return x + x;
    }
    k += (hx >> 20) - 1023;
    hx &= 0x000f_ffff;
    let i = (hx + 0x95f64) & 0x10_0000;
    x = with_hi(x, hx | (i ^ 0x3ff0_0000));
    k += i >> 20;
    let f = x - 1.0;
    if (0x000f_ffff & (2 + hx)) < 3 {
        if f == 0.0 {
            if k == 0 {
                return 0.0;
            }
            let dk = k as f64;
            return dk * LN2_HI + dk * LN2_LO;
        }
        let r = f * f * (0.5 - 0.33333333333333333 * f);
        if k == 0 {
            return f - r;
        }
        let dk = k as f64;
        return dk * LN2_HI - ((r - dk * LN2_LO) - f);
    }
    let s = f / (2.0 + f);
    let dk = k as f64;
    let z = s * s;
    let mut i = hx - 0x6147a;
    let w = z * z;
    let j = 0x6b851 - hx;
    let t1 = w * (LG2 + w * (LG4 + w * LG6));
    let t2 = z * (LG1 + w * (LG3 + w * (LG5 + w * LG7)));
    i |= j;
    let r = t2 + t1;
    if i > 0 {
        let hfsq = 0.5 * f * f;
        if k == 0 {
            f - (hfsq - s * (hfsq + r))
        } else {
            dk * LN2_HI - ((hfsq - (s * (hfsq + r) + dk * LN2_LO)) - f)
        }
    } else if k == 0 {
        f - s * (f - r)
    } else {
        dk * LN2_HI - ((s * (f - r) - dk * LN2_LO) - f)
    }
}
