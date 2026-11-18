//! Bit-exact reproductions of the Java / commons-math numerics that reach the VCF text:
//! `Math.round`, `String.format("%.Nf")`, `java.util.Random.nextGaussian` (with fdlibm
//! `StrictMath.log`), commons-math3 `Erf.erfInv`, `SaddlePointExpansion` and the
//! hypergeometric log-probability.

/// `Math.round(double)` (round half up, i.e. floor(x + 1/2) computed exactly), as i64.
#[inline]
pub fn java_round(x: f64) -> i64 {
    if x.is_nan() {
        return 0;
    }
    let f = x.floor();
    let r = if x - f >= 0.5 { f + 1.0 } else { f };
    if r >= i64::MAX as f64 {
        i64::MAX
    } else if r <= i64::MIN as f64 {
        i64::MIN
    } else {
        r as i64
    }
}

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

/// htsjdk `VCFEncoder.formatVCFDouble`
pub fn format_vcf_double(d: f64) -> String {
    if d < 1.0 {
        if d < 0.01 {
            if d.abs() >= 1e-20 {
                format_sci3(d)
            } else {
                "0.00".to_string()
            }
        } else {
            format_fixed(d, 3)
        }
    } else {
        format_fixed(d, 2)
    }
}

/// `String.format("%.3e")`: Java rounds the shortest decimal digits of the double HALF_UP to 4 significant
/// digits (`-10.125` -> `-1.013e+01`, where `{:.3e}` on the binary value gives `-1.012e+01`).
pub fn format_sci3(d: f64) -> String {
    if !d.is_finite() {
        return format_fixed(d, 3);
    }
    let neg = d.is_sign_negative();
    let s = format!("{:e}", d.abs()); // shortest round-trip digits, e.g. "1.0125e1"
    let (m, e) = s.split_once('e').unwrap();
    let mut exp: i32 = e.parse().unwrap();
    let mut digits: Vec<u8> = m.bytes().filter(|c| c.is_ascii_digit()).map(|c| c - b'0').collect();
    if digits.len() > 4 {
        let up = digits[4] >= 5;
        digits.truncate(4);
        if up {
            let mut i = 3;
            loop {
                if digits[i] == 9 {
                    digits[i] = 0;
                    if i == 0 {
                        digits.insert(0, 1);
                        digits.truncate(4);
                        exp += 1;
                        break;
                    }
                    i -= 1;
                } else {
                    digits[i] += 1;
                    break;
                }
            }
        }
    }
    while digits.len() < 4 {
        digits.push(0);
    }
    format!(
        "{}{}.{}{}{}e{}{:02}",
        if neg { "-" } else { "" },
        digits[0],
        digits[1],
        digits[2],
        digits[3],
        if exp < 0 { '-' } else { '+' },
        exp.abs()
    )
}

#[cfg(test)]
mod sci_tests {
    #[test]
    fn java_sci3() {
        for (x, s) in [(-0.431, "-4.310e-01"), (-1.514, "-1.514e+00"), (-12.345, "-1.235e+01"), (-10.125, "-1.013e+01"), (0.005, "5.000e-03"), (-9.9996, "-1.000e+01"), (1e-5, "1.000e-05")] {
            assert_eq!(super::format_sci3(x), s);
        }
    }
}

/// htsjdk `VCFEncoder.formatQualValue`: "%.2f" with a trailing ".00" removed.
pub fn format_qual(q: f64) -> String {
    let s = format_fixed(q, 2);
    match s.strip_suffix(".00") {
        Some(t) => t.to_string(),
        None => s,
    }
}

// ---------------------------------------------------------------------------------------------
// fdlibm __ieee754_log (StrictMath.log)

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

// ---------------------------------------------------------------------------------------------
// java.util.Random

pub struct JavaRandom {
    seed: u64,
    next_next_gaussian: f64,
    have_next_next_gaussian: bool,
}

const MULT: u64 = 0x5_DEEC_E66D;
const MASK: u64 = (1 << 48) - 1;

impl JavaRandom {
    pub fn new(seed: i64) -> Self {
        JavaRandom { seed: (seed as u64 ^ MULT) & MASK, next_next_gaussian: 0.0, have_next_next_gaussian: false }
    }
    fn next(&mut self, bits: u32) -> i32 {
        self.seed = (self.seed.wrapping_mul(MULT).wrapping_add(0xB)) & MASK;
        (self.seed >> (48 - bits)) as i64 as i32
    }
    /// Full generator state (48-bit seed, nextGaussian cache), to share one GATK
    /// `Utils.getRandomGenerator()` with hc_downsample (see variant-rs/crates/downsample/README.md).
    pub fn state(&self) -> (u64, f64, bool) {
        (self.seed, self.next_next_gaussian, self.have_next_next_gaussian)
    }
    pub fn set_state(&mut self, s: (u64, f64, bool)) {
        self.seed = s.0;
        self.next_next_gaussian = s.1;
        self.have_next_next_gaussian = s.2;
    }
    pub fn next_double(&mut self) -> f64 {
        let a = (self.next(26) as i64) << 27;
        let b = self.next(27) as i64;
        (a + b) as f64 * (1.0f64 / (1u64 << 53) as f64)
    }
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
}

/// `Utils.GATK_RANDOM_SEED`
pub const GATK_RANDOM_SEED: i64 = 47382911;

// ---------------------------------------------------------------------------------------------
// commons-math3

/// `FastMath.log`: the exact port of commons-math3 (`hc_assembly::fastmath`).
#[inline]
pub fn fast_log(x: f64) -> f64 {
    hc_assembly::fastmath::log(x)
}

/// commons-math3 `Erf.erfInv`
pub fn erf_inv(x: f64) -> f64 {
    const C1: [f64; 23] = [
        -3.64441206401782E-21, -1.6850591381820166E-19, 1.28584807152564E-18, 1.1157877678025181E-17,
        -1.333171662854621E-16, 2.0972767875968562E-17, 6.637638134358324E-15, -4.054566272975207E-14,
        -8.151934197605472E-14, 2.6335093153082323E-12, -1.2975133253453532E-11, -5.415412054294628E-11,
        1.0512122733215323E-9, -4.112633980346984E-9, -2.9070369957882005E-8, 4.2347877827932404E-7,
        -1.3654692000834679E-6, -1.3882523362786469E-5, 1.8673420803405714E-4, -7.40702534166267E-4,
        -0.006033670871430149, 0.24015818242558962, 1.6536545626831027,
    ];
    const C2: [f64; 19] = [
        2.2137376921775787E-9, 9.075656193888539E-8, -2.7517406297064545E-7, 1.8239629214389228E-8,
        1.5027403968909828E-6, -4.013867526981546E-6, 2.9234449089955446E-6, 1.2475304481671779E-5,
        -4.7318229009055734E-5, 6.828485145957318E-5, 2.4031110387097894E-5, -3.550375203628475E-4,
        9.532893797373805E-4, -0.0016882755560235047, 0.002491442096107851, -0.003751208507569241,
        0.005370914553590064, 1.0052589676941592, 3.0838856104922208,
    ];
    const C3: [f64; 17] = [
        -2.7109920616438573E-11, -2.555641816996525E-10, 1.5076572693500548E-9, -3.789465440126737E-9,
        7.61570120807834E-9, -1.496002662714924E-8, 2.914795345090108E-8, -6.771199775845234E-8,
        2.2900482228026655E-7, -9.9298272942317E-7, 4.526062597223154E-6, -1.968177810553167E-5,
        7.599527703001776E-5, -2.1503011930044477E-4, -1.3871931833623122E-4, 1.0103004648645344,
        4.849906401408584,
    ];
    let mut w = -fast_log((1.0 - x) * (1.0 + x));
    let p;
    let horner = |c: &[f64], w: f64| {
        let mut p = c[0];
        for &ci in &c[1..] {
            p = ci + p * w;
        }
        p
    };
    if w < 6.25 {
        w -= 3.125;
        p = horner(&C1, w);
    } else if w < 16.0 {
        w = w.sqrt() - 3.25;
        p = horner(&C2, w);
    } else if !w.is_infinite() {
        w = w.sqrt() - 5.0;
        p = horner(&C3, w);
    } else {
        p = f64::INFINITY;
    }
    p * x
}

/// `new NormalDistribution(0, 1).inverseCumulativeProbability(p)`
pub fn normal_inverse_cdf(p: f64) -> f64 {
    let sqrt2 = 2.0f64.sqrt();
    0.0 + (1.0 * sqrt2) * erf_inv(2.0 * p - 1.0)
}

const EXACT_STIRLING_ERRORS: [f64; 31] = [
    0.0, 0.15342640972002736, 0.08106146679532726, 0.05481412105191765, 0.0413406959554093,
    0.03316287351993629, 0.02767792568499834, 0.023746163656297496, 0.020790672103765093,
    0.018488450532673187, 0.016644691189821193, 0.015134973221917378, 0.013876128823070748,
    0.012810465242920227, 0.01189670994589177, 0.011104559758206917, 0.010411265261972096,
    0.009799416126158804, 0.009255462182712733, 0.008768700134139386, 0.00833056343336287,
    0.00793411456431402, 0.007573675487951841, 0.007244554301320383, 0.00694284010720953,
    0.006665247032707682, 0.006408994188004207, 0.006171712263039458, 0.0059513701127588475,
    0.0057462165130101155, 0.005554733551962801,
];

/// `SaddlePointExpansion.getStirlingError` (only reached with integer arguments)
fn stirling_error(z: f64) -> f64 {
    if z < 15.0 {
        let z2 = 2.0 * z;
        if z2.floor() == z2 {
            EXACT_STIRLING_ERRORS[z2 as usize]
        } else {
            panic!("non half-integer Stirling error argument {}", z)
        }
    } else {
        let z2 = z * z;
        (0.08333333333333333
            - (0.002777777777777778 - (7.936507936507937E-4 - (5.952380952380953E-4 - 8.417508417508417E-4 / z2) / z2) / z2)
                / z2)
            / z
    }
}

/// `SaddlePointExpansion.getDeviancePart`
fn deviance_part(x: f64, mu: f64) -> f64 {
    if (x - mu).abs() < 0.1 * (x + mu) {
        let d = x - mu;
        let mut v = d / (x + mu);
        let mut s1 = v * d;
        let mut s = f64::NAN;
        let mut ej = 2.0 * x * v;
        v *= v;
        let mut j: i32 = 1;
        while s1 != s {
            s = s1;
            ej *= v;
            s1 = s + ej / ((j * 2 + 1) as f64);
            j += 1;
        }
        s1
    } else {
        x * fast_log(x / mu) + mu - x
    }
}

/// `SaddlePointExpansion.logBinomialProbability`
fn log_binomial_probability(x: i32, n: i32, p: f64, q: f64) -> f64 {
    if x == 0 {
        if p < 0.1 {
            -deviance_part(n as f64, n as f64 * q) - n as f64 * p
        } else {
            n as f64 * fast_log(q)
        }
    } else if x == n {
        if q < 0.1 {
            -deviance_part(n as f64, n as f64 * p) - n as f64 * q
        } else {
            n as f64 * fast_log(p)
        }
    } else {
        let ret = stirling_error(n as f64)
            - stirling_error(x as f64)
            - stirling_error((n - x) as f64)
            - deviance_part(x as f64, n as f64 * p)
            - deviance_part((n - x) as f64, n as f64 * q);
        let f = (6.283185307179586 * x as f64 * (n - x) as f64) / n as f64;
        -0.5 * fast_log(f) + ret
    }
}

/// commons-math3 `BinomialDistribution(n, p).logProbability(k)` (GATK `MathUtils.logBinomialProbability`).
pub fn binomial_log_probability(n: i32, k: i32, p: f64) -> f64 {
    if n == 0 {
        return if k == 0 { 0.0 } else { f64::NEG_INFINITY };
    }
    if k < 0 || k > n {
        return f64::NEG_INFINITY;
    }
    log_binomial_probability(k, n, p, 1.0 - p)
}

/// `HypergeometricDistribution(populationSize, numberOfSuccesses, sampleSize).logProbability(x)`
pub fn hypergeometric_log_probability(population: i32, successes: i32, sample: i32, x: i32) -> f64 {
    let lo = 0.max(successes - (population - sample));
    let hi = successes.min(sample);
    if x < lo || x > hi {
        return f64::NEG_INFINITY;
    }
    let p = sample as f64 / population as f64;
    let q = (population - sample) as f64 / population as f64;
    let p1 = log_binomial_probability(x, successes, p, q);
    let p2 = log_binomial_probability(sample - x, population - successes, p, q);
    let p3 = log_binomial_probability(sample, population, p, q);
    p1 + p2 - p3
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn fixed() {
        assert_eq!(format_fixed(0.125, 2), "0.13");
        assert_eq!(format_fixed(1.005, 2), "1.01");
        assert_eq!(format_fixed(-0.726, 3), "-0.726");
        assert_eq!(format_fixed(1054.0, 2), "1054.00");
        assert_eq!(format_fixed(0.0004, 3), "0.000");
        assert_eq!(format_fixed(0.0005, 3), "0.001");
        assert_eq!(format_fixed(9.9996, 3), "10.000");
        assert_eq!(format_fixed(-0.0, 3), "-0.000");
        assert_eq!(format_qual(30.0), "30");
        assert_eq!(java_round(0.49999999999999994), 0);
        assert_eq!(java_round(-2.5), -2);
    }
}
