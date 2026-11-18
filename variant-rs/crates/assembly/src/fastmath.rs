//! commons-math 3.5 `FastMath.log` / `log1p` and `Gamma.logGamma`, ported operation by operation
//! (table-driven; they differ from libm in the last bit for some inputs). `LN_MANT` was dumped from
//! `FastMathLiteralArrays.loadLnMant()` by `somatic/javacheck/LnTables.java`.

use std::sync::OnceLock;

const LN_2_A: f64 = 0.693147063255310059;
const LN_2_B: f64 = 1.17304635250823482e-7;
const HEX_40000000: f64 = 1073741824.0;
const TWO_POWER_52: f64 = 4503599627370496.0;

const LN_QUICK_COEF: [[f64; 2]; 9] = [
    [1.0, 5.669184079525E-24],
    [-0.25, -0.25],
    [0.3333333134651184, 1.986821492305628E-8],
    [-0.25, -6.663542893624021E-14],
    [0.19999998807907104, 1.1921056801463227E-8],
    [-0.1666666567325592, -7.800414592973399E-9],
    [0.1428571343421936, 5.650007086920087E-9],
    [-0.12502530217170715, -7.44321345601866E-11],
    [0.11113807559013367, 9.219544613762692E-9],
];

const LN_HI_PREC_COEF: [[f64; 2]; 6] = [
    [1.0, -6.032174644509064E-23],
    [-0.25, -0.25],
    [0.3333333134651184, 1.9868161777724352E-8],
    [-0.2499999701976776, -2.957007209750105E-8],
    [0.19999954104423523, 1.5830993332061267E-10],
    [-0.16624879837036133, -2.6033824355191673E-8],
];

fn ln_mant() -> &'static [[f64; 2]] {
    static T: OnceLock<Vec<[f64; 2]>> = OnceLock::new();
    T.get_or_init(|| {
        let b = include_bytes!("fastmath_lnmant.bin");
        b.chunks_exact(16)
            .map(|c| [f64::from_le_bytes(c[0..8].try_into().unwrap()), f64::from_le_bytes(c[8..16].try_into().unwrap())])
            .collect()
    })
}

/// `FastMath.log(x)`
#[inline]
pub fn log(x: f64) -> f64 {
    log_impl(x, None)
}

fn log_impl(x: f64, hi_prec: Option<&mut [f64; 2]>) -> f64 {
    if x == 0.0 {
        return f64::NEG_INFINITY;
    }
    let mut bits = x.to_bits() as i64;
    if ((bits as u64 & 0x8000000000000000) != 0 || x.is_nan()) && x != 0.0 {
        if let Some(h) = hi_prec {
            h[0] = f64::NAN;
        }
        return f64::NAN;
    }
    if x == f64::INFINITY {
        if let Some(h) = hi_prec {
            h[0] = f64::INFINITY;
        }
        return f64::INFINITY;
    }
    let mut exp = (bits >> 52) as i32 - 1023;
    if (bits & 0x7ff0000000000000) == 0 {
        // subnormal (x == 0 handled above)
        bits <<= 1;
        while (bits & 0x0010000000000000) == 0 {
            exp -= 1;
            bits <<= 1;
        }
    }
    if (exp == -1 || exp == 0) && x < 1.01 && x > 0.99 && hi_prec.is_none() {
        let mut xa = x - 1.0;
        let mut tmp = xa * HEX_40000000;
        let mut aa = xa + tmp - tmp;
        let mut ab = xa - aa;
        xa = aa;
        let xb = ab;
        let last = LN_QUICK_COEF[LN_QUICK_COEF.len() - 1];
        let mut ya = last[0];
        let mut yb = last[1];
        for i in (0..LN_QUICK_COEF.len() - 1).rev() {
            aa = ya * xa;
            ab = ya * xb + yb * xa + yb * xb;
            tmp = aa * HEX_40000000;
            ya = aa + tmp - tmp;
            yb = aa - ya + ab;
            let c = LN_QUICK_COEF[i];
            aa = ya + c[0];
            ab = yb + c[1];
            tmp = aa * HEX_40000000;
            ya = aa + tmp - tmp;
            yb = aa - ya + ab;
        }
        aa = ya * xa;
        ab = ya * xb + yb * xa + yb * xb;
        tmp = aa * HEX_40000000;
        ya = aa + tmp - tmp;
        yb = aa - ya + ab;
        return ya + yb;
    }
    let lnm = ln_mant()[((bits & 0x000ffc0000000000) >> 42) as usize];
    let epsilon = (bits & 0x3ffffffffff) as f64 / (TWO_POWER_52 + (bits & 0x000ffc0000000000) as f64);
    let lnza;
    let mut lnzb = 0.0;
    if hi_prec.is_some() {
        let mut tmp = epsilon * HEX_40000000;
        let mut aa = epsilon + tmp - tmp;
        let mut ab = epsilon - aa;
        let xa = aa;
        let mut xb = ab;
        let numer = (bits & 0x3ffffffffff) as f64;
        let denom = TWO_POWER_52 + (bits & 0x000ffc0000000000) as f64;
        aa = numer - xa * denom - xb * denom;
        xb += aa / denom;
        let last = LN_HI_PREC_COEF[LN_HI_PREC_COEF.len() - 1];
        let mut ya = last[0];
        let mut yb = last[1];
        for i in (0..LN_HI_PREC_COEF.len() - 1).rev() {
            aa = ya * xa;
            ab = ya * xb + yb * xa + yb * xb;
            tmp = aa * HEX_40000000;
            ya = aa + tmp - tmp;
            yb = aa - ya + ab;
            let c = LN_HI_PREC_COEF[i];
            aa = ya + c[0];
            ab = yb + c[1];
            tmp = aa * HEX_40000000;
            ya = aa + tmp - tmp;
            yb = aa - ya + ab;
        }
        aa = ya * xa;
        ab = ya * xb + yb * xa + yb * xb;
        lnza = aa + ab;
        lnzb = -(lnza - aa - ab);
    } else {
        let mut z = -0.16624882440418567;
        z = z * epsilon + 0.19999954120254515;
        z = z * epsilon + -0.2499999997677497;
        z = z * epsilon + 0.3333333333332802;
        z = z * epsilon + -0.5;
        z = z * epsilon + 1.0;
        z *= epsilon;
        lnza = z;
    }
    let e = exp as f64;
    let mut a = LN_2_A * e;
    let mut b = 0.0;
    let mut c = a + lnm[0];
    let mut d = -(c - a - lnm[0]);
    a = c;
    b += d;
    c = a + lnza;
    d = -(c - a - lnza);
    a = c;
    b += d;
    c = a + LN_2_B * e;
    d = -(c - a - LN_2_B * e);
    a = c;
    b += d;
    c = a + lnm[1];
    d = -(c - a - lnm[1]);
    a = c;
    b += d;
    c = a + lnzb;
    d = -(c - a - lnzb);
    a = c;
    b += d;
    if let Some(h) = hi_prec {
        h[0] = a;
        h[1] = b;
    }
    a + b
}

/// `FastMath.log1p(x)`
pub fn log1p(x: f64) -> f64 {
    if x == -1.0 {
        return f64::NEG_INFINITY;
    }
    if x == f64::INFINITY {
        return f64::INFINITY;
    }
    if x > 1e-6 || x < -1e-6 {
        let xpa = 1.0 + x;
        let xpb = -(xpa - 1.0 - x);
        let mut hi = [0.0f64; 2];
        let lores = log_impl(xpa, Some(&mut hi));
        if lores.is_infinite() {
            return lores;
        }
        let fx1 = xpb / xpa;
        let epsilon = 0.5 * fx1 + 1.0;
        epsilon * fx1 + hi[1] + hi[0]
    } else {
        let y = (x * (1.0 / 3.0) - 0.5) * x + 1.0;
        y * x
    }
}

const LANCZOS_G: f64 = 607.0 / 128.0;
const LANCZOS: [f64; 15] = [
    0.99999999999999709182,
    57.156235665862923517,
    -59.597960355475491248,
    14.136097974741747174,
    -0.49191381609762019978,
    0.33994649984811888699e-4,
    0.46523628927048575665e-4,
    -0.98374475304879564677e-4,
    0.15808870322491248884e-3,
    -0.21026444172410488319e-3,
    0.21743961811521264320e-3,
    -0.16431810653676389022e-3,
    0.84418223983852743293e-4,
    -0.26190838401581408670e-4,
    0.36899182659531622704e-5,
];

fn half_log_2_pi() -> f64 {
    static V: OnceLock<f64> = OnceLock::new();
    // FastMath.PI = 105414357.0 / 33554432.0 + 1.984187159361080883e-9 (== Math.PI)
    *V.get_or_init(|| 0.5 * log(2.0 * (105414357.0 / 33554432.0 + 1.984187159361080883e-9)))
}

/// `Gamma.lanczos(x)`
fn lanczos(x: f64) -> f64 {
    let mut sum = 0.0;
    for i in (1..LANCZOS.len()).rev() {
        sum += LANCZOS[i] / (x + i as f64);
    }
    sum + LANCZOS[0]
}

/// `Gamma.invGamma1pm1(x)`, -0.5 <= x <= 1.5
pub fn inv_gamma1pm1(x: f64) -> f64 {
    const A0: f64 = 0.611609510448141581788E-08;
    const A1: f64 = 0.624730830116465516210E-08;
    const B1: f64 = 0.203610414066806987300E+00;
    const B2: f64 = 0.266205348428949217746E-01;
    const B3: f64 = 0.493944979382446875238E-03;
    const B4: f64 = -0.851419432440314906588E-05;
    const B5: f64 = -0.643045481779353022248E-05;
    const B6: f64 = 0.992641840672773722196E-06;
    const B7: f64 = -0.607761895722825260739E-07;
    const B8: f64 = 0.195755836614639731882E-09;
    const P0: f64 = 0.6116095104481415817861E-08;
    const P1: f64 = 0.6871674113067198736152E-08;
    const P2: f64 = 0.6820161668496170657918E-09;
    const P3: f64 = 0.4686843322948848031080E-10;
    const P4: f64 = 0.1572833027710446286995E-11;
    const P5: f64 = -0.1249441572276366213222E-12;
    const P6: f64 = 0.4343529937408594255178E-14;
    const Q1: f64 = 0.3056961078365221025009E+00;
    const Q2: f64 = 0.5464213086042296536016E-01;
    const Q3: f64 = 0.4956830093825887312020E-02;
    const Q4: f64 = 0.2692369466186361192876E-03;
    const C: f64 = -0.422784335098467139393487909917598E+00;
    const C0: f64 = 0.577215664901532860606512090082402E+00;
    const C1: f64 = -0.655878071520253881077019515145390E+00;
    const C2: f64 = -0.420026350340952355290039348754298E-01;
    const C3: f64 = 0.166538611382291489501700795102105E+00;
    const C4: f64 = -0.421977345555443367482083012891874E-01;
    const C5: f64 = -0.962197152787697356211492167234820E-02;
    const C6: f64 = 0.721894324666309954239501034044657E-02;
    const C7: f64 = -0.116516759185906511211397108401839E-02;
    const C8: f64 = -0.215241674114950972815729963053648E-03;
    const C9: f64 = 0.128050282388116186153198626328164E-03;
    const C10: f64 = -0.201348547807882386556893914210218E-04;
    const C11: f64 = -0.125049348214267065734535947383309E-05;
    const C12: f64 = 0.113302723198169588237412962033074E-05;
    const C13: f64 = -0.205633841697760710345015413002057E-06;
    assert!((-0.5..=1.5).contains(&x), "invGamma1pm1 out of range: {x}");
    let t = if x <= 0.5 { x } else { (x - 0.5) - 0.5 };
    if t < 0.0 {
        let a = A0 + t * A1;
        let mut b = B8;
        b = B7 + t * b;
        b = B6 + t * b;
        b = B5 + t * b;
        b = B4 + t * b;
        b = B3 + t * b;
        b = B2 + t * b;
        b = B1 + t * b;
        b = 1.0 + t * b;
        let mut c = C13 + t * (a / b);
        c = C12 + t * c;
        c = C11 + t * c;
        c = C10 + t * c;
        c = C9 + t * c;
        c = C8 + t * c;
        c = C7 + t * c;
        c = C6 + t * c;
        c = C5 + t * c;
        c = C4 + t * c;
        c = C3 + t * c;
        c = C2 + t * c;
        c = C1 + t * c;
        c = C + t * c;
        if x > 0.5 { t * c / x } else { x * ((c + 0.5) + 0.5) }
    } else {
        let mut p = P6;
        p = P5 + t * p;
        p = P4 + t * p;
        p = P3 + t * p;
        p = P2 + t * p;
        p = P1 + t * p;
        p = P0 + t * p;
        let mut q = Q4;
        q = Q3 + t * q;
        q = Q2 + t * q;
        q = Q1 + t * q;
        q = 1.0 + t * q;
        let mut c = C13 + (p / q) * t;
        c = C12 + t * c;
        c = C11 + t * c;
        c = C10 + t * c;
        c = C9 + t * c;
        c = C8 + t * c;
        c = C7 + t * c;
        c = C6 + t * c;
        c = C5 + t * c;
        c = C4 + t * c;
        c = C3 + t * c;
        c = C2 + t * c;
        c = C1 + t * c;
        c = C0 + t * c;
        if x > 0.5 { (t / x) * ((c - 0.5) - 0.5) } else { x * c }
    }
}

/// `Gamma.logGamma1p(x)`
pub fn log_gamma1p(x: f64) -> f64 {
    -log1p(inv_gamma1pm1(x))
}

/// commons-math 3.5 `Gamma.logGamma(x)`
pub fn log_gamma(x: f64) -> f64 {
    if x.is_nan() || x <= 0.0 {
        f64::NAN
    } else if x < 0.5 {
        log_gamma1p(x) - log(x)
    } else if x <= 2.5 {
        log_gamma1p((x - 0.5) - 0.5)
    } else if x <= 8.0 {
        let n = (x - 1.5).floor() as i32;
        let mut prod = 1.0;
        for i in 1..=n {
            prod *= x - i as f64;
        }
        log_gamma1p(x - (n + 1) as f64) + log(prod)
    } else {
        let sum = lanczos(x);
        let tmp = x + LANCZOS_G + 0.5;
        ((x + 0.5) * log(tmp)) - tmp + half_log_2_pi() + log(sum / x)
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn matches_commons_math() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../somatic/javacheck/fm.txt");
        let text = match std::fs::read_to_string(path) {
            Ok(t) => t,
            Err(_) => return,
        };
        let h = |s: &str| f64::from_bits(u64::from_str_radix(s, 16).unwrap());
        let mut bad = std::collections::BTreeMap::<&str, (u64, u64)>::new();
        for l in text.lines() {
            let f: Vec<&str> = l.split(' ').collect();
            let (x, want) = (h(f[1]), h(f[2]));
            let got = match f[0] {
                "log" => super::log(x),
                "log1p" => super::log1p(x),
                "lgamma" => super::log_gamma(x),
                "digamma" => crate::m2math::digamma(x),
                "jlog" => x.ln(),
                "jexp" => x.exp(),
                _ => unreachable!(),
            };
            let e = bad.entry(f[0]).or_default();
            e.0 += 1;
            if got.to_bits() != want.to_bits() && !(got.is_nan() && want.is_nan()) {
                e.1 += 1;
            }
        }
        eprintln!("{bad:?}");
        for k in ["log", "log1p", "lgamma", "digamma"] {
            assert_eq!(bad[k].1, 0, "{k} differs");
        }
    }
}
