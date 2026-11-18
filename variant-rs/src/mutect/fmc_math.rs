//! Numerics of FilterMutectCalls: commons-math 3.5 (`FastMath.log` / `log1p` / `exp` exact, via
//! `hc_assembly::fastmath` / `m2math::FastExp`) and the java.lang.Math / formatting helpers GATK uses.
//! java.lang.Math `exp` / `log` / `log10` / `pow` are HotSpot intrinsics (`crate::jmath` ports them bit
//! for bit; glibc differs in the last bit for some inputs); `Math.log1p` / `expm1` are fdlibm (`jmath::log1p` /
//! `expm1`). Every java.lang.Math call below goes through them.

use hc_assembly::fastmath;
use hc_assembly::m2math::FastExp;
use std::sync::OnceLock;

// ------------------------------------------------------------------------------------------------
// FastMath

#[inline]
pub fn flog(x: f64) -> f64 {
    fastmath::log(x)
}

#[inline]
pub fn flog1p(x: f64) -> f64 {
    fastmath::log1p(x)
}

pub fn fexp(x: f64) -> f64 {
    static E: OnceLock<FastExp> = OnceLock::new();
    E.get_or_init(FastExp::load).exp(x)
}

/// `MathUtils.TWO_PI` (FastMath.PI == Math.PI)
const TWO_PI: f64 = 2.0 * std::f64::consts::PI;

// ------------------------------------------------------------------------------------------------
// CombinatoricsUtils

fn binomial_exact(n: i64, k: i64) -> i64 {
    let k = k.min(n - k);
    let mut r: u128 = 1;
    for j in 1..=k {
        r = r * (n - k + j) as u128 / j as u128;
    }
    r as i64
}

fn binomial_double(n: i64, k: i64) -> f64 {
    if n == k || k == 0 {
        return 1.0;
    }
    if k == 1 || k == n - 1 {
        return n as f64;
    }
    if k > n / 2 {
        return binomial_double(n, n - k);
    }
    if n < 67 {
        return binomial_exact(n, k) as f64;
    }
    let mut r = 1.0f64;
    for i in 1..=k {
        r *= (n - k + i) as f64 / i as f64;
    }
    (r + 0.5).floor()
}

/// `CombinatoricsUtils.binomialCoefficientLog`
pub fn binomial_coefficient_log(n: i64, k: i64) -> f64 {
    assert!(n >= k && k >= 0, "binomial coefficient {n} choose {k}");
    if n == k || k == 0 {
        return 0.0;
    }
    if k == 1 || k == n - 1 {
        return flog(n as f64);
    }
    if n < 67 {
        return flog(binomial_exact(n, k) as f64);
    }
    if n < 1030 {
        return flog(binomial_double(n, k));
    }
    if k > n / 2 {
        return binomial_coefficient_log(n, n - k);
    }
    let mut s = 0.0;
    for i in n - k + 1..=n {
        s += flog(i as f64);
    }
    for i in 2..=k {
        s -= flog(i as f64);
    }
    s
}

// ------------------------------------------------------------------------------------------------
// SaddlePointExpansion / BinomialDistribution

const EXACT_STIRLING_ERRORS: [f64; 31] = [
    0.0,
    0.1534264097200273452913848,
    0.0810614667953272582196702,
    0.0548141210519176538961390,
    0.0413406959554092940938221,
    0.03316287351993628748511048,
    0.02767792568499833914878929,
    0.02374616365629749597132920,
    0.02079067210376509311152277,
    0.01848845053267318523077934,
    0.01664469118982119216319487,
    0.01513497322191737887351255,
    0.01387612882307074799874573,
    0.01281046524292022692424986,
    0.01189670994589177009505572,
    0.01110455975820691732662991,
    0.010411265261972096497478567,
    0.009799416126158803298389475,
    0.009255462182712732917728637,
    0.008768700134139385462952823,
    0.008330563433362871256469318,
    0.007934114564314020547248100,
    0.007573675487951840794972024,
    0.007244554301320383179543912,
    0.006942840107209529865664152,
    0.006665247032707682442354394,
    0.006408994188004207068439631,
    0.006171712263039457647532867,
    0.005951370112758847735624416,
    0.005746216513010115682023589,
    0.005554733551962801371038690,
];

fn stirling_error(z: f64) -> f64 {
    if z < 15.0 {
        let z2 = 2.0 * z;
        if z2.floor() == z2 {
            EXACT_STIRLING_ERRORS[z2 as usize]
        } else {
            let half_log_2_pi = 0.5 * flog(TWO_PI);
            fastmath::log_gamma(z + 1.0) - (z + 0.5) * flog(z) + z - half_log_2_pi
        }
    } else {
        let z2 = z * z;
        (0.083333333333333333333
            - (0.00277777777777777777778
                - (0.00079365079365079365079365 - (0.000595238095238095238095238 - 0.0008417508417508417508417508 / z2) / z2)
                    / z2)
                / z2)
            / z
    }
}

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
        x * flog(x / mu) + mu - x
    }
}

fn log_binomial_probability(x: i32, n: i32, p: f64, q: f64) -> f64 {
    if x == 0 {
        if p < 0.1 {
            -deviance_part(n as f64, n as f64 * q) - n as f64 * p
        } else {
            n as f64 * flog(q)
        }
    } else if x == n {
        if q < 0.1 {
            -deviance_part(n as f64, n as f64 * p) - n as f64 * q
        } else {
            n as f64 * flog(p)
        }
    } else {
        let ret = stirling_error(n as f64)
            - stirling_error(x as f64)
            - stirling_error((n - x) as f64)
            - deviance_part(x as f64, n as f64 * p)
            - deviance_part((n - x) as f64, n as f64 * q);
        let f = (TWO_PI * x as f64 * (n - x) as f64) / n as f64;
        -0.5 * flog(f) + ret
    }
}

/// `new BinomialDistribution(null, n, p).logProbability(k)`
pub fn binomial_log_probability(n: i32, k: i32, p: f64) -> f64 {
    // BinomialDistribution: OutOfRangeException if p < 0 || p > 1 (a NaN passes)
    assert!(!(p < 0.0 || p > 1.0), "binomial p {p}");
    if n == 0 {
        return if k == 0 { 0.0 } else { f64::NEG_INFINITY };
    }
    if k < 0 || k > n {
        return f64::NEG_INFINITY;
    }
    log_binomial_probability(k, n, p, 1.0 - p)
}

/// `new BinomialDistribution(null, n, p).probability(k)` (= `MathUtils.binomialProbability(n, k, p)`)
pub fn binomial_probability(n: i32, k: i32, p: f64) -> f64 {
    let l = binomial_log_probability(n, k, p);
    if l == f64::NEG_INFINITY { 0.0 } else { fexp(l) }
}

/// `new BinomialDistribution(null, n, p).cumulativeProbability(x)`
pub fn binomial_cdf(n: i32, p: f64, x: i32) -> f64 {
    if x < 0 {
        0.0
    } else if x >= n {
        1.0
    } else {
        1.0 - regularized_beta(p, x as f64 + 1.0, (n - x) as f64)
    }
}

// ------------------------------------------------------------------------------------------------
// Gamma / Beta

/// `Gamma.gamma(x)` for |x| <= 20 (the only range reached: logBeta with a < 1, b < 10).
fn gamma(x: f64) -> f64 {
    if x == x.round_ties_even() && x <= 0.0 {
        return f64::NAN;
    }
    assert!(x.abs() <= 20.0, "Gamma.gamma({x}) outside the ported range");
    if x >= 1.0 {
        let mut prod = 1.0;
        let mut t = x;
        while t > 2.5 {
            t -= 1.0;
            prod *= t;
        }
        prod / (1.0 + fastmath::inv_gamma1pm1(t - 1.0))
    } else {
        let mut prod = x;
        let mut t = x;
        while t < -0.5 {
            t += 1.0;
            prod *= t;
        }
        1.0 / (prod * (1.0 + fastmath::inv_gamma1pm1(t)))
    }
}

const HALF_LOG_TWO_PI: f64 = 0.9189385332046727;
const DELTA: [f64; 15] = [
    0.833333333333333333333333333333E-01,
    -0.277777777777777777777777752282E-04,
    0.793650793650793650791732130419E-07,
    -0.595238095238095232389839236182E-09,
    0.841750841750832853294451671990E-11,
    -0.191752691751854612334149171243E-12,
    0.641025640510325475730918472625E-14,
    -0.295506514125338232839867823991E-15,
    0.179643716359402238723287696452E-16,
    -0.139228964661627791231203060395E-17,
    0.133802855014020915603275339093E-18,
    -0.154246009867966094273710216533E-19,
    0.197701992980957427278370133333E-20,
    -0.234065664793997056856992426667E-21,
    0.171348014966398575409015466667E-22,
];

fn log_gamma_sum(a: f64, b: f64) -> f64 {
    let x = (a - 1.0) + (b - 1.0);
    if x <= 0.5 {
        fastmath::log_gamma1p(1.0 + x)
    } else if x <= 1.5 {
        fastmath::log_gamma1p(x) + flog1p(x)
    } else {
        fastmath::log_gamma1p(x - 1.0) + flog(x * (1.0 + x))
    }
}

fn log_gamma_minus_log_gamma_sum(a: f64, b: f64) -> f64 {
    let (d, w) = if a <= b { (b + (a - 0.5), delta_minus_delta_sum(a, b)) } else { (a + (b - 0.5), delta_minus_delta_sum(b, a)) };
    let u = d * flog1p(a / b);
    let v = a * (flog(b) - 1.0);
    if u <= v { (w - u) - v } else { (w - v) - u }
}

fn delta_minus_delta_sum(a: f64, b: f64) -> f64 {
    let h = a / b;
    let p = h / (1.0 + h);
    let q = 1.0 / (1.0 + h);
    let q2 = q * q;
    let mut s = [0.0f64; 15];
    s[0] = 1.0;
    for i in 1..s.len() {
        s[i] = 1.0 + (q + q2 * s[i - 1]);
    }
    let sqrt_t = 10.0 / b;
    let t = sqrt_t * sqrt_t;
    let mut w = DELTA[DELTA.len() - 1] * s[s.len() - 1];
    for i in (0..DELTA.len() - 1).rev() {
        w = t * w + DELTA[i] * s[i];
    }
    w * p / b
}

fn sum_delta_minus_delta_sum(p: f64, q: f64) -> f64 {
    let a = p.min(q);
    let b = p.max(q);
    let sqrt_t = 10.0 / a;
    let t = sqrt_t * sqrt_t;
    let mut z = DELTA[DELTA.len() - 1];
    for i in (0..DELTA.len() - 1).rev() {
        z = t * z + DELTA[i];
    }
    z / a + delta_minus_delta_sum(a, b)
}

/// `Beta.logBeta(p, q)`
pub fn log_beta(p: f64, q: f64) -> f64 {
    if p.is_nan() || q.is_nan() || p <= 0.0 || q <= 0.0 {
        return f64::NAN;
    }
    let a = p.min(q);
    let b = p.max(q);
    let lg = fastmath::log_gamma;
    if a >= 10.0 {
        let w = sum_delta_minus_delta_sum(a, b);
        let h = a / b;
        let c = h / (1.0 + h);
        let u = -(a - 0.5) * flog(c);
        let v = b * flog1p(h);
        if u <= v {
            (((-0.5 * flog(b) + HALF_LOG_TWO_PI) + w) - u) - v
        } else {
            (((-0.5 * flog(b) + HALF_LOG_TWO_PI) + w) - v) - u
        }
    } else if a > 2.0 {
        if b > 1000.0 {
            let n = (a - 1.0).floor() as i32;
            let mut prod = 1.0;
            let mut ared = a;
            for _ in 0..n {
                ared -= 1.0;
                prod *= ared / (1.0 + ared / b);
            }
            (flog(prod) - n as f64 * flog(b)) + (lg(ared) + log_gamma_minus_log_gamma_sum(ared, b))
        } else {
            let mut prod1 = 1.0;
            let mut ared = a;
            while ared > 2.0 {
                ared -= 1.0;
                let h = ared / b;
                prod1 *= h / (1.0 + h);
            }
            if b < 10.0 {
                let mut prod2 = 1.0;
                let mut bred = b;
                while bred > 2.0 {
                    bred -= 1.0;
                    prod2 *= bred / (ared + bred);
                }
                flog(prod1) + flog(prod2) + (lg(ared) + (lg(bred) - log_gamma_sum(ared, bred)))
            } else {
                flog(prod1) + lg(ared) + log_gamma_minus_log_gamma_sum(ared, b)
            }
        }
    } else if a >= 1.0 {
        if b > 2.0 {
            if b < 10.0 {
                let mut prod = 1.0;
                let mut bred = b;
                while bred > 2.0 {
                    bred -= 1.0;
                    prod *= bred / (a + bred);
                }
                flog(prod) + (lg(a) + (lg(bred) - log_gamma_sum(a, bred)))
            } else {
                lg(a) + log_gamma_minus_log_gamma_sum(a, b)
            }
        } else {
            lg(a) + lg(b) - log_gamma_sum(a, b)
        }
    } else if b >= 10.0 {
        lg(a) + log_gamma_minus_log_gamma_sum(a, b)
    } else {
        flog(gamma(a) * gamma(b) / gamma(a + b))
    }
}

/// `Precision.equals(x, y, 1)`
fn precision_equals(x: f64, y: f64) -> bool {
    let xi = x.to_bits() as i64;
    let yi = y.to_bits() as i64;
    let eq = if (xi ^ yi) & i64::MIN == 0 {
        (xi - yi).abs() <= 1
    } else {
        let nz = (-0.0f64).to_bits() as i64;
        let (dp, dm) = if xi < yi { (yi, xi.wrapping_sub(nz)) } else { (xi, yi.wrapping_sub(nz)) };
        if dp > 1 { false } else { dm <= 1 - dp }
    };
    eq && !x.is_nan() && !y.is_nan()
}

/// `Precision.equals(x, y, eps)`
fn precision_equals_eps(x: f64, y: f64, eps: f64) -> bool {
    precision_equals(x, y) || (y - x).abs() <= eps
}

/// `Beta.regularizedBeta(x, a, b)` (epsilon 1e-14, unbounded iterations)
pub fn regularized_beta(x: f64, a: f64, b: f64) -> f64 {
    if x.is_nan() || a.is_nan() || b.is_nan() || x < 0.0 || x > 1.0 || a <= 0.0 || b <= 0.0 {
        f64::NAN
    } else if x > (a + 1.0) / (2.0 + b + a) && 1.0 - x <= (b + 1.0) / (2.0 + b + a) {
        1.0 - regularized_beta(1.0 - x, b, a)
    } else {
        let get_b = |n: i32, x: f64| -> f64 {
            if n % 2 == 0 {
                let m = n as f64 / 2.0;
                (m * (b - m) * x) / ((a + (2.0 * m) - 1.0) * (a + (2.0 * m)))
            } else {
                let m = (n as f64 - 1.0) / 2.0;
                -((a + m) * (a + b + m) * x) / ((a + (2.0 * m)) * (a + (2.0 * m) + 1.0))
            }
        };
        // ContinuedFraction.evaluate(x, 1e-14, Integer.MAX_VALUE), getA == 1
        let small = 1e-50;
        let eps = 1e-14;
        let mut h_prev = 1.0;
        if precision_equals_eps(h_prev, 0.0, small) {
            h_prev = small;
        }
        let mut n = 1;
        let mut d_prev = 0.0;
        let mut c_prev = h_prev;
        let mut h_n = h_prev;
        while n < i32::MAX {
            let aa = 1.0;
            let bb = get_b(n, x);
            let mut d_n = aa + bb * d_prev;
            if precision_equals_eps(d_n, 0.0, small) {
                d_n = small;
            }
            let mut c_n = aa + bb / c_prev;
            if precision_equals_eps(c_n, 0.0, small) {
                c_n = small;
            }
            d_n = 1.0 / d_n;
            let delta_n = c_n * d_n;
            h_n = h_prev * delta_n;
            assert!(h_n.is_finite(), "continued fraction diverged at x = {x}");
            if (delta_n - 1.0).abs() < eps {
                break;
            }
            d_prev = d_n;
            c_prev = c_n;
            h_prev = h_n;
            n += 1;
        }
        fexp((a * flog(x)) + (b * flog1p(-x)) - flog(a) - log_beta(a, b)) * 1.0 / h_n
    }
}

/// GATK `BetaBinomialDistribution(alpha, beta, n).logProbability(k)`
pub fn beta_binomial_log_probability(alpha: f64, beta: f64, n: i32, k: i32) -> f64 {
    assert!(alpha > 0.0 && beta > 0.0 && n >= 0 && k >= 0, "beta-binomial({alpha}, {beta}, {n}) at {k}");
    if k > n {
        f64::NEG_INFINITY
    } else {
        binomial_coefficient_log(n as i64, k as i64) + log_beta(k as f64 + alpha, (n - k) as f64 + beta) - log_beta(alpha, beta)
    }
}

// ------------------------------------------------------------------------------------------------
// BrentOptimizer (no convergence checker), maximisation; returns the best point

pub fn brent_max(f: impl Fn(f64) -> f64, lo: f64, hi: f64, start: f64, rel: f64, abs_t: f64, max_eval: usize) -> f64 {
    brent_max_pv(f, lo, hi, start, rel, abs_t, max_eval).0
}

/// BrentOptimizer maximisation returning the `UnivariatePointValuePair` (point, objective value).
pub fn brent_max_pv(f: impl Fn(f64) -> f64, lo: f64, hi: f64, start: f64, rel: f64, abs_t: f64, max_eval: usize) -> (f64, f64) {
    let golden = 0.5 * (3.0 - 5.0f64.sqrt());
    let (mut a, mut b) = if lo < hi { (lo, hi) } else { (hi, lo) };
    let mut evals = 0usize;
    let mut eval = |x: f64| -> f64 {
        evals += 1;
        assert!(evals <= max_eval, "Brent optimizer: too many evaluations");
        -f(x)
    };
    let mut x = start;
    let mut v = x;
    let mut w = x;
    let mut d = 0.0f64;
    let mut e = 0.0f64;
    let mut fx = eval(x);
    let mut fv = fx;
    let mut fw = fx;
    // points as (x, value in maximisation sign)
    let mut previous: Option<(f64, f64)> = None;
    let mut current = (x, -fx);
    let mut best = current;
    let pick = |a: Option<(f64, f64)>, b: Option<(f64, f64)>| -> Option<(f64, f64)> {
        match (a, b) {
            (None, b) => b,
            (a, None) => a,
            (Some(a), Some(b)) => Some(if a.1 >= b.1 { a } else { b }),
        }
    };
    loop {
        let m = 0.5 * (a + b);
        let tol1 = rel * x.abs() + abs_t;
        let tol2 = 2.0 * tol1;
        let stop = (x - m).abs() <= tol2 - 0.5 * (b - a);
        if stop {
            return pick(Some(best), pick(previous, Some(current))).unwrap();
        }
        let mut p = 0.0;
        let mut q = 0.0;
        let mut r;
        let mut u;
        if e.abs() > tol1 {
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if q > 0.0 {
                p = -p;
            } else {
                q = -q;
            }
            r = e;
            e = d;
            if p > q * (a - x) && p < q * (b - x) && p.abs() < (0.5 * q * r).abs() {
                d = p / q;
                u = x + d;
                if u - a < tol2 || b - u < tol2 {
                    d = if x <= m { tol1 } else { -tol1 };
                }
            } else {
                e = if x < m { b - x } else { a - x };
                d = golden * e;
            }
        } else {
            e = if x < m { b - x } else { a - x };
            d = golden * e;
        }
        let _ = (p, q);
        u = if d.abs() < tol1 {
            if d >= 0.0 { x + tol1 } else { x - tol1 }
        } else {
            x + d
        };
        let fu = eval(u);
        previous = Some(current);
        current = (u, -fu);
        best = pick(Some(best), pick(previous, Some(current))).unwrap();
        if fu <= fx {
            if u < x {
                b = x;
            } else {
                a = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        } else {
            if u < x {
                a = u;
            } else {
                b = u;
            }
            if fu <= fw || precision_equals(w, x) {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            } else if fu <= fv || precision_equals(v, x) || precision_equals(v, w) {
                v = u;
                fv = fu;
            }
        }
    }
}

// ------------------------------------------------------------------------------------------------
// java.lang / GATK MathUtils helpers

/// `Math.exp` (HotSpot intrinsic)
#[inline]
pub fn jexp(x: f64) -> f64 {
    crate::jmath::exp(x)
}

/// `Math.log` (HotSpot intrinsic)
#[inline]
pub fn jlog(x: f64) -> f64 {
    crate::jmath::log(x)
}

/// `MathUtils.log10SumLog10`
pub fn log10_sum_log10(v: &[f64]) -> f64 {
    let mut mi = 0;
    for i in 1..v.len() {
        if v[i] > v[mi] {
            mi = i;
        }
    }
    let max = v[mi];
    if max == f64::NEG_INFINITY {
        return max;
    }
    let mut s = 1.0;
    for (i, &c) in v.iter().enumerate() {
        if i == mi || c == f64::NEG_INFINITY {
            continue;
        }
        s += crate::jmath::pow(10.0, c - max);
    }
    assert!(!(s.is_nan() || s == f64::INFINITY), "log10 p: Values must be non-infinite and non-NAN");
    max + if s != 1.0 { crate::jmath::log10(s) } else { 0.0 }
}

/// `MathUtils.normalizeLog10(array, false, true)`: linear-space probabilities
pub fn normalize_log10_to_linear(v: &mut [f64]) {
    let s = log10_sum_log10(v);
    for x in v.iter_mut() {
        *x -= s;
    }
    for x in v.iter_mut() {
        *x = crate::jmath::pow(10.0, *x);
    }
}

/// `NaturalLogUtils.log1mexp` (java.lang.Math)
pub fn log1mexp(a: f64) -> f64 {
    if a > 0.0 {
        return f64::NAN;
    }
    if a == 0.0 {
        return f64::NEG_INFINITY;
    }
    if a < jlog(0.5) { crate::jmath::log1p(-jexp(a)) } else { jlog(-crate::jmath::expm1(a)) }
}

/// `Math.ulp(double)`
fn java_ulp(x: f64) -> f64 {
    let x = x.abs();
    if x.is_nan() {
        return f64::NAN;
    }
    if x.is_infinite() {
        return f64::INFINITY;
    }
    if x == f64::MAX {
        return 2f64.powi(971);
    }
    let next = f64::from_bits(x.to_bits() + 1);
    next - x
}

/// `MathUtils.roundToNDecimalPlaces(in, n)`
pub fn round_to_n_decimal_places(x: f64, n: i32) -> f64 {
    let mult = crate::jmath::pow(10.0, n as f64); // Math.pow(10, n)
    hc_genotype::javamath::java_round((x + java_ulp(x)) * mult) as f64 / mult
}

/// `QualityUtils.errorProbToQual(p)` (max 93, min 1)
pub fn error_prob_to_qual(p: f64) -> i32 {
    assert!((0.0..=1.0).contains(&p), "errorRate must be good probability but got {p}");
    let d = hc_genotype::javamath::java_round(-10.0 * crate::jmath::log10(p));
    (d.clamp(i32::MIN as i64, i32::MAX as i64) as i32).min(93).max(1)
}

/// Java `Double.toString(double)` (shortest round-trip digits in Java's layout).
pub fn java_double_to_string(x: f64) -> String {
    if x.is_nan() {
        return "NaN".into();
    }
    if x.is_infinite() {
        return if x > 0.0 { "Infinity".into() } else { "-Infinity".into() };
    }
    if x == 0.0 {
        return if x.is_sign_negative() { "-0.0".into() } else { "0.0".into() };
    }
    // shortest digits and decimal exponent from Rust's `{:e}` (round-trip shortest)
    let s = format!("{:e}", x.abs());
    let (mant, exp) = s.split_once('e').unwrap();
    let exp: i32 = exp.parse().unwrap();
    let digits: String = mant.chars().filter(|c| *c != '.').collect();
    let sign = if x < 0.0 { "-" } else { "" };
    let ax = x.abs();
    if (1e-3..1e7).contains(&ax) {
        // plain notation: at least one digit after the point
        let point = exp + 1; // digits before the decimal point
        let mut out = String::new();
        if point <= 0 {
            out.push_str("0.");
            for _ in 0..-point {
                out.push('0');
            }
            out.push_str(&digits);
        } else {
            let p = point as usize;
            if digits.len() <= p {
                out.push_str(&digits);
                for _ in digits.len()..p {
                    out.push('0');
                }
                out.push_str(".0");
            } else {
                out.push_str(&digits[..p]);
                out.push('.');
                out.push_str(&digits[p..]);
            }
        }
        format!("{sign}{out}")
    } else {
        let frac = if digits.len() > 1 { &digits[1..] } else { "0" };
        format!("{sign}{}.{}E{}", &digits[..1], frac, exp)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn to_string() {
        assert_eq!(java_double_to_string(-11.252334427406248), "-11.252334427406248");
        assert_eq!(java_double_to_string(0.405), "0.405");
        assert_eq!(java_double_to_string(1.0e-8), "1.0E-8");
        assert_eq!(java_double_to_string(1.2), "1.2");
        assert_eq!(java_double_to_string(100.0), "100.0");
        assert_eq!(java_double_to_string(0.001), "0.001");
        assert_eq!(java_double_to_string(1.0e7), "1.0E7");
        assert_eq!(java_double_to_string(1.5e-4), "1.5E-4");
    }
    /// Against commons-math / GATK (`somatic/javacheck/CmCheck.java` -> cm.txt)
    #[test]
    fn matches_commons_math() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../somatic/javacheck/cm.txt");
        let text = match std::fs::read_to_string(path) {
            Ok(t) => t,
            Err(_) => return,
        };
        let mut bad: std::collections::BTreeMap<&str, (usize, usize, Vec<String>)> = Default::default();
        let f = |s: &str| s.parse::<f64>().unwrap();
        let i = |s: &str| s.parse::<i32>().unwrap();
        let same = |a: f64, b: f64| a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan());
        for l in text.lines() {
            let t: Vec<&str> = l.split(' ').collect();
            let checks: Vec<(&str, f64, f64)> = match t[0] {
                "logBeta" => vec![("logBeta", log_beta(f(t[1]), f(t[2])), f(t[3]))],
                "bb" => vec![("betaBinomial", beta_binomial_log_probability(f(t[1]), f(t[2]), i(t[3]), i(t[4])), f(t[5]))],
                "binom" => vec![
                    ("binomProb", binomial_probability(i(t[1]), i(t[2]), f(t[3])), f(t[4])),
                    ("binomLogProb", binomial_log_probability(i(t[1]), i(t[2]), f(t[3])), f(t[5])),
                    ("binomCdf", binomial_cdf(i(t[1]), f(t[3]), i(t[2]) - 1), f(t[6])),
                ],
                "bcl" => vec![("binomCoefLog", binomial_coefficient_log(t[1].parse().unwrap(), t[2].parse().unwrap()), f(t[3]))],
                "rbeta" => vec![("regularizedBeta", regularized_beta(f(t[1]), f(t[2]), f(t[3])), f(t[4]))],
                "digamma" => vec![("digamma", hc_assembly::m2math::digamma(f(t[1])), f(t[2]))],
                _ => vec![],
            };
            for (name, got, want) in checks {
                let e = bad.entry(name).or_default();
                e.0 += 1;
                if !same(got, want) {
                    e.1 += 1;
                    if e.2.len() < 3 {
                        e.2.push(format!("{l} -> {got:?}"));
                    }
                }
            }
        }
        let mut fail = false;
        for (k, (n, b, ex)) in &bad {
            eprintln!("{k}: {b} / {n} differ");
            for e in ex {
                eprintln!("   {e}");
            }
            fail |= *b > 0;
        }
        assert!(!fail);
    }

    #[test]
    fn rounding() {
        assert_eq!(round_to_n_decimal_places(0.4049999, 3), 0.405);
        assert_eq!(round_to_n_decimal_places(1.2, 2), 1.2);
    }
}
