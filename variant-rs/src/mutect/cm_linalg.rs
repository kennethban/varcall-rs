//! commons-math 3.5 pieces used by CalculateContamination's kernel segmenter and statistics, ported with
//! Java's evaluation order: `FastMath.hypot` / `scalb`, `SingularValueDecomposition` (square / tall
//! matrices), `Percentile` (LEGACY estimation), `Mean`, and GATK's `PersistenceOptimizer`.
//! Also java.lang helpers: `Double.compare`, `Math.max` / `Math.min` (NaN-propagating).

use std::cmp::Ordering;

/// `Double.compare` (NaN greater than everything, -0.0 < 0.0).
pub fn java_compare(a: f64, b: f64) -> Ordering {
    if a < b {
        return Ordering::Less;
    }
    if a > b {
        return Ordering::Greater;
    }
    let canon = |x: f64| if x.is_nan() { 0x7ff8000000000000i64 } else { x.to_bits() as i64 };
    canon(a).cmp(&canon(b))
}

/// `Math.max(double, double)`
pub fn jmax(a: f64, b: f64) -> f64 {
    if a.is_nan() {
        return a;
    }
    if a == 0.0 && b == 0.0 && a.to_bits() == (-0.0f64).to_bits() {
        return b;
    }
    if a >= b { a } else { b }
}

/// `Math.min(double, double)`
pub fn jmin(a: f64, b: f64) -> f64 {
    if a.is_nan() {
        return a;
    }
    if a == 0.0 && b == 0.0 && b.to_bits() == (-0.0f64).to_bits() {
        return b;
    }
    if a <= b { a } else { b }
}

/// `FastMath.getExponent(double)`
fn get_exponent(d: f64) -> i32 {
    (((d.to_bits() >> 52) & 0x7ff) as i32) - 1023
}

/// `FastMath.scalb(double, int)`
pub fn scalb(d: f64, n: i32) -> f64 {
    if n > -1023 && n < 1024 {
        return d * f64::from_bits(((n + 1023) as u64) << 52);
    }
    if d.is_nan() || d.is_infinite() || d == 0.0 {
        return d;
    }
    if n < -2098 {
        return if d > 0.0 { 0.0 } else { -0.0 };
    }
    if n > 2097 {
        return if d > 0.0 { f64::INFINITY } else { f64::NEG_INFINITY };
    }
    let bits = d.to_bits();
    let sign = bits & 0x8000000000000000;
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    let mut mantissa = bits & 0x000fffffffffffff;
    let mut scaled = exponent + n;
    if n < 0 {
        if scaled > 0 {
            f64::from_bits(sign | ((scaled as u64) << 52) | mantissa)
        } else if scaled > -53 {
            mantissa |= 1u64 << 52;
            let msb_lost = mantissa & (1u64 << (-scaled));
            mantissa >>= 1 - scaled;
            if msb_lost != 0 {
                mantissa += 1;
            }
            f64::from_bits(sign | mantissa)
        } else if sign == 0 {
            0.0
        } else {
            -0.0
        }
    } else {
        if exponent == 0 {
            while (mantissa >> 52) != 1 {
                mantissa <<= 1;
                scaled -= 1;
            }
            scaled += 1;
            mantissa &= 0x000fffffffffffff;
            if scaled < 2047 {
                f64::from_bits(sign | ((scaled as u64) << 52) | mantissa)
            } else if sign == 0 {
                f64::INFINITY
            } else {
                f64::NEG_INFINITY
            }
        } else if scaled < 2047 {
            f64::from_bits(sign | ((scaled as u64) << 52) | mantissa)
        } else if sign == 0 {
            f64::INFINITY
        } else {
            f64::NEG_INFINITY
        }
    }
}

/// `FastMath.hypot`
pub fn hypot(x: f64, y: f64) -> f64 {
    if x.is_infinite() || y.is_infinite() {
        f64::INFINITY
    } else if x.is_nan() || y.is_nan() {
        f64::NAN
    } else {
        let ex = get_exponent(x);
        let ey = get_exponent(y);
        if ex > ey + 27 {
            x.abs()
        } else if ey > ex + 27 {
            y.abs()
        } else {
            let middle = (ex + ey) / 2;
            let sx = scalb(x, -middle);
            let sy = scalb(y, -middle);
            let sh = (sx * sx + sy * sy).sqrt();
            scalb(sh, middle)
        }
    }
}

/// `new SingularValueDecomposition(A)` for an m x n matrix with m >= n: returns (singular values, U (m x n)).
pub fn svd(mut a: Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
    const EPS: f64 = 2.220446049250313e-16; // 0x1.0p-52
    let tiny: f64 = f64::from_bits(((1023 - 966) as u64) << 52); // 0x1.0p-966
    let m = a.len();
    let n = a[0].len();
    assert!(m >= n, "svd: transposed case not ported");
    let mut s = vec![0.0f64; n];
    let mut u = vec![vec![0.0f64; n]; m];
    let mut v = vec![vec![0.0f64; n]; n];
    let mut e = vec![0.0f64; n];
    let mut work = vec![0.0f64; m];
    let nct = (m as i64 - 1).min(n as i64).max(0) as usize;
    let nrt = (n as i64 - 2).max(0) as usize;
    for k in 0..nct.max(nrt) {
        if k < nct {
            s[k] = 0.0;
            for i in k..m {
                s[k] = hypot(s[k], a[i][k]);
            }
            if s[k] != 0.0 {
                if a[k][k] < 0.0 {
                    s[k] = -s[k];
                }
                for i in k..m {
                    a[i][k] /= s[k];
                }
                a[k][k] += 1.0;
            }
            s[k] = -s[k];
        }
        for j in k + 1..n {
            if k < nct && s[k] != 0.0 {
                let mut t = 0.0;
                for i in k..m {
                    t += a[i][k] * a[i][j];
                }
                t = -t / a[k][k];
                for i in k..m {
                    a[i][j] += t * a[i][k];
                }
            }
            e[j] = a[k][j];
        }
        if k < nct {
            for i in k..m {
                u[i][k] = a[i][k];
            }
        }
        if k < nrt {
            e[k] = 0.0;
            for i in k + 1..n {
                e[k] = hypot(e[k], e[i]);
            }
            if e[k] != 0.0 {
                if e[k + 1] < 0.0 {
                    e[k] = -e[k];
                }
                for i in k + 1..n {
                    e[i] /= e[k];
                }
                e[k + 1] += 1.0;
            }
            e[k] = -e[k];
            if k + 1 < m && e[k] != 0.0 {
                for i in k + 1..m {
                    work[i] = 0.0;
                }
                for j in k + 1..n {
                    for i in k + 1..m {
                        work[i] += e[j] * a[i][j];
                    }
                }
                for j in k + 1..n {
                    let t = -e[j] / e[k + 1];
                    for i in k + 1..m {
                        a[i][j] += t * work[i];
                    }
                }
            }
            for i in k + 1..n {
                v[i][k] = e[i];
            }
        }
    }
    let mut p = n;
    if nct < n {
        s[nct] = a[nct][nct];
    }
    if m < p {
        s[p - 1] = 0.0;
    }
    if nrt + 1 < p {
        e[nrt] = a[nrt][p - 1];
    }
    e[p - 1] = 0.0;

    for j in nct..n {
        for i in 0..m {
            u[i][j] = 0.0;
        }
        u[j][j] = 1.0;
    }
    for k in (0..nct).rev() {
        if s[k] != 0.0 {
            for j in k + 1..n {
                let mut t = 0.0;
                for i in k..m {
                    t += u[i][k] * u[i][j];
                }
                t = -t / u[k][k];
                for i in k..m {
                    u[i][j] += t * u[i][k];
                }
            }
            for i in k..m {
                u[i][k] = -u[i][k];
            }
            u[k][k] = 1.0 + u[k][k];
            for i in 0..(k as i64 - 1).max(0) as usize {
                u[i][k] = 0.0;
            }
        } else {
            for i in 0..m {
                u[i][k] = 0.0;
            }
            u[k][k] = 1.0;
        }
    }
    for k in (0..n).rev() {
        if k < nrt && e[k] != 0.0 {
            for j in k + 1..n {
                let mut t = 0.0;
                for i in k + 1..n {
                    t += v[i][k] * v[i][j];
                }
                t = -t / v[k + 1][k];
                for i in k + 1..n {
                    v[i][j] += t * v[i][k];
                }
            }
        }
        for i in 0..n {
            v[i][k] = 0.0;
        }
        v[k][k] = 1.0;
    }

    let pp = p - 1;
    while p > 0 {
        let mut k: i64 = p as i64 - 2;
        while k >= 0 {
            let ku = k as usize;
            let threshold = tiny + EPS * (s[ku].abs() + s[ku + 1].abs());
            if !(e[ku].abs() > threshold) {
                e[ku] = 0.0;
                break;
            }
            k -= 1;
        }
        let kase;
        if k == p as i64 - 2 {
            kase = 4;
        } else {
            let mut ks: i64 = p as i64 - 1;
            while ks >= k {
                if ks == k {
                    break;
                }
                let ksu = ks as usize;
                let t = (if ks != p as i64 { e[ksu].abs() } else { 0.0 }) + (if ks != k + 1 { e[ksu - 1].abs() } else { 0.0 });
                if s[ksu].abs() <= tiny + EPS * t {
                    s[ksu] = 0.0;
                    break;
                }
                ks -= 1;
            }
            if ks == k {
                kase = 3;
            } else if ks == p as i64 - 1 {
                kase = 1;
            } else {
                kase = 2;
                k = ks;
            }
        }
        k += 1;
        let k = k as usize;
        match kase {
            1 => {
                let mut f = e[p - 2];
                e[p - 2] = 0.0;
                let mut j = p as i64 - 2;
                while j >= k as i64 {
                    let ju = j as usize;
                    let mut t = hypot(s[ju], f);
                    let cs = s[ju] / t;
                    let sn = f / t;
                    s[ju] = t;
                    if ju != k {
                        f = -sn * e[ju - 1];
                        e[ju - 1] = cs * e[ju - 1];
                    }
                    for i in 0..n {
                        t = cs * v[i][ju] + sn * v[i][p - 1];
                        v[i][p - 1] = -sn * v[i][ju] + cs * v[i][p - 1];
                        v[i][ju] = t;
                    }
                    j -= 1;
                }
            }
            2 => {
                let mut f = e[k - 1];
                e[k - 1] = 0.0;
                for j in k..p {
                    let mut t = hypot(s[j], f);
                    let cs = s[j] / t;
                    let sn = f / t;
                    s[j] = t;
                    f = -sn * e[j];
                    e[j] = cs * e[j];
                    for i in 0..m {
                        t = cs * u[i][j] + sn * u[i][k - 1];
                        u[i][k - 1] = -sn * u[i][j] + cs * u[i][k - 1];
                        u[i][j] = t;
                    }
                }
            }
            3 => {
                let max_pm1_pm2 = jmax(s[p - 1].abs(), s[p - 2].abs());
                let scale = jmax(jmax(jmax(max_pm1_pm2, e[p - 2].abs()), s[k].abs()), e[k].abs());
                let sp = s[p - 1] / scale;
                let spm1 = s[p - 2] / scale;
                let epm1 = e[p - 2] / scale;
                let sk = s[k] / scale;
                let ek = e[k] / scale;
                let b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
                let c = (sp * epm1) * (sp * epm1);
                let mut shift = 0.0;
                if b != 0.0 || c != 0.0 {
                    shift = (b * b + c).sqrt();
                    if b < 0.0 {
                        shift = -shift;
                    }
                    shift = c / (b + shift);
                }
                let mut f = (sk + sp) * (sk - sp) + shift;
                let mut g = sk * ek;
                for j in k..p - 1 {
                    let mut t = hypot(f, g);
                    let mut cs = f / t;
                    let mut sn = g / t;
                    if j != k {
                        e[j - 1] = t;
                    }
                    f = cs * s[j] + sn * e[j];
                    e[j] = cs * e[j] - sn * s[j];
                    g = sn * s[j + 1];
                    s[j + 1] = cs * s[j + 1];
                    for i in 0..n {
                        t = cs * v[i][j] + sn * v[i][j + 1];
                        v[i][j + 1] = -sn * v[i][j] + cs * v[i][j + 1];
                        v[i][j] = t;
                    }
                    t = hypot(f, g);
                    cs = f / t;
                    sn = g / t;
                    s[j] = t;
                    f = cs * e[j] + sn * s[j + 1];
                    s[j + 1] = -sn * e[j] + cs * s[j + 1];
                    g = sn * e[j + 1];
                    e[j + 1] = cs * e[j + 1];
                    if j < m - 1 {
                        for i in 0..m {
                            t = cs * u[i][j] + sn * u[i][j + 1];
                            u[i][j + 1] = -sn * u[i][j] + cs * u[i][j + 1];
                            u[i][j] = t;
                        }
                    }
                }
                e[p - 2] = f;
            }
            _ => {
                let mut k = k;
                if s[k] <= 0.0 {
                    s[k] = if s[k] < 0.0 { -s[k] } else { 0.0 };
                    for i in 0..=pp {
                        v[i][k] = -v[i][k];
                    }
                }
                while k < pp {
                    if s[k] >= s[k + 1] {
                        break;
                    }
                    s.swap(k, k + 1);
                    if k < n - 1 {
                        for i in 0..n {
                            let t = v[i][k + 1];
                            v[i][k + 1] = v[i][k];
                            v[i][k] = t;
                        }
                    }
                    if k < m - 1 {
                        for i in 0..m {
                            let t = u[i][k + 1];
                            u[i][k + 1] = u[i][k];
                            u[i][k] = t;
                        }
                    }
                    k += 1;
                }
                p -= 1;
            }
        }
    }
    (s, u)
}

/// `new Percentile(p).evaluate(values)` (LEGACY estimation, NaNs removed; values here are never NaN).
pub fn percentile(values: &[f64], p: f64) -> f64 {
    let n = values.len();
    if n == 0 {
        return f64::NAN;
    }
    if n == 1 {
        return values[0];
    }
    let mut w = values.to_vec();
    w.sort_by(|a, b| java_compare(*a, *b));
    let q = p / 100.0;
    let pos = if java_compare(q, 0.0) == Ordering::Equal {
        0.0
    } else if java_compare(q, 1.0) == Ordering::Equal {
        n as f64
    } else {
        q * (n as f64 + 1.0)
    };
    let fpos = pos.floor();
    let int_pos = fpos as usize;
    let dif = pos - fpos;
    if pos < 1.0 {
        return w[0];
    }
    if pos >= n as f64 {
        return w[n - 1];
    }
    let lower = w[int_pos - 1];
    let upper = w[int_pos];
    lower + dif * (upper - lower)
}

/// `new Mean().evaluate(values)` (definitional mean plus a correction pass).
pub fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let n = values.len() as f64;
    let mut sum = 0.0;
    for &x in values {
        sum += x;
    }
    let xbar = sum / n;
    let mut corr = 0.0;
    for &x in values {
        corr += x - xbar;
    }
    xbar + corr / n
}

/// GATK `PersistenceOptimizer(data).getMinimaIndices()`: global minimum first, then the other local minima
/// by decreasing persistence.
pub fn persistence_minima(data: &[f64]) -> Vec<usize> {
    const NO_COLOR: i64 = -1;
    struct Comp {
        left: usize,
        right: usize,
        min_index: usize,
        min_value: f64,
    }
    let n = data.len();
    assert!(n > 0);
    let mut sorted: Vec<usize> = (0..n).collect();
    sorted.sort_by(|&a, &b| java_compare(data[a], data[b]));
    // (minIndex, persistence)
    let mut pairs: Vec<(usize, f64)> = Vec::new();
    let pair = |i1: usize, i2: usize| -> (usize, f64) {
        let (mn, mx) = if data[i1] > data[i2] {
            (i2, i1)
        } else if data[i1] < data[i2] {
            (i1, i2)
        } else if i1 < i2 {
            (i1, i2)
        } else {
            (i2, i1)
        };
        (mn, data[mx] - data[mn])
    };
    if n > 1 {
        let mut comps: Vec<Comp> = Vec::new();
        let mut colors = vec![NO_COLOR; n];
        let create = |comps: &mut Vec<Comp>, colors: &mut Vec<i64>, index: usize| {
            colors[index] = comps.len() as i64;
            comps.push(Comp { left: index, right: index, min_index: index, min_value: data[index] });
        };
        let extend = |comps: &mut Vec<Comp>, colors: &mut Vec<i64>, ci: i64, index: usize| {
            let c = &mut comps[ci as usize];
            if index + 1 == c.left {
                c.left = index;
            } else if index >= 1 && index - 1 == c.right {
                c.right = index;
            }
            colors[index] = ci;
        };
        for &index in &sorted {
            if index == 0 {
                if colors[1] == NO_COLOR {
                    create(&mut comps, &mut colors, index);
                } else {
                    let c = colors[1];
                    extend(&mut comps, &mut colors, c, index);
                }
                continue;
            } else if index == n - 1 {
                if colors[index - 1] == NO_COLOR {
                    create(&mut comps, &mut colors, index);
                } else {
                    let c = colors[index - 1];
                    extend(&mut comps, &mut colors, c, index);
                }
                continue;
            }
            let lc = colors[index - 1];
            let rc = colors[index + 1];
            if lc == NO_COLOR && rc == NO_COLOR {
                create(&mut comps, &mut colors, index);
            } else if lc != NO_COLOR && rc == NO_COLOR {
                extend(&mut comps, &mut colors, lc, index);
            } else if lc == NO_COLOR && rc != NO_COLOR {
                extend(&mut comps, &mut colors, rc, index);
            } else {
                let (l, r) = (lc as usize, rc as usize);
                if comps[r].min_value < comps[l].min_value {
                    pairs.push(pair(comps[l].min_index, index));
                } else {
                    pairs.push(pair(comps[r].min_index, index));
                }
                // mergeComponents
                let (keep, merge) = if comps[l].min_value < comps[r].min_value {
                    (l, r)
                } else if comps[l].min_value > comps[r].min_value {
                    (r, l)
                } else if l < r {
                    (l, r)
                } else {
                    (r, l)
                };
                let (ml, mr, mmin) = (comps[merge].left, comps[merge].right, comps[merge].min_index);
                colors[ml] = keep as i64;
                colors[mr] = keep as i64;
                if comps[keep].min_index > mmin {
                    comps[keep].left = ml;
                } else {
                    comps[keep].right = mr;
                }
                colors[index] = colors[index - 1];
            }
        }
        pairs.sort_by(|a, b| java_compare(b.1, a.1));
    }
    let mut minima = vec![sorted[0]];
    minima.extend(pairs.iter().map(|p| p.0));
    minima
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn svd_small() {
        let (s, u) = svd(vec![vec![4.0, 1.0], vec![1.0, 3.0]]);
        assert!((s[0] - 4.618033988749895).abs() < 1e-12 && (s[1] - 2.381966011250105).abs() < 1e-12, "{s:?}");
        assert!((u[0][0].abs() - 0.8506508083520399).abs() < 1e-12);
    }
    #[test]
    fn percentile_legacy() {
        assert_eq!(percentile(&[1.0, 2.0, 3.0, 4.0], 50.0), 2.5);
        assert_eq!(percentile(&[5.0, 1.0, 3.0], 50.0), 3.0);
        assert_eq!(persistence_minima(&[3.0, 1.0, 2.0, 0.0, 5.0]), vec![3, 1]);
    }
}
