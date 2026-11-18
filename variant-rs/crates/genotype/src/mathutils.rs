//! GATK `MathUtils` routines used by genotyping, ported operation by operation.

use std::sync::OnceLock;

static JACOBIAN: OnceLock<Vec<f64>> = OnceLock::new();

/// `MathUtils.JacobianLogTable.cache`, dumped from GATK itself (`scratch/probes/Tables.java`, a development probe not in the repository).
fn jacobian_table() -> &'static [f64] {
    JACOBIAN.get_or_init(|| {
        let bytes = include_bytes!("jacobian_java.bin");
        bytes.chunks_exact(8).map(|c| f64::from_le_bytes(c.try_into().unwrap())).collect()
    })
}

const JACOBIAN_MAX_TOLERANCE: f64 = 8.0;
const JACOBIAN_INV_STEP: f64 = 1.0 / 0.0001;

/// `MathUtils.fastRound`
#[inline]
fn fast_round(d: f64) -> i32 {
    if d > 0.0 {
        (d + 0.5) as i32
    } else {
        (d - 0.5) as i32
    }
}

/// `MathUtils.approximateLog10SumLog10(a, b)`
#[inline]
pub fn approximate_log10_sum_log10(a: f64, b: f64, table: &[f64]) -> f64 {
    let (a, b) = if a > b { (b, a) } else { (a, b) };
    if a == f64::NEG_INFINITY {
        return b;
    }
    let diff = b - a;
    b + if diff < JACOBIAN_MAX_TOLERANCE { table[fast_round(diff * JACOBIAN_INV_STEP) as usize] } else { 0.0 }
}

pub fn jacobian() -> &'static [f64] {
    jacobian_table()
}

/// `MathUtils.maxElementIndex(array, start, end)` (first maximum)
#[inline]
pub fn max_element_index(a: &[f64]) -> usize {
    let mut m = 0;
    for i in 1..a.len() {
        if a[i] > a[m] {
            m = i;
        }
    }
    m
}

/// `MathUtils.log10SumLog10(double[])`
pub fn log10_sum_log10(v: &[f64]) -> f64 {
    if v.is_empty() {
        return f64::NEG_INFINITY;
    }
    let mi = max_element_index(v);
    let max = v[mi];
    if max == f64::NEG_INFINITY {
        return max;
    }
    let mut sum = 1.0f64;
    for (i, &c) in v.iter().enumerate() {
        if i == mi || c == f64::NEG_INFINITY {
            continue;
        }
        sum += hc_assembly::jmath::pow(10.0, c - max);
    }
    assert!(!(sum.is_nan() || sum == f64::INFINITY), "log10 p: Values must be non-infinite and non-NAN");
    max + if sum != 1.0 { hc_assembly::jmath::log10(sum) } else { 0.0 }
}

/// `MathUtils.log10SumLog10(a, b)`
#[inline]
pub fn log10_sum_log10_2(a: f64, b: f64) -> f64 {
    if a > b {
        a + hc_assembly::jmath::log10(1.0 + hc_assembly::jmath::pow(10.0, b - a))
    } else {
        b + hc_assembly::jmath::log10(1.0 + hc_assembly::jmath::pow(10.0, a - b))
    }
}

/// `MathUtils.normalizeLog10(array)` (log10 output, in place)
pub fn normalize_log10_in_place(v: &mut [f64]) {
    let s = log10_sum_log10(v);
    for x in v.iter_mut() {
        *x -= s;
    }
}

/// `MathUtils.normalizeFromLog10ToLinearSpace`
pub fn normalize_from_log10_to_linear(v: &[f64]) -> Vec<f64> {
    let s = log10_sum_log10(v);
    v.iter().map(|x| hc_assembly::jmath::pow(10.0, x - s)).collect()
}

/// `MathUtils.sum(double[])`
#[inline]
pub fn sum(v: &[f64]) -> f64 {
    let mut s = 0.0;
    for &x in v {
        s += x;
    }
    s
}

/// `Math.log10(2)` as used for `log10(ploidy)` and het combination counts (0x3fd34413509f79ff).
pub const LOG10_2: f64 = 0.3010299956639812;
