//! HaplotypeCaller's default annotations (StandardAnnotation + StandardHCAnnotation) for a
//! single diploid sample.

use hc_assembly::cigar;
use hc_assembly::Read;

use crate::javamath::{self, format_fixed, hypergeometric_log_probability, normal_inverse_cdf};
use crate::mathutils;

#[inline]
pub fn is_reverse_strand(r: &Read) -> bool {
    r.flag & 0x10 != 0
}

/// `ReadUtils.getSoftStart`
pub fn soft_start(r: &Read) -> i64 {
    let mut s = r.start;
    for &(op, len) in &r.cigar {
        if op == b'S' {
            s -= len as i64;
        } else if op != b'H' {
            break;
        }
    }
    s
}

/// `ReadUtils.getSoftEnd`
pub fn soft_end(r: &Read) -> i64 {
    let mut found = false;
    let mut e = r.end;
    for &(op, len) in r.cigar.iter().rev() {
        if op == b'S' {
            e += len as i64;
        } else if op != b'H' {
            found = true;
            break;
        }
    }
    if !found {
        r.end
    } else {
        e
    }
}

/// `ReadUtils.getReadIndexForReferenceCoordinate(softStart, cigar, refCoord)` -> (index, operator)
fn read_index_for_ref_coord(start: i64, c: &[(u8, u32)], ref_coord: i64) -> Option<(i64, u8)> {
    if ref_coord < start {
        return None;
    }
    let mut last_read = 0i64;
    let mut last_ref = start;
    for &(op, len) in c {
        let first_read = last_read;
        let first_ref = last_ref;
        if cigar::consumes_read(op) {
            last_read += len as i64;
        }
        if cigar::consumes_ref(op) || op == b'S' {
            last_ref += len as i64;
        }
        if first_ref <= ref_coord && ref_coord < last_ref {
            let pos = first_read + if cigar::consumes_read(op) { ref_coord - first_ref } else { 0 };
            return Some((pos, op));
        }
    }
    None
}

/// `BaseQualityRankSumTest.getReadBaseQuality`
pub fn base_quality_at(r: &Read, ref_coord: i64) -> Option<f64> {
    if ref_coord < r.start || r.end < ref_coord {
        return None;
    }
    match read_index_for_ref_coord(soft_start(r), &r.cigar, ref_coord) {
        Some((i, op)) if cigar::consumes_read(op) => Some(r.qual[i as usize] as f64),
        _ => None,
    }
}

/// `ReadPosRankSumTest.getReadPosition`
pub fn read_position(r: &Read, vc_start: i64, vc_end: i64) -> Option<f64> {
    if r.start == vc_end + 1 {
        let first_non_clip = r.cigar.iter().map(|e| e.0).find(|&op| !cigar::is_clipping(op));
        if first_non_clip == Some(b'I') {
            return Some(0.0);
        }
    }
    let (offset, _) = read_index_for_ref_coord(soft_start(r), &r.cigar, vc_start)?;
    let first = r.cigar[0];
    let last = r.cigar[r.cigar.len() - 1];
    let lead_h = if first.0 == b'H' { first.1 as i64 } else { 0 };
    let trail_h = if last.0 == b'H' { last.1 as i64 } else { 0 };
    let left = lead_h + offset;
    let right = r.seq.len() as i64 - 1 - offset + trail_h;
    Some(left.min(right) as f64)
}

/// 2x2 strand table `[[ref fwd, ref rev], [alt fwd, alt rev]]` (`StrandBiasTest.getContingencyTable`).
/// `best[r]`: best allele of each evidence (None if not informative); `is_ref_allele` / `is_alt_allele`
/// classify allele indices.
pub fn strand_table(
    best: &[Option<usize>],
    reads: &[&Read],
    is_ref_allele: &dyn Fn(usize) -> bool,
    is_alt_allele: &dyn Fn(usize) -> bool,
    min_count: i32,
) -> [[i32; 2]; 2] {
    let mut t = [0i32; 4];
    for (b, r) in best.iter().zip(reads.iter()) {
        if let Some(a) = *b {
            let mr = is_ref_allele(a);
            if mr || is_alt_allele(a) {
                let off = if mr { 0 } else { 2 };
                t[off + if is_reverse_strand(r) { 1 } else { 0 }] += 1;
            }
        }
    }
    if t[0] + t[1] + t[2] + t[3] > min_count {
        [[t[0], t[1]], [t[2], t[3]]]
    } else {
        [[0, 0], [0, 0]]
    }
}

/// `FisherExactTest.twoSidedPValue`
fn fisher_two_sided(x: [[i32; 2]; 2]) -> f64 {
    let m = x[0][0] + x[0][1];
    let n = x[1][0] + x[1][1];
    let k = x[0][0] + x[1][0];
    let lo = 0.max(k - n);
    let hi = k.min(m);
    if hi + 1 - lo <= 1 {
        return 1.0;
    }
    let logds: Vec<f64> = (lo..=hi).map(|i| hypergeometric_log_probability(m + n, m, k, i)).collect();
    let rel_err = 1.0 - 10e-7;
    let threshold = logds[(x[0][0] - lo) as usize] * rel_err;
    let log10ds: Vec<f64> =
        logds.iter().filter(|&&d| d <= threshold).map(|&d| d * std::f64::consts::LOG10_E).collect();
    let p = hc_assembly::jmath::pow(10.0, mathutils::log10_sum_log10(&log10ds));
    p.min(1.0)
}

/// FisherStrand: "%.3f" of the phred-scaled two-sided p-value of the (normalised) table.
pub fn fisher_strand(table: [[i32; 2]; 2]) -> String {
    let sum = table[0][0] + table[0][1] + table[1][0] + table[1][1];
    let t = if (sum as f64) <= 200.0 * 2.0 {
        table
    } else {
        crate::stats::inc(&crate::stats::FS_TABLE_NORMALIZED);
        let f = sum as f64 / 200.0;
        [
            [(table[0][0] as f64 / f) as i32, (table[0][1] as f64 / f) as i32],
            [(table[1][0] as f64 / f) as i32, (table[1][1] as f64 / f) as i32],
        ]
    };
    let p = fisher_two_sided(t);
    let p = p.max(1E-320);
    // QualityUtils.phredScaleErrorRate -> phredScaleLog10ErrorRate
    let min_log10 = f64::from_bits(1).log10(); // Math.log10(Double.MIN_VALUE)
    let q = (-10.0 * hc_assembly::jmath::log10(p).max(min_log10)).abs();
    format_fixed(q, 3)
}

/// StrandOddsRatio: "%.3f"
pub fn strand_odds_ratio(table: [[i32; 2]; 2]) -> String {
    let t00 = table[0][0] as f64 + 1.0;
    let t01 = table[0][1] as f64 + 1.0;
    let t11 = table[1][1] as f64 + 1.0;
    let t10 = table[1][0] as f64 + 1.0;
    let ratio = (t00 / t01) * (t11 / t10) + (t01 / t00) * (t10 / t11);
    let ref_ratio = t00.min(t01) / t00.max(t01);
    let alt_ratio = t10.min(t11) / t10.max(t11);
    let sor = hc_assembly::jmath::log(ratio) + hc_assembly::jmath::log(ref_ratio) - hc_assembly::jmath::log(alt_ratio);
    format_fixed(sor, 3)
}

// ---------------------------------------------------------------------------------------------
// MannWhitneyU (one-sided, FIRST_DOMINATES)

struct Ranked {
    /// (value, rank, series) sorted by value
    ranks: Vec<(f64, f32, u8)>,
    ties: Vec<i32>,
}

fn calculate_rank(s1: &mut [f64], s2: &mut [f64]) -> Ranked {
    s1.sort_by(|a, b| hc_assembly::jmath::java_cmp(*a, *b));
    s2.sort_by(|a, b| hc_assembly::jmath::java_cmp(*a, *b));
    let total = s1.len() + s2.len();
    let mut ranks: Vec<(f64, f32, u8)> = Vec::with_capacity(total);
    let (mut i, mut j) = (0usize, 0usize);
    while ranks.len() < total {
        let r = (ranks.len() + 1) as f32;
        if i >= s1.len() {
            ranks.push((s2[j], r, 2));
            j += 1;
        } else if j >= s2.len() {
            ranks.push((s1[i], r, 1));
            i += 1;
        } else if s1[i] <= s2[j] {
            ranks.push((s1[i], r, 1));
            i += 1;
        } else {
            ranks.push((s2[j], r, 2));
            j += 1;
        }
    }
    let mut ties = Vec::new();
    let mut i = 0;
    while i < ranks.len() {
        let mut rank = ranks[i].1;
        let mut count = 1;
        let mut j = i + 1;
        while j < ranks.len() && ranks[j].0 == ranks[i].0 {
            rank += ranks[j].1;
            count += 1;
            j += 1;
        }
        if count > 1 {
            rank /= count as f32;
            for k in i..i + count {
                ranks[k].1 = rank;
            }
            ties.push(count as i32);
        }
        i += count;
    }
    Ranked { ranks, ties }
}

/// `MannWhitneyU.test(series1, series2, FIRST_DOMINATES).getZ()`
pub fn mann_whitney_z(series1: &[f64], series2: &[f64]) -> f64 {
    let n1 = series1.len();
    let n2 = series2.len();
    if n1 == 0 || n2 == 0 {
        return f64::NAN;
    }
    let mut a = series1.to_vec();
    let mut b = series2.to_vec();
    let ranked = calculate_rank(&mut a, &mut b);
    let n_ranks = ranked.ranks.len() as i32;
    let mut nties = 0.0f64;
    for &c in &ranked.ties {
        if c != n_ranks {
            nties += (c as f64).powi(3) - c as f64;
        }
    }
    let mut r1 = 0f32;
    let mut r2 = 0f32;
    for &(_, r, s) in &ranked.ranks {
        if s == 1 {
            r1 += r;
        } else {
            r2 += r;
        }
    }
    let _ = r2;
    let dn1 = n1 as f64;
    let u = r1 as f64 - (dn1 * (dn1 + 1.0)) / 2.0;
    if n1 >= 10 || n2 >= 10 {
        crate::stats::inc(&crate::stats::MWU_NORMAL);
        let (n1i, n2i) = (n1 as i64, n2 as i64);
        let m = (n1i * n2i) as f64 / 2.0;
        let correction = if nties == 0.0 { 0.0 } else { -0.5 };
        let sigma = (((n1i * n2i) as f64 / 12.0)
            * ((n1i + n2i + 1) as f64 - nties / ((n1i + n2i) * (n1i + n2i - 1)) as f64))
            .sqrt();
        (u - m - correction) / sigma
    } else {
        crate::stats::inc(&crate::stats::MWU_EXACT);
        let p = permutation_test(&ranked, n1, n2, u);
        normal_inverse_cdf(p)
    }
}

/// `MannWhitneyU.permutationTest`: exact distribution of U over all distinct assignments of the
/// ranks to the two series (computed by counting subsets instead of enumerating permutations).
fn permutation_test(ranked: &Ranked, n1: usize, _n2: usize, u: f64) -> f64 {
    // ranks are multiples of 1/2: work with 2*rank
    let r2: Vec<i64> = ranked.ranks.iter().map(|&(_, r, _)| (r as f64 * 2.0).round() as i64).collect();
    let max_sum: i64 = r2.iter().sum();
    let width = (max_sum + 1) as usize;
    // dp[k][s] = number of subsets of size k with doubled rank sum s
    let mut dp = vec![vec![0f64; width]; n1 + 1];
    dp[0][0] = 1.0;
    for &x in &r2 {
        for k in (1..=n1).rev() {
            let (lo, hi) = dp.split_at_mut(k);
            let prev = &lo[k - 1];
            let cur = &mut hi[0];
            for s in (x as usize..width).rev() {
                if prev[s - x as usize] != 0.0 {
                    cur[s] += prev[s - x as usize];
                }
            }
        }
    }
    let offset = (n1 * (n1 + 1)) as f64 / 2.0;
    let key = javamath::java_round(2.0 * u);
    let mut total = 0f64;
    let mut at = 0f64;
    let mut below = 0f64;
    for (s, &c) in dp[n1].iter().enumerate() {
        if c == 0.0 {
            continue;
        }
        let new_u = s as f64 / 2.0 - offset;
        let k = javamath::java_round(2.0 * new_u);
        total += c;
        if k == key {
            at += c;
        } else if k < key {
            below += c;
        }
    }
    (at / 2.0 + below) / total
}

/// RankSumTest: "%.3f" of the z-score, or None (missing) when NaN.
pub fn rank_sum(alt: &[f64], refv: &[f64]) -> Option<String> {
    if alt.is_empty() && refv.is_empty() {
        return None;
    }
    let z = mann_whitney_z(alt, refv);
    if z.is_nan() {
        None
    } else {
        Some(format_fixed(z, 3))
    }
}

/// Per-read values for the three rank sum tests.
pub fn base_qual_value(r: &Read, vc_start: i64) -> Option<f64> {
    base_quality_at(r, vc_start)
}
pub fn read_pos_value(r: &Read, vc_start: i64, vc_end: i64) -> Option<f64> {
    if !(soft_start(r) <= vc_end + 1 && soft_end(r) >= vc_start) {
        return None;
    }
    read_position(r, vc_start, vc_end)
}

// ---------------------------------------------------------------------------------------------
// ExcessHet (single sample)

pub fn exact_test(het: i32, refc: i32, hom: i32) -> f64 {
    let (obs_hom_r, obs_hom_c) = if refc < hom { (refc, hom) } else { (hom, refc) };
    let rare = 2 * obs_hom_r + het;
    let n = het + obs_hom_c + obs_hom_r;
    let mut probs = vec![0f64; (rare + 1) as usize];
    let mut mid = (rare as f64 * (2.0 * n as f64 - rare as f64) / (2.0 * n as f64 - 1.0)).floor() as i32;
    if mid % 2 != rare % 2 {
        mid += 1;
    }
    probs[mid as usize] = 1.0;
    let mut mysum = 1.0;
    let (mut ch, mut chr, mut chc) = (mid, (rare - mid) / 2, n - mid - (rare - mid) / 2);
    while ch >= 2 {
        let p = probs[ch as usize] * ch as f64 * (ch as f64 - 1.0) / (4.0 * (chr as f64 + 1.0) * (chc as f64 + 1.0));
        if p < 1.0E-16 {
            break;
        }
        probs[(ch - 2) as usize] = p;
        mysum += p;
        ch -= 2;
        chr += 1;
        chc += 1;
    }
    ch = mid;
    chr = (rare - mid) / 2;
    chc = n - ch - chr;
    while ch <= rare - 2 {
        let p = probs[ch as usize] * 4.0 * chr as f64 * chc as f64 / ((ch as f64 + 2.0) * (ch as f64 + 1.0));
        if p < 1.0E-16 {
            break;
        }
        probs[(ch + 2) as usize] = p;
        mysum += p;
        ch += 2;
        chr -= 1;
        chc -= 1;
    }
    let right = probs[het as usize] / mysum;
    if het == rare {
        return right.min(1.0).max(0.0);
    }
    let mut tail = 0.0;
    for &p in &probs[(het + 1) as usize..] {
        tail += p;
    }
    (right + tail / mysum).min(1.0).max(0.0)
}

/// ExcessHet ("%.4f") for one called diploid genotype with PLs over `n_alleles` alleles.
pub fn excess_het(pls: &[i32], n_alleles: usize) -> String {
    excess_het_samples(&[(true, pls)], n_alleles).unwrap()
}

/// `GenotypeUtils.computeDiploidGenotypeCounts` over diploid genotypes that all carry PLs (`(called, pls)`):
/// [refs, hets, homs], each genotype's contribution rounded (ExcessHet) or not (InbreedingCoeff).
pub fn diploid_genotype_counts(gs: &[(bool, &[i32])], n_alleles: usize, round: bool) -> [f64; 3] {
    let mut counts = [0f64; 3];
    for &(_, pls) in gs {
        let gls: Vec<f64> = pls.iter().map(|&p| p as f64 / -10.0).collect();
        let norm = mathutils::normalize_from_log10_to_linear(&gls);
        let (r, h, v);
        if n_alleles > 2 {
            let max_ind = {
                let mut m = 0;
                for i in 1..norm.len() {
                    if norm[i] > norm[m] {
                        m = i;
                    }
                }
                m
            };
            let gts = crate::genotyper::genotype_alleles(n_alleles);
            let (a1, a2) = gts[max_ind];
            if a1 != 0 && a2 != 0 {
                counts[2] += 1.0;
                continue;
            }
            let mut max_l = norm[1];
            let mut het_i = 1;
            let mut var_i = 2;
            for a in 1..n_alleles {
                let hi = crate::genotyper::genotype_index(0, a);
                if norm[hi] > max_l {
                    max_l = norm[hi];
                    het_i = hi;
                    var_i = crate::genotyper::genotype_index(a, a);
                }
            }
            let b = [norm[0], norm[het_i], norm[var_i]];
            let s = mathutils::sum(&b);
            r = b[0] / s;
            h = b[1] / s;
            v = b[2] / s;
        } else {
            r = norm[0];
            h = norm[1];
            v = norm[2];
        }
        if round {
            let fr = |d: f64| if d > 0.0 { (d + 0.5) as i32 } else { (d - 0.5) as i32 };
            counts[0] += fr(r) as f64;
            if r != h {
                counts[1] += fr(h) as f64;
            }
            if v != h {
                counts[2] += fr(v) as f64;
            }
        } else {
            counts[0] += r;
            counts[1] += h;
            counts[2] += v;
        }
    }
    counts
}

/// `ExcessHet.annotate` over the samples' genotypes (`(called, pls)`); None without a called genotype.
pub fn excess_het_samples(gs: &[(bool, &[i32])], n_alleles: usize) -> Option<String> {
    if !gs.iter().any(|g| g.0) {
        return None;
    }
    let counts = diploid_genotype_counts(gs, n_alleles, true);
    let p = exact_test(counts[1] as i32, counts[0] as i32, counts[2] as i32);
    let eh = if p < 10e-60 { -10.0 * hc_assembly::jmath::log10(1.0E-16) } else { -10.0 * hc_assembly::jmath::log10(p) + 0.0 };
    Some(format_fixed(eh, 4))
}

// ---------------------------------------------------------------------------------------------------------
// Allele-specific annotations (`-G AS_StandardAnnotation`). Alleles are positions in the call's allele list
// (reference first); `best[e]` is the call-allele position of evidence e's best allele over all likelihood
// alleles (None: not informative, or best allele not among the call's alleles).
// ---------------------------------------------------------------------------------------------------------

/// `Histogram.toString()` of a histogram holding one value (bin 0.1, "%.1f,count"); "NaN" when empty
/// (`AS_RankSumTest.outputSingletonValueAsHistogram`: NaN is not added).
pub fn histogram_singleton(z: f64) -> String {
    if z.is_nan() {
        return "NaN".to_string();
    }
    let key = crate::javamath::java_round(((z + 0.01 * 0.1) / 0.1).floor());
    assert!(key <= i32::MAX as i64 && key >= i32::MIN as i64, "Histogram values are suspiciously extreme.  Failed to add {z} to the Histogram.");
    format!("{},1", format_fixed(key as i32 as f64 * 0.1, 1))
}

/// `AS_RankSumTest.annotateRawData`: per-allele values (already `(int)`-truncated) -> `|`-joined singleton
/// histograms of each alt's Mann-Whitney z against the reference (nothing at all when the reference has no value).
pub fn as_raw_rank_sum(values: &[Vec<f64>]) -> String {
    let mut s = String::new();
    for (i, alt) in values.iter().enumerate().skip(1) {
        if i != 0 {
            s.push('|');
        }
        if !values[0].is_empty() {
            s.push_str(&histogram_singleton(mann_whitney_z(alt, &values[0])));
        }
    }
    s
}

/// `AS_RMSMappingQuality.makeRawAnnotationString`: "%.2f" of each allele's sum of squared MQs, `|`-joined.
pub fn as_raw_mq(sum_sq: &[Option<f64>]) -> String {
    sum_sq.iter().map(|v| format_fixed(v.unwrap_or(0.0), 2)).collect::<Vec<_>>().join("|")
}

/// `StrandBiasUtils.computeSBAnnotation` (AS_SB_TABLE): per allele "fwd,rev", `|`-joined; all zero unless more
/// than `MIN_COUNT` (2) reads count.
pub fn as_sb_table(best: &[Option<usize>], reads: &[&Read], n_alleles: usize) -> String {
    let mut t = vec![[0i32; 2]; n_alleles];
    for (b, r) in best.iter().zip(reads) {
        if let Some(a) = *b {
            t[a][if is_reverse_strand(r) { 1 } else { 0 }] += 1;
        }
    }
    let total: i32 = t.iter().map(|x| x[0] + x[1]).sum();
    if total <= 2 {
        t.iter_mut().for_each(|x| *x = [0, 0]);
    }
    t.iter().map(|x| format!("{},{}", x[0], x[1])).collect::<Vec<_>>().join("|")
}

#[cfg(test)]
mod as_tests {
    use super::*;
    #[test]
    fn histograms() {
        assert_eq!(histogram_singleton(f64::NAN), "NaN");
        assert_eq!(histogram_singleton(1.234), "1.2,1");
        assert_eq!(histogram_singleton(-0.35), "-0.4,1");
        assert_eq!(histogram_singleton(0.0), "0.0,1");
        assert_eq!(histogram_singleton(-0.0001), "0.0,1");
        assert_eq!(as_raw_rank_sum(&[vec![], vec![30.0], vec![]]), "||");
    }
}
