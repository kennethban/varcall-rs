//! `-contamination`: `AlleleLikelihoods.contaminationDownsampling` for one sample, applied by
//! `HaplotypeCallerGenotypingEngine` to the likelihoods used for genotyping (not to those used for annotation).

use crate::likelihoods::best_allele;
use crate::well::Well19937c;

/// `AlleleBiasedDownsamplingUtils.scoreAlleleCounts`
fn score(counts: &[i32]) -> i32 {
    if counts.len() < 2 {
        return 0;
    }
    let mut c = counts.to_vec();
    c.sort_unstable();
    let max = c[c.len() - 1];
    let next = c[c.len() - 2];
    let rest = c.iter().map(|&x| x as i64).sum::<i64>() as i32 - max - next;
    (max - next + rest).min((next + rest).abs())
}

/// `AlleleBiasedDownsamplingUtils.targetAlleleCounts`
fn target_allele_counts(counts: &[i32], remove: i32) -> Vec<i32> {
    let mut best_score = score(counts);
    let mut best = counts.to_vec();
    let per = remove / 2;
    for i in 0..counts.len() {
        for j in i..counts.len() {
            let mut c = counts.to_vec();
            if i == j {
                c[i] = (c[i] - remove).max(0);
            } else {
                c[i] = (c[i] - per).max(0);
                c[j] = (c[j] - per).max(0);
            }
            let s = score(&c);
            if s < best_score {
                best_score = s;
                best = c;
            }
        }
    }
    best
}

/// Evidence (columns of `lk`, rows = the event's alleles without `<NON_REF>`) that contamination downsampling
/// removes, sorted. `fraction` > 0. Draws from `rng` in GATK's order.
pub fn removals(lk: &[Vec<f64>], n_ev: usize, fraction: f64, rng: &mut Well19937c) -> Vec<usize> {
    if fraction >= 1.0 {
        return (0..n_ev).collect();
    }
    // evidenceByBestAlleleMap: informative evidence per best allele (canBeReference, no priorities)
    let n = lk.len();
    let mut bins: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut col = vec![0f64; n];
    for e in 0..n_ev {
        for a in 0..n {
            col[a] = lk[a][e];
        }
        let b = best_allele(&col, None);
        if b.is_informative() {
            bins[b.allele.unwrap()].push(e);
        }
    }
    // selectAlleleBiasedEvidence
    let total: i32 = bins.iter().map(|b| b.len() as i32).sum();
    let to_remove = (total as f64 * fraction) as i32;
    let counts: Vec<i32> = bins.iter().map(|b| b.len() as i32).collect();
    let target = target_allele_counts(&counts, to_remove);
    let mut out = Vec::new();
    for i in 0..n {
        if counts[i] > target[i] {
            // downsampleElements
            let k = (counts[i] - target[i]) as usize;
            if k >= bins[i].len() {
                out.extend_from_slice(&bins[i]);
            } else {
                out.extend(rng.next_permutation(bins[i].len(), k).into_iter().map(|x| bins[i][x]));
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

/// The columns of `lk` without `removed` (sorted), in order (`removeEvidenceByIndex`).
pub fn keep_columns(lk: &[Vec<f64>], removed: &[usize]) -> Vec<Vec<f64>> {
    lk.iter()
        .map(|row| {
            let mut k = 0;
            row.iter()
                .enumerate()
                .filter(|&(e, _)| {
                    if k < removed.len() && removed[k] == e {
                        k += 1;
                        false
                    } else {
                        true
                    }
                })
                .map(|(_, &v)| v)
                .collect()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn target_counts() {
        // worked through GATK's loop by hand: the first strictly better score wins
        assert_eq!(target_allele_counts(&[10, 10], 2), vec![10, 10]);
        assert_eq!(target_allele_counts(&[100, 2], 10), vec![95, 0]);
        assert_eq!(target_allele_counts(&[100, 20, 5], 10), vec![100, 10, 5]);
    }
}
