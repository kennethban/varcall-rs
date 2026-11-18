//! PairHMM input preparation of `PairHMMLikelihoodCalculationEngine.modifyReadQualities`
//! (PCR indel error model + `capMinimumReadQualities`) and the standard input-score imputator.

use crate::trim::find_number_of_repetitions;
use crate::Params;

pub const MAX_STR_UNIT_LENGTH: usize = 8;
pub const MAX_REPEAT_LENGTH: usize = 20;

/// `PairHMMLikelihoodCalculationEngine.PCRErrorModel`
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PcrErrorModel {
    None,
    Hostile,
    Aggressive,
    Conservative,
}

impl PcrErrorModel {
    pub fn rate_factor(self) -> f64 {
        match self {
            PcrErrorModel::None => 0.0,
            PcrErrorModel::Hostile => 1.0,
            PcrErrorModel::Aggressive => 2.0,
            PcrErrorModel::Conservative => 3.0,
        }
    }
}

/// `getErrorModelAdjustedQual(repeatLength, rateFactor)` =
/// `max(10, fastRound(40 - exp(repeatLength / (rateFactor * PI)) + 1))`.
pub fn pcr_indel_table(model: PcrErrorModel) -> Option<[u8; MAX_REPEAT_LENGTH + 1]> {
    if model == PcrErrorModel::None {
        return None;
    }
    let rf = model.rate_factor();
    let mut t = [0u8; MAX_REPEAT_LENGTH + 1];
    for (i, v) in t.iter_mut().enumerate() {
        let d = 40.0 - (i as f64 / (rf * std::f64::consts::PI)).exp() + 1.0;
        let r = if d > 0.0 { (d + 0.5) as i32 } else { (d - 0.5) as i32 };
        *v = r.max(10) as u8;
    }
    Some(t)
}

/// `ReadLikelihoodCalculationEngine.findTandemRepeatUnits(readBases, offset).getRight()` (the
/// repeat length, capped at MAX_REPEAT_LENGTH).  Direct transcription; see [`tandem_repeat_length`]
/// for the fast equivalent used in production.
pub fn tandem_repeat_length_reference(b: &[u8], offset: usize) -> usize {
    let mut max_bw: usize = 0;
    let mut bw_unit: &[u8] = &b[offset..offset + 1];
    for s in 1..=MAX_STR_UNIT_LENGTH {
        if offset + 1 < s {
            break;
        }
        // findNumberOfRepetitions(readBases, offset-str+1, str, readBases, 0, offset+1, false)
        max_bw = find_number_of_repetitions(&b[offset + 1 - s..offset + 1], &b[..offset + 1], false);
        if max_bw > 1 {
            bw_unit = &b[offset + 1 - s..offset + 1];
            break;
        }
    }
    let mut max_rl = max_bw;
    if offset < b.len() - 1 {
        let mut fw_unit: &[u8] = &b[offset + 1..offset + 2];
        let mut max_fw: usize = 0;
        for s in 1..=MAX_STR_UNIT_LENGTH {
            if offset + s + 1 > b.len() {
                break;
            }
            max_fw = find_number_of_repetitions(&b[offset + 1..offset + 1 + s], &b[offset + 1..], true);
            if max_fw > 1 {
                fw_unit = &b[offset + 1..offset + s + 1];
                break;
            }
        }
        if fw_unit == bw_unit {
            max_rl = max_bw + max_fw;
        } else {
            let bw2 = find_number_of_repetitions(fw_unit, &b[..offset + 1], false);
            max_rl = max_fw + bw2;
        }
    }
    max_rl.min(MAX_REPEAT_LENGTH)
}

/// `run / s` without a division in the common case `run < s`.
#[inline(always)]
fn div_small(run: usize, s: usize) -> usize {
    if run < s { 0 } else { run / s }
}

/// Number of consecutive `j = from, from-1, ..., 0` with `b[j] == b[j + s]`.
#[inline(always)]
fn run_back(b: &[u8], from: i64, s: usize) -> usize {
    let mut j = from;
    let mut m = 0;
    while j >= 0 && b[j as usize] == b[j as usize + s] {
        m += 1;
        j -= 1;
    }
    m
}

/// Number of consecutive `i = from, from+1, ...` (with `i + s < len`) with `b[i] == b[i + s]`.
#[inline(always)]
fn run_fwd(b: &[u8], from: usize, s: usize) -> usize {
    let mut i = from;
    while i + s < b.len() && b[i] == b[i + s] {
        i += 1;
    }
    i - from
}

/// Same result as [`tandem_repeat_length_reference`], computed with shifted self-comparisons:
/// `k` back-to-back copies of a unit of length `s` are exactly a run of `(k - 1) * s` positions
/// with `b[j] == b[j + s]`, so each `findNumberOfRepetitions` call becomes one tight byte loop.
pub fn tandem_repeat_length(b: &[u8], offset: usize) -> usize {
    let o = offset;
    // backward: unit b[o+1-s..=o], trailing repeats in b[..=o]
    let mut max_bw = 0usize;
    let mut bw_s = 1usize; // bestBWRepeatUnit = b[o..=o] unless a repeat is found
    for s in 1..=MAX_STR_UNIT_LENGTH {
        if o + 1 < s {
            break;
        }
        max_bw = 1 + div_small(run_back(b, o as i64 - s as i64, s), s);
        if max_bw > 1 {
            bw_s = s;
            break;
        }
    }
    let mut max_rl = max_bw;
    if o < b.len() - 1 {
        let mut fw_s = 1usize;
        let mut max_fw = 0usize;
        for s in 1..=MAX_STR_UNIT_LENGTH {
            if o + s + 1 > b.len() {
                break;
            }
            max_fw = 1 + div_small(run_fwd(b, o + 1, s), s);
            if max_fw > 1 {
                fw_s = s;
                break;
            }
        }
        let fw_unit = &b[o + 1..o + 1 + fw_s];
        let bw_unit = &b[o + 1 - bw_s..o + 1];
        if fw_unit == bw_unit {
            max_rl = max_bw + max_fw;
        } else {
            // trailing repeats of the forward unit in b[..=o]
            max_rl = max_fw + div_small(run_back(b, o as i64, fw_s), fw_s);
        }
    }
    max_rl.min(MAX_REPEAT_LENGTH)
}

/// Repeat lengths of every offset of a read at once (same values as [`tandem_repeat_length`]):
/// the runs of `b[j] == b[j + s]` are tabulated once per unit length, so each offset costs a
/// few table lookups.
pub struct RepeatTables {
    /// back[s-1][j] = run of matches at j, j-1, ...; fwd[s-1][i] = run at i, i+1, ...
    back: Vec<Vec<u16>>,
    fwd: Vec<Vec<u16>>,
}

impl Default for RepeatTables {
    fn default() -> Self {
        RepeatTables { back: vec![Vec::new(); MAX_STR_UNIT_LENGTH], fwd: vec![Vec::new(); MAX_STR_UNIT_LENGTH] }
    }
}

impl RepeatTables {
    /// Fills `out[o]` with the repeat length at offset `o` for `o in 0..b.len()-1`.
    pub fn repeat_lengths(&mut self, b: &[u8], out: &mut Vec<u8>) {
        let n = b.len();
        out.clear();
        if n < 2 {
            return;
        }
        for s in 1..=MAX_STR_UNIT_LENGTH {
            let back = &mut self.back[s - 1];
            let fwd = &mut self.fwd[s - 1];
            back.clear();
            back.resize(n, 0);
            fwd.clear();
            fwd.resize(n + 1, 0);
            if s < n {
                let mut run = 0u16;
                for j in 0..n - s {
                    run = if b[j] == b[j + s] { run + 1 } else { 0 };
                    back[j] = run;
                }
                let mut run = 0u16;
                for i in (0..n - s).rev() {
                    run = if b[i] == b[i + s] { run + 1 } else { 0 };
                    fwd[i] = run;
                }
            }
        }
        let rb = |from: i64, s: usize| -> usize { if from < 0 { 0 } else { self.back[s - 1][from as usize] as usize } };
        let rf = |from: usize, s: usize| -> usize { self.fwd[s - 1][from] as usize };
        for o in 0..n - 1 {
            let mut max_bw = 0usize;
            let mut bw_s = 1usize;
            for s in 1..=MAX_STR_UNIT_LENGTH {
                if o + 1 < s {
                    break;
                }
                max_bw = 1 + div_small(rb(o as i64 - s as i64, s), s);
                if max_bw > 1 {
                    bw_s = s;
                    break;
                }
            }
            // o < n - 1 always holds here
            let mut fw_s = 1usize;
            let mut max_fw = 0usize;
            for s in 1..=MAX_STR_UNIT_LENGTH {
                if o + s + 1 > n {
                    break;
                }
                max_fw = 1 + div_small(rf(o + 1, s), s);
                if max_fw > 1 {
                    fw_s = s;
                    break;
                }
            }
            let max_rl = if b[o + 1..o + 1 + fw_s] == b[o + 1 - bw_s..o + 1] { max_bw + max_fw } else { max_fw + div_small(rb(o as i64, fw_s), fw_s) };
            out.push(max_rl.min(MAX_REPEAT_LENGTH) as u8);
        }
    }
}

/// PairHMM inputs of one read, as GATK's `ReadDataHolder` (bases are the read's bases).
#[derive(Clone, Debug, Default)]
pub struct HmmQuals {
    pub quals: Vec<u8>,
    pub ins: Vec<u8>,
    pub del: Vec<u8>,
}

/// `modifyReadQualities` for one read (modifySoftclippedBases = true, i.e. no clipping):
/// insertion/deletion GOPs start at Q45 (no BI/BD tags), are capped by the PCR indel model
/// (all positions except the last), base qualities are capped at the mapping quality, and then
/// quals < 18 become 6 and ins/del < 6 become 6.
pub fn read_hmm_quals(seq: &[u8], qual: &[u8], mapq: u8, pcr: Option<&[u8; MAX_REPEAT_LENGTH + 1]>, p: &Params, out: &mut HmmQuals) {
    let mut rt = RepeatTables::default();
    let mut rl = Vec::new();
    read_hmm_quals_with(seq, qual, &[], &[], mapq, pcr, p, out, &mut rt, &mut rl)
}

/// [`read_hmm_quals`] with caller-owned scratch buffers.
#[allow(clippy::too_many_arguments)]
/// `bi` / `bd`: the read's BI / BD qualities (empty = no tag, Q45 per base).
pub fn read_hmm_quals_with(seq: &[u8], qual: &[u8], bi: &[u8], bd: &[u8], mapq: u8, pcr: Option<&[u8; MAX_REPEAT_LENGTH + 1]>, p: &Params, out: &mut HmmQuals, rt: &mut RepeatTables, rl: &mut Vec<u8>) {
    let n = seq.len();
    out.quals.clear();
    out.quals.extend_from_slice(qual);
    out.ins.clear();
    if bi.is_empty() {
        out.ins.resize(n, p.default_indel_qual);
    } else {
        out.ins.extend_from_slice(bi);
    }
    out.del.clear();
    if bd.is_empty() {
        out.del.resize(n, p.default_indel_qual);
    } else {
        out.del.extend_from_slice(bd);
    }
    if let Some(t) = pcr {
        rt.repeat_lengths(seq, rl);
        for (o, &len) in rl.iter().enumerate() {
            let c = t[len as usize];
            out.ins[o] = out.ins[o].min(c);
            out.del[o] = out.del[o].min(c);
        }
    }
    let minq = p.min_usable_q_score as i8;
    let bq_thresh = p.base_quality_score_threshold as i8;
    for i in 0..n {
        let mut q = out.quals[i];
        if !p.disable_cap_read_qualities_to_mapq {
            q = q.min(mapq);
        }
        // setToFixedValueIfTooLow compares signed Java bytes
        if (q as i8) < bq_thresh {
            q = p.min_usable_q_score;
        }
        out.quals[i] = q;
        if (out.ins[i] as i8) < minq {
            out.ins[i] = p.min_usable_q_score;
        }
        if (out.del[i] as i8) < minq {
            out.del[i] = p.min_usable_q_score;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fast_tandem_repeat_length_matches_reference() {
        let mut x: u64 = 0x9E3779B97F4A7C15;
        let mut rnd = || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        let mut rt = RepeatTables::default();
        let mut rl = Vec::new();
        for it in 0..200_000 {
            let len = 1 + (rnd() % 60) as usize;
            let alpha = 1 + (rnd() % 4) as usize;
            let mut b: Vec<u8> = (0..len).map(|_| b"ACGT"[(rnd() % alpha as u64) as usize]).collect();
            if it % 3 == 0 {
                // periodic sequences with a random period and some noise
                let per = 1 + (rnd() % 9) as usize;
                for i in per..len {
                    b[i] = b[i - per];
                }
                if rnd() % 2 == 0 {
                    let k = (rnd() % len as u64) as usize;
                    b[k] = b"ACGTN"[(rnd() % 5) as usize];
                }
            }
            for o in 0..len {
                assert_eq!(tandem_repeat_length(&b, o), tandem_repeat_length_reference(&b, o), "{:?} {}", String::from_utf8_lossy(&b), o);
            }
            rt.repeat_lengths(&b, &mut rl);
            assert_eq!(rl.len(), len - 1);
            for o in 0..len - 1 {
                assert_eq!(rl[o] as usize, tandem_repeat_length_reference(&b, o), "table {:?} {}", String::from_utf8_lossy(&b), o);
            }
        }
    }
}
