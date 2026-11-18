//! Active-region detection: HC read filters, per-locus ref-vs-any likelihoods (isActive), the band-pass
//! activity profile and region cutting. Ported 1:1 from GATK (see docs/notes/variant-rs.md).

use crate::bam::Rec;
use std::collections::VecDeque;

// ------------------------------------------------------------------------------------------
// Constants (GATK defaults)
// ------------------------------------------------------------------------------------------
pub const MIN_MAPQ: u8 = 20; // HaplotypeCallerEngine.DEFAULT_READ_QUALITY_FILTER_THRESHOLD
pub const MIN_BASE_QUAL: u8 = 10; // AssemblyBasedCallerArgumentCollection.minBaseQualityScore
pub const REF_MODEL_DELETION_QUAL: u8 = 30; // ReferenceConfidenceModel.REF_MODEL_DELETION_QUAL
pub const HQ_SOFTCLIP_QUAL: u8 = 28; // HQ_BASE_QUALITY_SOFTCLIP_THRESHOLD (strictly greater than)
pub const AVG_HQ_SOFTCLIPS_THRESHOLD: f64 = 6.0; // AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD
pub const ACTIVE_PROB_THRESHOLD: f64 = 0.002;
pub const MAX_PROB_PROPAGATION: i64 = 50;
pub const MIN_REGION: usize = 50;
pub const MAX_REGION: usize = 300;

// CIGAR ops (BAM encoding)
pub const OP_M: u8 = 0;
pub const OP_I: u8 = 1;
pub const OP_D: u8 = 2;
pub const OP_N: u8 = 3;
pub const OP_S: u8 = 4;
pub const OP_H: u8 = 5;
pub const OP_EQ: u8 = 7;
pub const OP_X: u8 = 8;

pub fn on_genome(op: u8) -> bool { matches!(op, OP_M | OP_EQ | OP_X | OP_D) }
pub fn consumes_ref(op: u8) -> bool { matches!(op, OP_M | OP_D | OP_N | OP_EQ | OP_X) }
pub fn consumes_read(op: u8) -> bool { matches!(op, OP_M | OP_I | OP_S | OP_EQ | OP_X) }

// ------------------------------------------------------------------------------------------
// GATK numeric tables (exact)
// ------------------------------------------------------------------------------------------
pub struct Tables {
    pub prior: [f64; 3],
    pub filter_size: i64,
    pub kernel: Vec<f64>,
    pub qual_to_prob_log10: [f64; 256],
    pub jacobian: Vec<f64>,
    pub log10_one_third: f64,
    pub log10_ploidy: f64, // log10(2)
    pub log10_one: f64,    // log10(1) = 0
    pub inv_step: f64,
    /// `add_element`'s additions to gl[0], gl[1], gl[2] for [is_alt][qual] (computed once with the same expressions,
    /// so the sums are bit-identical)
    pub elem: Box<[[[f64; 3]; 256]; 2]>,
}

impl Tables {
    pub fn load() -> Self {
        let b = include_bytes!("../tables/gatk_tables.bin");
        let mut it = b.chunks_exact(8).map(|c| f64::from_le_bytes(c.try_into().unwrap()));
        let mut next = || it.next().unwrap();
        let prior = [next(), next(), next()];
        let filter_size = next() as i64;
        let kernel: Vec<f64> = (0..2 * filter_size + 1).map(|_| next()).collect();
        let mut q = [0f64; 256];
        for v in q.iter_mut() { *v = next(); }
        let jsize = next() as usize;
        let jacobian: Vec<f64> = (0..jsize).map(|_| next()).collect();
        let mut t = Tables { prior, filter_size, kernel, qual_to_prob_log10: q, jacobian,
                 log10_one_third: next(), log10_ploidy: next(), log10_one: next(), inv_step: next(),
                 elem: Box::new([[[0.0; 3]; 256]; 2]) };
        let mut elem = Box::new([[[0.0; 3]; 256]; 2]);
        for (a, row) in elem.iter_mut().enumerate() {
            for (qual, e) in row.iter_mut().enumerate() {
                *e = t.element_terms(a == 1, qual as u8);
            }
        }
        t.elem = elem;
        t
    }

    /// ReferenceConfidenceModel.applyPileupElementRefVsNonRefLikelihoodAndCount: the terms added to gl[0..3]
    fn element_terms(&self, is_alt: bool, qual: u8) -> [f64; 3] {
        let (ref_lk, nonref_lk) = if is_alt {
            (qual_to_error_prob_log10(qual) + self.log10_one_third, self.qual_to_prob_log10[qual as usize])
        } else {
            (self.qual_to_prob_log10[qual as usize], qual_to_error_prob_log10(qual) + self.log10_one_third)
        };
        [
            1.0 * (ref_lk + self.log10_ploidy),
            1.0 * self.approximate_log10_sum_log10(ref_lk + self.log10_one, nonref_lk + self.log10_one),
            1.0 * (nonref_lk + self.log10_ploidy),
        ]
    }

    #[inline]
    pub fn approximate_log10_sum_log10(&self, a: f64, b: f64) -> f64 {
        let (a, b) = if a > b { (b, a) } else { (a, b) };
        if a == f64::NEG_INFINITY { return b; }
        let diff = b - a;
        b + if diff < 8.0 { self.jacobian[fast_round(diff * self.inv_step) as usize] } else { 0.0 }
    }
}

#[inline]
pub fn fast_round(d: f64) -> i32 { if d > 0.0 { (d + 0.5) as i32 } else { (d - 0.5) as i32 } }

#[inline]
pub fn qual_to_error_prob_log10(q: u8) -> f64 { (q as f64) * -0.1 }

pub fn max_element_index(a: &[f64]) -> usize {
    let mut m = 0;
    for i in 1..a.len() { if a[i] > a[m] { m = i; } }
    m
}

/// MathUtils.log10SumLog10(double[])
pub fn log10_sum_log10(v: &[f64]) -> f64 {
    let mi = max_element_index(v);
    let max = v[mi];
    if max == f64::NEG_INFINITY { return max; }
    let mut sum = 1.0;
    for (i, &x) in v.iter().enumerate() {
        if i == mi || x == f64::NEG_INFINITY { continue; }
        sum += crate::jmath::pow(10.0, x - max);
    }
    max + if sum != 1.0 { crate::jmath::log10(sum) } else { 0.0 }
}
// ------------------------------------------------------------------------------------------
// HaplotypeCaller read filters
// ------------------------------------------------------------------------------------------
#[derive(Default, Debug)]
pub struct FilterCounts { pub total: u64, pub passed: u64, pub mapq: u64, pub unmapped: u64, pub secondary: u64, pub duplicate: u64, pub qcfail: u64, pub zero_ref_len: u64, pub bad_cigar: u64, pub malformed: u64 }

impl FilterCounts {
    pub fn add(&mut self, o: &FilterCounts) {
        self.total += o.total; self.passed += o.passed; self.mapq += o.mapq; self.unmapped += o.unmapped; self.secondary += o.secondary;
        self.duplicate += o.duplicate; self.qcfail += o.qcfail; self.zero_ref_len += o.zero_ref_len; self.bad_cigar += o.bad_cigar; self.malformed += o.malformed;
    }
}

pub fn ref_length(c: &[(u8, u32)]) -> i64 { c.iter().filter(|e| consumes_ref(e.0)).map(|e| e.1 as i64).sum() }
pub fn read_length(c: &[(u8, u32)]) -> usize { c.iter().filter(|e| consumes_read(e.0)).map(|e| e.1 as usize).sum() }

/// CigarUtils.isGood (+ the structural parts of htsjdk Cigar.isValid)
pub fn good_cigar(c: &[(u8, u32)]) -> bool {
    if c.is_empty() { return true; }
    let n = c.len();
    for (i, &(op, len)) in c.iter().enumerate() {
        if len == 0 { return false; }
        if op == OP_H && i != 0 && i != n - 1 { return false; }
        if op == OP_S {
            // soft clips may only be preceded/followed by hard clips or the read end
            let left_ok = c[..i].iter().all(|e| e.0 == OP_H);
            let right_ok = c[i + 1..].iter().all(|e| e.0 == OP_H);
            if !left_ok && !right_ok { return false; }
        }
    }
    // consecutive indels
    let mut prev_indel = false;
    for &(op, _) in c {
        let ind = op == OP_I || op == OP_D;
        if prev_indel && ind { return false; }
        prev_indel = ind;
    }
    // starts or ends with deletion, ignoring clips
    for e in c.iter() { if e.0 == OP_D { return false; } if !(e.0 == OP_S || e.0 == OP_H) { break; } }
    for e in c.iter().rev() { if e.0 == OP_D { return false; } if !(e.0 == OP_S || e.0 == OP_H) { break; } }
    true
}

pub fn passes_filters(r: &Rec, ref_lens: &[i64], fc: &mut FilterCounts) -> bool {
    fc.total += 1;
    // order follows makeStandardHCReadFilters(); counts attribute each read to the first failing filter
    if r.mapq < MIN_MAPQ || r.mapq == 255 {
        // MappingQualityReadFilter / MappingQualityAvailable (unmapped reads usually also have mapq 0)
        if r.flag & 0x4 != 0 || r.ref_id < 0 { fc.unmapped += 1; } else { fc.mapq += 1; }
        return false;
    }
    if r.flag & 0x4 != 0 || r.ref_id < 0 || r.start < 1 { fc.unmapped += 1; return false; }
    if r.flag & 0x100 != 0 { fc.secondary += 1; return false; }
    if r.flag & 0x400 != 0 { fc.duplicate += 1; return false; }
    if r.flag & 0x200 != 0 { fc.qcfail += 1; return false; }
    if ref_length(&r.cigar) == 0 { fc.zero_ref_len += 1; return false; }
    if !good_cigar(&r.cigar) { fc.bad_cigar += 1; return false; }
    // WellformedReadFilter
    let wellformed = r.has_rg
        && !r.seq.is_empty()
        && !r.qual_missing
        && r.qual.len() == r.seq.len()
        && read_length(&r.cigar) == r.seq.len()
        && !r.cigar.iter().any(|e| e.0 == OP_N)
        && (r.ref_id as usize) < ref_lens.len()
        && r.start <= ref_lens[r.ref_id as usize];
    if !wellformed { fc.malformed += 1; return false; }
    fc.passed += 1;
    true
}

// ------------------------------------------------------------------------------------------
// Per-locus accumulators (the data-oriented replacement for pileups)
// ------------------------------------------------------------------------------------------
#[derive(Clone, Copy, Default)]
pub struct Acc {
    pub gl: [f64; 3],
    pub read_count: u32,
    pub n_elems: u32,      // pileup size (incl. low-qual bases)
    pub hq_n: u32,         // RunningAverage of HQ soft clips
    pub hq_mean: f64,
}

/// Contribution of one pileup element (ReferenceConfidenceModel.applyPileupElementRefVsNonRefLikelihoodAndCount)
#[inline]
pub fn add_element(acc: &mut Acc, t: &Tables, is_alt: bool, qual: u8, next_to_softclip: bool, hq_softclips: u32) {
    let e = &t.elem[is_alt as usize][qual as usize];
    acc.gl[0] += e[0];
    acc.gl[1] += e[1];
    acc.gl[2] += e[2];
    if is_alt && next_to_softclip {
        acc.hq_n += 1;
        acc.hq_mean += (hq_softclips as f64 - acc.hq_mean) / acc.hq_n as f64;
    }
}

/// Neighbour flags of a pileup element, which only depend on the cigar element index and on
/// whether the element is at the first / last base of its cigar element.
pub struct ElemFlags { pub alt_start: bool, pub alt_end: bool, pub alt_both: bool, pub sc_start: bool, pub sc_end: bool, pub sc_both: bool }

pub fn element_flags(c: &[(u8, u32)], i: usize) -> ElemFlags {
    let op = c[i].0;
    let is_del = op == OP_D;
    let adj_prev = if i > 0 { Some(c[i - 1].0) } else { None };
    let adj_next = if i + 1 < c.len() { Some(c[i + 1].0) } else { None };
    let near_prev = c[..i].iter().rev().find(|e| on_genome(e.0)).map(|e| e.0);
    let near_next = c[i + 1..].iter().find(|e| on_genome(e.0)).map(|e| e.0);
    // flags that apply when at start / at end of the cigar element
    let start = (!is_del && near_prev == Some(OP_D)) || adj_prev == Some(OP_I) || adj_prev == Some(OP_S);
    let end = (!is_del && near_next == Some(OP_D)) || adj_next == Some(OP_I) || adj_next == Some(OP_S);
    let sc_start = adj_prev == Some(OP_S);
    let sc_end = adj_next == Some(OP_S);
    ElemFlags { alt_start: start, alt_end: end, alt_both: start || end, sc_start, sc_end, sc_both: sc_start || sc_end }
}

// ------------------------------------------------------------------------------------------
// Activity profile (BandPassActivityProfile + ActivityProfile), ported 1:1
// ------------------------------------------------------------------------------------------
pub struct Profile {
    pub probs: VecDeque<f64>,
    pub start: i64,
    pub stop: i64,
    pub contig_len: i64,
    pub kernel: Vec<f64>,
    pub fsize: i64,
}

#[derive(Clone, Copy)]
pub struct Region { pub start: i64, pub end: i64, pub active: bool }

impl Profile {
    pub fn is_empty(&self) -> bool { self.probs.is_empty() }
    pub fn max_prob_propagation(&self) -> i64 { MAX_PROB_PROPAGATION + self.fsize }
    pub fn loc_for(&self, pos: i64, off: i64) -> Option<i64> {
        let s = pos + off;
        if s < 1 || s > self.contig_len { None } else { Some(s) }
    }
    pub fn incorporate(&mut self, loc: i64, p: f64) {
        let position = loc - self.start;
        assert!(position <= self.probs.len() as i64, "non-contiguous state");
        if position >= 0 {
            let position = position as usize;
            if position < self.probs.len() { self.probs[position] += p; } else { self.probs.push_back(p); }
        }
    }
    pub fn add(&mut self, pos: i64, prob: f64, hq_clips: Option<f64>) {
        if self.probs.is_empty() { self.start = pos; self.stop = pos; } else { assert_eq!(self.stop, pos - 1); self.stop = pos; }
        // ActivityProfile.processState: spread HQ soft-clip states
        let mut supers: Vec<f64> = Vec::new(); // only the probability matters (see quirk below)
        match hq_clips {
            Some(v) => {
                let n = std::cmp::min(v as i32 as i64, self.max_prob_propagation());
                for i in -n..=n { if self.loc_for(pos, i).is_some() { supers.push(prob); } }
            }
            None => supers.push(prob),
        }
        // BandPassActivityProfile.processState. NB: GATK centres every smeared super-state on the
        // *just added* locus (not the super-state's own locus); we replicate that exactly.
        for sp in supers {
            if sp > 0.0 {
                for i in -self.fsize..=self.fsize {
                    if let Some(l) = self.loc_for(pos, i) {
                        let np = sp * self.kernel[(i + self.fsize) as usize];
                        self.incorporate(l, np);
                    }
                }
            } else {
                self.incorporate(pos, prob);
            }
        }
    }
    pub fn is_minimum(&self, i: usize) -> bool {
        if i == self.probs.len() - 1 || i < 1 { return false; }
        let p = self.probs[i];
        p <= self.probs[i + 1] && p < self.probs[i - 1]
    }
    pub fn pop_ready(&mut self, force: bool, out: &mut Vec<Region>) {
        loop {
            if self.probs.is_empty() { return; }
            if force {
                let span = (self.stop - self.start + 1) as usize;
                self.probs.truncate(span);
            }
            let first = self.probs[0];
            let is_active = first > ACTIVE_PROB_THRESHOLD;
            if !force && (self.probs.len() as i64) < MAX_REGION as i64 + self.max_prob_propagation() { return; }
            let n = self.probs.len();
            let mut end = 0;
            while end < n && end < MAX_REGION {
                if (self.probs[end] > ACTIVE_PROB_THRESHOLD) != is_active { break; }
                end += 1;
            }
            if is_active && end == MAX_REGION {
                let mut min_i = end - 1;
                let mut min_p = f64::MAX;
                let mut i = end - 1;
                loop {
                    let cur = self.probs[i];
                    if cur < min_p && self.is_minimum(i) { min_p = cur; min_i = i; }
                    if i == MIN_REGION - 1 { break; }
                    i -= 1;
                }
                end = min_i + 1;
            }
            let offset = end - 1;
            out.push(Region { start: self.start, end: self.start + offset as i64, active: is_active });
            self.probs.drain(..=offset);
            if !self.probs.is_empty() { self.start += offset as i64 + 1; }
        }
    }
}
