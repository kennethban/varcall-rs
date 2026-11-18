//! Reference-confidence (`-ERC GVCF`) genotyping: `HaplotypeCallerGenotypingEngine.assignGenotypeLikelihoods`
//! with `emitReferenceConfidence = true`, `GenotypingEngine.calculateGenotypes` with the settings
//! `HaplotypeCallerEngine` forces in that mode (calling confidence -0.0, all-site PLs), the annotation
//! set of reference-confidence mode (raw MQ, StrandBiasBySample; no AC/AF/AN, FS, SOR, QD) and
//! physical phasing (`AssemblyBasedCallerUtils.phaseCalls`).

use hc_assembly::Read;

use crate::annotate;
use crate::engine::{alleles_to_keep_by_hap_scores, max_acceptable_allele_count, reverse_trim_alleles, Call, GenotypingEngine, RegionInput};
use crate::events::{create_allele_mapper, event_start_positions, merged_vc_at, Allele, EventMap};
use crate::genotyper;
use crate::likelihoods::{best_allele, ref_priorities};
use crate::stats;
use crate::vcf::{Phase, VcfRecord};

pub const NON_REF: &[u8] = b"<NON_REF>";

#[inline]
fn is_non_ref(a: &Allele) -> bool {
    !a.is_ref && a.bases == NON_REF
}

/// htsjdk `Allele.length()`: 0 for symbolic alleles.
#[inline]
fn allele_length(a: &Allele) -> usize {
    if is_non_ref(a) {
        0
    } else {
        a.len()
    }
}

/// `Math.min(double, double)` (NaN-propagating)
#[inline]
fn java_min(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.min(b)
    }
}

/// commons-math `Median` (`Percentile(50)`, LEGACY estimation) of `v` (no NaNs).
fn median(v: &mut [f64]) -> f64 {
    let n = v.len();
    if n == 0 {
        return f64::NAN;
    }
    if n == 1 {
        return v[0];
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let pos = 0.5 * (n as f64 + 1.0);
    let fpos = pos.floor();
    let int_pos = fpos as usize;
    let dif = pos - fpos;
    if pos < 1.0 {
        return v[0];
    }
    if pos >= n as f64 {
        return v[n - 1];
    }
    let lower = v[int_pos - 1];
    let upper = v[int_pos];
    lower + dif * (upper - lower)
}

/// `AlleleLikelihoods.updateNonRefAlleleLikelihoods(allelesToConsider)` for one read: `lks` holds all
/// alleles (the non-ref one at `nr`), `consider[a]` whether allele `a` is in `allelesToConsider`.
fn non_ref_likelihood(lks: &[f64], nr: usize, consider: &[bool], buf: &mut Vec<f64>) -> f64 {
    // searchBestAllele(canBeReference = true), no priorities
    let mut best_lk = lks[0];
    for &c in &lks[1..] {
        if c > best_lk {
            best_lk = c;
        }
    }
    buf.clear();
    for (i, &l) in lks.iter().enumerate() {
        if i != nr && l < best_lk && !l.is_nan() && consider[i] {
            buf.push(l);
        }
    }
    let m = median(buf);
    if !m.is_nan() {
        m
    } else if {
        stats::inc(&stats::ERC_NONREF_MEDIAN_NAN);
        lks.len() - 1 <= 1
    } {
        f64::NAN
    } else {
        best_lk
    }
}

/// One reference-confidence-mode call.
#[derive(Clone, Debug)]
pub struct ErcCall {
    pub rec: VcfRecord,
    /// VC end of the (trimmed) call
    pub end: i64,
    /// the genotype is called (not a no-call)
    pub called: bool,
}

struct PendingCall {
    rec: VcfRecord,
    /// haplotype indices supporting the call's alleles (for `calledHaplotypes`)
    haps: Vec<usize>,
}

/// `AssemblyBasedCallerUtils.isSiteSpecificAltAllele`
fn is_site_specific_alt(a: &[u8], is_ref: bool) -> bool {
    !(is_ref || a == NON_REF || a == b"*")
}

impl GenotypingEngine {
    /// `GenotypingEngine.calculateGenotypes` with `-ERC GVCF` settings (alleles include `<NON_REF>` last).
    fn calculate_genotypes_erc(&mut self, loc: i64, alleles: &[Allele], pls: &[i32]) -> Option<Call> {
        let n = alleles.len();
        if n > 50 {
            return None;
        }
        let nr = n - 1; // <NON_REF>
        let reduced: Vec<usize> = if self.cfg.max_alternate_alleles < n - 1 {
            // calculateMostLikelyAlleles: the non-ref allele does not count as a proper alt and is always kept
            let proper = n - 2;
            if self.cfg.max_alternate_alleles >= proper {
                (0..n).collect()
            } else {
                stats::inc(&stats::MAX_ALT_ALLELES_REDUCTION);
                let gls: Vec<f64> = pls.iter().map(|&p| p as f64 / -10.0).collect();
                let mut sums = vec![0f64; n];
                let best = crate::mathutils::max_element_index(&gls);
                let diff = (gls[best] - gls[0]).abs();
                let (i, j) = genotyper::genotype_alleles(n)[best];
                for a in 1..n {
                    if a == i || a == j {
                        sums[a] += diff;
                    }
                }
                let mut idx: Vec<usize> = (1..n).filter(|&a| a != nr).collect();
                idx.sort_by(|&a, &b| sums[b].partial_cmp(&sums[a]).unwrap());
                let kept: Vec<usize> = idx.into_iter().take(self.cfg.max_alternate_alleles).collect();
                (0..n).filter(|&a| a == 0 || a == nr || kept.contains(&a)).collect()
            }
        } else {
            (0..n).collect()
        };
        // the reduced genotype is re-derived through subsetAlleles whenever the reduction branch is taken
        let reduced_pls: Vec<i32> =
            if self.cfg.max_alternate_alleles < n - 1 { genotyper::subset_and_call(pls, &reduced).pl } else { pls.to_vec() };
        let lens: Vec<usize> = reduced.iter().map(|&a| allele_length(&alleles[a])).collect();
        let sd = reduced.iter().position(|&a| alleles[a].is_span_del());
        let af = self.af.calculate(&reduced_pls, &lens, sd);

        // qualToErrorProbLog10(-0.0)
        let threshold = self.cfg.standard_confidence_for_calling * -0.1;
        let mut out_alts: Vec<usize> = Vec::new(); // positions in `reduced`
        let mut monomorphic = true;
        for k in 1..reduced.len() {
            let a = &alleles[reduced[k]];
            let plausible = af.log10_p_ref_by_allele[k] + 1.0e-10 < threshold;
            let spurious = a.is_span_del() && !self.is_covered_by_deletion(loc);
            if spurious {
                stats::inc(&stats::SPAN_DEL_SPURIOUS);
            }
            // forceKeepAllele is true for every allele in reference-confidence mode
            if !spurious {
                out_alts.push(k);
            }
            monomorphic &= !(plausible && !spurious);
        }
        let log10_conf = af.log10_p_no_variant + 0.0;
        let phred = (-10.0 * log10_conf) + 0.0;
        let passes = !monomorphic && phred >= self.cfg.standard_confidence_for_calling;
        let first_is_non_ref = out_alts.first().map_or(false, |&k| reduced[k] == nr);
        if !passes && !first_is_non_ref {
            return None;
        }
        if !passes {
            stats::inc(&stats::ERC_MONOMORPHIC_EMITTED_NONREF_FIRST);
        }
        if out_alts.len() == 1 && alleles[reduced[out_alts[0]]].is_span_del() {
            stats::inc(&stats::ONLY_SPAN_DEL_DROPPED);
            return None;
        }
        let mut out_alleles = vec![0usize];
        out_alleles.extend(out_alts.iter().map(|&k| reduced[k]));
        // recordDeletions: <NON_REF> has length 0, so it records a "deletion" of the whole reference allele
        self.upstream_deletions.retain(|&(_, e)| e >= loc);
        let ref_len = alleles[0].len() as i64;
        for &a in &out_alleles {
            let del = ref_len - allele_length(&alleles[a]) as i64;
            if del > 0 {
                self.upstream_deletions.push((loc, loc + del));
            }
        }
        let genotype = genotyper::subset_and_call(pls, &out_alleles);
        let mleac = out_alts.iter().map(|&k| af.mle_alt_counts[k - 1]).collect();
        Some(Call { alleles: out_alleles, log10_p_error: log10_conf, genotype, others: Vec::new(), mleac })
    }

    /// Genotypes one region in reference-confidence mode. Returns the (phased) calls in order.
    pub fn genotype_region_erc(&mut self, input: &RegionInput, realigned: &[Read]) -> Vec<ErcCall> {
        let haps = input.haplotypes;
        let maps: Vec<EventMap> = haps
            .iter()
            .map(|h| EventMap::from_haplotype(h, input.ref_padded, input.ref_padded_start, self.cfg.max_mnp_distance))
            .collect();
        let starts = event_start_positions(&maps);
        let mut pending: Vec<PendingCall> = Vec::new();
        let n_reads = realigned.len();
        let mut buf: Vec<f64> = Vec::new();
        for loc in starts {
            if loc < input.span_start || loc > input.span_end {
                continue;
            }
            let mut merged = match merged_vc_at(loc, &maps, input.ref_padded, input.ref_padded_start) {
                Some(m) => m,
                None => continue,
            };
            let n_alleles_before_trim = merged.alleles.len() + 1;
            stats::inc(&stats::LOCI_EVALUATED);
            let (mapper_alleles, mut mapper_haps) = create_allele_mapper(&merged, loc, &maps);
            assert!(mapper_alleles == merged.alleles, "allele mapper differs from merged VC alleles at {}", loc);
            let practical = max_acceptable_allele_count(self.cfg.max_genotype_count);
            if mapper_alleles.len() > practical {
                stats::inc(&stats::TOO_MANY_GENOTYPES_REDUCTION);
                let keep = alleles_to_keep_by_hap_scores(&mapper_alleles, &mapper_haps, haps, practical);
                merged.alleles = keep.iter().map(|&a| mapper_alleles[a].clone()).collect();
                mapper_haps = keep.iter().map(|&a| mapper_haps[a].clone()).collect();
            }
            let vc_end = merged.end();
            // addNonRefSymbolicAllele
            let mut alleles = merged.alleles.clone();
            alleles.push(Allele::new(NON_REF, false));
            let n = alleles.len();
            let nr = n - 1;
            let lo = (loc - self.cfg.informative_read_overlap_margin).max(1);
            let hi = (vc_end + self.cfg.informative_read_overlap_margin).min(self.cfg.contig_length);
            let retained: Vec<usize> = (0..n_reads).filter(|&r| realigned[r].start <= hi && lo <= realigned[r].end).collect();
            // marginalize + addNonReferenceAllele (median of the non-best concrete alleles)
            let mut lk: Vec<Vec<f64>> = vec![vec![f64::NEG_INFINITY; retained.len()]; n];
            for a in 0..nr {
                for (k, &r) in retained.iter().enumerate() {
                    let mut m = f64::NEG_INFINITY;
                    for &h in &mapper_haps[a] {
                        m = m.max(input.matrix[h][r]);
                    }
                    lk[a][k] = m;
                }
            }
            let all = vec![true; n];
            let mut col = vec![0f64; n];
            for k in 0..retained.len() {
                for a in 0..n {
                    col[a] = lk[a][k];
                }
                lk[nr][k] = non_ref_likelihood(&col, nr, &all, &mut buf);
            }
            // contamination downsampling precedes addNonReferenceAllele: best alleles among the event's alleles
            let gls = match self.contamination_removals(&lk[..nr], retained.len()) {
                Some(rm) => genotyper::log10_genotype_likelihoods(&crate::contamination::keep_columns(&lk, &rm), retained.len() - rm.len()),
                None => genotyper::log10_genotype_likelihoods(&lk, retained.len()),
            };
            let pls = genotyper::gls_to_pls(&gls);
            let call = match self.calculate_genotypes_erc(loc, &alleles, &pls) {
                Some(c) => c,
                None => continue,
            };
            // prepareReadAlleleLikelihoodsForAnnotation: non-ref likelihoods re-derived from the call's alleles
            if call.alleles.len() != n {
                stats::inc(&stats::ERC_NONREF_RECOMPUTED_FOR_CALL_ALLELES);
                let mut consider = vec![false; n];
                for &a in &call.alleles {
                    consider[a] = true;
                }
                for k in 0..retained.len() {
                    for a in 0..n {
                        col[a] = lk[a][k];
                    }
                    lk[nr][k] = non_ref_likelihood(&col, nr, &consider, &mut buf);
                }
            }
            let filt: Vec<&Read> = input.filtered_reads.iter().filter(|r| r.start <= hi && lo <= r.end).collect();
            let mut ev: Vec<&Read> = retained.iter().map(|&r| &realigned[r]).collect();
            ev.extend(filt.iter().copied());
            let n_ret = retained.len();
            let lk_col = |a: usize, e: usize| -> f64 { if e < n_ret { lk[a][e] } else { 0.0 } };
            stats::inc(&stats::CALLS);
            match call.genotype.alleles {
                None => stats::inc(&stats::ERC_NO_CALL),
                Some((0, 0)) => stats::inc(&stats::ERC_HOM_REF_CALL),
                _ => {}
            }
            let rec = annotate_call_erc(&self.cfg.contig, loc, vc_end, &alleles, &call, &ev, &lk_col, n_alleles_before_trim, self.cfg.as_annotations);
            let mut call_haps: Vec<usize> = Vec::new();
            for &a in &call.alleles {
                if a < mapper_haps.len() {
                    call_haps.extend(mapper_haps[a].iter().copied());
                }
            }
            pending.push(PendingCall { rec, haps: call_haps });
        }
        // physical phasing
        let mut called: Vec<bool> = vec![false; haps.len()];
        for p in &pending {
            for &h in &p.haps {
                called[h] = true;
            }
        }
        let phases = phase_calls(&pending, &called, &maps);
        pending
            .into_iter()
            .zip(phases)
            .map(|(p, ph)| {
                let mut rec = p.rec;
                rec.phase = ph;
                let end = rec.pos + rec.alleles[0].len() as i64 - 1;
                let called = rec.gt.is_some();
                ErcCall { rec, end, called }
            })
            .collect()
    }
}

#[allow(clippy::too_many_arguments)]
fn annotate_call_erc(
    contig: &str,
    loc: i64,
    vc_end: i64,
    merged: &[Allele],
    call: &Call,
    ev: &[&Read],
    lk: &dyn Fn(usize, usize) -> f64,
    n_alleles_before_trim: usize,
    as_annotations: bool,
) -> VcfRecord {
    let n = merged.len();
    let n_ev = ev.len();
    let call_alleles = &call.alleles;
    let in_call = |a: usize| call_alleles.contains(&a);
    let g = &call.genotype;

    // genotype annotations (only for called genotypes): AD, DP, SB
    let mut ad: Option<Vec<i32>> = None;
    let mut dp: Option<i32> = None;
    let mut sb: Option<[i32; 4]> = None;
    // best allele of every evidence over all alleles (ref priority)
    let pr = ref_priorities(n, Some(0));
    let mut col = vec![0f64; n];
    let mut best: Vec<Option<usize>> = Vec::with_capacity(n_ev);
    for e in 0..n_ev {
        for (a, c) in col.iter_mut().enumerate() {
            *c = lk(a, e);
        }
        let b = best_allele(&col, Some(&pr));
        best.push(if b.is_informative() { b.allele } else { None });
    }
    if g.alleles.is_some() {
        let prc = ref_priorities(call_alleles.len(), Some(0));
        let mut counts = vec![0i32; call_alleles.len()];
        let mut colc = vec![0f64; call_alleles.len()];
        let mut informative = 0;
        for e in 0..n_ev {
            for (k, &a) in call_alleles.iter().enumerate() {
                colc[k] = lk(a, e);
            }
            let b = best_allele(&colc, Some(&prc));
            if b.is_informative() {
                counts[b.allele.unwrap()] += 1;
                informative += 1;
            }
        }
        ad = Some(counts);
        dp = Some(if n_ev == 0 { 0 } else { informative });
        let is_ref_a = |a: usize| a == 0;
        let is_alt_a = |a: usize| a != 0 && in_call(a);
        let t = annotate::strand_table(&best, ev, &is_ref_a, &is_alt_a, 0);
        sb = Some([t[0][0], t[0][1], t[1][0], t[1][1]]);
    }

    let mut info: Vec<(&'static str, String)> = Vec::new();
    // rank sum tests
    let (mut bq_ref, mut bq_alt, mut mq_ref, mut mq_alt, mut rp_ref, mut rp_alt) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for (e, b) in best.iter().enumerate() {
        let a = match b {
            Some(a) => *a,
            None => continue,
        };
        let r = ev[e];
        if !(r.mapq != 0 && r.mapq != 255) {
            continue;
        }
        let is_ref = a == 0;
        let is_alt = !is_ref && in_call(a);
        if !(is_ref || is_alt) {
            continue;
        }
        if let Some(v) = annotate::base_qual_value(r, loc) {
            if is_ref { bq_ref.push(v) } else { bq_alt.push(v) }
        }
        if is_ref { mq_ref.push(r.mapq as f64) } else { mq_alt.push(r.mapq as f64) }
        if let Some(v) = annotate::read_pos_value(r, loc, vc_end) {
            if is_ref { rp_ref.push(v) } else { rp_alt.push(v) }
        }
    }
    if let Some(s) = annotate::rank_sum(&bq_alt, &bq_ref) {
        info.push(("BaseQRankSum", s));
    }
    if n_ev > 0 {
        info.push(("DP", n_ev.to_string()));
    }
    if g.alleles.is_some() {
        info.push(("ExcessHet", annotate::excess_het(&g.pl, call_alleles.len())));
    }
    let an_called = if g.alleles.is_some() { 2.0 } else { 0.0 };
    info.push(("MLEAC", call.mleac.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(",")));
    info.push((
        "MLEAF",
        call.mleac
            .iter()
            .map(|&c| crate::javamath::format_vcf_double(java_min(1.0, c as f64 / an_called)))
            .collect::<Vec<_>>()
            .join(","),
    ));
    if let Some(s) = annotate::rank_sum(&mq_alt, &mq_ref) {
        info.push(("MQRankSum", s));
    }
    // RMSMappingQuality.annotateRawData
    if n_ev > 0 {
        let (mut ss, mut nmq) = (0i64, 0i64);
        for r in ev {
            let mq = r.mapq as i64;
            if mq != 255 {
                ss += mq * mq;
                nmq += 1;
            }
        }
        info.push(("RAW_MQandDP", format!("{},{}", ss, nmq)));
    }
    if let Some(s) = annotate::rank_sum(&rp_alt, &rp_ref) {
        info.push(("ReadPosRankSum", s));
    }
    if as_annotations {
        // AS_StandardAnnotation raw data (`annotateRawData`), per allele of the call (with <NON_REF>)
        let na = call_alleles.len();
        let pos: Vec<Option<usize>> = best.iter().map(|b| b.and_then(|a| call_alleles.iter().position(|&x| x == a))).collect();
        let (mut bq, mut mq, mut rp) = (vec![Vec::new(); na], vec![Vec::new(); na], vec![Vec::new(); na]);
        let mut sum_sq: Vec<Option<f64>> = vec![None; na];
        for (e, p) in pos.iter().enumerate() {
            let Some(k) = *p else { continue };
            let r = ev[e];
            if r.mapq != 255 {
                sum_sq[k] = Some(sum_sq[k].unwrap_or(0.0) + (r.mapq as i32 * r.mapq as i32) as f64);
            }
            if !(r.mapq != 0 && r.mapq != 255) {
                continue;
            }
            if let Some(v) = annotate::base_qual_value(r, loc) {
                bq[k].push(v as i32 as f64);
            }
            mq[k].push(r.mapq as i32 as f64);
            if let Some(v) = annotate::read_pos_value(r, loc, vc_end) {
                rp[k].push(v as i32 as f64);
            }
        }
        info.push(("AS_RAW_BaseQRankSum", annotate::as_raw_rank_sum(&bq)));
        info.push(("AS_RAW_MQ", annotate::as_raw_mq(&sum_sq)));
        info.push(("AS_RAW_MQRankSum", annotate::as_raw_rank_sum(&mq)));
        info.push(("AS_RAW_ReadPosRankSum", annotate::as_raw_rank_sum(&rp)));
        info.push(("AS_SB_TABLE", annotate::as_sb_table(&pos, ev, na)));
    }
    info.sort_by(|a, b| a.0.cmp(b.0));

    let mut out_alleles: Vec<Vec<u8>> = call_alleles.iter().map(|&a| merged[a].bases.clone()).collect();
    if call_alleles.len() != n_alleles_before_trim {
        let (_, trimmed) = reverse_trim_alleles(loc, &out_alleles);
        out_alleles = trimmed;
    }
    VcfRecord {
        contig: contig.to_string(),
        pos: loc,
        id: None,
        alleles: out_alleles,
        // CommonInfo.getPhredScaledQual: (log10PError * -10) + 0.0
        qual: (call.log10_p_error * -10.0) + 0.0,
        info,
        gt: g.alleles,
        ad,
        dp,
        gq: g.gq,
        pl: g.pl.clone(),
        sb,
        phase: None,
        more: Vec::new(),
    }
}

/// `AssemblyBasedCallerUtils.phaseCalls` (constructHaplotypeMapping, constructPhaseSetMapping, constructPhaseGroups).
fn phase_calls(calls: &[PendingCall], called_haps: &[bool], maps: &[EventMap]) -> Vec<Option<Phase>> {
    let nc = calls.len();
    let nh = called_haps.len();
    // haplotypes (as bitsets over haplotype indices) carrying each call's single site-specific alt allele
    let hap_sets: Vec<Vec<bool>> = calls
        .iter()
        .map(|c| {
            let r = &c.rec;
            let ss: Vec<&Vec<u8>> = r.alleles[1..].iter().filter(|a| is_site_specific_alt(a, false)).collect();
            if ss.len() != 1 {
                return vec![false; nh];
            }
            let alt = ss[0];
            (0..nh)
                .map(|h| called_haps[h] && maps[h].events.values().any(|e| e.start == r.pos && &e.alt_bases == alt))
                .collect()
        })
        .collect();
    let count = |s: &[bool]| s.iter().filter(|&&x| x).count();
    let contains_all = |a: &[bool], b: &[bool]| a.iter().zip(b).all(|(&x, &y)| x || !y);
    let mut union = vec![false; nh];
    for s in &hap_sets {
        for h in 0..nh {
            union[h] |= s[h];
        }
    }
    let total = count(&union);
    // phaseSetMapping: call -> (group id, PHASE_01 = true / PHASE_10 = false)
    let mut mapping: Vec<Option<(usize, bool)>> = vec![None; nc];
    let mut unique = 0usize;
    let mut aborted = false;
    'outer: for i in 0..nc.saturating_sub(1) {
        let hw = &hap_sets[i];
        let nw = count(hw);
        if nw == 0 {
            continue;
        }
        let call_on_all = nw == total;
        let mut available = hw.clone();
        for j in i + 1..nc {
            let hc = &hap_sets[j];
            let ncomp = count(hc);
            if ncomp == 0 {
                continue;
            }
            let comp_on_all = ncomp == total;
            if (nw == ncomp && contains_all(hw, hc)) || (call_on_all && contains_all(&available, hc)) || comp_on_all {
                if mapping[i].is_none() {
                    if mapping[j].is_some() {
                        aborted = true;
                        break 'outer;
                    }
                    mapping[i] = Some((unique, true));
                    mapping[j] = Some((unique, true));
                    for h in 0..nh {
                        available[h] &= hc[h];
                    }
                    unique += 1;
                } else if mapping[j].is_none() {
                    mapping[j] = mapping[i];
                }
            } else if nw + ncomp == total {
                let disjoint = (0..nh).all(|h| !(hw[h] && hc[h]));
                if disjoint {
                    if mapping[i].is_none() {
                        if mapping[j].is_some() {
                            aborted = true;
                            break 'outer;
                        }
                        mapping[i] = Some((unique, true));
                        mapping[j] = Some((unique, false));
                        unique += 1;
                    } else if mapping[j].is_none() {
                        let (g, p01) = mapping[i].unwrap();
                        mapping[j] = Some((g, !p01));
                    }
                }
            }
        }
    }
    let mut out: Vec<Option<Phase>> = vec![None; nc];
    if aborted {
        stats::inc(&stats::ERC_PHASING_ABORTED);
        return out;
    }
    for group in 0..unique {
        let idx: Vec<usize> = (0..nc).filter(|&k| mapping[k].map_or(false, |(g, _)| g == group)).collect();
        assert!(idx.len() >= 2, "phase group with fewer than 2 members");
        let first = &calls[idx[0]].rec;
        let pid = format!(
            "{}_{}_{}",
            first.pos,
            String::from_utf8_lossy(&first.alleles[0]),
            String::from_utf8_lossy(&first.alleles[1])
        );
        let ps = first.pos;
        for &k in &idx {
            let p01 = mapping[k].unwrap().1;
            let r = &calls[k].rec;
            // phaseVC: a het genotype is reversed when the allele at the phased index is not site-specific
            let reversed = match r.gt {
                Some((i, j)) if i != j => {
                    let at = if p01 { j } else { i };
                    !is_site_specific_alt(&r.alleles[at], at == 0)
                }
                _ => false,
            };
            stats::inc(&stats::ERC_PHASED_CALLS);
            if reversed {
                stats::inc(&stats::ERC_PHASED_GT_REVERSED);
            }
            out[k] = Some(Phase { pgt: if p01 { "0|1" } else { "1|0" }, pid: pid.clone(), ps, reversed });
        }
    }
    out
}
