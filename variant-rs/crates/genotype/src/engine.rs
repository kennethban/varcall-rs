//! `HaplotypeCallerGenotypingEngine.assignGenotypeLikelihoods` + `GenotypingEngine.calculateGenotypes`
//! + annotation + `RMSMappingQuality.finalizeRawMQ`, producing VCF records.

use hc_assembly::sw::SwAligner;
use hc_assembly::{Haplotype, Read};

use crate::annotate;
use crate::events::{create_allele_mapper, event_start_positions, merged_vc_at, Allele, EventMap};
use crate::genotyper::{self, AfCalculator};
use crate::javamath::{format_fixed, JavaRandom, GATK_RANDOM_SEED};
use crate::likelihoods::{best_allele, ref_priorities};
use crate::realign;
use crate::stats;
use crate::vcf::VcfRecord;

/// Default arguments that matter here (all checked against HaplotypeCallerArgumentCollection).
#[derive(Clone, Debug)]
pub struct EngineConfig {
    pub contig: String,
    /// contig length, for `SimpleInterval.expandWithinContig`
    pub contig_length: i64,
    /// `--max-mnp-distance` (0)
    pub max_mnp_distance: i64,
    /// `--allele-informative-reads-overlap-margin` (2)
    pub informative_read_overlap_margin: i64,
    /// `--standard-min-confidence-threshold-for-calling` (30)
    pub standard_confidence_for_calling: f64,
    /// `--max-alternate-alleles` (6)
    pub max_alternate_alleles: usize,
    /// `--max-genotype-count` (1024)
    pub max_genotype_count: usize,
    /// Do not draw the QD jitter here: write [`QD_JITTER_PLACEHOLDER`] as the QD value instead, and let
    /// the caller replace the placeholders, in output order, with `format_fixed(30 + 3 * nextGaussian(), 2)`
    /// drawn from its own run-wide generator (see `resolve_qd_jitter`). Default false.
    pub defer_qd_jitter: bool,
    /// `-contamination` (0: off). Genotyping then draws from [`GenotypingEngine::contam_rng`].
    pub contamination: f64,
    /// `-G AS_StandardAnnotation` (GVCF mode only: GATK refuses allele-specific annotations in VCF mode)
    pub as_annotations: bool,
    /// number of samples (VCF mode only when > 1); reads carry their sample index, samples are in sorted-name
    /// order (GATK's `SampleList`)
    pub n_samples: usize,
}

/// QD value written for a jittered QD when `EngineConfig::defer_qd_jitter` is set (a byte that never occurs
/// in VCF text).
pub const QD_JITTER_PLACEHOLDER: &str = "\u{1}";

/// Replaces every [`QD_JITTER_PLACEHOLDER`] in `text` (in order) with the jittered QD that
/// `QualByDepth` draws from `rng` (`30 + nextGaussian() * 3`, formatted `%.2f`). Returns the number of draws.
pub fn resolve_qd_jitter(text: &[u8], rng: &mut JavaRandom, out: &mut Vec<u8>) -> u64 {
    let mut n = 0;
    let mut last = 0;
    for (i, &b) in text.iter().enumerate() {
        if b == 1 {
            out.extend_from_slice(&text[last..i]);
            let qd = 30.0 + rng.next_gaussian() * 3.0;
            out.extend_from_slice(format_fixed(qd, 2).as_bytes());
            last = i + 1;
            n += 1;
        }
    }
    out.extend_from_slice(&text[last..]);
    n
}

impl EngineConfig {
    pub fn new(contig: &str, contig_length: i64) -> Self {
        EngineConfig {
            contig: contig.to_string(),
            contig_length,
            max_mnp_distance: 0,
            informative_read_overlap_margin: 2,
            standard_confidence_for_calling: 30.0,
            max_alternate_alleles: 6,
            max_genotype_count: 1024,
            defer_qd_jitter: false,
            contamination: 0.0,
            as_annotations: false,
            n_samples: 1,
        }
    }

    /// Settings `HaplotypeCallerEngine` forces in reference-confidence mode (`-ERC GVCF`): the
    /// calling confidence becomes -0.0 (genotyping must then go through `genotype_region_erc`).
    pub fn new_erc(contig: &str, contig_length: i64) -> Self {
        EngineConfig { standard_confidence_for_calling: -0.0, ..Self::new(contig, contig_length) }
    }
}

/// Inputs of one region, as `callRegion` has them right after `computeReadLikelihoods`.
pub struct RegionInput<'a> {
    /// `readLikelihoods.alleles()` (contains the reference haplotype)
    pub haplotypes: &'a [Haplotype],
    /// `readLikelihoods.sampleEvidence(0)`: reads in evidence order, before realignment
    pub reads: &'a [Read],
    /// `[haplotype][read]` log10 likelihoods
    pub matrix: &'a [Vec<f64>],
    /// reads removed by `filterNonPassingReads` (perSampleFilteredReadList): added to the annotation
    /// likelihoods with likelihood 0 for every allele when they overlap the variant (+/- 2 bp)
    pub filtered_reads: &'a [Read],
    /// `readLikelihoods.filteredSampleEvidence(0)` (evidence dropped by the likelihood engine). GATK
    /// carries it through `marginalize`/`retainEvidence`, but no default annotation reads it, so it
    /// has no effect on the output; accepted for completeness.
    pub filtered_evidence: &'a [Read],
    /// `assemblyResult.getFullReferenceWithPadding()` and the start of `getPaddedReferenceLoc()`
    pub ref_padded: &'a [u8],
    pub ref_padded_start: i64,
    /// `regionForGenotyping.getSpan()` (1-based inclusive)
    pub span_start: i64,
    pub span_end: i64,
}

/// Cross-region state of GATK's genotyping engine plus a reusable aligner. Regions must be
/// passed in genomic order: `upstreamDeletionsLoc` (spanning-deletion bookkeeping) and the QD
/// jitter random generator (`Utils.getRandomGenerator()`) carry over from region to region.
pub struct GenotypingEngine {
    pub cfg: EngineConfig,
    pub(crate) af: AfCalculator,
    pub(crate) upstream_deletions: Vec<(i64, i64)>,
    pub(crate) rng: JavaRandom,
    pub(crate) sw: SwAligner,
    /// number of QD values that were replaced by the random jitter
    pub n_qd_jitter: u64,
    /// GATK's `Utils.getRandomDataGenerator()` while contamination downsampling is on (the driver hands it from
    /// contig to contig in GATK's order)
    pub contam_rng: Option<crate::well::Well19937c>,
}

pub(crate) struct Call {
    /// indices into the merged allele list: ref first, then emitted alts
    pub(crate) alleles: Vec<usize>,
    pub(crate) log10_p_error: f64,
    /// the first sample's genotype
    pub(crate) genotype: genotyper::CalledGenotype,
    /// the other samples' genotypes, in sample order
    pub(crate) others: Vec<genotyper::CalledGenotype>,
    pub(crate) mleac: Vec<i32>,
}

impl Call {
    /// every sample's genotype, in sample order
    pub(crate) fn genotypes(&self) -> impl Iterator<Item = &genotyper::CalledGenotype> {
        std::iter::once(&self.genotype).chain(self.others.iter())
    }
}

impl GenotypingEngine {
    pub fn new(cfg: EngineConfig) -> Self {
        GenotypingEngine {
            cfg,
            af: AfCalculator::default(),
            upstream_deletions: Vec::new(),
            rng: JavaRandom::new(GATK_RANDOM_SEED),
            sw: SwAligner::new(),
            n_qd_jitter: 0,
            contam_rng: None,
        }
    }

    /// State of the QD-jitter generator (GATK's global `Utils.getRandomGenerator()`). With read
    /// downsampling on, the driver keeps ONE generator for the downsampler and the QD jitter: it
    /// calls `set_rng_state` before and `rng_state` after each `call_region` / `genotype_region`.
    pub fn rng_state(&self) -> (u64, f64, bool) {
        self.rng.state()
    }
    pub fn set_rng_state(&mut self, s: (u64, f64, bool)) {
        self.rng.set_state(s)
    }

    /// Realigns the reads (step 1) then genotypes the region (step 2). Returns the VCF records.
    pub fn call_region(&mut self, input: &RegionInput) -> Vec<VcfRecord> {
        let mut reads = input.reads.to_vec();
        realign::realign_reads_to_best_haplotype(
            &mut self.sw,
            &mut reads,
            input.haplotypes,
            input.matrix,
            input.ref_padded_start,
        );
        self.genotype_region(input, &reads)
    }

    /// Step 2 alone: `realigned` are the evidence reads after `realignReadsToTheirBestHaplotype`.
    pub fn genotype_region(&mut self, input: &RegionInput, realigned: &[Read]) -> Vec<VcfRecord> {
        let haps = input.haplotypes;
        let maps: Vec<EventMap> = haps
            .iter()
            .map(|h| EventMap::from_haplotype(h, input.ref_padded, input.ref_padded_start, self.cfg.max_mnp_distance))
            .collect();
        let starts = event_start_positions(&maps);
        let mut out = Vec::new();
        let n_reads = realigned.len();
        for loc in starts {
            if loc < input.span_start || loc > input.span_end {
                continue;
            }
            let mut merged = match merged_vc_at(loc, &maps, input.ref_padded, input.ref_padded_start) {
                Some(m) => m,
                None => continue,
            };
            let n_alleles_before_trim = merged.alleles.len();
            stats::inc(&stats::LOCI_EVALUATED);
            if merged.alleles.len() > 2 {
                stats::inc(&stats::MULTIALLELIC_LOCI);
            }
            if merged.alleles.iter().any(|a| a.is_span_del()) {
                stats::inc(&stats::SPAN_DEL_LOCI);
            }
            let (mapper_alleles, mut mapper_haps) = create_allele_mapper(&merged, loc, &maps);
            assert!(mapper_alleles == merged.alleles, "allele mapper differs from merged VC alleles at {}", loc);
            // removeAltAllelesIfTooManyGenotypes
            let practical = max_acceptable_allele_count(self.cfg.max_genotype_count);
            if mapper_alleles.len() > practical {
                stats::inc(&stats::TOO_MANY_GENOTYPES_REDUCTION);
                let keep = alleles_to_keep_by_hap_scores(&mapper_alleles, &mapper_haps, haps, practical);
                merged.alleles = keep.iter().map(|&a| mapper_alleles[a].clone()).collect();
                mapper_haps = keep.iter().map(|&a| mapper_haps[a].clone()).collect();
            }
            let alleles = &merged.alleles;
            let n = alleles.len();
            let vc_end = merged.end();
            // readAlleleLikelihoods.retainEvidence(read overlaps [start-margin, end+margin])
            let lo = (loc - self.cfg.informative_read_overlap_margin).max(1);
            let hi = (vc_end + self.cfg.informative_read_overlap_margin).min(self.cfg.contig_length);
            let retained: Vec<usize> = (0..n_reads).filter(|&r| realigned[r].start <= hi && lo <= realigned[r].end).collect();
            // marginalize: max over the haplotypes of each allele
            let mut lk: Vec<Vec<f64>> = Vec::with_capacity(n);
            for a in 0..n {
                let mut row = vec![f64::NEG_INFINITY; retained.len()];
                for (k, &r) in retained.iter().enumerate() {
                    let mut m = f64::NEG_INFINITY;
                    for &h in &mapper_haps[a] {
                        m = m.max(input.matrix[h][r]);
                    }
                    row[k] = m;
                }
                lk.push(row);
            }
            // calculateGLsForThisEvent: one genotype per sample (its evidence in evidence order; no evidence: PLs 0)
            let ns = self.cfg.n_samples.max(1);
            let mut pls: Vec<Vec<i32>> = Vec::with_capacity(ns);
            if ns == 1 {
                let gls = match self.contamination_removals(&lk, retained.len()) {
                    Some(rm) => genotyper::log10_genotype_likelihoods(&crate::contamination::keep_columns(&lk, &rm), retained.len() - rm.len()),
                    None => genotyper::log10_genotype_likelihoods(&lk, retained.len()),
                };
                pls.push(genotyper::gls_to_pls(&gls));
            } else {
                let mut cols: Vec<Vec<usize>> = vec![Vec::new(); ns];
                for (k, &r) in retained.iter().enumerate() {
                    cols[realigned[r].sample as usize].push(k);
                }
                for c in &cols {
                    let lk_s: Vec<Vec<f64>> = lk.iter().map(|row| c.iter().map(|&k| row[k]).collect()).collect();
                    let gls = match self.contamination_removals(&lk_s, c.len()) {
                        Some(rm) => genotyper::log10_genotype_likelihoods(&crate::contamination::keep_columns(&lk_s, &rm), c.len() - rm.len()),
                        None => genotyper::log10_genotype_likelihoods(&lk_s, c.len()),
                    };
                    pls.push(genotyper::gls_to_pls(&gls));
                }
            }
            let call = match self.calculate_genotypes(loc, alleles, &pls) {
                Some(c) => c,
                None => continue,
            };

            // annotation likelihoods: retained evidence + overlapping filtered reads (likelihood 0)
            let filt: Vec<&Read> = input.filtered_reads.iter().filter(|r| r.start <= hi && lo <= r.end).collect();
            let mut ev: Vec<&Read> = retained.iter().map(|&r| &realigned[r]).collect();
            ev.extend(filt.iter().copied());
            let n_ev = ev.len();
            let lk_col = |a: usize, e: usize| -> f64 { if e < retained.len() { lk[a][e] } else { 0.0 } };
            stats::inc(&stats::CALLS);
            if !filt.is_empty() {
                stats::inc(&stats::FILTERED_READS_IN_ANNOTATION);
            }
            if call.alleles.len() != n_alleles_before_trim {
                stats::inc(&stats::ALLELES_DROPPED);
            }
            if call.alleles.iter().any(|&a| alleles[a].is_span_del()) {
                stats::inc(&stats::SPAN_DEL_EMITTED);
            }
            if call.genotype.alleles.is_none() {
                stats::inc(&stats::NO_CALL_GENOTYPE);
            }
            let rec = self.annotate_call(loc, vc_end, alleles, &call, &ev, n_ev, &lk_col, n_alleles_before_trim);
            out.push(rec);
        }
        out
    }

    /// `readAlleleLikelihoods.contaminationDownsampling(...)` when `-contamination` is on: the evidence removed for
    /// genotyping (sorted), or None when nothing is removed. `lk` rows: the event's alleles (no `<NON_REF>`).
    pub(crate) fn contamination_removals(&mut self, lk: &[Vec<f64>], n_ev: usize) -> Option<Vec<usize>> {
        if !(self.cfg.contamination > 0.0) {
            return None;
        }
        let rng = self.contam_rng.as_mut().expect("contamination downsampling without its generator");
        let rm = crate::contamination::removals(lk, n_ev, self.cfg.contamination, rng);
        (!rm.is_empty()).then_some(rm)
    }

    pub(crate) fn is_covered_by_deletion(&self, start: i64) -> bool {
        self.upstream_deletions.iter().any(|&(s, e)| s < start && start <= e)
    }

    /// `GenotypingEngine.calculateGenotypes` (USE_PLS_TO_ASSIGN, no forced alleles); `pls[s]`: sample s's PLs.
    fn calculate_genotypes(&mut self, loc: i64, alleles: &[Allele], pls: &[Vec<i32>]) -> Option<Call> {
        let n = alleles.len();
        if n > 50 {
            return None; // GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED
        }
        // reduce to the most likely max-alternate-alleles
        let reduced: Vec<usize> = if n - 1 > self.cfg.max_alternate_alleles {
            stats::inc(&stats::MAX_ALT_ALLELES_REDUCTION);
            most_likely_alleles(pls, n, self.cfg.max_alternate_alleles)
        } else {
            (0..n).collect()
        };
        let reduced_pls: Vec<Vec<i32>> =
            if reduced.len() == n { pls.to_vec() } else { pls.iter().map(|p| genotyper::subset_and_call(p, &reduced).pl).collect() };
        let lens: Vec<usize> = reduced.iter().map(|&a| alleles[a].len()).collect();
        let sd = reduced.iter().position(|&a| alleles[a].is_span_del());
        let refs: Vec<&[i32]> = reduced_pls.iter().map(|p| &p[..]).collect();
        let af = self.af.calculate_multi(&refs, &lens, sd);

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
            if plausible && !spurious {
                out_alts.push(k);
            }
            monomorphic &= !(plausible && !spurious);
        }
        if monomorphic {
            return None; // passesEmitThreshold requires a non-monomorphic site
        }
        let log10_conf = af.log10_p_no_variant + 0.0;
        let phred = (-10.0 * log10_conf) + 0.0;
        if !(phred >= self.cfg.standard_confidence_for_calling) {
            return None;
        }
        if out_alts.len() == 1 && alleles[reduced[out_alts[0]]].is_span_del() {
            stats::inc(&stats::ONLY_SPAN_DEL_DROPPED);
            return None;
        }
        let mut out_alleles = vec![0usize];
        out_alleles.extend(out_alts.iter().map(|&k| reduced[k]));
        // recordDeletions
        self.upstream_deletions.retain(|&(_, e)| e >= loc);
        let ref_len = alleles[0].len() as i64;
        for &a in &out_alleles {
            let del = ref_len - alleles[a].len() as i64;
            if del > 0 {
                self.upstream_deletions.push((loc, loc + del));
            }
        }
        let mut gts = pls.iter().map(|p| genotyper::subset_and_call(p, &out_alleles));
        let genotype = gts.next().unwrap();
        let others = gts.collect();
        let mleac = out_alts.iter().map(|&k| af.mle_alt_counts[k - 1]).collect();
        Some(Call { alleles: out_alleles, log10_p_error: log10_conf, genotype, others, mleac })
    }

    #[allow(clippy::too_many_arguments)]
    fn annotate_call(
        &mut self,
        loc: i64,
        vc_end: i64,
        merged: &[Allele],
        call: &Call,
        ev: &[&Read],
        n_ev: usize,
        lk: &dyn Fn(usize, usize) -> f64,
        n_alleles_before_trim: usize,
    ) -> VcfRecord {
        let n = merged.len();
        let ns = self.cfg.n_samples.max(1);
        let call_alleles = &call.alleles;
        let in_call = |a: usize| call_alleles.contains(&a);
        let gts: Vec<&genotyper::CalledGenotype> = call.genotypes().collect();
        // the samples' evidence (annotation likelihoods: retained evidence, then the overlapping filtered reads)
        let sample_of = |e: usize| if ns == 1 { 0 } else { ev[e].sample as usize };
        let mut ev_count = vec![0usize; ns];
        for e in 0..n_ev {
            ev_count[sample_of(e)] += 1;
        }

        // genotype annotations of the called genotypes: AD and DP over the call's alleles
        let mut ad: Vec<Option<Vec<i32>>> = vec![None; ns];
        let mut dp: Vec<Option<i32>> = vec![None; ns];
        if gts.iter().any(|g| g.alleles.is_some()) {
            let pr = ref_priorities(call_alleles.len(), Some(0));
            let mut counts = vec![vec![0i32; call_alleles.len()]; ns];
            let mut informative = vec![0i32; ns];
            let mut col = vec![0f64; call_alleles.len()];
            for e in 0..n_ev {
                let s = sample_of(e);
                if gts[s].alleles.is_none() {
                    continue;
                }
                for (k, &a) in call_alleles.iter().enumerate() {
                    col[k] = lk(a, e);
                }
                let b = best_allele(&col, Some(&pr));
                if b.is_informative() {
                    counts[s][b.allele.unwrap()] += 1;
                    informative[s] += 1;
                }
            }
            for s in 0..ns {
                if gts[s].alleles.is_some() {
                    ad[s] = Some(std::mem::take(&mut counts[s]));
                    dp[s] = Some(if ev_count[s] == 0 { 0 } else { informative[s] });
                }
            }
        }

        // best allele of every evidence over all merged alleles (ref priority)
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

        let mut info: Vec<(&'static str, String)> = Vec::new();
        // ChromosomeCounts (no called chromosome: none)
        let an: i32 = gts.iter().filter(|g| g.alleles.is_some()).count() as i32 * 2;
        if an > 0 {
            let mut ac = Vec::new();
            let mut afs = Vec::new();
            for k in 1..call_alleles.len() {
                let c: i32 = gts.iter().filter_map(|g| g.alleles).map(|(i, j)| (i == k) as i32 + (j == k) as i32).sum();
                ac.push(c.to_string());
                afs.push(crate::javamath::format_vcf_double(c as f64 / an as f64));
            }
            info.push(("AC", ac.join(",")));
            info.push(("AF", afs.join(",")));
            info.push(("AN", an.to_string()));
        }
        // rank sum tests
        let (mut bq_ref, mut bq_alt, mut mq_ref, mut mq_alt, mut rp_ref, mut rp_alt) =
            (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for (e, b) in best.iter().enumerate() {
            let a = match b {
                Some(a) => *a,
                None => continue,
            };
            let r = ev[e];
            let usable = r.mapq != 0 && r.mapq != 255;
            if !usable {
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
        // Coverage
        if n_ev > 0 {
            info.push(("DP", n_ev.to_string()));
        }
        // ExcessHet (samples with a called genotype)
        let het_in: Vec<(bool, &[i32])> = gts.iter().map(|g| (g.alleles.is_some(), &g.pl[..])).collect();
        if let Some(v) = annotate::excess_het_samples(&het_in, call_alleles.len()) {
            info.push(("ExcessHet", v));
        }
        // FisherStrand / StrandOddsRatio: per-sample tables, each kept if its total exceeds the minimum count
        let is_ref_a = |a: usize| a == 0;
        let is_alt_a = |a: usize| a != 0 && call_alleles.contains(&a);
        let sample_tables = |min_count: i32| -> [[i32; 2]; 2] {
            if ns == 1 {
                return annotate::strand_table(&best, ev, &is_ref_a, &is_alt_a, min_count);
            }
            let mut t = [[0i32; 2]; 2];
            for s in 0..ns {
                let (b_s, ev_s): (Vec<Option<usize>>, Vec<&Read>) = (0..n_ev).filter(|&e| sample_of(e) == s).map(|e| (best[e], ev[e])).unzip();
                let x = annotate::strand_table(&b_s, &ev_s, &is_ref_a, &is_alt_a, min_count);
                for i in 0..2 {
                    for j in 0..2 {
                        t[i][j] += x[i][j];
                    }
                }
            }
            t
        };
        let fs_table = sample_tables(2);
        info.push(("FS", annotate::fisher_strand(fs_table)));
        // InbreedingCoeff (10 samples or more, unrounded genotype counts)
        if ns >= 10 {
            let [refs, hets, homs] = annotate::diploid_genotype_counts(&het_in, call_alleles.len(), false);
            let p = (2.0 * refs + hets) / (2.0 * (refs + hets + homs));
            let q = 1.0 - p;
            let ic = 1.0 - (hets / (2.0 * p * q * ns as f64));
            info.push(("InbreedingCoeff", format_fixed(ic, 4)));
        }
        // MLEAC / MLEAF
        let an_called = an as f64;
        info.push(("MLEAC", call.mleac.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(",")));
        info.push((
            "MLEAF",
            call.mleac
                .iter()
                .map(|&c| crate::javamath::format_vcf_double((1.0f64).min(c as f64 / an_called)))
                .collect::<Vec<_>>()
                .join(","),
        ));
        // RMSMappingQuality (raw, then finalizeRawMQ)
        let (mut ss, mut nmq) = (0i64, 0i64);
        for r in ev {
            let mq = r.mapq as i64;
            if mq != 255 {
                ss += mq * mq;
                nmq += 1;
            }
        }
        if nmq > 0 {
            info.push(("MQ", format_fixed((ss as f64 / nmq as f64).sqrt(), 2)));
        }
        if let Some(s) = annotate::rank_sum(&mq_alt, &mq_ref) {
            info.push(("MQRankSum", s));
        }
        // QualByDepth: depth of the het / hom-var genotypes
        let qual = -10.0 * call.log10_p_error;
        let mut depth = 0i32;
        let mut ad_restricted = 0i32;
        let mut any_variant = false;
        for (s, g) in gts.iter().enumerate() {
            let Some((i, j)) = g.alleles else { continue };
            if i == 0 && j == 0 {
                continue;
            }
            any_variant = true;
            let ad_v = ad[s].as_ref().unwrap();
            let total: i32 = ad_v.iter().sum();
            if total != 0 {
                if total - ad_v[0] > 1 {
                    ad_restricted += total;
                }
                depth += total;
            } else {
                depth += ev_count[s] as i32;
            }
        }
        if ad_restricted > 0 {
            depth = ad_restricted;
        }
        if any_variant && depth != 0 {
            let mut qd = qual / depth as f64;
            if !(qd < 35.0) {
                stats::inc(&stats::QD_JITTER);
                self.n_qd_jitter += 1;
                if self.cfg.defer_qd_jitter {
                    info.push(("QD", QD_JITTER_PLACEHOLDER.to_string()));
                } else {
                    qd = 30.0 + self.rng.next_gaussian() * 3.0;
                    info.push(("QD", format_fixed(qd, 2)));
                }
            } else {
                info.push(("QD", format_fixed(qd, 2)));
            }
        }
        if let Some(s) = annotate::rank_sum(&rp_alt, &rp_ref) {
            info.push(("ReadPosRankSum", s));
        }
        let sor_table = sample_tables(0);
        info.push(("SOR", annotate::strand_odds_ratio(sor_table)));
        info.sort_by(|a, b| a.0.cmp(b.0));

        // reverse trimming when alleles were dropped
        let mut out_alleles: Vec<Vec<u8>> = call_alleles.iter().map(|&a| merged[a].bases.clone()).collect();
        let mut pos = loc;
        if call_alleles.len() != n_alleles_before_trim {
            let (p, trimmed) = reverse_trim_alleles(loc, &out_alleles);
            if trimmed != out_alleles {
                stats::inc(&stats::REVERSE_TRIMMED);
            }
            pos = p;
            out_alleles = trimmed;
        }
        let g = &call.genotype;
        let mut ad_it = ad.into_iter();
        let mut dp_it = dp.into_iter();
        let (ad0, dp0) = (ad_it.next().unwrap(), dp_it.next().unwrap());
        let more = call
            .others
            .iter()
            .zip(ad_it.zip(dp_it))
            .map(|(g, (ad, dp))| crate::vcf::SampleGenotype { gt: g.alleles, ad, dp, gq: g.gq, pl: g.pl.clone() })
            .collect();
        VcfRecord {
            contig: self.cfg.contig.clone(),
            pos,
            id: None,
            alleles: out_alleles,
            qual,
            info,
            gt: g.alleles,
            ad: ad0,
            dp: dp0,
            gq: g.gq,
            pl: g.pl.clone(),
            sb: None,
            phase: None,
            more,
        }
    }
}

/// `GenotypeIndexCalculator.computeMaxAcceptableAlleleCount(2, maxGenotypeCount)`
pub(crate) fn max_acceptable_allele_count(max_genotype_count: usize) -> usize {
    let mut n = 1;
    while genotyper::genotype_count(n + 1) <= max_genotype_count {
        n += 1;
    }
    n
}

/// `HaplotypeCallerGenotypingEngine.whichAllelesToKeepBasedonHapScores`: keeps the `desired`
/// alleles with the best (then second best) supporting haplotype scores; returns kept indices in
/// their original order. Uses the same `java.util.PriorityQueue` sift operations as GATK.
pub(crate) fn alleles_to_keep_by_hap_scores(alleles: &[Allele], haps_of: &[Vec<usize>], haps: &[Haplotype], desired: usize) -> Vec<usize> {
    if alleles.len() <= desired {
        return (0..alleles.len()).collect();
    }
    // Double.compareTo ordering (NaN largest, -0.0 < 0.0)
    let dcmp = |a: f64, b: f64| -> std::cmp::Ordering {
        if a < b {
            std::cmp::Ordering::Less
        } else if a > b {
            std::cmp::Ordering::Greater
        } else {
            let (x, y) = (if a.is_nan() { 0x7ff8000000000000u64 as i64 } else { a.to_bits() as i64 }, if b.is_nan() { 0x7ff8000000000000u64 as i64 } else { b.to_bits() as i64 });
            x.cmp(&y)
        }
    };
    struct Scored {
        allele: usize,
        best: f64,
        second: f64,
    }
    let items: Vec<Scored> = (0..alleles.len())
        .map(|a| {
            let mut s: Vec<f64> = haps_of[a].iter().map(|&h| haps[h].score).collect();
            s.sort_by(|x, y| dcmp(*x, *y));
            let best = if !s.is_empty() { s[s.len() - 1] } else { f64::NEG_INFINITY };
            let second = if s.len() > 1 { s[s.len() - 2] } else { f64::NEG_INFINITY };
            Scored { allele: a, best, second }
        })
        .collect();
    let bits_eq = |a: f64, b: f64| (a.is_nan() && b.is_nan()) || a.to_bits() == b.to_bits();
    let cmp = |x: &Scored, y: &Scored| -> i32 {
        let (ax, ay) = (&alleles[x.allele], &alleles[y.allele]);
        if ax.is_ref && !ay.is_ref {
            -1
        } else if !ax.is_ref && ay.is_ref {
            1
        } else if x.best > y.best {
            -1
        } else if x.best < y.best {
            1
        } else if !bits_eq(x.second, y.second) {
            if x.second > y.second { -1 } else { 1 }
        } else {
            // Allele.compareTo: bases lexicographically (signed bytes), then length
            let (a, b) = (&ax.bases, &ay.bases);
            for i in 0..a.len().min(b.len()) {
                let d = a[i] as i8 as i32 - b[i] as i8 as i32;
                if d != 0 {
                    return if d < 0 { -1 } else { 1 };
                }
            }
            (a.len() as i64 - b.len() as i64).signum() as i32
        }
    };
    // java.util.PriorityQueue: offer (siftUp) then poll (siftDown)
    let mut q: Vec<usize> = Vec::new();
    for i in 0..items.len() {
        let mut k = q.len();
        q.push(i);
        while k > 0 {
            let parent = (k - 1) >> 1;
            if cmp(&items[q[k]], &items[q[parent]]) >= 0 {
                break;
            }
            q.swap(k, parent);
            k = parent;
        }
    }
    let mut retain: Vec<usize> = Vec::new();
    while retain.len() < desired {
        let n = q.len() - 1;
        q.swap(0, n);
        let top = q.pop().unwrap();
        if n > 0 {
            let mut k = 0usize;
            let half = n >> 1;
            while k < half {
                let mut child = 2 * k + 1;
                let right = child + 1;
                if right < n && cmp(&items[q[child]], &items[q[right]]) > 0 {
                    child = right;
                }
                if cmp(&items[q[k]], &items[q[child]]) <= 0 {
                    break;
                }
                q.swap(k, child);
                k = child;
            }
        }
        if !retain.contains(&items[top].allele) {
            retain.push(items[top].allele);
        }
    }
    (0..alleles.len()).filter(|a| retain.contains(a)).collect()
}

/// `AlleleSubsettingUtils.calculateMostLikelyAlleles(vc, 2, numAltAllelesToKeep, false)`: likelihood sums over
/// the samples (in sample-name order).
fn most_likely_alleles(pls: &[Vec<i32>], n: usize, keep: usize) -> Vec<usize> {
    let mut sums = vec![0f64; n];
    for p in pls {
        let gls: Vec<f64> = p.iter().map(|&p| p as f64 / -10.0).collect();
        let best = crate::mathutils::max_element_index(&gls);
        let diff = (gls[best] - gls[0]).abs();
        let (i, j) = genotyper::genotype_alleles(n)[best];
        for a in 1..n {
            if a == i || a == j {
                sums[a] += diff;
            }
        }
    }
    let mut idx: Vec<usize> = (1..n).collect();
    // stable sort by decreasing likelihood sum
    idx.sort_by(|&a, &b| sums[b].partial_cmp(&sums[a]).unwrap());
    let kept: Vec<usize> = idx.into_iter().take(keep).collect();
    (0..n).filter(|&a| a == 0 || kept.contains(&a)).collect()
}

/// `GATKVariantContextUtils.reverseTrimAlleles` (trimAlleles(vc, false, true)) on allele bases (ref first).
pub fn reverse_trim_alleles(start: i64, alleles: &[Vec<u8>]) -> (i64, Vec<Vec<u8>>) {
    // symbolic alleles (<NON_REF>) have length 0 and, like `*`, are never trimmed
    let exempt = |a: &[u8]| a == b"*" || a.first() == Some(&b'<');
    if alleles.len() <= 1 || alleles.iter().any(|a| a.len() == 1 && a.as_slice() != b"*") {
        return (start, alleles.to_vec());
    }
    let seqs: Vec<&Vec<u8>> = alleles.iter().filter(|a| !exempt(a)).collect();
    let mut from: Vec<i64> = vec![0; seqs.len()];
    let mut to: Vec<i64> = seqs.iter().map(|s| s.len() as i64).collect();
    // AlignmentUtils.normalizeAlleles(seqs, ranges, maxShift=0, trim=true)
    let mut start_shift = 0i64;
    let mut end_shift = 0i64;
    let mut min_size = (0..seqs.len()).map(|i| to[i] - from[i]).min().unwrap();
    let last_same = |from: &[i64], to: &[i64]| {
        let _ = from;
        let b = seqs[0][(to[0] - 1) as usize];
        (0..seqs.len()).all(|i| seqs[i][(to[i] - 1) as usize] == b)
    };
    let first_same = |from: &[i64]| {
        let b = seqs[0][from[0] as usize];
        (0..seqs.len()).all(|i| seqs[i][from[i] as usize] == b)
    };
    let next_left_same = |from: &[i64]| {
        let b = seqs[0][(from[0] - 1) as usize];
        (0..seqs.len()).all(|i| seqs[i][(from[i] - 1) as usize] == b)
    };
    while min_size > 0 && last_same(&from, &to) {
        for t in to.iter_mut() {
            *t -= 1;
        }
        min_size -= 1;
        end_shift += 1;
    }
    while min_size > 0 && first_same(&from) {
        for f in from.iter_mut() {
            *f += 1;
        }
        min_size -= 1;
        start_shift -= 1;
    }
    while start_shift < 0 && next_left_same(&from) && last_same(&from, &to) {
        for i in 0..seqs.len() {
            from[i] -= 1;
            to[i] -= 1;
        }
        start_shift += 1;
        end_shift += 1;
    }
    let end_trim = end_shift;
    let start_trim = -start_shift;
    let empty = (0..seqs.len()).any(|i| to[i] - from[i] == 0);
    let restore_end = empty && start_trim == 0;
    let end_clip = if restore_end { end_trim - 1 } else { end_trim };
    // trimAlleles(vc, fwdTrimEnd = -1, revTrim = end_clip)
    if end_clip == 0 {
        return (start, alleles.to_vec());
    }
    let trimmed = alleles
        .iter()
        .map(|a| if exempt(a) { a.clone() } else { a[..a.len() - end_clip as usize].to_vec() })
        .collect();
    (start, trimmed)
}
