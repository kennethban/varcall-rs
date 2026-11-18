//! Rust port of the part of GATK 4.7 `HaplotypeCallerEngine.callRegion` that runs between local
//! assembly and genotyping, with default HaplotypeCaller arguments:
//!
//! 1. variation events of the assembled haplotypes (`AssemblyResultSet.getVariationEvents(0)`);
//! 2. `AssemblyRegionTrimmer.trim` (stop with "no variation" when no event overlaps the active span);
//! 3. `AssemblyResultSet.trimTo(variantRegion)` (trimmed, merged and re-sorted haplotypes);
//! 4. `getRegionForGenotyping()` reads (hard-clipped to the trimmed padded span), stub-read removal,
//!    `filterNonPassingReads`, and the `isVariationPresent` / empty-region early exits;
//! 5. `PairHMMLikelihoodCalculationEngine.computeReadLikelihoods`: read quality modification,
//!    the PairHMM (CUDA, through `hc_cudahmm`), likelihood normalisation and
//!    `filterPoorlyModeledEvidence`.
//!
//! The API is two-phase so that the PairHMM of many regions can run in one GPU launch:
//! [`prepare`] -> ([`Prepared`], PairHMM inputs via [`Prepared::hmm_group`]) -> caller batches
//! groups into `hc_cudahmm::compute` (or [`compute_likelihoods_gpu`]) -> [`finish`].
//!
//! See README.md for the replicated GATK behaviours.

pub mod clip;
pub mod cpu;
pub mod events;
pub mod hmm;
pub mod trim;

pub use hc_assembly::{Haplotype, Read};
pub use hc_cudahmm::{HmmGroup, HmmRead};
pub use hmm::PcrErrorModel;

use std::time::{Duration, Instant};

/// 1-based inclusive interval on the region's contig.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct Span {
    pub start: i64,
    pub end: i64,
}

/// Arguments that affect this stage.  `Default` = GATK 4.7 HaplotypeCaller defaults.
#[derive(Clone, Debug)]
pub struct Params {
    /// `--max-mnp-distance` (HaplotypeCaller: 0)
    pub max_mnp_distance: usize,
    /// `--padding-around-snps` (20)
    pub snp_padding: i64,
    /// `--padding-around-indels` (75)
    pub indel_padding: i64,
    /// `--padding-around-strs` (75)
    pub str_padding: i64,
    /// `AssemblyBasedCallerUtils.MINIMUM_READ_LENGTH_AFTER_TRIMMING` (10)
    pub min_read_length_after_trimming: usize,
    /// `HaplotypeCallerEngine.READ_LENGTH_FILTER_THRESHOLD` (10)
    pub read_length_filter_threshold: usize,
    /// `--mapping-quality-threshold-for-genotyping` (20)
    pub mapping_quality_threshold: u8,
    /// `--disable-optimizations` (false)
    pub disable_optimizations: bool,
    /// `--pair-hmm-gap-continuation-penalty` (10)
    pub gcp: u8,
    /// `--pcr-indel-model` (CONSERVATIVE)
    pub pcr_error_model: PcrErrorModel,
    /// `--base-quality-score-threshold` (`PairHMM.BASE_QUALITY_SCORE_THRESHOLD` = 18)
    pub base_quality_score_threshold: u8,
    /// `QualityUtils.MIN_USABLE_Q_SCORE` (6)
    pub min_usable_q_score: u8,
    /// `ReadUtils.DEFAULT_INSERTION_DELETION_QUAL` (45; reads carry no BI/BD tags)
    pub default_indel_qual: u8,
    /// `--disable-cap-base-qualities-to-map-quality` (false)
    pub disable_cap_read_qualities_to_mapq: bool,
    /// `--phred-scaled-global-read-mismapping-rate` (45; negative disables the cap)
    pub phred_scaled_global_read_mismapping_rate: i32,
    /// `!--disable-symmetric-hmm-normalizing` (true)
    pub symmetrically_normalize_alleles_to_reference: bool,
    /// `--expected-mismatch-rate-for-read-disqualification` (0.02)
    pub expected_error_rate_per_base: f64,
    /// `--enable-dynamic-read-disqualification-for-genotyping` (false)
    pub dynamic_read_disqualification: bool,
    /// `--dynamic-read-disqualification-threshold` (1.0)
    pub read_disqualification_scale: f64,
    /// Mutect2's genotyping reads: `removeReadStubs` (fewer than 10 bases) only, no `filterNonPassingReads`
    pub mutect: bool,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            max_mnp_distance: 0,
            snp_padding: 20,
            indel_padding: 75,
            str_padding: 75,
            min_read_length_after_trimming: 10,
            read_length_filter_threshold: 10,
            mapping_quality_threshold: 20,
            disable_optimizations: false,
            gcp: 10,
            pcr_error_model: PcrErrorModel::Conservative,
            base_quality_score_threshold: 18,
            min_usable_q_score: 6,
            default_indel_qual: 45,
            disable_cap_read_qualities_to_mapq: false,
            phred_scaled_global_read_mismapping_rate: 45,
            symmetrically_normalize_alleles_to_reference: true,
            expected_error_rate_per_base: 0.02,
            dynamic_read_disqualification: false,
            read_disqualification_scale: 1.0,
            mutect: false,
        }
    }
}

impl Params {
    /// `AssemblyBasedCallerUtils.getGlobalMismatchingRateFromArgs`: `-inf` disables the cap.
    pub fn log10_global_read_mismapping_rate(&self) -> f64 {
        if self.phred_scaled_global_read_mismapping_rate < 0 {
            f64::NEG_INFINITY
        } else {
            self.phred_scaled_global_read_mismapping_rate as f64 * -0.1
        }
    }
}

/// One assembled active region (everything `callRegion` has right after `assembleReads`).
pub struct RegionInput<'a> {
    /// active span and padded span of the assembly region
    pub span: Span,
    pub padded_span: Span,
    /// `region.getReads()` after `finalizeRegion` (the reads given to the assembler), in order
    pub reads: &'a [Read],
    /// per read: paired, mate mapped and mate on another contig (`MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE`
    /// fails).  Empty slice = take it from the reads (`Read::mate_on_other_contig`, i.e. `ref_id` /
    /// `mate_ref_id`).
    pub mate_on_other_contig: &'a [bool],
    /// untrimmed `AssemblyResultSet.getHaplotypeList()` (includes the reference haplotype)
    pub haplotypes: &'a [Haplotype],
    /// `getFullReferenceWithPadding()` and the start of `getPaddedReferenceLoc()`
    pub ref_padded: &'a [u8],
    pub ref_padded_start: i64,
}

/// Which way `callRegion` leaves (or continues) after this stage.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Exit {
    /// trimmer found no event overlapping the active span ("no_variation_after_trim")
    NoVariationAfterTrim,
    /// `!assemblyResult.isVariationPresent()` after trimming ("no_assembly_variation")
    NoAssemblyVariation,
    /// no reads left for genotyping ("no_reads_for_genotyping")
    NoReadsForGenotyping,
    /// read likelihoods are computed and genotyping follows
    Genotyped,
}

impl Exit {
    pub fn as_str(self) -> &'static str {
        match self {
            Exit::NoVariationAfterTrim => "no_variation_after_trim",
            Exit::NoAssemblyVariation => "no_assembly_variation",
            Exit::NoReadsForGenotyping => "no_reads_for_genotyping",
            Exit::Genotyped => "genotyped",
        }
    }
}

/// Phase-1 result: trimming and genotyping reads, plus the PairHMM inputs when `exit == Genotyped`.
#[derive(Clone, Debug)]
pub struct Prepared {
    pub exit: Exit,
    /// variant span / padded variant span = regionForGenotyping span / padded span
    /// (`None` iff `exit == NoVariationAfterTrim`)
    pub variant_span: Option<Span>,
    pub variant_padded_span: Option<Span>,
    /// `assemblyResult.isVariationPresent()` of the trimmed result set
    pub assembly_variation_present: bool,
    /// trimmed haplotypes in `AssemblyResultSet` order (= likelihood allele order)
    pub haplotypes: Vec<Haplotype>,
    /// regionForGenotyping reads after stub removal and `filterNonPassingReads` (= PairHMM batch order)
    pub reads: Vec<Read>,
    /// reads removed by `filterNonPassingReads`
    pub filtered_reads: Vec<Read>,
    /// PairHMM qualities per read (empty unless `exit == Genotyped`)
    pub hmm: Vec<hmm::HmmQuals>,
    gcp: Vec<u8>,
}

impl Prepared {
    /// Whether this region needs PairHMM likelihoods (`exit == Genotyped`).
    pub fn needs_hmm(&self) -> bool {
        self.exit == Exit::Genotyped
    }
    /// The PairHMM batch of this region (reads in evidence order, haplotypes in allele order).
    pub fn hmm_group(&self) -> HmmGroup<'_> {
        let reads = self
            .reads
            .iter()
            .zip(&self.hmm)
            .map(|(r, q)| HmmRead { bases: &r.seq, quals: &q.quals, ins: &q.ins, del: &q.del, gcp: &self.gcp[..r.seq.len()] })
            .collect();
        HmmGroup { reads, haps: self.haplotypes.iter().map(|h| &h.bases[..]).collect() }
    }
    /// The reads with `keep[i]` (with their PairHMM inputs), e.g. one sample's evidence.
    pub fn subset(&self, keep: &[bool]) -> Prepared {
        let pick = |v: &[Read]| v.iter().zip(keep).filter(|x| *x.1).map(|x| x.0.clone()).collect::<Vec<_>>();
        Prepared {
            exit: self.exit,
            variant_span: self.variant_span,
            variant_padded_span: self.variant_padded_span,
            assembly_variation_present: self.assembly_variation_present,
            haplotypes: self.haplotypes.clone(),
            reads: pick(&self.reads),
            filtered_reads: Vec::new(),
            hmm: self.hmm.iter().zip(keep).filter(|x| *x.1).map(|x| x.0.clone()).collect(),
            gcp: self.gcp.clone(),
        }
    }

    /// The gap continuation penalties of read `i` (constant GCP).
    pub fn gcp(&self, i: usize) -> &[u8] {
        &self.gcp[..self.reads[i].seq.len()]
    }
}

/// `AlleleLikelihoods<GATKRead, Haplotype>` for the single sample.
#[derive(Clone, Debug)]
pub struct ReadLikelihoods {
    /// allele order (= `Prepared::haplotypes`)
    pub haplotypes: Vec<Haplotype>,
    /// evidence that passed `filterPoorlyModeledEvidence`, in order
    pub reads: Vec<Read>,
    /// `filteredSampleEvidence(0)`: poorly modelled reads, in order
    pub filtered_reads: Vec<Read>,
    /// log10 likelihoods, row-major `[haplotype][read]` (`values[h * reads.len() + r]`)
    pub values: Vec<f64>,
}

impl ReadLikelihoods {
    #[inline]
    pub fn get(&self, hap: usize, read: usize) -> f64 {
        self.values[hap * self.reads.len() + read]
    }
}

/// Phase-2 result: everything of [`Prepared`] that downstream needs, plus the likelihoods.
#[derive(Clone, Debug)]
pub struct RegionResult {
    pub exit: Exit,
    pub variant_span: Option<Span>,
    pub variant_padded_span: Option<Span>,
    pub assembly_variation_present: bool,
    /// reads removed by `filterNonPassingReads`
    pub filtered_reads: Vec<Read>,
    /// `Some` iff `exit == Genotyped`
    pub likelihoods: Option<ReadLikelihoods>,
    /// trimmed haplotypes / genotyping reads when no likelihoods were computed
    pub haplotypes: Vec<Haplotype>,
    pub reads: Vec<Read>,
}

/// Accumulated wall time per stage (for profiling; pass `None` to skip timing).
#[derive(Clone, Copy, Debug, Default)]
pub struct StageTimes {
    pub events: Duration,
    pub trim: Duration,
    pub trim_haplotypes: Duration,
    pub genotyping_reads: Duration,
    pub hmm_inputs: Duration,
    pub finish: Duration,
}

impl std::ops::AddAssign for StageTimes {
    fn add_assign(&mut self, o: Self) {
        self.events += o.events;
        self.trim += o.trim;
        self.trim_haplotypes += o.trim_haplotypes;
        self.genotyping_reads += o.genotyping_reads;
        self.hmm_inputs += o.hmm_inputs;
        self.finish += o.finish;
    }
}

struct Timer<'a> {
    t: Option<&'a mut StageTimes>,
    last: Instant,
}
impl<'a> Timer<'a> {
    fn new(t: Option<&'a mut StageTimes>) -> Self {
        Timer { t, last: Instant::now() }
    }
    fn lap(&mut self, f: impl FnOnce(&mut StageTimes) -> &mut Duration) {
        if let Some(t) = self.t.as_deref_mut() {
            let now = Instant::now();
            *f(t) += now - self.last;
            self.last = now;
        }
    }
}

/// Phase 1 with default parameters.
pub fn prepare(input: &RegionInput) -> Prepared {
    prepare_with(input, &Params::default(), None)
}

/// Phase 1: events, trimming, trimmed haplotypes, genotyping reads and PairHMM inputs.
pub fn prepare_with(input: &RegionInput, p: &Params, times: Option<&mut StageTimes>) -> Prepared {
    let mut tm = Timer::new(times);
    let ev = events::variation_events(input.haplotypes, input.ref_padded, input.ref_padded_start, p.max_mnp_distance);
    tm.lap(|t| &mut t.events);

    // ReferenceContext of AssemblyRegionWalker = the region's padded span
    let off = (input.padded_span.start - input.ref_padded_start) as usize;
    let ctx = &input.ref_padded[off..off + (input.padded_span.end - input.padded_span.start + 1) as usize];
    let trimmed = trim::trim_region(input.span, input.padded_span, &ev, ctx, input.padded_span.start, p);
    tm.lap(|t| &mut t.trim);

    let mut out = Prepared {
        exit: Exit::NoVariationAfterTrim,
        variant_span: None,
        variant_padded_span: None,
        assembly_variation_present: false,
        haplotypes: Vec::new(),
        reads: Vec::new(),
        filtered_reads: Vec::new(),
        hmm: Vec::new(),
        gcp: Vec::new(),
    };
    let (vspan, vpadded) = match trimmed {
        Some(x) => x,
        None => {
            if !p.disable_optimizations {
                return out;
            }
            // GATK continues to trimTo(trimmingResult.getVariantRegion()), which throws
            // "There is no variation present."
            panic!("--disable-optimizations with no variation: GATK throws here as well");
        }
    };
    // AssemblyRegion.trim(span, paddedSpan): intersect with the original spans
    let vspan = trim::intersect(input.span, vspan);
    let vpadded = trim::intersect(input.padded_span, vpadded);
    out.variant_span = Some(vspan);
    out.variant_padded_span = Some(vpadded);

    out.haplotypes = trim::trim_haplotypes(input.haplotypes, vpadded);
    // trimTo: variationPresent comes from the *untrimmed* set; isVariationPresent also needs > 1 haplotype
    out.assembly_variation_present = input.haplotypes.iter().any(|h| !h.is_ref) && out.haplotypes.len() > 1;
    tm.lap(|t| &mut t.trim_haplotypes);

    let mut scratch = Vec::new();
    let mut reads: Vec<(Read, bool)> = Vec::with_capacity(input.reads.len());
    for (i, r) in input.reads.iter().enumerate() {
        let mut c = r.clone();
        if clip::hard_clip_to_region(&mut c, vpadded.start, vpadded.end, &mut scratch)
            && !c.seq.is_empty()
            && clip::coord_overlaps(vpadded.start, vpadded.end, c.start, c.end)
        {
            let other = input.mate_on_other_contig.get(i).copied().unwrap_or_else(|| r.mate_on_other_contig());
            reads.push((c, other));
        }
    }
    reads.sort_by(|a, b| clip::read_coordinate_cmp(&a.0, &b.0)); // stable, as List.sort
    for (r, mate_other) in reads {
        if p.mutect {
            if r.seq.len() >= p.min_read_length_after_trimming {
                out.reads.push(r);
            }
            continue;
        }
        let unclipped_len = r.seq.len() - r.cigar.iter().filter(|e| e.0 == b'S').map(|e| e.1 as usize).sum::<usize>();
        if unclipped_len < p.min_read_length_after_trimming {
            continue; // read stub
        }
        let paired_mate_mapped = r.flag & 0x1 != 0 && !clip::mate_is_unmapped(&r);
        if unclipped_len < p.read_length_filter_threshold || r.mapq < p.mapping_quality_threshold || (paired_mate_mapped && mate_other) {
            out.filtered_reads.push(r);
        } else {
            out.reads.push(r);
        }
    }
    tm.lap(|t| &mut t.genotyping_reads);

    out.exit = if !out.assembly_variation_present && !p.disable_optimizations {
        Exit::NoAssemblyVariation
    } else if out.reads.is_empty() && !p.disable_optimizations && !p.mutect {
        Exit::NoReadsForGenotyping
    } else {
        Exit::Genotyped
    };
    if out.exit == Exit::Genotyped {
        let pcr = hmm::pcr_indel_table(p.pcr_error_model);
        let mut rt = hmm::RepeatTables::default();
        let mut rl = Vec::new();
        out.hmm = out
            .reads
            .iter()
            .map(|r| {
                let mut q = hmm::HmmQuals::default();
                hmm::read_hmm_quals_with(&r.seq, &r.qual, &r.bi, &r.bd, r.mapq, pcr.as_ref(), p, &mut q, &mut rt, &mut rl);
                q
            })
            .collect();
        let maxlen = out.reads.iter().map(|r| r.seq.len()).max().unwrap_or(0);
        out.gcp = vec![p.gcp; maxlen];
    }
    tm.lap(|t| &mut t.hmm_inputs);
    out
}

/// Runs the CUDA PairHMM for many prepared regions in one launch (regions that do not need the
/// PairHMM get an empty matrix).  `hc_cudahmm::init` must have succeeded.
///
/// A region with a haplotype longer than the GPU's shared-memory limit (`hc_cudahmm::max_hap_len`,
/// 4,223 bases) goes to the bit-identical CPU PairHMM instead.
pub fn compute_likelihoods_gpu(preps: &[&Prepared]) -> Vec<Vec<f64>> {
    let limit = hc_cudahmm::max_hap_len();
    let fits = |g: &HmmGroup| g.haps.iter().all(|h| h.len() <= limit);
    let groups: Vec<HmmGroup> = preps.iter().filter(|p| p.needs_hmm()).map(|p| p.hmm_group()).collect();
    let (gpu_idx, cpu_idx): (Vec<usize>, Vec<usize>) = (0..groups.len()).partition(|&i| fits(&groups[i]));
    let mut out: Vec<Vec<f64>> = vec![Vec::new(); groups.len()];
    if cpu_idx.is_empty() {
        out = hc_cudahmm::compute(&groups);
    } else {
        let (g, c): (Vec<HmmGroup>, Vec<HmmGroup>) = {
            let mut v: Vec<Option<HmmGroup>> = groups.into_iter().map(Some).collect();
            (gpu_idx.iter().map(|&i| v[i].take().unwrap()).collect(), cpu_idx.iter().map(|&i| v[i].take().unwrap()).collect())
        };
        for (i, r) in gpu_idx.iter().zip(hc_cudahmm::compute(&g)) {
            out[*i] = r;
        }
        for g in &c {
            let h = g.haps.iter().map(|h| h.len()).max().unwrap_or(0);
            eprintln!("[pairhmm] region with a {h}-base haplotype ({} reads x {} haplotypes) computed on the CPU", g.reads.len(), g.haps.len());
        }
        for (i, r) in cpu_idx.iter().zip(cpu::compute_batched(&c, cpu::Flavor::Avx512)) {
            out[*i] = r;
        }
    }
    let mut res = out.into_iter();
    preps.iter().map(|p| if p.needs_hmm() { res.next().unwrap() } else { Vec::new() }).collect()
}

/// CPU fallback for [`compute_likelihoods_gpu`]: GKL-compatible scalar PairHMM (see [`cpu`]).
pub fn compute_likelihoods_cpu(preps: &[&Prepared], flavor: cpu::Flavor) -> Vec<Vec<f64>> {
    let groups: Vec<HmmGroup> = preps.iter().filter(|p| p.needs_hmm()).map(|p| p.hmm_group()).collect();
    let mut res = cpu::compute_batched(&groups, flavor).into_iter();
    preps.iter().map(|p| if p.needs_hmm() { res.next().unwrap() } else { Vec::new() }).collect()
}

/// Phase 2 with default parameters.
pub fn finish(prep: Prepared, raw: &[f64]) -> RegionResult {
    finish_with(prep, raw, &Params::default(), None)
}

/// Phase 2: turns the raw PairHMM matrix (row-major read x haplotype, as returned by
/// `hc_cudahmm::compute` for `prep.hmm_group()`) into the engine's `AlleleLikelihoods`:
/// transposition to `[haplotype][read]`, `normalizeLikelihoods` and `filterPoorlyModeledEvidence`.
pub fn finish_with(prep: Prepared, raw: &[f64], p: &Params, times: Option<&mut StageTimes>) -> RegionResult {
    let mut tm = Timer::new(times);
    let Prepared { exit, variant_span, variant_padded_span, assembly_variation_present, haplotypes, reads, filtered_reads, hmm, .. } = prep;
    if exit != Exit::Genotyped {
        return RegionResult { exit, variant_span, variant_padded_span, assembly_variation_present, filtered_reads, likelihoods: None, haplotypes, reads };
    }
    let nh = haplotypes.len();
    let nr = reads.len();
    assert_eq!(raw.len(), nh * nr, "PairHMM matrix has the wrong size");

    // read-major working copy: v[r * nh + h]
    let mut v = raw.to_vec();
    // normalizeLikelihoods(log10GlobalReadMismappingRate, symmetricallyNormalizeAllelesToReference)
    let cap = p.log10_global_read_mismapping_rate();
    let ref_idx = haplotypes.iter().position(|h| h.is_ref);
    if cap != f64::NEG_INFINITY && nh > 1 {
        let can_be_ref = p.symmetrically_normalize_alleles_to_reference;
        for r in 0..nr {
            let row = &mut v[r * nh..(r + 1) * nh];
            // searchBestAllele(canBeReference)
            let mut best_i = if can_be_ref || ref_idx != Some(0) { 0 } else { 1 };
            let mut best = row[best_i];
            for (a, &x) in row.iter().enumerate().skip(best_i + 1) {
                if !can_be_ref && ref_idx == Some(a) {
                    continue;
                }
                if x > best {
                    best = x;
                    best_i = a;
                }
            }
            let _ = best_i;
            let worst = best + cap;
            for x in row.iter_mut() {
                if *x < worst {
                    *x = worst;
                }
            }
        }
    }

    // filterPoorlyModeledEvidence(log10MinTrueLikelihood(expectedErrorRatePerBase, capLikelihoods=true))
    let mut keep = vec![true; nr];
    for r in 0..nr {
        let row = &v[r * nh..(r + 1) * nh];
        let mut maxv = f64::NEG_INFINITY;
        for &x in row {
            if x > maxv {
                maxv = x;
            }
        }
        let thr = min_true_likelihood(&hmm[r].quals, p);
        if maxv < thr {
            keep[r] = false;
        }
    }
    let n_keep = keep.iter().filter(|&&k| k).count();
    let mut values = vec![0f64; nh * n_keep];
    let mut kept = Vec::with_capacity(n_keep);
    let mut poor = Vec::new();
    let mut j = 0;
    for (r, read) in reads.into_iter().enumerate() {
        if keep[r] {
            for h in 0..nh {
                values[h * n_keep + j] = v[r * nh + h];
            }
            kept.push(read);
            j += 1;
        } else {
            poor.push(read);
        }
    }
    tm.lap(|t| &mut t.finish);
    RegionResult {
        exit,
        variant_span,
        variant_padded_span,
        assembly_variation_present,
        filtered_reads,
        likelihoods: Some(ReadLikelihoods { haplotypes, reads: kept, filtered_reads: poor, values }),
        haplotypes: Vec::new(),
        reads: Vec::new(),
    }
}

/// The per-read threshold of `filterPoorlyModeledEvidence`: `log10MinTrueLikelihood` (static:
/// `min(2, ceil(len * rate)) * -4`) or, with dynamic disqualification, the smaller of that
/// (uncapped) and `calculateLog10DynamicReadQualThreshold`.  `len` is the length of the HMM
/// base-quality array (= read length).
pub fn min_true_likelihood(hmm_quals: &[u8], p: &Params) -> f64 {
    let len = hmm_quals.len() as f64;
    if !p.dynamic_read_disqualification {
        let max_err = (2.0f64).min((len * p.expected_error_rate_per_base).ceil());
        max_err * -4.0
    } else {
        let max_err = (len * p.expected_error_rate_per_base).ceil();
        let static_thr = max_err * -4.0;
        let dynamic = dynamic_read_qual_threshold(hmm_quals, p.read_disqualification_scale);
        if dynamic < static_thr {
            dynamic
        } else {
            static_thr
        }
    }
}

/// `ReadLikelihoodCalculationEngine.calculateLog10DynamicReadQualThreshold`.
fn dynamic_read_qual_threshold(quals: &[u8], k: f64) -> f64 {
    #[rustfmt::skip]
    const T: [f64; 80] = [
        5.996842844, 0.196616587, 5.870018422, 1.388545569, 5.401558531, 5.641990128, 4.818940919, 10.33176216,
        4.218758304, 14.25799688, 3.646319832, 17.02880749, 3.122346753, 18.64537883, 2.654731979, 19.27521677,
        2.244479156, 19.13584613, 1.88893867, 18.43922003, 1.583645342, 17.36842261, 1.3233807, 16.07088712,
        1.102785365, 14.65952563, 0.916703025, 13.21718577, 0.760361881, 11.80207947, 0.629457387, 10.45304833,
        0.520175654, 9.194183767, 0.42918208, 8.038657241, 0.353590663, 6.991779595, 0.290923699, 6.053379213,
        0.23906788, 5.219610436, 0.196230431, 4.484302033, 0.160897421, 3.839943445, 0.131795374, 3.27839108,
        0.1078567, 2.791361596, 0.088189063, 2.370765375, 0.072048567, 2.008921719, 0.058816518, 1.698687797,
        0.047979438, 1.433525748, 0.039111985, 1.207526336, 0.031862437, 1.015402928, 0.025940415, 0.852465956,
        0.021106532, 0.714585285, 0.017163711, 0.598145851, 0.013949904, 0.500000349, 0.011332027, 0.41742159,
        0.009200898, 0.348056286, 0.007467036, 0.289881373, 0.006057179, 0.241163527, 0.004911394, 0.200422214,
    ];
    let mut sum_mean = 0.0;
    let mut sum_var = 0.0;
    for &q in quals {
        let bq = q as usize;
        let e = if bq <= 1 { 0 } else { bq.min(40) - 1 };
        sum_mean += T[e * 2];
        sum_var += T[e * 2 + 1];
    }
    let threshold = sum_mean + k * sum_var.sqrt();
    threshold * -0.1
}

/// CIGAR to SAM string.
pub fn cigar_string(c: &[(u8, u32)]) -> String {
    hc_assembly::cigar::to_string(c)
}
