//! Rust port of GATK 4.7 HaplotypeCaller local assembly:
//! `AssemblyBasedCallerUtils.assembleReads` (after `finalizeRegion`) ->
//! `ReadThreadingAssembler.runLocalAssembly` with default HaplotypeCaller arguments.
//!
//! See README.md for the list of ported components and replicated GATK quirks.

pub mod cigar;
pub mod cigar_utils;
pub mod graph;
pub mod fastmath;
pub mod hotspot;
pub mod jmath;
mod fdlibm; // Math.log1p / expm1, re-exported by jmath
mod hotspot_tables;
pub mod m2math;
pub mod kbest;
pub mod kmertable;
pub mod log10;
pub mod prof;
pub mod rtgraph;
pub mod seqgraph;
pub mod sw;

use cigar::Cigar;
use rtgraph::{ReadThreadingGraph, SeqForKmers};
use seqgraph::SeqGraph;
use sw::{Overhang, SwAligner, SwParameters};

/// A read as handed to the assembler (i.e. after `finalizeRegion`).
#[derive(Clone, Debug, Default)]
pub struct Read {
    pub name: String,
    pub flag: u16,
    /// 1-based alignment start
    pub start: i64,
    /// 1-based inclusive alignment end
    pub end: i64,
    pub cigar: Cigar,
    pub seq: Vec<u8>,
    /// base qualities as raw phred values (NOT +33)
    pub qual: Vec<u8>,
    pub mapq: u8,
    /// mate alignment start, 0 if unpaired or the mate is unmapped (GATK `getMateStart()`)
    pub mate_start: i64,
    pub tlen: i64,
    /// BAM reference index of the read and of its mate (RNEXT). Only used for
    /// `ReadCoordinateComparator` ties and `MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE`; the default 0 / 0
    /// means "mate on the read's contig".
    pub ref_id: i32,
    pub mate_ref_id: i32,
    /// BI / BD base insertion / deletion qualities (raw phred); empty = absent (GATK default 45 per base)
    pub bi: Vec<u8>,
    pub bd: Vec<u8>,
    /// sample index (Mutect2; 0 for single-sample callers)
    pub sample: u16,
    /// Mutect2 `OriginalAlignment`: the read's OA tag names another contig
    pub oa_other_contig: bool,
    /// the original BAM record (bytes after block_size; Mutect2 `--bam-output` only)
    pub bam_rec: Option<std::sync::Arc<[u8]>>,
}

impl Read {
    /// GATK `mateIsUnmapped()` (meaningful for paired reads).
    #[inline]
    pub fn mate_is_unmapped(&self) -> bool {
        self.flag & 0x8 != 0 || self.mate_start == 0
    }
    /// `ReadUtils.getMateReferenceIndex`: -1 when the mate is unmapped.
    #[inline]
    pub fn mate_reference_index(&self) -> i32 {
        if self.mate_is_unmapped() { -1 } else { self.mate_ref_id }
    }
    /// Fails `ReadFilterLibrary.MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE` (paired, mate mapped, mate on
    /// another contig; the read itself is mapped).
    #[inline]
    pub fn mate_on_other_contig(&self) -> bool {
        self.flag & 0x1 != 0 && !self.mate_is_unmapped() && self.mate_ref_id != self.ref_id
    }
    /// `ReadUtils.getSAMFlagsForRead`: the SAM flags GATK rebuilds from the read (used by
    /// `ReadCoordinateComparator`).
    #[inline]
    pub fn gatk_sam_flags(&self) -> u16 {
        let flag = self.flag;
        let mut f = 0u16;
        if flag & 0x1 != 0 {
            f |= 0x1 | (flag & 0x2) | (flag & 0xc0);
            if self.mate_is_unmapped() { f |= 0x8; } else { f |= flag & 0x20; }
        }
        if flag & 0x4 != 0 { f |= 0x4; } else { f |= flag & 0x10; }
        f | (flag & 0xf00)
    }
}

/// The assembly region: active span and padded span (1-based, inclusive).
#[derive(Clone, Copy, Debug, Default)]
pub struct RegionInfo {
    pub start: i64,
    pub end: i64,
    pub padded_start: i64,
    pub padded_end: i64,
}

#[derive(Clone, Debug)]
pub struct Haplotype {
    pub bases: Vec<u8>,
    pub cigar: Cigar,
    /// getAlignmentStartHapwrtRef
    pub align_start_wrt_ref: i64,
    /// getGenomeLocation() start/end (1-based inclusive); -1 if null
    pub loc_start: i64,
    pub loc_end: i64,
    pub is_ref: bool,
    /// NaN when GATK's score is NaN (never set); may be -inf
    pub score: f64,
    pub kmer: u32,
}

#[derive(Clone, Debug)]
pub struct AssemblyResult {
    pub ref_haplotype: Haplotype,
    /// AssemblyResultSet.getHaplotypeList() in GATK order (contains the reference haplotype)
    pub haplotypes: Vec<Haplotype>,
}

/// Assembler parameters. `Default` gives GATK HaplotypeCaller defaults; only the default code
/// path is implemented (non-default flags are accepted but some combinations are unsupported).
#[derive(Clone, Debug)]
pub struct AssemblerConfig {
    pub kmer_sizes: Vec<usize>,
    pub dont_increase_kmer_sizes_for_cycles: bool,
    pub allow_non_unique_kmers_in_ref: bool,
    pub max_num_haplotypes_in_population: usize,
    pub min_prune_factor: i32,
    pub min_base_quality_to_use_in_assembly: u8,
    pub min_matching_bases_to_dangling_end_recovery: i32,
    pub recover_dangling_branches: bool,
    pub dangling_end_sw: SwParameters,
    pub haplotype_to_reference_sw: SwParameters,
    /// Mutect2: `AdaptiveChainPruner` instead of `LowWeightChainPruner(min_prune_factor)`
    pub adaptive_pruning: Option<rtgraph::AdaptivePruning>,
    /// `--recover-all-dangling-branches`
    pub recover_all_dangling_branches: bool,
}

impl Default for AssemblerConfig {
    fn default() -> Self {
        AssemblerConfig {
            kmer_sizes: vec![10, 25],
            dont_increase_kmer_sizes_for_cycles: false,
            allow_non_unique_kmers_in_ref: false,
            max_num_haplotypes_in_population: 128,
            min_prune_factor: 2,
            min_base_quality_to_use_in_assembly: 10,
            min_matching_bases_to_dangling_end_recovery: -1,
            recover_dangling_branches: true,
            dangling_end_sw: sw::STANDARD_NGS,
            haplotype_to_reference_sw: sw::NEW_SW_PARAMETERS,
            adaptive_pruning: None,
            recover_all_dangling_branches: false,
        }
    }
}

const KMER_SIZE_ITERATION_INCREASE: usize = 10;
const MAX_KMER_ITERATIONS_TO_ATTEMPT: usize = 6;
const MIN_HAPLOTYPE_REFERENCE_LENGTH: i64 = 30;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum AssemblyStatus {
    Failed,
    JustAssembledReference,
    AssembledSomeVariation,
}

struct GraphResult {
    status: AssemblyStatus,
    seq_graph: Option<SeqGraph>,
    kmer: usize,
}

/// Reusable assembler (holds SW buffers). Not thread-safe; use one per thread.
pub struct Assembler {
    pub config: AssemblerConfig,
    aligner: SwAligner,
    sw_cache: sw::FixedRefAligner,
    cache_ref_set: bool,
    pad_buf: (Vec<u8>, Vec<u8>),
    m2math: Option<Box<m2math::M2Math>>,
}

impl Default for Assembler {
    fn default() -> Self {
        Self::new(AssemblerConfig::default())
    }
}

/// Convenience wrapper using default parameters.
pub fn assemble_region(
    reads: &[Read],
    region: RegionInfo,
    ref_padded: &[u8],
    ref_padded_start: i64,
    ref_padded_end: i64,
) -> AssemblyResult {
    Assembler::default().assemble_region(reads, region, ref_padded, ref_padded_start, ref_padded_end)
}

/// ReadClipper.hardClipSoftClippedBases: returns (bases, quals) without soft-clipped bases.
fn hard_clip_soft_clipped<'a>(r: &'a Read) -> (&'a [u8], &'a [u8]) {
    if r.seq.is_empty() {
        return (&r.seq, &r.qual);
    }
    let mut read_index: i64 = 0;
    let mut cut_left: i64 = -1;
    let mut cut_right: i64 = -1;
    let mut right_tail = false;
    for &(op, len) in &r.cigar {
        if op == b'S' {
            if right_tail {
                cut_right = read_index;
            } else {
                cut_left = read_index + len as i64 - 1;
            }
        } else if op != b'H' {
            right_tail = true;
        }
        if cigar::consumes_read(op) {
            read_index += len as i64;
        }
    }
    let mut end = r.seq.len();
    let mut start = 0usize;
    if cut_right >= 0 {
        end = cut_right as usize;
    }
    if cut_left >= 0 {
        start = (cut_left + 1) as usize;
    }
    if start >= end {
        return (&r.seq[0..0], &r.qual[0..0]);
    }
    (&r.seq[start..end], &r.qual[start..end])
}

impl Assembler {
    pub fn new(config: AssemblerConfig) -> Self {
        Assembler {
            config,
            aligner: SwAligner::new(),
            sw_cache: sw::FixedRefAligner::new(),
            cache_ref_set: false,
            pad_buf: (Vec::new(), Vec::new()),
            m2math: None,
        }
    }

    /// Like `assemble_region`, but returns `Err(message)` in the situations where GATK throws an
    /// exception (e.g. "Graph must have ref source and sink vertices"), which aborts a GATK run.
    /// While it runs, [`panic_is_caught`] is true on this thread.
    pub fn try_assemble_region(
        &mut self,
        reads: &[Read],
        region: RegionInfo,
        ref_padded: &[u8],
        ref_padded_start: i64,
        ref_padded_end: i64,
    ) -> Result<AssemblyResult, String> {
        GUARDED.with(|g| g.set(true));
        let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            // test hook (VARIANT_RS_TEST_ASSEMBLY_FAIL=POS): assembly fails for every region starting at or after
            // POS, to test what a run does with a failed region (it must fail, as GATK does)
            static FAIL_AT: std::sync::OnceLock<Option<i64>> = std::sync::OnceLock::new();
            if let Some(p) = *FAIL_AT.get_or_init(|| std::env::var("VARIANT_RS_TEST_ASSEMBLY_FAIL").ok().and_then(|v| v.parse().ok())) {
                if region.start >= p {
                    panic!("test assembly failure (VARIANT_RS_TEST_ASSEMBLY_FAIL={p})");
                }
            }
            self.assemble_region(reads, region, ref_padded, ref_padded_start, ref_padded_end)
        }));
        GUARDED.with(|g| g.set(false));
        r.map_err(|e| {
            if let Some(s) = e.downcast_ref::<&str>() {
                s.to_string()
            } else if let Some(s) = e.downcast_ref::<String>() {
                s.clone()
            } else {
                "assembly failed".to_string()
            }
        })
    }

    /// AssemblyBasedCallerUtils.assembleReads (after finalizeRegion) + ReadThreadingAssembler.runLocalAssembly.
    /// Panics where GATK would throw; see `try_assemble_region`.
    pub fn assemble_region(
        &mut self,
        reads: &[Read],
        region: RegionInfo,
        ref_padded: &[u8],
        ref_padded_start: i64,
        ref_padded_end: i64,
    ) -> AssemblyResult {
        assert_eq!((ref_padded_end - ref_padded_start + 1) as usize, ref_padded.len(), "Reference bases and reference loc must be the same size.");
        // createReferenceHaplotype
        let a = (region.padded_start - ref_padded_start) as usize;
        let b = (region.padded_end - ref_padded_start + 1) as usize;
        let ref_bases = ref_padded[a..b].to_vec();
        let ref_hap = Haplotype {
            cigar: vec![(b'M', ref_bases.len() as u32)],
            bases: ref_bases,
            align_start_wrt_ref: region.padded_start - ref_padded_start,
            loc_start: region.padded_start,
            loc_end: region.padded_end,
            is_ref: true,
            score: f64::NAN,
            kmer: 0,
        };
        self.cache_ref_set = false;
        let _p = prof::scope("assembly total");
        let haplotypes = self.run_local_assembly(reads, &ref_hap, region);
        drop(_p);
        AssemblyResult { ref_haplotype: ref_hap, haplotypes }
    }

    fn run_local_assembly(&mut self, reads: &[Read], ref_hap: &Haplotype, region: RegionInfo) -> Vec<Haplotype> {
        // Build the shared arena: ref bases followed by all read bases (after hard-clipping soft clips)
        let _pa = prof::scope("arena");
        let mut arena: Vec<u8> = Vec::with_capacity(ref_hap.bases.len() + reads.iter().map(|r| r.seq.len()).sum::<usize>());
        arena.extend_from_slice(&ref_hap.bases);
        let mut read_spans: Vec<(u32, u32, usize)> = Vec::with_capacity(reads.len()); // (arena off, len, sample group)
        let mut quals: Vec<&[u8]> = Vec::with_capacity(reads.len());
        for r in reads {
            let (bases, q) = hard_clip_soft_clipped(r);
            read_spans.push((arena.len() as u32, bases.len() as u32, r.sample as usize));
            quals.push(q);
            arena.extend_from_slice(bases);
        }

        drop(_pa);
        // assemble(): requested kmer sizes, then larger ones if none worked
        let _pc = prof::scope("create_graph (all k)");
        let mut kmer_sizes = self.config.kmer_sizes.clone();
        kmer_sizes.sort();
        let mut results: Vec<GraphResult> = Vec::new();
        for &k in &kmer_sizes {
            let dont_inc = self.config.dont_increase_kmer_sizes_for_cycles;
            let allow_nu = self.config.allow_non_unique_kmers_in_ref;
            if let Some(r) = self.create_graph(&arena, &read_spans, &quals, ref_hap, k, dont_inc, allow_nu) {
                results.push(r);
            }
        }
        if results.is_empty() && !self.config.dont_increase_kmer_sizes_for_cycles {
            let mut k = *kmer_sizes.iter().max().unwrap() + KMER_SIZE_ITERATION_INCREASE;
            let mut num_iterations = 1;
            while results.is_empty() && num_iterations <= MAX_KMER_ITERATIONS_TO_ATTEMPT {
                let last_attempt = num_iterations == MAX_KMER_ITERATIONS_TO_ATTEMPT;
                if let Some(r) = self.create_graph(&arena, &read_spans, &quals, ref_hap, k, last_attempt, last_attempt) {
                    results.push(r);
                }
                k += KMER_SIZE_ITERATION_INCREASE;
                num_iterations += 1;
            }
        }

        drop(_pc);
        let _pf = prof::scope("find paths + cigars (all graphs)");
        // findBestPaths over all graphs with variation
        let active_region_start = ref_hap.align_start_wrt_ref;
        let ref_len = cigar::reference_length(&ref_hap.cigar);
        // returnHaplotypes: LinkedHashSet<Haplotype> keyed by (bases, isRef)
        let mut ret: Vec<Haplotype> = Vec::new();
        let mut ret_index: std::collections::HashMap<(Vec<u8>, bool), usize> = std::collections::HashMap::new();
        let mut result_set: Option<Vec<Haplotype>> = None;
        for r in results.iter() {
            if r.status != AssemblyStatus::AssembledSomeVariation {
                continue;
            }
            let sg = r.seq_graph.as_ref().unwrap();
            let source = sg.g.reference_source_vertex().expect("source");
            let sink = sg.g.reference_sink_vertex().expect("sink");
            let _gk = prof::scope("kbest");
            let best = kbest::find_best_haplotypes(sg, source, sink, self.config.max_num_haplotypes_in_population);
            drop(_gk);
            let _gs = prof::scope("hap_cigar");
            // calculateCigar is a pure function of the bases: compute all needed SOFTCLIP cigars
            // up front, sorted (duplicates adjacent for the dedup; the aligner reuses the columns of shared prefixes).
            let mut todo: Vec<usize> = (0..best.len())
                .filter(|&i| !ret_index.contains_key(&(best[i].bases.clone(), false)))
                .collect();
            todo.sort_by(|&a, &b| best[a].bases.cmp(&best[b].bases));
            todo.dedup_by(|a, b| best[*a].bases == best[*b].bases);
            if !todo.is_empty() && prof::enabled() {
                let cells: u64 = todo.iter().map(|&i| (best[i].bases.len() * ref_hap.bases.len()) as u64).sum();
                let b = match todo.len() { 1 => "1", 2..=4 => "2-4", 5..=8 => "5-8", 9..=16 => "9-16", 17..=32 => "17-32", 33..=64 => "33-64", _ => "65+" };
                let name: &'static str = match b { "1" => "hapSW cells by haps per graph: 1", "2-4" => "hapSW cells by haps per graph: 2-4", "5-8" => "hapSW cells by haps per graph: 5-8", "9-16" => "hapSW cells by haps per graph: 9-16", "17-32" => "hapSW cells by haps per graph: 17-32", "33-64" => "hapSW cells by haps per graph: 33-64", _ => "hapSW cells by haps per graph: 65+" };
                prof::add(name, cells);
            }
            let mut softclip_cigars: std::collections::HashMap<Vec<u8>, Option<Cigar>> =
                std::collections::HashMap::with_capacity(todo.len());
            if !todo.is_empty() {
                if !self.cache_ref_set {
                    cigar_utils::pad(&ref_hap.bases, &mut self.pad_buf.0);
                    self.sw_cache.set_reference(&self.pad_buf.0, &self.config.haplotype_to_reference_sw);
                    self.cache_ref_set = true;
                }
                self.sw_cache.begin_batch();
                let alts: Vec<&[u8]> = todo.iter().map(|&i| &best[i].bases[..]).collect();
                let cigars = cigar_utils::calculate_cigars_softclip(&mut self.sw_cache, &ref_hap.bases, &alts);
                for (&i, c) in todo.iter().zip(cigars) {
                    softclip_cigars.insert(best[i].bases.clone(), c);
                }
            }
            for kb in best {
                // KBestHaplotype.isReference() is always false in GATK (bug), so haplotypes are non-ref
                let key = (kb.bases, false);
                if ret_index.contains_key(&key) {
                    continue;
                }
                let bases = key.0;
                let cig = softclip_cigars.get(&bases).expect("precomputed cigar").clone();
                let cig = match cig {
                    None => {
                        prof::hit("hap cigar SW failure");
                        continue; // failedCigars++
                    }
                    Some(c) => c,
                };
                assert!(!cig.is_empty(), "Smith-Waterman alignment failure");
                let rl = cigar::reference_length(&cig);
                if cig.iter().any(|e| e.0 == b'N') || rl < MIN_HAPLOTYPE_REFERENCE_LENGTH {
                    prof::hit("hap cigar too short");
                    continue;
                } else if rl != ref_len {
                    let c2 = cigar_utils::calculate_cigar(
                        &mut self.aligner,
                        &ref_hap.bases,
                        &bases,
                        &self.config.haplotype_to_reference_sw,
                        Overhang::Indel,
                        &mut self.pad_buf,
                    );
                    let ok = c2.map(|c| cigar::reference_length(&c) == ref_len).unwrap_or(false);
                    if ok {
                        prof::hit("hap cigar ref-length mismatch (INDEL ok)");
                        continue;
                    } else {
                        panic!("Smith-Waterman alignment failure (GATK throws IllegalStateException)");
                    }
                }
                let norm = cigar::haplotype_normalize(&cig);
                assert_eq!(cigar::read_length(&norm), bases.len() as i64);
                let h = Haplotype {
                    bases: bases.clone(),
                    cigar: norm,
                    align_start_wrt_ref: active_region_start,
                    loc_start: region.padded_start,
                    loc_end: region.padded_end,
                    is_ref: false,
                    score: kb.score,
                    kmer: r.kmer as u32,
                };
                ret_index.insert((bases, false), ret.len());
                ret.push(h);
            }
            if !ret.is_empty() {
                // remove the non-ref haplotype with reference bases, then add the reference haplotype
                let tmp_key = (ref_hap.bases.clone(), false);
                if let Some(idx) = ret_index.remove(&tmp_key) {
                    ret.remove(idx);
                    for v in ret_index.values_mut() {
                        if *v > idx {
                            *v -= 1;
                        }
                    }
                }
                let ref_key = (ref_hap.bases.clone(), true);
                if !ret_index.contains_key(&ref_key) {
                    ret_index.insert(ref_key, ret.len());
                    ret.push(ref_hap.clone());
                }
                result_set = Some(ret.clone());
            }
        }
        match result_set {
            Some(v) => v,
            None => vec![ref_hap.clone()],
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn create_graph(
        &mut self,
        arena: &[u8],
        read_spans: &[(u32, u32, usize)],
        quals: &[&[u8]],
        ref_hap: &Haplotype,
        k: usize,
        allow_low_complexity: bool,
        allow_non_unique_in_ref: bool,
    ) -> Option<GraphResult> {
        let ref_len = ref_hap.bases.len();
        if ref_len < k {
            prof::hit("graph FAILED: ref shorter than k");
            return Some(GraphResult { status: AssemblyStatus::Failed, seq_graph: None, kmer: k });
        }
        let ref_seq = SeqForKmers { seq_off: 0, seq_len: ref_len as u32, start: 0, stop: ref_len as u32, is_ref: true };
        let _p0 = prof::scope("ref non-unique check + read segments");
        if !allow_non_unique_in_ref {
            let mut seen = kmertable::KmerTable::new(k, ref_len);
            let mut any = false;
            rtgraph::for_each_non_unique_kmer(arena, &ref_seq, k, &mut seen, |_| any = true);
            if any {
                prof::hit("graph rejected: non-unique ref kmer");
                return None;
            }
        }
        // addSequence("ref") then addRead for each read
        // ReadThreadingGraph.pending: one list per sample, in the order of each sample's first added sequence
        let mut group_seqs: Vec<(usize, Vec<SeqForKmers>)> = Vec::new();
        let minq = self.config.min_base_quality_to_use_in_assembly;
        let mut segs: Vec<(usize, SeqForKmers)> = Vec::new();
        for (i, &(off, len, sample)) in read_spans.iter().enumerate() {
            segs.clear();
            let q = quals[i];
            let seq = &arena[off as usize..(off + len) as usize];
            let mut last_good: i64 = -1;
            for end in 0..=len as usize {
                if end == len as usize || seq[end] == b'N' || q[end] < minq {
                    let start = last_good;
                    let l = end as i64 - start;
                    if start != -1 && l >= k as i64 {
                        segs.push((sample, SeqForKmers {
                            seq_off: off,
                            seq_len: len,
                            start: start as u32,
                            stop: end as u32,
                            is_ref: false,
                        }));
                    }
                    last_good = -1;
                } else if last_good == -1 {
                    last_good = end as i64;
                }
            }
            for (sample, sq) in segs.drain(..) {
                match group_seqs.iter_mut().find(|g| g.0 == sample) {
                    Some(g) => g.1.push(sq),
                    None => group_seqs.push((sample, vec![sq])),
                }
            }
        }
        let mut pending = vec![vec![ref_seq]];
        for (_, read_seqs) in group_seqs {
            if !read_seqs.is_empty() {
                pending.push(read_seqs);
            }
        }
        drop(_p0);
        let _g1 = prof::scope("rt_build");
        let mut rt = ReadThreadingGraph::new(k, arena.to_vec(), pending, self.config.min_matching_bases_to_dangling_end_recovery);
        rt.build();
        drop(_g1);
        let _g2 = prof::scope("rt_prune_cycles");
        // prune before cycle counting (default: legacy cycle detection disabled)
        rt.prune_factor = self.config.min_prune_factor;
        rt.recover_all = self.config.recover_all_dangling_branches;
        match &self.config.adaptive_pruning {
            Some(ap) => {
                let math = self.m2math.get_or_insert_with(|| Box::new(m2math::M2Math::new()));
                rt.prune_adaptive(ap, math);
            }
            None => rt.prune_low_weight_chains(self.config.min_prune_factor),
        }
        if rt.g.has_cycles() {
            prof::hit("graph rejected: cycles");
            return None;
        }
        if !allow_low_complexity && rt.is_low_quality_graph() {
            prof::hit("graph rejected: low complexity");
            return None;
        }
        if k > 64 {
            prof::hit("graph built with k > 64");
        } else if k > 32 {
            prof::hit("graph built with 32 < k <= 64");
        }
        drop(_g2);
        // getAssemblyResult
        let _g3 = prof::scope("dangling_recovery");
        if self.config.recover_dangling_branches {
            let p = self.config.dangling_end_sw;
            rt.recover_dangling_tails(&mut self.aligner, &p);
            rt.recover_dangling_heads(&mut self.aligner, &p);
        }
        drop(_g3);
        let _g4 = prof::scope("seqgraph");
        rt.remove_paths_not_connected_to_ref();
        // ReadThreadingAssembler.createGraph: "check whether recovering dangling ends created cycles"
        if self.config.recover_all_dangling_branches && rt.g.has_cycles() {
            prof::hit("graph rejected: cycles after recovering all dangling branches");
            return None;
        }
        let mut sg = rt.to_sequence_graph();
        sg.clean_non_ref_paths();
        let (status, sg) = cleanup_seq_graph(sg);
        if status == AssemblyStatus::AssembledSomeVariation {
            // sanityCheckReferenceGraph
            let s = sg.g.reference_source_vertex().expect("All reference graphs must have a reference source vertex.");
            let t = sg.g.reference_sink_vertex().expect("All reference graphs must have a reference sink vertex.");
            let rb = sg.reference_bytes(s, t);
            assert!(rb == ref_hap.bases, "Mismatch between the reference haplotype and the reference assembly graph path");
        }
        Some(GraphResult { status, seq_graph: Some(sg), kmer: k })
    }
}

/// ReadThreadingAssembler.cleanupSeqGraph
pub fn cleanup_seq_graph(mut sg: SeqGraph) -> (AssemblyStatus, SeqGraph) {
    sg.zip_linear_chains();
    sg.remove_singleton_orphan_vertices();
    sg.remove_vertices_not_connected_to_ref_regardless_of_edge_direction();
    sg.simplify_graph();
    if sg.g.reference_source_vertex().is_none() || sg.g.reference_sink_vertex().is_none() {
        prof::hit("seqgraph: JUST_ASSEMBLED_REFERENCE");
        return (AssemblyStatus::JustAssembledReference, sg);
    }
    sg.remove_paths_not_connected_to_ref();
    sg.simplify_graph();
    if sg.g.vertex_count() == 1 {
        prof::hit("seqgraph: single vertex (dummy added)");
        let complete = sg.g.vertex_iter().next().unwrap();
        let dummy = sg.g.add_vertex(Vec::new());
        sg.g.add_edge(complete, dummy, seqgraph::SeqEdge { mult: 0, is_ref: true });
    }
    (AssemblyStatus::AssembledSomeVariation, sg)
}

thread_local! {
    static GUARDED: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// True while this thread is inside `try_assemble_region`, whose panics are expected and caught
/// (GATK's assembly exceptions). A panic hook can use it to treat every other panic as fatal.
pub fn panic_is_caught() -> bool {
    GUARDED.with(|g| g.get())
}
