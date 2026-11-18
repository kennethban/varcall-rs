# hc-likelihood: GATK 4.7 HaplotypeCaller from assembly to read likelihoods, in Rust

This crate ports the part of `HaplotypeCallerEngine.callRegion` that runs between
`AssemblyBasedCallerUtils.assembleReads(...)` and `likelihoodCalculationEngine.computeReadLikelihoods(...)`,
with the default HaplotypeCaller arguments:

1. variation events of the assembled haplotypes (`AssemblyResultSet.getVariationEvents(maxMnpDistance=0)`, `EventMap`, `Event`);
2. `AssemblyRegionTrimmer.trim` (the non-legacy trimmer), including the "no variation" exit;
3. `AssemblyResultSet.trimTo(variantRegion)`: `Haplotype.trim`, merging of haplotypes that become equal, and the sort order;
4. `getRegionForGenotyping()` reads (`AssemblyRegion.trim` → `ReadClipper.hardClipToRegion`, sorted), stub-read removal,
   `filterNonPassingReads`, and the `isVariationPresent` / empty-region exits;
5. `PairHMMLikelihoodCalculationEngine.computeReadLikelihoods`: `modifyReadQualities` (PCR indel model and quality caps),
   the PairHMM batch, the PairHMM itself (CUDA through `hc_cudahmm`, or an optional CPU port of GKL),
   `AlleleLikelihoods.normalizeLikelihoods`, and `filterPoorlyModeledEvidence`.

Everything matches GATK exactly, including the PairHMM likelihoods: both the CUDA kernel and the optional CPU
PairHMM match GATK's GKL AVX-512 output bit for bit.

## Results (NA12878 30x, chr20:10–15 Mb, `$VARCALL_DATA/asm/gatk_assembly.jsonl` + `$VARCALL_DATA/geno/gatk_genotyping.jsonl`)

```
$ target/release/likelihood-validate $VARCALL_DATA/asm/gatk_assembly.jsonl $VARCALL_DATA/geno/gatk_genotyping.jsonl
regions compared: 6345
GATK exits: {"genotyped": 6270, "no_variation_after_trim": 75}
  hmm_reads (bases, qual, ins, del, gcp; batch order)       6270 / 6270
  exit                                                      6345 / 6345
  trim_variation_present                                    6345 / 6345
  variant spans                                             6270 / 6270
  geno spans                                                6270 / 6270
  assembly_variation_present                                6270 / 6270
  trimmed haplotypes (all fields, order)                    6270 / 6270
  filtered_reads (filterNonPassingReads)                    6270 / 6270
  genotyping_reads                                          6270 / 6270
  lk_haplotypes                                             6270 / 6270
  lk_reads (evidence order)                                 6270 / 6270
  lk_filtered_reads                                         6270 / 6270
  lk_matrix shape                                           6270 / 6270
  lk_matrix within 1e-4                                     6270 / 6270

reads: 271620 kept as evidence, 11799 poorly modelled (filterPoorlyModeledEvidence), 0 removed by filterNonPassingReads
lk_matrix: 1004565 values, max |diff| = 0.000e0, > 1e-5: 0, > 1e-4: 0, bit-identical: 1004565
reads kept whose best likelihood is within 1e-4 of the poorly-modelled threshold: 0
```

* All reads, haplotypes, spans, decisions and PairHMM inputs are compared field by field (name, flag,
  start, end, cigar, bases, quals, mapq, mate, tlen; haplotype bases, cigar, alignment start, location,
  ref flag, score, kmer). Order is checked too.
* The tolerance never changes a discrete outcome: the evidence and poorly-modelled sets are identical, and no
  kept read has its best likelihood within 1e-4 of the filtering threshold.
* With `--cpu avx512`, the CPU PairHMM gives `max |diff| = 0`, with all 1,004,565 values bit-identical.
  With `--cpu avx2` (GKL's non-FMA build), max |diff| is 1.9e-6, with 1,004,325 values bit-identical.
* The validator parses the dumps with serde_json's `float_roundtrip` feature. The default serde_json float
  parser is sometimes off by 1 ulp, which would look like a PairHMM mismatch.

### Differential fuzzing against GATK (`Fuzz.java`)

The real data never reaches several code paths: `filterNonPassingReads` removes nothing, the
"no_assembly_variation" and "no_reads_for_genotyping" exits never happen, and there are no compound
SNP+indel events. `scratch/Fuzz.java` (a development probe, not in the repository) generates synthetic regions to cover them:

* a random reference with planted STRs, and regions near both contig ends;
* raw reads run through the real `finalizeRegion`, with MAPQ < 20, short reads, unmapped mates and mates
  on another contig;
* haplotypes with SNPs, insertions, deletions, STR expansions and contractions, SNP+insertion and
  SNP+deletion at the same anchor, and edits near the padded ends.

It then runs GATK's own code for the whole stage, from `getVariationEvents` to `computeReadLikelihoods`
with the AVX PairHMM, and writes both dumps with GATK's dumpers.

```
GATK_ASSEMBLY_DUMP=a.jsonl GATK_GENOTYPING_DUMP=g.jsonl java -cp $GATK_JAR scratch/Fuzz.java <seed> <n> m.jsonl
target/release/likelihood-validate a.jsonl g.jsonl --mates m.jsonl [--cpu avx512]
```

Three seeds were run, 12,000 regions in total. GATK exits: genotyped 9,791, no_variation_after_trim 1,661,
no_reads_for_genotyping 519, no_assembly_variation 29. The reads included 25,555 removed by
`filterNonPassingReads`, 49,313 poorly modelled and 49,720 kept as evidence. Every check matches in every
region. The CUDA max |diff| was 3.8e-6 before the kernel was made bit-exact; the CPU PairHMM is bit-identical. GATK threw no exceptions.

## Timings (full dump, RTX 5070 Ti, Ryzen 9 9950X3D)

| Stage | Time |
|---|---:|
| events (`getVariationEvents`), 1 thread | 0.017 s |
| trimmer, 1 thread | 0.001 s |
| `trimTo` (haplotype trimming), 1 thread | 0.009 s |
| genotyping reads (clip, sort, stub removal, filter), 1 thread | 0.069 s |
| PairHMM inputs (PCR model, caps), 1 thread | 0.29 s |
| **prepare total**: 1 thread / 8 threads | 0.38 s / 0.064 s |
| **PairHMM GPU**: 1,087,302 pairs, 1.39e10 cells, 256 regions per launch (25 launches) | 0.16 s (kernel 0.079 s) |
| finish (normalisation + filtering), 1 thread | 0.024 s |
| PairHMM CPU fallback (`compute_batched`): 1 / 8 / 16 threads | 6.8 s / 1.0 s / 0.79 s |

GPU time by regions per launch: 32 → 0.20 s, 128 → 0.15 s, 256 → 0.16 s, 1024 → 0.16 s, everything in
one launch → 0.18 s. For comparison, GKL's AVX-512 PairHMM takes 4.1 s on one thread for the same batches
(measured while developing the [CUDA PairHMM](../../cuda-pairhmm/native/cuda_pairhmm.cu)). The whole validator run, including 0.7 s of JSON parsing, takes 2.3 s.

## API

```rust
use hc_likelihood::{prepare_with, finish_with, compute_likelihoods_gpu, Params, RegionInput, Span, Exit};
// re-exported: hc_assembly::{Read, Haplotype}, hc_cudahmm::{HmmGroup, HmmRead}

let input = RegionInput {
    span, padded_span,                       // assembly region (1-based inclusive)
    reads: &final_reads,                     // region.getReads() after finalizeRegion, in order
    mate_on_other_contig: &[],               // per read; empty = none (see "Inputs" below)
    haplotypes: &assembly.haplotypes,        // untrimmed getHaplotypeList(), includes the reference
    ref_padded: &ref_padded, ref_padded_start,  // getFullReferenceWithPadding / getPaddedReferenceLoc
};
let p = Params::default();                   // GATK defaults
// phase 1, per region, parallelisable
let prep: Prepared = prepare_with(&input, &p, None);
// prep.exit, prep.variant_span / variant_padded_span (= regionForGenotyping spans),
// prep.haplotypes (trimmed, likelihood allele order), prep.reads (genotyping reads, PairHMM order),
// prep.filtered_reads (filterNonPassingReads), prep.assembly_variation_present
// phase 2: batch many regions per GPU launch (regions with exit != Genotyped get an empty matrix)
hc_cudahmm::init(false);
let raws: Vec<Vec<f64>> = compute_likelihoods_gpu(&[&prep, /* ... */]);
//   or build the groups yourself: prep.hmm_group() -> HmmGroup, then hc_cudahmm::compute(&groups)
//   or the CPU fallback: compute_likelihoods_cpu(&preps, cpu::Flavor::Avx512)
// phase 3, per region
let res: RegionResult = finish_with(prep, &raws[0], &p, None);
// res.likelihoods: Option<ReadLikelihoods { haplotypes, reads (evidence), filtered_reads (poorly modelled),
//                                           values: [hap][read] row-major }>
// res.filtered_reads (filterNonPassingReads); res.haplotypes/res.reads for exits without likelihoods
```

Lower-level functions are public too: `events::{variation_events, haplotype_events}`,
`trim::{trim_region, trim_haplotypes, trim_haplotype, event_tandem_repeat}`, `clip::hard_clip_to_region`,
`hmm::{read_hmm_quals, pcr_indel_table, tandem_repeat_length}`, `min_true_likelihood`, and
`cpu::{compute, compute_batched, pair_log10}`. `StageTimes` collects per-stage timings.

### Inputs

* `Read` is `hc_assembly::Read`: raw phred quals, 1-based `start`/`end`, and `mate_start = 0` when unpaired
  or the mate is unmapped. The reads must be `finalizeRegion`'s output, so they have no soft clips.
* `MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE`: `Read` carries `ref_id` / `mate_ref_id` (from the BAM's
  `refID` / `next_refID`); with `mate_on_other_contig: &[]` the filter uses them (`Read::mate_on_other_contig`).
  A non-empty slice overrides them per read (the fuzz validator's `--mates` input does this). variant-rs relies on
  the read fields; the whole-genome tests in [docs/notes/variant-rs.md](../../../docs/notes/variant-rs.md) exercise it.
* Reads are assumed to have no BQSR `BI`/`BD` tags, so insertion/deletion GOPs start at Q45. This matches
  the dump: every GOP is ≤ 45. With such tags GATK would use them.
* There is one sample. GATK splits reads by sample, then normalises and filters per sample, which is the
  same thing for one sample.

## Defaults checked in the argument collections

| Argument | Default | Effect |
|---|---|---|
| `--max-mnp-distance` (`HaplotypeCallerArgumentCollection.getDefaultMaxMnpDistance`) | 0 | only SNPs, no MNPs |
| `--padding-around-snps` / `-indels` / `-strs` (`AssemblyRegionArgumentCollection`) | 20 / 75 / 75 | trimming |
| `--enable-legacy-assembly-region-trimming` | false | non-legacy `trim` |
| `MINIMUM_READ_LENGTH_AFTER_TRIMMING`, `READ_LENGTH_FILTER_THRESHOLD` | 10, 10 | stub removal and length filter |
| `--mapping-quality-threshold-for-genotyping` | 20 | `filterNonPassingReads` |
| `--keep-rg`, `--sample-name` | null | no read-group or sample filter |
| `--disable-optimizations`, `--dont-genotype` | false | early exits active |
| `--pileup-detection`, `--generate-pd-haplotypes` (pileup args) | false | `goodAndBadPileupEvents` returns empty sets, so `removeHaplotypesWithBadAlleles` and `injectPileupEvents` are no-ops |
| `--alleles` | none | `addGivenAlleles` no-op |
| `--stepwise-filtering`, `--filter-alleles`, `--flow-mode` | false / false / NONE | plain `PairHMMLikelihoodCalculationEngine`, no uncollapsing |
| `--likelihood-calculation-engine` | PairHMM | |
| `--pair-hmm-implementation` | FASTEST_AVAILABLE (GKL AVX-512 here) | |
| `--pair-hmm-gap-continuation-penalty` | 10 | constant GCP |
| `--pcr-indel-model` | CONSERVATIVE (rate factor 3) | PCR indel GOP caps |
| `--base-quality-score-threshold` (`PairHMM.BASE_QUALITY_SCORE_THRESHOLD`) | 18 | quals < 18 become 6 |
| `QualityUtils.MIN_USABLE_Q_SCORE` | 6 | |
| `--disable-cap-base-qualities-to-map-quality` | false | base quals capped at MAPQ |
| `--dragstr-params-path` / `--dont-use-dragstr-pair-hmm-scores` | null / false | `DragstrParamUtils.parse(null)` is null, so the standard imputator is used and the PCR model stays on |
| `--phred-scaled-global-read-mismapping-rate` | 45 | normalisation cap is best + (45 × -0.1) |
| `--disable-symmetric-hmm-normalizing` | false | best allele may be the reference |
| `--expected-mismatch-rate-for-read-disqualification` | 0.02 | poorly-modelled threshold |
| `--enable-dynamic-read-disqualification-for-genotyping` | false | static threshold (the dynamic one is ported but unvalidated) |
| `softClipLowQualityEnds` | false | `handleSoftclips = true`, so reads are not clipped before the PairHMM |

## GATK behaviours replicated

**Events**
1. `EventMap.getEvents`: an insertion is an event only when `refPos > 0` and the `I` element is neither
   the first nor the last cigar element. A deletion is an event only when `refPos > 0`. Both are kept only
   when all their bases are "regular" (`ACGTacgt*`). A mismatch counts only when both bases are regular.
   Event positions use the full padded reference (`ref_padded_start`) and `alignStartHapwrtRef`.
2. Events at the same start are compounded (`makeCompoundEvents`): SNP + insertion gives alt = SNP alt +
   inserted bases; SNP + deletion gives ref = deletion ref and alt = SNP alt. `Event`'s constructor then
   applies the minimal representation (it strips a common suffix when both alleles are longer than 1).
3. The variation-event set is a `TreeSet` ordered by (start, ref length, alt bases).

**Trimming (`AssemblyRegionTrimmer.trim`)**
4. Only events overlapping the *active* span count (`CoordMath.overlaps`). With none, the result is
   "no variation" and `callRegion` exits.
5. Variant span = [min start, max end] ∩ active span. Padding: 20 for SNPs, 75 for indels. For an indel with
   a tandem repeat, padding is 75 + (max repeat count over ref/alt) × unit length. The start is clamped at 1,
   and the padded span is ∩ the region's padded span. The min/max that get padded are the unintersected ones.
6. The tandem repeat comes from `TandemRepeat.getNumTandemRepeatUnits(ReferenceContext, Event)`. The
   `ReferenceContext` window is the region's padded span (`AssemblyRegionWalker`), and the remaining context
   starts one base after the event start. `findRepeatedSubstring` compares zero-padded pieces, so only whole
   repeats match. The alt repeat count can be negative, as in Java, and it still counts as "not 0".
7. `AssemblyRegion.trim(span, paddedSpan)` intersects both spans with the original ones. The genotyping
   spans equal the variant spans.

**`trimTo` / `Haplotype.trim`**
8. `getBasesCoveringRefInterval` returns null, which drops the haplotype, when either end of the new span
   falls in a deletion. `trimCigarByReference` uses a `CigarBuilder(true)` that removes deletions at the
   ends. Then a leading or trailing insertion element is removed. A cigar that is only insertion drops the
   haplotype. `setCigar` normalises with `CigarBuilder(false)`. The new alignment start is
   `newStart + alignStartHapwrtRef`, the location is the padded span, and score and kmer are copied.
9. All haplotypes are trimmed with `ignoreRefState = true`, so they are keyed by bases only. A duplicate
   keeps the first haplotype, unless the duplicate is the reference, which replaces it. The kept haplotype
   is then rebuilt as a reference haplotype: `new Haplotype(bases, true)` copies the cigar, location, score
   and alignment start but **not the kmer size, which becomes 0**.
10. The trimmed list is sorted by `SIZE_AND_BASE_ORDER` (length, then bases). The reference is not
    necessarily first. This is the likelihood allele order.
11. `isVariationPresent()` is the `variationPresent` flag of the *untrimmed* set (any non-reference
    haplotype) **and** more than one trimmed haplotype. When it is false the exit is "no_assembly_variation".

**Genotyping reads**
12. Each read is `ReadClipper.hardClipToRegion(read, paddedStart, paddedEnd)` (code shared with `hc_preproc`,
    which was fuzz-validated). Empty reads and reads not overlapping the padded span are dropped. The rest
    are stably sorted with `ReadCoordinateComparator`.
13. Stub reads are removed **silently**: unclipped length (bases minus soft clips) < 10, so they are not
    "filtered reads". Then `filterNonPassingReads` collects the remaining reads with length < 10, MAPQ < 20,
    or a mapped mate on another contig, in read order. These are the "filtered reads".
14. If no reads are left, the exit is "no_reads_for_genotyping". It is tested after the variation check.

**Read likelihoods**
15. `modifyReadQualities` uses GOPs of 45 (no `BI`/`BD`) and the PCR model on every position except the
    last: `min(gop, table[repeatLength(i)])`, with `table[k] = max(10, fastRound(40 - exp(k / (3π)) + 1))` =
    40,40,40,40,39,…,33. `findTandemRepeatUnits` is reproduced exactly, including the backward repeat count
    left over from the last unit length tried, the forward/backward unit comparison and the
    `MAX_REPEAT_LENGTH = 20` cap. The fast table-driven version is checked against a direct transcription
    on 200,000 random sequences.
16. `capMinimumReadQualities`: base qual = min(qual, MAPQ). Then signed-byte comparisons set qual < 18 to 6
    and ins/del < 6 to 6. GCP is 10 everywhere.
17. The PairHMM batch is genotyping reads × trimmed haplotypes, in that order. The matrix is transposed
    to `[haplotype][read]`.
18. `normalizeLikelihoods(-4.5, symmetric = true)`: per read, best = max over *all* alleles (strict `>`
    scan from allele 0). Every value below best − 4.5 is raised to it. Nothing happens with a single
    haplotype.
19. `filterPoorlyModeledEvidence`: a read is removed when its max likelihood (after normalisation) is
    below `min(2, ceil(len × 0.02)) × -4`, where `len` is the HMM base-quality array length. The removed
    reads keep their order as the filtered evidence.

**PairHMM numerics**
20. The CUDA kernel (`hc_cudahmm`) follows GKL's algorithm in float (2^120 scale), with a double fallback
    below 1e-28. It is bit-exact: separate `sumM`/`sumX`, GKL's FMA placement, `__fdiv_rn` (nvcc otherwise
    turns `distm / 3.0f` into a multiply by a rounded 1/3), and explicit flush-to-zero of double results (GKL's
    MXCSR.FTZ covers doubles; nvcc `-ftz` does not). `src/bin/hmmcmp.rs` compares it with the CPU port pair by
    pair (`WILD=1` for arbitrary qualities, `CELLS=1` to locate the first differing DP cell).
21. The CPU fallback (`cpu.rs`) is an operation-for-operation port of GKL's `compute_full_prob`. It builds
    the same tables from the same libm calls (`pow`, `powf`, `log1p`, `log10`, the truncated
    `INV_LN10 = 0.434294`, the Jacobian table and `approximateLog10SumLog10` in `NUMBER` precision). It
    sets flush-to-zero as `initNative` does, uses 'N' as a wildcard, maps other non-ACGTN bases to 'A', and
    masks quals with `& 127`. It sums the last row's M and X separately and returns
    `(double)(log10f(sum) − log10f(2^120))`. GKL's AVX-512 build is FMA-contracted: from its disassembly,
    `M = fma(Y2, gapm, fma(X2, gapm, M2*pMM)) * prior`, `X = fma(X1, pXX, M1*pMX)` and
    `Y = fma(Y1, pYY, M1*pMY)`. `Flavor::Avx512` reproduces that; `Flavor::Avx2` uses plain mul/add, like GKL's
    AVX2 build. `compute_batched` runs 16 pairs in SIMD lanes, sorted by size, and gives the same bits.

## Not ported / limitations

* Non-default modes: pileup detection and PD haplotypes (PDHMM), force-calling alleles, DRAGEN/dragstr
  PairHMM scores, flow mode, stepwise filtering, allele filtering, legacy trimming, `--disable-optimizations`,
  dynamic read disqualification (implemented but not validated), more than one sample, and BQSR `BI`/`BD`
  indel qualities.
* GATK removes stubs and filtered reads with `List.removeAll`, which uses `SAMRecord.equals`, so it would also
  remove a different read that is identical in every field. The port removes by identity. This can only
  differ for two fully identical records; duplicates are filtered upstream.
* The read comparator uses `ReadUtils.getSAMFlagsForRead` (rebuilt flags) and `getMateReferenceIndex`
  (`mate_ref_id`, -1 for an unmapped mate), as GATK does.

## Files

* `src/lib.rs`: API, `prepare`/`finish`, normalisation and filtering.
* `src/events.rs`: `Event`, `EventMap`, `getVariationEvents`.
* `src/trim.rs`: tandem repeats, `AssemblyRegionTrimmer.trim`, `Haplotype.trim`, `trimTo`.
* `src/clip.rs`: `hardClipToRegion` and `ReadCoordinateComparator` (from `hc_preproc`).
* `src/hmm.rs`: PCR indel error model and quality caps.
* `src/cpu.rs`: GKL-exact CPU PairHMM (scalar and 16-lane).
* `src/bin/validate.rs`: the validator (`--batch N`, `--threads N`, `--limit N`, `--cpu avx512|avx2`,
  `--mates file`). It prints per-check counts, code-path coverage and stage timings. Set `SHOWDIFF=1` to
  list every non-identical likelihood.
* `build.rs`: adds rpaths for `libcudapairhmm.so` and for `libstdc++` in `$CONDA_PREFIX/lib`, because cargo does not propagate
  `hc-cudahmm`'s link arguments to this crate's binaries.
* `scratch/Fuzz.java` (development probe, not in the repository): the differential fuzzer. The generated dumps (about 130 MB per 4,000 regions) were deleted; they were made with the commands above (seeds 2, 3 and 4, 4,000 regions each).

Build (from the repository root, in an environment with htslib; `environment.yml` at the repo root):
`RUSTFLAGS="-C target-cpu=native" cargo build --release -p hc-likelihood`.
Tests: `cargo test --release -p hc-likelihood`.
