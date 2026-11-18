# hc-refconf: GATK 4.7 HaplotypeCaller `-ERC GVCF` in Rust

This crate ports the reference-confidence part of HaplotypeCaller's GVCF mode (`--emit-ref-confidence GVCF`,
default arguments otherwise):

* `ReferenceConfidenceModel.calculateRefConfidence`: the pileups, the ref-vs-any likelihoods, the indel
  informativeness and the per-position GQ/PL (`model.rs`);
* `GVCFWriter` / `GVCFBlockCombiner` / `HomRefBlock` with the default GQ bands (`gvcf.rs`), plus the writer of
  the unbanded per-position records that GATK's `GATK_RAW_CALLS_VCF` debug option produces.

The ERC branches of genotyping and annotation live in [`../genotype`](../genotype) (`src/erc.rs`). The pipeline wiring
(inactive regions, early exits, flanks, the writer thread) lives in `variant-rs` (`genotype_erc_batch` in [`../../src/backend.rs`](../../src/backend.rs)).

## Results (NA12878 30x, chr20:10,000,000–15,000,000)

```
$ target/release/variant-rs -I $VARCALL_DATA/data/NA12878.chr20_10-15M.bam -R $VARCALL_DATA/data/chr20.fa -L chr20:10000000-15000000 \
      -ERC GVCF -O rust.g.vcf.gz --gvcf-raw rust_raw.vcf --pairhmm cuda --threads 16
GVCF: regions 24567 (active 6345, genotyped 6270), calls 8905, GVCF lines 798047, raw lines 5000001, assembly failures 0

$ target/release/gvcfcmp $VARCALL_DATA/gvcf/gatk_raw_calls.vcf.gz rust_raw.vcf
lines compared (same POS): 5000001   identical: 5000001   differing: 0
  of which with a real ALT on either side: 8905   identical: 8905
only in GATK: 0   only in rust: 0
RESULT: IDENTICAL

$ target/release/gvcfcmp $VARCALL_DATA/gvcf/gatk.g.vcf.gz rust.g.vcf
lines compared (same POS): 798047   identical: 798047   differing: 0
  of which with a real ALT on either side: 8905   identical: 8905
only in GATK: 0   only in rust: 0
RESULT: IDENTICAL
```

The whole files are byte-identical too, header included except the command-line line (variant-rs now generates the
header: `##variant-rsCommandLine` instead of `##GATKCommandLine`, same position; the development script `validate.sh`, not in the repository, checked this). This holds with `--threads 1` and `--threads 16` and with the CPU PairHMM. VCF mode is
unchanged: its 8,064 records are still byte-identical, and every crate's validator still passes.

### Differential tests on synthetic data (GATK's own classes from gatk.jar)

The Java fuzzers and `scratch/fuzz.sh` below are development probes, not in the repository.

| Test | Inputs | Result |
|---|---|---|
| `scratch/fuzz/RefconfFuzz.java` + `refconf_fuzzcheck`: GATK `ReferenceConfidenceModel.calculateRefConfidence` vs `calculate_ref_confidence`. The inputs are random references with homopolymers, STRs and Ns; reads with mismatches, N/IUPAC bases, insertions, deletions, soft and hard clips and low or very low qualities; short fragments (adaptor exclusion); up to 350 reads per region; overlapping calls and calls inside deletions | 6,300 regions, 512,588 positions | all identical |
| `scratch/fuzz/GvcfFuzz.java` + `gvcf_fuzzcheck`: GATK `GVCFWriter` (default bands) plus HaplotypeCaller's raw writer, vs `GvcfWriter`. The streams contain sites across all GQ bands, deletions, het/hom/hom-ref/no-call calls, biallelic `<NON_REF>` calls (hom-ref and no-call, which start or extend blocks) and missing DP | 620,000 records | banded and raw identical |
| genotype crate: `ERC=1 scratch/fuzz.sh SEED K` (seed 101 with K = 5, seeds 102–105 with K = 3) runs GATK `assignGenotypeLikelihoods` in ERC mode on mutated regions: more than 6 alt alleles, more than 44 alleles, compound events, filtered reads, ties, and so on | 9,549 regions, 14,470 calls | all identical |
| `genotype-validate --erc` on `$VARCALL_DATA/gvcf/gatk_genotyping_erc.jsonl` | 6,270 regions, 8,905 calls | all identical |

Code paths that the real data never reaches but the fuzzers do: adaptor exclusion in the reference model, the
"reference ends first" branch of the indel informativeness, reduction to 6 alt alleles in ERC mode,
`removeAltAllelesIfTooManyGenotypes`, and blocks started by a call (they keep the call's QUAL).

One command ran the GVCF mode and both comparisons: `validate.sh [OUTDIR] [THREADS]` (a development script, not in the
repository: the `variant-rs` and `gvcfcmp` commands above, the whole-file header check and a tabix check). It exited
with 0 only if everything was identical.

## Usage

```
variant-rs ... --gvcf out.g.vcf [--gvcf-raw raw.vcf] [--vcf-header header.txt]
```

* `--gvcf`: banded GVCF (plain text). Without `--vcf-header`, the header is GATK's header for this run
  (`tables/gvcf_header.txt`, copied verbatim from `$VARCALL_DATA/gvcf/gatk.g.vcf.gz`).
* `--gvcf-raw`: the records `callRegion` returns, one per position and before banding, exactly as GATK writes
  them with `GATK_RAW_CALLS_VCF` (header `tables/raw_header.txt`).
* `gvcfcmp <gatk.vcf[.gz]> <rust.vcf[.gz]> [--show N] [--only-variants]` joins two files on POS and prints
  identical and differing counts, mismatches by column and by INFO/FORMAT key, and examples. It reads raw text,
  gzip included.

Library API:

```rust
// one region: reference haplotype = reference bases of the padded span
let entries: Vec<Entry> = hc_refconf::calculate_ref_confidence(ref_padded_span, padded_start, span_start, span_end, &reads, &calls);
let sites: Vec<SiteRecord> = hc_refconf::ref_model_for_no_variation(&finalized_reads, span_start, span_end, padded_start, ref_padded_span);
let mut w = hc_refconf::GvcfWriter::new("chr20", out, Some(raw_out), &hc_refconf::default_gq_bands());
w.add_region(items)?;   // the records of one region, in order (Item::Site / Item::Call)
w.finish()?;
```

Build (from the repository root, in an environment with htslib; `environment.yml` at the repo root): `RUSTFLAGS="-C target-cpu=native" cargo build --release -p hc-refconf`

## What GVCF mode changes, and where it is ported

| GATK 4.7 | Rust |
|---|---|
| `HaplotypeCallerEngine.validateAndInitializeArgs`: calling confidence -0.0, `annotateAllSitesWithPLs`, physical phasing on | `EngineConfig::new_erc`, `GenotypingEngine::genotype_region_erc` (genotype crate) |
| `HaplotypeCaller.makeVariantAnnotations` → `filterReferenceConfidenceAnnotations`, `VariantAnnotatorEngine(useRaw = true)` | `erc.rs::annotate_call_erc` |
| `HaplotypeCallerGenotypingEngine.assignGenotypeLikelihoods(emitReferenceConfidence = true)`, `addNonRefSymbolicAllele`, `AlleleLikelihoods.addNonReferenceAllele` / `updateNonRefAlleleLikelihoods`, `prepareReadAlleleLikelihoodsForAnnotation` | `erc.rs` |
| `GenotypingEngine.calculateGenotypes` with `<NON_REF>` (`calculateMostLikelyAlleles`, `forceKeepAllele`, `noAllelesOrFirstAlleleIsNotNonRef`, `recordDeletions`) | `erc.rs::calculate_genotypes_erc` |
| `AssemblyBasedCallerUtils.phaseCalls` (constructHaplotypeMapping, constructPhaseSetMapping, constructPhaseGroups, phaseVC) | `erc.rs::phase_calls` |
| `HaplotypeCallerEngine.callRegion` ERC exits, `referenceModelForNoVariation`, `containsCalls`, the flanks from `AssemblyRegionTrimmer.Result` | `variant-rs/src/backend.rs::genotype_erc_batch` |
| `ReferenceConfidenceModel.calculateRefConfidence`, `makeReferenceConfidenceVariantContext`, `calcGenotypeLikelihoodsOfRefVsAny`, `doIndelRefConfCalc`, `calcNReadsWithNoPlausibleIndelsReads`, `readHasNoPlausibleIdealsOfSize`, `getIndelPLs` | `model.rs` |
| `AssemblyBasedCallerUtils.getPileupsOverReference` (`LocusIteratorByState`, `AlignmentStateMachine`, `PileupElement`) | `model.rs` (streamed pileup) |
| `GVCFWriter`, `GVCFBlockCombiner`, `HomRefBlock`, `GVCFBlock` | `gvcf.rs` |

## GATK behaviours replicated

### Engine settings and genotyping (`variant-rs/crates/genotype/src/erc.rs`)

1. **Calling confidence -0.0.** `passesThreshold(allele, -0.0)` compares `log10P(absent) + 1e-10` with
   `qualToErrorProbLog10(-0.0) = 0.0`. Every non-monomorphic site is emitted, whatever its QUAL, and FILTER stays `.`.
2. **`forceKeepAllele` is always true in reference-confidence mode** (`... || referenceConfidenceMode != NONE`).
   Every allele is output, including implausible alts and `<NON_REF>`, except spurious spanning deletions.
   This is why hom-ref calls such as `A G,<NON_REF> ... 0/0` appear (710 here). A monomorphic site is still
   dropped, unless its first output alt is `<NON_REF>`.
3. **QUAL.** `CommonInfo.getPhredScaledQual()` is `log10PError * -10 + 0.0`, so a site with P(no variant) = 1
   prints `0`, not `-0` (13 spanning-deletion sites here).
4. **`<NON_REF>` likelihood.** `addNonReferenceAllele` adds it with -∞ and then sets, per read, the
   commons-math `Median` (LEGACY percentile: `lower + (pos - floor(pos)) * (upper - lower)`, pos = (n+1)/2) of
   the likelihoods of the concrete alleles that are strictly below the read's best. If there are none, the median
   is NaN and the read gets the best likelihood. Before annotation, when the call dropped alleles, it is recomputed
   from the call's alleles only, but the best is still searched over *all* alleles, the current `<NON_REF>`
   included. Filtered reads are added afterwards with likelihood 0 for every allele, `<NON_REF>` included.
5. **`<NON_REF>` has `Allele.length() == 0`.**
   * The AF calculator gives it the indel pseudocount.
   * `recordDeletions` sees `ref.length() - 0 > 0` and records a "deletion" `[start, start + refLength]` for
     *every* call. This makes a `*` allele at the next position(s) look covered by an upstream deletion.
6. **Max-alt-allele reduction.** `calculateMostLikelyAlleles` does not count `<NON_REF>` as a proper alt and
   always keeps it. The reduction branch (`subsetAlleles`, which re-derives the PLs) runs whenever the alt count
   *including* `<NON_REF>` exceeds 6, even when nothing is removed.
7. **Annotations.**
   * The set is the HC defaults minus ChromosomeCounts, FisherStrand, StrandOddsRatio and QualByDepth, plus
     StrandBiasBySample.
   * `RAW_MQandDP=sumMQ²,n` covers all annotation evidence (filtered reads included) with MQ ≠ 255 and is
     written whenever there is evidence. There is no `finalizeRawMQ`.
   * `SB` uses `getContingencyTable(min count 0)`. AD and SB treat `<NON_REF>` as an alt allele; so do the rank
     sums (`vc.hasAllele`).
   * AD, DP and SB are written only for called genotypes, so a no-call is `GT:GQ:PL`. ExcessHet is likewise only
     written for called genotypes.
   * `MLEAF = Math.min(1.0, AC/AN)` prints `NaN` for a no-call (Rust's `f64::min` would give 1.0).
8. **Reverse trimming** skips symbolic alleles (`<NON_REF>`), as it skips `*`.
9. **Physical phasing.**
   * Only calls with exactly one site-specific alt (not ref, not `<NON_REF>`, not `*`) take part.
   * The haplotypes considered are those of `calledHaplotypes` (the union over calls of the haplotypes of every
     output allele, ref included) whose event map has an event at the *call's* start with the call's alt bases.
   * Groups follow the same "same haplotypes / on all alt haplotypes / complementary" rules, in call order. A
     conflict clears all phasing for the region (2 regions here).
   * PID = `start_REF_ALT0` of the group's first call, PS = its start, and PGT = `0|1` / `1|0`.
   * A het genotype is reversed when the allele at the phased index is not site-specific, which gives `1|0` or `2|1`.
   * Non-het genotypes are only marked phased, so they print `0|0`, `1|1` or `.|.`.
   * FORMAT order: `GT:AD:DP:GQ:PGT:PID:PL:PS:SB`.

### Which records a region produces (`variant-rs`, `run_erc_batch`)

10. `callRegion` returns the records of `referenceModelForNoVariation(region)` for:
    * an inactive region, which is first finalized (`finalizeRegion`, the `hc_preproc` port);
    * an active region with no reads;
    * the exits after assembly (no variation after trimming, no assembly variation, no reads for genotyping);
    * a genotyped region whose calls contain no called genotype (`containsCalls`).

    In all of these, the reads are the region's finalized reads minus `filterNonPassingReads` (unclipped length < 10
    or MAPQ < 20), and the reference haplotype is the padded span.
11. When the region does contain called genotypes, the output is:
    * the left flank, if the region starts before the genotyping span;
    * `calculateRefConfidence` over `regionForGenotyping`, with the trimmed reference haplotype, the realigned
      evidence reads (poorly modelled reads excluded) and the calls;
    * the right flank.

    A flank is `originalRegion.trim(flank, 100)`: its padded span is the flank ± 100, clipped to the contig and to
    the region's padded span. Its reads are the finalized reads after `hardClipToRegion`, with empty and
    non-overlapping reads dropped, sorted, then `filterNonPassingReads`.

### Reference confidence model (`model.rs`)

12. **Pileups** are `LocusIteratorByState` over the reads sorted with `ReadCoordinateComparator`. Reads join the
    active list at their start and leave after their end; elements come in that order. Deletions are included.
    Bases inside the adaptor are dropped by `isBaseInsideAdaptor`: its boundary is computed on the current
    (realigned or clipped) read, and its `fragmentLength > 100` test uses the *signed* TLEN.
13. **Ref-vs-any likelihoods.**
    * `readsWereRealigned` is always true here, so alt = `base != ref || deletion`, with no neighbour or
      soft-clip rules.
    * Deletions have qual 30 and are never skipped; other bases with qual ≤ 6 are skipped.
    * The original-soft-clip check (`os`/`oe` tags) is skipped, because the model is built with
      `useSoftClippedBases = !overrideSoftclipFragmentCheck = true`.
    * Het = Σ `approximateLog10SumLog10` (GATK's Jacobian table) - n·log10(2). The het and hom-var values are
      capped by hom-ref, then turned into PLs with `GLsToPLs` (Java `Math.round`).
14. **Indel likelihoods.**
    * The cached diploid indel GLs for n informative reads (capped at 40) are dumped from GATK
      (`tables/RefconfTables.java`).
    * The final PLs are the indel PLs only if their hom-ref GQ (`getGQLog10FromLikelihoods`) is strictly greater
      than the SNP one; ties go to the SNP PLs.
    * GQ = `calculateGQFromPLs`. AD = [ref, non-ref] counted over the bases used; DP = their sum.
15. **Indel informativeness.**
    * The per-read `BitSet` (`IDL` transient attribute) is built the *first* time a read is queried, using that
      position's read and reference offsets. The global mismatch cost (`baselineMMSums[0]`), the "last relevant"
      index and the reference-ends-first correction all depend on that first position, and every later position
      reads the cached bits.
    * The count returns early at 41 informative reads, so the reads after that point do not get their cache
      built at that position. This happened 169,412 times here.
    * The cache is cleared at the end of every `calculateRefConfidence` call.
16. **Read offsets for indel informativeness.**
    * `getCigarModifiedOffset` counts soft clips as reference-consuming, but
      `getBasesAndBaseQualitiesAlignedOneToOne` does not copy soft-clipped bases when the read has indels. The
      array keeps zeros there; this is replicated.
    * Mismatches use `Nucleotide.intersect` (IUPAC masks; `N` matches everything, unknown bytes nothing) and
      exclude the gap `-`.
    * Base qualities are added as signed bytes.
17. **Positions covered by calls.** At a position where a call starts, the call is emitted. A position covered by
    a call that starts earlier (inside a deletion) still gets a reference-confidence record (3,602 here). The
    GVCF writer drops it later. `getOverlappingVariantContext` keeps the overlapping call with the largest start.

### GVCF writer (`gvcf.rs`)

18. **Bands** are `[0,1) … [59,60) [60,70) [70,80) [80,90) [90,99) [99,100)`.
    * A site joins the current block if `min(GQ, 99)` (the *site's* GQ) is in the band of the block's first site
      and the PL arrays have the same length.
    * The block's GQ is `calculateGQFromPLs` of the element-wise minimum PLs, so it can fall outside the band.
    * DP is `Math.round` of the commons-math median of the DPs, and MIN_DP is the minimum.
    * The FORMAT is `GT:DP:GQ:MIN_DP:PL`, or `GT:GQ:PL` when no record had DP. INFO is `END=…`.
19. **Hom-ref sites under a variant.**
    * After a variant record, `nextAvailableStart = vc.getEnd()`. A following hom-ref site with
      `start ≤ nextAvailableStart && end ≤ nextAvailableStart` is dropped (3,463 here), and the first site that
      does not satisfy this resets it.
    * A non-contiguous record flushes the block.
    * A variant record flushes the block and is written unchanged.
20. **Calls treated as hom-ref sites.** The combiner also treats a call as a hom-ref site when it is
    `[REF, <NON_REF>]` with GT 0/0, or a no-call with PL[0] = 0. A block started by such a call keeps the call's
    QUAL (`new VariantContextBuilder(startingVC)` only clears INFO).
21. **PL aliasing in the unbanded records.**
    * htsjdk's `FastGenotype.getPL()` returns its internal array, and `HomRefBlock` keeps it as `minPLs`,
      updating it in place as sites are added.
    * HaplotypeCaller writes the raw records only after the whole region went through the GVCF writer. The raw
      line of a record that *starts* a block therefore shows the block's running minimum PL (up to the end of that
      region), while its GQ is the original. For example, `chr20 10000000 … 0/0:25,0:25:75:0,72,834`, where the
      record's own PL was `0,75,…`.
    * The port keeps the block's starting record in the region buffer and applies every min update to it
      (3,868,597 updates here). The banded GVCF is not affected.

## Timings

Same machine (Ryzen 9 9950X3D, 16 cores; RTX 5070 Ti), 3 runs each, wall clock including start-up. GATK writes
`.g.vcf.gz` plus index (GKL AVX-512 PairHMM, 4 OpenMP threads). Rust writes a plain `.g.vcf` with the CUDA PairHMM.
The script is `scratch/timing.py` (a development script, not in the repository).

| Run (monotonic wall clock, start to exit) | Runs | Median | vs GATK |
|---|---|---:|---:|
| GATK HaplotypeCaller `-ERC GVCF` | 57.34 / 53.99 / 54.21 s (GATK's own log: 0.83 / 0.86 / 0.82 min) | 54.2 s | 1.0x |
| **Rust GVCF mode, CUDA PairHMM, `--threads 1`** | 13.01 / 13.27 / 12.88 s | 13.0 s | **4.2x** |
| **Rust GVCF mode, CUDA PairHMM, `--threads 16`** | 4.30 / 4.30 / 4.28 s | 4.3 s | **12.6x** |

All six Rust outputs were byte-identical to `$VARCALL_DATA/gvcf/gatk.g.vcf.gz`. In the earlier timing set, GATK's three outputs
had the same data lines as that file. That earlier set of runs timed with `/usr/bin/time` (GATK 51.1 / 48.8 / 49.2 s) was discarded for timing purposes: under
WSL2 the realtime clock jumped during the Rust runs (one reported 2.03 s while the program's own monotonic clock
measured 4.23 s), so `scratch/timing.py` uses Python's monotonic clock.

Stage breakdown for Rust with 16 threads (GVCF mode adds roughly 1 s over VCF mode):
* streaming (BAM decode, pileup scan, regions, read copies): about 2 s, single thread;
* finalizeRegion of all 24,567 regions: 0.34 s, serial (serial is faster than the rayon version, which suffered
  cross-thread free contention);
* assembly: 0.55 s;
* PairHMM: 0.15 s;
* realignment: 0.12 s;
* ERC genotyping: 0.17 s;
* reference confidence: 0.24 s (3.5 s of single-thread work: 158 M pileup elements and 2.07 M indel caches);
* banding and writing: 0.45 s on a separate writer thread, overlapped (inline with `--threads 1`).

## Files

* `src/model.rs`: reference confidence model.
* `src/gvcf.rs`: GVCF writer and raw writer.
* `src/stats.rs`: path counters (printed by `variant-rs`).
* `src/lib.rs`: API and `ref_model_for_no_variation`.
* `tables/RefconfTables.java` produces `refconf_tables.bin`: `qualToProbLog10`, `LOG10_ONE_THIRD`, `log10(2)` and
  the indel GL cache, all from gatk.jar.
* `tables/gvcf_header.txt` and `tables/raw_header.txt`: GATK's headers.
* `src/bin/gvcfcmp.rs`: the comparison tool.
* `src/bin/refconf_fuzzcheck.rs` with `scratch/fuzz/RefconfFuzz.java`, and `src/bin/gvcf_fuzzcheck.rs` with
  `scratch/fuzz/GvcfFuzz.java` (the Java halves are development probes, not in the repository): the differential fuzzers. Their inputs were made with
  `java -cp gatk.jar scratch/fuzz/RefconfFuzz.java SEED N f.txt` and
  `java -cp gatk.jar scratch/fuzz/GvcfFuzz.java SEED N prefix`; the generated files were deleted.
* `validate.sh` (development script, not in the repository): the end-to-end validation (raw records and banded GVCF against GATK).
* `scratch/timing.py` (not in the repository): the timing script.
* `scratch/erc/gatk_erc_calls.vcf` (not in the repository): the variant lines of `$VARCALL_DATA/gvcf/gatk_raw_calls.vcf.gz`, used by the genotype
  crate's `genotype-validate --erc`.

## Limitations

* Same scope as the rest of the port: single sample, diploid, one contig, default arguments.
* Not ported: `BP_RESOLUTION`, `--floor-blocks`, custom `--gvcf-gq-bands` (the writer takes bands, but there is
  no command-line option), `--population-callset` / priors, DRAGEN and flow modes, and the flow-based deletion
  quality.
* `filterNonPassingReads` (`passes_hc_read_filter`) includes `MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE`, from the reads'
  `ref_id` / `mate_ref_id`; the pileup comparator uses GATK's rebuilt SAM flags and mate reference index.
* GATK aborts when assembly fails. The port falls back to the reference model for that region (it never happens
  here).
* Output: variant-rs writes bgzip + tabix (`.g.vcf.gz`) or plain text, with a generated header (see [docs/notes/variant-rs.md](../../../docs/notes/variant-rs.md)).
* With the CUDA PairHMM, GVCF output can differ from GATK in rare deep or `<NON_REF>`-sensitive sites (the kernel's
  likelihoods are within 3.8e-6 of GKL's, which can flip a `<NON_REF>` median or a poorly-modelled read); `--pairhmm cpu`
  (bit-exact GKL port) removes it. Seen in 1 of 182,673 lines on alt/decoy contigs and 4 of 868 lines on chrM.
