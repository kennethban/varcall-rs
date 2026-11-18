# hc-genotype — GATK 4.7 HaplotypeCaller genotyping and VCF output in Rust

This crate ports the last stage of `HaplotypeCallerEngine.callRegion`, from the point after
`computeReadLikelihoods(...)` to the returned calls, using **default arguments** and normal VCF
output (no `-ERC`). It covers:

1. `AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype` and `changeEvidence`
2. `HaplotypeCallerGenotypingEngine.assignGenotypeLikelihoods`: event maps, merged VCs, the allele
   mapper, marginalisation, genotype likelihoods, the AF calculator and QUAL, genotype calls, all
   default annotations and reverse trimming
3. `RMSMappingQuality.finalizeRawMQ`
4. htsjdk `VCFEncoder` formatting of the data lines

## Status

`target/release/genotype-validate` runs on GATK's dump of NA12878 30x, chr20:10–15 Mb
(`$VARCALL_DATA/geno/gatk_genotyping.jsonl`, 6,270 genotyped regions) and compares with `$VARCALL_DATA/geno/gatk.vcf.gz`:

```
== realignment (6270 genotyped regions)
reads: 271620  (changed by realignment: 4278)  start match: 271620  cigar match: 271620  both: 271620
== calls per region
regions with rust count == n_calls == GATK VCF count: 6270   differing: 0
records: rust 8064  GATK VCF 8064
== VCF lines (paired by region and POS)
compared: 8064  identical: 8064  unpaired GATK: 0  unpaired rust: 0
QD jitter draws: 540
== time: realignment 0.92s, genotyping+annotation+formatting 0.11s, total wall 2.0s
```

All 8,064 data lines are byte-identical across all 10 columns. Time is single-threaded with
`-C target-cpu=native`; JSON loading is excluded. GATK's own code for the same stage takes 3.9 s
in-JVM on the same inputs (measured with `RunGeno.java`, a development probe not in the repository; one JVM, JIT warm-up
included). In realignment, 258,446 reads use the exact-substring shortcut and 13,174 need a full
Smith-Waterman alignment (counters from `stats::report()`).

### Differential tests on synthetic regions

`scratch/fuzz.sh SEED K` takes every K-th genotyped region of the dump and applies random mutations
with `scratch/mutate.py` (both development scripts, like the Java probes below; none of them is in the repository):

* filtered reads, some with MQ < 20;
* random MAPQ values (0, 255, …);
* strand flips;
* likelihood noise;
* exact ties and near ties at the 0.2 informative threshold;
* reads duplicated up to 10× (tables larger than 400 reads);
* 6–14 extra haplotypes carrying SNPs, insertions, deletions and SNP+indel compound events at one locus;
* 46–60 haplotypes at one locus.

It then runs **GATK's own classes** on the result (`RunGeno.java`: GATK realignment,
`HaplotypeCallerGenotypingEngine`, the default annotations, `finalizeRawMQ` and htsjdk
`VCFEncoder`). Finally `genotype-validate` compares the Rust output with GATK's. `RunGeno` on the unmodified
dump reproduces `gatk.vcf.gz` byte for byte. Seeds 1–5 gave 14,545 synthetic records, all
identical, and every realigned read matched. Paths this exercised that the real data never reaches:

* reduction to 6 alternate alleles (about 5,000 loci);
* `removeAltAllelesIfTooManyGenotypes` (more than 44 alleles, about 1,600 loci);
* about 8,700 compound events;
* filtered reads in the annotations;
* FS table normalisation;
* about 4,600 QD jitter draws.

`target/release/numcheck` checks the numeric kernels against Java test vectors
(`NumProbe.java`, a development probe not in the repository):

| Check | Values | Mismatches |
|---|---:|---:|
| `String.format` `%.2f`/`%.3f`/`%.4f` | 2,000,000 | 0 |
| FisherStrand and SOR strings (random 2×2 tables) | 300,000 | 0 |
| Mann-Whitney z (exact and normal paths), bitwise | 100,000 | 0 |

## API

```rust
use hc_genotype::{EngineConfig, GenotypingEngine, RegionInput, VcfRecord};

let mut engine = GenotypingEngine::new(EngineConfig::new("chr20", 64444167)); // contig, contig length
for region in regions_in_genomic_order {
    let input = RegionInput {
        haplotypes: &haps,          // readLikelihoods.alleles() (hc_assembly::Haplotype, includes the ref haplotype)
        reads: &reads,              // readLikelihoods.sampleEvidence(0), BEFORE realignment (hc_assembly::Read)
        matrix: &matrix,            // Vec<Vec<f64>> [haplotype][read], log10
        filtered_reads: &filtered,  // filterNonPassingReads(regionForGenotyping)
        filtered_evidence: &lkfilt, // readLikelihoods.filteredSampleEvidence(0) (has no effect, see below)
        ref_padded: &ref_bytes,     // assemblyResult.getFullReferenceWithPadding()
        ref_padded_start,           // assemblyResult.getPaddedReferenceLoc().getStart()
        span_start, span_end,       // regionForGenotyping.getSpan()
    };
    let records: Vec<VcfRecord> = engine.call_region(&input);   // realign + genotype
    for r in &records { println!("{}", r.to_line()); }          // exact VCF data line
}
```

* `GenotypingEngine::call_region` realigns the reads, then genotypes the region. The two steps can
  be split:
  * `realign::realign_reads_to_best_haplotype(&mut SwAligner, &mut reads, haps, matrix, ref_padded_start)`
    is pure. It can run in parallel across regions, with one `SwAligner` per thread.
  * `engine.genotype_region(&input, &realigned)` must then be called **sequentially in genomic
    order**. The engine carries two pieces of cross-region state, exactly as GATK does: the
    upstream-deletion list that decides whether a `*` allele is emitted, and the global random
    generator used for QD jitter. This step costs about 0.1 s per 5 Mb.
* `VcfRecord` has `pos`, `alleles` (reference first), `qual`, the pre-formatted `info` pairs (sorted),
  `gt`, `ad`, `dp`, `gq` (uncapped) and `pl`. `to_line()` produces the 10-column line. The header
  can be copied verbatim from GATK's VCF.
* `realign::create_read_aligned_to_ref` aligns a single read.
* `stats::report()` returns code-path hit counters.
* `EngineConfig::defer_qd_jitter` (default false): the engine writes `QD_JITTER_PLACEHOLDER` (byte 1) instead of
  drawing the jitter; the caller replaces the placeholders, in output order, with `30 + 3·nextGaussian()` (`%.2f`)
  from its run-wide generator (`resolve_qd_jitter`, or variant-rs's `rng::resolve_qd` with `hc_downsample::JavaRandom`).
  variant-rs uses this so that contigs can be genotyped in parallel while the draws keep GATK's order, interleaved with
  the downsampler's draws.
* `GenotypingEngine::{rng_state, set_rng_state}` (added by a development patch, see [`../downsample/README.md`](../downsample/README.md))
  expose the jitter generator's state.

Binaries:

* `genotype-validate [--geno F] [--asm F|-] [--vcf F] [--show N] [--out F]` validates against the dump (default: under `$VARCALL_DATA`). With
  `--asm -`, the reference must be inlined in the regions (as `fuzz.sh` does).
* `numcheck vectors.tsv` checks the numeric kernels against Java's.

Build (from the repository root, in an environment with htslib; `environment.yml` at the repo root): `RUSTFLAGS="-C target-cpu=native" cargo build --release -p hc-genotype`

## Defaults (checked in the 4.7.0.0 sources and the run's command line)

| Setting | Default | Effect here |
|---|---|---|
| `--flow-filter-alleles` (`filterAlleles`), `stepwiseFiltering`, `--use-pdhmm` / `--pileup-detection` | false | Paths not taken |
| `--emit-ref-confidence` | NONE | Calls go through `finalizeRawMQ`; no GVCF |
| Physical phasing | off | `HaplotypeCallerEngine` turns `doNotRunPhysicalPhasing` on when ERC is NONE (logged "Disabling physical phasing…"), so there are no PGT/PID/PS fields |
| `--apply-bqd` / `--apply-frd` | false | `IndependentSampleGenotypesModel`; the read-overlap predicate is `target.overlaps(read)` |
| DRAGstr params | none | Priors from `GenotypePriorCalculator.assumingHW`, unused by `USE_PLS_TO_ASSIGN` |
| `--max-mnp-distance` | 0 | |
| `--allele-informative-reads-overlap-margin` | 2 | |
| `--disable-spanning-event-genotyping` | false | `*` alleles are genotyped |
| `--standard-min-confidence-threshold-for-calling` | 30 | |
| `--max-alternate-alleles` | 6 | |
| `--max-genotype-count` | 1024 | Gives 44 alleles at ploidy 2 |
| Heterozygosity, indel heterozygosity, stdev | 0.001, 1.25e-4, 0.01 | |
| `--genotype-assignment-method` | USE_PLS_TO_ASSIGN | |
| `--use-posteriors-to-calculate-qual` | false | |
| `--output-mode` | EMIT_VARIANTS_ONLY | |
| Sample contamination | none | |
| `--use-filtered-reads-for-annotations` | false | |
| Read-to-haplotype SW | match 10, mismatch -15, gap open -30, gap extend -5 | SOFTCLIP overhang |
| `--smith-waterman` | FASTEST_AVAILABLE (GKL) | Exact-substring shortcut first |

**Annotations.** The groups are `StandardAnnotation` + `StandardHCAnnotation`, and the annotation
engine runs with `useRaw=false`.

* INFO: AC, AF, AN (ChromosomeCounts), BaseQRankSum, DP (Coverage), ExcessHet, FS, InbreedingCoeff,
  MQ (RMSMappingQuality), MQRankSum, QD, ReadPosRankSum, SOR.
* FORMAT: AD (DepthPerAlleleBySample), DP (DepthPerSampleHC).
* MLEAC and MLEAF come from the genotyping engine. GT, GQ and PL come from genotype assignment.
* InbreedingCoeff needs at least 10 samples, so it is never written.

## Which reads each annotation sees

For a locus at `[start, end]` (untrimmed merged-VC coordinates):

* **Genotyping evidence.** These are the realigned evidence reads whose new alignment overlaps
  `[start-2, end+2]`, clamped to the contig (`retainEvidence`). Their likelihoods are marginalised
  to the merged VC's alleles.
* **Annotation likelihoods.** These are the genotyping evidence plus the `filterNonPassingReads`
  reads whose original alignment overlaps the same window. The filtered reads get likelihood 0 for
  every allele, so they are never informative, but they count in DP and MQ. The alleles are all the
  merged VC's alleles, including those the call dropped.
* `filteredSampleEvidence` (the reads the likelihood engine dropped) is carried through
  `marginalize` but no default annotation reads it.
* An *informative* read is one whose best allele beats the second best by more than 0.2 log10.
  Ties within 0.2 are broken towards the reference allele.

| Annotation | Reads used |
|---|---|
| INFO DP | `evidenceCount()`: every read above, informative or not |
| MQ | All reads above with MQ ≠ 255: `sqrt(Σmq²/n)` formatted `%.2f` |
| FORMAT AD / DP | Informative reads after re-marginalising to the **call's** alleles. A read that preferred a dropped allele can count for a kept one |
| FS / SOR | Informative reads (best allele over all merged alleles) whose best is the reference or any emitted alt. 2×2 strand table: FS zeroes it unless the total is > 2; SOR uses its own table with min count 0 |
| BaseQRankSum, MQRankSum, ReadPosRankSum | Informative reads (all merged alleles) with MQ ∉ {0, 255}. Reference reads vs reads whose best allele is an emitted alt; reads preferring dropped alleles are ignored. ReadPos also needs the soft-clipped span to touch `[start, end+1]` |
| QD | `QUAL / depth`. The depth is the sum of AD, or the AD-restricted depth when more than 1 read supports an alt; it falls back to the evidence count when AD sums to 0 |
| AC / AF / AN, ExcessHet | The genotype only |

## GATK / htsjdk behaviours replicated (non-obvious)

### Realignment

1. **Best haplotype.** `bestAllelesBreakingTies(HAPLOTYPE_ALIGNMENT_TIEBREAKING_PRIORITY)` picks the
   best haplotype, with priority = isRef + (1 − #cigar elements). The priorities only matter when
   best − second < 0.2; then every haplotype within 0.2 of the best competes on priority. The port
   copies `searchBestAllele` exactly, including its initial `secondBestIndex = 0` quirk.
2. **SW input.** Soft-clipped bases are removed before SW (`hardClipSoftClippedBases`). The exact
   substring shortcut (last occurrence) of `SWNativeAlignerWrapper` is applied before the full SW.
   The SW cigar passes through `new CigarBuilder()`, which drops leading and trailing deletions.
3. **Haplotype cigar.** `getConsolidatedPaddedCigar(1000)` is the haplotype cigar + `1000M` through
   `CigarBuilder(true)`, so a leading deletion in the haplotype cigar is dropped. The read's
   reference start uses `readStartOnReferenceHaplotype`: an element that does not consume the
   reference contributes no excess.
4. **Projection to the reference.** `trimCigarByBases` and then `applyCigarToCigar`, where the
   read-D / haplotype-I pair emits nothing. `leftAlignIndels` runs against the **trimmed reference
   haplotype's** bases, not the padded reference. The start is shifted by the leading deletion
   bases it removes.
5. **Clips.** The original clip elements (S and H) are re-attached at both ends. The original read is
   kept when SW reports offset −1. The evidence order does not change.

### Events and merged alleles

6. **EventMap.** Events come from the haplotype cigar against the padded reference, starting at
   `alignmentStartHapwrtRef`.
   * Insertions in the first or last cigar element are ignored.
   * Deletions need `refPos > 0`.
   * "Regular" bases include lowercase and `*`.
   * With max MNP distance 0, adjacent SNPs stay separate.
   * Two events with the same start merge into a compound event (SNP+insertion → ref with longer
     alt; SNP+deletion → deletion ref with SNP alt; insertion+deletion). Every `Event` is reduced
     to its minimal representation (common suffix trimmed).
7. **Overlapping events.** `getOverlappingEvents(loc)` drops deletions ending at `loc` when an
   insertion overlaps `loc`.
8. **Locus set.** Event starts are the union over all haplotypes, restricted to the **genotyping
   span** (`regionForGenotyping`), not the active span.
9. **Event list per locus.** Events at a locus include those spanning it: they are de-duplicated as
   `Event` objects in haplotype order, and spanning ones become `ref-base → *` VCs.
10. **`simpleMerge`.**
    * The longest reference allele wins.
    * VCs with a shorter reference get their alts extended with the extra reference bases.
    * The allele `LinkedHashSet` receives such a VC's mapped alleles as **alts first, then the
      reference**. `VariantContext` then moves the reference to the front, so ALT order is
      first-appearance order over haplotypes.
11. **`*` is not symbolic.** In htsjdk, `Allele.SPAN_DEL.isSymbolic()` is false. So `*` is an ordinary
    allele of the merged VC and of the allele mapper, in merge order, and the genotype-likelihood
    allele order equals the VC order.
12. **`createAlleleMapper`.**
    * A haplotype with no overlapping event supports the reference.
    * An event starting at the locus with a shorter reference is remapped by extension.
    * A spanning event maps the haplotype to `*` and **breaks**, so the first spanning event wins.
13. **`removeAltAllelesIfTooManyGenotypes`** (more than 44 alleles).
    * Alleles are ranked by best, then second-best supporting haplotype score.
    * The reference haplotype's score is NaN: it sorts last with `Double.compareTo`, and `>`/`<`
      comparisons with it are false.
    * Final ties use `Allele.compareTo` (signed bytes, then length).
    * The ranking runs through a ported `java.util.PriorityQueue` (sift-up / sift-down), so even an
      inconsistent comparator gives the same result.

### Genotyping and QUAL

14. **Marginalisation and evidence window.** An allele's likelihood is the max over its haplotypes,
    or −∞ if it has none. The overlap window of ±2 uses the realigned start/end.
15. **Genotype likelihoods.** Hom = plain sum. Het = Σ `approximateLog10SumLog10` − n·log10(2).
    `approximateLog10SumLog10` uses the Jacobian table, which is **dumped from GATK**
    (`src/jacobian_java.bin`, 80,001 entries, via reflection), with `fastRound` indexing and an
    8.0 cutoff. PL = `Math.round(-10·(gl − max))`, with Java's round-half-up implemented exactly.
16. **AF calculator.**
    * It works from **PL/−10**, not from the raw likelihoods.
    * Pseudocounts: ref 0.001/0.01² = 10, SNP ×0.001, indel ×1.25e-4. Any allele whose length equals
      the reference's gets the SNP pseudocount, including `*` at a 1-bp reference.
    * The first EM iteration uses a flat prior; iteration stops once the change in counts is ≤ 0.1.
    * Diploid combination counts are 0 and log10(2) (the same bits as GATK's FastMath-based value).
    * With `*`, P(no variant) sums the 0/0, 0/\* and \*/\* posteriors, capped at 0.
    * Multi-allelic sites get a per-allele absent probability; biallelic sites take a shortcut.
    * MLE counts use `Math.round` of the EM counts.
17. **Emission.**
    * An allele is emitted iff log10 P(absent) + 1e-10 < 30·(−0.1).
    * A monomorphic site is not called.
    * QUAL = −10·log10 P(no variant), and the site needs QUAL ≥ 30. FILTER is therefore always `.`
      with default thresholds.
    * A site whose only emitted alt is `*` is dropped.
18. **`*` bookkeeping across regions.** A `*` allele is spurious, and dropped, unless an
    **upstream deletion recorded by an earlier emitted call** covers the locus (`start < loc ≤ end`).
    The recorded length is ref length − allele length for each emitted allele. This includes `*`
    itself when the reference is longer than 1 bp. The state persists across regions, as GATK's
    `upstreamDeletionsLoc` does.
19. **More than 6 alt alleles.** They are reduced with `calculateMostLikelyAlleles`: a likelihood
    sum for the alleles in the best genotype, then a stable descending sort. The reduced PLs are
    re-derived by subsetting and rounding before the AF calculation. The final genotype is subset
    from the original PLs.
20. **Genotype assignment.**
    * The original PLs are divided by −10, subset to the output alleles, rescaled to max 0 and
      rounded back to PLs.
    * GT = first maximum.
    * The call is no-call when it is hom-ref with GQ log10 > −0.1.
    * GQ = `Math.round(-10·gqLog10)` (= second-smallest PL). The cap at 99 is applied only by the
      VCF writer (`GQAccessor`).
21. **MLEAF** = min(1, MLEAC/AN), and AF = AC/AN, both formatted by `formatVCFDouble`.

### Annotation details

22. **FisherStrand.**
    * The table is normalised to 200 when its total exceeds 400 (int truncation).
    * The two-sided p-value uses commons-math's saddle-point hypergeometric `logProbability`, with
      Stirling-error table and deviance-part series ported from the bytecode.
    * Tables with probability ≤ observed·(1−1e-6) are summed via `log10SumLog10`.
    * p is floored at 1e-320; the output is `|−10·log10 p|` formatted `%.3f`.
23. **Mann-Whitney U.**
    * Ranks are `float`, and the rank sum is a float accumulation.
    * Tie correction excludes a tie group that spans all values.
    * The normal approximation (with a −0.5 continuity correction when there are ties) applies when
      either series has ≥ 10 values. Otherwise an exact permutation distribution is used and z =
      `NormalDistribution.inverseCumulativeProbability(p)` = √2·`Erf.erfInv(2p−1)`; the erfInv
      coefficients are taken from the jar's bytecode. The permutation distribution is computed by
      subset-sum counting instead of enumerating permutations; the result is identical.
    * A NaN z (an empty series) means the key is omitted.
24. **QD jitter.** QD ≥ 35 is replaced by 30 + 3·`nextGaussian()` from GATK's **global**
    `java.util.Random`, seeded with 47382911. The port reproduces the polar method, including
    `StrictMath.log`, which is ported from fdlibm `e_log.c`. Draws happen in call order across the
    whole run: 540 draws on this dataset, all matching. Positional downsampling, the only other
    consumer of that generator in HaplotypeCaller, never triggers at 30x (see `variant-rs`).
25. **ReadPosRankSum.**
    * A read starting right after the variant whose first non-clip element is an insertion gets 0.
    * Positions are counted from the soft-clip-inclusive start (`getReadIndexForReferenceCoordinate`
      treats S as reference-consuming).
    * Hard clips are added to both distances.
26. **BaseQRankSum.** The base quality at `vc.getStart()`; reads with a deletion there are skipped.
27. **ExcessHet.** The single-sample exact test (`ROUND_GENOTYPE_COUNTS`) always gives `0.0000`, but
    it is implemented in full.

### Trimming and VCF text

28. **`reverseTrimAlleles`.** It runs only when the call dropped alleles, and **after** annotation,
    so every annotation sees untrimmed alleles and positions. `*` is exempt. Nothing is trimmed if
    any allele has length 1. It uses `normalizeAlleles` semantics (forward trim computed, only
    reverse applied, one base restored if an allele would become empty).
29. **Number formatting.** Java `%.Nf` rounds **HALF_UP on the shortest decimal representation**
    (`FormattedFloatingDecimal`), so `1.005` → `1.01` and `0.125` → `0.13`, unlike C's `printf`.
    `-0.0` prints as `-0.000`.
    * `formatVCFDouble`: `%.3f` below 1, `%.2f` from 1 up, `%.3e` below 0.01.
    * QUAL: `%.2f` with a trailing `.00` removed.
    * GATK formats most INFO values itself as strings: FS, SOR and the rank sums `%.3f`; MQ and QD
      `%.2f`; ExcessHet `%.4f`.
30. **Record layout.**
    * INFO keys are written through a `TreeMap`, i.e. ASCII-sorted (`BaseQRankSum` < `DP` …
      `ReadPosRankSum` < `SOR`).
    * FORMAT is GT first, then the other keys sorted: `GT:AD:DP:GQ:PL`.
    * ID and FILTER are `.`.

## Not ported / differences

* **Off-by-default modes are not ported:** DRAGEN BQD/FRD
  genotyping, DRAGstr priors, contamination downsampling, allele filtering, stepwise filtering,
  PDHMM, given/force-called alleles, `--use-filtered-reads-for-annotations`,
  `--annotate-with-num-discovered-alleles`, flow-mode annotations, the bamout writer, and ploidy
  other than 2. Only a single sample is supported.
* **Elementary functions are not bit-identical to the JVM's.** The port uses Rust/glibc
  `log10`, `pow` and `ln`, but HotSpot's `Math.log10` / `Math.pow` intrinsics and commons-math's
  `FastMath.log` sometimes differ in the last bit. `numcheck` measured, over 1,000,000 random
  arguments each:

  | JVM function vs Rust | Last-bit differences |
  |---|---:|
  | `Math.log10` | 3.3% |
  | `Math.pow(10, y)` | 0.08% |
  | `FastMath.log` | 0.011% |
  | `Math.log` | 0.008% |

  These functions feed internal values (AF-calculator posteriors, FS, SOR). A 1-ulp difference
  reaches the text only if a value falls within about 1e-13 of a rounding boundary of
  `%.2f`/`%.3f`, of the QUAL ≥ 30 or 1e-10 thresholds, or of an MLE `.5`. That never happened on the
  8,064 real or 14,545 synthetic records. Matching it bit for bit would need a port of HotSpot's
  Intel-LIBM `log10`/`pow` stubs and of commons-math's `FastMath.log`.

## Reference-confidence mode (`-ERC GVCF`)

`src/erc.rs` adds `GenotypingEngine::genotype_region_erc` (with `EngineConfig::new_erc`). It ports the ERC
branches of `assignGenotypeLikelihoods` and `calculateGenotypes` (`<NON_REF>` allele and likelihoods, confidence
-0.0, all-site PLs, `forceKeepAllele`), the ERC annotation set (raw `RAW_MQandDP`, `SB`; no AC/AF/AN, FS, SOR, QD)
and physical phasing (PGT/PID/PS). The quirks are listed in [`../refconf/README.md`](../refconf/README.md). The reference model and the
GVCF writer are in [`../refconf`](../refconf).

* `genotype-validate --erc --geno $VARCALL_DATA/gvcf/gatk_genotyping_erc.jsonl --asm $VARCALL_DATA/gvcf/gatk_assembly_erc.jsonl --vcf <variant
  lines of $VARCALL_DATA/gvcf/gatk_raw_calls.vcf.gz>` gives 8,905 / 8,905 identical call lines.
* `ERC=1 scratch/fuzz.sh SEED K` (development script, not in the repository) runs `RunGeno.java` with GATK's ERC settings (`ERC=1`) on mutated regions. Seeds
  101–105 with K = 3 or 5 gave 14,470 identical calls, including max-alt reduction, `removeAltAllelesIfTooManyGenotypes`
  and compound events. Their VCFs were in `scratch/runs/erc_*.vcf` (not in the repository).
* `VcfRecord` gained `sb` and `phase` (both `None` in VCF mode). `reverse_trim_alleles` leaves symbolic alleles
  untouched. Neither change affects VCF mode.

## Files

* `src/realign.rs`: realignment.
* `src/likelihoods.rs`: `searchBestAllele`.
* `src/events.rs`: EventMap, merged VC, allele mapper.
* `src/genotyper.rs`: genotype likelihoods, AF calculator, subsetting, genotype call.
* `src/engine.rs`: per-locus driver, cross-region state, annotation assembly, reverse trim.
* `src/erc.rs`: reference-confidence (GVCF) genotyping, annotations and physical phasing.
* `src/annotate.rs`: FS, SOR, Mann-Whitney U, read helpers, ExcessHet.
* `src/javamath.rs`: Java `Math.round`, `%.Nf` formatting, `java.util.Random`, fdlibm log, erfInv,
  saddle-point hypergeometric.
* `src/mathutils.rs`: GATK `MathUtils`.
* `src/vcf.rs`: record formatting.

Development-only, not in the repository:

* `scratch/probes/`: Java probes.
  * `Tables.java` dumps the Jacobian table and checks htsjdk behaviours.
  * `RunGeno.java` runs GATK's genotyping on a dump.
  * `NumProbe.java` writes the numeric test vectors.
* `scratch/mutate.py`, `scratch/fuzz.sh`: synthetic differential tests. The large generated JSONL
  files were deleted; they were regenerated with `scratch/fuzz.sh SEED K`. Their VCFs were kept in
  `scratch/runs/`.
* `scratch/*.javap`: disassembly of the htsjdk and commons-math methods that were ported.
