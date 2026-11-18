# hc_downsample: GATK HaplotypeCaller read downsampling, exact

This crate ports GATK 4.7 HaplotypeCaller's `--max-reads-per-alignment-start` downsampling (default 50) to Rust.
It also ports the way the downsampler shares GATK's global `java.util.Random` with the QD jitter.

The port has been validated against instrumented GATK runs. Kept-read set and order, RNG draw sequence and
jittered QD values are all identical, on real 30x WGS data at the default cap and on stress data (caps 1–3,
synthetic depth, two samples, interval gaps, coverage holes).

| File | Content |
|---|---|
| `src/javarandom.rs` | `JavaRandom`: bit-exact `java.util.Random` (`nextInt(bound)` with its rejection loop, `nextGaussian` with its cache). `strict_log` and `format_fixed` are copied from `hc_genotype::javamath`. |
| `src/downsampler.rs` | `PositionalDownsampler<T>` (one per contig shard), `reservoir_downsample`, and `horizon()` (GATK's read-pull lookahead). |
| `src/schedule.rs` | The interleaving rule: `ReadyTracker`, `lookahead_locus`, `pileup_span`, `LibsCoverage`. |
| `src/bin/validate.rs` | Validator against an instrumented GATK run. |
| `src/bin/rngcheck.rs` + `probes/RngProbe.java` | `JavaRandom` vs JDK 17 `java.util.Random`. |
| `probes/patched/` | GATK sources with logging added. They are compiled to `probes/classes` and put ahead of gatk.jar. |
| `probes/run_gatk_probe.sh` | Runs HaplotypeCaller with the instrumented classes. |
| `probes/make_synthetic.sh` | Builds the stress BAMs. |
| `integration/hc_genotype_rng_hook.patch` | Patch to `hc_genotype` exposing the QD generator's state (applied). |
| `scratch/` | Test BAMs, GATK logs and VCFs from the validation runs. |

Only `src/` is in the repository: `probes/`, `integration/` and `scratch/` were development-only and are not included.

## 1. GATK semantics (source: GATK 4.7.0.0)

**Creation.** `AssemblyRegionWalker.createDownsampler()` returns
`maxReadsPerAlignmentStart > 0 ? new PositionalDownsampler(max, header, nonRandomDownsamplingMode) : null`.
- The default is `AssemblyRegionArgumentCollection.DEFAULT_MAX_READS_PER_ALIGNMENT = 50`. A value of 0 disables downsampling.
- `nonRandomDownsamplingMode()` is `false` in `AssemblyRegionWalker`, and HaplotypeCaller does not override it.
  Random mode is therefore used. (The non-random mode would use `|name.hashCode()| % seen` instead of the RNG.)

**Attachment.** One shard is a `MultiIntervalLocalReadShard` per contig, holding all of that contig's intervals.
- `traverse()` sets up each shard with:
  - an identity pre-filter transformer;
  - the counted HaplotypeCaller read filter;
  - a **new** downsampler, so one per contig;
  - an identity post-filter transformer.
- `MultiIntervalLocalReadShard.iterator()` wraps the reads in this order:
  1. the source: reads overlapping the intervals padded by 100 (`assemblyRegionPadding`), queried together, no duplicates;
  2. the pre-filter transformer;
  3. `ReadFilteringIterator`;
  4. the post-filter transformer;
  5. `ReadsDownsamplingIterator`.

  **So the downsampler sees the post-filter reads, in BAM order.** Filtered reads (MQ < 20, duplicates,
  secondary, QC-fail, unmapped including mate-placed unmapped reads, bad CIGAR, malformed) never count toward the cap.
  Supplementary alignments (0x800) are not filtered, so they do count: flags 2131/2147 appear in the kept streams.
- `ReadsDownsamplingIterator` is a `PushToPullIterator`. Downsampling is the last stage before
  `ReadCachingIterator` → `LocusIteratorByState`. Everything downstream sees only the kept reads, in kept order:
  - isActive pileups;
  - region read lists, including the reads carried over from the previous region.

  The driver must therefore feed its pileup accumulators and its region read cache in kept order, not BAM order.

**Semantics** (`PositionalDownsampler`, `ReservoirDownsampler`):
- **Not per sample, not per read group.** There is one reservoir per alignment start for all reads.
  Validated with a 2-sample BAM (`syn.bam`, a development test BAM not in the repository; read group SYN2 with SM SYNTH2).
- **Position** = `ReadCoordinateComparator.compareCoordinates`, which compares (assigned contig index, assigned start).
  For mapped reads this is BAM POS, 1-based, with soft clips not counted. Reads with no assigned position
  would pass through, but `MappedReadFilter` removes them first.
- **Reservoir** = Vitter's Algorithm R:
  - the first `target` reads of a position fill slots 0..target-1;
  - read number `k > target` draws `Utils.getRandomGenerator().nextInt(k)` and replaces slot `s` if `s < target`.

  A position with `n` reads therefore makes `n − target` `nextInt` calls if `n > target`, and none otherwise.
  This README calls such a position a *drawing group*.
- **Output order.** Positions are output in coordinate order. Within a position, reads come out in **reservoir slot
  order**: a replacing read takes the evicted read's list index. This is not BAM order once a position overflows.
- **Finalization.** A position is finalized when the first read of the next position is submitted
  (`handlePositionalChange` → `finalizeReservoir`: `signalEndOfInput` + `consumeFinalizedItems`, which resets
  `totalReadsSeen`). End of the shard: `PushToPullIterator.fillCache` calls `signalEndOfInput` (called again on an
  empty reservoir; harmless) and the last position is finalized.
- **Unsorted input** throws `IllegalStateException`. The Rust port asserts instead.
- **Shard and interval boundaries.** The downsampler is recreated for each contig. Within a contig, the reads from all
  padded intervals form one stream: an interval gap resets nothing. Reads that start in the padding or before the
  interval are downsampled like any other read.
- **Draw timing.** `ReservoirDownsampler.submit` draws at submit time. The draws of one position all happen inside a
  single `fillCache` call, as one contiguous block. That call submits the rest of position P plus the first read of P+.

**No other downsampling on the default path** (checked in the code and by the probes):
- `AssemblyRegionIterator` and `PileupReadErrorCorrector` build `LocusIteratorByState` with `DownsamplingMethod.NONE`.
- Contamination downsampling (`AlleleBiasedDownsamplingUtils`) is skipped. `HaplotypeCallerGenotypingEngine` calls
  `contaminationDownsampling` only if `isSampleContaminationPresent()`, which is false for the default
  `CONTAMINATION_FRACTION = 0.0` with no file. Even when on, it draws from `Utils.getRandomDataGenerator()`
  (a separate Well19937c), not from the `Random`.
- `LevelingDownsampler` and `maxDepthPerSample` belong to LocusWalker and are not used here.

## 2. The random generator

- `Utils.randomGenerator = new Random(47382911)` is a JVM-wide static. It is reset only by
  `FilterMutectCalls` (`Utils.resetRandomGenerator`), never on HaplotypeCaller's path.
- The probe replaced it with a subclass that logs every LCG step `next(bits)` together with its first non-JDK caller.
  Across all runs, the **only** callers are:
  - `ReservoirDownsampler.submit`: `nextInt`;
  - `QualByDepth.fixTooHighQD`: `nextGaussian`, used when raw QD ≥ 35 to give `30 + 3·N(0,1)`.
- **GVCF mode (`-ERC GVCF`) makes no QD draws** (run g1: 593 steps, all from the reservoir). Only the downsampler uses the RNG there.
- **The generator is global across contigs.** Shards run sequentially, in sequence-dictionary order, and the state carries over.
  The draws on chr1 shift every later contig, for downsampling and for QD jitter alike. Runs w1 and w3 cover
  multiple contigs in one run and validate this continuity.
  This already matters without downsampling: QD values on chr20 in a whole-genome run depend on how many jitter
  draws chr1–chr19 made.
- `nextGaussian` caches its second value. A cached value survives intervening `nextInt` calls, so it must be
  modelled. `JavaRandom` does this. `rngcheck` compared 1.2 M mixed calls against JDK 17 (4 seeds; bounds that
  force `nextInt`'s rejection loop) and found 0 mismatches.

## 3. When the draws happen relative to `apply()` (QD jitter)

GATK pulls reads lazily. The chain is:

```
PushToPullIterator(downsampler) -> ReadCachingIterator -> LIBS (ReadStateManager: PeekableIterator)
 -> IntervalAlignmentContextIterator -> AssemblyRegionIterator.loadNextAssemblyRegion -> apply()
```

1. `AssemblyRegionIterator.next()` returns R_i only after it has loaded R_{i+1} (`readyRegion = loadNextAssemblyRegion()`).
   When `apply(R_i)` runs, the locus loop therefore stands at A_{i+1}, the locus where R_{i+1} became ready.
   - A region becomes ready at the first locus L, with at most one region per locus, such that it heads
     `pendingRegions` and `L > paddedSpan.end`.
   - Pops happen at the top of the loop body, before the locus's own state is added.
2. After returning the pileup at A, `IntervalAlignmentContextIterator` advances to N, the next interval locus.
   It pulls LIBS contexts until it holds one at or after N. That context is at C, the first locus ≥ N with a
   **non-empty LIBS pileup**: deletions count, adaptor-excluded bases do not.
   - If no interval locus is left (A is the last one, or R_{i+1} came from the out-of-loci branch), the iterator
     drains LIBS (`while (libs.hasNext()) libs.next()`), so the whole shard has been pulled.
   - If R_{i+1} does not exist, the whole shard has also been pulled.
3. To produce C, LIBS consumed every kept read with start ≤ C, say m of them.
   - Its `PeekableIterator` already holds kept read k_{m+1}.
   - `PushToPullIterator` already holds k_{m+2}. Getting it required finalizing k_{m+2}'s position.

**Rule.** The draws of drawing group P happen before `apply(R_i)` if and only if P ≤ X_i, where:
- X_i = ∞ if R_{i+1} is absent, ready at End, ready at the last interval locus, or no pileup exists at or after N;
- otherwise X_i = the start of k_{m+2}. That is P1, the first start after C, if P1 keeps ≥ 2 reads
  (kept = min(n, target), which does not depend on the RNG). Otherwise it is P2, the next start after P1.

X_i is non-decreasing. Within `apply(R_i)`, QD draws happen in the genotyping engine's call order, which hc_genotype already reproduces.

Evidence that the P1/P2 lookahead is exact and not approximate. Deliberately wrong rules were run through the same
validator; "FAIL" means the kept stream or the QD values differ from GATK:

| run | data | exact | X=C | X=P1 always | X=P2 always | no interleave |
|---|---|---|---|---|---|---|
| t2 | chr20 5 Mb, cap 1 (P1 always keeps 1 read) | PASS | FAIL | FAIL | PASS | FAIL |
| s3 | synthetic 6x copies, cap 3 (P1 keeps ≥ 2) | PASS | FAIL | PASS | FAIL | — |
| w1 | WGS chrM:1-1500 + decoy, cap 50 | PASS | FAIL | PASS | FAIL | FAIL |
| w2 | WGS chr20 2.0–2.5 + 30.9–31.3 Mb, cap 50 | PASS | FAIL | PASS | PASS | FAIL |
| w3 | WGS chrM + chrUn_KI270467v1 + JTFH01001899 decoy, cap 50 | PASS | FAIL | PASS | FAIL | FAIL |

Only the exact rule passes every run.

## 4. Integration recipe for the Rust driver

**One generator for the whole run.**
- Create one `JavaRandom::gatk()`, shared by the downsampler and the QD jitter, for all contigs in GATK order
  (sequence-dictionary order, all contigs when there is no `-L`).
- `hc_genotype::GenotypingEngine` keeps its own private `rng`. `integration/hc_genotype_rng_hook.patch` (development-only,
  not in the repository) added `rng_state` / `set_rng_state` to it; the change is applied here (`../genotype/src/engine.rs`).
- Around every `genotype_region` / `call_region`, write the shared state into the engine and read it back:

  ```rust
  engine.set_rng_state(rng.state().to_tuple());
  /* genotype */
  rng.set_state(JavaRandomState::from_tuple(engine.rng_state()));
  ```

**Per contig shard:**
```
ds = PositionalDownsampler::new(50); cov = LibsCoverage::new(); tracker = ReadyTracker::new()
for each read passing the HC filters, in BAM order (reads overlapping the intervals padded by 100):
    ds.push(read.start, read)
at end of the shard's reads: ds.finish()

To process loci <= L (all reads with start <= L are needed):
    loop:
        match ds.finalize_free_upto(L, &mut kept):         // groups with n <= target: no RNG use
            None => break
            Some(P) => {                                     // drawing group P
                process all loci < P (pileups are complete for them): isActive, pops, readiness, horizons
                genotype, in order, every pending region i whose X_i is known and X_i < P
                ds.finalize_head(&mut rng, &mut kept)          // P's draws
            }
    for r in kept, in order: cov.add(pileup_span(r)); feed pileup accumulators and the region read cache

Readiness: per interval locus l, ReadyTracker::on_locus(l, padded ends popped at l) -> region id ready at l.
After the last locus: on_end(final forced pops) -> those regions are ReadyAt::End.

Horizon of region i, once ready_at[i+1] is known (or R_i is known to be the last region):
    N = lookahead_locus(ready_at[i+1], intervals)
    C = N.and_then(|n| cov.first_covered(n, ds.first_pending_start() - 1))
        // None while groups are still pending means: not known yet. Release more groups, then retry.
        // None with nothing pending and ds finished means LIBS is exhausted: X = All.
    X_i = ds.horizon(C)          // NeedInput: push more reads first

Genotyping region i (requires X_i known, and regions < i already genotyped):
    ds.finalize_upto(X_i, &mut rng, &mut kept)   // or finalize_all for Horizon::All; feed `kept` as above
    set engine rng state; genotype; read back
```

Why this can neither deadlock nor reorder draws:
- Whether X_i < P can be decided from loci < P and from raw read starts.
  - A_{i+1} < C_i < X_i < P, and pops and readiness at a locus depend only on earlier loci.
  - The P1/P2 kept counts are min(n, target).
- So if a horizon is still unknown when group P is reached, it is ≥ P, and that region is correctly genotyped after P's draws.

**Cost for a batching driver.** It must flush, i.e. genotype the pending regions, before each drawing group that
follows unflushed regions. Neighbouring drawing groups need no flush between them. Drawing groups are rare and
clustered on 30x WGS at cap 50:

| contig | positions with > 50 post-filter reads | `nextInt` draws | clusters (gap > 2 kb) |
|---|---:|---:|---:|
| chr20 | 79 | 4,499 | 23 |
| chr1 | 3,763 | 236,048 | 42 (mostly pericentromeric, 125 Mb and 143 Mb) |
| chrM | 16,297 of 16,530 starts | 590,411 | 1 (whole contig) |

These counts are approximate: `samtools view -F 0x704 -q 20` with `awk`, without the CIGAR, wellformed and MQ255 filters.
Most regions can still be batched; chrM and high-copy decoys are effectively sequential.

**Sharded parallelism.** The draw chain is genome-wide.
- In VCF mode, the state at a drawing group on contig k depends on every QD draw and every group draw before it,
  on all earlier contigs.
- In GVCF mode it depends on group draws only.
  - Their count is fixed by the input: Σ(n − target).
  - The LCG step count is not quite fixed: `nextInt` rejects with probability (2³¹ mod bound)/2³¹, about 10⁻⁶ per call here.
    The runs here saw no rejections: 662,173 `nextInt` calls took 662,173 steps.
- Exact options:
  - (a) Contigs in RNG order, with parallel work inside each contig between flushes.
  - (b) Speculate on contigs in parallel, then rerun each contig from its first drawing group once the true
    incoming state is known. Contigs without drawing groups only need their QD values recomputed in a final
    sequential pass, because QD jitter values do not feed anything else.

  A per-contig fresh generator is not GATK-identical once any contig before it drew.

### How variant-rs integrates it

`variant-rs` (see its README) follows this recipe with one change: the genotyping engine does not draw the QD jitter
itself (`EngineConfig::defer_qd_jitter`); it writes a placeholder per jittered QD. The placeholders of region i are
drawn from the shared generator when the schedule says so: before drawing group P for every region with X_i < P
(the front end asks the back end to genotype and draw those first), the rest at the end of the contig. Contigs of
one RNG domain (the run, or one shard file) hand the generator on in dictionary order (`variant-rs/src/rng.rs`).
`LibsCoverage` is replaced by an equivalent heap-based structure (`variant-rs/src/frontend.rs`, `Coverage`): the
`Vec`-based one rescans every stored span per query, which is quadratic over a whole contig.
variant-rs reproduces w1, w2, w3, t1, g1 and t2 byte for byte (`ds_regress.sh`, a development script not in the repository); s3 has two samples, which
variant-rs does not support.

## 5. GATK quirks worth knowing

- Kept order within an overflowing position is slot order, not BAM order. This reorders LIBS pileup elements,
  which changes the floating-point sum order of isActive, and it reorders region read lists.
- The lookahead depends on whether the next start keeps one read or several: `PeekableIterator` plus
  `PushToPullIterator`, two reads ahead. On real data at cap 50 the P1/P2 distinction decides the outcome
  (w1 and w3 fail with "always P2").
- C counts deletions in pileups but not adaptor-clipped bases. It uses the signed-TLEN `> 100` quirk, which never
  applies to reverse reads, as in variant-rs (`pileup_span`).
- The readiness check is one region per locus. The out-of-loci branch drains LIBS, so the last regions of a shard see everything drawn.
- Supplementary alignments count toward the cap. Duplicates and MQ < 20 reads do not.

## 6. Validation

All runs use `probes/run_gatk_probe.sh` (development probe, not in the repository), i.e. gatk.jar with 6 patched classes that only add logging. Plain gatk.jar
and the instrumented run gave byte-identical VCF data lines on t1.

The validator does four things:
1. It re-derives the read stream from the BAM itself: `samtools view -M` over the padded intervals, then the HC filters.
2. It checks that this stream equals GATK's downsampler input.
3. ORACLE: it replays GATK's logged interleaving through the Rust downsampler and RNG.
4. RULE: it computes the interleaving with §3's rule from the BAM, from GATK's region pops (the driver's own
   activity profile provides these in production) and from the number of QD draws per region.

It compares the kept stream (read name, flag, start, in order), every jittered QD value (bit-exact, and as the
VCF's `%.2f` text), the total LCG step count and `ReadyTracker` against GATK's READY loci.

| run | input | cap | reads in → kept | drawing groups | regions | jittered QD | result |
|---|---|---:|---|---:|---:|---:|---|
| t1 | chr20:10.0–10.3M (30x) | 2 | 63,848 → 63,255 | 480 | 1,511 | 35 | ALL PASS |
| g1 | same, `-ERC GVCF` | 2 | 63,848 → 63,255 | 480 | 1,512 | 0 (no QD in GVCF) | ALL PASS |
| t2 | chr20:10–15M (30x) | 1 | 1,058,307 → 948,581 | 100,517 | 24,541 | 537 | ALL PASS |
| s3 | syn.bam: 6 copies per read, 2 samples, 4 intervals with gaps | 3 | 248,070 → 111,288 | 37,096 | 1,019 | 31 | ALL PASS |
| s2_1 / s2_2 | syn2.bam: syn.bam with 400 bp start-free holes every 7.9 kb, interval past the end of coverage | 1 / 2 | 243,222 → 36,373 / 72,746 | 36,373 | 1,253 / 1,246 | 32 / 32 | ALL PASS |
| w1 | WGS chrM:1-1500 + chrUn_JTFH01001929v1_decoy (2 shards) | 50 | 197,335 → 104,408 | 1,812 | 28 | 2 | ALL PASS |
| w2 | WGS chr20:2.0–2.5M + 30.9–31.3M | 50 | 260,297 → 256,521 | 65 | 4,789 | 201 | ALL PASS |
| w3 | WGS chrM (full) + chrUn_KI270467v1 + chrUn_JTFH01001899v1_decoy | 50 | 1,584,628 → 922,455 | 17,777 | 115 | 36 | ALL PASS |

The horizon cases exercised were:
- C = N;
- C > N (coverage holes; 368 times in w2);
- no pileup after N (199 times in s2);
- stream drained;
- the largest group: 5,161 reads at one chrM start.

Reproduce:

```bash
# from the repository root, in an environment with htslib (environment.yml at the repo root)
RUSTFLAGS="-C target-cpu=native" cargo build --release -p hc-downsample && cargo test --release -p hc-downsample
# probes/ is a development directory, not in the repository
java probes/RngProbe.java 300000 | target/release/rngcheck
probes/run_gatk_probe.sh t2.log $VARCALL_DATA/data/chr20.fa $VARCALL_DATA/data/NA12878.chr20_10-15M.bam t2.vcf.gz \
    -L chr20:10000000-15000000 --max-reads-per-alignment-start 1
target/release/downsample-validate --log t2.log --bam $VARCALL_DATA/data/NA12878.chr20_10-15M.bam -L chr20:10000000-15000000 \
    --target 1 --vcf t2.vcf.gz [--variant C|P1|P2|none] [--dump-kept kept.tsv]
```

The synthetic BAMs were built by `probes/make_synthetic.sh`, run from `scratch/` (development files, not in the repository).
- `syn.bam`: each record of chr20:9.999–10.201M copied 6 times; odd copies moved to read group SYN2 (sample SYNTH2).
- `syn2.bam`: `syn.bam` with coverage holes cut into it.

The interval lists were `scratch/syn.intervals` and `scratch/syn2.intervals` (not in the repository).
