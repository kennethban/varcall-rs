# hc-assembly — GATK 4.7 HaplotypeCaller local assembly in Rust

Port of `AssemblyBasedCallerUtils.assembleReads` (everything after `finalizeRegion`) and
`ReadThreadingAssembler.runLocalAssembly` with **default HaplotypeCaller arguments**.
Input: the finalized reads of one active region plus the padded reference window.
Output: GATK's `AssemblyResultSet.getHaplotypeList()` — same haplotypes, cigars, order, kmer
sizes and scores.

## Status

Validator on `$VARCALL_DATA/asm/gatk_assembly.jsonl` (6,345 regions of NA12878 chr20:10–15 Mb):

```
regions:                         6345
haplotype set matches exactly:   6345
  ...and order matches:          6345
  ...and kmer+score match too:   6345
reference haplotype mismatches:  0
total haplotypes GATK / Rust:    29642 / 29642
missing / extra:                 0 / 0
kmer / score mismatches:         0 / 0
assembly time (excl. JSON):      ~4.6 s  (single thread, -C target-cpu=native)
```

Scores are bit-identical; the validator's 1e-9 tolerance is never needed.

Additional differential tests against GATK, run on the same jar, cover code paths the real data
never reaches (see "Testing" below): 3,500 synthetic regions (all identical, including one
region where GATK throws and the port throws too), 120,000 Smith-Waterman alignments against
GKL AVX-512 plus 40,000 `calculateCigar` calls (0 mismatches), and 20,000 random seq-graphs
through cleanup, simplification and k-best (0 mismatches).

## API

```rust
use hc_assembly::{assemble_region, Assembler, AssemblerConfig, Read, RegionInfo, AssemblyResult, Haplotype};

// one-shot, default parameters
let res: AssemblyResult = assemble_region(&reads, region, ref_padded, ref_padded_start, ref_padded_end);

// reusable (keeps Smith-Waterman buffers); one Assembler per thread
let mut asm = Assembler::new(AssemblerConfig::default());
let res = asm.assemble_region(&reads, region, ref_padded, ref_padded_start, ref_padded_end);
let res: Result<AssemblyResult, String> = asm.try_assemble_region(&reads, region, ref_padded, s, e);
```

* `Read { name, flag, start, end, cigar: Vec<(u8 op, u32 len)>, seq, qual (raw phred, not +33), mapq, mate_start, tlen }`.
  Pass the reads exactly as `finalizeRegion` leaves them (`region.getReads()`), in that order.
  Only `cigar` (for soft-clip removal), `seq` and `qual` affect the result.
* `RegionInfo { start, end, padded_start, padded_end }` (1-based, inclusive). Only the padded span is used.
* `ref_padded`, `ref_padded_start`, `ref_padded_end`: `getFullReferenceWithPadding()` and
  `getPaddedReferenceLoc()`, i.e. padded span ± `REFERENCE_PADDING_FOR_ASSEMBLY` (500), clipped to the contig.
* `AssemblyResult { ref_haplotype, haplotypes }`. `Haplotype { bases, cigar, align_start_wrt_ref,
  loc_start, loc_end, is_ref, score (NaN = GATK NaN, may be -inf), kmer }`.
  `haplotypes` is `getHaplotypeList()` and contains the reference haplotype.
* `try_assemble_region` returns `Err` where GATK throws (for example "Graph must have ref source and sink
  vertices"). An exception there aborts a GATK run. `assemble_region` panics instead.

Binaries:

* `assembly-validate <jsonl> [--limit N] [--detail N] [--only IDX] [--verbose]` — compares against an
  AssemblyDumper JSONL. Set `VARIANT_RS_ASM_PROFILE=1` for per-phase timings and code-path hit counters.
* `swcheck <tsv>` and `graphcheck <jsonl>` — component checks against the outputs of the Java probes listed under "Testing" (development probes, not in the repository).

Build (from the repository root, in an environment with htslib; `environment.yml` at the repo root): `RUSTFLAGS="-C target-cpu=native" cargo build --release -p hc-assembly`

## What was ported (default path only)

| GATK | Rust |
|---|---|
| `AssemblyBasedCallerUtils.assembleReads` (ref haplotype via `ReferenceConfidenceModel.createReferenceHaplotype`, padded ref loc) | `Assembler::assemble_region` |
| `ReadThreadingAssembler.runLocalAssembly`, `assemble`, `createGraph`, `getAssemblyResult`, `cleanupSeqGraph`, `findBestPaths`, `assembleKmerGraphsAndHaplotypeCall` | `lib.rs` |
| `ReadClipper.hardClipSoftClippedBases` (bases and quals only) | `hard_clip_soft_clipped` |
| `AbstractReadThreadingGraph`, `ReadThreadingGraph`, `MultiDeBruijnVertex`, `MultiSampleEdge` (numPruningSamples = 1) | `rtgraph.rs` |
| `ChainPruner` + `LowWeightChainPruner` | `ReadThreadingGraph::prune_low_weight_chains` |
| dangling tail recovery, legacy dangling head recovery (`mergeDanglingHeadLegacy`, `extendDanglingPathAgainstReference`) | `rtgraph.rs` |
| `BaseGraph` (ref source/sink, `getNextReferenceVertex`, `cleanNonRefPaths`, `removePathsNotConnectedToRef`, …) on top of jgrapht 1.1.0 `DefaultDirectedGraph` semantics | `graph.rs` (generic arena graph), `seqgraph.rs` |
| `SeqGraph`: `zipLinearChains`, `simplifyGraph` (incl. `graphEquals` loop guard), `MergeDiamonds`, `MergeTails`, `SplitCommonSuffices`/`CommonSuffixSplitter`, `MergeCommonSuffices`/`SharedSequenceMerger`, `SharedVertexSequenceSplitter` | `seqgraph.rs` |
| `GraphBasedKBestHaplotypeFinder` + `KBestHaplotypeFinder.removeCyclesAndVerticesThatDontLeadToSinks`, `KBestHaplotype` scores | `kbest.rs` |
| `SmithWatermanJavaAligner` (all overhang strategies), the `lastIndexOf` exact-match shortcut of `SWNativeAlignerWrapper` | `sw.rs` |
| `CigarUtils.calculateCigar`, `AlignmentUtils.trimCigarByBases`, `leftAlignIndels`, `normalizeAlleles`, `removeTrailingDeletions`, `CigarBuilder`, `Haplotype.setCigar` | `cigar.rs`, `cigar_utils.rs` |

Not ported (off by default, verified in the argument collections): read error correction
(`--error-correct-reads`, pileup error correction), linked de Bruijn graphs / JunctionTree
k-best finder, adaptive pruning, flow-based hmer collapsing, PDHMM/pileup haplotypes,
DRAGEN modes, `--recover-all-dangling-branches`, non-legacy dangling head merging
(`--num-matching-bases-in-dangling-end-to-recover >= 0`), debug graph output. Only one
sample is supported: reads go into a single pending-sequence group.

## Defaults assumed (checked in the 4.7.0.0 sources)

* kmer sizes `[10, 25]`; `dontIncreaseKmerSizesForCycles=false` → retry 35, 45, …, 85 (6 attempts) if
  **no** graph was produced; the last attempt allows low-complexity graphs and non-unique ref kmers.
* `allowNonUniqueKmersInRef=false`, `maxNumHaplotypesInPopulation=128`, `numPruningSamples=1`.
* `useAdaptivePruning=false` → `LowWeightChainPruner(minPruneFactor=2)`; `pruneFactor=2` also for dangling ends.
* `enableLegacyGraphCycleDetection=false` → prune *before* the cycle check.
* `minDanglingBranchLength=4`, `recoverDanglingBranches=true`, `recoverAllDanglingBranches=false`,
  `minMatchingBasesToDanglingEndRecovery=-1` (legacy head merging).
* `minBaseQualityToUseInAssembly = --min-base-quality-score = 10` (set by `createReadThreadingAssembler`,
  not the assembler's internal default of 10, which is the same value).
* `useLinkedDeBruijnGraph=false` → SeqGraph + `GraphBasedKBestHaplotypeFinder`.
* Error correction off (`errorCorrectReads=false`, `pileupErrorCorrectionLogOdds=-inf`);
  `flowAssemblyCollapseHKerSize=0` (no hmer collapsing).
* SW parameters: dangling ends `STANDARD_NGS (25,-50,-110,-6)`, haplotype→ref `NEW_SW_PARAMETERS (200,-150,-260,-11)`;
  overhang strategies: dangling ends `LEADING_INDEL`, haplotypes `SOFTCLIP` with `INDEL` fallback check.
* `REFERENCE_PADDING_FOR_ASSEMBLY=500`, `MIN_HAPLOTYPE_REFERENCE_LENGTH=30`, `SW_PAD="NNNNNNNNNN"`.

## GATK behaviours replicated (non-obvious)

Ordering and data structures

1. **jgrapht 1.1.0 iteration order.** `vertexSet()`/`edgeSet()` are `LinkedHashMap`s and per-vertex
   edge sets are `ArrayUnenforcedSet`s, so everything iterates in insertion order and removal keeps the
   order of what is left. The arena graph assigns ids at insertion, so id order is insertion order. Edge
   order decides the vertex order in `outgoingVerticesOf`, the splitter's `middles`, zip starts, ties
   in "heaviest edge", and so on.
2. **No multi-edges.** `addEdge(u,v,e)` silently returns false when `u→v` already exists. This matters in
   `SharedSequenceMerger` (two predecessors with the same source drop one edge and its multiplicity),
   in splitter reconnection and in dangling-end merges. `addEdge` with an edge object already in the graph
   also returns false (the shared prefix→suffix edge of `SharedVertexSequenceSplitter`).
3. **New vertices go to the end of `vertexSet()`.** `VertexBasedTransformer` restarts from the first
   vertex after every successful transform. `SplitCommonSuffices` keeps its `alreadySplit` set for one
   `transformUntilComplete` call only.
4. **`java.util.PriorityQueue` is ported exactly** (siftUp/siftDown), so k-best ties resolve the same way.
   The comparator is score descending, then `BASES_COMPARATOR` (signed bytes, then length) **reversed**.
5. **`Math.log10` differs from glibc** in the last ulp for about 1.6% of integers. k-best scores use a table of
   Java's own `Math.log10(0..131071)` (`src/log10_java.bin`, generated by `Log10Dump.java`, a development probe not in the repository).
   Without it, scores are off by 1 ulp and tied paths can swap order.

Read threading graph

6. `determineNonUniqueKmers` scans each sequence **from index 0** to `stop-k`, not from its `start`.
   Kmers in the low-quality or N prefix of a read therefore count, and every quality-split segment of a read rescans the read
   from 0. Only kmers repeated *within one sequence* are non-unique.
7. `findStart` uses `i < stop - k`, so the last kmer is never a threading start. `isThreadingStart` = "not non-unique".
8. The reference-source kmer is never a merge target (`getKmerVertex(kmer, allowRefSource=false)`). A read
   containing it creates a duplicate vertex.
9. `increaseCountsInMatchedKmers` walks back up to `k-1` edges from the start vertex while in-degree is 1,
   comparing against the start vertex's own kmer.
10. `MultiSampleEdge` with one pruning sample: the pruning multiplicity is `max` over per-sample flushes.
    The reference is its own pseudo-sample (`XXX_UNNAMED_XXX`) flushed before the reads, so a ref edge has
    pruning multiplicity `max(1, read support)` while its multiplicity is `1 + read support`. Edges created
    by dangling-end recovery keep pruning multiplicity 1.
11. `AbstractReadThreadingGraph.removeVertex` also removes the vertex's **kmer** from `kmerToVertexMap`,
    even when the map points to another vertex with that kmer. `isLowQualityGraph` (`nonUnique*4 > map.size()`)
    is evaluated *after* chain pruning, so pruned orphan vertices shrink the map.
    Its `removeSingletonOrphanVertices` has no ref-source exception (the `BaseGraph` one does).
12. Chain pruning only visits chains reachable from sources. Chains inside a source-less cycle are never pruned.
13. Dangling tails: `findPath` with `giveUpAtBranch`; the path is cleared, not aborted, when an edge has
    pruning multiplicity < 2. The reference path is followed down with the heaviest incoming edge of
    `altPath[1]` blacklisted. Then `removeTrailingDeletions`, `cigarIsOkayToMerge(≤3 elements, last M)`,
    `longestSuffixMatch`, the leading-deletion `+1` case and the `refIndexToMerge==0` guard.
14. Dangling heads (legacy): `bestPrefixMatchLegacy` returns the index of the **last mismatch** (or -1
    if there is none, which means *no* merge). Max mismatches is `max(1, firstM/k)`. If needed,
    `extendDanglingPathAgainstReference` grows the path with new untracked vertices, and the old source stays
    behind as an orphan. `getBasesForPath(expandSource=true)` reverses the full kmer of source vertices.
15. `toSequenceGraph`: a source vertex contributes its whole kmer and every other vertex its last base.

Seq graph and haplotypes

16. `cleanNonRefPaths` runs before `cleanupSeqGraph`. `removeVerticesNotConnectedToRefRegardlessOfEdgeDirection`
    removes **all** vertices when there is no ref source.
17. `simplifyGraph`: zip once, then rounds of diamonds → tails (min common suffix 10) → split suffixes →
    merge suffixes → zip. After round 5 it compares against a snapshot with `graphEquals` (sequence-based
    vertex and edge sets).
18. A graph that collapses to one vertex gets an empty dummy sink with a ref edge of multiplicity 0.
    Its k-best score is `log10(0) - log10(0) = NaN`; this haplotype is the reference, which is removed anyway.
19. **`KBestHaplotype.isReference()` is always false** (`isReference &= edge.isRef()` on a field that is
    initialized false). Consequences:
    * `refHaplotype.setScore` is never called, so the reference keeps score NaN and kmer 0;
    * the path with the reference bases is added as a *non-ref* haplotype, then removed at the end of
      each graph as `tmpRefHaplotype`. The real reference haplotype is appended after **the first graph
      that produced anything**, and it keeps that position. Final order:
      `[graph-1 haplotypes…, REF, new graph-2 haplotypes…]`;
    * `Haplotype.equals` includes `isRef`, so in later graphs the reference-bases path is not deduplicated
      against REF: it is re-aligned, added and removed again.
20. Deduplication across kmer sizes is by bases (and isRef), first one wins. A haplotype's `kmer` is
    therefore the smallest kmer size that produced it, and the score comes from that graph.
21. `findBestPaths` filters: SW failure (`calculateCigar` returns null when the offset is > 0 or there is a soft clip),
    reference length < 30, N operators, and reference-length mismatch, where GATK throws unless the `INDEL` strategy
    gives the right reference length (then the haplotype is silently dropped).
22. `Haplotype.setCigar` re-normalizes the cigar with `new CigarBuilder(false)` (merges ops and puts D before I).
    All non-ref haplotypes get `alignmentStartHapwrtRef = refHaplotype's` (padded_start - ref_padded_start)
    and genome location = padded span.
23. `CigarUtils.calculateCigar`: identical bases give `NM`. Otherwise both sequences are padded with 10 Ns,
    aligned with the SOFTCLIP strategy, then `trimCigarByBases` (whose `CigarBuilder` removes leading and
    trailing deletions and counts them), `leftAlignIndels`, and the removed leading and trailing deletions
    are re-added.
24. **GKL vs Java SW.** GATK uses the GKL AVX aligner by default. Reading its source (`PairWiseSW.h`)
    shows Gotoh with the same tie rules as `SmithWatermanJavaAligner`: match ≥ insertion ≥ deletion,
    gap extension preferred on ties, and the same end-cell selection for each overhang strategy. The probe
    found 0 differences between them in 120k alignments.

## Implementation notes (speed)

* Arena graphs with `u32` ids and ordered `Vec` adjacency lists; vertex k-mers are offsets into one byte
  arena (reference + reads), so there is no allocation per k-mer.
* Non-unique k-mer detection uses 2-bit packed k-mers (`u64` for k ≤ 32, `u128` for k ≤ 64), with a byte-slice
  table for k > 64 or for k-mers containing non-ACGT bytes. Sets clear in O(1) with generation stamps.
* Smith-Waterman (haplotype→reference and read→haplotype) fills the Java recurrence by anti-diagonals, as GKL's
  AVX aligner does, with one traceback byte per cell (move + gap-extension flags); the inner loop is plain Rust
  that the compiler vectorises. The haplotype aligner (`FixedRefAligner::align_batch`) also skips the columns of
  the prefix a haplotype shares with the previous one in the sorted batch: traceback cells keep fixed positions
  for a fixed reference, and the fill resumes from a saved column (a stack of saved columns, chosen by looking
  ahead in the batch); about 42% of the cells are computed. Results are bit-identical to the row-major Java port:
  `swbench` replays an `VARIANT_RS_SW_DUMP` trace (1.07 M alignments of HG002 shard023) and `swbench --random N` all four
  overhang strategies plus random prefix-sharing batches against `src/bin/swbench/old_sw.rs`. 2026-11: 6x
  (haplotypes) and 8x (reads) faster than the scalar port. (The figures below predate it.)
* GATK's own `runLocalAssembly` on the same 6,345 inputs (`RunAsm.java`, a development probe not in the repository; one JVM, JIT warm-up
  included, GKL AVX-512 SW) takes 15.0 s. The port takes about 4.6 s on one thread.
* Phase profile on the full dump: SW for haplotype cigars about 2.7 s, graph construction about 1.1 s,
  seq graph about 0.4 s, pruning about 0.16 s, everything else < 0.1 s.

## Testing

The Java probes and scripts below were development-only and are not in the repository.

* `RunAsm.java` runs GATK's real `runLocalAssembly` (default `HaplotypeCallerArgumentCollection`,
  FASTEST_AVAILABLE aligner = GKL AVX-512) on any JSONL in the dump schema and writes the same schema.
  Re-running it on `$VARCALL_DATA/asm/gatk_assembly_small.jsonl` reproduces the dump exactly.
* `mkstress.py N seed out [extreme]` generates synthetic regions with tandem repeats, duplicated
  segments, homopolymers, large indels near the ends, paralogous reads, 0–8% errors, low qualities, Ns,
  up to 500x coverage and tiny regions. Runs: 2,000 regions (seed 1) plus 1,500 "extreme" regions (seed 2).
  All 3,499 regions where GATK succeeds are identical, including order, kmer and score. In the one
  region where GATK throws, the port fails the same way.
  These exercise k up to 85, cycle rejection, low-complexity rejection, the FAILED status, head path
  extension and dropped duplicate edges.
* `SwProbe.java` + `swcheck`: 20,000 random pairs × 3 strategies × 2 parameter sets against GKL,
  plus `calculateCigar` (SOFTCLIP and INDEL), with 291 null results. 0 mismatches.
* `GraphProbe.java` + `graphcheck`: 20,000 random seq graphs (bubbles, tails, heads, cycles) through
  cleanNonRefPaths → zip → … → simplify → removePathsNotConnectedToRef → simplify → k-best. The final
  vertex order, edge order, multiplicities and k-best bases and scores all match. This covers MergeTails, k-best cycle
  removal, the `graphEquals` loop guard and JUST_ASSEMBLED_REFERENCE.

Paths not reached by any test: the `findPath` loop guard in dangling recovery and the INDEL-fallback
drop in `findBestPaths`. Both are direct transcriptions.

The large generated inputs and outputs were deleted to save space. They were generated with the development probes: `python3 mkstress.py 2000 1 s.jsonl && java -cp $GATK_JAR RunAsm.java s.jsonl s_gatk.jsonl && target/release/validate s_gatk.jsonl`.
