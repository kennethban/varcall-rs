# hc-preproc

A Rust port of GATK 4.7.0.0 HaplotypeCaller's per-region read preprocessing,
`AssemblyBasedCallerUtils.finalizeRegion(...)`, as called from `assembleReads(...)` with
HaplotypeCaller's default arguments. It takes the reads GATK puts into an active region and
returns the reads GATK hands to the local assembler, in the same order.

## API

```rust
pub type CigarElem = (u8, u32);              // (b'M', 150)
pub struct Read {
    pub name: String, pub flag: u16,
    pub start: i64, pub end: i64,            // 1-based inclusive; `end` is recomputed from start+cigar
    pub cigar: Vec<CigarElem>,
    pub seq: Vec<u8>, pub qual: Vec<u8>,     // qual = raw phred (not +33)
    pub mapq: u8, pub mate_start: i64,       // mate_start = 0 if unpaired / mate unmapped
    pub tlen: i64,
    pub ref_id: i32, pub mate_ref_id: i32,   // BAM reference index of the read and of its mate (RNEXT)
}
pub fn finalize_region(raw_reads: &[Read], padded_start: i64, padded_end: i64) -> Vec<Read>;   // clones the input
pub fn finalize_region_owned(reads: Vec<Read>, padded_start: i64, padded_end: i64) -> Vec<Read>; // no clones, in place
pub fn finalize_region_with(..., p: &FinalizeParams) / finalize_region_owned_with(..., p: &FinalizeParams)
pub struct FinalizeParams { min_tail_quality: u8 /*9*/, dont_use_soft_clipped_bases: bool /*false*/,
                            override_softclip_fragment_check: bool /*false*/, correct_overlapping_base_qualities: bool /*true*/ }
pub fn has_well_defined_fragment_size(&Read) -> bool;  pub fn parse_cigar(&str); pub fn cigar_to_string(&[CigarElem]);
```

`raw_reads` must be in `region.getReads()` order, which is the order GATK added them in.
The library has no dependencies. serde/serde_json are only pulled in by the `validate` feature
(on by default) for the validator binary.

## Validation

```
RUSTFLAGS="-C target-cpu=native" cargo build --release -p hc-preproc    # from the repository root
./target/release/preproc-validate $VARCALL_DATA/asm/gatk_assembly.jsonl [n_details]
```
On the full dump (6,345 regions, 586,027 raw reads, 586,022 expected final reads), every region
and every read matches exactly: same reads, same order, and the same start, end, cigar, seq and
qual. The run takes about 0.87 s in total: about 0.38 s parsing JSON and about 0.054 s in
`finalize_region` (roughly 90 ns per read).

`fuzz/Fuzz.java` (a development probe, not in the repository) is a differential fuzzer. It builds random regions and reads (soft and hard
clips, I/D next to each other, low-quality tails, odd TLENs, overlapping mates, reads near
contig position 1), runs the real GATK `finalizeRegion` on them, and writes JSONL in the dump
format:
`java -cp $GATK_JAR fuzz/Fuzz.java <seed> <nRegions> out.jsonl`.
Seeds 1 to 5 (23,000 regions, about 497k final reads) all match 100%. The fuzzer covers paths
that the real data never reaches: clipping both ends to the region, reads that are entirely
low quality, reverted soft clips that run past position 1, and the CigarBuilder D+I+clip quirk.

## What was ported

* `finalizeRegion` main loop and `HardClipAndPossiblyAddToCollection`. The `trackHardclippedReads`
  branch is skipped because it does not affect the returned reads.
* `ReadClipper`: `hardClipSoftClippedBases`, `revertSoftClippedBases`, `hardClipLowQualEnds`,
  `hardClipAdaptorSequence`, `hardClipToRegion`, `hardClipBothEndsByReferenceCoordinates`,
  `clipByReferenceCoordinates`, `clipRead`.
* `ClippingOp.applyHardClipBases` and `applyRevertSoftClippedBases`.
* `CigarUtils.clipCigar`, `alignmentStartShift`, `revertSoftClips`, and `CigarBuilder` in full.
* `ReadUtils`: `hasWellDefinedFragmentSize`, `getAdaptorBoundary`, `getSoftStart`,
  `getReadIndexForReferenceCoordinate`, `isInsideRead`.
* `ReadCoordinateComparator`, with a stable sort.
* `cleanOverlappingReadPairs`, `FragmentCollection.create` and
  `FragmentUtils.adjustQualsOfOverlappingPairedFragments`.
* htsjdk: `SAMRecord.getAlignmentEnd` (start + reference length - 1), `Locatable.overlaps` /
  `CoordMath.overlaps`, and `Cigar.getReadLength`.

## Defaults checked in the argument collections

`dontUseSoftClippedBases=false`, `overrideSoftclipFragmentCheck=false`,
`softClipLowQualityEnds=false` (so the code path is `hardClipLowQualEnds`),
`errorCorrectReads=false`, `minBaseQualityScore=10` (so minTailQuality is 9),
`doNotCorrectOverlappingBaseQualities=false` (so overlap correction is on), and
`pileupErrorCorrectionLogOdds=-inf` (no read error corrector, so the assembler leaves the reads
unchanged before the dump).

## GATK behaviours replicated

1. **Revert or hard clip soft clips.** The choice uses `hasWellDefinedFragmentSize` on the
   *original* read. It returns true when TLEN != 0, the read is paired, both mates are mapped,
   the strands differ, and then:
   - reverse read: `end > mateStart`;
   - forward read: `start <= mateStart + TLEN`. TLEN is signed here; no `abs` is taken.
2. **Reverting soft clips** happens only if the first *or* last cigar element is a clip. That
   includes `H`, so a read with only hard clips still has its cigar rebuilt by CigarBuilder.
   S becomes M and merges with neighbouring M elements, and `start = softStart`. If
   softStart <= 0, the read is set to start 1, the first `1 - softStart` bases are hard-clipped,
   and the start is set back to 1.
3. **Low-quality tails.** Bases with qual `<= 9` are clipped. The right tail is clipped first,
   then the left. If every base is low quality, the read is dropped.
4. **Adaptor boundary** is computed on the *already reverted and clipped* read, and
   `hasWellDefinedFragmentSize` is evaluated again on that read:
   - reverse read: `mateStart - 1`;
   - forward read: `start + |TLEN|`. This uses the current start, which reverting soft clips may
     have moved.

   Clipping happens only if the boundary lies within [start, end]. A reverse read is clipped on
   the left tail, a forward read on the right tail.
5. **`getReadIndexForReferenceCoordinate`** starts at the soft start and counts S as consuming
   reference. If the coordinate falls in a deletion, it returns the index *after* the deletion.
   The left-tail clip then subtracts 1, so it stops before the deletion. A right-tail clip
   starts at the base after the deletion.
6. **`clipCigar`** drops deletions that sit exactly on the clip edge. CigarBuilder then removes
   leading deletions (including D right after a leading I), trailing deletions, and a D sitting
   before a trailing I. It merges equal neighbouring operators and moves D to the left of I.
   - **Stale `lastOperator` quirk:** `[.., D, I] + clip` becomes `[.., I, clip]` without
     updating `lastOperator`, which stays `I`. A trailing `H` added next is therefore *not*
     merged with the new H. For example, `50M2D3I45M5H` clipped from read index 53 becomes
     `50M3I45H5H`. This is replicated.
7. **`alignmentStartShift`** counts reference bases of elements that lie entirely within the
   clipped prefix, including deletions immediately after it, plus the clipped part of the
   element the clip ends in.
8. **`applyHardClipBases`** for a right-tail clip keeps the first `newLength` bases. When the
   clip starts at 0, it shifts the start by `alignmentStartShift`.
9. **Read checks** in `HardClipAndPossiblyAddToCollection`, in order:
   - `start <= end` and the read is mapped;
   - after adaptor clipping: the read is not empty and the cigar read length is > 0;
   - after region clipping: `start <= end`, length > 0, and the *adaptor-clipped* (not
     region-clipped) read overlaps the padded span. The overlap test uses htsjdk
     `CoordMath.overlaps`: `(s2>=s1&&s2<=e1)||(e2>=s1&&e2<=e1)||(s1>=s2&&e1<=e2)`.
10. **`hardClipToRegion`** compares against the padded span:
    - a read sticking out on both sides is clipped at `paddedStart-1` and `paddedEnd+1`, right
      tail first;
    - if, after the right clip, `paddedStart-1 > end`, the read is dropped;
    - a read entirely outside the padded span is dropped.
11. **Empty reads** (GATK `emptyRead`, unmapped) can come out of any step and are always
    dropped.
12. **Sort order** (`ReadCoordinateComparator`): start, then forward before reverse, then name
    (String.compareTo), SAM flags, MAPQ, (if both reads are paired) mate reference index and
    mate start, and finally TLEN. Java's `List.sort` is stable, so `slice::sort_by` is used.
    * The flags are `ReadUtils.getSAMFlagsForRead`, which GATK *rebuilds* from the read: 0x8 is set
      whenever `mateIsUnmapped()` (also for a mate without position), 0x20 only for a mapped mate,
      0x10 only for a mapped read, 0x2/0x40/0x80 only for paired reads (`Read::gatk_sam_flags`).
    * The mate reference index is `ReadUtils.getMateReferenceIndex`: -1 for an unmapped mate, else
      the mate contig's index (`mate_ref_id`).
13. **Overlapping mates** (`FragmentCollection`). A read is a singleton if it is unpaired, its
    mate is unmapped, `mateStart == 0`, or `mateStart > end`. Here `mateStart` is the mate's
    *original* position, while `end` is the read's end *after* clipping. Otherwise the read is
    paired with an earlier read of the same name from a name map.
14. **Quality adjustment** (`adjustQualsOfOverlappingPairedFragments`):
    - "first" is the read with the strictly smaller soft start. On a tie, the *later* read in
      sort order is "first".
    - No adjustment if `first.end < second.start`, or if the second read's start falls in a
      clip of the first read.
    - Overlap length = `min(firstEndIdx - offset, secondEndIdx - secondOffset) + 1`. It is
      walked base by base with no indel awareness.
    - Equal bases (N == N counts as equal): both quals are capped at 20, which is
      `round(-10*log10(1e-4)) / 2`.
    - Different bases: both quals are set to 0 (`setConflictingToZero=true`).
    - Indel quals are not touched (`halfOfPcrIndelQual` is empty).
15. **`end`** is always recomputed as `start + refLen(cigar) - 1`. A read whose cigar has no
    reference-consuming operators gets `end = start - 1`.

## Not ported / assumptions

* `softClipLowQualityEnds=true` (the SOFTCLIP_BASES representation, DRAGEN BQD mode),
  `errorCorrectReads`, and the `trackHardclippedReads` (pileup-detection) output are not ported.
* A single sample is assumed. GATK splits reads by sample before pairing mates, but with one
  sample the result is the same.
* Mates on another contig: `Read` carries `ref_id` and `mate_ref_id` (the validator's dumps have no
  mate contig and use the default 0 / 0, i.e. "same contig"). GATK's `hasWellDefinedFragmentSize`,
  `getAdaptorBoundary` and `FragmentCollection` have **no** contig test: they use TLEN (0 for such pairs
  in BWA output) and the mate *position* even when it is on another contig, and so does the port.
  The mate contig only enters the comparator (above) and `Read::mate_on_other_contig`
  (`MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE`, used downstream by `filterNonPassingReads`).
* Inputs on which GATK would throw (an NPE or GATKException in `clipByReferenceCoordinates`,
  CigarBuilder validation errors, or `setPosition(<1)`) are handled leniently, by leaving the
  read unclipped, instead of aborting.
