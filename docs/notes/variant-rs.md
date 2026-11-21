# variant-rs development notes

The detailed record of variant-rs: every GATK behaviour reproduced, how each tool was validated (with the data
and the record counts), and performance measurements, mostly from November 2026. Paths are relative to
`variant-rs/` unless they start with `tests/harness/` (the harness scripts) or `$VARCALL_DATA` (their data, see
`tests/harness/README.md`). For the current picture see [reference.md](../reference.md),
[architecture.md](../architecture.md) and [validation.md](../validation.md).

`variant-rs` runs HaplotypeCaller, BaseRecalibrator, ApplyBQSR, MarkDuplicates and the Mutect2 tools with outputs
byte-identical to GATK 4.7. It takes GATK's own command lines (`variant-rs HaplotypeCaller -R ... -I ... -O ...`, see
"GATK command lines" below) as well as its original options. As a HaplotypeCaller it is the Rust port of GATK 4.7 HaplotypeCaller (default arguments; one sample
in VCF or `-ERC GVCF` mode, several samples called jointly in VCF mode). It reads a coordinate-sorted, indexed BAM or CRAM and writes bgzipped, tabix-indexed VCF/GVCF whose
data lines are byte-identical to GATK's. The stages after region detection live in the crates of `variant-rs/crates`
(preproc, assembly, likelihood, cudahmm, genotype, refconf, downsample).

```
per contig (one GATK read shard), several contigs at a time:
  BAI chunks -> BGZF inflate (reader thread) -> htsjdk interval filter + HC read filters (decode thread)
    -> positional downsampling (hc_downsample) -> pileup scatter / isActive / activity profile / regions
    -> region filling -> batches of regions  ==> back end thread:
       finalizeRegion, assembly, trimming, PairHMM (one CUDA launch per batch, GPU shared by all contigs),
       realignment, genotyping (in order) -> VCF text   |   + reference confidence -> GVCF banding (writer thread)
  -> per-contig part -> main thread: output in dictionary order, bgzip blocks + tabix index
```

## GATK command lines (`src/gatk_cli.rs`)

```bash
variant-rs HaplotypeCaller -R ref.fa -I x.bam -L chr20 -O x.g.vcf.gz -ERC GVCF      # as `gatk HaplotypeCaller ...`
variant-rs BaseRecalibrator -R ref.fa -I x.bam --known-sites dbsnp.vcf.gz -O recal.table
variant-rs ApplyBQSR -R ref.fa -I x.bam --bqsr-recal-file recal.table -O x.bqsr.bam
variant-rs MarkDuplicates -I x.bam -O x.md.bam -M metrics.txt [--TMP_DIR t]          # or Picard's I=x.bam O=... M=...
variant-rs GatherBQSRReports -I scatter1.table -I scatter2.table -O recal.table
variant-rs GatherBamFiles INPUT=scatter1.bam INPUT=scatter2.bam OUTPUT=x.bam CREATE_INDEX=true CREATE_MD5_FILE=true
variant-rs MergeVcfs INPUT=scatter1.g.vcf.gz INPUT=scatter2.g.vcf.gz OUTPUT=x.g.vcf.gz
variant-rs ReblockGVCF -R ref.fa -V x.g.vcf.gz -do-qual-approx --floor-blocks -GQB 20 -GQB 30 -GQB 40 -O x.rb.g.vcf.gz
variant-rs GenotypeGVCFs -R ref.fa -V x.g.vcf.gz [-D dbsnp.vcf.gz -G StandardAnnotation -G AS_StandardAnnotation -L chr8 \
        --only-output-calls-starting-in-intervals --merge-input-intervals] -O x.vcf.gz     # a GVCF file or gendb://ws
variant-rs SortSam INPUT=x.md.bam OUTPUT=x.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true CREATE_MD5_FILE=true
variant-rs MarkDuplicates INPUT=a.bam INPUT=b.bam OUTPUT=x.md.bam METRICS_FILE=m.txt ASSUME_SORT_ORDER=queryname \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CLEAR_DT=false ADD_PG_TAG_TO_READS=false   # WARP's command line
variant-rs Mutect2 -R ref.fa -I tumor.bam [-I normal.bam -normal N] --germline-resource g.vcf.gz [-pon pon.vcf.gz] \
        [--f1r2-tar-gz f1r2.tar.gz] [--genotype-germline-sites] [--genotype-pon-sites] [--bam-output bamout.bam] -O m2.vcf
variant-rs Mutect2 -R ref.fa -I m.bam -L chrM --mitochondria-mode [--read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --read-filter MateUnmappedAndUnmappedReadFilter --annotation StrandBiasBySample --max-mnp-distance 0] -O m2.vcf
variant-rs LearnReadOrientationModel -I f1r2.tar.gz -O priors.tar.gz
variant-rs GetPileupSummaries -I tumor.bam -V common.vcf.gz -L common.vcf.gz -O tumor.pileups.table   # or -L piece.bed (sarek)
variant-rs GatherPileupSummaries --I p1.pileups.table --I p2.pileups.table --O tumor.pileups.table --sequence-dictionary ref.dict
variant-rs MergeMutectStats --stats p1.vcf.gz.stats --stats p2.vcf.gz.stats -O m2.vcf.gz.stats   # scattered Mutect2 (sarek)
variant-rs CalculateContamination -I tumor.pileups.table [-matched normal.pileups.table] -O contamination.table [-segments segments.table]
variant-rs FilterMutectCalls -R ref.fa -V m2.vcf [--contamination-table ...] [--tumor-segmentation ...] [--ob-priors priors.tar.gz] -O filtered.vcf
variant-rs FilterMutectCalls -R ref.fa -V m2.vcf --stats m2.vcf.stats --mitochondria-mode [--max-alt-allele-count 4 \
        --min-allele-fraction 0.03 --f-score-beta 1.0 --contamination-estimate 0.02] -O filtered.vcf   # AlignAndCall.wdl
variant-rs <tool> --help                                       # what is mapped, ignored, and accepted at default only
GATK_FALLBACK=/path/to/gatk variant-rs <any tool> ...           # anything variant-rs refuses runs the real GATK, unchanged
```

A first argument that is a GATK tool name (CamelCase) selects this mode; lowercase subcommands and the bare
HaplotypeCaller options below keep working. GATK's grammar is accepted: `--name value`, `--name=value`, short
names (with either `-` or `--`, as GATK allows), repeated list arguments, booleans with or without a value,
`--java-options "..."` (ignored). Every GATK
argument of the tool (tables in `gatkargs/`, generated from `gatk <Tool> --help` of GATK 4.7 by `gatkargs/gen.py`)
is either mapped to the tool's option, ignored because it cannot change the output (logging, temporary
directories, JDK codecs, cloud buffers, index creation, `--native-pair-hmm-threads`, `--smith-waterman`, Mutect2's
deprecated and unused `--tumor-sample`), or
accepted only at GATK's default value (`-ploidy 2`, `--kmer-size 10 --kmer-size 25`, ...); anything else,
including a second `-I` where variant-rs takes one, VCF interval files (GetPileupSummaries instead also takes a
sites VCF as `-L`), and unknown names, is refused with the reason. The thread count defaults to all CPUs (outputs
do not depend on it). Without `GATK_FALLBACK` a refused command exits with status
2. Checked by `tests/gatk_cli_check.sh` (GATK's command lines from `tests/harness/hg002/run_gatk.sh`, `tests/harness/bqsr/truth.sh`,
`tests/harness/markdup/truth.sh`, `tests/harness/somatic/orient/truth.sh`, `tests/harness/somatic/contam/truth.sh`: outputs identical to GATK's;
fallback with a stub) and `cargo test gatk_cli`.

## Usage

```bash
variant-rs -I NA12878.bam -R GRCh38_full_analysis_set_plus_decoy_hla.fa \
       [-L chr20:10000000-15000000 | -L file.intervals]...   # repeatable; none = every contig of the BAM
       -O out.vcf.gz [-ERC GVCF] --threads 16 [--pairhmm cuda|cpu|auto]

# GATK-per-shard mode: one RNG/state domain per shard file (= one GATK JVM per shard)
variant-rs -I NA12878.bam -R ref.fa --shard shards/shard00.intervals ... --shard shards/shard22.intervals \
       -O all.vcf.gz [--shard-output-dir per_shard/] [-ERC GVCF] --threads 16

# shard planner: cut the genome into ~96 load-balanced shard files (contigs split where they are large)
variant-rs plan -I NA12878.bam -R ref.fa -o plan96 [--pieces 96] [--min-gap 1000] [--base-cost 0.05]
variant-rs -I NA12878.bam -R ref.fa $(for f in plan96/shard*.intervals; do echo --shard $f; done) -O all.vcf.gz

# base quality score recalibration (GATK BaseRecalibrator + ApplyBQSR, default arguments)
variant-rs bqsr -I in.bam -R ref.fa --known-sites dbsnp138.vcf.gz --known-sites Mills.vcf.gz \
       --known-sites known_indels.vcf.gz -O recal.table [-L ...] [--threads 16]
variant-rs applybqsr -I in.bam -R ref.fa --bqsr-recal-file recal.table -O out.bam [-L ...] [--threads 16]
# duplicate marking (Picard MarkDuplicates, default arguments; coordinate-sorted input)
variant-rs markdup -I sorted.bam -O dedup.bam -M dup_metrics.txt [--tmp-dir DIR] [--threads 16]
# hard filters (GATK Best Practices: SelectVariants SNP/INDEL -> VariantFiltration -> MergeVcfs, in one pass)
variant-rs filter -V calls.vcf.gz -O filtered.vcf.gz [--keep-other]
# or recalibrate on the fly inside HaplotypeCaller (no recalibrated BAM written)
variant-rs -I in.bam -R ref.fa --bqsr-recal-file recal.table -O out.vcf.gz [-ERC GVCF] ...
```

| Option | Meaning |
|---|---|
| `-I/--input` | BAM (`X.bam.bai` or `X.bai` must exist) or CRAM (`X.cram.crai` or `X.crai`; decoded with `-R`), detected by content |
| `-R/--reference` | FASTA with `.fai` (the `.dict`, if present, gives the dictionary order) |
| `-L/--intervals` | `contig`, `contig:start-end`, `contig:start`, `contig:start+`, or a file: `.intervals`/`.list` (one such string per line), `.bed`, `.interval_list`. Contig names containing `:` (HLA alleles) are resolved against the dictionary, as GATK does. Intervals are sorted in dictionary order and overlapping or abutting ones merged (`--interval-merging-rule ALL`). |
| `--shard FILE` | repeatable; per-shard mode (below). Not combinable with `-L`. |
| `-O/--output` | `.vcf.gz`: bgzip + `.tbi`; anything else: plain text. (A `-O x.igv` is still accepted as the old region output.) |
| `-ERC GVCF` | reference-confidence mode (`--emit-ref-confidence GVCF`); `--gvcf FILE` is the old spelling |
| `--shard-output-dir DIR` | per-shard mode: also write `DIR/<shard name>.vcf.gz` (or `.g.vcf.gz`) + `.tbi` per shard file |
| `--threads N` | rayon threads (assembly, likelihoods, realignment, reference confidence, deflate/inflate) |
| `--contig-workers N` | contigs processed at once (default `threads/2`, at most 8) |
| `--max-workers-per-shard N` | at most N workers on one shard at a time (default 2; not applied with a single domain) |
| `--max-reads-per-alignment-start N` | positional downsampling cap (GATK default 50; 0 = off) |
| `--pairhmm cuda/cpu/auto` | CUDA PairHMM (default), the bit-exact CPU port of GKL's AVX-512 PairHMM, or CUDA if available |
| `--compression-level N` | BGZF level (libdeflate, default 5) |
| `--assembly-region-out F` | GATK's `--assembly-region-out` IGV file |
| `--gvcf-raw F` | GVCF mode: the unbanded per-position records (GATK `GATK_RAW_CALLS_VCF` debug output), plain text |
| `--tmp-dir D` | where GVCF parts are written (default: next to the output) |
| `--vcf-header F` | copy this header verbatim instead of generating it |
| `--dump-activity F` | per-locus isActive dump (single contig only) |

## Two modes: which GATK run they equal

GATK's only state that crosses contigs is the JVM-wide `java.util.Random` (`Utils.getRandomGenerator()`,
seed 47382911). Two things draw from it, in traversal order: positional downsampling (`nextInt`, once per read
beyond 50 starting at one position) and the QD jitter (`nextGaussian`, every call with QD ≥ 35; VCF mode only).
Everything else is per contig: GATK makes one read shard per contig, with a fresh activity profile and region
iterator; the spanning-deletion list is effectively per contig (it checks the contig); the GVCF block combiner
flushes at contig changes.

* **Single-run mode** (`-L ...` or no `-L`): equals **one GATK JVM over all the intervals**. One generator for
  the whole run, handed from contig to contig in dictionary order.
* **Per-shard mode** (`--shard f1 --shard f2 ...`): equals **one GATK JVM per shard file**: a fresh generator at
  the start of each shard, carried across that shard's contigs in order. `-O` is the concatenation of the shards in
  the order given (one header); `--shard-output-dir` additionally writes each shard on its own, identical to that
  shard's GATK output. This is the mode that matches `wgs/truth/` (one JVM per `shards/shardNN.intervals`).
* **Per-interval mode** (`-L file --scatter-intervals`): every `-L` interval is a shard of its own, = one GATK JVM
  per interval, as nf-core/sarek scatters HaplotypeCaller. On NA12878 30x with sarek's 356 calling regions:
  4,954,925 records identical to sarek's stock VCF; 407-419 s against 672 s in single-run mode (2026-11-16).

How the draws keep GATK's order while contigs run in parallel (`src/rng.rs`, `src/frontend.rs`):
* The genotyping engine does not draw the QD jitter; it writes a placeholder (`EngineConfig::defer_qd_jitter`).
* The downsampler (`hc_downsample`) needs the true generator state at a *drawing group* (a start position with
  more than 50 reads). Following `variant-rs/crates/downsample/README.md` §3–4, the front end computes each region's
  read-pull horizon X_i (readiness tracking, lookahead locus, first covered locus, P1/P2). Before drawing group P,
  the back end genotypes every region with X_i < P and draws their QD placeholders; then P draws.
* The generator is a baton per domain: contig k gets it when contig k−1 of the same domain is done with it.
  A contig that never needs it (no drawing group) has its QD placeholders drawn by the main thread when its turn
  comes; a contig without any draw is skipped. So contigs without drawing groups never wait; a contig with
  drawing groups waits at its first one until the previous contigs of its domain are finished.
* GVCF mode makes no QD draws, so only the downsampler uses the baton.

Scheduling: longest estimated time first (LPT). A domain's estimate is its BAI cost (below) divided by the number
of workers it can use at once (its contig count, capped by `--max-workers-per-shard`, default 2), so long baton
chains of small contigs (shard22, shard21) still start early. At most `--max-workers-per-shard` workers per domain.
A contig may be split across consecutive shard files (as `variant-rs plan` does); shards must be in dictionary order
and disjoint so their outputs concatenate into a sorted file.

## Shard planner (`variant-rs plan`, `src/plan.rs`)

Writes `OUT/shardNNN.intervals` (one `contig:start-end` per line) and `OUT/plan.tsv` (shard, contig count, first and
last interval, bases, estimated cost fraction). Each shard file is one GATK JVM's worth of work, so `--shard` output
equals GATK run with `-L shardNNN.intervals` on the same files.

* **Cost model, from the BAI only** (no BAM read): the linear index holds, per 16 kb window, the smallest file
  offset of the reads overlapping it, so compressed bytes between consecutive windows estimate the read data per
  window. Plus `--base-cost` (default 0.05 of the mean per-base read cost) so read-free sequence is not free.
* **Partition:** contigs in dictionary order; target = total cost / `--pieces`. A contig worth ≥ 2 targets is split
  into that many near-equal-cost pieces. Each cut moves to the middle of an N run ≥ `--min-gap` (1 kb) within ±1/8
  of a piece length of the ideal cut (no read or region crosses an N gap there); otherwise to the 16 kb window with
  the fewest reads within ±1/32 of a piece. Smaller contigs are grouped, in order, into shards of ≥ 0.9 × target.
  A gapless (T2T) reference just gets low-coverage cuts.
* On NA12878 / GRCh38 + decoys + HLA (`--pieces 96`): 97 shards, 22 contigs split (chr1 into 8), 15 cuts in N
  gaps, largest shard 1.37% of the estimated total.
* **Cuts change nothing but RNG domains**, exactly as in GATK: with 97 instead of 23 shards, GATK's own output
  differs in ~317 k records, 311,774 of them only QD (the jitter generator restarts per shard) and 5,731 at very
  high depth (downsampling draws); none within 2 kb of a cut differ in anything but QD.
* On one 16-core machine it does not make the run faster (see Performance): the whole genome is throughput-bound,
  not bound by chr1. Its use is many workers (several GPUs or nodes), where 23 chromosome shards are too coarse.

## CRAM input (`src/cram.rs`, `src/input.rs`, `csrc/cram_shim.c`)

CRAM is decoded by htslib (the conda env's `libhts`, 1.25), the same library that converts CRAM to BAM
(`samtools view -b`). A small C shim runs htslib's multi-interval iterator (`sam_itr_regions`) over the padded
intervals of a contig and serializes every record exactly as `bam_write1` would (BAM record bytes, htslib's
query-name padding removed), so the front end and its htsjdk interval filter are unchanged. One reader thread per
contig decodes; open handles (each with its own index, iterator and reference cache) are pooled across contigs.
The scheduler and `variant-rs plan` take their cost model from the `.crai` (compressed slice bytes over each slice's
span) instead of the BAI.

* Do **not** share htslib's reference cache between handles (`CRAM_OPT_SHARED_REF`): it keeps only the last-used
  contig, so workers on different contigs evict each other and reload whole chromosomes for almost every slice
  (5.4 TB read from the page cache, ~2.6 cores busy, in a first whole-genome attempt). Each handle keeps its own
  current contig (≤ ~250 MB).
* **GATK 4.7 cannot open the 1000 Genomes CRAM with its own `.crai`**: htsjdk converts the CRAI into a BAI in memory
  and rejects the 11 unmapped-slice entries (`-1 0 150 ...`; it expects span 0): `CRAMException: Attempt to create a
  bai entry for an unmapped slice with unexpected alignment start (0) or span (150)`. A copy of the index with those
  spans set to 0 works (`wgs/cramtest/NA12878.cram.crai`). variant-rs (htslib) reads the original index.

## BQSR (`variant-rs bqsr`, `variant-rs applybqsr`; `src/bqsr/`)

GATK 4.7 `BaseRecalibrator` and `ApplyBQSR` with default arguments: standard covariates (read group, reported
quality, context 2, cycle), substitution tables only (GATK4 builds no indel tables by default), no BAQ, 16-level
quantization in the report, no quantization when applying (`--quantize-quals 0`), qualities < 6 preserved.

* **`bqsr`** (`recal.rs`, `covariates.rs`, `known.rs`, `datum.rs`, `report.rs`): BQSR read filters (MQ 0/255,
  unmapped, secondary, duplicate, QC fail, WellformedReadFilter), then the engine's read transform (cigar
  consolidation, `hardClipAdaptorSequence`, `hardClipSoftClippedBases`; from `hc_preproc`), mismatch and
  known-site arrays, covariates. Without BAQ every error is 0 or 1, so a table cell is an integer
  (observations, errors) pair: chunks of 2 Mb run in parallel and their tables are summed. Finalization follows
  GATK exactly: the read-group table is collapsed from the quality table in NestedIntegerArray leaf order
  (`RecalDatum.combine`, floating point), values rounded (`roundToNDecimalPlaces`), the Bayesian empirical quality
  (commons-math `SaddlePointExpansion` from `hc_genotype::javamath`), `QualQuantizer`, and GATKReport v1.1
  formatting (column widths, alignment, SORT_BY_COLUMN row order, Java `%.Nf` rounding).
* **Known sites**: all records of the VCFs (any FILTER), span `POS..POS+len(REF)-1` (or `INFO/END`), loaded with
  parallel BGZF inflate (64 M records of dbSNP 138 + Mills + known_indels in ~4 s, ~1 GB).
* **WARP's arguments** (S3 group 1): `--use-original-qualities` (`-OQ`) in both tools (the OQ tag's qualities
  replace the read's after the read filters, as `ReadUtils.resetOriginalBaseQualities`; the OQ tag stays);
  ApplyBQSR `--static-quantized-quals` / `--round-down-quantized` (`constructStaticQuantizedMapping`, applied after
  the report's quantization as GATK does; checked against GATK's own method in `cargo test`) and
  `--create-output-bam-md5` (`<out>.bam.md5`, 32 hex digits, no newline, as htsjdk). `tests/harness/bqsr/cmp_warp.sh`: WARP's
  command lines on NA12878 chr20:10-20M with OQ tags on half the reads: table and 2,427,612 records identical.
* **`applybqsr`** (`apply.rs`): reads the report (ours or GATK's), precomputes the posterior quality per
  (read group, quality) and the context/cycle deltas (`hierarchicalBayesianQualityEstimate`), then rewrites the
  quality bytes of each record in place. GATK's htsjdk re-encodes every record's tags when writing a modified
  read: tags sorted by binary key (`(c1 << 8) | c0`), integers retyped (smallest of c, C, s, S, i, I, signed
  first), BI/BD cleared. `applybqsr` does the same, so records are byte-identical to GATK's. Reads failing
  WellformedReadFilter are dropped, as in GATK. Output: BGZF (libdeflate level 5, parallel), `.bai` written
  alongside, header = input header + one `@PG` line. Without `-L`, unplaced unmapped reads are copied at the end
  (BAM and, since 2026-11-17, CRAM input: before, a CRAM's unplaced reads were dropped from the output without a
  warning, here and in `markdup`; `tests/harness/markdup/cram/run.sh` checks CRAM against BAM input).

| Test (GATK 4.7 on the same input, same known sites) | Recal table | Recalibrated BAM records |
|---|---|---|
| HG002 PCR+ NovaSeq 30x, `-L chr20` (1 read group) | **byte-identical** (14,556,032 reads) | **16,489,891 byte-identical** |
| **HG002 PCR+ NovaSeq 30x, whole genome** | **byte-identical** (642,394,755 reads); `bqsr` 95 s vs GATK 99 min | **721,883,698 byte-identical**; `applybqsr` 6 min vs GATK 1 h 47 min |
| NA12878 1000G/NYGC 30x, `-L chr20` (16 read groups) | **byte-identical** (14,388,910 reads) | **17,705,654 byte-identical** (after htsjdk tag re-encoding) |

**On the fly** (`variant-rs ... --bqsr-recal-file recal.table`): the HaplotypeCaller decode thread recalibrates each
record's qualities (`apply::Recalibrator`) before parsing it, so the caller sees what it would read from the
ApplyBQSR output. chr20 of HG002 PCR+: VCF (130,089 lines) and GVCF (11,686,333 lines) identical to `variant-rs` on
GATK's recalibrated BAM; run time +1 s.

`.bai`: `samtools idxstats` and region queries identical to GATK's. Times (chr20): `bqsr` 6.6 s vs GATK 130 s;
`applybqsr` 9.6 s vs GATK 185 s.

## MarkDuplicates (`variant-rs markdup`; `src/markdup/`, `src/bamwrite.rs`)

Picard 3.5.0 `MarkDuplicates` (as bundled with GATK 4.7) with default arguments on a coordinate-sorted BAM:
SUM_OF_BASE_QUALITIES scoring, optical duplicates from the default read-name parsing (last three `:` fields,
100 pixels, graph clustering for sets of 4+), no DT tagging, CLEAR_DT, `PG` tags chained to a new `@PG`.

* **Pass 1** (2 Mb chunks in parallel): one fragment end per primary mapped read (library, 5' unclipped
  position, orientation, score, tile/x/y, record index) and, for pairs with both mates mapped, a pair end once
  both mates are seen, joined by read group + name exactly like Picard's read-ends map (filed under the mate's
  contig); mates in other chunks are joined afterwards in file order. Record indexes are `chunk << 32 | local`,
  which orders like the file index. Ends are written to on-disk buckets by (contig, 5' position) (56 bytes each).
* **Grouping** (buckets in parallel): each bucket is sorted with `ReadEndsMDComparator` (library, positions,
  orientation, tile, x, y, indexes) and duplicate sets are marked as in `markDuplicatePairs` /
  `markDuplicateFragments` (best = first maximum score); optical duplicates (`OpticalDuplicateFinder`, including
  the Java HashMap iteration order of its cluster map) and the set-size histograms. Duplicates go into a bitset.
* **Pass 2** (chunks in parallel, written in order): duplicate flag set or cleared, `PG` remapped (Picard's
  `HashSet` order and `PgIdGenerator` base-36 suffixes), DT removed, tags re-encoded the htsjdk way, BGZF + `.bai`.
* **Scheduling** (2026-11-16): chunks whose index-estimated read data is well above a 2 Mb chunk at mean depth are
  cut into equal pieces (on GRCh38 + decoys a 33 kb decoy, chrEBV and chrM each hold ~6x a normal chunk; one thread
  each finished both passes alone); mates across chunks are joined in parallel, partitioned by read name; buckets
  are grouped largest first, with duplicates set straight into the bitset; pass 2 runs as a pipeline (the pool
  keeps compressing up to 4 x threads chunks ahead while one thread appends them in order) instead of rounds of
  128 chunks each followed by a single-threaded write. NA12878 30x (41 GB BAM, 32 threads): 247 s -> 196 s
  (pass 2 153 -> 103 s), mean busy threads 17.7 -> 24.9, time under 16 busy 45% -> 20%; all 768,580,569 records,
  header, idxstats, index queries and metrics identical to the previous version's.
* **Metrics**: `DuplicationMetrics` per library (Lander-Waterman library size by bisection, ROI histogram for
  one library) in htsjdk's MetricsFile format (FormatUtil: 6 fraction digits, HALF_DOWN).

| Test (Picard 3.5.0 via GATK 4.7, same input) | Records | Metrics (from `## METRICS CLASS`) |
|---|---|---|
| NA12878 smoke BAM (bwa-mem2 + samtools sort, 1.19 M reads, no PG tags) | **1,192,285 byte-identical** (128,457 duplicates) | **identical** |
| NA12878 1000G chr20 slice, duplicate flags cleared (16 read groups, existing PG tags, 378,597 reads whose mates are on other contigs) | **17,705,654 byte-identical**; `@PG` IDs identical (`MarkDuplicates.G`) | **identical** |

| **HG002 PCR+ NovaSeq 30x whole genome**, duplicate flags and PG tags cleared | **721,883,698 byte-identical** (37,487,487 duplicates) | **identical** |

**WARP's arguments** (S3 group 1; `MarkDuplicates INPUT=a INPUT=b ... VALIDATION_STRINGENCY=SILENT
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 ASSUME_SORT_ORDER=queryname CLEAR_DT=false ADD_PG_TAG_TO_READS=false`):

* **Query-grouped input** (`src/markdup/qname.rs`; `ASSUME_SORT_ORDER=queryname` or a `SO:queryname` header):
  read ends carry the index of the first record of their name group; the same buckets and duplicate sets as
  above; pass 2 is Picard's walk over the sorted duplicate indexes (`nextIndexIfNeeded`), so every record of a
  duplicate group is flagged (secondary, supplementary and unmapped mates too). Output in input order (records
  decided in one thread, BGZF in parallel), `@HD SO:unknown GO:query` when the order was assumed, no index.
* **Several inputs**: htsjdk's `SamFileHeaderMerger` (identical dictionaries required; @PG trees merged level by
  level and @RG merged, colliding IDs renamed `ID.<base36>`, both sorted by ID; comments concatenated) and
  `MergingSamRecordIterator` (a `java.util.PriorityQueue` replica ordered by `SAMRecordQueryNameComparator`, RG / PG
  tags translated on collisions). Several coordinate-sorted inputs (nf-core/sarek gives MarkDuplicates one BAM per
  lane / FASTQ chunk of a sample): each 2 Mb chunk reads every input that has data there and merges them by
  `SAMRecordCoordinateComparator` (ties to the earlier input), so the rest of the algorithm sees Picard's merged
  stream. `tests/harness/markdup/multi/cmp.sh` (inputs and truth: `make.sh`): NA12878 chr20:10-15M as 12 lane BAMs (in
  order and reversed), one lane as 3 chunks sharing a read group (colliding `@PG` IDs), chunks + lanes, and all of
  chr1 as 12 lane BAMs: records, metrics and headers identical to Picard.
* `OPTICAL_DUPLICATE_PIXEL_DISTANCE`, `CLEAR_DT`, `ADD_PG_TAG_TO_READS` (with neither CLEAR_DT nor PG tagging,
  Picard sets no attribute and htsjdk writes each record's tag bytes as read), `CREATE_MD5_FILE`;
  `VALIDATION_STRINGENCY` is ignored (variant-rs does not validate; it only decides whether Picard stops on malformed
  input).
* `tests/harness/markdup/cmp_warp.sh` (truth: `tests/harness/markdup/truth_warp.sh`): NA12878 chr20 query-grouped (WARP's arguments),
  the same split into two inputs with different read groups and a colliding `@PG` ID, a `SO:queryname` input with
  defaults, and a coordinate BAM with DT tags (the four options; and default re-marking): records (17,705,654 and
  1,192,285) byte-identical, headers identical apart from our own `@PG` line, metrics identical, `.md5` correct.
  Query-grouped chr20 13 s, two inputs 19 s.
* The output header is written as htsjdk does: @HD, @SQ, @RG, @PG, @CO, whatever order the input had them in
  (samtools collate writes @RG after @PG; this applies to ApplyBQSR too).

Times: smoke 3.2 s vs Picard 9.5 s; chr20 20 s (8 threads, machine shared) vs 1.7 min; whole genome **6.0 min**
(16 threads, 7.8 GB peak, ~50 GB of temporary bucket files; machine shared with two GATK jobs) vs Picard 42 min
(14 GB).

## Hard filters (`variant-rs filter`; `src/filter.rs`)

The output of GATK's single-sample hard-filter workflow in one pass: `SelectVariants -select-type SNP` and
`-select-type INDEL`, `VariantFiltration` with one `--filter-name` per expression (SNPs: QD2, QUAL30, SOR3, FS60,
MQ40, MQRankSum-12.5, ReadPosRankSum-8; indels: QD2, QUAL30, FS200, ReadPosRankSum-20), `MergeVcfs`.
Record type as htsjdk's `getType()` (`*` counts as a 1-base allele, so `*` with an insertion is MIXED): MNP,
MIXED and symbolic records are dropped exactly as the GATK workflow drops them (`--keep-other` keeps them
unfiltered). JEXL semantics: INFO values compared as doubles, QUAL from the column, a missing annotation never
fails. FILTER = existing filters + failing names, sorted; PASS if none. Header in htsjdk's merged order.

| VCF (from `variant-rs`) | Records in / out (dropped) | vs GATK workflow |
|---|---|---|
| HG002 35x WGS | 5,108,955 / 5,101,025 (7,930) | **records and header identical** (except command lines); 5.9 s vs 86 s |
| HG002 PCR+ 30x WGS | 5,083,936 / 5,076,883 (7,053) | **identical** |
| HG002 exome | 84,873 / 84,803 (70) | **identical** |
| NA12878 30x WGS (from CRAM) | 5,064,948 / 5,056,847 (8,101) | **identical** |

End to end (`pipeline/germline.sh --rust` vs the all-GATK pipeline, NA12878 smoke FASTQs, chr20:10-15 Mb):
duplicate-marked BAM records, recal table and final filtered VCF (8,023 records) identical.

## Output and header

* bgzip (BGZF blocks of 0xff00 bytes, libdeflate, compressed in parallel) and a tabix `.tbi` built while writing
  (htslib's binning/linear index, `INFO/END` honoured, as `tabix -p vcf`). Checked with `bgzip -t`, `tabix -l` and
  region queries against `zcat | awk`.
* The header is generated: GATK's HaplotypeCaller header lines (`tables/vcf_header_template.txt`,
  `tables/gvcf_header_template.txt`: INFO/FORMAT/FILTER/ALT/GVCFBlock/source, taken from GATK's own output), the
  BAM's sequence dictionary as `##contig=<ID=..,length=..[,assembly=AS]>` lines (htsjdk's
  `VCFHeader.setSequenceDictionary` on the reads' dictionary, all contigs), and the `@RG SM` sample.
  **Differences from GATK's header:** GATK's `##GATKCommandLine=<ID=HaplotypeCaller,CommandLine="...",Version=...,Date=...>`
  is replaced, at the same position, by `##variant-rsCommandLine=<ID=variant-rs,CommandLine="<variant-rs argv>",Version="0.1.0",Date="<UTC>">`.
  Every other header line is byte-identical (checked on all tests below, 3,366 contig lines included).

## HaplotypeCaller arguments beyond the defaults (S3 group 2, `tests/harness/hg002/group2/`)

WARP's HaplotypeCaller command lines (`tasks/wdl/GermlineVariantDiscovery.wdl`, `HaplotypeCaller_GATK4_VCF`:
`-contamination X -G StandardAnnotation -G StandardHCAnnotation [-G AS_StandardAnnotation] -GQB 10 ... -GQB 90
[-ERC GVCF]`) and a few related arguments:
* **`-G StandardAnnotation`, `-G StandardHCAnnotation`**: HaplotypeCaller's default groups; naming them changes
  nothing (`GATKAnnotationPluginDescriptor` merges annotations by class name). Other groups are refused.
* **`-G AS_StandardAnnotation`** (GVCF mode): `AS_RAW_BaseQRankSum`, `AS_RAW_MQRankSum`, `AS_RAW_ReadPosRankSum`
  (per alt allele, the Mann-Whitney z against the reference as a one-value histogram, `NaN` for an allele without
  reads, nothing when the reference has none: `||`), `AS_RAW_MQ` (sum of squared MQs per allele), `AS_SB_TABLE`
  (strand counts per allele, zero unless more than 2 reads), plus GATK's 7 header lines (AS_QD and
  AS_InbreedingCoeff are declared but never written for one sample). In VCF mode GATK refuses AS annotations
  ("not yet supported in the VCF mode", exit 2); so does variant-rs.
* **`-GQB`** (`--gvcf-gq-bands`): the given bands replace the default ones (`GVCFBlockCombiner.parsePartitions`,
  same checks and messages), `##GVCFBlock` header lines regenerated; header lines are sorted as htsjdk's writer does.
* **`-contamination F`**: `AlleleLikelihoods.contaminationDownsampling` before genotyping (best allele per read,
  `AlleleBiasedDownsamplingUtils.targetAlleleCounts`, random reads removed per allele); annotations still use every
  read, as GATK. The reads are drawn from GATK's second global generator, commons-math `Well19937c` seeded
  47382911 (`Utils.getRandomDataGenerator()`, `RandomDataGenerator.nextPermutation`, ported in
  `variant-rs/crates/genotype/src/well.rs` and checked against commons-math through jshell), which only this step uses. It runs
  through genotyping in GATK's order: contigs of a domain hand it on one after the other (a contig's genotyping
  waits for the previous contig's to finish). `-contamination 0` (WARP's default) changes nothing.
* **`--dbsnp`**: the ID column and the `DB` flag (`VariantOverlapAnnotator`: dbSNP records starting at the call,
  not filtered, sharing a biallelic event after `splitVariantContextToEvents` + left trimming + `Event`'s minimal
  representation; several IDs joined with `;`), read through the tabix index (`src/dbsnp.rs`).
* **`--interval-padding`, `-XL`, `--interval-exclusion-padding`**: `IntervalArgumentCollection.parseIntervals`
  (padded `-L`, minus padded `-XL`; `-XL` alone subtracts from every contig).

`tests/harness/hg002/group2/cmp.sh` (truth: `truth.sh`, HG002 chr8:1-15M unless noted) — all identical to GATK, data lines
and headers:

| case | GATK arguments | records |
|---|---|---|
| warp_gvcf | WARP GVCF command (AS, `-GQB 10..90`, `-contamination 0`) | 1,699,395 |
| warp_gvcf_shard023 | WARP GVCF command on the whole shard023 (chr8:1-74,211,328) | 8,234,095 |
| warp_vcf | WARP VCF command | 48,066 |
| dbsnp_vcf / dbsnp_gvcf | `--dbsnp` dbsnp138 (43,730 IDs, 1,317 with several) | 48,066 / 2,314,612 |
| contam_vcf / contam_gvcf | `-contamination 0.05` | 47,574 / 2,314,626 |
| contam3_vcf / contam3_gvcf | `-contamination 0.05`, chr8+chr9+chr10 1-3M | 20,668 / 1,437,831 |
| iv_vcf | two `-L`, `-ip 500`, `-XL`, `-ixp 50` | 4,756 |
| as_vcf | `-G AS_StandardAnnotation` without `-ERC GVCF` | refused by both |

## Mutect2 arguments beyond the defaults (S3 group 3, `tests/harness/somatic/group3/`)

* **Tumor-only** (no `-normal`): isActive's germline test (a germline-resource record covering the locus with
  AF >= 0.01 makes it inactive when the tumor's likeliest alt is an indel and the record's alt is an indel, or the
  likeliest alt is a base equal to the record's *reference* first base, as GATK compares it), the tumor-only default
  AF 5e-8, no NLOD / NALOD.
* **`--genotype-germline-sites`** (skips the normal / germline test in isActive and genotypes alleles the normal
  supports; the fraction argument stays at 1.0, its draws come from an unseeded generator in GATK),
  **`--genotype-pon-sites`**, **`--af-of-alleles-not-in-resource`**.
* **`--mitochondria-mode`** (`ModeArgumentUtils`, only for arguments not given): default AF 4e-3, emission and
  initial LOD 0, `--recover-all-dangling-branches` (ported in `variant-rs/crates/assembly`: heaviest-edge paths through
  forks, a k-mer graph with cycles after recovery is rejected), adaptive pruning LOD ln(1e-4), and the
  `OriginalAlignment` annotation (OCM: alt reads whose `OA` tag names another contig) unless `--annotation` is
  given. WARP's mitochondria command (`AlignAndCall.wdl`) adds `--read-filter MateOnSameContigOrNoMappedMateReadFilter
  --read-filter MateUnmappedAndUnmappedReadFilter` (applied after the tool's filters), `--annotation
  StrandBiasBySample` (a default annotation; it only stops the mode adding OriginalAlignment) and
  `--max-mnp-distance 0`. Giving a mode-set argument explicitly (`-emit-lod`, `-init-lod`, ...) is refused.
* **`--bam-output`** (`src/mutect/bamout.rs`, HaplotypeBAMWriter CALLED_HAPLOTYPES): per region with calls, the
  called haplotypes as reads (`HC_<n>_<i>`, `i` = position in the iteration order of the `HashSet` GATK copies
  them into), then every likelihood read realigned to its best haplotype, with the original record's other fields
  and tags plus `HC` (informative reads: `BestAllele.isInformative` compares with 0.2 even in natural log), `AR`,
  `CR`, `XA` (supported alleles appended to bwa's own `XA`), `BI` / `BD` (overlap cleaning); tags re-encoded as
  htsjdk does, records sorted stably with `SAMRecordCoordinateComparator`, header = the reads' (merged for two
  inputs, RG / PG tags translated) plus `ArtificialHaplotypeRG` and `@PG ID:HaplotypeBAMWriter`, `.bai` beside it.
* Found on the way, shared with ApplyBQSR / MarkDuplicates: htsjdk rewrites `@RG DT` in the machine's time zone
  (`Iso8601Date`); `bamwrite::htsjdk_rg_line` does the same (glibc `localtime_r`).

`tests/harness/somatic/group3/cmp.sh` (truth: `truth.sh`): all identical to GATK (VCF, `.stats`, bamout records and header):

| case | GATK arguments | records |
|---|---|---|
| gg / gp / af | HCC1143 tumor-normal + PoN + gnomAD with `--genotype-germline-sites` / `--genotype-pon-sites` / `--af-of-alleles-not-in-resource 0.0001` | 3,274 / 345 / 181 |
| to / to_gg | tumor-only (+ both genotype flags) | 3,159 / 9,646 |
| bamout | tumor-normal with `--bam-output` | 181 (bamout 34,328) |
| mito / mito_warp | NA12878 chrM `--mitochondria-mode` / AlignAndCall.wdl's command | 32 / 36 |

**FilterMutectCalls** `--mitochondria-mode` (ln priors of a somatic SNV / indel -2.5·ln10 / -3.75·ln10 unless
given; no clustered_events, multiallelic, fragment, slippage, haplotype or germline filter), `--min-allele-fraction`,
`--max-alt-allele-count`, `--f-score-beta`. WARP runs it twice (`AlignAndCall.wdl`: InitialFilter, then
FilterContamination on that output with haplochecker's estimate). `tests/harness/somatic/group3/fmc_cmp.sh` (truth:
`fmc_truth.sh`), VCF and filtering stats identical to GATK:

| case | GATK arguments | records (PASS) |
|---|---|---|
| warp1 / warp2 | AlignAndCall.wdl's two Filter calls on mito_warp (warp2 re-filters warp1: `--min-allele-fraction 0.03 --f-score-beta 1.5 --contamination-estimate 0.02`) | 36 (22) / 36 (17) |
| mito | mito with `--mitochondria-mode` | 32 (20) |
| to_mito / to_opts | tumor-only calls, mito mode / `--max-alt-allele-count 2 --min-allele-fraction 0.05 --f-score-beta 2.0` | 3,159 (592) / 3,159 (233) |
| tn_opts | gg, `--max-alt-allele-count 2 --min-allele-fraction 0.02 --f-score-beta 0.7` | 3,274 (49) |

## GATK semantics reproduced (new in the whole-genome driver)

* **Shards:** one read shard per contig with all of the contig's intervals (`AssemblyRegionWalker.makeReadShards`,
  `MultiIntervalLocalReadShard`). Reads are queried over the intervals padded by 100 and merged
  (`IntervalUtils.getIntervalsWithFlanks`), through the BAI, each record once, in file order, with htsjdk's
  `BAMQueryMultipleIntervalsIteratorFilter` overlap test (alignment end = start + reference length − 1; a placed
  unmapped read uses its start). Loci are only the (unpadded) interval loci; a non-contiguous locus forces the
  activity profile out; reads in interval gaps still fill regions whose padded span reaches them.
* **Reference:** upper-cased and IUPAC codes → `N` (`CachingIndexedFastaSequenceFile`); 87 FASTA lines of this hg38 analysis set contain IUPAC codes (chr1–chr22, chrX).
* **Mates on another contig** (grep of GATK for `getMateContig` / `mateIsUnmapped` on HC's path):
  * `filterNonPassingReads` (`MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE`) in genotyping and in the reference model:
    reads carry `ref_id` / `mate_ref_id` from the BAM (`hc_likelihood`, `hc_refconf::passes_hc_read_filter`).
  * `ReadCoordinateComparator` ties (preproc, likelihood, refconf): `getSAMFlagsForRead` (GATK *rebuilds* the flags)
    and `getMateReferenceIndex` (-1 if the mate is unmapped, else its contig index).
  * `hasWellDefinedFragmentSize`, `getAdaptorBoundary` (pileup adaptor exclusion here and in preproc/refconf) and
    `FragmentCollection`: GATK has **no** contig test there, only `TLEN == 0` and the mate *position*; reproduced as is.
* **Positional downsampling** (`hc_downsample`): after the read filters, per contig, kept reads in reservoir order
  feed both the pileups and the region read cache.

## Validation (NA12878 30x WGS BAM, GRCh38 + decoys + HLA; GATK 4.7 gatk.jar)

All comparisons are raw text (`zcat`, data lines; header separately). `scratch/cmp.sh NAME vcf|gvcf|both -L ...` (development script, not in the repository)
runs GATK and variant-rs on the same intervals; `scratch/ds_regress.sh` (development script, not in the repository) replays the downsampler's GATK runs.

| Test | Intervals | VCF data lines | GVCF data lines |
|---|---|---|---|
| chr20 slice (`$VARCALL_DATA/data`, old BAM) | chr20:10–15 Mb | 8,064 identical | 798,047 identical (+ 5,000,001 raw records) |
| t0 | chr20:10.0–10.2 Mb | 454 identical | 30,810 identical |
| t2 alts / decoys / HLA / unplaced | 13 whole contigs (4 HLA with `:` in the name, 4 decoys, 1 alt, 4 unplaced; 6.7 Mb). 342,093 reads pass MQ/flags, **34,369 with the mate on another contig** | 6,295 identical | 182,673: 1 line differs with CUDA (2 `<NON_REF>` PLs); **identical with `--pairhmm cpu`** |
| t3 downsampling | chrM + chrUn_KI270467v1 + chrUn_JTFH01001899v1_decoy: 1,584,628 reads → 922,455 kept, 17,777 drawing groups, 662,173 `nextInt` | 172 identical (= w3) | 868: 4 lines differ with CUDA (high-depth chrM PLs, one MIN_DP); **identical with `--pairhmm cpu`** |
| t4 contig ends + several intervals per contig | chr1:1-250000, chr20:1-120000, chr20:64334167-end, chr21:5010000-5150000, chr21:10000000-10050000 + chr21:10050150-10100000 (gap 150: separate loci, merged read query), chrX:155900000-end | 1,112 identical | 73,862 identical |
| t5 shard semantics | shard A = chr21:20–21 Mb + chr22:30–31 Mb; shard B = chr19:1.0–1.6 Mb + chrM:1-3000 | per-shard: A 3,856 and B 1,101 identical; `-O` concatenation identical; **single-run over A+B = one GATK JVM over A+B: 4,957 identical** (GATK's own per-shard and single-JVM outputs differ in 411 lines: RNG carry-over) | — |
| downsampler's GATK runs (`scratch/ds_regress.sh` (development script, not in the repository)) | w1, w2, w3 (WGS, cap 50); t1 (cap 2), g1 (cap 2, GVCF), t2 (cap 1) on the chr20 slice | 6 / 3,838 / 172 / 573 / 8,030 identical | g1: 47,007 identical |
| **whole genome, per-shard mode vs `wgs/truth/vcf`** | shard00–shard13 (chr1–chr14, the truth shards finished at the time) | 3,815,005 lines: shards 02, 04, 06, 09, 10 identical; 15 records differ in the other nine (≤ 3 per shard). Re-running each of the 15 sites on a ±250 kb window with `--pairhmm cpu` (no drawing group within 300 kb) gives GATK's QUAL/PL/DP/GQ/rank sums: **all 15 are CUDA PairHMM rounding** (QUAL/PL ±1, a GQ, a ReadPosRankSum, one read in or out of DP/MQ). Headers identical except the command line in all 14. | not compared at the time |
| **whole genome, bit-exact CUDA kernel (2026-11-12)** | all 23 shards, VCF and GVCF | **23/23 identical (5,064,948 lines)** | **23/23 identical (461,802,035 lines)** |
| **whole genome from the CRAM** (`wgs/cram/NA12878.final.cram`, 1000G/NYGC, CRAM 3.0) | all 23 shards, vs GATK on the BAM | **23/23 identical (5,064,948 lines)**, same read counts | chr20 slice and the t4 intervals + chrM/HLA/alt/decoy contigs: VCF and GVCF identical to the BAM run |
| **GATK on the CRAM** (`wgs/cramtest`) | htsjdk decoding (PrintReads) of chr20, chrM, alt/decoy/HLA contigs and slices of chr21/chrX vs htslib: 20,020,307 records; HaplotypeCaller on shard22 (2,968 alt/decoy/HLA contigs) and shard20 (chr21) | records identical (sorted tags); **GATK-CRAM = GATK-BAM: 8,871 and 158,940 records identical** | — |
| **whole genome, `variant-rs plan` 97 shards** (`wgs/plan96`, GATK truth `wgs/truth96`) | 97 shards, 22 contigs split mid-contig | **97/97 identical (5,064,997 lines)**, headers identical except the command line | GATK not run |

Headers: identical to GATK's except the command-line line in every test (3,366 `##contig` lines on the WGS BAM).
Tabix indexes: `tabix -l` and region queries agree with `zcat | awk` (chr20, GVCF with `END`, multi-contig), and on
HG002 (split contigs) with an index built by `tabix -p vcf` for all 195 contigs plus region queries, VCF and GVCF.

**Second sample: HG002** (`tests/harness/hg002`; NovaSeq 2x151 PCR-free BAM from the DeepVariant case study, ~35x, GRCh38 no-alt analysis set, 834 M reads; plan:
47 shards from 195 contigs, 14 contigs split):

| Test | Result |
|---|---|
| Rust vs GATK 4.7, 47 shards, VCF | **47/47 identical (5,108,955 lines)**, headers identical except the command line |
| Rust vs GATK 4.7, 47 shards, GVCF | **47/47 identical (456,525,374 lines)** |
| GIAB HG002 v4.2.1, chr1–22, `rtg vcfeval --ref-overlap` | Rust = GATK: precision 0.9887, recall 0.9923, F 0.9905 (TP 3,860,537, FP 44,313, FN 29,903) |

**HG002 exome** (`$VARCALL_DATA/hg002_exome`; NovaSeq IDT capture, ~100x, DeepVariant exome case study; targets + 100 bp padding,
merged, 178,603 intervals / 77.5 Mb in 16 shards built by `make_shards.py`):

| Test | Result |
|---|---|
| Rust vs GATK 4.7, 16 shards (~11,000 intervals each), VCF | **16/16 identical (84,873 lines)** |
| Rust vs GATK 4.7, GVCF | **16/16 identical (3,385,734 lines)** |
| GIAB HG002 v4.2.1 within targets ∩ benchmark (36.1 Mb) | Rust = GATK: precision 0.9600, recall 0.9847, F 0.9722 (TP 26,095, FP 1,087, FN 406) |

The downsampler barely triggers here (4,865 of 53.9 M reads removed, 217 drawing groups): 100x of 150 bp reads is
< 1 read per start position, so a deduplicated exome only exceeds the cap of 50 in pileups; high-depth panels or
amplicon data would stress it more.

**HG002 PCR+** (`$VARCALL_DATA/hg002_pcrplus`; NovaSeq TruSeq PCR+ library, 30x, hs38DH; Google's "gold-standard samples" dataset,
`gs://brain-genomics-public/research/sequencing/grch38/bam/novaseq/wgs_pcr_plus/30x/`; `variant-rs plan --pieces 48`: 46
shards, 13 contigs split):

| Test | Result |
|---|---|
| Rust vs GATK 4.7, 46 shards, VCF | **46/46 identical (5,083,936 lines)** |
| Rust vs GATK 4.7, 46 shards, GVCF | **46/46 identical (583,697,036 lines)** |
| GIAB HG002 v4.2.1, chr1–22 | Rust = GATK: precision 0.9761, recall 0.9841, F 0.9801 (TP 3,828,531, FP 93,848, FN 61,930); PCR-free HG002 above: F 0.9905 |

Downsampling is busier than in the PCR-free genome: 3.7 M of 630 M reads removed, 44,394 drawing groups.

HG002 found two bugs that NA12878 did not exercise: haplotypes over 1,055 bp exceeded the CUDA kernel's shared memory
(now fewer warps per block up to 4,223 bp, CPU PairHMM beyond; none fell back on HG002), and the merged `.tbi` listed
a split contig once per piece, so readers saw only its first piece (now one entry per contig). The VCF data was never
affected.

The crate validators still pass after the mate-contig changes (preproc 100%, assembly 6,345/6,345, likelihood all
checks 6,270/6,270, genotype 8,064/8,064 and `--erc` 8,905/8,905, refconf `validate.sh`).

## Performance

Machine: 16 cores / 32 threads, RTX 5070 Ti.

**Bit-exact kernel, whole genome (2026-11-12)**, wall clock from `/usr/bin/time`:

| Run | Machine state | Wall | CPU (user) | Peak RSS |
|---|---|---:|---:|---:|
| 23 shards, VCF, `--threads 16` (8 contig workers) | idle | 12:51 | 11,044 s | 9.3 GB |
| 23 shards, GVCF, `--threads 16` | idle | 17:04 | 14,689 s | 11.2 GB |
| 23 shards, VCF, `--threads 32` (16 workers) | idle | 11:58 | 15,223 s | 13.8 GB |
| 23 shards, VCF, **from CRAM**, `--threads 16` | idle | 13:54 | 13,550 s | 7.5 GB |
| 97 shards (`plan96`), VCF, `--threads 32` | CPU thermal throttling | 12:41 | 15,963 s | 9.4 GB |
| 97 shards (`plan96`), GVCF, `--threads 32` | CPU thermal throttling | 18:41 | 22,711 s | 14.5 GB |
| GATK 4.7, 23 shards, 5 JVMs × 2 threads, VCF / GVCF | — | 2 h 11 m / 2 h 23 m | 9.7 h / 10.7 h summed | |
| GATK 4.7, 97 shards, 4 JVMs × 2 threads, VCF | CPU thermal throttling | 1 h 57 m | 7.6 h summed | |
| **HG002**, 47 shards, VCF, `--threads 16` | idle | 15:24 | 13,343 s | 6.5 GB |
| **HG002**, 47 shards, GVCF, `--threads 16` | idle | 20:00 | 17,239 s | 10.4 GB |
| GATK 4.7 on HG002, 47 shards, 5 JVMs × 2 threads, VCF / GVCF | — | 2 h 19 m / 2 h 41 m | 10.6 h / 12.4 h summed | |
| **HG002 exome**, 16 shards, VCF / GVCF, `--threads 16` | idle | 0:41 / 0:52 | 605 s / 857 s | 3.3 / 5.9 GB |
| GATK 4.7 on the exome, 16 shards, 5 JVMs × 2 threads, VCF / GVCF | — | 6 m 18 s / 8 m 27 s | 1,773 s / 2,371 s summed | |
| **HG002 PCR+**, 46 shards, VCF / GVCF, `--threads 16` | idle | 15:06 / 20:24 | 13,504 s / 17,553 s | 6.0 / 10.1 GB |
| GATK 4.7 on HG002 PCR+, 46 shards, 5 JVMs × 2 threads, VCF / GVCF | — | 1 h 52 m / 2 h 36 m | 9.2 h / 12.7 h summed | |

CRAM decoding takes ~6,200 s summed over the reader threads (BGZF inflate of the BAM: ~580 s) and costs ~8% wall
clock. Doubling the threads gains 7%, and 97 shards gain nothing over 23 here: the run is bound by total CPU work (front
end above all), not by its longest contig. The GPU is busy ~220 s of a ~12 min run.

**Earlier runs, approximate-kernel build.** The numbers below were measured while 5 GATK JVMs
(the truth jobs, 2 threads each) ran on the same machine** (load average 20–27), so they are upper bounds.
Times are variant-rs's own monotonic clock (`time: total`).

| Run | Wall | CPU (user) | Peak RSS | Notes |
|---|---:|---:|---:|---|
| chr20:10–15 Mb VCF (`$VARCALL_DATA/data`), 3 runs | 2.00–2.12 s | | | was 3.1 s with the old driver; GATK 33.8 s |
| chr20:10–15 Mb GVCF, 3 runs | 2.92–3.00 s | | | was 4.3 s; GATK 54.2 s |
| chr21, VCF (46.7 Mb) | 38.4 s | | 1.3 GB | single contig: its front end and its sequential batches bound it |
| shards 01–04 (chr2–chr5, 812 Mb), VCF | 242 s | 3,006 s | 4.8 GB | GATK: 1,649–2,211 s per shard |
| **whole genome, 23 shards, VCF** | **1,040 s** | 13,584 s | 9.1 GB | 5,064,948 records; GATK needed 25,013 JVM-s for chr1–chr14 alone (1,405–2,887 s per shard) |
| **whole genome, 23 shards, GVCF** | **1,385 s** | 18,288 s | 11.4 GB | 461,802,035 lines, 5.5 GB `.g.vcf.gz` + 4.1 MB `.tbi`, written twice (merged + per shard) |

Where the time goes (whole genome; summed per-contig stage wall times, which overlap across the 8 contig workers;
CPU-bound overall at ~13 cores of the ~22 left by the GATK jobs):

| Stage | VCF | GVCF |
|---|---:|---:|
| front end (decode, filters, downsampling, pileups, regions; includes waiting on its back end) | 7,263 s (BGZF inflate 471 s) | 9,768 s (536 s) |
| finalizeRegion | 998 s | 1,420 s |
| assembly | 2,110 s | 1,822 s |
| trimming + likelihood prep / finish | 676 / 633 s | 626 / 832 s |
| PairHMM (GPU, one global lock) | 841 s; GPU busy 221 s (kernel 175 s), 1.57 G pairs, 3.2e13 cells | 1,188 s; GPU 239 s (kernel 187 s) |
| realignment | 1,094 s | 909 s |
| genotyping | 176 s | 273 s |
| reference confidence | — | 2,078 s |
| GVCF banding + bgzip (writer threads) | — | 1,059 s (deflate 527 s) |

* The longest contigs set the floor: chr1 took 615 s (VCF) / 819 s (GVCF) from start to finish. Each contig's front
  end is one thread (plus its inflate and decode threads) and its back end processes one batch at a time.
* The baton chains are short: shard21 (375 contigs incl. chrY, chrM) finished at 180 s and shard22 (2,968 contigs)
  at 272 s in the GVCF run; the end of the run is the bulk of shard13–shard19 on CPU.
* Downsampling over the genome: 635,682,919 reads → 632,847,114 kept, 40,524 drawing groups (VCF run).

**GPU PairHMM engine pool (2026-11-19).** HaplotypeCaller's contig workers shared one PairHMM engine behind one
lock (one CUDA stream, one pinned buffer), with the packing of a batch on the CPU inside the lock: on the L40S, 30x
whole genome, 64 threads, the callers waited 2,543 s in all for it, and more threads made the run slower (212 s at
64, 250 s at 384; GPU busy ~40%). Now the device setup is shared and a pool of engines (`VARIANT_RS_GPU_STREAMS`,
default 8), each with its own stream and buffers, serves the callers; each pair is computed on its own, so the
results cannot change (VCF identical in every run below, selftest identical with the CUDA PairHMM). With the queue
gone, more contig workers pay off where threads are plentiful:

| Machine | Contig workers | Engines | Wall | GPU busy | Peak RSS |
|---|---:|---:|---:|---:|---:|
| L40S, 64 threads | 8 | 1 (before) | 189.9 s | 48% | 8.4 GB |
| L40S, 64 threads | 8 | 8 | 159.7 s | 48% | |
| L40S, 64 threads | **16 (new default)** | 8 | **121.0 s** (defaults, 2 runs: 121.8, 122.7) | 61% | 15.0 GB |
| L40S, 128 threads | 32 (new default) | 8 | 113.6 s | 64% | |
| RTX 5070 Ti, 32 threads | 8 (default) | 1 (before) | 396.7 s | 34% | 8.6 GB |
| RTX 5070 Ti, 32 threads | 8 (default) | 8 | 401.2 s | 33% | 8.2 GB |
| RTX 5070 Ti, 32 threads | 16 | 8 | 403.3 s | 34% | 14.8 GB |

On 32 threads the CPU stages set the pace (front end, assembly), so the defaults there are as before: contig
workers are half the threads up to 8 and, beyond 32 threads, a quarter of them up to 32 (`--contig-workers`
overrides).

**Fewer passes per batch (2026-11-21).** On 32 threads (RTX 5070 Ti) the run was CPU-bound with only ~23 threads
busy, and adding parallelism made it slower (more batches in flight, deeper queues, 12 or 16 contig workers, a
48-thread pool: 449-476 s against 445): the rayon pool woke and slept ~38,000 times a second, one parallel pass per
stage per batch. Each region's work before the PairHMM (finalizeRegion, assembly, PairHMM inputs) is now one task,
and its work after it (likelihoods, realignment) another: 30x genome, 456 -> 400 s (2 runs each), ~15,000 switches a
second, VCF identical. Assembly on chr20 (`VARIANT_RS_ASM_PROFILE=1`): graph building 60% (read threading 52 of 54 s),
paths and haplotype cigars 39% (the haplotype Smith-Waterman 32 of 35 s; 43 G cells with prefix reuse). Workstation
runs: `$VARCALL_DATA/runs/stage6/hc`.

**Several GPUs in one process (2026-11-19).** The engines are spread over every visible GPU (`VARIANT_RS_GPUS` for
at most N): each GPU gets its own copy of the tables (constant memory) and kernel attributes, each engine belongs to
one GPU, and a call takes a free engine on the GPU with the fewest calls running. HaplotypeCaller, 30x genome, L40S,
128 threads (32 contig workers), VCF identical in every run:

| GPUs | Wall | GPU busy (mean) | Callers' wait for an engine (summed) |
|---:|---:|---:|---:|
| 1 | 114.9 s | 63% | 1,312 s |
| 2 | 97.9 s | 41% | 165 s |
| 4 | 97.3 s | 21% | 4 s |
| 6 | 99.5 s | 15% | 6 s |

Beyond 2 GPUs the longest calling interval sets the time: with `--scatter-intervals` each interval is its own GATK
run, sequential, and the calling-regions list has one of 141 Mb (chr2:97,489,619-238,903,659, ~93 s). The pipeline
gives HaplotypeCaller up to 2 GPUs (`hc_gpus`).

## Earlier single-contig results (chr20:10–15 Mb, `$VARCALL_DATA/data`)

Byte-identical to GATK: 8,064 VCF records; the 798,047-line GVCF and the 5,000,001 unbanded records
(`bash variant-rs/crates/refconf/validate.sh (development script, not in the repository)`, header compared without the command-line line). Region detection:
5,000,001 per-locus isActive values and 24,567 regions identical (`--assembly-region-out` vs GATK's IGV file).
16 threads + CUDA: 2.0–2.1 s (VCF) and 2.9–3.0 s (GVCF) wall, vs GATK's 33.8 s and 54.2 s.

## GATK behaviours of region detection

* **Integer PL quantization:** the ref-vs-any likelihoods go through `GenotypeBuilder.PL(double[])`,
  i.e. they are rounded to integer Phred PLs (normalized to the best genotype) before the posterior is computed.
* **Adaptor exclusion in pileups:** LocusIteratorByState drops bases past `ReadUtils.getAdaptorBoundary`
  (reverse-strand read: positions ≤ mate start − 1). The `fragmentLength > 100` shortcut compares the *signed*
  TLEN, so it never applies to reverse reads.
* **Band-pass quirk:** when an HQ-soft-clip state is spread over ±N loci, GATK centres every resulting Gaussian
  on the *just-added* locus rather than on each spread state's own locus.
* Deletions count with fixed quality 30 and are never quality-filtered; bases with Q ≤ 10 are skipped; the
  HQ-soft-clip running average uses quality > 28.
* Numeric tables (Phred→log10 cache, Jacobian log-sum table, Gaussian kernel, genotype priors) are generated by
  GATK itself (`tables/GenTables.java`) and embedded.

Why it is fast: no pileup objects (each read's CIGAR is walked once and scattered into a ring of per-locus
accumulators, in file order, so sums are bit-identical); per-read precomputation; batched GPU PairHMM.

## Several `-I` in HaplotypeCaller (S3 group 4, `tests/harness/hg002/group4/`)

Several inputs of one sample (lanes, batches): each input is read, interval-filtered and read-filtered on its own
thread, then merged as htsjdk's `MergingSamRecordIterator` does (`SAMRecordCoordinateComparator` over every record
the interval filter returns, filtered or not, since GATK merges before its read filters; ties to the earlier input).
Inputs holding several samples are called jointly (next section). Refused, with the reason (or run by GATK under
`GATK_FALLBACK`, now also for refusals found after opening the inputs): different sequence dictionaries,
`--bqsr-recal-file` with several inputs, and the same read in two inputs. GATK's reads are `equals()` by their
fields, so byte-identical copies collapse in some of its hash-keyed read maps but not others (one BAM given twice:
INFO DP / MQ off by one read at 9,374 sites); that is not reproduced.

`tests/harness/hg002/group4/cmp.sh` (truth: `truth.sh`; inputs: `make_inputs.sh`, HG002 chr8:1-15M split by read name so
mates stay together) — all identical to GATK, data lines and headers:

| case | inputs, GATK arguments | records |
|---|---|---|
| lanes_vcf / lanes_gvcf | L1 + L2 (read groups L1 / L2, sample HG002), WARP VCF / GVCF command | 48,066 / 1,699,401 |
| same_rg_vcf | S0 + S1, both with the original read group | 48,066 |
| three_vcf / three_rev_gvcf | T0 + T1 + T2 / T2 + T0 + T1, GVCF | 48,066 / 2,314,614 |
| contam_vcf | L1 + L2, `-contamination 0.05` | 47,574 |
| other_vcf | L1 (HG002) + O (L2's reads as sample OTHER), two samples, WARP VCF command | 48,148 |

## Several samples: joint calling (S3 group 4, `tests/harness/hg002/group4/ms/`)

Inputs whose read groups name several samples (in one file or several) are called jointly, as GATK does in VCF
mode (GATK refuses `-ERC GVCF` with several samples, and so does variant-rs). Samples are GATK's sorted sample list;
each read's sample comes from its read group (per input, so read-group IDs may repeat across inputs). What changes:

* isActive: one ref-vs-any likelihood accumulator per sample and locus. One sample with pileup elements: the
  single-sample posterior; several: `MinimalGenotypingEngine.calculateGenotypes` on the biallelic `N` /
  `<FAKE_ALT>` site (the symbolic alt has length 0, so the indel prior) with stand-call-conf 4, QUAL to
  probability. The HQ soft-clip running average runs over the samples in `splitContextBySampleName` order (a Java
  `HashMap`: bucket of the spread `String.hashCode`, ties in sorted order), so the values are kept per element.
* Overlapping mates are paired within their sample; assembly prunes on the largest per-sample multiplicity (already
  ported for Mutect2).
* Genotyping: one genotype per sample (no evidence: PLs 0,0,0), `AlleleFrequencyCalculator` over all samples
  (`calculate_multi`; one sample gives the old result bit for bit), `calculateMostLikelyAlleles` over the sum,
  USE_PLS_TO_ASSIGN per sample; contamination downsampling per sample in sample order.
* Annotations: AD / DP for called genotypes only; AC / AF / AN over called chromosomes; ExcessHet over all
  genotypes (rounded counts); InbreedingCoeff from 10 samples (unrounded counts); FS / SOR from per-sample strand
  tables each kept when its total exceeds the minimum; QD from the het / hom-var genotypes; DP, MQ and rank sums
  over all evidence. FORMAT is the union of the genotypes' keys, missing values `.`, trailing ones dropped.

`tests/harness/hg002/group4/ms/cmp.sh` (inputs and truth: `make.sh`; chr8 of HG002, NA12878 and HG002's PCR+ library renamed
HG002pcr, all on the reference's dictionary) — all identical to GATK, data lines and headers:

| case | inputs, GATK arguments | records | variant-rs | GATK |
|---|---|---|---|---|
| ms2 | HG002 + NA12878, chr8:1-15M | 63,689 | 12.5 s | 8.9 min |
| ms3 | HG002 + NA12878 + HG002pcr, chr8:1-15M | 64,449 | 21.5 s | 14.6 min |
| ms11 | the three split by read name into 4 + 4 + 3 samples (InbreedingCoeff), chr8:1-5M | 20,794 | 27.7 s | 8.2 min |
| onebam | one BAM holding HG002 and NA12878, chr8:5M-10M | 22,000 | 5.1 s | 3.1 min |
| ms3_opts | three samples, chr8:1-5M, `--dbsnp`, `-contamination 0.05` | 20,623 | 8.0 s | 4.4 min |

## GatherBQSRReports (S3 group 4, `src/bqsr/gather.rs`, `tests/harness/bqsr/gather/`)

`RecalibrationReport.gatherReports`: each input report is parsed back into datums (observations, the printed
2-decimal errors, reported quality = `EstimatedQReported` in the read-group table, the quality score elsewhere),
read groups re-keyed by the sorted union of all inputs' names, tables combined left to right (`RecalDatum.combine`:
expected errors recompute the reported quality), no rounding afterwards, the quantization recomputed from the
combined quality table, the first input's Arguments table printed again; base-substitution tables only
(BaseRecalibrator's default). `tests/harness/bqsr/gather/cmp.sh` (truth: `truth.sh`): byte-identical to GATK on three chr20
scatter tables of NA12878 (as WARP scatters), three HG002 lane tables out of order, HG002 and NA12878 tables mixed
(five read groups), a single table and one table twice; ApplyBQSR with the gathered table is identical too.

## GatherBamFiles (S3 group 4, `src/gatherbam.rs`, `tests/harness/gatherbam/`)

Picard `GatherBamFiles` on BAM inputs is htsjdk's block copy (`BamFileIoUtils.gatherWithBlockCopying`): the first
input's blocks unchanged (EOF block dropped), then per further input the rest of the block holding its first record
deflated again into one block and its later blocks unchanged, one EOF block; no record is decoded. variant-rs does the
same, so the output differs from Picard's only in those re-deflated blocks (Picard's deflater vs libdeflate);
`CREATE_INDEX`: the `.bai` built by inflating the output (as htsjdk's `StreamInflatingIndexingOutputStream`),
`CREATE_MD5_FILE` as elsewhere. `tests/harness/gatherbam/cmp.sh` (truth: `truth.sh`): WARP's GatherSortedBamFiles (index, md5)
and GatherUnsortedBamFiles on ApplyBQSR scatter outputs (three contig stretches + unmapped reads), one input, and three
samtools-written BAMs: records byte for byte and headers identical, every BGZF block identical except the one
re-deflated block per further input (`blocks.py`), `.bai` answers identical, `.md5` right.

## MergeVcfs (S3 group 4, `src/mergevcfs.rs`, `tests/harness/mergevcfs/`)

Picard `MergeVcfs`: records merged by `VariantContextComparator` (contig index of the first header, then start)
through htsjdk's `MergingIterator`, a `java.util.PriorityQueue` keyed by each input's next record, simulated operation
for operation so records at equal positions in several inputs come out in Picard's order; an unsorted input is an
error, as in Picard. Header: `smartMergeHeaders` (each input's lines in sorted order, the first line per key wins:
`KEY-ID` for INFO / FILTER / FORMAT / contig / ALT / PEDIGREE / META / SAMPLE, the bare key otherwise, so one
`##GATKCommandLine` is kept), written sorted. Records re-encoded as htsjdk's `VCFEncoder` does (genotypes as read,
INFO keys sorted, `KEY=` -> `KEY`, QUAL `%.2f` without `.00`, FILTER sorted, alleles upper-cased). `.tbi` for
`.vcf.gz` (`CREATE_INDEX`, default true). Refused: INFO / FORMAT / FILTER header lines htsjdk would rewrite,
conflicting INFO / FORMAT definitions (Picard promotes them), an uncompressed output with an index (Tribble `.idx`),
inputs without contig lines. `tests/harness/mergevcfs/cmp.sh` (truth: `truth.sh`): text byte-identical and `.tbi` answers
identical on WARP's MergeVCFs of three HaplotypeCaller GVCF / VCF shards (HG002 chr8, WARP's commands), the shards out
of order, one shard twice, a VCF + GVCF + `-contamination` VCF of one stretch (15,323 positions in several inputs),
two FilterMutectCalls VCFs, and an uncompressed output. `SEQUENCE_DICTIONARY` (nf-core/sarek passes the reference's
`.dict` in every gather): Picard only fills contig-less inputs from it (refused here) and hands it to the index, so
the output is unchanged; the file is checked as Picard opens it. Its cases (sarek's command on the VCF and GVCF
shards, and a 3,366-contig dictionary over 195-contig headers) are identical to Picard too.

## GenotypeGVCFs (S3 group 4, `src/genotypegvcfs.rs`, `tests/harness/genotypegvcfs/`)

GATK `GenotypeGVCFs` on a GVCF file: a single-sample HaplotypeCaller or ReblockGVCF GVCF, or a multi-sample one from
CombineGVCFs; diploid; the default annotations or WARP's (`JointGenotypingTasks.wdl`: `-D dbsnp -G StandardAnnotation
-G AS_StandardAnnotation --only-output-calls-starting-in-intervals -L ... --merge-input-intervals`). WARP reads a
GenomicsDB workspace (`-V gendb://`, made by GATK's GenomicsDBImport): see below. By-variant
traversal, each record on its own: `ReferenceConfidenceVariantContextMerger.merge` (NON_REF dropped, PL / AD remapped,
genotypes no-called or best-matched, raw annotations combined, other numeric INFO typed by the input header),
`MinimalGenotypingEngine.calculateGenotypes` (the allele-frequency EM over all samples, GQ-only hom-refs of reblocked
GVCFs approximated from GQ; QUAL; output alleles incl. spurious-`*` tracking; PREFER_PLS genotypes; MLEAC / MLEAF /
AS_QUAL), `finalizeAnnotations` (AS_*RankSum histogram medians, AS_FS / AS_SOR from AS_SB_TABLE, AS_QD, AS_MQ, MQ from
RAW_MQandDP), reverse trim, `annotateContext` without reads (AC / AF / AN, ExcessHet, FS / SOR from FORMAT SB, QD,
InbreedingCoeff and AS_InbreedingCoeff from 10 samples, dbSNP ID and DB), genotype cleanup (SB and MIN_DP dropped,
PGT `1|1` for hom-var, AD filled). QD / AS_QD at or above 35 are jittered with GATK's global `java.util.Random(47382911)`
(AS_QD's draws first, per alt, then QD's), as in HaplotypeCaller. Annotation layer specified in
[specs/genotypegvcfs_annotations.md](../specs/genotypegvcfs_annotations.md) (checked against the truth). Found on the way: the non-AS rank sums are not
recomputed but written back through htsjdk's `formatVCFDouble` (negative values as `-4.310e-01`), and Java's `%.3e`
rounds the shortest digits HALF_UP (`hc_genotype::javamath::format_sci3` fixed: it used Rust's `{:.3e}`);
InbreedingCoeff uses unrounded genotype counts (no-PL hom-refs as 1/3 each at GQ 0, else `qualToProb(GQ)`), ExcessHet
rounded ones. `tests/harness/genotypegvcfs/cmp.sh` (truth: `truth.sh`, `multi/make.sh`, `multi/make11.sh`): byte-identical to
GATK on 11 cases: HG002 chr8:1-15M with WARP's arguments on the HaplotypeCaller and on the reblocked GVCF, the
defaults, a GVCF without AS annotations; three samples (HG002, NA12878, HG002 PCR+; raw and reblocked, CombineGVCFs)
with WARP's arguments and the defaults; eleven samples (copies under other names: InbreedingCoeff, AS_InbreedingCoeff).
`wgs_cmp.sh` (truth: `wgs_truth.sh`): WARP's arguments on the reblocked whole-genome HG002 GVCF (47 shards from
`tests/harness/reblock/wgs`, 5.1M calls): all identical, 26 s with the comparison.

### GenotypeGVCFs on a GenomicsDB workspace (`-V gendb://`, `src/gdb_store.rs`, `src/gdb_merge.rs`, `tests/harness/genomicsdb/`)

Workspaces written by GATK's GenomicsDBImport (GenomicsDB 1.5.5) are read directly; GenomicsDBImport itself stays in
GATK. Storage ([specs/genomicsdb_tiledb.md](../specs/genomicsdb_tiledb.md)): GenomicsDB's TileDB fork, one sparse array per imported interval
(`contig$start$end`, rows = samples by `callset.json`, columns = global positions by `vidmap.json`), live fragments
by marker file in timestamp order, `__book_keeping.tdb.gz`, zlib tiles with delta-encoded coordinates and var
offsets, cells merged over fragments (the newest wins), values decoded with GenomicsDB's sentinels, GT as
allele/phase pairs, allele-specific fields as packed blobs. Query / merge ([specs/genomicsdb_merge.md](../specs/genomicsdb_merge.md)): what GATK's
GenomicsDB reader hands the walker: queries in 200 kb windows (shard 100 kb + feature-cache lookahead 100 kb), each
an independent `scan_and_operate` (one call per sample, an overlapped call replaced, 1-bp records while a deletion or
MNV is active, records at every breakpoint and window start), `BroadCombinedGVCFOperator` (spanning deletions reduced
to `N` -> `*`, REF / ALT merge, AD / PL / AS remapping from NON_REF, INFO combine operations from the vid map: sums,
medians printed with htslib's `kputd`, element-wise and per-allele sums, histogram sums; INFO DP rule; GT called,
the GQ0 / PL0 rule; a sample without a call as a haploid `.`), then the same GenotypeGVCFs engine as for files.
Records are streamed (per window: the window's cells and each sample's last record before it). The storage and the
merge were specified and validated in Python first (every cell of the test workspaces re-derived from the input
GVCFs; the received stream rebuilt 100 %), then ported. `tests/harness/genomicsdb/cmp.sh` (truth: `truth.sh`, `truth2.sh`;
`RECEIVED=1` also compares the received stream with the harness dump in `received/`): byte-identical to GATK on 11
cases: 3 samples (reblocked, raw, imported in 3 fragments) and 11 samples with WARP's and the default arguments, the
whole reblocked HG002 chr20-22 (three arrays, 309K calls, 12 s vs GATK 44 s), an array starting mid-contig queried
mid-array and whole. Not supported (refused): several IDs at a site (GenomicsDB's order is not deterministic), PL
remapping for ploidy > 2, the `mean` / `concatenate` combine operations, `-V gendb://` without `-L`.

## ReblockGVCF (S3 group 4, `src/reblock.rs`, `tests/harness/reblock/`)

GATK `ReblockGVCF` on a single-sample, diploid HaplotypeCaller GVCF, with its defaults or WARP's arguments
(`GermlineVariantDiscovery.wdl` Reblock: `-do-qual-approx --floor-blocks -GQB 20 -GQB 30 -GQB 40`); options that need
the genotyping engine (`--drop-low-quals`), dbSNP, kept / removed annotations, tree scores or site filters are refused.
A line-by-line port of `regenotypeVC` and its writer: reference blocks passed on; a variant whose most likely genotype
has no concrete alt becomes a hom-ref block (GQ0 when PL[0] != 0, else PLs subset to ref + the best alt by
`calculateLikelihoodSums`); a called variant loses its uncalled alts (`subsetAlleles`, genotype from the subset PLs),
is reverse-trimmed with a reference block for the trimmed bases, gets its NON_REF AD zeroed, its INFO copied (the
StandardAnnotation / AS_StandardAnnotation keys; AS raw values subset, the NON_REF entry set to the empty raw value)
plus QUALapprox, VarDP, AS_QUALapprox, AS_VarDP, RAW_GT_COUNT; `ReblockingGVCFBlockCombiner`'s buffer trims and splits
reference blocks around variants; `GVCFBlockCombiner` / `HomRefBlock` merge them into the new GQ bands (floor blocks:
GQ = the band's lower bound, median DP, no PL). Header: the input's without GVCFBlock / MLEAC / MLEAF / HaplotypeScore /
DS lines, plus the lines ReblockGVCF adds (`tables/reblock_header_lines.txt`, taken from GATK's output) and the new
bands. GATK details reproduced: GQ is written capped at 99 (`IntGenotypeFieldAccessors`), not computed capped; the
trimmed-bases block's PLs come from an in-place loop (`L[i] - L[0]` with `L[0]` already zeroed); VarDP / AS_VarDP see
the zeroed NON_REF AD (the AD array is shared). `tests/harness/reblock/cmp.sh` (truth: `truth.sh`): byte-identical to GATK on
WARP's command (HG002 chr8:1-15M and all of shard023: 8.2M records -> 967,421 in 9 s), a GVCF without AS
annotations, one with dbSNP IDs, and the default arguments (no floor blocks); `wgs_cmp.sh` (truth: `wgs_truth.sh`):
the whole-genome HG002 GVCF (47 shards, 451M records -> 49.9M), all identical. `VARIANT_RS_REBLOCK_STATS=1` prints which code paths a run
took; never taken on these inputs: re-calling a genotype from its PLs (HaplotypeCaller's genotypes agree with them) and
the buffer's drop rules for overlapping deletions.

## SortSam (S3 group 4, `src/sortsam.rs`, `tests/harness/sortsam/`)

Picard `SortSam SORT_ORDER=coordinate` (other orders refused). htsjdk's sorting writer is a stable sort by
`SAMRecordCoordinateComparator` (TimSort chunks of `MAX_RECORDS_IN_RAM`, spill files merged with ties to the earlier
file), so any stable sort gives its order and `MAX_RECORDS_IN_RAM` is ignored. variant-rs sorts chunks of a memory
budget (`--sort-buffer-mb`, default 4096) in parallel, spills them as BGZF level 1 to `TMP_DIR` (default: the output's
directory), merges, and deflates the output in parallel; records keep their bytes (htsjdk re-encodes only records
whose attributes were set), the header is htsjdk's rewrite with `SO:coordinate`, `CREATE_INDEX` / `CREATE_MD5_FILE`
as in MarkDuplicates. `tests/harness/sortsam/cmp.sh` (truth: `truth.sh`): 17.7M records byte for byte, headers identical, the
`.bai` answers identically and the `.md5` is right, on WARP's command (BamProcessing.wdl) after WARP's MarkDuplicates
(query-grouped), a queryname-sorted input, an already sorted one and Picard spilling every 100,000 records; also with
variant-rs spilling 113 files. 13-14 s vs Picard 62-75 s.

## Limits / not done

* Diploid, default arguments apart from those listed under "HaplotypeCaller arguments beyond the defaults";
  several samples in VCF mode only (as GATK); no `--sample-name`; no BQSR `BI`/`BD` indel qualities.
* CRAM: needs the reference given with `-R` (no `REF_PATH`/EBI lookup is configured); multi-reference slices are
  read correctly but ignored by the `.crai` cost model.
* One machine: no gather step across machines, no multi-GPU (`variant-rs plan` shards can be run as separate
  processes, one per GPU/node, and their per-shard outputs concatenated in shard order).
* ~~CUDA PairHMM numerics~~ fixed 2026-11-12: the kernel is now bit-exact; whole-genome VCF and GVCF are
  byte-identical to GATK on all 23 shards with `--pairhmm cuda`.
* Single-run mode over many contigs with drawing groups is close to sequential (each contig waits at its first
  drawing group for the previous contigs); per-shard mode is the fast one.
* Errors: a panic on any thread (GPU, reader, decoder, back end) ends the run with exit code 101, before any
  output is finalized. That includes an assembly failure (where GATK throws, e.g. "Graph must have ref source and
  sink vertices"): GATK rethrows it and fails; until 2026-11-17 this port skipped the region and wrote a
  complete-looking VCF without its calls (`variant-rs/tests/assembly_fail.sh`, hook `VARIANT_RS_TEST_ASSEMBLY_FAIL`). A BAM that ends before the blocks its index points to is an error
  ("BAM file truncated"); before 2026-11-12 the reader stopped there silently and a crashed GPU thread left the
  run hanging.
* `--dump-activity` needs a single contig; the old `--dump-region-reads/--dump-final-reads/--dump-haplotypes` were dropped.

## Files

* `src/main.rs`: CLI, domains/jobs, scheduling (LPT), result ordering, output.
* `src/plan.rs`: `variant-rs plan` (BAI cost model, contig splitting, shard files) and the scheduler's cost model.
* `src/frontend.rs`: per-contig front end (decode thread, downsampling + RNG schedule, pileups, regions, filling).
* `src/active.rs`: read filters, isActive accumulators, activity profile (ported 1:1).
* `src/backend.rs`: batches through preproc/assembly/likelihood/genotype/refconf, QD placeholder bookkeeping.
* `src/rng.rs`: the per-domain generator baton, QD placeholder resolution.
* `src/filter.rs`: hard filters.
* `src/markdup/`: MarkDuplicates (`mod.rs`, `metrics.rs`); `src/bamwrite.rs`: BAM output shared with applybqsr.
* `src/bqsr/`: BQSR (`mod.rs` driver and record decoding, `recal.rs`, `covariates.rs`, `known.rs`, `datum.rs`,
  `report.rs`, `apply.rs`).
* `src/input.rs`: BAM or CRAM input; `src/cram.rs` + `csrc/cram_shim.c`: htslib CRAM reader and CRAI cost model.
* `src/bam.rs`, `src/bgzf.rs`, `src/tabix.rs`, `src/reference.rs`, `src/intervals.rs`, `src/output.rs`: I/O.
* `tests/harness/`: the comparisons with GATK (`regress.sh` runs them all).
