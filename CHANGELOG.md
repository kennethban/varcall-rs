# Changelog

All notable changes are listed here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and versions follow [Semantic Versioning](https://semver.org/). A version's outputs are tied to the reference tools
named in its entry; a change of reference version is a new minor version at least.

## [Unreleased]

### Added

* align-rs 0.1.0: `bwa-mem2 mem` 2.2.1 with seeding, SA lookups, extension and mate rescue on the GPU; sorted BAM and
  `.bai` output identical to `samtools sort` / `samtools index` (1.25); several lanes per process (`--lanes`); `tune`.
* variant-rs 0.1.0: GATK 4.7.0.0 / Picard 3.5.0 HaplotypeCaller (VCF, GVCF), Mutect2, FilterMutectCalls,
  LearnReadOrientationModel, GetPileupSummaries, GatherPileupSummaries, CalculateContamination, MergeMutectStats,
  BaseRecalibrator, ApplyBQSR, GatherBQSRReports, MarkDuplicates, GenotypeGVCFs (GVCF and GenomicsDB input),
  ReblockGVCF, MergeVcfs, SortSam, GatherBamFiles; CUDA PairHMM; GATK command lines with fallback to GATK.
* A Nextflow FASTQ → VCF pipeline; Docker images for align-rs and for variant-rs as a drop-in `gatk` in nf-core/sarek.
* Tests: unit tests, a self-contained test on synthetic data (CI), and a full harness on public data.
* Pipeline: several GPUs. With `CUDA_VISIBLE_DEVICES` listing them at launch, each GPU step gets its own GPU, and a
  sample's lanes are split over one align-rs task per GPU.
* Pipeline: `pipeline/run.sh` launches it on the idle GPUs (`VARCALL_GPUS` for at most N), unless
  `CUDA_VISIBLE_DEVICES` is set. Threads by default from the machine and the GPUs: each alignment task the machine's
  threads / GPUs (at most 64, `align_cpus`), each other heavy step the machine's (at most 128, `cpus`; was 16).
* variant-rs: the GPU PairHMM serves its callers from a pool of engines (`VARIANT_RS_GPU_STREAMS`, default 8), and
  HaplotypeCaller runs more contig workers beyond 32 threads: a 30x genome on an L40S at 64 threads 190 -> 122 s.
* variant-rs: the GPU PairHMM's engines spread over every visible GPU (`VARIANT_RS_GPUS` for at most N); the pipeline
  gives HaplotypeCaller up to 2 GPUs (`hc_gpus`): a 30x genome on L40S at 128 threads 115 -> 98 s.
* variant-rs: VariantRecalibrator and ApplyVQSR (GATK 4.7.0.0), identical to GATK on a 30x genome and a
  1000 Genomes cohort with WARP's, sarek's and best-practice settings; a 30x genome's SNP model in ~1 min (GATK ~11).
* variant-rs: VariantFiltration (filter expressions, masks) and Picard MakeSitesOnlyVcf, identical to GATK on WARP's
  joint-calling, exome-mask and DRAGEN filters, GATK's hard-filtering examples and the expression language's edge
  cases.
* variant-rs: GenomicsDBImport (GenomicsDB 1.5.5 workspaces): the same cells, `vidmap.json`, `callset.json`, header
  and arrays as GATK's workspace, which GATK's GenotypeGVCFs reads too; WARP's and sarek's commands (sample map or
  `-V`, intervals, batches, `--consolidate`, `--bypass-feature-reader`).
* variant-rs: CreateSomaticPanelOfNormals, with GenomicsDBImport of Mutect2 VCFs (flag INFO fields): GATK's
  `mutect2_pon.wdl` chain identical to GATK on the DREAM normals.
* Pipeline: joint germline calling (`joint_germline`, as nf-core/sarek's): GVCFs, GenomicsDBImport, GenotypeGVCFs and
  VQSR. GenotypeGVCFs genotypes each interval as its own GATK run (`scatter_intervals`), as sarek's per-interval tasks
  do; the VQSR recal files are published with their indexes.
* variant-rs GenomicsDBImport: the intervals (arrays) imported side by side (`--threads`, default all CPUs); the
  workspace's cells are unchanged. A whole-genome trio (356 intervals): ~2.5 h -> 6.5 min.
* variant-rs GenotypeGVCFs `--scatter-intervals` (with `--threads`): each `-L` interval genotyped as its own GATK run,
  side by side, records in interval order: identical to GATK run once per interval and concatenated (nf-core/sarek's
  scatter). GATK's single run differs from that (QD's random jitter is one stream across all intervals), so without
  the flag variant-rs stays one run.
* align-rs: the sorted BAM's final pass in parallel (a 30x lane 25 -> 5 s); extension on its own thread pool.
* variant-rs MarkDuplicates: smaller on-disk read ends (buckets 60 -> 43 GB), read straight into their arrays, and a
  cheaper first pass: a 30x genome on 32 threads 281 -> 255 s. `VARIANT_RS_TRACE_CPU=1` prints its CPU by step
  (and slows it). On the server at 128 threads, with level 2: 57 -> 48 s.
* variant-rs HaplotypeCaller: each region's work before the PairHMM, and after it, runs as one task (two parallel
  passes over a batch instead of five): a 30x genome on 32 threads 456 -> 400 s. The stage-time line reports
  `finalize + assembly + prep` and `finish + realign` together.
* Pipeline: MarkDuplicates writes its BAM at compression level 2, as nf-core/sarek's (the gatk launcher's
  `-Dsamjdk.compression_level=2`): the same records, the BAM ~7% larger, the step 255 -> 246 s.
* align-rs: `ALIGN_RS_MEMLOG=SECS` prints the process's memory and its queues every SECS seconds.
* variant-rs CUDA PairHMM: a tiled float kernel (several threads per pair, read rows in registers, pairs binned by
  read length) and the double-precision redo spread over the SMs; results still bit-identical to GKL AVX-512. Kernel
  282 -> 663 GCUPS on an L40S; HaplotypeCaller on a 30x genome, 1 L40S at 128 threads 106 -> 94 s, the workstation
  400 -> 355 s. `VARIANT_RS_HMM_KERNEL=warp` selects the old one-warp-per-pair kernel; `VARIANT_RS_HMM_BINTIME=1`
  prints the kernel time per read-length bin.
* align-rs: mate SW as a wavefront over the query's 16 segments (no lazy-F loop, no per-row vote), the same results
  as bwa-mem2's ksw_u8: kernel 21 -> 12 ms per 65,536 jobs on an L40S, 31 -> 18 on an RTX 5070 Ti; a 30x lane on
  one L40S 86.5 -> 82.4 s. `ALIGN_RS_MSW=striped` runs the previous kernel (`ALIGN_RS_MSW_MINB` and `tune` now
  apply to it only).
* align-rs: extension's ungapped shortcut (from bwa-mem3): a job whose seed diagonal loses less than a gap open
  (default scoring: at most one mismatch) takes its result from a walk along the diagonal on the GPU instead of the
  banded DP, identical results; a quarter of the jobs, the extension kernel -12% on an L40S. `ALIGN_RS_UNGAPPED=0`
  turns it off, `=check` compares it with the DP.
* align-rs: the CIGAR's global alignment runs equal-length regions in the narrower band the ungapped score certifies
  (from bwa-mem3), identical score, CIGAR, NM and MD: half the band cells, CIGAR CPU -31% on 5 M pairs.
  `ALIGN_RS_CIGAR_BAND=0` turns it off, `=check` compares it with the usual band.
* align-rs: mate rescue's region dedup (`mem_sort_dedup_patch` after each of the four directions) adds the rescued
  region to the already deduplicated regions without sorting when they are clean (distinct ends, nothing the overlap
  loop would drop, distinct keys), and skips the call when nothing was rescued; region sorts use one packed integer
  key and permute in place. Same records; a 30x lane on the workstation 149.5 -> 138.8 s, its CPU time -21%.
  `ALIGN_RS_DEDUP_CHECK=1` / `ALIGN_RS_SORT_CHECK=1` compare with the full function / the previous sort.

### Fixed

* align-rs: memory grew without bound on long runs (a workstation's 12 lanes filled 48 GB in ~6 minutes and took
  WSL down): the thread freeing each batch's structures fell behind. Four threads now share a bounded queue. On
  the 6-GPU server the alignment tasks drop from 54-169 GB to 37 GB each and the germline pipeline from 460 to 409 s.
* variant-rs HaplotypeCaller `-ERC GVCF`: a site whose most likely genotype includes `<NON_REF>` is written `0/0` with
  all PLs 0, as GATK's `makeGenotypeCall` does (its GQ stays). Once in three 30x genomes (M3: NA12892,
  chr10:128583486, a tandem repeat with no informative reads).

* variant-rs GenotypeGVCFs: the input records' IDs (a GVCF's rsIDs from HaplotypeCaller `--dbsnp`) are dropped as
  GATK drops them; the ID is `--dbsnp`'s matches only. Seen in the pipeline's joint calling test.

* variant-rs ApplyVQSR, MergeVcfs: INFO decoded and written as htsjdk does (`src/vcfinfo.rs`): an empty value
  (`KEY=`) is written `KEY=.`, a declared flag given as `FLAG=0` is dropped, and a key the header does not declare
  fails the run, as in GATK and Picard. Inputs whose samples are not in sorted order are refused (htsjdk re-encodes
  their genotypes; ApplyVQSR reordered them as read).

* variant-rs MergeVcfs: standard INFO / FORMAT header lines with a non-standard Number or Type (e.g. `AD` with
  `Number=.`, as in the 1000 Genomes calls) are repaired as htsjdk does on reading, as Picard's output has them.

[Unreleased]: https://github.com/kennethban/varcall-rs/commits/main
