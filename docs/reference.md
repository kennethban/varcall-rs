# Reference

This page documents the commands, options and environment variables of align-rs and variant-rs.
[getting-started.md](getting-started.md) gives an introduction, and [how-to.md](how-to.md) gives procedures for
common tasks.

* [align-rs](#align-rs): [environment variables](#align-rs-environment-variables)
* [variant-rs](#variant-rs): [GATK tools](#gatk-command-lines), [native options](#native-options-and-subcommands),
  [environment variables](#variant-rs-environment-variables), [exit status](#exit-status)

## align-rs

```bash
align-rs mem [-t N] [-K N] [-Y] [-R RG] [-o out.sam | -o out.bam [-l level]] [--gpu] ref.fa R1.fq[.gz] [R2.fq[.gz]]
align-rs mem [-t N] [-K N] [-Y] [-l level] [--gpu] --lanes LANES.tsv ref.fa
align-rs tune [-t N] [-K N] ref.fa R1.fq[.gz] [R2.fq[.gz]]
```

* **`mem`** aligns reads as `bwa-mem2 mem` does and produces the same records. It accepts four options of
  bwa-mem2 (`-t`, `-K`, `-Y` and `-R`). Every other option remains at the default of bwa-mem2, and specifying one
  stops the run ([unsupported features](#unsupported-features)). The reference index is the one built by
  `bwa-mem2 index ref.fa`.
  * `-o x.bam` writes a coordinate-sorted BAM and its `.bai`. The records are in the order produced by
    `samtools sort`, and the index is byte-identical to that of `samtools index`. `-l` sets the compression level
    (default 1).
  * `-o x.sam`, or the absence of `-o`, writes the SAM output of bwa-mem2.
  * `--gpu` runs seeding, SA lookups, extension and mate rescue on the GPU, or on the CPU when no GPU is present.
* **`-K`** is the batch size of bwa-mem2 in bases. The output of bwa-mem2 depends on this value, because
  insert-size statistics are estimated for each batch, and does not depend on the number of threads. The value
  should equal that of the run being reproduced. The pipelines of GATK and nf-core use `-K 100000000`.
* **`--lanes FILE`** aligns several lanes of a sample in one process, which loads the index once. The file has one
  line per lane, of the form `out.bam<TAB>@RG\tID:...<TAB>R1.fq.gz[<TAB>R2.fq.gz]`. The BAM file of each lane is
  byte-identical to the file obtained by aligning that lane on its own. Two lanes cannot share an output file.
* **`tune`** measures the GPU settings listed below on the GPU at hand, which takes about 35 s in addition to
  loading the index, and saves them for `mem`. None of these settings changes the output.

### align-rs environment variables

| environment | |
|---|---|
| `ALIGN_RS_SMEM_G` | seeding: GPU lanes per read (1, 2, 4, 8, 16, 32; default 4) |
| `ALIGN_RS_EXT_BPSM` | extension: blocks of 128 threads per SM (default from the L1 size) |
| `ALIGN_RS_MSW=striped` | mate SW: the striped kernel (16 threads per job emulating bwa's SSE lanes) instead of the wavefront one; same results, slower |
| `ALIGN_RS_UNGAPPED` | extension's ungapped shortcut: `0` off (every job runs the banded DP; same results), `check` runs both and stops on a difference (default on) |
| `ALIGN_RS_CIGAR_BAND` | the CIGAR's global alignment in the band certified by the ungapped score (equal lengths): `0` the usual band (same results), `check` runs both and stops on a difference (default on) |
| `ALIGN_RS_MSW_MINB` | mate SW, striped kernel only: blocks per SM the registers are capped for (default from the thread slots) |
| `ALIGN_RS_TUNE_CACHE` | where `tune` saves (default `$XDG_CACHE_HOME/align-rs/gpu-settings.tsv` or `~/.cache/align-rs/...`) |
| `ALIGN_RS_PROFILE=1` / `2` | per-stage timeline and counters at the end / plus per-item CPU timers |
| `ALIGN_RS_GPU_SERIAL=1` | one GPU call at a time (kernel-only timings) |
| `ALIGN_RS_TRACE=FILE` | a line per stage run, extension device call (upload / sort / kernel / download) and mate-SW chunk (CPU preparation, upload / kernel / download), on one time origin (not free: a traced 30x lane runs ~10% slower) |
| `ALIGN_RS_MEMLOG=SECS` | every SECS seconds, the process's memory (anonymous, file-backed, swap) and how many batches wait to be dropped or finished |
| `ALIGN_RS_EXT_POOL=0` | extension's CPU steps in the global thread pool instead of their own (for comparison) |
| `ALIGN_RS_REPLAY_CHECK=1` | runs bwa's linear scan beside the optimized replay lookup and stops if they disagree |
| `ALIGN_RS_GPU_FAIL_GROW=N` | test hook: device allocations after the first N fail (exercises the CPU fallback) |

The remaining variables (`ALIGN_RS_DUMP_*`, `ALIGN_RS_KSW_EMULATE`, `ALIGN_RS_SLACK*`, `ALIGN_RS_MI_*`,
`ALIGN_RS_PTRACER`) are switches for debugging and measurement. They are documented in `align-rs/src` at the points
where they are read.

## variant-rs

### GATK command lines

```bash
variant-rs <Tool> <GATK arguments>         # as `gatk <Tool> ...`
variant-rs <Tool> --help                   # which arguments are mapped, ignored, or accepted at the default only
```

| Tool | Notes |
|---|---|
| HaplotypeCaller | VCF and `-ERC GVCF`; one sample, or several jointly (VCF mode); WARP's arguments (`-G` groups incl. AS, `-GQB`, `-contamination`, `--dbsnp`, interval padding / exclusion); on-the-fly BQSR (`--bqsr-recal-file`) |
| Mutect2 | tumor / normal and tumor-only; PoN, germline resource, `--f1r2-tar-gz`, `--genotype-germline-sites`, `--genotype-pon-sites`, `--mitochondria-mode`, `--bam-output` |
| CreateSomaticPanelOfNormals | from a GenomicsDB workspace of the normals' Mutect2 calls, as GATK's `mutect2_pon.wdl` runs it: `--germline-resource`, `--min-sample-count`, `--max-germline-probability` |
| FilterMutectCalls, LearnReadOrientationModel, GetPileupSummaries, GatherPileupSummaries, CalculateContamination, MergeMutectStats | as GATK's somatic Best Practices and sarek run them |
| BaseRecalibrator, ApplyBQSR, GatherBQSRReports | default covariates; WARP's `-OQ`, `--static-quantized-quals` |
| MarkDuplicates | coordinate or query-grouped input, several inputs, WARP's and sarek's arguments, Picard `KEY=value` syntax |
| GenotypeGVCFs | a GVCF file or a GenomicsDB workspace (`-V gendb://`); default or WARP annotations; `--scatter-intervals`: each `-L` interval its own GATK run (sarek's scatter), side by side |
| GenomicsDBImport | `--sample-name-map` or `-V` per GVCF, `-L` (strings or files), `--merge-input-intervals`, `--batch-size`, `--consolidate`, `--bypass-feature-reader` (WARP's and sarek's commands): the same cells, `vidmap.json`, `callset.json` and header as GATK's workspace; the intervals imported side by side (`--threads`) |
| VariantFiltration, MakeSitesOnlyVcf | INFO-level filter expressions (GATK's JEXL: comparisons, `&&` / `||` / `!`, strings, flags, `QUAL`, `TYPE`, `CHROM`, `POS`, `N_ALLELES`, `FILTER`), `--missing-values-evaluate-as-failing`, `--invert-filter-expression`, `--invalidate-previous-filters`, masks (BED, interval list, VCF); WARP's joint-calling and DRAGEN hard filters |
| VariantRecalibrator, ApplyVQSR | SNP, INDEL and BOTH modes; WARP's and sarek's arguments (`--resource:name,known=...,training=...,truth=...,prior=...`, tranches, `--max-gaussians`, `--trust-all-polymorphic`) |
| ReblockGVCF, MergeVcfs, SortSam (coordinate), GatherBamFiles | |

Every GATK argument of a tool belongs to one of three groups:

* **mapped:** variant-rs implements the argument;
* **ignored:** the argument cannot change the output (logging, temporary directories, `--java-options`,
  `--native-pair-hmm-threads` and similar);
* **default only:** the argument is accepted only at the default value of GATK.

Any other value, any unknown argument and any input that variant-rs cannot reproduce is refused with an
explanation and exit status 2. The argument tables are in `variant-rs/gatkargs/` and were generated from the
`--help` output of GATK 4.7. The three groups are listed for every tool in [gatk-arguments.md](gatk-arguments.md).

* `GATK_FALLBACK=/path/to/gatk` causes a refused command to be run by GATK instead, with its original arguments.
  The variant-rs image uses this mechanism to replace `gatk`
  ([containers/variant-rs](../containers/variant-rs/README.md)).
* `VARIANT_RS_DRY_RUN=1` prints the command that would run (variant-rs with its translated options, or GATK with
  the reason) and exits.
* variant-rs uses all CPUs that the process is permitted to use, including cgroup quotas. The outputs do not depend
  on the number of threads.

### Native options and subcommands

The GATK command lines are translated into the native options of variant-rs, which can also be given directly.
Two functions are available only in this form:

```bash
# HaplotypeCaller over many shard files in one process, each equal to its own GATK run
variant-rs -I x.bam -R ref.fa --shard s/shard00.intervals ... --shard s/shard22.intervals -O all.vcf.gz [-ERC GVCF]
# cut a genome into load-balanced shard files from the BAM index alone
variant-rs plan -I x.bam -R ref.fa -o plan96 --pieces 96
```

* **`--pairhmm cuda|cpu|auto`** selects the PairHMM implementation. Both implementations give the results of GATK.
* **`--threads`, `--contig-workers`** set the degree of parallelism. The default number of contig workers is half
  the number of threads up to 8, and a quarter of the number of threads beyond that, up to 32.
* **`-L file --scatter-intervals`** calls every interval as a separate GATK run, as sarek scatters HaplotypeCaller.

[notes/variant-rs.md](notes/variant-rs.md) contains the complete option table and describes the relation between
shards and the random-number streams of GATK.

### variant-rs environment variables

| environment | |
|---|---|
| `VARIANT_RS_THREADS` | Mutect2 / contamination: worker threads (default all CPUs) |
| `VARIANT_RS_CONTIG_THREADS`, `VARIANT_RS_DECODE_THREADS`, `VARIANT_RS_TWO_STAGE` | Mutect2 / HaplotypeCaller front-end sizing |
| `VARIANT_RS_GPU_STREAMS` | GPU PairHMM: engines (CUDA streams with their own buffers) per GPU, shared by the callers (default 8; 1 on one GPU: one call at a time) |
| `VARIANT_RS_GPUS` | GPU PairHMM: use at most this many of the visible GPUs (default all of `CUDA_VISIBLE_DEVICES`) |
| `VARIANT_RS_HMM_KERNEL=warp`, `VARIANT_RS_HMM_BINTIME` | GPU PairHMM: the one-warp-per-pair float kernel for every pair (same output, slower); kernel time per read-length bin at exit |
| `VARIANT_RS_TRACE_STAGES`, `VARIANT_RS_TRACE_CPU`, `VARIANT_RS_TRACE_FRONT`, `VARIANT_RS_TRACE_REGIONS` | timing and tracing output |
| `VARIANT_RS_ASM_PROFILE` (with `VARIANT_RS_TRACE_CPU`), `VARIANT_RS_REBLOCK_STATS` | assembly phase timings and code-path counters; ReblockGVCF code-path counters |
| `VARIANT_RS_NO_AVX512` | CPU PairHMM: take the path for CPUs without AVX-512 (same output) |
| `VARIANT_RS_TEST_ASSEMBLY_FAIL` | test hook: makes an assembly fail (the run must fail, as GATK's does) |

### Exit status

| status | |
|---|---|
| 0 | success |
| 2 | refused: a tool, argument value or input variant-rs does not reproduce exactly (the message says which); with `GATK_FALLBACK` set, GATK runs instead |
| other | an error: unreadable or malformed input, a failed assembly (as in GATK), ... (101 for errors raised as Rust panics); no output is finalized |

## Unsupported features

Each item in this section is refused with an explanation and is never run with a different result. A refused
variant-rs command exits with status 2, or is run by GATK when `GATK_FALLBACK` is set. Only the reference versions
are reproduced: bwa-mem2 2.2.1 (and 2.3, which gives the same alignments), GATK 4.7.0.0, Picard 3.5.0 and htsjdk
5.0.0.

### align-rs

* **Other options of `bwa-mem2 mem`.** align-rs accepts `-t`, `-K`, `-R`, `-Y` and its own output options. Any
  other option (the scoring options `-A -B -O -E -L -U`, `-k`, `-w`, `-r`, `-c`, `-T`, `-M`, `-a`, `-C`, `-p`, `-I`,
  the `-x` presets and others) stops the run with the message "unsupported option". Without these options, align-rs
  uses the defaults of bwa-mem2, as the pipelines of GATK and nf-core do.
* **Other commands of bwa-mem2.** The index must be built with `bwa-mem2 index`.
* **Reads longer than 65,535 bases.** The run stops and reports the name of the read.

### variant-rs: unsupported GATK tools

variant-rs provides the 23 tools listed in the [table above](#gatk-command-lines). The following tools are among
those that it does not provide:

* **Model-based filtering other than VQSR:** CNNScoreVariants / NVScoreVariants, FilterVariantTranches, VETS
  (ExtractVariantAnnotations, TrainVariantAnnotationsModel, ScoreVariantAnnotations).
* **SelectVariants.** The command `variant-rs filter` (the HARDFILTER step of the pipeline) gives, in one pass, the
  result of the single-sample hard filtering of the GATK Best Practices (SelectVariants → VariantFiltration →
  MergeVcfs).
* **RNA-seq** (SplitNCigarReads), **annotation** (Funcotator), **copy number and structural variants**, and
  **Picard's metrics** tools (CollectWgsMetrics, ...).

### variant-rs: unsupported modes of supported tools

Each tool supports the default arguments of GATK and the arguments that WARP and nf-core/sarek pass. The list for
every tool is in [gatk-arguments.md](gatk-arguments.md). The following modes are among those that are refused:

* **HaplotypeCaller:** ploidy other than 2 (`-ploidy`); `--dragen-mode`; flow-based calling (`--flow-mode`);
  `--linked-de-bruijn-graph`; force-calling given alleles (`--alleles`); a non-default `--pcr-indel-model` or
  `--min-base-quality-score`; extra read filters (`-RF`) or single annotations (`-A`; annotation groups `-G` are
  supported); `-bamout` (Mutect2 supports it); the same read in two inputs.
* **BaseRecalibrator:** covariates other than the default ones (`-cov`).
* **GenomicsDBImport:** adding samples to a workspace (`--genomicsdb-update-workspace-path`),
  `--merge-contigs-into-num-partitions`, `--header`, `--avoid-nio`, sample maps with index files; GVCFs with several
  samples or with FILTER values (GATK's have none); `-V` headers that disagree on a field's Number or Type (GATK
  promotes the line).
* **CreateSomaticPanelOfNormals:** intervals (`-L`: the workspace is read whole), a VCF as input instead of a
  workspace, `--call-genotypes`, `--genomicsdb-use-bcf-codec`.
* **VariantFiltration:** genotype filters (`--genotype-filter-expression`, `--set-filtered-genotype-to-no-call`), SNP
  clusters (`--cluster-window-size`), `--apply-allele-specific-filters`, intervals (`-L`); in expressions, arithmetic,
  method calls (`vc.isSNP()`), regular expressions, the `vc` / `g` / `ALLELES` / genotype-count variables, and
  decimal literals with more than 6 significant digits. An expression GATK itself cannot evaluate (an integer
  literal against a decimal value, as in `QD > 2`; a non-boolean value) fails the run, as in GATK.
* **VCF inputs whose samples are not in sorted order** (VariantFiltration, ApplyVQSR, MergeVcfs): htsjdk then
  re-encodes every genotype, which is not reproduced. GATK writes its VCFs with sorted samples.
* **VariantRecalibrator / ApplyVQSR:** allele-specific mode (`-AS`), `--aggregate` inputs, saving or reading the
  model (`--output-model`, `--input-model`), scattered tranches (`--output-tranches-for-scatter`), the R plots
  (`--rscript-file`), intervals (`-L`) and more than one `-V`.
* **SortSam:** sort orders other than `coordinate`.
