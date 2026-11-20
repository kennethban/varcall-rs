# Reference

Commands, options and environment variables of align-rs and variant-rs. For a walk-through see
[getting-started.md](getting-started.md); for recipes, [how-to.md](how-to.md).

* [align-rs](#align-rs): [environment variables](#align-rs-environment-variables)
* [variant-rs](#variant-rs): [GATK tools](#gatk-command-lines), [native options](#native-options-and-subcommands),
  [environment variables](#variant-rs-environment-variables), [exit status](#exit-status)

## align-rs

```bash
align-rs mem [-t N] [-K N] [-Y] [-R RG] [-o out.sam | -o out.bam [-l level]] [--gpu] ref.fa R1.fq[.gz] [R2.fq[.gz]]
align-rs mem [-t N] [-K N] [-Y] [-l level] [--gpu] --lanes LANES.tsv ref.fa
align-rs tune [-t N] [-K N] ref.fa R1.fq[.gz] [R2.fq[.gz]]
```

* **`mem`** takes these of bwa-mem2 `mem`'s options (the others stay at bwa-mem2's defaults: see
  [what it does not cover](#what-it-does-not-cover)) and gives its records. `-o x.bam` writes a coordinate-sorted BAM and its
  `.bai`: records as `samtools sort` orders them, the index byte-identical to `samtools index`. `-o x.sam` (or no
  `-o`) writes bwa-mem2's SAM. `--gpu` runs seeding, SA lookups, extension and mate rescue on the GPU (the CPU when
  there is none). The index is bwa-mem2's (`bwa-mem2 index ref.fa`).
* **`-K`** is bwa-mem2's batch size in bases. bwa-mem2's output depends on it (insert-size statistics are per batch),
  not on the thread count. Use the value the reference run used. GATK's and nf-core's pipelines use `-K 100000000`.
* **`--lanes FILE`**: several lanes of a sample in one process (the index loads once). One line per lane:
  `out.bam<TAB>@RG\tID:...<TAB>R1.fq.gz[<TAB>R2.fq.gz]`. Each lane's BAM is byte-identical to aligning that lane on
  its own. Two lanes cannot share an output file.
* **`tune`** measures the GPU settings below on the GPU at hand (~35 s plus the index load) and saves them for
  `mem`. No setting changes the output.

### align-rs environment variables

| environment | |
|---|---|
| `ALIGN_RS_SMEM_G` | seeding: GPU lanes per read (1, 2, 4, 8, 16, 32; default 4) |
| `ALIGN_RS_EXT_BPSM` | extension: blocks of 128 threads per SM (default from the L1 size) |
| `ALIGN_RS_MSW_MINB` | mate SW: blocks per SM the registers are capped for (default from the thread slots) |
| `ALIGN_RS_TUNE_CACHE` | where `tune` saves (default `$XDG_CACHE_HOME/align-rs/gpu-settings.tsv` or `~/.cache/align-rs/...`) |
| `ALIGN_RS_PROFILE=1` / `2` | per-stage timeline and counters at the end / plus per-item CPU timers |
| `ALIGN_RS_GPU_SERIAL=1` | one GPU call at a time (kernel-only timings) |
| `ALIGN_RS_TRACE=FILE` | a line per stage run, extension device call (upload / sort / kernel / download) and mate-SW chunk (CPU preparation, upload / kernel / download), on one time origin (not free: a traced 30x lane runs ~10% slower) |
| `ALIGN_RS_EXT_POOL=0` | extension's CPU steps in the global thread pool instead of their own (for comparison) |
| `ALIGN_RS_REPLAY_CHECK=1` | runs bwa's linear scan beside the optimized replay lookup and stops if they disagree |
| `ALIGN_RS_GPU_FAIL_GROW=N` | test hook: device allocations after the first N fail (exercises the CPU fallback) |

The rest (`ALIGN_RS_DUMP_*`, `ALIGN_RS_KSW_EMULATE`, `ALIGN_RS_SLACK*`, `ALIGN_RS_MI_*`, `ALIGN_RS_PTRACER`) are
debugging and measurement switches, described where they are read in `align-rs/src`.

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
| FilterMutectCalls, LearnReadOrientationModel, GetPileupSummaries, GatherPileupSummaries, CalculateContamination, MergeMutectStats | as GATK's somatic Best Practices and sarek run them |
| BaseRecalibrator, ApplyBQSR, GatherBQSRReports | default covariates; WARP's `-OQ`, `--static-quantized-quals` |
| MarkDuplicates | coordinate or query-grouped input, several inputs, WARP's and sarek's arguments, Picard `KEY=value` syntax |
| GenotypeGVCFs | a GVCF file or a GenomicsDB workspace (`-V gendb://`); default or WARP annotations |
| GenomicsDBImport | `--sample-name-map` or `-V` per GVCF, `-L` (strings or files), `--merge-input-intervals`, `--batch-size`, `--consolidate`, `--bypass-feature-reader` (WARP's and sarek's commands): the same cells, `vidmap.json`, `callset.json` and header as GATK's workspace |
| VariantFiltration, MakeSitesOnlyVcf | INFO-level filter expressions (GATK's JEXL: comparisons, `&&` / `||` / `!`, strings, flags, `QUAL`, `TYPE`, `CHROM`, `POS`, `N_ALLELES`, `FILTER`), `--missing-values-evaluate-as-failing`, `--invert-filter-expression`, `--invalidate-previous-filters`, masks (BED, interval list, VCF); WARP's joint-calling and DRAGEN hard filters |
| VariantRecalibrator, ApplyVQSR | SNP, INDEL and BOTH modes; WARP's and sarek's arguments (`--resource:name,known=...,training=...,truth=...,prior=...`, tranches, `--max-gaussians`, `--trust-all-polymorphic`) |
| ReblockGVCF, MergeVcfs, SortSam (coordinate), GatherBamFiles | |

Every GATK argument of a tool is in one of three groups:

* **mapped:** variant-rs implements it;
* **ignored:** it cannot change the output (logging, temporary directories, `--java-options`,
  `--native-pair-hmm-threads`, ...);
* **accepted at GATK's default value only.**

Everything else is refused with the reason, exit status 2: other values, unknown arguments, inputs variant-rs
cannot reproduce. The tables are in `variant-rs/gatkargs/`, generated from GATK 4.7's `--help`; the three groups of
every tool are listed in [gatk-arguments.md](gatk-arguments.md).

* `GATK_FALLBACK=/path/to/gatk`: a refused command runs the real GATK instead, with the original arguments. This is
  what makes the variant-rs image a drop-in `gatk` ([containers/variant-rs](../containers/variant-rs/README.md)).
* `VARIANT_RS_DRY_RUN=1`: prints what would run (variant-rs, with its translated options, or GATK and why) and exits.
* Threads: all CPUs the process may use (cgroup quotas included). Outputs do not depend on the thread count.

### Native options and subcommands

The GATK tools translate into variant-rs's own options, which can also be given directly. Two things exist only
there:

```bash
# HaplotypeCaller over many shard files in one process, each equal to its own GATK run
variant-rs -I x.bam -R ref.fa --shard s/shard00.intervals ... --shard s/shard22.intervals -O all.vcf.gz [-ERC GVCF]
# cut a genome into load-balanced shard files from the BAM index alone
variant-rs plan -I x.bam -R ref.fa -o plan96 --pieces 96
```

* **`--pairhmm cuda|cpu|auto`** chooses the PairHMM. Both give GATK's results.
* **`--threads`, `--contig-workers`** size the parallelism (contig workers by default: half the threads up to 8, a
  quarter of them beyond, up to 32).
* **`-L file --scatter-intervals`** calls every interval as its own GATK run, as sarek scatters HaplotypeCaller.

See [notes/variant-rs.md](notes/variant-rs.md) for the full option table and for how shards relate to GATK's
random-number streams.

### variant-rs environment variables

| environment | |
|---|---|
| `VARIANT_RS_THREADS` | Mutect2 / contamination: worker threads (default all CPUs) |
| `VARIANT_RS_CONTIG_THREADS`, `VARIANT_RS_DECODE_THREADS`, `VARIANT_RS_TWO_STAGE` | Mutect2 / HaplotypeCaller front-end sizing |
| `VARIANT_RS_GPU_STREAMS` | GPU PairHMM: engines (CUDA streams with their own buffers) per GPU, shared by the callers (default 8; 1 on one GPU: one call at a time) |
| `VARIANT_RS_GPUS` | GPU PairHMM: use at most this many of the visible GPUs (default all of `CUDA_VISIBLE_DEVICES`) |
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

## What it does not cover

Each item below is refused with the reason, never run differently. A refused variant-rs command exits with status
2, or runs the real GATK under `GATK_FALLBACK`. Only the reference versions are reproduced: bwa-mem2 2.2.1, GATK
4.7.0.0, Picard 3.5.0, htsjdk 5.0.0.

### align-rs

* **Other `bwa-mem2 mem` options.** align-rs takes `-t`, `-K`, `-R`, `-Y` and its own output options; any other
  (scoring `-A -B -O -E -L -U`, `-k`, `-w`, `-r`, `-c`, `-T`, `-M`, `-a`, `-C`, `-p`, `-I`, `-x` presets, ...)
  stops it with "unsupported option". Without them it runs bwa-mem2's defaults, as GATK's and nf-core's pipelines
  do.
* **Other bwa-mem2 commands.** Build the index with `bwa-mem2 index`.
* **Reads over 65,535 bases** (the run stops with the read's name).

### variant-rs: GATK tools it does not have

It has the 22 tools in the [table above](#gatk-command-lines). Notable tools it does not have:

* **Model-based filtering other than VQSR:** CNNScoreVariants / NVScoreVariants, FilterVariantTranches, VETS
  (ExtractVariantAnnotations, TrainVariantAnnotationsModel, ScoreVariantAnnotations).
* **SelectVariants.** `variant-rs filter` (the pipeline's HARDFILTER step) gives the result of GATK's Best Practices
  single-sample hard filtering (SelectVariants → VariantFiltration → MergeVcfs) in one pass.
* **CreateSomaticPanelOfNormals.** Mutect2 takes a panel of normals; making one needs GATK.
* **RNA-seq** (SplitNCigarReads), **annotation** (Funcotator), **copy number and structural variants**, and
  **Picard's metrics** tools (CollectWgsMetrics, ...).

### variant-rs: modes of the tools it has

Each tool covers GATK's defaults plus the arguments that WARP and nf-core/sarek pass; every tool's list is in
[gatk-arguments.md](gatk-arguments.md). Notable modes it refuses:

* **HaplotypeCaller:** ploidy other than 2 (`-ploidy`); `--dragen-mode`; flow-based calling (`--flow-mode`);
  `--linked-de-bruijn-graph`; force-calling given alleles (`--alleles`); a non-default `--pcr-indel-model` or
  `--min-base-quality-score`; extra read filters (`-RF`) or single annotations (`-A`; annotation groups `-G` are
  supported); `-bamout` (Mutect2 supports it); the same read in two inputs.
* **BaseRecalibrator:** covariates other than the default ones (`-cov`).
* **GenomicsDBImport:** adding samples to a workspace (`--genomicsdb-update-workspace-path`),
  `--merge-contigs-into-num-partitions`, `--header`, `--avoid-nio`, sample maps with index files; GVCFs with several
  samples, with FILTER values (GATK's have none) or with flag INFO fields; `-V` headers that disagree on a field's
  Number or Type (GATK promotes the line).
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
