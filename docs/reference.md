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

* **`mem`** takes bwa-mem2 `mem`'s options and gives its records. `-o x.bam` writes a coordinate-sorted BAM and its
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
| `ALIGN_RS_TRACE=FILE` | a line per stage run, extension device call (upload / sort / kernel / download), extension CPU step, mate-SW chunk (CPU preparation, upload / kernel / download) and final-pass step per window, on one time origin |
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
| ReblockGVCF, MergeVcfs, SortSam (coordinate), GatherBamFiles | |

Every GATK argument of a tool is in one of three groups:

* **mapped:** variant-rs implements it;
* **ignored:** it cannot change the output (logging, temporary directories, `--java-options`,
  `--native-pair-hmm-threads`, ...);
* **accepted at GATK's default value only.**

Everything else is refused with the reason, exit status 2: other values, unknown arguments, inputs variant-rs
cannot reproduce. The tables are in `variant-rs/gatkargs/`, generated from GATK 4.7's `--help`.

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
* **`--threads`, `--contig-workers`** size the parallelism.
* **`-L file --scatter-intervals`** calls every interval as its own GATK run, as sarek scatters HaplotypeCaller.

See [notes/variant-rs.md](notes/variant-rs.md) for the full option table and for how shards relate to GATK's
random-number streams.

### variant-rs environment variables

| environment | |
|---|---|
| `VARIANT_RS_THREADS` | Mutect2 / contamination: worker threads (default all CPUs) |
| `VARIANT_RS_CONTIG_THREADS`, `VARIANT_RS_DECODE_THREADS`, `VARIANT_RS_TWO_STAGE` | Mutect2 / HaplotypeCaller front-end sizing |
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
