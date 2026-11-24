# variant-rs container: a drop-in `gatk` for nf-core/sarek

The image is sarek's own GATK 4.7.0.0 + samtools image with variant-rs added. The base is
`community.wave.seqera.io/library/gatk4-lite_samtools:9b6806e26dd32afe`, which sarek's MarkDuplicates and ApplyBQSR
modules use because they convert to CRAM with samtools in the same task. Its other GATK modules use the same GATK
without samtools, so one image serves them all.

A `gatk` wrapper takes the place of the image's `gatk`, and the real launcher becomes `gatk.real` beside it. Putting
the wrapper earlier on `PATH` would not work: the image's entrypoint (`/shell-hook.sh`, which Nextflow keeps) resets
`PATH`, and every command would then quietly run in GATK. GATK and its Java stay in the image only as the fallback;
the commands that variant-rs runs use no Java.

* `gatk <Tool> ...`: variant-rs runs every command it reproduces byte for byte. Everything else goes to the image's
  real `gatk` launcher (`GATK_FALLBACK`), with the original arguments (`--java-options` included): unsupported
  tools and arguments, and inputs it refuses.
* `gatk --version` prints `The Genome Analysis Toolkit (GATK) v4.7.0.0`, which nf-core's
  `sed -n '/GATK.*v/s/.*v//p'` reads as `4.7.0.0`; a second line names variant-rs.
* `gatk --list` / `--help` are answered by GATK.
* `VARIANT_RS_DRY_RUN=1 gatk <Tool> ...` prints what variant-rs would run, or why GATK would run it, and exits.

## Build

From the repository root (the context is limited to the sources by `.dockerignore`):

```bash
docker build -f containers/variant-rs/Dockerfile -t variant-rs:0.1.0-gatk4.7.0.0 .
```

* **CPU code:** built for `x86-64-v3` (AVX2). The CPU PairHMM picks AVX-512 at run time. Measured on HG002
  chr8:1-15M, HaplotypeCaller:
  * GPU: 4.0 s, the same as a `-C target-cpu=native` build;
  * CPU with AVX-512: 7.3 s (6.2 s native);
  * CPU without AVX-512: 16.6 s;
  * output identical in all of these cases.
* **GPU code:** a fat binary for sm_70 … sm_120 plus PTX, built with `-fmad=false -ftz=true` (IEEE operations, no
  fused multiply-add). It has been validated on sm_120 only (RTX 5070 Ti). Without a GPU, or without the
  container's GPU access, the CPU PairHMM runs and the output is the same.

## Use with sarek

```bash
nextflow run nf-core/sarek -r 3.12.0 -profile docker -c containers/variant-rs/sarek_variant_rs.config ...
```

`sarek_variant_rs.config` gives the image to these GATK processes:

* MarkDuplicates, BaseRecalibrator, GatherBQSRReports, ApplyBQSR;
* HaplotypeCaller, GenotypeGVCFs, and the MergeVcfs steps after them;
* Mutect2 and its gathers (MergeMutectStats, GatherPileupSummaries), GetPileupSummaries, CalculateContamination,
  LearnReadOrientationModel, FilterMutectCalls.

**CPUs.** By default these tasks keep sarek's CPU allocation: 2 per task, MarkDuplicates 6, GenotypeGVCFs 1. Sarek
sizes it for GATK, which uses little more than one core per task. variant-rs is multi-threaded and uses whatever the
task is given, so on a machine with cores to spare, `--variant_rs_cpus N` gives every variant-rs task N CPUs. It is still
capped by `process.resourceLimits` (`--max_cpus`), and a later `-c` config overrides it. Checked on sarek's test
scenarios: unset, every task gets sarek's CPUs; set to 3, every variant-rs task gets 3; outputs identical either way.

HaplotypeCaller and Mutect2 also get `--gpus all` (Docker) or `--nv` (Singularity / Apptainer); set
`params.variant_rs_gpu = false` without NVIDIA's container runtime. For Singularity, push the image to a registry or
convert it (`singularity build variant-rs.sif docker-daemon://variant-rs:0.1.0-gatk4.7.0.0`) and point
`params.variant_rs_container` at it.

## Checks

`containers/test_images.sh` (see `containers/README.md`) runs every check through the image's entrypoint and reports
which engine ran it:

* `tests/selftest.sh`'s variant-rs steps with the image's binary, CPU and GPU PairHMM: identical to GATK 4.7 / Picard;
* all libraries resolve; `gatk --version` gives `4.7.0.0` to nf-core's parser; an unsupported tool
  (IntervalListToBed) is handed to the image's GATK;
* with `VARCALL_DATA`: HaplotypeCaller on HG002 chr8:1-15M, GPU and CPU (identical to GATK 4.7.0.0, the release
  sarek pins, run in this image: 48,066 records, 4.1 s GPU / 6.9 s CPU vs GATK 2.6 min); MarkDuplicates over 12 lane
  BAMs, MergeVcfs and GetPileupSummaries in sarek's forms, and Mutect2 on DREAM 1: identical to their GATK truths.

## Open points

* **Threads:** variant-rs uses all CPUs the container can see (Rust's `available_parallelism` follows cgroup CPU quotas).
  With Docker the config caps each container at the task's CPUs (`--cpus`); Singularity / Apptainer has no such cap
  by default.
* **GPUs other than sm_120:** not run (only an RTX 5070 Ti here).
