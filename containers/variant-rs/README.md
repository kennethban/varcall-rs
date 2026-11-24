# variant-rs container: a replacement for `gatk` in nf-core/sarek

The image consists of the GATK 4.7.0.0 and samtools image used by sarek, with variant-rs added. The base image is
`community.wave.seqera.io/library/gatk4-lite_samtools:9b6806e26dd32afe`. The MarkDuplicates and ApplyBQSR modules of
sarek use this image because they convert to CRAM with samtools in the same task, and the other GATK modules of
sarek use the same GATK without samtools, so a single image serves all of them.

A wrapper replaces the `gatk` launcher of the image, and the original launcher is kept beside it as `gatk.real`.
Placing the wrapper earlier on `PATH` would not be sufficient, because the entrypoint of the image
(`/shell-hook.sh`, which Nextflow retains) resets `PATH`, and every command would then be run by GATK without
notice. GATK and its Java runtime remain in the image only as the fallback, and the commands run by variant-rs do
not use Java.

* `gatk <Tool> ...`: variant-rs runs every command that it reproduces byte for byte. Unsupported tools, unsupported
  arguments and refused inputs are passed to the original `gatk` launcher (`GATK_FALLBACK`) with the original
  arguments, including `--java-options`.
* `gatk --version` prints `The Genome Analysis Toolkit (GATK) v4.7.0.0`, which nf-core's
  `sed -n '/GATK.*v/s/.*v//p'` reads as `4.7.0.0`; a second line names variant-rs.
* `gatk --list` / `--help` are answered by GATK.
* `VARIANT_RS_DRY_RUN=1 gatk <Tool> ...` prints what variant-rs would run, or why GATK would run it, and exits.

## Build

The image is built from the repository root, where `.dockerignore` limits the build context to the sources:

```bash
docker build -f containers/variant-rs/Dockerfile -t variant-rs:0.1.0-gatk4.7.0.0 .
```

* **CPU code:** built for `x86-64-v3` (AVX2). The CPU PairHMM picks AVX-512 at run time. Measured on HG002
  chr8:1-15M, HaplotypeCaller:
  * GPU: 4.0 s, the same as a `-C target-cpu=native` build;
  * CPU with AVX-512: 7.3 s (6.2 s native);
  * CPU without AVX-512: 16.6 s;
  * output identical in all of these cases.
* **GPU code:** a binary for sm_70 to sm_120 together with PTX, built with `-fmad=false -ftz=true` (IEEE operations, no
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

**CPUs.** By default, these tasks retain the CPU allocation of sarek: 2 per task, 6 for MarkDuplicates and 1 for
GenotypeGVCFs. Sarek chooses these values for GATK, which uses little more than one core per task. variant-rs is
multi-threaded and uses all the CPUs that a task is given, so on a machine with unused cores, `--variant_rs_cpus N`
gives every variant-rs task N CPUs. This value remains limited by `process.resourceLimits` (`--max_cpus`), and a
later `-c` configuration overrides it. The option was checked on the test scenarios of sarek: when it is unset,
every task receives the CPUs assigned by sarek; when it is set to 3, every variant-rs task receives 3; and the
outputs are identical in both cases.

HaplotypeCaller and Mutect2 are also given `--gpus all` (Docker) or `--nv` (Singularity or Apptainer). Without the
NVIDIA container runtime, `params.variant_rs_gpu = false` must be set. For Singularity, push the image to a registry
or convert it (`singularity build variant-rs.sif docker-daemon://variant-rs:0.1.0-gatk4.7.0.0`) and point
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

## Limitations

* **Threads:** variant-rs uses all CPUs the container can see (Rust's `available_parallelism` follows cgroup CPU quotas).
  With Docker the config caps each container at the task's CPUs (`--cpus`); Singularity / Apptainer has no such cap
  by default.
* **GPUs other than sm_120:** the image has been run only on an RTX 5070 Ti.
