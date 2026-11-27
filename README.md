# varcall-rs

[![CI](https://github.com/kennethban/varcall-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/kennethban/varcall-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

varcall-rs is a GPU-accelerated reimplementation of bwa-mem2 2.2.1 and GATK 4.7.0.0 for short-read alignment and
germline and somatic variant calling, written in Rust and CUDA. It produces the same BAM records, VCF records,
metrics and recalibration tables as the reference tools when given the same arguments. On a workstation with one
GPU, it processes a 30x human genome from FASTQ files to a filtered VCF in 34 minutes, compared with about 5 hours
for the reference tools.

It consists of two programs:

* **align-rs** reimplements `bwa-mem2 mem` 2.2.1 together with `samtools sort` and `samtools index`. Seeding,
  suffix-array lookups, seed extension and mate rescue run on the GPU.
* **variant-rs** reimplements the GATK 4.7 and Picard tools of the GATK Best Practices, including HaplotypeCaller,
  Mutect2 and its filters, MarkDuplicates, BaseRecalibrator, ApplyBQSR and GenotypeGVCFs
  ([complete list](docs/reference.md#variant-rs)). Its PairHMM runs on the GPU and is bit-exact with the AVX-512
  kernel that GATK uses.

## Features

* **Identical output.** Outputs are compared byte for byte with those of the reference tools, on synthetic data in
  continuous integration and on public whole-genome, exome and tumor-normal data ([validation](docs/validation.md)).
* **GATK command-line interface.** variant-rs accepts the arguments of `gatk`, for example
  `variant-rs HaplotypeCaller -R ref.fa -I x.bam -O x.vcf.gz`. An argument or input that it cannot reproduce exactly
  is refused with an explanation, or is passed to an installed GATK.
* **Reproduction of scattered runs.** The output of GATK depends on how a genome is divided among separate runs,
  because each run has its own random-number stream. variant-rs reproduces a given division (a single run, the
  shards of WARP, or the intervals of nf-core/sarek) within one process
  ([how-to](docs/how-to.md#call-a-whole-genome-in-shards-in-one-process)).
* **Use within nf-core/sarek.** The variant-rs container image replaces the GATK image of sarek for the steps that
  variant-rs reproduces.
* **Optional GPU.** Without a GPU, all steps run on the CPU and give the same output.
* **A FASTQ-to-VCF pipeline.** A Nextflow pipeline covers germline calling (per sample, or joint calling with VQSR)
  and tumor-normal somatic calling, and accepts the samplesheet format of sarek.

## Performance

The first two tables were measured on a workstation with a 16-core Ryzen 9 9950X3D (32 threads), one RTX 5070 Ti
(16 GB), and Linux under WSL2 with 48 GB of RAM. In every row, the output of varcall-rs is identical to that of the
reference tools.

| Whole pipeline | Data | Reference tools | varcall-rs | Speed-up |
|---|---|---|---:|---:|
| Germline, FASTQ to filtered VCF | NA12878 30x whole genome, 12 lanes | bwa-mem2 + samtools + GATK: ~5.2 h | **34.0 min** | **9.2×** |
| Somatic, BAMs to filtered Mutect2 calls (with HaplotypeCaller on the normal) | SEQC2 exome tumor-normal pair | nf-core/sarek with GATK: 46.5 min | **2.8 min** | **17×** |

| Single step | Data | Reference tool | varcall-rs | Speed-up |
|---|---|---|---:|---:|
| Alignment, sorting and indexing | NA12878 30x, one lane (31.5 M pairs) | bwa-mem2 + samtools: 587 s | **124 s** | **4.7×** |
| HaplotypeCaller, VCF | NA12878 30x whole genome, 23 shards | GATK, 5 JVMs × 2 threads: 2 h 11 min | **6.2 min** | **21×** |
| HaplotypeCaller, GVCF | as above | GATK: 2 h 23 min | **9.3 min** | **15×** |
| MarkDuplicates | HG002 30x whole genome | Picard: 42 min | **3.1 min** | **13×** |
| BaseRecalibrator | as above | GATK: 99 min | **65 s** | **91×** |
| ApplyBQSR | as above | GATK: 1 h 47 min | **4.4 min** | **24×** |

The germline time is the mean of three runs (34.0 ± 0.2 min); with memory limited to 30 GB, the run takes 41.6 min.
The germline baseline was not run end to end. It is the sum of measured steps: bwa-mem2 and samtools on one lane,
multiplied by 12, and the GATK steps of sarek on the same genome.

The comparisons with other pipelines were run on a server with two EPYC 9654 processors (384 threads), 2.2 TiB of
RAM, NVMe storage and L40S GPUs (48 GB each). Each run was held to fixed resources on one NUMA node: **S** is 1 GPU,
32 threads and 64 GB; **L** is 4 GPUs, 128 threads and 256 GB; **C** is 32 threads and 64 GB without a GPU.
NVIDIA documents at least 100 GB of host memory for Parabricks, also on a single GPU, so S is below its stated
requirement and L meets it; L is the main Parabricks comparison. The Mutect2-alone comparison was run on S only.

| Germline, NA12878 30x, 12 FASTQ lanes to filtered VCF | S | L | C |
|---|---:|---:|---:|
| varcall-rs (three runs each on S and L) | **29.1 min** | **9.2 min** | 3.1 h |
| NVIDIA Parabricks 4.8.0-1 `germline` (three runs each) | 52.7 min | 18.3 min | |
| nf-core/sarek 3.12.0 with stock GATK | 7.7 h | | |
| bwa-mem2 + samtools + GATK, one process per step, run one after another | | | 9.6 h |

| Somatic, SEQC2 HCC1395 whole-genome tumor-normal BAMs to filtered Mutect2 calls | S |
|---|---:|
| varcall-rs (three runs) | **45.2 min** |
| nf-core/sarek 3.12.0 with the variant-rs image | 2.9 h |
| nf-core/sarek 3.12.0 with stock GATK | 9.4 h |
| Mutect2 alone on the same recalibrated BAMs: varcall-rs / Parabricks `mutectcaller` / GATK as one process | 8.9 min / 21.5 min / 12.8 h |

With the CPU threads and memory of L, the germline run takes 16.5 min on 1 GPU, 10.4 min on 2 and 9.2 min on 4.
Every output compared with GATK's was identical: all outputs of the germline genome; the Mutect2 VCF, its statistics
and the filtered VCF of the somatic genome (against GATK Mutect2 run as one process); and all published outputs of
sarek with the variant-rs image (against stock sarek). The outputs of Parabricks are concordant with those of GATK but
not identical to them (99.85% of germline sites and 99.97% of Mutect2 sites with the same alleles). The conditions of
each measurement, including which runs were repeated, are in [pipeline/README.md](pipeline/README.md#benchmarks).

## Quick start

```bash
conda env create -f environment.yml && conda activate varcall-rs   # Rust, htslib, CUDA compiler, bwa-mem2, samtools
make -C align-rs/cuda && make -C variant-rs/cuda-pairhmm           # GPU code (omit for a CPU-only build)
RUSTFLAGS="-C target-cpu=native" cargo build --release              # target/release/align-rs, target/release/variant-rs
tests/selftest.sh                                                   # ~1 min: synthetic data against bwa-mem2 and GATK

# Align one lane with the options of bwa-mem2; the output is a sorted and indexed BAM
align-rs mem -t 32 -K 100000000 -Y --gpu -R '@RG\tID:L1\tSM:S1' -o S1.L1.bam ref.fa R1.fq.gz R2.fq.gz
# Run the GATK steps with the arguments of GATK
variant-rs MarkDuplicates -I S1.L1.bam -O S1.md.bam -M S1.md.metrics
variant-rs BaseRecalibrator -R ref.fa -I S1.md.bam --known-sites dbsnp.vcf.gz -O S1.recal.table
variant-rs HaplotypeCaller -R ref.fa -I S1.md.bam --bqsr-recal-file S1.recal.table -O S1.vcf.gz
```

The [tutorial](docs/getting-started.md) covers germline and somatic calling on the test data and shows how to
compare the results with those of bwa-mem2 and GATK.

## Installation

* **From source.** varcall-rs builds on Linux x86-64 with Rust 1.99 and htslib, and requires the CUDA toolkit for
  the GPU code ([building.md](docs/building.md)). The file `environment.yml` installs all of these.
* **Hardware.** For a 30x genome, 32 threads, 48 GB of RAM and an NVIDIA GPU with 16 GB of memory give the times
  reported above. A GPU is optional, and 30 GB of RAM is sufficient at a lower speed
  ([memory and GPU requirements per step](docs/building.md#hardware)).
* **Docker.** Two images are built from this repository: `align-rs`, and `variant-rs`, which replaces `gatk` in
  sarek ([containers/](containers/README.md)).

## Documentation

| | |
|---|---|
| **Tutorial** | [Getting started](docs/getting-started.md): germline and somatic calling on the test data, with a comparison against bwa-mem2 and GATK |
| **How-to guides** | [How-to](docs/how-to.md): aligning the lanes of a sample, using variant-rs as `gatk`, running it in sarek, running the pipeline, calling a genome in shards, joint calling · [Building](docs/building.md) · [Containers](containers/README.md) · [Pipeline](pipeline/README.md) |
| **Reference** | [Reference](docs/reference.md): commands, options, supported GATK arguments, environment variables, exit status · [Changelog](CHANGELOG.md) |
| **Explanation** | [Architecture](docs/architecture.md) · [Optimizations](docs/optimizations.md): how each optimization preserves the output · [Validation](docs/validation.md): the definition of identity, the tests and the known differences · [Development notes](docs/README.md#background) |

## Limitations

* **Versions.** The outputs are those of GATK 4.7.0.0 (Picard 3.5.0, htsjdk 5.0.0) and bwa-mem2 2.2.1, whose source
  align-rs was ported from. bwa-mem2 2.3 gives the same alignments; its release still reports itself as 2.2.1, and it
  is the build that `environment.yml` installs and the server benchmarks used. Other versions of the reference tools may give
  different outputs.
* **Coverage of GATK.** variant-rs supports the default arguments of each tool and the arguments that WARP and
  nf-core/sarek pass ([per-tool list](docs/gatk-arguments.md)). HaplotypeCaller is limited to diploid calling, and
  CNN filtering, copy-number and structural variant calling and RNA-seq are not supported
  ([unsupported features](docs/reference.md#unsupported-features)). Any other tool or argument is refused, or is
  passed to GATK when `GATK_FALLBACK` is set.
* **Platform.** Only Linux on x86-64 is supported, and the GPU code requires CUDA. The requirements of a port to
  AMD, Intel or Apple GPUs are discussed in [architecture.md](docs/architecture.md#other-gpu-vendors).
* **Multiple GPUs.** One align-rs process uses one GPU. The [pipeline](docs/how-to.md#use-several-gpus) distributes
  the lanes of a sample over several GPUs, and the PairHMM of variant-rs uses every visible GPU.
* **Tested hardware.** varcall-rs has been tested on two GPU models, the RTX 5070 Ti (sm_120) and the L40S (sm_89).
  The GPU code is built for sm_70 (sm_75 for align-rs) to sm_120 but has not been run on other models.

## Contributing

Bug reports and pull requests are welcome. Any output that differs from that of bwa-mem2 or GATK is treated as a
bug and should be [reported](CONTRIBUTING.md#reporting-a-difference) with an input that reproduces it.
[CONTRIBUTING.md](CONTRIBUTING.md) describes the tests that a change requires.

## Authors

varcall-rs is written by Kenneth Ban, with Claude (the AI model of Anthropic, through Claude Code) as a
contributor. Much of the code, tests and documentation was written with the assistance of Claude, and the commits
to which it contributed carry a `Co-Authored-By: Claude` trailer. All outputs are checked against bwa-mem2 and GATK
as described in [validation.md](docs/validation.md).

## Citation

If you use varcall-rs, please cite this repository ([CITATION.cff](CITATION.cff); the "Cite this repository" button
on GitHub formats the entry) and the tools whose results it reproduces: bwa-mem2 (Vasimuddin et al., IPDPS 2019)
and GATK (Van der Auwera & O'Connor, *Genomics in the Cloud*, O'Reilly 2020).

## Acknowledgements

varcall-rs reimplements algorithms of [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2),
[GATK](https://github.com/broadinstitute/gatk), [Picard](https://github.com/broadinstitute/picard),
[htsjdk](https://github.com/samtools/htsjdk), [GKL](https://github.com/IntelLabs/GKL),
[GenomicsDB](https://github.com/GenomicsDB/GenomicsDB) and Apache Commons Math, and links
[htslib](https://github.com/samtools/htslib). Their licenses are listed in [NOTICE](NOTICE). The test data come
from Genome in a Bottle, the 1000 Genomes Project, the SEQC2 consortium, the ICGC-TCGA DREAM challenge and the GATK
resource bundle ([tests/harness/README.md](tests/harness/README.md)).

## License

varcall-rs is distributed under the MIT license ([LICENSE](LICENSE)).

## Disclaimer

varcall-rs is not affiliated with or endorsed by the Broad Institute, Intel or the authors of bwa-mem2. GATK and
Picard are named only to identify the outputs that varcall-rs reproduces. varcall-rs is research software. It is
not a medical device and is not intended for clinical use.
