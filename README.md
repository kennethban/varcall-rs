# varcall-rs

[![CI](https://github.com/kennethban/varcall-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/kennethban/varcall-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Short-read alignment and germline / somatic variant calling in Rust and CUDA, with **outputs byte-identical to
bwa-mem2 and GATK 4.7**: the same BAM records, VCF lines, metrics and tables. A 30x whole genome goes from FASTQ to a
filtered VCF in 34 minutes on one workstation, instead of about 5 hours.

* **align-rs** reimplements `bwa-mem2 mem` 2.2.1, plus `samtools sort` and `samtools index`. It runs seeding, SA
  lookups, seed extension and mate rescue on the GPU.
* **variant-rs** reimplements the GATK 4.7 / Picard tools of the GATK Best Practices: HaplotypeCaller, Mutect2 and
  its filters, MarkDuplicates, BaseRecalibrator, ApplyBQSR, GenotypeGVCFs and [more](docs/reference.md#variant-rs).
  Its PairHMM runs on the GPU, bit-exact with GATK's AVX-512 kernel.

## Highlights

* **Identical, not approximately equal.** Every output is checked byte for byte against the reference tools: on
  synthetic data in CI, and on whole genomes, exomes and tumor / normal pairs ([validation](docs/validation.md)).
* **GATK's command lines.** `variant-rs HaplotypeCaller -R ref.fa -I x.bam -O x.vcf.gz` takes the same arguments as
  `gatk`. An argument or input it does not reproduce exactly is refused with the reason, or passed to a real GATK.
* **Drop-in for nf-core/sarek.** The variant-rs image replaces sarek's GATK image for the steps it reproduces.
* **GPU optional.** Without a GPU everything runs on the CPU, with the same output.
* **A FASTQ → VCF pipeline** in Nextflow, germline (per sample, or joint calling with VQSR) and tumor / normal, that
  takes sarek's samplesheet.

## Speed

One workstation: 16-core Ryzen 9 9950X3D (32 threads), RTX 5070 Ti (16 GB), Linux (WSL2 with 48 GB; held to 30 GB the
germline run takes 41.6 min). The output is identical in every row.

| Whole pipeline | Data | Reference tools | varcall-rs | |
|---|---|---|---:|---:|
| germline: FASTQ → filtered VCF | NA12878 30x whole genome, 12 lanes | bwa-mem2 + samtools + GATK: ~5.2 h | **34.3 min** | **9.1×** |
| somatic: BAMs → filtered Mutect2 calls (+ HaplotypeCaller on the normal) | SEQC2 exome tumor / normal | nf-core/sarek with GATK: 46.5 min | **2.8 min** | **17×** |

| Single step | Data | Reference tool | varcall-rs | |
|---|---|---|---:|---:|
| align, sort, index | NA12878 30x, one lane (31.5 M pairs) | bwa-mem2 + samtools: 587 s | **124 s** | **4.7×** |
| HaplotypeCaller, VCF | NA12878 30x whole genome, 23 shards | GATK, 5 JVMs × 2 threads: 2 h 11 min | **12.9 min** | **10×** |
| HaplotypeCaller, GVCF | same | GATK: 2 h 23 min | **17.1 min** | **8.4×** |
| MarkDuplicates | HG002 30x whole genome | Picard: 42 min | **6 min** | **7×** |
| BaseRecalibrator | same | GATK: 99 min | **95 s** | **63×** |
| ApplyBQSR | same | GATK: 1 h 47 min | **6 min** | **18×** |

The germline baseline adds up measured steps: bwa-mem2 + samtools on one lane times 12, and sarek's GATK steps on the
same genome. On a server (2 × EPYC 9654, 384 threads, 2.2 TiB of RAM, NVMe; L40S GPUs with 48 GB) with 6 L40S and 384
threads, the germline pipeline takes **6.8 min**; with 4 L40S and 128 threads, **9.3 min** (NVIDIA Parabricks
`germline`: 16.8 min on the same 4 GPUs and 128 threads, and its own outputs, not GATK's). How each pipeline number
was measured: [pipeline/README.md](pipeline/README.md#benchmarks).

## Getting started

```bash
conda env create -f environment.yml && conda activate varcall-rs   # Rust, htslib, CUDA compiler, bwa-mem2, samtools
make -C align-rs/cuda && make -C variant-rs/cuda-pairhmm           # GPU code (skip for a CPU-only build)
RUSTFLAGS="-C target-cpu=native" cargo build --release              # target/release/align-rs, target/release/variant-rs
tests/selftest.sh                                                   # ~1 min: synthetic data vs bwa-mem2 and GATK

# align a lane: bwa-mem2's options, a sorted and indexed BAM out
align-rs mem -t 32 -K 100000000 -Y --gpu -R '@RG\tID:L1\tSM:S1' -o S1.L1.bam ref.fa R1.fq.gz R2.fq.gz
# GATK's command lines
variant-rs MarkDuplicates -I S1.L1.bam -O S1.md.bam -M S1.md.metrics
variant-rs BaseRecalibrator -R ref.fa -I S1.md.bam --known-sites dbsnp.vcf.gz -O S1.recal.table
variant-rs HaplotypeCaller -R ref.fa -I S1.md.bam --bqsr-recal-file S1.recal.table -O S1.vcf.gz
```

The [tutorial](docs/getting-started.md) walks through germline and somatic calling on the test data and compares
the results with bwa-mem2 and GATK.

## Installation

* **From source:** Linux x86-64, Rust 1.99, htslib, and the CUDA toolkit for the GPU code
  ([building.md](docs/building.md)). `environment.yml` installs all of them.
* **Hardware:** for a 30x genome, 32 threads, 48 GB of RAM and an NVIDIA GPU with 16 GB give the times above; a GPU
  is optional, and 30 GB of RAM works, more slowly ([memory and GPU per step](docs/building.md#hardware)).
* **Docker:** two images built from this repository: `align-rs`, and `variant-rs`, a drop-in `gatk` for sarek
  ([containers/](containers/README.md)).

## Documentation

| | |
|---|---|
| **Tutorial** | [Getting started](docs/getting-started.md): from a build to germline and somatic calls, checked against bwa-mem2 and GATK |
| **How-to guides** | [How-to](docs/how-to.md): align a sample's lanes, use variant-rs as `gatk`, run it in sarek, run the pipeline, shard a genome, joint-call · [Building](docs/building.md) · [Containers](containers/README.md) · [Pipeline](pipeline/README.md) |
| **Reference** | [Reference](docs/reference.md): commands, options, supported GATK arguments, environment variables, exit status · [Changelog](CHANGELOG.md) |
| **Explanation** | [Architecture](docs/architecture.md) · [Optimizations](docs/optimizations.md): fast without changing an output bit · [Validation](docs/validation.md): how "identical" is checked, known differences · [Development notes](docs/README.md#background) |

## Limitations

varcall-rs is the wrong tool if:

* **you need other GATK or bwa-mem2 versions.** The outputs are those of GATK 4.7.0.0 (Picard 3.5.0, htsjdk 5.0.0)
  and bwa-mem2 2.2.1. Other versions can give other outputs.
* **you need GATK tools or arguments it does not cover.** It covers each tool's defaults plus the arguments that WARP
  and nf-core/sarek pass ([tool by tool](docs/gatk-arguments.md)). HaplotypeCaller is diploid only; there is no
  CNN filtering, CNV / SV calling or RNA-seq ([what it does not cover](docs/reference.md#what-it-does-not-cover)).
  Everything else is refused, or passed to GATK with `GATK_FALLBACK`.
* **you run something other than Linux on x86-64**, or want one align-rs process to use several GPUs (the
  [pipeline](docs/how-to.md#use-several-gpus) spreads a sample's lanes over them; variant-rs's PairHMM uses every
  visible GPU).

It has been tested on two GPU models: RTX 5070 Ti (sm_120) and L40S (sm_89). The GPU code is built for sm_70 (sm_75
for align-rs) through sm_120, but has not run on other GPUs yet; reports are welcome. It is CUDA only; what a port to
AMD, Intel or Apple GPUs would need is in [architecture.md](docs/architecture.md#other-gpu-vendors).

## Contributing

Bug reports and pull requests are welcome. An output that differs from bwa-mem2's or GATK's is a bug: please
[report it](CONTRIBUTING.md#reporting-a-difference) with an input that shows it. [CONTRIBUTING.md](CONTRIBUTING.md)
describes the tests a change needs.

## Authors

varcall-rs is written by Kenneth Ban, with Claude (Anthropic's AI model, through Claude Code) as a contributor.
Much of the code, tests and documentation was written with Claude's assistance; commits it helped with carry a
`Co-Authored-By: Claude` trailer. All outputs are checked against bwa-mem2 and GATK as described in
[validation.md](docs/validation.md).

## Citing

If you use varcall-rs, please cite this repository ([CITATION.cff](CITATION.cff); GitHub's "Cite this repository"
formats it), and the tools whose results it reproduces: bwa-mem2 (Vasimuddin et al., IPDPS 2019) and GATK
(Van der Auwera & O'Connor, *Genomics in the Cloud*, O'Reilly 2020).

## Acknowledgements

varcall-rs reimplements algorithms of [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2),
[GATK](https://github.com/broadinstitute/gatk), [Picard](https://github.com/broadinstitute/picard),
[htsjdk](https://github.com/samtools/htsjdk), [GKL](https://github.com/Intel-HLS/GKL),
[GenomicsDB](https://github.com/GenomicsDB/GenomicsDB) and Apache Commons Math, and links
[htslib](https://github.com/samtools/htslib); their licenses are listed in [NOTICE](NOTICE). The test data come
from Genome in a Bottle, the 1000 Genomes Project, the SEQC2 consortium, the ICGC-TCGA DREAM challenge and the GATK
resource bundle ([tests/harness/README.md](tests/harness/README.md)).

## License

MIT, see [LICENSE](LICENSE).

## Disclaimer

varcall-rs is not affiliated with or endorsed by the Broad Institute, Intel or the bwa-mem2 authors. GATK and
Picard are named only to say which outputs it reproduces. It is research software, not a medical device, and is not
intended for clinical use.
