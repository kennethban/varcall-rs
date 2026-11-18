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

[Unreleased]: https://github.com/kennethban/varcall-rs/commits/main
