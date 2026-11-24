# Documentation

The documentation is organised by purpose: a tutorial, task-oriented guides, reference material and explanations
of the design.

## Tutorial

* [Getting started](getting-started.md): building, alignment, and germline and somatic calling on the synthetic
  test data, with a comparison of the results against bwa-mem2 and GATK (about 10 minutes).

## How-to guides

* [How-to](how-to.md): aligning all lanes of a sample, tuning align-rs for a GPU, using variant-rs in place of
  `gatk` and within nf-core/sarek, running the pipeline, calling a genome in shards, joint calling, running
  without a GPU, and checking a result.
* [Building](building.md): requirements, GPU architectures, CPU-only builds, tests and Docker.
* [Containers](../containers/README.md): the two images, and the
  [use of the variant-rs image in sarek](../containers/variant-rs/README.md).
* [Pipeline](../pipeline/README.md): the Nextflow FASTQ-to-VCF pipeline, its parameters and its benchmarks.
* [Full test harness](../tests/harness/README.md): the public data and the procedure for running every comparison.

## Reference

* [Reference](reference.md): the commands of align-rs and variant-rs, the supported GATK tools and arguments,
  environment variables, exit status and [unsupported features](reference.md#unsupported-features).
* [GATK arguments by tool](gatk-arguments.md): the GATK arguments of each tool that are implemented, ignored, or
  accepted only at their default value (generated from `variant-rs <Tool> --help`).
* [Changelog](../CHANGELOG.md).

## Explanation

* [Architecture](architecture.md): the structure of each tool, the location of each component in the source code,
  and the requirements of a port to [other GPU vendors](architecture.md#other-gpu-vendors).
* [Optimizations](optimizations.md): the optimizations and the argument that each preserves the output, with code.
* [Validation](validation.md): the definition of identity, the levels of testing, the adversarial review and the
  known differences.

## Background

These documents were written during development. They are more detailed than the pages above and are partly
historical.

* Development notes for [align-rs](notes/align-rs.md) and [variant-rs](notes/variant-rs.md): how each behaviour
  of the reference tools was matched, the measurements, and the approaches that were tried and abandoned.
* [Pipeline benchmarks: measurement log](notes/pipeline-benchmarks.md): how each pipeline timing was taken, and
  earlier configurations.
* The stages of HaplotypeCaller, one crate each, with their own validation:
  [variant-rs/crates](../variant-rs/crates/README.md).
* Specifications written for the reimplementation of GATK components: [GenomicsDB merge](specs/genomicsdb_merge.md),
  [GenomicsDB / TileDB storage](specs/genomicsdb_tiledb.md),
  [GenotypeGVCFs annotations](specs/genotypegvcfs_annotations.md).
