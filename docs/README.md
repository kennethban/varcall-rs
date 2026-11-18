# Documentation

The pages are sorted by what you want to do: learn, get a task done, look something up, or understand how it works.

## Tutorial

* [Getting started](getting-started.md): build, align and call germline and somatic variants on the synthetic test
  data, and check the results against bwa-mem2 and GATK (~10 minutes).

## How-to guides

* [How-to](how-to.md): align all lanes of a sample, tune align-rs for a GPU, use variant-rs in place of `gatk`, use it
  in nf-core/sarek, run the pipeline, call a genome in shards, joint-call, run without a GPU, check a result.
* [Building](building.md): requirements, GPU architectures, CPU-only builds, tests, Docker.
* [Containers](../containers/README.md): the two images, and [the sarek drop-in](../containers/variant-rs/README.md).
* [Pipeline](../pipeline/README.md): the Nextflow FASTQ → VCF pipeline, its parameters and its benchmarks.
* [Full test harness](../tests/harness/README.md): the public data and how to run every comparison.

## Reference

* [Reference](reference.md): align-rs and variant-rs commands, supported GATK tools and arguments, environment
  variables, exit status.
* [Changelog](../CHANGELOG.md).

## Explanation

* [Architecture](architecture.md): how each tool is put together, and where things are in the code.
* [Optimizations](optimizations.md): what makes it fast without changing an output bit, with code.
* [Validation](validation.md): what "identical" means, the test tiers, the adversarial review, known differences.

## Background

Written during development; more detailed than the pages above, and partly historical.

* Development notes: [align-rs](notes/align-rs.md), [variant-rs](notes/variant-rs.md): how each behaviour was
  matched, measurements, what was tried and dropped.
* The HaplotypeCaller stages, one crate each, with their own validation: [variant-rs/crates](../variant-rs/crates/README.md).
* Specifications written to reimplement GATK components: [GenomicsDB merge](specs/genomicsdb_merge.md),
  [GenomicsDB / TileDB storage](specs/genomicsdb_tiledb.md), [GenotypeGVCFs annotations](specs/genotypegvcfs_annotations.md).
