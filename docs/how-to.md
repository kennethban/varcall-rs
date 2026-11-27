# How-to guides

This page gives short procedures for common tasks. Each assumes that varcall-rs has been built
([building.md](building.md)) or that the [container images](../containers/README.md) are used.
[getting-started.md](getting-started.md) presents the basic commands in sequence, and [reference.md](reference.md)
lists every option.

* [Align all lanes of a sample](#align-all-lanes-of-a-sample)
* [Tune align-rs for your GPU](#tune-align-rs-for-your-gpu)
* [Use variant-rs in place of `gatk`](#use-variant-rs-in-place-of-gatk)
* [Use variant-rs in nf-core/sarek](#use-variant-rs-in-nf-coresarek)
* [Run the FASTQ → VCF pipeline](#run-the-fastq--vcf-pipeline)
* [Use several GPUs](#use-several-gpus)
* [Call a whole genome in shards, in one process](#call-a-whole-genome-in-shards-in-one-process)
* [Joint-call several samples](#joint-call-several-samples)
* [Run without a GPU](#run-without-a-gpu)
* [Check a result against GATK](#check-a-result-against-gatk)

## Align all lanes of a sample

A single align-rs process can align every lane of a sample. The index is loaded once, and the final BAM pass of
each lane overlaps the alignment of the next. The BAM file of each lane is identical to the file obtained by
aligning that lane on its own.

```bash
# lanes.tsv: output BAM <TAB> read group <TAB> R1 [<TAB> R2], one line per lane
align-rs mem -t 32 -K 100000000 -Y --gpu --lanes lanes.tsv ref.fa
variant-rs MarkDuplicates -I S1.L1.bam -I S1.L2.bam ... -O S1.md.bam -M S1.md.metrics
```

The read group is written as for `-R`, with `\t` between its fields. MarkDuplicates merges the lanes, as Picard
does.

## Tune align-rs for your GPU

The launch settings of the GPU kernels were tuned on one GPU model. The `tune` command measures them on the GPU at
hand, which takes about 35 s in addition to loading the index, and saves the result for later `mem` runs. None of
these settings changes the output.

```bash
align-rs tune -t 32 -K 100000000 ref.fa R1.fq.gz R2.fq.gz
```

## Use variant-rs in place of `gatk`

variant-rs accepts the command lines of GATK. When `GATK_FALLBACK` is set to the path of a GATK launcher, any
command that variant-rs cannot reproduce exactly is passed to GATK with its original arguments:

```bash
export GATK_FALLBACK=/opt/gatk-4.7.0.0/gatk
variant-rs HaplotypeCaller -R ref.fa -I S1.bam -O S1.vcf.gz          # run by variant-rs
variant-rs CollectWgsMetrics -I S1.bam -O S1.wgs.txt -R ref.fa      # not a variant-rs tool: run by GATK
VARIANT_RS_DRY_RUN=1 variant-rs HaplotypeCaller ...                 # print what would run, or why it is refused
```

`variant-rs <Tool> --help` lists the arguments that variant-rs supports for a tool. The
[variant-rs image](../containers/variant-rs/README.md) installs this arrangement as the `gatk` command.

## Use variant-rs in nf-core/sarek

Build the variant-rs image and add its configuration file to the sarek run. To use the published image without a
build, add `--variant_rs_container ghcr.io/kennethban/variant-rs:0.1.0-gatk4.7.0.0` to the sarek command in place
of the first line:

```bash
docker build -f containers/variant-rs/Dockerfile -t variant-rs:0.1.0-gatk4.7.0.0 .
nextflow run nf-core/sarek -r 3.12.0 -profile docker -c containers/variant-rs/sarek_variant_rs.config ...
```

The GATK steps that variant-rs reproduces then run in this image, and the other steps of sarek are unchanged. The
option `--variant_rs_cpus N` gives the variant-rs tasks more CPUs, and `params.variant_rs_gpu = false` is required
when the NVIDIA container runtime is not installed. Singularity, Apptainer and the remaining settings are described
in [containers/variant-rs/README.md](../containers/variant-rs/README.md).

## Run the FASTQ → VCF pipeline

The Nextflow pipeline in [pipeline/](../pipeline/README.md) runs the same steps as the germline and somatic paths
of sarek, using align-rs and variant-rs, and does not divide the genome into the approximately 350 interval tasks
that sarek creates. It accepts the samplesheet format of sarek. With `joint_germline: true`, the germline samples
are called jointly ([see below](#joint-call-several-samples)).

```bash
pipeline/run.sh -params-file my.yaml    # runs pipeline/main.nf on the idle GPUs, with threads set from the machine
```

```yaml
# my.yaml
input: samples.csv           # patient,sample,status,lane,fastq_1,fastq_2 (status 1 = tumor)
outdir: results
tools: haplotypecaller,mutect2
fasta: /ref/Homo_sapiens_assembly38.fasta
known_sites: /ref/dbsnp.vcf.gz,/ref/mills.vcf.gz,/ref/known_indels.vcf.gz
germline_resource: /ref/af-only-gnomad.hg38.vcf.gz
intervals: /ref/wgs_calling_regions.hg38.interval_list
```

All parameters are listed in [pipeline/README.md](../pipeline/README.md).

## Use several GPUs

The two programs use GPUs differently:

* An **align-rs** process uses one GPU.
* A **variant-rs** process distributes its PairHMM over every visible GPU (`VARIANT_RS_GPUS=N` limits it to N).
  HaplotypeCaller on a single sample gains little from more than 2 GPUs.

The pipeline uses both. `pipeline/run.sh` selects the GPUs that no process is using (`VARCALL_GPUS=N` limits it to
N); alternatively, the GPUs can be listed explicitly:

```bash
VARCALL_GPUS=4 pipeline/run.sh -params-file my.yaml
CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7 pipeline/run.sh -params-file my.yaml
```

With several GPUs, the pipeline divides the work as follows. The outputs are the same as with one GPU.

| Step | GPUs | Threads |
|---|---|---|
| Alignment | the lanes of each sample are divided among one align-rs task per GPU | the threads of the machine divided by the number of GPUs (`--align_cpus`) |
| HaplotypeCaller | up to 2 per sample (`--hc_gpus`) | the threads of the machine, at most 128 (`--cpus`) |
| Mutect2 | 1 per tumor | as for HaplotypeCaller |

Memory requirements and further details are given in [pipeline/README.md](../pipeline/README.md#several-gpus).

Outside the pipeline, one align-rs process can be run per GPU by setting `CUDA_VISIBLE_DEVICES=N` for each lane or
sample, and variant-rs can be given the GPUs that it should share.

## Call a whole genome in shards, in one process

Pipelines such as WARP and sarek call a genome as many separate GATK runs, one per shard (a set of intervals), and
join the results. The output of GATK depends on this division, because each run has its own random-number stream.
The same genome called in one run and in 23 shards gives VCF files that agree on nearly every call but differ in
many records, mostly in the QD annotation alone.

variant-rs can reproduce either result in a single process. By default it reproduces a single GATK run over the
whole input. With the options below it treats each shard as a separate GATK run, processes the shards in parallel,
and writes one output:

```bash
# As in WARP: one output, equal to the concatenation of 23 separate GATK runs over the shard files
variant-rs -I S1.bam -R ref.fa --shard s/shard00.intervals ... --shard s/shard22.intervals -O S1.vcf.gz
# As in sarek: one interval list, each interval treated as a separate GATK run
variant-rs HaplotypeCaller -R ref.fa -I S1.bam -L intervals.bed --scatter-intervals -O S1.vcf.gz
```

To reproduce an existing GATK result, use the shard files or the interval list with which it was produced. For a
new analysis, the `plan` command writes load-balanced shard files from the BAM index alone:

```bash
variant-rs plan -I S1.bam -R ref.fa -o plan96 --pieces 96
```

## Joint-call several samples

Joint calling follows the procedure of GATK: a GVCF is produced for each sample, and GenotypeGVCFs is then run on a
combined GVCF or on a GenomicsDB workspace.

```bash
variant-rs HaplotypeCaller -R ref.fa -I S1.bam -O S1.g.vcf.gz -ERC GVCF
variant-rs GenomicsDBImport --genomicsdb-workspace-path my_workspace --sample-name-map samples.map -L chr20 \
  --batch-size 50 --consolidate
variant-rs GenotypeGVCFs -R ref.fa -V gendb://my_workspace -O cohort.vcf.gz
```

The workspace contains the same cells and metadata as a workspace written by GATK, and GenotypeGVCFs of GATK can
read it.

VariantRecalibrator and ApplyVQSR, or VariantFiltration, are then applied as in GATK. The pipeline performs all of
these steps from FASTQ files or BAM files when `joint_germline: true` is set, which corresponds to the option of
sarek ([pipeline/README.md](../pipeline/README.md)).

## Run without a GPU

All steps run on the CPU and give the same output. Either build without the two `make` steps
([building.md](building.md)), or run a GPU build on a machine without a GPU. In the second case, `--gpu` in
align-rs and `--pairhmm auto` in variant-rs (the default) fall back to the CPU, whereas `--pairhmm cuda` stops with
an error.

## Check a result against GATK

Run the same command with GATK 4.7.0.0 and compare the outputs. For a VCF file, compare the data lines, because the
header differs in the line that records the command. For a BAM file, compare the records:

```bash
cmp <(zgrep -v '^#' variant-rs.vcf.gz) <(zgrep -v '^#' gatk.vcf.gz) && echo identical
cmp <(samtools view variant-rs.bam) <(samtools view gatk.bam) && echo identical
```

Any difference is a bug and should be [reported](../CONTRIBUTING.md#reporting-a-difference).
