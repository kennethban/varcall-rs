# How-to guides

Short recipes for common tasks. They assume you have built varcall-rs ([building.md](building.md)) or use the
[images](../containers/README.md). [getting-started.md](getting-started.md) shows the basic commands end to end, and
[reference.md](reference.md) lists every option.

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

One align-rs process can align every lane of a sample. The index is loaded once, and each lane's final BAM pass
overlaps the next lane's alignment. Each lane's BAM is the one aligning it alone would give.

```bash
# lanes.tsv: output BAM <TAB> read group <TAB> R1 [<TAB> R2], one line per lane
align-rs mem -t 32 -K 100000000 -Y --gpu --lanes lanes.tsv ref.fa
variant-rs MarkDuplicates -I S1.L1.bam -I S1.L2.bam ... -O S1.md.bam -M S1.md.metrics
```

The read group is written as in `-R`, with `\t` between fields. MarkDuplicates merges the lanes, as Picard does.

## Tune align-rs for your GPU

The GPU kernels' launch settings were tuned on one GPU. `tune` measures them on yours (about 35 s plus the index
load) and saves them; later `mem` runs pick them up. No setting changes the output.

```bash
align-rs tune -t 32 -K 100000000 ref.fa R1.fq.gz R2.fq.gz
```

## Use variant-rs in place of `gatk`

variant-rs takes GATK's command lines. Point `GATK_FALLBACK` at a real GATK launcher, and anything variant-rs does
not reproduce exactly is passed to GATK with the original arguments:

```bash
export GATK_FALLBACK=/opt/gatk-4.7.0.0/gatk
variant-rs HaplotypeCaller -R ref.fa -I S1.bam -O S1.vcf.gz          # variant-rs runs it
variant-rs CollectWgsMetrics -I S1.bam -O S1.wgs.txt -R ref.fa      # not a variant-rs tool: GATK runs it
VARIANT_RS_DRY_RUN=1 variant-rs HaplotypeCaller ...                 # print what would run, or why it is refused
```

`variant-rs <Tool> --help` lists the arguments variant-rs supports for that tool. The
[variant-rs image](../containers/variant-rs/README.md) installs this as the `gatk` command.

## Use variant-rs in nf-core/sarek

Build the variant-rs image, then add its config to your sarek run:

```bash
docker build -f containers/variant-rs/Dockerfile -t variant-rs:0.1.0-gatk4.7.0.0 .
nextflow run nf-core/sarek -r 3.12.0 -profile docker -c containers/variant-rs/sarek_variant_rs.config ...
```

The GATK steps that variant-rs reproduces run in the image; sarek's other steps are unchanged. Give the variant-rs
tasks more CPUs with `--variant_rs_cpus N`, and set `params.variant_rs_gpu = false` without NVIDIA's container
runtime. Singularity / Apptainer and the other settings:
[containers/variant-rs/README.md](../containers/variant-rs/README.md).

## Run the FASTQ → VCF pipeline

The Nextflow pipeline in [pipeline/](../pipeline/README.md) runs the same steps as sarek's germline and somatic
paths, with align-rs and variant-rs, and without sarek's splitting into ~350 interval tasks. It takes sarek's
samplesheet format.

```bash
pipeline/run.sh -params-file my.yaml    # nextflow run pipeline/main.nf on the idle GPUs, threads from the machine
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

All parameters: [pipeline/README.md](../pipeline/README.md).

## Use several GPUs

An align-rs process uses one GPU. A variant-rs process spreads its PairHMM over every visible GPU
(`VARIANT_RS_GPUS=N` for at most N), though HaplotypeCaller on one sample gains little beyond 2. The pipeline uses
several: `pipeline/run.sh` takes the GPUs no
process is using (`VARCALL_GPUS=N` for at most N), or list them yourself when you launch. Each alignment task gets
the machine's threads divided by the GPUs (`--align_cpus`), the steps after it the machine's (`--cpus`, at most 128).

```bash
VARCALL_GPUS=4 pipeline/run.sh -params-file my.yaml
CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7 pipeline/run.sh -params-file my.yaml
```

Each sample's lanes are split over one align-rs task per GPU; HaplotypeCaller runs each sample on up to 2 GPUs
(`--hc_gpus`) and Mutect2 each tumor on one. The outputs are the same as with one GPU. Memory and other details:
[pipeline/README.md](../pipeline/README.md#several-gpus). Outside the pipeline, run one align-rs process per GPU
with `CUDA_VISIBLE_DEVICES=N` (per lane or per sample), and give variant-rs the GPUs it should share.

## Call a whole genome in shards, in one process

Pipelines such as WARP call a genome as many separate GATK runs, one per shard. Their output depends on the
sharding, because each GATK run has its own random-number stream. variant-rs reproduces this in one process:

```bash
# one output, equal to concatenating 23 separate GATK runs over the shard files
variant-rs -I S1.bam -R ref.fa --shard s/shard00.intervals ... --shard s/shard22.intervals -O S1.vcf.gz
# sarek's form: one interval list, each interval its own GATK run
variant-rs HaplotypeCaller -R ref.fa -I S1.bam -L intervals.bed --scatter-intervals -O S1.vcf.gz
# new shard files, balanced from the BAM index alone
variant-rs plan -I S1.bam -R ref.fa -o plan96 --pieces 96
```

## Joint-call several samples

As with GATK: a GVCF per sample, then GenotypeGVCFs on a combined GVCF or a GenomicsDB workspace.

```bash
variant-rs HaplotypeCaller -R ref.fa -I S1.bam -O S1.g.vcf.gz -ERC GVCF
variant-rs GenotypeGVCFs -R ref.fa -V gendb://my_workspace -O cohort.vcf.gz
```

The GenomicsDB workspace is made by GATK's `GenomicsDBImport`, which variant-rs does not reimplement.

## Run without a GPU

Everything runs on the CPU and gives the same output. Build without the two `make` steps ([building.md](building.md)),
or run a GPU build on a machine without a GPU: align-rs `--gpu` and variant-rs `--pairhmm auto` (the default) fall
back to the CPU. `--pairhmm cuda` stops instead of falling back.

## Check a result against GATK

Run the same command with GATK 4.7.0.0 and compare. For a VCF, compare the data lines (the header differs only in the
command-line line); for a BAM, the records:

```bash
cmp <(zgrep -v '^#' variant-rs.vcf.gz) <(zgrep -v '^#' gatk.vcf.gz) && echo identical
cmp <(samtools view variant-rs.bam) <(samtools view gatk.bam) && echo identical
```

If they differ, please [open an issue](../CONTRIBUTING.md#reporting-a-difference): an output difference is a bug.
