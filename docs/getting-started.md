# Getting started

This tutorial takes you from a fresh build to germline and somatic calls in about ten minutes. It uses the small
synthetic data set that the tests use, so there is nothing to download. Along the way you check the two promises
varcall-rs makes: align-rs writes bwa-mem2's records, and variant-rs writes GATK's output.

You need:

* varcall-rs built as in [building.md](building.md) (a GPU is optional);
* the `varcall-rs` conda environment from `environment.yml` (it has bwa-mem2, samtools, bgzip and tabix);
* for step 5 only, GATK 4.7.0.0 (the release zip from the
  [GATK releases page](https://github.com/broadinstitute/gatk/releases/tag/4.7.0.0)).

## 1. Make the data

```bash
REPO=$PWD                                    # run from the varcall-rs checkout
export PATH=$REPO/target/release:$PATH       # align-rs and variant-rs
python3 tests/synth/simulate.py ~/varcall-tutorial
cd ~/varcall-tutorial
```

`simulate.py` writes a 1.9 Mb genome with repeats, STRs and gaps; a germline sample sequenced on two lanes; and a
tumor / normal pair with three somatic clones. It is seeded, so you get the same files as everyone else.

Index the reference for bwa-mem2 (align-rs reads the same index), and compress and index the two VCFs that BQSR and
Mutect2 use:

```bash
bwa-mem2 index ref.fa
for v in known_sites germline_resource; do bgzip -k $v.vcf && tabix -p vcf $v.vcf.gz; done
```

## 2. Align

```bash
RG='@RG\tID:GERM.L1\tSM:GERM\tPL:ILLUMINA\tLB:GERM\tPU:FCGERM1.1'
align-rs mem -t 8 -K 100000000 -Y --gpu -R "$RG" -o GERM.L1.bam ref.fa germline_L1_R1.fq.gz germline_L1_R2.fq.gz
```

These are bwa-mem2's options. `-K` matters: bwa-mem2's output depends on the batch size, so use the value your
reference run used (GATK's and nf-core's pipelines use `-K 100000000`). `--gpu` uses the GPU when there is one and
the CPU otherwise. With `-o x.bam`, align-rs writes the coordinate-sorted BAM and its index itself, as
`bwa-mem2 mem | samtools sort` and `samtools index` would.

A sample's other lanes can go through one process, which loads the index once. Each line of the lanes file gives an
output BAM, a read group and the FASTQ files:

```bash
rg() { printf '@RG\\tID:%s\\tSM:%s\\tPL:ILLUMINA\\tLB:%s\\tPU:%s' $1 $2 $2 $3; }
printf '%s\t%s\t%s\t%s\n' \
  GERM.L2.bam "$(rg GERM.L2 GERM FCGERM1.2)" germline_L2_R1.fq.gz germline_L2_R2.fq.gz \
  TUMOR.bam   "$(rg TUMOR TUMOR FCTUMR1.1)"  tumor_R1.fq.gz       tumor_R2.fq.gz \
  NORMAL.bam  "$(rg NORMAL NORMAL FCNORM1.1)" normal_R1.fq.gz     normal_R2.fq.gz > lanes.tsv
align-rs mem -t 8 -K 100000000 -Y --gpu --lanes lanes.tsv ref.fa
```

**Check it against bwa-mem2.** Align the first lane with bwa-mem2 and samtools, and compare the records:

```bash
bwa-mem2 mem -t 8 -K 100000000 -Y -R "$RG" ref.fa germline_L1_R1.fq.gz germline_L1_R2.fq.gz 2> /dev/null \
  | samtools sort -o bwa.GERM.L1.bam - && samtools index bwa.GERM.L1.bam
cmp <(samtools view GERM.L1.bam) <(samtools view bwa.GERM.L1.bam) && echo identical
```

## 3. Call germline variants

variant-rs takes GATK's command lines. Mark duplicates over both lanes, recalibrate base qualities, and call:

```bash
variant-rs MarkDuplicates -I GERM.L1.bam -I GERM.L2.bam -O GERM.md.bam -M GERM.md.metrics
variant-rs BaseRecalibrator -R ref.fa -I GERM.md.bam --known-sites known_sites.vcf.gz -O GERM.recal.table
variant-rs ApplyBQSR -R ref.fa -I GERM.md.bam --bqsr-recal-file GERM.recal.table -O GERM.recal.bam
variant-rs HaplotypeCaller -R ref.fa -I GERM.recal.bam -O GERM.vcf.gz
```

Add `-ERC GVCF` to HaplotypeCaller for a GVCF, then `variant-rs GenotypeGVCFs` as with GATK. HaplotypeCaller can
also apply the recalibration table while it reads (`--bqsr-recal-file GERM.recal.table` on `GERM.md.bam`), which
saves writing the recalibrated BAM; the pipeline does that.

## 4. Call somatic variants

The same preprocessing for the tumor and the normal, then GATK's somatic Best Practices:

```bash
for s in TUMOR NORMAL; do
  variant-rs MarkDuplicates -I $s.bam -O $s.md.bam -M $s.md.metrics
  variant-rs BaseRecalibrator -R ref.fa -I $s.md.bam --known-sites known_sites.vcf.gz -O $s.recal.table
  variant-rs ApplyBQSR -R ref.fa -I $s.md.bam --bqsr-recal-file $s.recal.table -O $s.recal.bam
  variant-rs GetPileupSummaries -I $s.recal.bam -V germline_resource.vcf.gz -L germline_resource.vcf.gz \
    -O $s.pileups.table
done
variant-rs Mutect2 -R ref.fa -I TUMOR.recal.bam -I NORMAL.recal.bam -normal NORMAL \
  --germline-resource germline_resource.vcf.gz --f1r2-tar-gz f1r2.tar.gz -O somatic.vcf.gz
variant-rs LearnReadOrientationModel -I f1r2.tar.gz -O priors.tar.gz
variant-rs CalculateContamination -I TUMOR.pileups.table -matched NORMAL.pileups.table \
  -O contamination.table --tumor-segmentation segments.table
variant-rs FilterMutectCalls -R ref.fa -V somatic.vcf.gz --contamination-table contamination.table \
  --tumor-segmentation segments.table --ob-priors priors.tar.gz -O somatic.filtered.vcf.gz
```

## 5. Compare with GATK

With GATK 4.7.0.0 installed, run the same command with `gatk` and compare the records. A VCF's header differs only
in the command-line line, so compare the data lines:

```bash
gatk HaplotypeCaller -R ref.fa -I GERM.recal.bam -O gatk.GERM.vcf.gz 2> /dev/null
cmp <(zgrep -v '^#' GERM.vcf.gz) <(zgrep -v '^#' gatk.GERM.vcf.gz) && echo identical
```

Any step works the same way. `tests/selftest.sh` runs this whole tutorial, compares 29 outputs with GATK's, and does
it once with the CPU PairHMM and once with the GPU one.

## What next

* Run a whole sample from FASTQ with the [Nextflow pipeline](../pipeline/README.md), or put variant-rs into
  nf-core/sarek: [how-to.md](how-to.md).
* Every tool, option and environment variable: [reference.md](reference.md).
* How the outputs are checked at scale: [validation.md](validation.md).
