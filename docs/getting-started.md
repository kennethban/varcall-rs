# Getting started

This tutorial covers alignment, germline calling and somatic calling with varcall-rs, starting from a completed
build, and takes about ten minutes. It uses the small synthetic data set of the test suite, so no data need to be
downloaded. Two comparisons are included: the records of align-rs against those of bwa-mem2, and the output of
variant-rs against that of GATK.

Requirements:

* varcall-rs, built as described in [building.md](building.md) (a GPU is optional);
* the `varcall-rs` conda environment from `environment.yml`, which provides bwa-mem2, samtools, bgzip and tabix;
* for step 5 only, GATK 4.7.0.0 (the release archive from the
  [GATK releases page](https://github.com/broadinstitute/gatk/releases/tag/4.7.0.0)).

## 1. Generate the data

```bash
REPO=$PWD                                    # run from the varcall-rs checkout
export PATH=$REPO/target/release:$PATH       # align-rs and variant-rs
python3 tests/synth/simulate.py ~/varcall-tutorial
cd ~/varcall-tutorial
```

`simulate.py` writes a 1.9 Mb genome that contains repeats, short tandem repeats and gaps, a germline sample
sequenced on two lanes, and a tumor-normal pair with three somatic clones. The simulation is seeded and therefore
produces the same files on every machine.

The reference is then indexed for bwa-mem2 (align-rs reads the same index), and the two VCF files used by BQSR and
Mutect2 are compressed and indexed:

```bash
bwa-mem2 index ref.fa
for v in known_sites germline_resource; do bgzip -k $v.vcf && tabix -p vcf $v.vcf.gz; done
```

## 2. Align

```bash
RG='@RG\tID:GERM.L1\tSM:GERM\tPL:ILLUMINA\tLB:GERM\tPU:FCGERM1.1'
align-rs mem -t 8 -K 100000000 -Y --gpu -R "$RG" -o GERM.L1.bam ref.fa germline_L1_R1.fq.gz germline_L1_R2.fq.gz
```

The options are those of bwa-mem2:

* `-K` sets the batch size. The output of bwa-mem2 depends on this value, so it should equal the value used in the
  run that is being reproduced. The pipelines of GATK and nf-core use `-K 100000000`.
* `--gpu` uses the GPU when one is present and the CPU otherwise.
* `-o x.bam` makes align-rs write the coordinate-sorted BAM and its index directly, with the result that
  `bwa-mem2 mem | samtools sort` followed by `samtools index` would give.

The remaining lanes of a sample can be aligned by a single process, which loads the index once. Each line of the
lanes file specifies an output BAM, a read group and the FASTQ files:

```bash
rg() { printf '@RG\\tID:%s\\tSM:%s\\tPL:ILLUMINA\\tLB:%s\\tPU:%s' $1 $2 $2 $3; }
printf '%s\t%s\t%s\t%s\n' \
  GERM.L2.bam "$(rg GERM.L2 GERM FCGERM1.2)" germline_L2_R1.fq.gz germline_L2_R2.fq.gz \
  TUMOR.bam   "$(rg TUMOR TUMOR FCTUMR1.1)"  tumor_R1.fq.gz       tumor_R2.fq.gz \
  NORMAL.bam  "$(rg NORMAL NORMAL FCNORM1.1)" normal_R1.fq.gz     normal_R2.fq.gz > lanes.tsv
align-rs mem -t 8 -K 100000000 -Y --gpu --lanes lanes.tsv ref.fa
```

**Comparison with bwa-mem2.** The first lane is aligned with bwa-mem2 and samtools, and the records of the two BAM
files are compared:

```bash
bwa-mem2 mem -t 8 -K 100000000 -Y -R "$RG" ref.fa germline_L1_R1.fq.gz germline_L1_R2.fq.gz 2> /dev/null \
  | samtools sort -o bwa.GERM.L1.bam - && samtools index bwa.GERM.L1.bam
cmp <(samtools view GERM.L1.bam) <(samtools view bwa.GERM.L1.bam) && echo identical
```

## 3. Call germline variants

variant-rs accepts the command lines of GATK. The following commands mark duplicates across both lanes,
recalibrate the base qualities and call variants:

```bash
variant-rs MarkDuplicates -I GERM.L1.bam -I GERM.L2.bam -O GERM.md.bam -M GERM.md.metrics
variant-rs BaseRecalibrator -R ref.fa -I GERM.md.bam --known-sites known_sites.vcf.gz -O GERM.recal.table
variant-rs ApplyBQSR -R ref.fa -I GERM.md.bam --bqsr-recal-file GERM.recal.table -O GERM.recal.bam
variant-rs HaplotypeCaller -R ref.fa -I GERM.recal.bam -O GERM.vcf.gz
```

Adding `-ERC GVCF` to HaplotypeCaller produces a GVCF, which `variant-rs GenotypeGVCFs` then genotypes as GATK
does. HaplotypeCaller can also apply the recalibration table while it reads the input
(`--bqsr-recal-file GERM.recal.table` with `GERM.md.bam`), which avoids writing the recalibrated BAM. The pipeline
uses this form.

## 4. Call somatic variants

The tumor and the normal sample are preprocessed in the same way, and the somatic Best Practices of GATK follow:

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

With GATK 4.7.0.0 installed, the same command can be run with `gatk` and the records compared. The header of a VCF
file differs only in the line that records the command, so the comparison is made on the data lines:

```bash
gatk HaplotypeCaller -R ref.fa -I GERM.recal.bam -O gatk.GERM.vcf.gz 2> /dev/null
cmp <(zgrep -v '^#' GERM.vcf.gz) <(zgrep -v '^#' gatk.GERM.vcf.gz) && echo identical
```

Any other step can be compared in the same way. `tests/selftest.sh` runs the whole of this tutorial and compares 34
outputs with those of GATK, once with the CPU PairHMM and once with the GPU PairHMM.

## Next steps

* The [Nextflow pipeline](../pipeline/README.md) processes a whole sample from FASTQ files, and
  [how-to.md](how-to.md) describes the use of variant-rs within nf-core/sarek.
* [reference.md](reference.md) lists every tool, option and environment variable.
* [validation.md](validation.md) describes how the outputs are checked on complete genomes.
