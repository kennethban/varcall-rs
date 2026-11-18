#!/usr/bin/env bash
# Test FASTQs and samplesheets (data/):
#   germline: NA12878 chr20:10-15M ($VD/pipeline/r1.fq.gz, r2.fq.gz: 1.19M pairs) cut into two lanes
#   somatic:  HCC1143 tumor / normal (GATK somatic workshop, hg38 chr17 subset), FASTQs per read group (the
#             normal is paired-end, the workshop's tumor BAM single-end)
set -eo pipefail
cd "$(dirname "$0")"
mkdir -p data; D=$PWD/data
# germline: first and second half of the pairs (4 FASTQ lines per read)
n=$(( $(zcat $VD/pipeline/r1.fq.gz | wc -l) / 8 * 4 ))
for r in 1 2; do
  zcat $VD/pipeline/r$r.fq.gz | awk -v n=$n -v a=data/NA12878_L1_R$r.fq -v b=data/NA12878_L2_R$r.fq '{ print > (NR <= n ? a : b) }'
  gzip -1 -f data/NA12878_L1_R$r.fq data/NA12878_L2_R$r.fq
done
{ echo patient,sample,status,lane,fastq_1,fastq_2
  for l in L1 L2; do echo NA12878,NA12878,0,$l,$D/NA12878_${l}_R1.fq.gz,$D/NA12878_${l}_R2.fq.gz; done; } > germline.csv
# somatic: per read group
VD=${VARCALL_DATA:?set VARCALL_DATA to the test data directory}; H=$VD/somatic/hcc1143/data
echo patient,sample,status,lane,fastq_1,fastq_2 > somatic.csv
for s in normal tumor; do
  st=$([ $s = tumor ] && echo 1 || echo 0)
  for rg in $(samtools view -H $H/$s.bam | awk -F'\t' '/^@RG/{sub("ID:","",$2); print $2}'); do
    f=data/HCC1143_${s}_$rg
    if [ "$(samtools view -c -f 1 -r $rg $H/$s.bam)" = 0 ]; then
      # single-end reads (the workshop's tumor BAM): one FASTQ, no fastq_2
      samtools view -u -r $rg $H/$s.bam | samtools collate -u -O - data/coll.$s | \
        samtools fastq -F 0x900 -0 $f.R1.fq.gz -n - 2>/dev/null
      echo HCC1143,HCC1143_$s,$st,$rg,$D/HCC1143_${s}_$rg.R1.fq.gz, >> somatic.csv
    else
      samtools view -u -r $rg $H/$s.bam | samtools collate -u -O - data/coll.$s | \
        samtools fastq -F 0x900 -1 $f.R1.fq.gz -2 $f.R2.fq.gz -0 /dev/null -s /dev/null -n - 2>/dev/null
      echo HCC1143,HCC1143_$s,$st,$rg,$D/HCC1143_${s}_$rg.R1.fq.gz,$D/HCC1143_${s}_$rg.R2.fq.gz >> somatic.csv
    fi
  done
done
cat germline.csv somatic.csv
