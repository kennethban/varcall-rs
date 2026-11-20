#!/usr/bin/env bash
# Inputs for the VQSR checks, in $VARCALL_DATA/vqsr:
#   resources/   HapMap 3.3, Omni 2.5, 1000G phase 1 SNPs, Axiom Exome Plus (Broad hg38 bundle; dbSNP and Mills come
#                from ../bqsr/resources)
#   in/NA12878.vcf.gz   single sample: GATK HaplotypeCaller on the NA12878 30x WGS (../wgs/truth/vcf, 23 shards joined)
#   in/cohort.vcf.gz    joint-called cohort: the 1000 Genomes 30x GATK joint calls (3202 samples), chr20:5-17 Mb,
#                first 50 samples kept (INFO not recomputed), FILTER cleared and INFO cut to GATK's annotations, so it
#                reads like GenotypeGVCFs output. Streamed with bcftools in six parallel pieces (~8 GB read, ~20 min).
# Needs curl, bcftools (with libcurl), tabix. Resumable: finished files are kept.
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
mkdir -p "$VARCALL_DATA/vqsr" && cd "$VARCALL_DATA/vqsr" || exit 1
mkdir -p resources in

B=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0
for f in hapmap_3.3.hg38.vcf.gz 1000G_omni2.5.hg38.vcf.gz 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
         Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz; do
  for x in "" .tbi; do
    [ -s resources/$f$x ] && continue
    curl -sSf -o resources/$f$x.part $B/$f$x && mv resources/$f$x.part resources/$f$x
  done
done

if [ ! -s in/NA12878.vcf.gz.tbi ]; then
  bcftools concat --no-version $(ls ../wgs/truth/vcf/shard*.vcf.gz) -Ou | bcftools sort -Oz -o in/NA12878.vcf.gz && tabix -p vcf in/NA12878.vcf.gz
fi

if [ ! -s in/cohort.vcf.gz.tbi ]; then
  K=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.recalibrated_variants.vcf.gz
  [ -s in/chr20.tbi ] || curl -sSf -o in/chr20.tbi $K.tbi
  bcftools query -l "$K##idx##in/chr20.tbi" | head -50 > in/samples50.txt
  keep=INFO/AC,INFO/AF,INFO/AN,INFO/BaseQRankSum,INFO/ClippingRankSum,INFO/DP,INFO/DS,INFO/END,INFO/ExcessHet,INFO/FS,INFO/InbreedingCoeff,INFO/MLEAC,INFO/MLEAF,INFO/MQ,INFO/MQ0,INFO/MQRankSum,INFO/QD,INFO/RAW_MQ,INFO/ReadPosRankSum,INFO/SOR
  piece() { # start end n: records starting in [start, end]
    [ -s in/piece$3.vcf.gz ] && return
    bcftools view -r chr20:$1-$2 --regions-overlap 0 -S in/samples50.txt -I -Ou "$K##idx##in/chr20.tbi" \
      | bcftools annotate -x FILTER -Ou | bcftools annotate -x "^$keep" -Oz -o in/piece$3.tmp.vcf.gz && mv in/piece$3.tmp.vcf.gz in/piece$3.vcf.gz
  }
  for i in 0 1 2 3 4 5; do
    s=$((5000000 + i * 2000000)); e=$((s + 1999999)); [ $i = 5 ] && e=17000000
    piece $s $e $i &
  done
  wait
  bcftools concat -n in/piece{0,1,2,3,4,5}.vcf.gz -o in/cohort.vcf.gz && tabix -p vcf in/cohort.vcf.gz && rm in/piece*.vcf.gz
fi
ls -la resources in
