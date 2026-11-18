#!/usr/bin/env bash
# GATK 4.7 GetPileupSummaries + CalculateContamination truth tables (and FilterMutectCalls with them).
# usage: truth.sh [dream] [deep] [wgs]
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/contam" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
G="java -Xmx8g -jar $J"
D=../dream/data; P=$VARCALL_DATA
mkdir -p truth
gps() { # name bam vcf
  t0=$(date +%s); $G GetPileupSummaries -I $2 -V $3 -L $3 -O truth/$1.pileups.table > truth/$1.gps.log 2>&1
  echo "GPS $1 rc=$? $(( $(date +%s) - t0 ))s $(grep -vc '^#\|^contig' truth/$1.pileups.table) sites"; }
cc() { # out tumor [normal]
  local m=""; [ -n "$3" ] && m="-matched truth/$3.pileups.table"
  $G CalculateContamination -I truth/$2.pileups.table $m -O truth/$1.contamination.table --tumor-segmentation truth/$1.segments.table > truth/$1.cc.log 2>&1
  echo "CC $1 rc=$? $(tail -1 truth/$1.contamination.table) $(grep -vc '^#\|^contig' truth/$1.segments.table) segments"; }
for what in ${@:-dream deep wgs}; do case $what in
dream)
  for i in 1 2 3 4; do
    gps tumor_$i $D/tumor_$i.bam res/exac_20.b37.vcf.gz; gps normal_$i $D/normal_$i.bam res/exac_20.b37.vcf.gz
    cc dream_$i tumor_$i normal_$i; cc dream_${i}_to tumor_$i
    $G FilterMutectCalls -R $D/human_g1k_v37.20.21.fasta -V ../dream/truth/m2_$i.vcf --contamination-table truth/dream_$i.contamination.table \
      --tumor-segmentation truth/dream_$i.segments.table -O truth/fc_$i.vcf > truth/fc_$i.log 2>&1
    echo "FMC $i rc=$? PASS $(grep -v '^#' truth/fc_$i.vcf | awk '$7=="PASS"' | wc -l)"
  done ;;
deep)
  gps na12878_chr20 $P/markdup/input/na12878.chr20.bam res/exac_chr20.hg38.vcf.gz
  gps hg002_chr20 $P/hg002_pcrplus/bam/HG002.pcrplus.bam res/exac_chr20.hg38.vcf.gz
  cc deep na12878_chr20 hg002_chr20; cc deep_to na12878_chr20
  $G FilterMutectCalls -R ../dsmerge/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa -V ../dsmerge/gatk_ds.vcf --contamination-table truth/deep.contamination.table \
    --tumor-segmentation truth/deep.segments.table -O truth/fc_deep.vcf > truth/fc_deep.log 2>&1
  echo "FMC deep rc=$? PASS $(grep -v '^#' truth/fc_deep.vcf | awk '$7=="PASS"' | wc -l)" ;;
wgs)
  gps na12878_wgs $P/wgs/bam/NA12878.bam res/small_exac_common_3.hg38.vcf.gz
  gps hg002_wgs $P/hg002_pcrplus/bam/HG002.pcrplus.bam res/small_exac_common_3.hg38.vcf.gz
  cc wgs na12878_wgs hg002_wgs; cc wgs_to na12878_wgs; cc wgs_hg002_to hg002_wgs ;;
esac; done
