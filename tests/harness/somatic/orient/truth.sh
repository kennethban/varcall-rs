#!/usr/bin/env bash
# GATK 4.7 Mutect2 --f1r2-tar-gz, LearnReadOrientationModel, FilterMutectCalls --ob-priors. usage: truth.sh [1 2 3 4 deep]
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/orient" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
G="java -Xmx16g -jar $J"
D=../dream/data; P=$VARCALL_DATA
for x in ${@:-1 2 3 4 deep}; do
  t0=$(date +%s)
  if [ $x = deep ]; then
    R=../dsmerge/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
    $G Mutect2 -R $R -I $P/markdup/input/na12878.chr20.bam -I $P/hg002_pcrplus/bam/HG002.pcrplus.bam -normal HG002 -L ../dsmerge/ds.intervals \
      --f1r2-tar-gz truth/f1r2_deep.tar.gz -O truth/m2_deep.vcf > truth/m2_deep.log 2>&1
  else
    R=$D/human_g1k_v37.20.21.fasta
    N=$(samtools view -H $D/normal_$x.bam | grep -m1 '^@RG' | tr '\t' '\n' | grep '^SM:' | cut -c4-)
    $G Mutect2 -R $R -I $D/tumor_$x.bam -I $D/normal_$x.bam -normal "$N" -L 20 --germline-resource $D/very-small-gnomad.vcf --f1r2-tar-gz truth/f1r2_$x.tar.gz -O truth/m2_$x.vcf > truth/m2_$x.log 2>&1
  fi
  echo "M2 $x rc=$? $(( $(date +%s) - t0 ))s"
  $G LearnReadOrientationModel -I truth/f1r2_$x.tar.gz -O truth/priors_$x.tar.gz > truth/lrom_$x.log 2>&1; echo "LROM $x rc=$?"
  $G FilterMutectCalls -R $R -V truth/m2_$x.vcf --ob-priors truth/priors_$x.tar.gz -O truth/fo_$x.vcf > truth/fo_$x.log 2>&1
  echo "FMC $x rc=$? PASS $(grep -v '^#' truth/fo_$x.vcf | awk '$7=="PASS"' | wc -l) orientation $(grep -v '^#' truth/fo_$x.vcf | grep -c orientation)"
done
