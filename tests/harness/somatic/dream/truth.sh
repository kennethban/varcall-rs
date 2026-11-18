#!/usr/bin/env bash
# GATK 4.7 Mutect2 tumor-normal reference runs on the DREAM pairs (default arguments, -L 20, small gnomAD).
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/dream" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
for i in ${@:-1 2 3 4}; do
  N=$(samtools view -H data/normal_$i.bam | grep -m1 '^@RG' | grep -o 'SM:[^[:space:]]*' | cut -c4-)
  t0=$(date +%s)
  java -Xmx8g -jar $J Mutect2 -R data/human_g1k_v37.20.21.fasta -I data/tumor_$i.bam -I data/normal_$i.bam -normal "$N" \
    -L 20 --germline-resource data/very-small-gnomad.vcf -O truth/m2_$i.vcf > truth/m2_$i.log 2>&1
  echo "pair $i rc=$? $(( $(date +%s) - t0 ))s records $(grep -vc '^#' truth/m2_$i.vcf)"
done
