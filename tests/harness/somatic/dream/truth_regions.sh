#!/usr/bin/env bash
# GATK Mutect2 assembly-region / activity-profile outputs on the DREAM pairs (front-end truth).
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/dream" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
for i in ${@:-1 2 3 4}; do
  N=$(samtools view -H data/normal_$i.bam | grep -m1 '^@RG' | grep -o 'SM:[^[:space:]]*' | cut -c4-)
  java -Xmx8g -jar $J Mutect2 -R data/human_g1k_v37.20.21.fasta -I data/tumor_$i.bam -I data/normal_$i.bam -normal "$N" \
    -L 20 --germline-resource data/very-small-gnomad.vcf -O truth_regions/m2_$i.vcf \
    --assembly-region-out truth_regions/regions_$i.igv > truth_regions/m2_$i.log 2>&1
  echo "pair $i rc=$? regions $(grep -vc '^#' truth_regions/regions_$i.igv)"
done
