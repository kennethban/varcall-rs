#!/usr/bin/env bash
# FilterMutectCalls with one sample's --contamination-table (or --tumor-segmentation) given twice (review finding 35):
# GATK fails (Collectors.toMap: Duplicate key); so must variant-rs, instead of using one of them.
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/fmcdup" || exit 1

J="java -Xmx4g -jar ${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}"; X=${VARIANT_RS:-$REPO/target/release/variant-rs}
O=${1:-/tmp/fmc-dup}; mkdir -p $O; C=../contam/truth; V=../dream/truth/m2_1.vcf; R="-R ../dream/data/human_g1k_v37.20.21.fasta"
cp $C/dream_1.contamination.table $O/c2.table; cp $C/dream_1.segments.table $O/s2.table
for c in contamination segmentation; do
  case $c in contamination) a="--contamination-table $C/dream_1.contamination.table --contamination-table $O/c2.table";;
    segmentation) a="--tumor-segmentation $C/dream_1.segments.table --tumor-segmentation $O/s2.table";; esac
  $J FilterMutectCalls $R -V $V $a -O $O/$c.gatk.vcf > $O/$c.gatk.log 2>&1; g=$?
  $X FilterMutectCalls $R -V $V $a -O $O/$c.rs.vcf > $O/$c.rs.log 2>&1; r=$?
  why=$(grep -h -o "Duplicate key [^ ]*" $O/$c.gatk.log $O/$c.rs.log | sort -u | paste -sd" " -)
  [ $g != 0 ] && [ $r != 0 ] && echo "  $c table twice: both fail ($why): same" || echo "  $c table twice: GATK exit $g, variant-rs exit $r: DIFFERENT"
done
