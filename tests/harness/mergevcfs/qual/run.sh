#!/usr/bin/env bash
# MergeVcfs on QUALs at %.2f rounding edges and the VCF 3 missing value (review finding 36): Picard writes QUAL after
# htsjdk's QUAL / -10 * -10 round trip (0.105 -> 0.10), and '.' for -1 / -1.0 (and -10). variant-rs vs Picard.
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/mergevcfs/qual" || exit 1

J="java -Xmx2g -jar ${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}"; X=${VARIANT_RS:-$REPO/target/release/variant-rs}
O=${1:-/tmp/mergevcfs-qual}; mkdir -p $O
zcat ../in/v1.vcf.gz | awk -F'\t' 'BEGIN{OFS="\t"; split("0.055 0.105 0.225 0.115 1.005 2.675 -1 -1.0 -10 12.345 0.125 3.0001", q, " ")}
  /^#/ {print; next} { if (++n <= 12) $6 = q[n]; print }' | bgzip > $O/q.vcf.gz && tabix -f $O/q.vcf.gz
$J MergeVcfs -I $O/q.vcf.gz -O $O/picard.vcf.gz > $O/picard.log 2>&1; p=$?
$X MergeVcfs -I $O/q.vcf.gz -O $O/rs.vcf.gz > $O/rs.log 2>&1; r=$?
d=$(cmp -s <(zcat $O/picard.vcf.gz | grep -v '^#') <(zcat $O/rs.vcf.gz | grep -v '^#') && echo identical || echo DIFFERENT)
echo "  QUAL edges: Picard rc=$p variant-rs rc=$r, records $d ($(zcat $O/picard.vcf.gz | grep -v '^#' | head -12 | cut -f6 | paste -sd' ' -))"
