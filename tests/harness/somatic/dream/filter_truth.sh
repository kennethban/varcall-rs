#!/usr/bin/env bash
# GATK 4.7 FilterMutectCalls (default arguments, no contamination / orientation tables) on the Mutect2 truth VCFs.
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/dream" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
for i in ${@:-1 2 3 4}; do
  t0=$(date +%s)
  java -Xmx4g -jar $J FilterMutectCalls -R data/human_g1k_v37.20.21.fasta -V truth/m2_$i.vcf -O truth/f_$i.vcf > truth/f_$i.log 2>&1
  echo "pair $i rc=$? $(( $(date +%s) - t0 ))s PASS $(grep -v '^#' truth/f_$i.vcf | awk '$7=="PASS"' | wc -l) / $(grep -vc '^#' truth/f_$i.vcf)"
done
