#!/usr/bin/env bash
# variant-rs vs GATK truth (truth.sh), with the front end on one thread and on two (VARIANT_RS_TWO_STAGE)
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/hclong" || exit 1

X=${VARIANT_RS:-$REPO/target/release/variant-rs}; O=${1:-out}; mkdir -p $O
REF=../hg002/ref/GRCh38_no_alt_analysis_set.fa
for ts in 0 1; do for m in vcf gvcf; do
  e=""; [ $m = gvcf ] && e="-ERC GVCF"
  VARIANT_RS_TWO_STAGE=$ts $X HaplotypeCaller -R $REF -I long.bam -L chr20:10000000-10400000 -O $O/out.$m.gz $e > $O/$m.$ts.log 2>&1; rc=$?
  a=$(zcat truth/out.$m.gz | grep -v '^##' | md5sum | cut -c1-32); b=$(zcat $O/out.$m.gz 2>/dev/null | grep -v '^##' | md5sum | cut -c1-32)
  echo "  $m two-stage=$ts: rc=$rc records $( [ $a = $b ] && echo "identical ($(zcat truth/out.$m.gz | grep -vc '^#'))" || echo DIFFERENT)"
done; done
