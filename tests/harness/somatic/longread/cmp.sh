#!/usr/bin/env bash
# variant-rs vs GATK truth (truth.sh): contig split into chunks (VARIANT_RS_CONTIG_THREADS=4) or not, front end on one or two threads
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/longread" || exit 1

X=${VARIANT_RS:-$REPO/target/release/variant-rs}; O=${1:-out}; mkdir -p $O
REF=$VARCALL_DATA/ref_hg38/Homo_sapiens_assembly38.fasta
for ct in 1 4; do for ts in 0 1; do
  VARIANT_RS_CONTIG_THREADS=$ct VARIANT_RS_TWO_STAGE=$ts $X Mutect2 -R $REF -I tumor.bam -L ivs.bed -O $O/out.vcf > $O/m2.$ct.$ts.log 2>&1; rc=$?
  a=$(grep -v '^##' truth/out.vcf | md5sum | cut -c1-32); b=$(grep -v '^##' $O/out.vcf 2>/dev/null | md5sum | cut -c1-32)
  s=$(cmp -s truth/out.vcf.stats $O/out.vcf.stats && echo identical || echo DIFFERENT)
  echo "  chunks-threads=$ct two-stage=$ts: rc=$rc records $( [ $a = $b ] && echo "identical ($(grep -vc '^#' truth/out.vcf))" || echo DIFFERENT), stats $s"
done; done
