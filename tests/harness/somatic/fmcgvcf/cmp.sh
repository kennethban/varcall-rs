#!/usr/bin/env bash
# variant-rs FilterMutectCalls vs GATK truth (truth.sh): VCF records and filteringStats.tsv
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/fmcgvcf" || exit 1

X=${VARIANT_RS:-$REPO/target/release/variant-rs}; O=${1:-out}; mkdir -p $O
REF=$VARCALL_DATA/ref_hg38/Homo_sapiens_assembly38.fasta
for c in plain ob; do
  e=""; [ $c = ob ] && e="--ob-priors priors.tar.gz"
  $X FilterMutectCalls -R $REF -V in.g.vcf $e -O $O/$c.vcf > $O/$c.log 2>&1; rc=$?
  n=$(grep -v '^#' truth/$c.vcf | paste - <(grep -v '^#' $O/$c.vcf 2>/dev/null) | awk -F'\t' '{h=NF/2; for(i=1;i<=h;i++) if($i!=$(i+h)) {d++; break}} END{print d+0}')
  echo "  $c: rc=$rc records $(grep -vc '^#' truth/$c.vcf) vs $(grep -vc '^#' $O/$c.vcf 2>/dev/null), $n differ; filteringStats $(cmp -s truth/$c.vcf.filteringStats.tsv $O/$c.vcf.filteringStats.tsv && echo identical || echo DIFFERENT)"
done
