#!/usr/bin/env bash
# variant-rs FilterMutectCalls vs GATK (truth.sh): records and filteringStats.tsv
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/fmcnan" || exit 1

X=${VARIANT_RS:-$REPO/target/release/variant-rs}; REF=$VARCALL_DATA/ref_hg38/Homo_sapiens_assembly38.fasta
O=${1:-/tmp/fmc-nan}; mkdir -p $O
$X FilterMutectCalls -R $REF -V in.vcf -O $O/out.vcf > $O/fmc.log 2>&1; rc=$?
echo "  AD 0,0 record: rc=$rc records $(cmp -s <(grep -v '^##' truth/out.vcf) <(grep -v '^##' $O/out.vcf 2>/dev/null) && echo identical || echo DIFFERENT), filteringStats $(cmp -s truth/out.vcf.filteringStats.tsv $O/out.vcf.filteringStats.tsv && echo identical || echo DIFFERENT)"
