#!/usr/bin/env bash
# variant-rs GenotypeGVCFs (WARP's arguments) on the 47 reblocked whole-genome HG002 shards vs GATK (wgs_truth.sh).
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/genotypegvcfs" || exit 1; D=$PWD
G=${VARIANT_RS:-$REPO/target/release/variant-rs}
REF=$VARCALL_DATA/hg002/ref/GRCh38_no_alt_analysis_set.fa
DBSNP=$VARCALL_DATA/bqsr/resources/Homo_sapiens_assembly38.dbsnp138.vcf.gz
O=${1:-$D/rust_wgs}; mkdir -p $O; O=$(realpath $O)
t0=$(date +%s)
ls $VARCALL_DATA/reblock/wgs/shard*.rb.g.vcf.gz | xargs -P 8 -I{} sh -c '
  s=$(basename {} .rb.g.vcf.gz)
  '$G' GenotypeGVCFs -R '$REF' -V {} -D '$DBSNP' -G StandardAnnotation -G AS_StandardAnnotation -O '$O'/$s.vcf.gz > '$O'/$s.log 2>&1 || echo "$s rc=$?"'
n=0; bad=0; recs=0
for t in wgs/shard*.vcf.gz; do
  s=$(basename $t .vcf.gz); n=$((n+1))
  cmp -s <(zcat $t | grep -v '^##GATKCommandLine=<ID=GenotypeGVCFs') <(zcat $O/$s.vcf.gz | grep -v '^##variant-rsCommandLine') || { bad=$((bad+1)); echo "  $s DIFFERENT"; }
  recs=$((recs + $(zcat $t | grep -vc '^#')))
done
echo "  whole genome, $n shards: $((n-bad)) identical, $recs records ($(( $(date +%s) - t0 ))s incl. comparison)"
