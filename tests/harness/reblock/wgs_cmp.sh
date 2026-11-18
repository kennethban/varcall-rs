#!/usr/bin/env bash
# variant-rs ReblockGVCF (WARP's command) on the 47 whole-genome HG002 GVCF shards vs GATK (wgs_truth.sh): data lines
# and headers (but ReblockGVCF's own command line).  usage: wgs_cmp.sh [OUTDIR]
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/reblock" || exit 1; D=$PWD
G=${VARIANT_RS:-$REPO/target/release/variant-rs}
REF=$VARCALL_DATA/hg002/ref/GRCh38_no_alt_analysis_set.fa
O=${1:-$D/rust_wgs}; mkdir -p $O; O=$(realpath $O)
t0=$(date +%s)
ls $VARCALL_DATA/hg002/gatk/gvcf/shard*.g.vcf.gz | xargs -P 8 -I{} sh -c '
  s=$(basename {} .g.vcf.gz)
  '$G' ReblockGVCF -R '$REF' -V {} -do-qual-approx --floor-blocks -GQB 20 -GQB 30 -GQB 40 -O '$O'/$s.rb.g.vcf.gz > '$O'/$s.log 2>&1 || echo "$s rc=$?"'
n=0; bad=0; recs=0
for t in wgs/shard*.rb.g.vcf.gz; do
  s=$(basename $t .rb.g.vcf.gz); o=$O/$s.rb.g.vcf.gz; n=$((n+1))
  if cmp -s <(zcat $t | grep -v '^##GATKCommandLine=<ID=ReblockGVCF') <(zcat $o | grep -v '^##variant-rsCommandLine'); then :; else bad=$((bad+1)); echo "  $s DIFFERENT"; fi
  recs=$((recs + $(zcat $t | grep -vc '^#')))
done
echo "  whole genome, $n shards: $((n-bad)) identical, $recs records ($(( $(date +%s) - t0 ))s incl. comparison)"
