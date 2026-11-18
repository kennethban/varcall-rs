#!/usr/bin/env bash
# GATK ReblockGVCF with WARP's command on the whole-genome HG002 GVCF (hg002/gatk/gvcf: GATK HaplotypeCaller
# -ERC GVCF, 47 shards). Outputs wgs/<shard>.rb.g.vcf.gz; compare with wgs_cmp.sh.
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/reblock" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
REF=$VARCALL_DATA/hg002/ref/GRCh38_no_alt_analysis_set.fa
mkdir -p wgs
ls $VARCALL_DATA/hg002/gatk/gvcf/shard*.g.vcf.gz | xargs -P 8 -I{} sh -c '
  s=$(basename {} .g.vcf.gz); [ -f wgs/$s.rb.g.vcf.gz.tbi ] && exit 0
  java -Xmx2g -XX:ParallelGCThreads=1 -jar '$J' ReblockGVCF -R '$REF' -V {} -do-qual-approx --floor-blocks -GQB 20 -GQB 30 -GQB 40 -O wgs/$s.rb.g.vcf.gz > wgs/$s.log 2>&1 || echo "$s FAILED"'
ls wgs/*.rb.g.vcf.gz | wc -l
