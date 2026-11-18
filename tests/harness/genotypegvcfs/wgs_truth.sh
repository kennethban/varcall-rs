#!/usr/bin/env bash
# GATK GenotypeGVCFs with WARP's arguments (JointGenotypingTasks.wdl, a GVCF file in place of the GenomicsDB
# workspace; no -L) on the reblocked whole-genome HG002 GVCF (reblock/wgs: 47 shards). Outputs wgs/<shard>.vcf.gz.
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/genotypegvcfs" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
REF=$VARCALL_DATA/hg002/ref/GRCh38_no_alt_analysis_set.fa
DBSNP=$VARCALL_DATA/bqsr/resources/Homo_sapiens_assembly38.dbsnp138.vcf.gz
mkdir -p wgs
ls $VARCALL_DATA/reblock/wgs/shard*.rb.g.vcf.gz | xargs -P 8 -I{} sh -c '
  s=$(basename {} .rb.g.vcf.gz); [ -f wgs/$s.vcf.gz.tbi ] && exit 0
  java -Xmx2g -XX:ParallelGCThreads=1 -jar '$J' GenotypeGVCFs -R '$REF' -V {} -D '$DBSNP' -G StandardAnnotation -G AS_StandardAnnotation -O wgs/$s.vcf.gz > wgs/$s.log 2>&1 || echo "$s FAILED"'
ls wgs/*.vcf.gz | wc -l
