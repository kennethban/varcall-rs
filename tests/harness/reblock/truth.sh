#!/usr/bin/env bash
# GATK ReblockGVCF truth (S3 group 4). Inputs: hg002/group2 HaplotypeCaller GVCFs (WARP's GVCF command = WARP's
# Reblock input after MergeVCFs). warp / warp_shard: GermlineVariantDiscovery.wdl's Reblock command on chr8:1-15M /
# all of shard023; noas: a GVCF without AS annotations (-contamination 0.05); dbsnp: one with dbSNP IDs; default:
# ReblockGVCF's default arguments (no floor blocks, bands 20 / 100, no QUALapprox).
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/reblock" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
REF=$VARCALL_DATA/hg002/ref/GRCh38_no_alt_analysis_set.fa
T=$VARCALL_DATA/hg002/group2/truth
W="-do-qual-approx --floor-blocks -GQB 20 -GQB 30 -GQB 40"
mkdir -p truth
declare -A C
C[warp]="-V $T/warp_gvcf.g.vcf.gz $W"
C[warp_shard]="-V $T/warp_gvcf_shard023.g.vcf.gz $W"
C[noas]="-V $T/contam_gvcf.g.vcf.gz $W"
C[dbsnp]="-V $T/dbsnp_gvcf.g.vcf.gz $W"
C[default]="-V $T/warp_gvcf.g.vcf.gz"
run() {
  local c=$1
  [ -f truth/$c.rb.g.vcf.gz.tbi ] && return 0
  echo "${C[$c]}" > truth/$c.cmd
  local t0=$(date +%s)
  java -Xmx3g -XX:ParallelGCThreads=2 -jar $J ReblockGVCF -R $REF ${C[$c]} -O truth/$c.rb.g.vcf.gz > truth/$c.log 2>&1
  echo "$c rc=$? $(( $(date +%s) - t0 ))s" >> truth/times.txt
}
for c in "${!C[@]}"; do run $c & done
wait
cat truth/times.txt
