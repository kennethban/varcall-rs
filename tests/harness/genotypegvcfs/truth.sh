#!/usr/bin/env bash
# GATK GenotypeGVCFs truth (S3 group 4), single-sample GVCF inputs (HG002 chr8:1-15M, hg002/group2 + reblock truth).
# WARP (JointGenotypingTasks.wdl GenotypeGVCFs) reads a GenomicsDB workspace; here the same arguments on GVCF files:
# warp: the HaplotypeCaller GVCF; warp_rb: the reblocked GVCF (what WARP's joint calling reads); default: no
# arguments; noas: a GVCF without AS annotations.
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/genotypegvcfs" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
REF=$VARCALL_DATA/hg002/ref/GRCh38_no_alt_analysis_set.fa
DBSNP=$VARCALL_DATA/bqsr/resources/Homo_sapiens_assembly38.dbsnp138.vcf.gz
T=$VARCALL_DATA/hg002/group2/truth; RB=$VARCALL_DATA/reblock/truth
W="-D $DBSNP -G StandardAnnotation -G AS_StandardAnnotation --only-output-calls-starting-in-intervals -L chr8:1-15000000 --merge-input-intervals"
mkdir -p truth
declare -A C
C[warp]="-V $T/warp_gvcf.g.vcf.gz $W"
C[warp_rb]="-V $RB/warp.rb.g.vcf.gz $W"
C[default]="-V $T/warp_gvcf.g.vcf.gz"
C[noas]="-V $T/contam_gvcf.g.vcf.gz"
run() {
  local c=$1
  [ -f truth/$c.vcf.gz.tbi ] && return 0
  echo "${C[$c]}" | $HARNESS/cmdargs.sh -w > truth/$c.cmd
  local t0=$(date +%s)
  java -Xmx4g -XX:ParallelGCThreads=2 -jar $J GenotypeGVCFs -R $REF ${C[$c]} -O truth/$c.vcf.gz > truth/$c.log 2>&1
  echo "$c rc=$? $(( $(date +%s) - t0 ))s" >> truth/times.txt
}
for c in "${!C[@]}"; do run $c & done
wait
cat truth/times.txt
