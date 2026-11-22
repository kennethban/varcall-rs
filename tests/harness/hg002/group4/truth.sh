#!/usr/bin/env bash
# GATK truth for S3 group 4, several -I in HaplotypeCaller (inputs: make_inputs.sh), HG002 chr8:1-15M (more multi-sample cases: ms/make.sh).
# Outputs truth/<case>.{vcf.gz,g.vcf.gz}; each case's arguments in truth/<case>.cmd (input paths relative to here).
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/hg002/group4" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
REF=../ref/GRCh38_no_alt_analysis_set.fa
L="-L chr8:1-15000000"
GQB="-GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90"
WV="-contamination 0 -G StandardAnnotation -G StandardHCAnnotation $GQB"
WG="-contamination 0 -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation $GQB -ERC GVCF"
mkdir -p truth
declare -A C
C[lanes_vcf]="-I in/L1.bam -I in/L2.bam $L $WV"
C[lanes_gvcf]="-I in/L1.bam -I in/L2.bam $L $WG"
C[same_rg_vcf]="-I in/S0.bam -I in/S1.bam $L"
C[three_vcf]="-I in/T0.bam -I in/T1.bam -I in/T2.bam $L"
C[three_rev_gvcf]="-I in/T2.bam -I in/T0.bam -I in/T1.bam $L -ERC GVCF"
C[contam_vcf]="-I in/L1.bam -I in/L2.bam $L -contamination 0.05"
# two samples (HG002 and OTHER: L2's reads renamed), called jointly
C[other_vcf]="-I in/L1.bam -I in/O.bam $L $WV"
run() {
  local c=$1 args=$2 ext=vcf.gz
  [[ $args == *"-ERC GVCF"* ]] && ext=g.vcf.gz
  [ -f truth/$c.$ext.tbi ] && return 0
  echo "$args" | $HARNESS/cmdargs.sh -w > truth/$c.cmd
  local t0=$(date +%s)
  java -Xmx3g -XX:ParallelGCThreads=2 -jar $J HaplotypeCaller -R $REF -O truth/$c.$ext $args \
       --native-pair-hmm-threads 2 > truth/$c.log 2>&1
  echo "$c rc=$? $(( $(date +%s) - t0 ))s" >> truth/times.txt
}
for c in "${!C[@]}"; do run $c "${C[$c]}" & done
wait
cat truth/times.txt
