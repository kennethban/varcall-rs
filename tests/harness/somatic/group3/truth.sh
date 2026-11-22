#!/usr/bin/env bash
# GATK truth for S3 group 3 (Mutect2 arguments). HCC1143 (hg38 chr17plus, tumor/normal, PoN, gnomAD) for the
# tumor-normal and tumor-only cases; NA12878 WGS chrM for mitochondria mode (WARP/GATK AlignAndCall.wdl command).
# Outputs truth/<case>.vcf (+ .stats, bamout); each case's arguments in truth/<case>.cmd.
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/group3" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
H=$VARCALL_DATA/somatic/hcc1143/data; HR=$VARCALL_DATA/ref_hg38/Homo_sapiens_assembly38.fasta
W=$VARCALL_DATA/wgs/bam/NA12878.bam; WR=$VARCALL_DATA/bqsr/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
TN="-R $HR -I $H/tumor.bam -I $H/normal.bam -normal HCC1143_normal -pon $H/chr17_m2pon.vcf.gz --germline-resource $H/chr17_af-only-gnomad_grch38.vcf.gz -L $H/chr17plus.interval_list"
TO="-R $HR -I $H/tumor.bam -pon $H/chr17_m2pon.vcf.gz --germline-resource $H/chr17_af-only-gnomad_grch38.vcf.gz -L $H/chr17plus.interval_list"
mkdir -p truth
declare -A C
C[gg]="$TN --genotype-germline-sites"
C[gp]="$TN --genotype-pon-sites"
C[af]="$TN --af-of-alleles-not-in-resource 0.0001"
C[to]="$TO"
C[to_gg]="$TO --genotype-germline-sites --genotype-pon-sites"
C[bamout]="$TN --bam-output truth/bamout.bam"
C[mito]="-R $WR -I $W -L chrM --mitochondria-mode"
# AlignAndCall.wdl M2 task (m2_extra_args empty, max_reads_per_alignment_start 75)
C[mito_warp]="-R $WR -I $W -L chrM --read-filter MateOnSameContigOrNoMappedMateReadFilter --read-filter MateUnmappedAndUnmappedReadFilter --annotation StrandBiasBySample --mitochondria-mode --max-reads-per-alignment-start 75 --max-mnp-distance 0"
run() {
  local c=$1 args=$2
  [ -f truth/$c.vcf.stats ] && return 0
  echo "$args" | $HARNESS/cmdargs.sh -w > truth/$c.cmd
  local t0=$(date +%s)
  java -Xmx4g -XX:ParallelGCThreads=2 -jar $J Mutect2 $args -O truth/$c.vcf --native-pair-hmm-threads 2 > truth/$c.log 2>&1
  echo "$c rc=$? $(( $(date +%s) - t0 ))s" >> truth/times.txt
}
for c in "${!C[@]}"; do run $c "${C[$c]}" & done
wait
cat truth/times.txt
