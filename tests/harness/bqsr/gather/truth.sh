#!/usr/bin/env bash
# GATK truth for GatherBQSRReports (S3 group 4). Inputs as WARP makes them (BamProcessing.wdl: BaseRecalibrator per
# scatter interval, then GatherBqsrReports): scatter3 = NA12878 chr20 in three intervals; lanes = HG002 chr8:1-15M
# lane BAMs of hg002/group4 (read groups L1, L2, T0: each table has its own) given out of order; single = one table;
# twice = the same table twice.
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/bqsr/gather" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
R=../resources
KS="--known-sites $R/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites $R/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --known-sites $R/Homo_sapiens_assembly38.known_indels.vcf.gz"
NA=$VARCALL_DATA/wgs/bam/NA12878.bam; NR=../ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
G4=$VARCALL_DATA/hg002/group4/in; HR=$VARCALL_DATA/hg002/ref/GRCh38_no_alt_analysis_set.fa
mkdir -p in truth
br() { # out args...
  local o=$1; shift
  [ -f in/$o.table ] || java -Xmx4g -XX:ParallelGCThreads=2 -jar $J BaseRecalibrator "$@" $KS -O in/$o.table > in/$o.log 2>&1 || echo "$o FAILED"
}
br na_a -R $NR -I $NA -L chr20:1-21000000 &
br na_b -R $NR -I $NA -L chr20:21000001-42000000 &
br na_c -R $NR -I $NA -L chr20:42000001-64444167 &
br l1 -R $HR -I $G4/L1.bam -L chr8:1-15000000 &
br l2 -R $HR -I $G4/L2.bam -L chr8:1-15000000 &
br t0 -R $HR -I $G4/T0.bam -L chr8:1-15000000 &
wait
declare -A C
C[scatter3]="-I in/na_a.table -I in/na_b.table -I in/na_c.table"
C[lanes]="-I in/t0.table -I in/l2.table -I in/l1.table"
C[single]="-I in/na_b.table"
C[mixed]="-I in/l1.table -I in/na_b.table -I in/t0.table -I in/na_c.table"
C[twice]="-I in/l1.table -I in/l1.table"
for c in "${!C[@]}"; do
  [ -f truth/$c.table ] && continue
  echo "${C[$c]}" | $HARNESS/cmdargs.sh -w > truth/$c.cmd
  java -Xmx2g -jar $J GatherBQSRReports ${C[$c]} -O truth/$c.table > truth/$c.log 2>&1; echo "$c rc=$?"
done
# WARP's next step, ApplyBQSR with a gathered table (cmp.sh runs variant-rs with truth/apply.applycmd)
NR=../ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
echo "-R $NR -I $NA -L chr20:10000000-12000000 --bqsr-recal-file truth/scatter3.table" > truth/apply.applycmd
[ -f truth/apply.bam ] || java -Xmx4g -jar $J ApplyBQSR $(cat truth/apply.applycmd) -O truth/apply.bam > truth/apply.log 2>&1
