#!/usr/bin/env bash
# GenotypeGVCFs --scatter-intervals truth: GATK run once per (merged) -L interval, as nf-core/sarek runs it per
# interval, the runs' records in interval order under the first run's header. On the two workspaces with several
# arrays: ws_sarek (sarek's BED: 4 intervals, 3 after merging) and ws_wg (WARP's arguments, chr20-22).
# The case files truth/<ws>.scatter.cmd are cmp.sh's (variant-rs with --scatter-intervals).
HARNESS=$(cd "$(dirname "$0")/.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/genomicsdb" || exit 1
J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
mkdir -p truth/scatter
# case workspace reference "other arguments" interval...
one() {
  local c=$1 ws=$2 ref=$3 args=$4; shift 4
  [ -f truth/$c.vcf.gz.tbi ] && return
  local i=0 parts=()
  for iv in "$@"; do
    i=$((i + 1))
    java -Xmx8g -jar $J GenotypeGVCFs -R $ref -V gendb://$ws $args -L $iv -O truth/scatter/$c.$i.vcf.gz > truth/scatter/$c.$i.log 2>&1 \
      || { echo "$c interval $iv FAILED"; return 1; }
    parts+=(truth/scatter/$c.$i.vcf.gz)
  done
  { zcat ${parts[0]} | grep '^#'; for p in "${parts[@]}"; do zcat $p | grep -v '^#'; done; } | bgzip > truth/$c.vcf.gz
  tabix -p vcf truth/$c.vcf.gz
  echo "$c: $i GATK runs, $(zcat truth/$c.vcf.gz | grep -vc '^#') records"
}
B=$VARCALL_DATA/bqsr/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
H=$VARCALL_DATA/hg002/ref/GRCh38_no_alt_analysis_set.fa
WARP="-D $VARCALL_DATA/bqsr/resources/Homo_sapiens_assembly38.dbsnp138.vcf.gz -G StandardAnnotation -G AS_StandardAnnotation --only-output-calls-starting-in-intervals --merge-input-intervals"
echo "-V gendb://ws_sarek -L sarek.bed --scatter-intervals" > truth/ws_sarek.scatter.cmd
echo "$WARP -V gendb://ws_wg -L chr20 -L chr21 -L chr22 --scatter-intervals" > truth/ws_wg.scatter.cmd
# sarek.bed (0-based) merged: chr8:100001-2500000 (two overlapping lines), chr8:5000001-5600000, chr8:9000001-9100000
one ws_sarek.scatter ws_sarek $B "" chr8:100001-2500000 chr8:5000001-5600000 chr8:9000001-9100000 &
one ws_wg.scatter ws_wg $H "$WARP" chr20 chr21 chr22 &
wait
