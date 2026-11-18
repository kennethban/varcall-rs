#!/usr/bin/env bash
# Picard SortSam (GATK 4.7's bundled Picard) truth. warp: BamProcessing.wdl's command on the query-grouped
# MarkDuplicates output of markdup/truth_warp.sh (WARP's chain); qn: same on the queryname-sorted one; sorted: an
# already coordinate-sorted BAM (NA12878 chr20, with unmapped reads); spill: warp with MAX_RECORDS_IN_RAM=100000.
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/sortsam" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
M=$VARCALL_DATA/markdup
mkdir -p truth/tmp
declare -A C
C[warp]="INPUT=$M/truth/warp/qg.md.bam SORT_ORDER=coordinate CREATE_INDEX=true CREATE_MD5_FILE=true MAX_RECORDS_IN_RAM=300000"
C[qn]="INPUT=$M/truth/warp/qn.md.bam SORT_ORDER=coordinate CREATE_INDEX=true"
C[sorted]="INPUT=$M/input/na12878.chr20.bam SORT_ORDER=coordinate"
C[spill]="INPUT=$M/truth/warp/qg.md.bam SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=100000"
# an htslib-written input with unplaced reads (bin 4680 there; htsjdk writes 0): REVIEW_FINDINGS.md, 2026-11-17
C[unplaced]="INPUT=$M/edge/in/unmapped.bam SORT_ORDER=coordinate"
run() {
  local c=$1 args=$2
  [ -f truth/$c.bam ] && return 0
  echo "$args" > truth/$c.cmd
  local t0=$(date +%s)
  java -Dsamjdk.compression_level=2 -Xmx4g -jar $J SortSam $args OUTPUT=truth/$c.bam TMP_DIR=truth/tmp > truth/$c.log 2>&1
  echo "$c rc=$? $(( $(date +%s) - t0 ))s" >> truth/times.txt
}
for c in "${!C[@]}"; do run $c "${C[$c]}" & done
wait
cat truth/times.txt
