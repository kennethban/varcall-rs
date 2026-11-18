#!/usr/bin/env bash
# variant-rs MarkDuplicates with the command lines of truth_warp.sh vs Picard: records (byte for byte), header (the
# VN / CL of this run's own MarkDuplicates @PG lines aside), metrics, .md5, and which index files exist.
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/markdup" || exit 1

G=${VARIANT_RS:-$REPO/target/release/variant-rs}
I=input/warp; T=truth/warp; O=${1:-rust/warp}; mkdir -p $O/tmp
hdr() { samtools view -H --no-PG $1 | awk 'BEGIN{FS=OFS="\t"} /^@PG/ && /\tPN:MarkDuplicates/ {o=$1; for (i=2; i<=NF; i++) if ($i !~ /^(VN|CL):/) o=o OFS $i; print o; next} {print}'; }
cmp1() { # name
  local n=$1
  local r=$(python3 $HARNESS/bqsr/cmpbam.py $T/$n.md.bam $O/$n.md.bam | tail -1)
  local h=$(diff <(hdr $T/$n.md.bam) <(hdr $O/$n.md.bam) > /dev/null && echo identical || echo DIFFERENT)
  local m=$(diff <(grep -v '^#' $T/$n.metrics.txt) <(grep -v '^#' $O/$n.metrics.txt) > /dev/null && echo identical || echo DIFFERENT)
  local ix="$(ls $T/$n.md.bai $T/$n.md.bam.bai 2>/dev/null | wc -l)/$(ls $O/$n.md.bai $O/$n.md.bam.bai 2>/dev/null | wc -l)"
  local md5=""; [ -f $T/$n.md.bam.md5 ] && md5="; md5 file $([ "$(cat $O/$n.md.bam.md5)" = "$(md5sum < $O/$n.md.bam | cut -c1-32)" ] && echo correct || echo WRONG)"
  echo "  $n: records $r; header $h; metrics $m; index files Picard/variant-rs $ix$md5 ($(grep -o 'total [0-9.]*s' $O/$n.log))"
}
W="VALIDATION_STRINGENCY=SILENT OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 ASSUME_SORT_ORDER=queryname CLEAR_DT=false ADD_PG_TAG_TO_READS=false"
$G MarkDuplicates INPUT=$I/qg.bam OUTPUT=$O/qg.md.bam METRICS_FILE=$O/qg.metrics.txt $W TMP_DIR=$O/tmp > $O/qg.log 2>&1; cmp1 qg
$G MarkDuplicates INPUT=$I/a.bam INPUT=$I/b.bam OUTPUT=$O/multi.md.bam METRICS_FILE=$O/multi.metrics.txt $W TMP_DIR=$O/tmp > $O/multi.log 2>&1; cmp1 multi
$G MarkDuplicates INPUT=$I/qn.bam OUTPUT=$O/qn.md.bam METRICS_FILE=$O/qn.metrics.txt TMP_DIR=$O/tmp > $O/qn.log 2>&1; cmp1 qn
$G MarkDuplicates -I $I/dt.bam -O $O/opts.md.bam -M $O/opts.metrics.txt --CLEAR_DT false --ADD_PG_TAG_TO_READS false \
  --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --CREATE_MD5_FILE true --CREATE_INDEX true --TMP_DIR $O/tmp > $O/opts.log 2>&1; cmp1 opts
$G MarkDuplicates -I $I/dt.bam -O $O/remark.md.bam -M $O/remark.metrics.txt --CREATE_INDEX true --TMP_DIR $O/tmp > $O/remark.log 2>&1; cmp1 remark
