#!/usr/bin/env bash
# Picard 3.5.0 (GATK 4.7 bundle) MarkDuplicates truth for WARP's arguments (tasks/wdl/BamProcessing.wdl) and the
# other newly supported options. Inputs (input/warp/):
#   qg.bam      NA12878 chr20, query-grouped (samtools collate), as MergeBamAlignment output
#   a.bam b.bam qg.bam split by read name; b's read group renamed NA12878.b (same library) and an extra @PG
#               ID:samtools given another command line (ID collision with a's)
#   qn.bam      qg.bam sorted by Picard SortSam (header SO:queryname, htsjdk's order)
#   dt.bam      smoke.bam marked with TAGGING_POLICY=All (DT tags, duplicate flags, PG tags)
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/markdup" || exit 1

G="java -Xmx16g -jar ${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}"
I=input/warp; T=truth/warp; mkdir -p $I $T $T/tmp
t() { local n=$1; shift; local t0=$(date +%s); "$@" > $T/$n.log 2>&1; echo "$n rc=$? $(( $(date +%s) - t0 ))s"; }
if [ ! -f $I/qn.bam ]; then
  [ -f $I/qg.bam ] || samtools collate -@ 8 -o $I/qg.bam input/na12878.chr20.bam $I/collate_tmp
  samtools view -h $I/qg.bam | awk '/^@/ || $1 ~ /[02468]$/' | samtools view -b -o $I/a.bam -
  samtools view -h $I/qg.bam | awk 'BEGIN{OFS="\t"}
    /^@RG/ {sub(/ID:NA12878/, "ID:NA12878.b")} /^@PG\tID:samtools\t/ {sub(/CL:.*/, "CL:samtools view (split b)")} /^@/ {print; next}
    $1 ~ /[13579]$/ {for (i=12; i<=NF; i++) if ($i == "RG:Z:NA12878") $i = "RG:Z:NA12878.b"; print}' | samtools view -b -o $I/b.bam -
  $G SortSam -I $I/qg.bam -O $I/qn.bam --SORT_ORDER queryname --TMP_DIR $T/tmp > $I/qn.log 2>&1 # htsjdk's queryname order
  [ -f $I/dt.bam ] || $G MarkDuplicates -I input/smoke.bam -O $I/dt.bam -M $I/dt.metrics.txt --TAGGING_POLICY All --CREATE_INDEX true > $I/dt.log 2>&1
fi
for f in qg a b qn; do echo "$f: $(samtools view -H $I/$f.bam | grep -m1 '^@HD') $(samtools view -c $I/$f.bam) records"; done
W="VALIDATION_STRINGENCY=SILENT OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 ASSUME_SORT_ORDER=queryname CLEAR_DT=false ADD_PG_TAG_TO_READS=false"
[ -f $T/qg.md.bam ] || t qg $G MarkDuplicates INPUT=$I/qg.bam OUTPUT=$T/qg.md.bam METRICS_FILE=$T/qg.metrics.txt $W TMP_DIR=$T/tmp
t multi $G MarkDuplicates INPUT=$I/a.bam INPUT=$I/b.bam OUTPUT=$T/multi.md.bam METRICS_FILE=$T/multi.metrics.txt $W TMP_DIR=$T/tmp
t qn $G MarkDuplicates INPUT=$I/qn.bam OUTPUT=$T/qn.md.bam METRICS_FILE=$T/qn.metrics.txt TMP_DIR=$T/tmp
[ -f $T/opts.md.bam ] || t opts $G MarkDuplicates -I $I/dt.bam -O $T/opts.md.bam -M $T/opts.metrics.txt --CLEAR_DT false --ADD_PG_TAG_TO_READS false \
  --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --CREATE_MD5_FILE true --CREATE_INDEX true --TMP_DIR $T/tmp
[ -f $T/remark.md.bam ] || t remark $G MarkDuplicates -I $I/dt.bam -O $T/remark.md.bam -M $T/remark.metrics.txt --CREATE_INDEX true --TMP_DIR $T/tmp
