#!/bin/bash
# Pipeline edge cases (review findings 12, 13): each run must fail, early and with the reason.
#  dup:      a samplesheet with sample NA12878 lane L1 on two rows (lane BAM, @RG ID and PU would collide)
#  names:    aligner bwa-mem2, lane L2's R2 reads renamed after the first 10,000: bwa-mem2 exits ("paired reads have
#            different names") inside `bwa-mem2 | samtools sort`, whose partial BAM must not pass as the lane's
#            (a truncated .gz does not do: bwa-mem2 warns "the 2nd file has fewer sequences" and exits 0, as for a
#            shorter R2)
# Usage: test/edge.sh [main.nf]
REPO=$(cd "$(dirname "$0")/../.." && pwd); R=${VARCALL_DATA:?set VARCALL_DATA to the test data directory}
cd $REPO/pipeline; export NXF_ANSI_LOG=false; bash test/render.sh || exit 1
M=${1:-main.nf}; E=test/work/edge; rm -rf $E; mkdir -p $E
head -1 test/germline.csv > $E/dup.csv; sed -n 2p test/germline.csv >> $E/dup.csv; sed -n 2p test/germline.csv >> $E/dup.csv
zcat test/data/NA12878_L2_R2.fq.gz | awk 'NR % 4 == 1 && NR > 40000 { $1 = $1 "x" } { print }' | gzip > $E/L2_R2.fq.gz
sed "s|$R/pipeline/nf/test/data/NA12878_L2_R2.fq.gz|$REPO/pipeline/$E/L2_R2.fq.gz|" test/germline.csv > $E/names.csv
for c in dup names; do
  case $c in dup) a="--input $E/dup.csv"; want='lane L1 on more than one row';; names) a="--input $E/names.csv --aligner bwa-mem2"; want='Process `ALIGN';; esac
  ${NEXTFLOW:-nextflow} run $M -params-file test/germline.yaml $a --outdir $E/out.$c -w $E/work.$c > $E/$c.log 2>&1; rc=$?
  if [ $rc != 0 ] && grep -q "$want" $E/$c.log; then echo "$c: failed as it should ($(grep -m1 -oE "$want[^\"]*" $E/$c.log | cut -c1-90))"
  else echo "$c: NOT FAILED (exit $rc; completed: $(grep -oE '\[[^]]+\] Submitted process > [A-Z_]+' $E/$c.log | awk '{print $NF}' | sort -u | tr '\n' ' '))"; fi
done
