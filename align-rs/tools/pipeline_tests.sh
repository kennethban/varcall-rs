#!/bin/bash
# Pipeline tests (pipeline/test) with align-rs: germline + somatic outputs vs GATK (test/compare.py), and every
# lane BAM's records + header vs bwa-mem2's digests (test/log/lanes.bwamem2.md5). Fresh run each time.
REPO=$(cd "$(dirname "$0")/../.." && pwd); R=${VARCALL_DATA:?set VARCALL_DATA to the test data directory}
cd $REPO/pipeline; export NXF_ANSI_LOG=false; bash test/render.sh || exit 1
for t in germline somatic; do
  echo "== $t"
  rm -rf test/out/$t test/work/$t
  s=$(date +%s)
  ${NEXTFLOW:-nextflow} run main.nf -params-file test/$t.yaml -w test/work/$t -with-trace test/log/$t.trace.tsv | tail -1
  echo "wall $(( $(date +%s)-s )) s"
  python3 test/compare.py $t | tail -1
  for f in test/out/$t/lanes/*.bam; do
    b=$(basename $f); exp=$(grep "^$b " test/log/lanes.bwamem2.md5)
    got="$b $(samtools view $f | md5sum | cut -c1-32) $(samtools view -H --no-PG $f | grep -v '^@PG' | md5sum | cut -c1-32)"
    [ "$got" = "$exp" ] && echo "lane BAMs identical to bwa-mem2" || echo "LANE DIFFERS: $b"
  done | sort | uniq -c
done
