#!/usr/bin/env bash
# variant-rs GatherBamFiles with each truth case's arguments vs Picard: records byte for byte, header, which bytes differ
# (only the re-deflated first block of each further input may), .bai answers, .md5.  usage: cmp.sh [OUTDIR] [case...]
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/gatherbam" || exit 1; D=$PWD

G=${VARIANT_RS:-$REPO/target/release/variant-rs}
O=${1:-$D/rust}; shift; mkdir -p $O; O=$(realpath $O)
for c in ${@:-$(ls truth/*.cmd | xargs -n1 basename | sed 's/\.cmd$//' | sort)}; do
  $G GatherBamFiles $(cat truth/$c.cmd) OUTPUT=$O/$c.bam > $O/$c.log 2>&1; rc=$?
  r=$(python3 $HARNESS/bqsr/cmpbam.py truth/$c.bam $O/$c.bam | tail -1)
  h=$(diff <(samtools view -H --no-PG truth/$c.bam) <(samtools view -H --no-PG $O/$c.bam) > /dev/null && echo identical || echo DIFFERENT)
  bytes=$(cmp -s truth/$c.bam $O/$c.bam && echo "file identical" || echo "file differs: $(python3 $HARNESS/gatherbam/blocks.py truth/$c.bam $O/$c.bam), $(( $(grep -o 'INPUT=' truth/$c.cmd | wc -l) - 1 )) re-deflated")
  ix="index Picard/variant-rs $(ls truth/$c.bai 2>/dev/null | wc -l)/$(ls $O/$c.bai 2>/dev/null | wc -l)"
  if [ -f truth/$c.bai ]; then
    q() { samtools idxstats -X $1 $2; for r in chr20:10500000-10600000 chr21:21000000-21000500 chr22 chr20:11999000-12100000; do samtools view -c -X $1 $2 $r; done; }
    ix="$ix ($(diff <(q truth/$c.bam truth/$c.bai) <(q $O/$c.bam $O/$c.bai) > /dev/null && echo "same answers" || echo "DIFFERENT answers"))"
  fi
  md5=""; [ -f truth/$c.bam.md5 ] && md5="; md5 file $([ "$(cat $O/$c.bam.md5)" = "$(md5sum < $O/$c.bam | cut -c1-32)" ] && echo correct || echo WRONG)"
  echo "  $c: rc=$rc records $r; header $h; $bytes; $ix$md5 $(grep -o '([0-9.]*s)' $O/$c.log)"
done
