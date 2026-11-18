#!/usr/bin/env bash
# variant-rs SortSam with each truth case's arguments (truth/<case>.cmd) vs Picard: records byte for byte, header,
# .bai (same idxstats and region counts through each file's own index; the bytes differ with the BGZF blocks), .md5.
# usage: cmp.sh [OUTDIR] [case...]
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/sortsam" || exit 1; D=$PWD

G=${VARIANT_RS:-$REPO/target/release/variant-rs}
O=${1:-$D/rust}; shift; mkdir -p $O/tmp; O=$(realpath $O)
cases=${@:-$(ls truth/*.cmd | xargs -n1 basename | sed 's/\.cmd$//' | sort)}
for c in $cases; do
  /usr/bin/time -f "%e s" $G SortSam $(cat truth/$c.cmd) OUTPUT=$O/$c.bam TMP_DIR=$O/tmp > $O/$c.log 2>&1
  rc=$?
  r=$(python3 $HARNESS/bqsr/cmpbam.py truth/$c.bam $O/$c.bam | tail -1)
  h=$(diff <(samtools view -H --no-PG truth/$c.bam) <(samtools view -H --no-PG $O/$c.bam) > /dev/null && echo identical || echo DIFFERENT)
  ix="index files Picard/variant-rs $(ls truth/$c.bai 2>/dev/null | wc -l)/$(ls $O/$c.bai 2>/dev/null | wc -l)"
  if [ -f truth/$c.bai ]; then
    q() { samtools idxstats -X $1 $2; for r in chr20:1-100000 chr20:30000000-30100000 chr1:1000000-2000000 chrX chr20:60000000-; do samtools view -c -X $1 $2 $r 2>/dev/null; done; }
    ix="$ix ($(diff <(q truth/$c.bam truth/$c.bai) <(q $O/$c.bam $O/$c.bai) > /dev/null && echo "same answers" || echo "DIFFERENT answers"))"
  fi
  md5=""; [ -f truth/$c.bam.md5 ] && md5="; md5 file $([ "$(cat $O/$c.bam.md5)" = "$(md5sum < $O/$c.bam | cut -c1-32)" ] && echo correct || echo WRONG)"
  echo "  $c: rc=$rc records $r; header $h; $ix$md5 ($(tail -1 $O/$c.log))"
done
