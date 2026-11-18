#!/usr/bin/env bash
# variant-rs MarkDuplicates with each case's command (truth/<case>.cmd, see make.sh) vs Picard: records, metrics, header.
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/markdup/multi" || exit 1

G=${VARIANT_RS:-$REPO/target/release/variant-rs}
O=${1:-rust}; shift; mkdir -p $O
cases=${@:-$(ls truth/*.cmd | xargs -n1 basename | sed 's/\.cmd$//' | sort)}
for c in $cases; do
  t0=$(date +%s.%N)
  $G --java-options "-Xmx4g -XX:-UsePerfData" MarkDuplicates $(cat truth/$c.cmd) --OUTPUT $O/$c.bam --METRICS_FILE $O/$c.metrics > $O/$c.log 2>&1
  rc=$?
  t=$(echo "$(date +%s.%N) - $t0" | bc | xargs printf '%.1f')
  rec=$(python3 $HARNESS/bqsr/cmpbam.py truth/$c.bam $O/$c.bam 2>&1 | tail -1)
  met=$(diff -q <(grep -v '^#' truth/$c.metrics) <(grep -v '^#' $O/$c.metrics) > /dev/null && echo identical || echo DIFFERENT)
  hdr=$(diff <(samtools view --no-PG -H truth/$c.bam | grep -v "ID:MarkDuplicates") <(samtools view --no-PG -H $O/$c.bam | grep -v "ID:MarkDuplicates") > /dev/null && echo identical || echo DIFFERENT)
  echo "  $c: rc=$rc records: $rec; metrics: $met; header: $hdr ($t s)"
done
