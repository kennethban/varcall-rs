#!/usr/bin/env bash
# variant-rs GatherBQSRReports with each truth case's arguments (truth/<case>.cmd) vs GATK (byte for byte).
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/bqsr/gather" || exit 1; D=$PWD

G=${VARIANT_RS:-$REPO/target/release/variant-rs}
O=${1:-$D/rust}; shift; mkdir -p $O
for c in ${@:-$(ls truth/*.cmd | xargs -n1 basename | sed 's/\.cmd$//' | sort)}; do
  $G GatherBQSRReports $($HARNESS/cmdargs.sh truth/$c.cmd) -O $O/$c.table > $O/$c.log 2>&1; rc=$?
  echo "  $c: rc=$rc $(cmp -s truth/$c.table $O/$c.table && echo identical || echo "DIFFERENT ($(diff truth/$c.table $O/$c.table | grep -c '^[<>]') lines)") ($(grep -c . truth/$c.table) lines)"
done
# WARP's next step: ApplyBQSR with the gathered table
$G ApplyBQSR $(cat truth/apply.applycmd) -O $O/apply.bam > $O/apply.log 2>&1; rc=$?
echo "  ApplyBQSR with scatter3's gathered table: rc=$rc records $(python3 $HARNESS/bqsr/cmpbam.py truth/apply.bam $O/apply.bam | tail -1)"
