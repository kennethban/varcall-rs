#!/usr/bin/env bash
# variant-rs ApplyBQSR vs GATK (truth.sh): records (bytes), the POS 0 read included
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/bqsr/pos0" || exit 1

X=${VARIANT_RS:-$REPO/target/release/variant-rs}; REF=$VARCALL_DATA/ref_hg38/Homo_sapiens_assembly38.fasta
O=${1:-/tmp/bqsr-pos0}; mkdir -p $O
$X ApplyBQSR -R $REF -I ../../markdup/edge/in/pos0.bam --bqsr-recal-file truth/pos0.table -O $O/pos0.bam > $O/apply.log 2>&1
echo "  ApplyBQSR, POS 0 read: rc=$? records $(python3 $HARNESS/bqsr/cmpbam.py truth/pos0.bam $O/pos0.bam 2>&1 | tail -1)"
