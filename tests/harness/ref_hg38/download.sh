#!/usr/bin/env bash
# GATK GRCh38 reference (Homo_sapiens_assembly38) + prebuilt bwa-mem2 index from AWS iGenomes. Resumable.
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/ref_hg38" || exit 1
B=https://ngi-igenomes.s3.amazonaws.com/igenomes/Homo_sapiens/GATK/GRCh38/Sequence
get() { for i in 1 2 3 4 5; do curl -sSL --retry 5 -C - -o "$2" "$1" && return 0; sleep 10; done; echo "FAILED $1"; return 1; }
date "+start %H:%M"
for f in Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38.fasta.fai Homo_sapiens_assembly38.dict; do get $B/WholeGenomeFasta/$f $f; done
for f in amb ann alt pac 0123 bwt.2bit.64; do get $B/BWAmem2Index/Homo_sapiens_assembly38.fasta.$f Homo_sapiens_assembly38.fasta.$f; done
date "+done %H:%M"
ls -la
