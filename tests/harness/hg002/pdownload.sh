#!/usr/bin/env bash
# parallel ranged download of the HG002 BAM (8 connections), resumable per part
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/hg002" || exit 1
U=https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam
S=46033519192; N=8; C=$(( (S + N - 1) / N ))
date "+start %H:%M"
for i in $(seq 0 $((N-1))); do
  a=$((i*C)); b=$(( (i+1)*C - 1 )); [ $b -ge $S ] && b=$((S-1))
  ( want=$((b-a+1)); for k in 1 2 3 4 5 6 7 8; do
      have=$(stat -c %s bam/parts/p$i 2>/dev/null || echo 0); [ $have -ge $want ] && break
      curl -sS --retry 5 -r $((a+have))-$b $U >> bam/parts/p$i; done ) &
done
wait
tot=0; for i in $(seq 0 $((N-1))); do tot=$((tot + $(stat -c %s bam/parts/p$i))); done
[ $tot -eq $S ] || { echo "size mismatch $tot != $S"; exit 1; }
cat $(for i in $(seq 0 $((N-1))); do echo bam/parts/p$i; done) > bam/HG002.bam.tmp && mv bam/HG002.bam.tmp bam/HG002.bam && rm -r bam/parts

samtools quickcheck -v bam/HG002.bam && echo "quickcheck ok"
date "+bam done %H:%M"
