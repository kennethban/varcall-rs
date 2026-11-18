#!/bin/bash
# .alt parsing as bns_restore (review finding 26): tests/data/small.fa linked under prefixes whose .alt names the
# contig (a) with a final newline: ALT (@SQ AH:*), (b) on a last line without one: not looked up, not ALT. align-rs
# vs bwa-mem2 on the smallref reads, header and records.
# Usage: tests/alt.sh [align-rs binary] [out dir]
REPO=$(cd "$(dirname "$0")/../.." && pwd); G=${VARCALL_DATA:?set VARCALL_DATA to the test data directory}; S=$G/wgs/fastq/sub; D=$REPO/align-rs/tests/data
A=$(realpath ${1:-$REPO/target/release/align-rs}); O=${2:-/tmp/align-rs-alt}; mkdir -p $O
for c in nl nonl; do
  for x in fa fa.0123 fa.amb fa.ann fa.bwt.2bit.64 fa.pac; do ln -sf $D/small.$x $O/$c.$x; done
  [ $c = nl ] && printf '@SQ\tSN:x\nchr20_10M_15M\n' > $O/$c.fa.alt || printf '@SQ\tSN:x\nchr20_10M_15M' > $O/$c.fa.alt
  bwa-mem2 mem -t 8 -K 100000000 -Y $O/$c.fa $S/smallref.R1.fq.gz $S/smallref.R2.fq.gz 2> $O/$c.bwa.log | grep -v '^@PG' > $O/$c.bwa.sam
  $A mem -t 8 -K 100000000 -Y $O/$c.fa $S/smallref.R1.fq.gz $S/smallref.R2.fq.gz 2> $O/$c.aln.log | grep -v '^@PG' > $O/$c.aln.sam
  echo "  .alt $c: bwa-mem2 @SQ $(grep -c 'AH:\*' $O/$c.bwa.sam) ALT; align-rs $(cmp -s $O/$c.bwa.sam $O/$c.aln.sam && echo identical || echo DIFFERENT)"
done
