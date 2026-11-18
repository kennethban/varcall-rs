#!/bin/bash
# Malformed FASTQ records (review finding 10) vs bwa-mem2 on tests/data/small.fa: kseq returns -2 for a record whose
# quality line is missing or of another length; bwa-mem2 drops it, ends the batch there and goes on after it.
#  trunc_se  single-end, the last record cut inside its quality line
#  plus_se   single-end, the last record cut right after '+' (no newline)
#  long_se   single-end, record 101 with one extra quality character (reading goes on after it)
#  trunc_pe  paired, R2 cut inside its last quality line ("the 2nd file has fewer sequences")
#  emptyq_se single-end, an empty read (record 1) with a 3-character quality line: kseq reads at least one quality
#            line, finds the lengths differ (-2) and drops the record
#  bad2_se   single-end, records 1 and 2 both malformed (empty read with qualities, a short quality line)
#  short_pe  paired, R2 without its last 50 records ("the 2nd file has fewer sequences", repeatedly)
#  gzpad_pe  paired .gz files with bytes after the last gzip member: R1 zero padding, R2 a stray newline (zlib
#            ignores them; review finding 37)
#  bom_se    single-end, a UTF-8 BOM before the first '@' (kseq skips bytes up to '@'; review finding 38)
#  ws_se     single-end, a space inside record 4's sequence line (and one more quality): kseq keeps it, an N
#  crlf_se   single-end, CRLF line ends: '\r' dropped from lines longer than 1, the empty read becomes '\r'
#  junk_se   single-end, bytes before record 6's '@' on its header line
#  latin1_pe paired, a Latin-1 byte (0xE9) in every 10th read name and quality bytes >= 0x80 (review finding 27:
#            SAM text is bytes; bwa-mem2 copies them verbatim)
# Usage: tests/badfq.sh [align-rs binary] [out dir]
REPO=$(cd "$(dirname "$0")/../.." && pwd); G=${VARCALL_DATA:?set VARCALL_DATA to the test data directory}; S=$G/wgs/fastq/sub; REF=$REPO/align-rs/tests/data/small.fa
A=$(realpath ${1:-$REPO/target/release/align-rs}); O=${2:-/tmp/align-rs-badfq}; mkdir -p $O
zcat $S/smallref.R1.fq.gz | head -800 > $O/r1.fq; zcat $S/smallref.R2.fq.gz | head -800 > $O/r2.fq
head -c -40 $O/r1.fq > $O/trunc_se.fq
{ head -798 $O/r1.fq; printf '+'; } > $O/plus_se.fq
awk 'NR == 404 { $0 = $0 "I" } { print }' $O/r1.fq > $O/long_se.fq
head -c -40 $O/r2.fq > $O/trunc_pe.R2.fq
awk 'NR == 4 { $0 = "III" } { print }' $O/r1.fq > $O/emptyq_se.fq
awk 'NR == 4 { $0 = "III" } NR == 8 { $0 = substr($0, 1, 30) } { print }' $O/r1.fq > $O/bad2_se.fq
head -600 $O/r2.fq > $O/short_pe.R2.fq
{ head -400 $O/r1.fq | gzip; tail -n +401 $O/r1.fq | gzip; head -c 512 /dev/zero; } > $O/gzpad_pe.R1.fq.gz
{ gzip < $O/r2.fq; echo; } > $O/gzpad_pe.R2.fq.gz
{ printf '\357\273\277'; cat $O/r1.fq; } > $O/bom_se.fq
awk 'NR == 14 { $0 = substr($0, 1, 10) " " substr($0, 11) } NR == 16 { $0 = $0 "I" } { print }' $O/r1.fq > $O/ws_se.fq
sed 's/$/\r/' $O/r1.fq > $O/crlf_se.fq
awk 'NR == 21 { $0 = "xx" $0 } { print }' $O/r1.fq > $O/junk_se.fq
for r in 1 2; do LC_ALL=C awk 'NR % 4 == 2 { l = length($0) } NR % 40 == 5 { sub(/^@[^ \/]*/, "&\351") } NR % 40 == 8 && l > 22 { $0 = substr($0, 1, 20) "\351\200" substr($0, 23) } { print }' $O/r$r.fq > $O/latin1_pe.R$r.fq; done
for c in trunc_se plus_se long_se emptyq_se bad2_se bom_se ws_se crlf_se junk_se trunc_pe short_pe gzpad_pe latin1_pe; do
  in=$O/$c.fq; [ $c = trunc_pe ] && in="$O/r1.fq $O/trunc_pe.R2.fq"; [ $c = latin1_pe ] && in="$O/latin1_pe.R1.fq $O/latin1_pe.R2.fq"; [ $c = short_pe ] && in="$O/r1.fq $O/short_pe.R2.fq"; [ $c = gzpad_pe ] && in="$O/gzpad_pe.R1.fq.gz $O/gzpad_pe.R2.fq.gz"
  bwa-mem2 mem -t 4 -K 20000 -Y $REF $in 2> $O/$c.bwa.log | grep -av '^@PG' > $O/$c.bwa.sam
  $A mem -t 4 -K 20000 -Y $REF $in 2> $O/$c.aln.log | grep -av '^@PG' > $O/$c.aln.sam; rc=${PIPESTATUS[0]}
  n=$(grep -avc '^@' $O/$c.bwa.sam)
  echo "  $c: align-rs exit $rc, $(cmp -s $O/$c.bwa.sam $O/$c.aln.sam && echo "identical to bwa-mem2 ($n records)" || echo "DIFFERENT ($n vs $(grep -avc '^@' $O/$c.aln.sam) records)")"
done
