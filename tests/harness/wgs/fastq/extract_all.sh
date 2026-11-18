#!/usr/bin/env bash
# FASTQ pairs of every NA12878 lane from the 30x BAM (12 read groups): one pass splits the BAM by read group, then
# each lane is collated by name and written as R1/R2 (primary records only, as the two lanes extracted earlier).

set -euo pipefail
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/wgs/fastq" || exit 1
B=$VARCALL_DATA/wgs/bam/NA12878.bam
mkdir -p split
t0=$(date +%s)
[ -e split/.done ] || { samtools split -@ 16 --output-fmt BAM,level=1 -f 'split/%!.%.' $B && touch split/.done; }
echo "split $(( $(date +%s) - t0 ))s"
one() {  # read group ID
  local fl=${1##*_TTGCCTAG-ACCACTTA_}
  [ -s NA12878_${fl}_R2.fq.gz ] && { rm -f split/$1.bam; return; }
  # .part.fq.gz, not .fq.gz.part: samtools fastq gzips only names ending in .gz
  samtools collate -@ 8 -u -O -T split/$fl.tmp split/$1.bam \
    | samtools fastq -@ 8 -F 0x900 -1 NA12878_${fl}_R1.part.fq.gz -2 NA12878_${fl}_R2.part.fq.gz -0 /dev/null -s /dev/null -n - 2> $fl.log
  mv NA12878_${fl}_R1.part.fq.gz NA12878_${fl}_R1.fq.gz; mv NA12878_${fl}_R2.part.fq.gz NA12878_${fl}_R2.fq.gz
  rm -f split/$1.bam
  echo "$fl done $(( $(date +%s) - t0 ))s"
}
export -f one; export t0
ls split/*.bam | xargs -n1 basename | sed 's/\.bam$//' | xargs -P 3 -I{} bash -c 'one {}'
rmdir split/ 2>/dev/null || true
echo "all lanes $(( $(date +%s) - t0 ))s"
