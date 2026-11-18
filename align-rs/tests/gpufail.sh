#!/bin/bash
# GPU calls that fail (review finding 9: device memory runs out) are redone on the CPU: ALIGN_RS_GPU_FAIL_GROW=N makes
# every device allocation after the first N fail, so N = 0, 1, 2, ... fails each stage in turn (read upload,
# seeding, SA, extension, mate SW). The `small` set (hg38) must still give bwa-mem2's records.
# Usage: tests/gpufail.sh [align-rs binary] [out dir]
REPO=$(cd "$(dirname "$0")/../.." && pwd); G=${VARCALL_DATA:?set VARCALL_DATA to the test data directory}; S=$G/wgs/fastq/sub; REF=$G/ref_hg38/Homo_sapiens_assembly38.fasta
A=$(realpath ${1:-$REPO/target/release/align-rs}); O=${2:-/tmp/align-rs-gpufail}; mkdir -p $O
rg='@RG\tID:NA12878.L\tSM:NA12878\tPL:ILLUMINA\tLB:NA12878\tPU:L'
for n in 0 1 2 3 4 6 8 12 16; do
  ALIGN_RS_GPU_FAIL_GROW=$n $A mem -t 32 -K 100000000 -Y --gpu -R "$rg" $REF $S/small.R1.fq.gz $S/small.R2.fq.gz -o $O/f$n.sam 2> $O/f$n.log; rc=$?
  w=$(grep -o "GPU [a-zA-Z ]* failed" $O/f$n.log | sort -u | sed 's/GPU //; s/ failed//' | paste -sd, -)
  echo "  fail after $n: exit $rc, $(cmp -s <(grep -v '^@PG' $S/small.bwamem2.sam) <(grep -v '^@PG' $O/f$n.sam) && echo "identical to bwa-mem2" || echo DIFFERENT) (CPU fallback: ${w:-none})"
done
