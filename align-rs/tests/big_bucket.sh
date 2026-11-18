#!/usr/bin/env bash
# BAM writer with one sort bucket over 4 GiB (REVIEW_FINDINGS.md: 32-bit in-bucket offsets): ~3.6M random 2x250
# pairs with long names, all unmapped, so all in the unplaced bucket. Truth: align-rs's SAM output (the path
# compared with bwa-mem2 elsewhere) through samtools sort.   big_bucket.sh [OUTDIR]  (~25 GB of disk, a few min)
REPO=$(cd "$(dirname "$0")/../.." && pwd); G=${VARCALL_DATA:?set VARCALL_DATA to the test data directory}; O=${1:-/tmp/align-rs-big}; mkdir -p $O
A=$REPO/target/release/align-rs; REF=$G/ref_hg38/Homo_sapiens_assembly38.fasta; rg='@RG\tID:x\tSM:x'
[ -s $O/u.R2.fq.gz ] || python3 $REPO/align-rs/tests/synth.py unmapped $O/u 3600000
$A mem -t 32 -K 100000000 -Y --gpu -R "$rg" -l 1 $REF $O/u.R1.fq.gz $O/u.R2.fq.gz -o $O/u.bam 2> $O/u.bam.log; rc=$?
[ $rc = 0 ] || { echo "big bucket: align-rs -o BAM exit $rc (FAILED)"; tail -3 $O/u.bam.log; exit 1; }
$A mem -t 32 -K 100000000 -Y --gpu -R "$rg" $REF $O/u.R1.fq.gz $O/u.R2.fq.gz 2> $O/u.sam.log | samtools sort -@ 16 -m 1G -l 1 -T $O/srt -o $O/truth.bam - 2>/dev/null
raw=$(samtools view -u --no-PG $O/truth.bam | wc -c)
a=$(samtools view $O/u.bam | md5sum | cut -c1-32); t=$(samtools view $O/truth.bam | md5sum | cut -c1-32)
echo "big bucket: $(samtools view -c $O/truth.bam) records, $((raw >> 20)) MiB uncompressed: $([ "$a" = "$t" ] && echo "records identical" || echo "DIFFERENT")"
samtools quickcheck $O/u.bam && echo "big bucket: BAM passes quickcheck" || echo "big bucket: quickcheck FAILED"
