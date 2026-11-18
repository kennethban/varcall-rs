#!/usr/bin/env bash
# MarkDuplicates on CRAM vs the same reads as BAM (REVIEW_FINDINGS.md: unplaced unmapped reads dropped from CRAM
# input). Input: chr20:10-11 Mb + the unplaced reads of one 30x lane (align-rs output). The BAM path is compared with
# Picard elsewhere (markdup/cmp.sh), so here the CRAM run must give the same records and metrics.  run.sh
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/markdup/cram" || exit 1; G=$VARCALL_DATA

R=$G/ref_hg38/Homo_sapiens_assembly38.fasta; X=${VARIANT_RS:-$REPO/target/release/variant-rs}
if [ ! -s in/t.cram ]; then
  mkdir -p in; L=$(ls $G/pipeline/nf/work/na12878_fastq2/*/*/NA12878.HCLHLDSXX.1.bam | head -1)
  samtools view -b -o in/t.bam $L chr20:10000000-11000000 '*' && samtools index in/t.bam
  samtools view -C -T $R -o in/t.cram in/t.bam && samtools index in/t.cram
fi
mkdir -p out
for f in bam cram; do $X markdup -I in/t.$f -R $R -O out/t.$f.md.bam -M out/t.$f.metrics 2> out/t.$f.log || echo "$f run FAILED"; done
echo "input: $(samtools view -c in/t.bam) records, unplaced $(samtools view -c in/t.bam '*')"
for f in bam cram; do echo "$f input -> output $(samtools view -c out/t.$f.md.bam) records, unplaced $(samtools view out/t.$f.md.bam | awk '$3 == "*"' | wc -l)"; done
a=$(samtools view out/t.bam.md.bam | md5sum | cut -c1-32); b=$(samtools view out/t.cram.md.bam | md5sum | cut -c1-32)
m=$(cmp -s <(grep -v '^#' out/t.bam.metrics) <(grep -v '^#' out/t.cram.metrics) && echo identical || echo DIFFERENT)
echo "CRAM vs BAM input: records $([ "$a" = "$b" ] && echo identical || echo DIFFERENT), metrics $m"
