#!/usr/bin/env bash
# HaplotypeCaller without -L on a BAM whose @SQ order differs from the reference .dict (chrM, chr21, chr20 first).
# GATK traverses in the reads' dictionary order (AssemblyRegionWalker.onStartup); review finding 16.
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/hcorder" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
REF=../hg002/ref/GRCh38_no_alt_analysis_set.fa; SRC=../hg002/bam/HG002.bam
if [ ! -s reordered.bam ]; then
  # chrM subsampled to ~1000x (still downsampled at 50 reads per start; GATK's GVCF mode runs out of heap on all of it)
  samtools merge -f -o src.bam <(samtools view -b $SRC chr20:10000000-10300000 chr21:20000000-20300000) <(samtools view -b -s 0.05 $SRC chrM)
  samtools view -H src.bam > h.sam
  { grep -v '^@SQ' h.sam | grep '^@HD'
    for c in chrM chr21 chr20; do grep -P "^@SQ\tSN:$c\t" h.sam; done
    grep '^@SQ' h.sam | grep -vP "^@SQ\tSN:(chrM|chr21|chr20)\t"
    grep -v '^@SQ' h.sam | grep -v '^@HD'; } > h2.sam
  # reheader keeps the reference ids, so rebuild from SAM text (names) under the new header
  { cat h2.sam; samtools view src.bam; } | samtools sort -o reordered.bam - && samtools index reordered.bam
  rm -f src.bam h.sam h2.sam
fi
mkdir -p truth
for m in vcf gvcf; do
  e=""; [ $m = gvcf ] && e="-ERC GVCF"
  java -Xmx16g -jar $J HaplotypeCaller -R $REF -I reordered.bam -O truth/out.$m.gz $e > truth/$m.log 2>&1; echo "GATK $m rc=$?"
done
