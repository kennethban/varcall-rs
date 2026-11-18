#!/usr/bin/env bash
# align-rs vs bwa-mem2 2.2.1 on the subsampled lane sets (wgs/fastq/sub: tiny 20k pairs, small 1M, indel 3.4k
# synthetic) and synthetic edge-case sets from tests/synth.py (longdel: 2x250 with a 75-100 bp deletion; smallref:
# against tests/data/small.fa, a 5 Mb reference, with zero-length reads; longread: 700-3000 bp pairs): SAM records identical (header apart from @PG), and -o x.bam vs bwa-mem2 | samtools sort: records,
# header apart from @PG, .bai bytes vs samtools index. Then --lanes (several lanes in one run, one without reads):
# each lane's BAM and .bai byte-identical to the single-lane command in its @PG line.   cmp.sh [OUTDIR]
REPO=$(cd "$(dirname "$0")/../.." && pwd); G=${VARCALL_DATA:?set VARCALL_DATA to the test data directory}; S=$G/wgs/fastq/sub; O=${1:-/tmp/align-rs-cmp}; mkdir -p $O
command -v samtools > /dev/null || { echo "samtools: not found (FAILED)"; exit 1; }
A=$REPO/target/release/align-rs; REF=$G/ref_hg38/Homo_sapiens_assembly38.fasta
for n in tiny small indel longdel smallref longread; do
  rg='@RG\tID:NA12878.L\tSM:NA12878\tPL:ILLUMINA\tLB:NA12878\tPU:L'; case $n in indel|longdel|smallref|longread) rg='@RG\tID:x\tSM:x';; esac
  R=$REF; [ $n = smallref ] && R=$REPO/align-rs/tests/data/small.fa
  $A mem -t 32 -K 100000000 -Y --gpu -R "$rg" $R $S/$n.R1.fq.gz $S/$n.R2.fq.gz -o $O/$n.sam 2> $O/$n.sam.log
  cmp -s <(grep -v "^@PG" $S/$n.bwamem2.sam) <(grep -v "^@PG" $O/$n.sam) && echo "$n SAM: identical" || echo "$n SAM: DIFFERENT"
  [ -f $O/$n.truth.bam ] || samtools sort -@16 -l 1 -o $O/$n.truth.bam $S/$n.bwamem2.sam 2>/dev/null
  $A mem -t 32 -K 100000000 -Y --gpu -R "$rg" $R $S/$n.R1.fq.gz $S/$n.R2.fq.gz -o $O/$n.bam 2> $O/$n.bam.log
  [ -s $O/$n.bam ] && [ -s $O/$n.truth.bam ] || { echo "$n BAM: missing (FAILED)"; continue; }
  r=$( [ "$(samtools view $O/$n.bam | md5sum)" = "$(samtools view $O/$n.truth.bam | md5sum)" ] && echo identical || echo DIFFERENT)
  h=$(cmp -s <(samtools view -H --no-PG $O/$n.bam | grep -v "^@PG") <(samtools view -H --no-PG $O/$n.truth.bam | grep -v "^@PG") && echo identical || echo DIFFERENT)
  cp $O/$n.bam.bai $O/$n.own.bai; samtools index $O/$n.bam
  i=$(cmp -s $O/$n.own.bai $O/$n.bam.bai && echo identical || echo DIFFERENT)
  echo "$n BAM: records $r, header $h, bai $i"
  rm -f $O/$n.sam
done
# --lanes: the single-lane runs first (each writes the path its lane will), then all lanes in one run
D=$O/lanes; mkdir -p $D; : | gzip > $D/empty.fq.gz; L=$D/lanes.tsv; : > $L
for n in tiny empty small indel; do
  rg='@RG\tID:'$n'\tSM:NA12878'; f1=$S/$n.R1.fq.gz; f2=$S/$n.R2.fq.gz; [ $n = empty ] && f1=$D/empty.fq.gz && f2=$f1
  printf '%s\t%s\t%s\t%s\n' $D/$n.bam "$rg" $f1 $f2 >> $L
  $A mem -t 32 -K 100000000 -Y --gpu -R "$rg" -l 1 $REF $f1 $f2 -o $D/$n.bam 2> $D/$n.single.log
  mv $D/$n.bam $D/$n.single.bam; mv $D/$n.bam.bai $D/$n.single.bam.bai
done
$A mem -t 32 -K 100000000 -Y --gpu -l 1 --lanes $L $REF 2> $D/lanes.log || echo "--lanes run: FAILED (exit $?)"
for n in tiny empty small indel; do
  [ -s $D/$n.bam ] && [ -s $D/$n.single.bam ] || { echo "lanes $n: missing (FAILED)"; continue; }
  b=$(cmp -s $D/$n.bam $D/$n.single.bam && echo identical || echo DIFFERENT)
  i=$(cmp -s $D/$n.bam.bai $D/$n.single.bam.bai && echo identical || echo DIFFERENT)
  echo "lanes $n ($(samtools view -c $D/$n.bam) records): BAM $b, bai $i (vs the single-lane run)"
done
# two lanes with one output (spelled differently): refused before any work, nothing written
L2=$D/dup.tsv; rm -f $D/dup.bam*
printf '%s\t%s\t%s\t%s\n' $D/dup.bam '@RG\tID:a\tSM:x' $S/smallref.R1.fq.gz $S/smallref.R2.fq.gz > $L2
printf '%s\t%s\t%s\t%s\n' $D/../lanes/dup.bam '@RG\tID:b\tSM:x' $S/smallref.R1.fq.gz $S/smallref.R2.fq.gz >> $L2
$A mem -t 32 -K 100000000 -Y --gpu -l 1 --lanes $L2 $REPO/align-rs/tests/data/small.fa > $D/dup.log 2>&1; rc=$?
[ $rc = 1 ] && grep -q "output of more than one lane" $D/dup.log && [ ! -e $D/dup.bam ] && echo "lanes duplicate output: refused" || echo "lanes duplicate output: NOT REFUSED (exit $rc, FAILED)"
