#!/bin/bash
# Times align-rs (env passed through) on a FASTQ pair and checks its records against a digest.
#   timerun.sh OUTDIR NAME DIGEST_NAME [R1 R2]   (default: wgs/fastq/sub/mid)
# Prints the steady batch rate, peak anonymous / file RSS, minor faults, and whether the records match
# digests.txt's entry (md5 of `samtools view` of bwa-mem2's records, coordinate-sorted).
T=$(cd "$(dirname "$0")" && pwd); REPO=$(cd $T/../.. && pwd); G=${VARCALL_DATA:?set VARCALL_DATA to the test data directory}; O=$1; NAME=$2; DG=$3
R1=${4:-$G/wgs/fastq/sub/mid.R1.fq.gz}; R2=${5:-$G/wgs/fastq/sub/mid.R2.fq.gz}
mkdir -p $O; cd $O
A=$REPO/target/release/align-rs; R=$G/ref_hg38/Homo_sapiens_assembly38.fasta
/usr/bin/time -f "%R %U %S" -o $NAME.time $A mem -t 32 -K 100000000 -Y --gpu -R '@RG\tID:NA12878.L\tSM:NA12878\tPL:ILLUMINA\tLB:NA12878\tPU:L' -l 1 $R $R1 $R2 -o $NAME.bam 2> $NAME.log &
sleep 0.5; p=$(pgrep -n -x align-rs); ma=0; mf=0
while kill -0 $p 2>/dev/null; do
  a=$(awk '/RssAnon/{print $2}' /proc/$p/status 2>/dev/null); f=$(awk '/RssFile/{print $2}' /proc/$p/status 2>/dev/null)
  [ -n "$a" ] && [ $a -gt $ma ] && ma=$a; [ -n "$f" ] && [ $f -gt $mf ] && mf=$f; sleep 0.5
done; wait
read flt usr sys < $NAME.time
got=$(samtools view $NAME.bam | md5sum | cut -c1-32); want=$(awk -v n=$DG '$1 == n {print $2}' $T/digests.txt)
[ "$got" = "$want" ] && ok="records identical to bwa-mem2" || ok="RECORDS DIFFER ($got vs $want)"
echo "$NAME: $($T/steady.sh $NAME.log); user ${usr} s, sys ${sys} s; peak anon $((ma/1024)) MB, file $((mf/1024)) MB, faults $flt; $ok"
