#!/usr/bin/env bash
# Picard GatherBamFiles (GATK 4.7's bundled Picard) truth. Inputs as WARP makes them: ApplyBQSR per scatter interval
# (NA12878, chr20/21/22 stretches and the unmapped reads, with the gathered table of bqsr/gather), then
# BamProcessing.wdl's GatherSortedBamFiles (index + md5) / GatherUnsortedBamFiles. lanes: samtools-written BAMs
# (hg002/group4: first records in the header's block).
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/gatherbam" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
NR=$VARCALL_DATA/bqsr/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa; NA=$VARCALL_DATA/wgs/bam/NA12878.bam
T=$VARCALL_DATA/bqsr/gather/truth/scatter3.table
G4=$VARCALL_DATA/hg002/group4/in
mkdir -p in truth
ap() { # name interval
  [ -f in/$1.bam ] || java -Xmx4g -XX:ParallelGCThreads=2 -jar $J ApplyBQSR -R $NR -I $NA -L $2 --bqsr-recal-file $T -O in/$1.bam > in/$1.log 2>&1 || echo "$1 FAILED"
}
ap a chr20:10000000-12000000 & ap b chr21:20000000-22000000 & ap c chr22:30000000-32000000 & ap u unmapped &
wait
declare -A C
C[sorted]="INPUT=in/a.bam INPUT=in/b.bam INPUT=in/c.bam INPUT=in/u.bam CREATE_INDEX=true CREATE_MD5_FILE=true"
C[unsorted]="INPUT=in/a.bam INPUT=in/b.bam INPUT=in/c.bam CREATE_INDEX=false CREATE_MD5_FILE=false"
C[single]="INPUT=in/b.bam CREATE_INDEX=true"
C[lanes]="INPUT=$G4/L1.bam INPUT=$G4/L2.bam INPUT=$G4/T0.bam"
# a header-only input in the middle (an interval group without reads), and inputs given as a list file (Picard
# unrolls any input without a BAM/SAM/CRAM extension): REVIEW_FINDINGS.md, 2026-11-17
[ -f in/e.bam ] || samtools view -H -b -o in/e.bam in/a.bam
printf '%s\n' $PWD/in/a.bam $PWD/in/b.bam > in/ab.list
C[empty]="INPUT=in/a.bam INPUT=in/e.bam INPUT=in/b.bam"
C[list]="INPUT=in/ab.list"
for c in "${!C[@]}"; do
  [ -f truth/$c.bam ] && continue
  echo "${C[$c]}" | $HARNESS/cmdargs.sh -w > truth/$c.cmd
  java -Dsamjdk.compression_level=2 -Xmx2g -jar $J GatherBamFiles ${C[$c]} OUTPUT=truth/$c.bam > truth/$c.log 2>&1; echo "$c rc=$?"
done
