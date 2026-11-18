#!/usr/bin/env bash
# Multi-sample HaplotypeCaller test inputs and GATK truth (S3 group 4). chr8 extracts of HG002, NA12878 and HG002's
# PCR+ library (renamed HG002pcr), all with the reference's 3366-contig dictionary; an 11-sample set (each source split
# by the last digit of the read name, so mates stay together: HG002 and NA12878 in 4, HG002pcr in 3); one BAM holding
# two samples. Truth: GATK HaplotypeCaller in VCF mode (multi-sample calling is VCF-only in GATK).
HARNESS=$(cd "$(dirname "$0")/../../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/hg002/group4/ms" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
REF=$VARCALL_DATA/bqsr/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
DBSNP=$VARCALL_DATA/bqsr/resources/Homo_sapiens_assembly38.dbsnp138.vcf.gz
SQ=in/sq.txt
mkdir -p in truth
grep "^@SQ" ${REF%.fa}.dict | cut -f1-3 > $SQ
R=chr8:1-15000000
# extract <out> <bam> <sample>: chr8 reads, the common dictionary, read groups renamed to <sample>
extract() {
  [ -f in/$1.bam.bai ] && return
  samtools view -h $2 $R 2>/dev/null | awk -v sm=$3 -v sq=$SQ 'BEGIN{OFS="\t"}
    /^@SQ/ { next }
    /^@HD/ { print; while ((getline l < sq) > 0) print l; next }
    /^@RG/ { sub(/SM:[^\t]*/, "SM:" sm); print; next }
    /^@/ { print; next }
    { print }' | samtools view -b -o in/$1.bam - && samtools index in/$1.bam
}
extract HG002 $VARCALL_DATA/hg002/bam/HG002.bam HG002 &
extract NA12878 $VARCALL_DATA/wgs/bam/NA12878.bam NA12878 &
extract HG002pcr $VARCALL_DATA/hg002_pcrplus/bam/HG002.pcrplus.bam HG002pcr &
wait
# split <src> <n> <k> <sample>: reads whose name's last digit % n == k, as sample <sample> (read groups suffixed)
split() {
  [ -f in/$4.bam.bai ] && return
  samtools view -h in/$1.bam chr8:1-5000000 | awk -v n=$2 -v k=$3 -v sm=$4 'BEGIN{OFS="\t"}
    /^@RG/ { match($0, /ID:[^\t]*/); id = substr($0, RSTART + 3, RLENGTH - 3); sub(/ID:[^\t]*/, "ID:" id "." sm); sub(/SM:[^\t]*/, "SM:" sm); print; next }
    /^@/ { print; next }
    { if (substr($1, length($1), 1) % n != k) next
      for (i = 12; i <= NF; i++) if ($i ~ /^RG:Z:/) $i = $i "." sm
      print }' | samtools view -b -o in/$4.bam - && samtools index in/$4.bam
}
for k in 0 1 2 3; do split HG002 4 $k HG002_$k & split NA12878 4 $k NA12878_$k & done
for k in 0 1 2; do split HG002pcr 3 $k HG002pcr_$k & done
wait
[ -f in/two.bam.bai ] || { samtools merge -f -o in/two.bam in/HG002.bam in/NA12878.bam && samtools index in/two.bam; }
I3="-I in/HG002.bam -I in/NA12878.bam -I in/HG002pcr.bam"
I11=$(for s in HG002_0 HG002_1 HG002_2 HG002_3 NA12878_0 NA12878_1 NA12878_2 NA12878_3 HG002pcr_0 HG002pcr_1 HG002pcr_2; do echo -n "-I in/$s.bam "; done)
declare -A C
C[ms2]="-I in/HG002.bam -I in/NA12878.bam -L chr8:1-15000000"
C[ms3]="$I3 -L chr8:1-15000000"
C[ms11]="$I11 -L chr8:1-5000000"
C[onebam]="-I in/two.bam -L chr8:5000000-10000000"
C[ms3_opts]="$I3 -L chr8:1-5000000 --dbsnp $DBSNP -contamination 0.05"
run() {
  echo "${C[$1]}" > truth/$1.cmd
  [ -f truth/$1.vcf.gz.tbi ] || java -Xmx5g -XX:ParallelGCThreads=2 -jar $J HaplotypeCaller -R $REF ${C[$1]} -O truth/$1.vcf.gz --native-pair-hmm-threads 2 > truth/$1.log 2>&1
}
run ms2 & run ms3 & run ms11 & wait
run onebam & run ms3_opts & wait
for c in "${!C[@]}"; do echo "$c $(grep -c -v '^#' <(zcat truth/$c.vcf.gz)) $(grep 'Elapsed time' truth/$c.log)"; done
