#!/usr/bin/env bash
# GATK 4.7 truth for WARP's BQSR arguments (tasks/wdl/BamProcessing.wdl): --use-original-qualities in
# BaseRecalibrator and ApplyBQSR, --static-quantized-quals 10 20 30, --create-output-bam-md5,
# --add-output-sam-program-record; plus --round-down-quantized with 10..50. Input: NA12878 chr20:10-20M with OQ tags
# (GATK ApplyBQSR --emit-original-quals), the OQ tag then removed from reads whose name ends in an even digit.
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/bqsr" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
G="java -Xmx8g -XX:ParallelGCThreads=2 -jar $J"
REF=ref/GRCh38_full_analysis_set_plus_decoy_hla.fa; R=resources; W=truth/warp; L=chr20:10000000-20000000
KS="--known-sites $R/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites $R/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --known-sites $R/Homo_sapiens_assembly38.known_indels.vcf.gz"
mkdir -p $W
t() { local n=$1; shift; local t0=$(date +%s); "$@" > $W/$n.log 2>&1; echo "$n rc=$? $(( $(date +%s) - t0 ))s"; }
t emit $G ApplyBQSR -R $REF -I $VARCALL_DATA/wgs/bam/NA12878.bam -L $L --bqsr-recal-file truth/na12878.chr20.recal.table \
  --emit-original-quals -O $W/na.oq_all.bam
samtools view -h $W/na.oq_all.bam | awk 'BEGIN{OFS="\t"} /^@/{print; next}
  $1 ~ /[02468]$/ {o=$1; for (i=2; i<=NF; i++) if ($i !~ /^OQ:Z:/) o=o OFS $i; print o; next} {print}' \
  | samtools view -b -o $W/na.oq.bam - && samtools index $W/na.oq.bam
echo "input: $(samtools view -c $W/na.oq.bam) reads, $(samtools view $W/na.oq.bam | grep -c 'OQ:Z:') with OQ"
t table $G BaseRecalibrator -R $REF -I $W/na.oq.bam --use-original-qualities -O $W/oq.recal.table $KS -L $L
t apply $G ApplyBQSR --create-output-bam-md5 --add-output-sam-program-record -R $REF -I $W/na.oq.bam --use-original-qualities \
  -O $W/oq.bqsr.bam -bqsr $W/oq.recal.table --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 -L $L
t rounddown $G ApplyBQSR -R $REF -I $W/na.oq.bam -O $W/rd.bqsr.bam -bqsr $W/oq.recal.table --round-down-quantized \
  --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 --static-quantized-quals 40 --static-quantized-quals 50 -L $L
