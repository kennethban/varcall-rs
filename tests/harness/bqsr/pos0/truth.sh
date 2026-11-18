#!/usr/bin/env bash
# ApplyBQSR without -L on a BAM with an unmapped read that keeps its contig at POS 0 (markdup/edge/pos0.py; review
# finding 32): GATK copies it through. Table: GATK BaseRecalibrator on the same BAM.
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/bqsr/pos0" || exit 1

J="java -Xmx4g -jar ${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}"; REF=$VARCALL_DATA/ref_hg38/Homo_sapiens_assembly38.fasta; R=../resources
B=../../markdup/edge/in/pos0.bam; mkdir -p truth
$J BaseRecalibrator -R $REF -I $B --known-sites $R/Homo_sapiens_assembly38.dbsnp138.vcf.gz -O truth/pos0.table > truth/bqsr.log 2>&1; echo "BaseRecalibrator rc=$?"
$J ApplyBQSR -R $REF -I $B --bqsr-recal-file truth/pos0.table -O truth/pos0.bam > truth/apply.log 2>&1; echo "ApplyBQSR rc=$?"
