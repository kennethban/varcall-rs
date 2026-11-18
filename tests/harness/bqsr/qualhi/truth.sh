#!/usr/bin/env bash
# GATK 4.7 truth for base qualities outside 0-93 (REVIEW_FINDINGS.md). Bytes >= 128 are negative in Java: GATK fails
# wherever it sets such qualities on a read (BaseRecalibrator when it clips the read, ApplyBQSR always), and skips
# them otherwise. Bytes 94-127: BaseRecalibrator fails (no table cell), ApplyBQSR uses the read-group estimate.
# Expected: hi -> both fail; noclip -> table OK, ApplyBQSR fails; mid -> table fails, ApplyBQSR OK (with noclip's table).
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/bqsr/qualhi" || exit 1
J="java -Xmx4g -jar ${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}"; REF=../ref/GRCh38_full_analysis_set_plus_decoy_hla.fa; R=../resources
KS="--known-sites $R/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites $R/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --known-sites $R/Homo_sapiens_assembly38.known_indels.vcf.gz"
mkdir -p truth
for m in hi noclip mid; do [ -s $m.bam ] || { python3 $HARNESS/bqsr/qualhi/make.py src.bam $m.bam $m && samtools index $m.bam; }; done
# nocigar: every 100th mapped read with CIGAR '*' (WellformedReadFilter drops it; it is not an out-of-range case)
[ -s nocigar.bam ] || { samtools view -h src.bam | awk 'BEGIN{OFS="\t"} /^@/ {print; next} {if (and($2, 4) == 0 && ++n % 100 == 0) $6 = "*"; print}' | samtools view -b -o nocigar.bam - && samtools index nocigar.bam; }
$J BaseRecalibrator -R $REF -I nocigar.bam $KS -O truth/nocigar.table > truth/nocigar.table.log 2>&1; echo "BaseRecalibrator nocigar rc=$?"
for m in hi noclip mid; do $J BaseRecalibrator -R $REF -I $m.bam $KS -O truth/$m.table > truth/$m.table.log 2>&1; echo "BaseRecalibrator $m rc=$?"; done
for m in hi noclip mid; do $J ApplyBQSR -R $REF -I $m.bam -bqsr truth/noclip.table -O truth/$m.bqsr.bam > truth/$m.apply.log 2>&1; echo "ApplyBQSR $m rc=$?"; done
