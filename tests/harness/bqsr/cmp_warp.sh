#!/usr/bin/env bash
# variant-rs with WARP's BQSR command lines (truth_warp.sh) vs GATK: recalibration table with --use-original-qualities;
# ApplyBQSR records with --use-original-qualities --static-quantized-quals 10 20 30 (+ .md5 file), and with
# --round-down-quantized 10..50.
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/bqsr" || exit 1

G=${VARIANT_RS:-$REPO/target/release/variant-rs}
REF=ref/GRCh38_full_analysis_set_plus_decoy_hla.fa; R=resources; W=truth/warp; O=${1:-rust/warp}; L=chr20:10000000-20000000; mkdir -p $O
KS="--known-sites $R/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites $R/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --known-sites $R/Homo_sapiens_assembly38.known_indels.vcf.gz"
$G BaseRecalibrator -R $REF -I $W/na.oq.bam --use-original-qualities -O $O/oq.recal.table $KS -L $L 2> $O/table.log
echo "  BaseRecalibrator --use-original-qualities: table $(cmp -s $O/oq.recal.table $W/oq.recal.table && echo identical || echo DIFFERENT)"
$G ApplyBQSR --create-output-bam-md5 --add-output-sam-program-record -R $REF -I $W/na.oq.bam --use-original-qualities \
  -O $O/oq.bqsr.bam -bqsr $W/oq.recal.table --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 -L $L 2> $O/apply.log
echo "  ApplyBQSR (WARP): records $(python3 $HARNESS/bqsr/cmpbam.py $W/oq.bqsr.bam $O/oq.bqsr.bam | tail -1); md5 file $([ "$(cat $O/oq.bqsr.bam.md5)" = "$(md5sum < $O/oq.bqsr.bam | cut -c1-32)" ] && echo correct || echo WRONG)"
$G ApplyBQSR -R $REF -I $W/na.oq.bam -O $O/rd.bqsr.bam -bqsr $W/oq.recal.table --round-down-quantized \
  --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 --static-quantized-quals 40 --static-quantized-quals 50 -L $L 2> $O/rd.log
echo "  ApplyBQSR --round-down-quantized: records $(python3 $HARNESS/bqsr/cmpbam.py $W/rd.bqsr.bam $O/rd.bqsr.bam | tail -1)"
