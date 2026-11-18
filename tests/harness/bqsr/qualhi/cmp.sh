#!/usr/bin/env bash
# variant-rs vs GATK on base qualities outside 0-93 (truth.sh): the same runs succeed or fail, tables and records identical.
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/bqsr/qualhi" || exit 1
X=${VARIANT_RS:-$REPO/target/release/variant-rs}; REF=../ref/GRCh38_full_analysis_set_plus_decoy_hla.fa; R=../resources
KS="--known-sites $R/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites $R/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --known-sites $R/Homo_sapiens_assembly38.known_indels.vcf.gz"
O=${1:-/tmp/qualhi}; rm -rf $O; mkdir -p $O
st() { [ $1 = 0 ] && echo ok || echo fails; }
$X BaseRecalibrator -R $REF -I nocigar.bam $KS -O $O/nocigar.table > $O/nocigar.table.log 2>&1
echo "BaseRecalibrator nocigar (mapped reads with CIGAR '*'): table $(cmp -s truth/nocigar.table $O/nocigar.table && echo identical || echo DIFFERENT)"
for m in hi noclip mid; do
  $X BaseRecalibrator -R $REF -I $m.bam $KS -O $O/$m.table > $O/$m.table.log 2>&1; r=$?
  t=$( [ -s truth/$m.table ] && echo 0 || echo 1 )
  c=""; [ $r = 0 ] && [ $t = 0 ] && c=", table $(cmp -s truth/$m.table $O/$m.table && echo identical || echo DIFFERENT)"
  echo "BaseRecalibrator $m: GATK $(st $t), variant-rs $(st $r)$c"
  $X ApplyBQSR -R $REF -I $m.bam -bqsr truth/noclip.table -O $O/$m.bqsr.bam > $O/$m.apply.log 2>&1; r=$?
  t=$( [ -s truth/$m.bqsr.bam ] && grep -q "^.*rc=0\|SUCCESS\|Shutting down engine" truth/$m.apply.log && ! grep -q "ERROR\|Exception" truth/$m.apply.log && echo 0 || echo 1 )
  c=""; [ $r = 0 ] && [ $t = 0 ] && c=", records $(python3 $HARNESS/bqsr/cmpbam.py truth/$m.bqsr.bam $O/$m.bqsr.bam | tail -1)"
  echo "ApplyBQSR $m: GATK $(st $t), variant-rs $(st $r)$c"
done
