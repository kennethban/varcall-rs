#!/usr/bin/env bash
# More GenomicsDBImport truth (truth.sh and truth2.sh make the WARP-style workspaces): the rb.map / raw.map imports
# with sarek's --bypass-feature-reader --genomicsdb-shared-posixfs-optimizations, and sarek's form: one --variant per
# GVCF and a BED file of intervals (two of them overlapping: merged into one array), with and without
# --merge-input-intervals (its WES setting). GenotypeGVCFs truth on ws_sarek.
HARNESS=$(cd "$(dirname "$0")/.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/genomicsdb" || exit 1
J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
REF=$VARCALL_DATA/bqsr/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
M=$VARCALL_DATA/genotypegvcfs/multi/in
SAREK="--genomicsdb-shared-posixfs-optimizations true --bypass-feature-reader"
printf 'chr8\t100000\t2000000\nchr8\t1999990\t2500000\nchr8\t5000000\t5600000\nchr8\t9000000\t9100000\n' > sarek.bed
imp() { [ -f $1/callset.json ] || { rm -rf $1; java -Xmx8g -jar $J GenomicsDBImport --genomicsdb-workspace-path $1 $2 > $1.import.log 2>&1 || echo "$1 import FAILED"; }; }
imp ws_rb_bypass "-L chr8:1-15000000 --sample-name-map rb.map --reader-threads 5 --merge-input-intervals --batch-size 50 --consolidate $SAREK" &
imp ws_raw_bypass "-L chr8:1-15000000 --sample-name-map raw.map --reader-threads 5 --merge-input-intervals --batch-size 50 --consolidate $SAREK" &
V="--variant $M/NA12878.g.vcf.gz --variant $M/HG002.g.vcf.gz --variant $M/HG002pcr.g.vcf.gz --intervals sarek.bed --tmp-dir ."
imp ws_sarek "$V $SAREK" &
imp ws_sarek_wes "$V --merge-input-intervals $SAREK" &
wait
mkdir -p truth
echo "-V gendb://ws_sarek -L sarek.bed" | $HARNESS/cmdargs.sh -w > truth/ws_sarek.cmd
[ -f truth/ws_sarek.vcf.gz.tbi ] || java -Xmx6g -jar $J GenotypeGVCFs -R $REF -V gendb://ws_sarek -L sarek.bed -O truth/ws_sarek.vcf.gz > truth/ws_sarek.log 2>&1 || echo "ws_sarek FAILED"
