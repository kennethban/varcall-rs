#!/usr/bin/env bash
# Larger GenomicsDB truth: the whole-genome reblocked HG002 GVCF (reblock/wgs shards merged with MergeVcfs), imported
# with WARP's import command. ws_wg: -L chr20 -L chr21 -L chr22 (three arrays); ws_mid: -L chr21:20000000-30000000 (an
# array starting mid-contig: records that begin before it are stored before its lower bound). GenotypeGVCFs truth:
# ws_wg over the three whole chromosomes (WARP arguments), ws_mid over chr21:22500000-27000000 (windows starting
# mid-array) and over the whole array.
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/genomicsdb" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
REF=$VARCALL_DATA/hg002/ref/GRCh38_no_alt_analysis_set.fa
DBSNP=$VARCALL_DATA/bqsr/resources/Homo_sapiens_assembly38.dbsnp138.vcf.gz
G=${VARIANT_RS:-$REPO/target/release/variant-rs}
mkdir -p in2 truth
[ -f in2/HG002.rb.g.vcf.gz.tbi ] || $G MergeVcfs $(ls $VARCALL_DATA/reblock/wgs/shard*.rb.g.vcf.gz | sed 's/^/INPUT=/') OUTPUT=in2/HG002.rb.g.vcf.gz > in2/merge.log 2>&1
printf "HG002\t$PWD/in2/HG002.rb.g.vcf.gz\n" > wg.map
imp() { [ -f $1/callset.json ] || { rm -rf $1; java -Xmx8g -jar $J GenomicsDBImport --genomicsdb-workspace-path $1 $2 --sample-name-map wg.map --batch-size 50 --reader-threads 5 --merge-input-intervals --consolidate > $1.import.log 2>&1 || echo "$1 import FAILED"; }; }
imp ws_wg "-L chr20 -L chr21 -L chr22" &
imp ws_mid "-L chr21:20000000-30000000" &
wait
W="-D $DBSNP -G StandardAnnotation -G AS_StandardAnnotation --only-output-calls-starting-in-intervals --merge-input-intervals"
gg() { echo "$2" | $HARNESS/cmdargs.sh -w > truth/$1.cmd; [ -f truth/$1.vcf.gz.tbi ] || java -Xmx8g -jar $J GenotypeGVCFs -R $REF $2 -O truth/$1.vcf.gz > truth/$1.log 2>&1 || echo "$1 FAILED"; }
gg ws_wg "$W -V gendb://ws_wg -L chr20 -L chr21 -L chr22" &
gg ws_mid "$W -V gendb://ws_mid -L chr21:22500000-27000000" &
gg ws_mid.whole "-V gendb://ws_mid -L chr21:20000000-30000000" &
wait
for c in ws_wg ws_mid ws_mid.whole; do echo "$c $(zcat truth/$c.vcf.gz | grep -vc '^#')"; done
