#!/usr/bin/env bash
# GATK 4.7 panel-of-normals truth, as GATK's mutect2_pon.wdl makes one: Mutect2 tumor-only on the four DREAM normals
# (b37 chr20, --max-mnp-distance 0), GenomicsDBImport -R -V ... -L (one interval: one array; an interval list of three:
# three arrays), CreateSomaticPanelOfNormals -R --germline-resource -V gendb:// (defaults; --min-sample-count 1; no
# germline resource; the three-array workspace: intervals with gaps, as abutting ones merge), and SelectVariants on the workspaces (the merged records).
HARNESS=$(cd "$(dirname "$0")/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
mkdir -p "$VARCALL_DATA/somatic/pon" && cd "$VARCALL_DATA/somatic/pon" || exit 1
J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
D=../dream/data; REF=$D/human_g1k_v37.20.21.fasta; GN=$D/very-small-gnomad.vcf
mkdir -p in truth
for i in 1 2 3 4; do
  [ -f in/normal_$i.vcf.gz.tbi ] && continue
  java -Xmx8g -jar $J Mutect2 -R $REF -I $D/normal_$i.bam -L 20 --max-mnp-distance 0 -O in/normal_$i.vcf.gz > in/normal_$i.log 2>&1 || echo "normal_$i FAILED"
done
printf '20:1-20000000\n20:20000101-45000000\n20:45000101-63025520\n' > in/three.intervals
V="-V in/normal_1.vcf.gz -V in/normal_2.vcf.gz -V in/normal_3.vcf.gz -V in/normal_4.vcf.gz"
imp() { [ -f $1/callset.json ] || { rm -rf $1; java -Xmx8g -jar $J GenomicsDBImport --genomicsdb-workspace-path $1 -R $REF $V -L $2 > $1.import.log 2>&1 || echo "$1 import FAILED"; }; }
imp ws_pon 20
imp ws_pon3 in/three.intervals
# what a VariantWalker receives from each workspace (SelectVariants reads gendb:// with the same options as
# CreateSomaticPanelOfNormals): the target for the GenomicsDB merge on Mutect2 VCFs
for w in ws_pon ws_pon3; do
  [ -f truth/$w.received.vcf ] || java -Xmx8g -jar $J SelectVariants -R $REF -V gendb://$w -O truth/$w.received.vcf > truth/$w.received.log 2>&1 || echo "$w SelectVariants FAILED"
done
pon() { # name args
  echo "$2" > truth/$1.cmd
  [ -f truth/$1.vcf ] || java -Xmx8g -jar $J CreateSomaticPanelOfNormals -R $REF $2 -O truth/$1.vcf > truth/$1.log 2>&1 || echo "$1 FAILED"
  echo "$1 $(grep -vc '^#' truth/$1.vcf) records"
}
pon default "--germline-resource $GN -V gendb://ws_pon"
pon min1 "--germline-resource $GN -V gendb://ws_pon --min-sample-count 1"
pon nogermline "-V gendb://ws_pon"
pon three "--germline-resource $GN -V gendb://ws_pon3"
pon strict "--germline-resource $GN -V gendb://ws_pon --min-sample-count 3 --max-germline-probability 0.05"
