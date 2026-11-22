#!/usr/bin/env bash
# GATK truth for S3 group 2 (HaplotypeCaller arguments from WARP's command lines), HG002 chr8:1-15M.
# Outputs: truth/<case>.{vcf.gz,g.vcf.gz}; each case's command line is in truth/<case>.cmd.
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/hg002/group2" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
REF=../ref/GRCh38_no_alt_analysis_set.fa; BAM=../bam/HG002.bam
DBSNP=$VARCALL_DATA/bqsr/resources/Homo_sapiens_assembly38.dbsnp138.vcf.gz
L="-L chr8:1-15000000"
GQB="-GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90"
mkdir -p truth
declare -A C
# WARP GermlineVariantDiscovery.wdl HaplotypeCaller_GATK4_VCF, make_gvcf true / false
C[warp_gvcf]="$L -contamination 0 -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation $GQB -ERC GVCF"
C[warp_vcf]="$L -contamination 0 -G StandardAnnotation -G StandardHCAnnotation $GQB"
C[as_vcf]="$L -G AS_StandardAnnotation"
C[dbsnp_vcf]="$L --dbsnp $DBSNP"
C[dbsnp_gvcf]="$L --dbsnp $DBSNP -ERC GVCF"
C[contam_vcf]="$L -contamination 0.05"
C[contam_gvcf]="$L -contamination 0.05 -ERC GVCF"
# contamination across contigs: the downsampling generator carries over from contig to contig
C[contam3_vcf]="-L chr8:1-3000000 -L chr9:1-3000000 -L chr10:1-3000000 -contamination 0.05"
C[contam3_gvcf]="-L chr8:1-3000000 -L chr9:1-3000000 -L chr10:1-3000000 -contamination 0.05 -ERC GVCF"
# WARP's GVCF command on a whole HaplotypeCaller shard (hg002/plan/shard023.intervals, chr8:1-74211328)
C[warp_gvcf_shard023]="-L $VARCALL_DATA/hg002/plan/shard023.intervals -contamination 0 -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation $GQB -ERC GVCF"
C[iv_vcf]="-L chr8:2000000-3000000 -L chr8:6000000-6500000 -ip 500 -XL chr8:2400000-2500000 -ixp 50"
run() {
  local c=$1 args=$2 ext=vcf.gz
  [[ $args == *"-ERC GVCF"* ]] && ext=g.vcf.gz
  [ -f truth/$c.$ext.tbi ] && return 0
  echo "$args" | $HARNESS/cmdargs.sh -w > truth/$c.cmd
  local t0=$(date +%s)
  java -Xmx3g -XX:ParallelGCThreads=2 -jar $J HaplotypeCaller -R $REF -I $BAM -O truth/$c.$ext $args \
       --native-pair-hmm-threads 2 > truth/$c.log 2>&1
  echo "$c rc=$? $(( $(date +%s) - t0 ))s" >> truth/times.txt
}
for c in "${!C[@]}"; do run $c "${C[$c]}" & done
wait
cat truth/times.txt
