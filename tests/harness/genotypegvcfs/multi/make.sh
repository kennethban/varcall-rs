#!/usr/bin/env bash
# Multi-sample GenotypeGVCFs inputs and GATK truth (S3 group 4). Three samples over chr8:1-15M with WARP's
# HaplotypeCaller GVCF command (reference: GRCh38 full analysis set, which both BAMs' dictionaries fit): HG002,
# NA12878, and HG002's PCR+ library renamed HG002pcr. Each reblocked with WARP's command; GATK CombineGVCFs over the
# raw and over the reblocked GVCFs; then GenotypeGVCFs truth on both (WARP's arguments and the defaults).
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/genotypegvcfs/multi" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
REF=$VARCALL_DATA/bqsr/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
DBSNP=$VARCALL_DATA/bqsr/resources/Homo_sapiens_assembly38.dbsnp138.vcf.gz
L="-L chr8:1-15000000"
GQB="-GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90"
WG="-contamination 0 -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation $GQB -ERC GVCF"
mkdir -p in truth
hc() { # name bam
  [ -f in/$1.g.vcf.gz.tbi ] || java -Xmx3g -XX:ParallelGCThreads=2 -jar $J HaplotypeCaller -R $REF -I $2 $L -O in/$1.g.vcf.gz $WG --native-pair-hmm-threads 2 > in/$1.hc.log 2>&1 || echo "$1 HC FAILED"
}
hc HG002 $VARCALL_DATA/hg002/bam/HG002.bam &
hc NA12878 $VARCALL_DATA/wgs/bam/NA12878.bam &
hc pcr $VARCALL_DATA/hg002_pcrplus/bam/HG002.pcrplus.bam &
wait
# rename the PCR+ sample
if [ ! -f in/HG002pcr.g.vcf.gz.tbi ]; then
  zcat in/pcr.g.vcf.gz | sed '/^#CHROM/s/\tHG002$/\tHG002pcr/' | bgzip > in/HG002pcr.g.vcf.gz && tabix -p vcf in/HG002pcr.g.vcf.gz
fi
for s in HG002 NA12878 HG002pcr; do
  [ -f in/$s.rb.g.vcf.gz.tbi ] || java -Xmx3g -jar $J ReblockGVCF -R $REF -V in/$s.g.vcf.gz -do-qual-approx --floor-blocks -GQB 20 -GQB 30 -GQB 40 -O in/$s.rb.g.vcf.gz > in/$s.rb.log 2>&1 &
done
wait
[ -f in/combined.g.vcf.gz.tbi ] || java -Xmx4g -jar $J CombineGVCFs -R $REF -V in/HG002.g.vcf.gz -V in/NA12878.g.vcf.gz -V in/HG002pcr.g.vcf.gz -G StandardAnnotation -G AS_StandardAnnotation -O in/combined.g.vcf.gz > in/combined.log 2>&1 &
[ -f in/combined_rb.g.vcf.gz.tbi ] || java -Xmx4g -jar $J CombineGVCFs -R $REF -V in/HG002.rb.g.vcf.gz -V in/NA12878.rb.g.vcf.gz -V in/HG002pcr.rb.g.vcf.gz -G StandardAnnotation -G AS_StandardAnnotation -O in/combined_rb.g.vcf.gz > in/combined_rb.log 2>&1 &
wait
W="-D $DBSNP -G StandardAnnotation -G AS_StandardAnnotation --only-output-calls-starting-in-intervals $L --merge-input-intervals"
declare -A C
C[m_warp]="-V in/combined.g.vcf.gz $W"
C[m_warp_rb]="-V in/combined_rb.g.vcf.gz $W"
C[m_default]="-V in/combined.g.vcf.gz"
C[m_default_rb]="-V in/combined_rb.g.vcf.gz"
for c in "${!C[@]}"; do
  echo "${C[$c]}" > truth/$c.cmd
  [ -f truth/$c.vcf.gz.tbi ] || java -Xmx4g -jar $J GenotypeGVCFs -R $REF ${C[$c]} -O truth/$c.vcf.gz > truth/$c.log 2>&1 &
done
wait
for c in "${!C[@]}"; do echo "$c $(zcat truth/$c.vcf.gz | grep -vc '^#')"; done
