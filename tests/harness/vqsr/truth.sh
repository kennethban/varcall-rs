#!/usr/bin/env bash
# GATK 4.7 VariantRecalibrator / ApplyVQSR truth (inputs: download.sh). Each case's arguments go to truth/<case>.cmd
# with @O@ for the output folder; truth/cases.txt lists the cases in run order (ApplyVQSR steps read earlier outputs).
#   s_*: single sample (NA12878 30x, 5.06M records)   c_*: cohort (1000 Genomes, 50 of 3202 samples, chr20:5-17 Mb)
#   warp_*: WARP JointGenotyping's settings (the WGS one downsamples 3.1M training sites to 2.5M)
#   sarek_*: nf-core/sarek joint_germline's (default Gaussians and tranches, resource genotypes checked, -R)
#   inb_*: GATK best practice for cohorts, with InbreedingCoeff (absent at 18% of the sites)
#   both: -mode BOTH;  *apply*: ApplyVQSR, chained SNP -> INDEL or INDEL -> SNP as the pipelines run it
HARNESS=$(cd "$(dirname "$0")/.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/vqsr" || exit 1
J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}

REF=../wgs/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
S=in/NA12878.vcf.gz C=in/cohort.vcf.gz
HAPMAP=resources/hapmap_3.3.hg38.vcf.gz OMNI=resources/1000G_omni2.5.hg38.vcf.gz
G1K=resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz AXIOM=resources/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
MILLS=../bqsr/resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz DBSNP=../bqsr/resources/Homo_sapiens_assembly38.dbsnp138.vcf.gz
snpres() { echo "--resource:hapmap,known=false,training=true,truth=true,prior=15$1 $HAPMAP --resource:omni,known=false,training=true,truth=true,prior=12$1 $OMNI --resource:1000G,known=false,training=true,truth=false,prior=10$1 $G1K --resource:dbsnp,known=true,training=false,truth=false,prior=$2 $DBSNP"; }
WARP_SNP_T="-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0"
WARP_INDEL_T="-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0"
WARP_INDELRES="--resource:mills,known=false,training=true,truth=true,prior=12 $MILLS --resource:axiomPoly,known=false,training=true,truth=false,prior=10 $AXIOM --resource:dbsnp,known=true,training=false,truth=false,prior=2 $DBSNP"
vr() { echo "VariantRecalibrator -V $2 -O @O@/$1.recal.vcf.gz --tranches-file @O@/$1.tranches"; }
warp_snp() { echo "$(vr $1 $2) --trust-all-polymorphic $WARP_SNP_T -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP -mode SNP --max-gaussians 6 $(snpres '' 7)"; }
warp_indel() { echo "$(vr $1 $2) --trust-all-polymorphic $WARP_INDEL_T -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP -mode INDEL --max-gaussians 4 $WARP_INDELRES"; }
apply() { echo "ApplyVQSR $1 --recal-file @O@/$2.recal.vcf.gz --tranches-file @O@/$2.tranches --truth-sensitivity-filter-level $3 -mode $4 -O @O@/$5.vcf.gz"; }

CASES=(); declare -A CMD
add() { CASES+=($1); CMD[$1]=$2; }
add s_warp_snp "$(warp_snp s_warp_snp $S)"
add s_warp_indel "$(warp_indel s_warp_indel $S)"
add s_sarek_snp "$(vr s_sarek_snp $S) -R $REF -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP $(snpres .0 2.0)"
add s_sarek_indel "$(vr s_sarek_indel $S) -R $REF -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -mode INDEL --resource:mills,known=false,training=true,truth=true,prior=12.0 $MILLS --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP"
add s_sarek_apply_snp "$(apply "-R $REF -V $S" s_sarek_snp 99.9 SNP s_sarek_apply_snp)"
add s_sarek_apply_indel "$(apply "-R $REF -V @O@/s_sarek_apply_snp.vcf.gz" s_sarek_indel 99.9 INDEL s_sarek_apply_indel)"
add s_warp_apply_indel "$(apply "-V $S" s_warp_indel 99.0 INDEL s_warp_apply_indel) --create-output-variant-index true"
add s_warp_apply_snp "$(apply "-V @O@/s_warp_apply_indel.vcf.gz" s_warp_snp 99.7 SNP s_warp_apply_snp) --create-output-variant-index true"
add c_warp_snp "$(warp_snp c_warp_snp $C)"
add c_warp_indel "$(warp_indel c_warp_indel $C)"
add c_inb_snp "$(vr c_inb_snp $C) -R $REF -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff -mode SNP $(snpres '' 7)"
add c_inb_indel "$(vr c_inb_indel $C) -R $REF -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff -mode INDEL --max-gaussians 4 $WARP_INDELRES"
add c_both "$(vr c_both $C) -an QD -an FS -an SOR -an ReadPosRankSum -mode BOTH --max-gaussians 4 --resource:hapmap,known=false,training=true,truth=true,prior=15 $HAPMAP --resource:mills,known=false,training=true,truth=true,prior=12 $MILLS"
add c_warp_apply_indel "$(apply "-V $C" c_warp_indel 99.0 INDEL c_warp_apply_indel) --create-output-variant-index true"
add c_warp_apply_snp "$(apply "-V @O@/c_warp_apply_indel.vcf.gz" c_warp_snp 99.7 SNP c_warp_apply_snp) --create-output-variant-index true"
add c_inb_apply_snp "$(apply "-R $REF -V $C" c_inb_snp 99.5 SNP c_inb_apply_snp)"
add c_inb_apply_indel "$(apply "-R $REF -V @O@/c_inb_apply_snp.vcf.gz" c_inb_indel 99.0 INDEL c_inb_apply_indel)"
add c_both_apply "$(apply "-V $C" c_both 99.0 BOTH c_both_apply)"

mkdir -p truth; printf '%s\n' "${CASES[@]}" > truth/cases.txt
for c in "${CASES[@]}"; do
  echo "${CMD[$c]}" | $HARNESS/cmdargs.sh -w > truth/$c.cmd
  [ -f truth/$c.log ] && grep -q 'rc=0' truth/$c.log && continue
  t=$(date +%s)
  java -Xmx16g -jar $J ${CMD[$c]//@O@/truth} > truth/$c.log 2>&1; rc=$?
  echo "rc=$rc $(( $(date +%s) - t ))s" | tee -a truth/$c.log | sed "s/^/$c /"
done
