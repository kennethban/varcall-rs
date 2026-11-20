#!/usr/bin/env bash
# GATK 4.7 VariantFiltration / Picard MakeSitesOnlyVcf truth. Inputs from other areas: the NA12878 30x and 1000
# Genomes cohort VCFs (vqsr/in, vqsr/download.sh), an ApplyVQSR output (vqsr/truth), GenotypeGVCFs' 11-sample joint
# calls (genotypegvcfs/multi/truth), the SEQC2 exome targets (somatic/seqc2/data). Each case's arguments go to
# truth/<case>.cmd with @O@ for the output folder; truth/cases.txt lists them in run order. Cases where GATK fails
# (an expression it cannot evaluate) are kept: variant-rs must fail too.
HARNESS=$(cd "$(dirname "$0")/.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
mkdir -p "$VARCALL_DATA/vfilter" && cd "$VARCALL_DATA/vfilter" || exit 1
J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}

S=../vqsr/in/NA12878.vcf.gz C=../vqsr/in/cohort.vcf.gz F=../vqsr/truth/c_inb_apply_indel.vcf.gz
JG=../genotypegvcfs/multi/truth/m11_warp.vcf.gz
mkdir -p in truth
# masks: the exome targets as BED and as a Picard interval list (WARP's targets_interval_list)
[ -s in/targets.bed ] || sort -k1,1 -k2,2n ../somatic/seqc2/data/targets.bed > in/targets.bed
if [ ! -s in/targets.interval_list ]; then
  (grep '^@' ../wgs/ref/GRCh38_full_analysis_set_plus_decoy_hla.dict; awk -v OFS='\t' '{print $1, $2 + 1, $3, "+", "t" NR}' in/targets.bed) > in/targets.interval_list
fi
# GATK queries a mask by interval: it needs an index (in/noidx.bed: a copy without one, which GATK refuses)
for m in targets.bed targets.interval_list; do [ -s in/$m.idx ] || java -jar $J IndexFeatureFile -I in/$m > /dev/null 2>&1; done
[ -s in/noidx.bed ] || head -1000 in/targets.bed > in/noidx.bed
# edge: the cohort's first records plus one with repeated and empty INFO
# keys, a declared flag given as DS=0 (dropped) and odd genotype formatting; edge_undeclared.vcf.gz has an INFO key
# the header does not declare (DB), which htsjdk's writer refuses
if [ ! -s in/edge.vcf.gz ]; then
  zcat $C | head -3500 > in/edge.vcf
  grep -v '^#' in/edge.vcf | sed -n 3p | awk -F'\t' -v OFS='\t' '{ $2 = $2 + 1; $8 = $8 ";QD=;FS=2;DS=0"; $10 = "0/1:0.50000:3,4:7:10:1,0,5"; print }' > in/extra.txt
  (grep '^#' in/edge.vcf; grep -v '^#' in/edge.vcf | head -3; cat in/extra.txt; grep -v '^#' in/edge.vcf | sed -n '4,$p') | bgzip > in/edge.vcf.gz
  tabix -p vcf in/edge.vcf.gz && rm in/edge.vcf in/extra.txt
fi
if [ ! -s in/edge_undeclared.vcf.gz ]; then
  zcat in/edge.vcf.gz | awk -F'\t' -v OFS='\t' '!/^#/ && ++n == 5 { $8 = $8 ";DB" } { print }' | bgzip > in/edge_undeclared.vcf.gz
  tabix -p vcf in/edge_undeclared.vcf.gz
fi

vf() { echo "VariantFiltration -V $1 -O @O@/$2"; }
CASES=(); declare -A CMD
add() { CASES+=($1); CMD[$1]=$2; }
# WARP JointGenotyping (HardFilterAndMakeSitesOnlyVcf), WGS and exome
add warp_jg "$(vf $JG warp_jg.vcf.gz) --filter-expression \"ExcessHet > 54.69\" --filter-name ExcessHet"
add warp_jg_sites "MakeSitesOnlyVcf -I @O@/warp_jg.vcf.gz -O @O@/warp_jg_sites.vcf.gz"
add warp_jg_exome "$(vf $JG warp_jg_exome.vcf.gz) --filter-expression \"ExcessHet > 10.0\" --filter-name ExcessHet --filter-not-in-mask --mask-description \"Does not overlap the targets interval list.\" --mask-name OUTSIDE_OF_TARGETS --mask in/targets.interval_list"
add warp_jg_exome_sites "MakeSitesOnlyVcf -I @O@/warp_jg_exome.vcf.gz -O @O@/warp_jg_exome_sites.vcf.gz"
# WARP single sample: DRAGEN mode's hard filter, and the HardFilterVcf task's expression (with and without
# missing values failing: hom-var sites have no MQRankSum)
add dragen "$(vf $S dragen.vcf.gz) --filter-expression \"QUAL < 10.4139\" --filter-name DRAGENHardQUAL"
add warp_hard "$(vf $S warp_hard.vcf.gz) --filter-expression \"QD < 2.0 || FS > 30.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -3.0 || ReadPosRankSum < -3.0\" --filter-name HardFiltered"
add warp_hard_missing "$(vf $S warp_hard_missing.vcf) --filter-expression \"QD < 2.0 || FS > 30.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -3.0 || ReadPosRankSum < -3.0\" --filter-name HardFiltered --missing-values-evaluate-as-failing"
# GATK's hard-filtering example: one name per expression
add bp_snp "$(vf $S bp_snp.vcf.gz) -filter \"QD < 2.0\" --filter-name QD2 -filter \"QUAL < 30.0\" --filter-name QUAL30 -filter \"SOR > 3.0\" --filter-name SOR3 -filter \"FS > 60.0\" --filter-name FS60 -filter \"MQ < 40.0\" --filter-name MQ40 -filter \"MQRankSum < -12.5\" --filter-name MQRankSum-12.5 -filter \"ReadPosRankSum < -8.0\" --filter-name ReadPosRankSum-8"
# the expression language on the cohort: strings, word operators, integers, lists (missing), negation of a string
# (missing), null, CHROM / POS / TYPE / N_ALLELES
add cohort "$(vf $C cohort.vcf.gz) -filter 'TYPE == \"SNP\" && FS > 2.0' --filter-name A -filter \"InbreedingCoeff lt -0.001 or not (MQ ge 40.0)\" --filter-name B -filter \"AF > 0.5\" --filter-name C -filter \"CHROM == 'chr20' && POS >= 10000000 && N_ALLELES > 2\" --filter-name D -filter \"DP > 100000\" --filter-name E -filter \"-ReadPosRankSum > 1.0\" --filter-name F -filter \"QUAL > 100000.0 || MQ0 != 0\" --filter-name G -filter \"null == QD || TYPE eq 'MIXED'\" --filter-name H"
add cohort_missing "$(vf $C cohort_missing.vcf.gz) -filter \"AF > 0.5\" --filter-name C -filter \"-ReadPosRankSum > 1.0\" --filter-name F --missing-values-evaluate-as-failing"
add cohort_invert "$(vf $C cohort_invert.vcf.gz) -filter \"QD > 10.0\" --filter-name Q -filter \"FS < 5.0 && SOR < 1.0\" --filter-name S --invert-filter-expression"
# an input with filters: kept and visible to expressions (FILTER, a filter name), or invalidated
add prev "$(vf $F prev.vcf.gz) -filter 'FILTER == \"1\" && QD < 10.0' --filter-name F1 -filter \"QD < 5.0 && LowQual\" --filter-name F2"
add invalidate "$(vf $F invalidate.vcf.gz) -filter \"QD < 5.0\" --filter-name LowQD --invalidate-previous-filters"
# masks: BED with an extension, a VCF, an interval list
add mask_bed "$(vf $S mask_bed.vcf.gz) --mask in/targets.bed --mask-extension 50"
add mask_vcf "$(vf $C mask_vcf.vcf.gz) --mask ../vqsr/resources/hapmap_3.3.hg38.vcf.gz --mask-name HAPMAP -filter \"QD < 2.0\" --filter-name QD2"
# repeated / empty INFO keys, a FLAG=0, unsorted samples, genotypes as read; a name given twice (the last wins)
add edge "$(vf in/edge.vcf.gz edge.vcf.gz) -filter \"FS > 1.0\" --filter-name X -filter \"FS > 3.0\" --filter-name X"
# GATK fails: an integer literal against a decimal string, a non-boolean value, an empty QD ('.')
add fail_int "$(vf $S fail_int.vcf.gz) -filter \"QD > 2\" --filter-name Q"
add fail_value "$(vf $S fail_value.vcf.gz) -filter \"FILTER\" --filter-name F"
add fail_dot "$(vf in/edge.vcf.gz fail_dot.vcf.gz) -filter \"QD > 2.0\" --filter-name Q"
# GATK fails: an INFO key the header does not declare; a mask without an index
add fail_undeclared "$(vf in/edge_undeclared.vcf.gz fail_undeclared.vcf.gz) -filter \"FS > 1.0\" --filter-name X"
add fail_mask_noindex "$(vf $S fail_mask_noindex.vcf.gz) --mask in/noidx.bed"

printf '%s\n' "${CASES[@]}" > truth/cases.txt
for c in "${CASES[@]}"; do
  echo "${CMD[$c]}" > truth/$c.cmd
  [ -f truth/$c.log ] && grep -q '^rc=' truth/$c.log && continue
  t=$(date +%s)
  eval "java -Xmx8g -jar $J ${CMD[$c]//@O@/truth}" > truth/$c.log 2>&1; rc=$?
  echo "rc=$rc $(( $(date +%s) - t ))s" | tee -a truth/$c.log | sed "s/^/$c /"
done
