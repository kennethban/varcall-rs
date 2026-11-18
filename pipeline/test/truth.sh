#!/usr/bin/env bash
# The pipeline's chain with GATK 4.7 (Picard MarkDuplicates included).  usage: truth.sh germline|somatic|somatic_bam
#   germline, somatic: from the pipeline's own lane BAMs (out/SCEN/lanes, --save_lane_bams: alignment is not a GATK step)
#   somatic_bam: from data/wes (SEQC2 exome pair, chr22; read-group SM WES_FD_*, not the samplesheet's FD_*), as the
#   pipeline's --step markduplicates; Mutect2 on the pair, HaplotypeCaller on the normal, over the chr22 target BED
set -eo pipefail
cd "$(dirname "$0")"
scen=$1; L=out/$scen/lanes; O=truth/$scen; mkdir -p $O
VD=${VARCALL_DATA:?set VARCALL_DATA to the test data directory}
G="java -Xmx16g -jar ${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}"
R=$VD/ref_hg38/Homo_sapiens_assembly38.fasta
B=$VD/bqsr/resources; H=$VD/somatic/hcc1143/data; W=$VD/somatic/seqc2/data
if [ $scen = somatic_bam ]; then L=data/wes; R=$W/GRCh38.d1.vd1.fa; fi
KS="--known-sites $B/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites $B/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --known-sites $B/Homo_sapiens_assembly38.known_indels.vcf.gz"
t() { local name=$1; shift; local s=$(date +%s); "$@" > $O/$name.log 2>&1 || { echo "FAILED $name"; exit 1; }; echo "$name $(( $(date +%s) - s ))s"; }
# hard filters as the pipeline's HARDFILTER (GATK's germline short-variant recommendations): hf SAMPLE
hf() {
  local s=$1
  t sel_snp $G SelectVariants -V $O/$s.vcf.gz -select-type SNP -O $O/snps.vcf.gz
  t sel_indel $G SelectVariants -V $O/$s.vcf.gz -select-type INDEL -O $O/indels.vcf.gz
  t filt_snp $G VariantFiltration -V $O/snps.vcf.gz -O $O/snps.f.vcf.gz \
    -filter "QD < 2.0" --filter-name QD2 -filter "QUAL < 30.0" --filter-name QUAL30 -filter "SOR > 3.0" --filter-name SOR3 \
    -filter "FS > 60.0" --filter-name FS60 -filter "MQ < 40.0" --filter-name MQ40 \
    -filter "MQRankSum < -12.5" --filter-name MQRankSum-12.5 -filter "ReadPosRankSum < -8.0" --filter-name ReadPosRankSum-8
  t filt_indel $G VariantFiltration -V $O/indels.vcf.gz -O $O/indels.f.vcf.gz \
    -filter "QD < 2.0" --filter-name QD2 -filter "QUAL < 30.0" --filter-name QUAL30 -filter "FS > 200.0" --filter-name FS200 \
    -filter "ReadPosRankSum < -20.0" --filter-name ReadPosRankSum-20
  t merge $G MergeVcfs -I $O/snps.f.vcf.gz -I $O/indels.f.vcf.gz -O $O/$s.filtered.vcf.gz
}
samples=$(ls $L/*.bam | xargs -n1 basename | cut -d. -f1 | sort -u)
for s in $samples; do
  ins=$(ls $L/$s.*.bam | sort | sed 's/^/-I /' | tr '\n' ' ')
  t md_$s $G MarkDuplicates $ins -O $O/$s.md.bam -M $O/$s.md.metrics --CREATE_INDEX true
  t bqsr_$s $G BaseRecalibrator -R $R -I $O/$s.md.bam $KS -O $O/$s.recal.table
  t apply_$s $G ApplyBQSR -R $R -I $O/$s.md.bam --bqsr-recal-file $O/$s.recal.table -O $O/$s.recal.bam
done
if [ $scen = germline ]; then
  s=NA12878
  # one GATK run per interval (the pipeline's --scatter_intervals, sarek's HaplotypeCaller scatter); the runs'
  # records in interval order under the first run's header
  i=0; parts=()
  while read -r iv; do
    i=$((i+1)); t hc_$i $G HaplotypeCaller -R $R -I $O/$s.recal.bam -L $iv -O $O/hc_$i.vcf.gz; parts+=($O/hc_$i.vcf.gz)
  done < germline.intervals
  { zcat ${parts[0]} | grep "^#"; for p in "${parts[@]}"; do zcat $p | grep -v "^#"; done; } | bgzip > $O/$s.vcf.gz
  tabix -f -p vcf $O/$s.vcf.gz
  hf $s
elif [ $scen = somatic ]; then
  T=HCC1143_tumor; N=HCC1143_normal; id=${T}_vs_$N; SITES=$H/chr17_small_exac_common_3_grch38.vcf.gz
  t m2 $G Mutect2 -R $R -I $O/$T.recal.bam -I $O/$N.recal.bam -normal $N -L $H/chr17plus.interval_list \
    --germline-resource $H/chr17_af-only-gnomad_grch38.vcf.gz -pon $H/chr17_m2pon.vcf.gz --f1r2-tar-gz $O/$id.f1r2.tar.gz -O $O/$id.mutect2.vcf.gz
  for s in $T $N; do t gps_$s $G GetPileupSummaries -I $O/$s.recal.bam -V $SITES -L $SITES -O $O/$s.pileups.table; done
  t cc $G CalculateContamination -I $O/$T.pileups.table -matched $O/$N.pileups.table -O $O/$id.contamination.table --tumor-segmentation $O/$id.segments.table
  t lrom $G LearnReadOrientationModel -I $O/$id.f1r2.tar.gz -O $O/$id.priors.tar.gz
  t fmc $G FilterMutectCalls -R $R -V $O/$id.mutect2.vcf.gz --contamination-table $O/$id.contamination.table \
    --tumor-segmentation $O/$id.segments.table --ob-priors $O/$id.priors.tar.gz -O $O/$id.filtered.vcf.gz
else
  T=FD_T_1; N=FD_N_1; id=${T}_vs_$N; BED=data/wes/targets.chr22.bed; GN=$W/af-only-gnomad.hg38.vcf.gz
  # -normal: the read groups' SM (sarek's patient_sample convention), not the samplesheet name
  t m2 $G Mutect2 -R $R -I $O/$T.recal.bam -I $O/$N.recal.bam -normal WES_FD_N_1 -L $BED \
    --germline-resource $GN -pon $W/1000g_pon.hg38.vcf.gz --f1r2-tar-gz $O/$id.f1r2.tar.gz -O $O/$id.mutect2.vcf.gz
  for s in $T $N; do t gps_$s $G GetPileupSummaries -I $O/$s.recal.bam -V $GN -L $BED -O $O/$s.pileups.table; done
  t cc $G CalculateContamination -I $O/$T.pileups.table -matched $O/$N.pileups.table -O $O/$id.contamination.table --tumor-segmentation $O/$id.segments.table
  t lrom $G LearnReadOrientationModel -I $O/$id.f1r2.tar.gz -O $O/$id.priors.tar.gz
  t fmc $G FilterMutectCalls -R $R -V $O/$id.mutect2.vcf.gz --contamination-table $O/$id.contamination.table \
    --tumor-segmentation $O/$id.segments.table --ob-priors $O/$id.priors.tar.gz -O $O/$id.filtered.vcf.gz
  # HaplotypeCaller on the normal: one GATK run over the targets (scatter_intervals auto: 6,644 regions > 1,000)
  t hc $G HaplotypeCaller -R $R -I $O/$N.recal.bam -L $BED -O $O/$N.vcf.gz
  hf $N
fi
