# The variant-rs part of the self-contained test, as GATK command lines: sourced by tests/selftest.sh (X =
# variant-rs) and tests/make_truth.sh (X = java -jar gatk.jar), and run in the data directory with relative paths
# (GATK writes its arguments into some outputs). The caller defines run NAME COMMAND... and sets HMM (variant-rs:
# --pairhmm cpu or cuda; GATK: empty).
variant_steps() {
  run md_GERM   $X MarkDuplicates -I GERM.L1.bam -I GERM.L2.bam -O GERM.md.bam -M GERM.md.metrics
  run md_TUMOR  $X MarkDuplicates -I TUMOR.bam -O TUMOR.md.bam -M TUMOR.md.metrics
  run md_NORMAL $X MarkDuplicates -I NORMAL.bam -O NORMAL.md.bam -M NORMAL.md.metrics
  for s in GERM TUMOR NORMAL; do
    run bqsr_$s  $X BaseRecalibrator -R ref.fa -I $s.md.bam --known-sites known_sites.vcf.gz -O $s.recal.table
    run apply_$s $X ApplyBQSR -R ref.fa -I $s.md.bam --bqsr-recal-file $s.recal.table -O $s.recal.bam
  done
  run hc      $X HaplotypeCaller -R ref.fa -I GERM.recal.bam -O GERM.vcf.gz $HMM
  run hc_gvcf $X HaplotypeCaller -R ref.fa -I GERM.recal.bam -O GERM.g.vcf.gz -ERC GVCF $HMM
  run gg      $X GenotypeGVCFs -R ref.fa -V GERM.g.vcf.gz -O GERM.gg.vcf.gz
  run m2      $X Mutect2 -R ref.fa -I TUMOR.recal.bam -I NORMAL.recal.bam -normal NORMAL \
                 --germline-resource germline_resource.vcf.gz --f1r2-tar-gz f1r2.tar.gz -O somatic.vcf.gz $HMM
  run lrom    $X LearnReadOrientationModel -I f1r2.tar.gz -O priors.tar.gz
  for s in TUMOR NORMAL; do
    run gps_$s $X GetPileupSummaries -I $s.recal.bam -V germline_resource.vcf.gz -L germline_resource.vcf.gz -O $s.pileups.table
  done
  run cc      $X CalculateContamination -I TUMOR.pileups.table -matched NORMAL.pileups.table \
                 -O contamination.table --tumor-segmentation segments.table
  run fmc     $X FilterMutectCalls -R ref.fa -V somatic.vcf.gz --contamination-table contamination.table \
                 --tumor-segmentation segments.table --ob-priors priors.tar.gz -O somatic.filtered.vcf.gz
  # VQSR on the germline calls, the known sites as the training / truth resource (small: 2 / 1 Gaussians)
  run vqsr_snp   $X VariantRecalibrator -R ref.fa -V GERM.vcf.gz -O GERM.snp.recal.vcf.gz --tranches-file GERM.snp.tranches \
                    -an QD -an FS -an SOR -an ReadPosRankSum -mode SNP --max-gaussians 2 \
                    --resource:known,known=true,training=true,truth=true,prior=12 known_sites.vcf.gz
  run vqsr_indel $X VariantRecalibrator -R ref.fa -V GERM.vcf.gz -O GERM.indel.recal.vcf.gz --tranches-file GERM.indel.tranches \
                    -an QD -an FS -an SOR -mode INDEL --max-gaussians 1 \
                    --resource:known,known=true,training=true,truth=true,prior=12 known_sites.vcf.gz
  run applyvqsr_snp   $X ApplyVQSR -R ref.fa -V GERM.vcf.gz --recal-file GERM.snp.recal.vcf.gz --tranches-file GERM.snp.tranches \
                         --truth-sensitivity-filter-level 99.0 -mode SNP -O GERM.vqsr_snp.vcf.gz
  run applyvqsr_indel $X ApplyVQSR -R ref.fa -V GERM.vqsr_snp.vcf.gz --recal-file GERM.indel.recal.vcf.gz \
                         --tranches-file GERM.indel.tranches --truth-sensitivity-filter-level 99.0 -mode INDEL -O GERM.vqsr.vcf.gz
  # hard filters as GATK's examples and WARP write them, then WARP's sites-only step
  run vfilter $X VariantFiltration -V GERM.vcf.gz -O GERM.vf.vcf.gz \
                 --filter-expression "QD < 2.0 || FS > 30.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -3.0 || ReadPosRankSum < -3.0" \
                 --filter-name HardFiltered --filter-expression "QUAL < 30.0" --filter-name QUAL30 \
                 --filter-expression 'TYPE == "INDEL" && FS > 10.0' --filter-name IndelFS
  run sitesonly $X MakeSitesOnlyVcf -I GERM.vf.vcf.gz -O GERM.sites.vcf.gz
}

# What is compared, and how: BAM records (bamdigest.py), VCF data lines, metrics without '#' lines, files as they are
VARIANT_BAMS="GERM.md.bam TUMOR.md.bam NORMAL.md.bam GERM.recal.bam TUMOR.recal.bam NORMAL.recal.bam"
VARIANT_VCFS="GERM.vcf.gz GERM.g.vcf.gz GERM.gg.vcf.gz somatic.vcf.gz somatic.filtered.vcf.gz GERM.snp.recal.vcf.gz
  GERM.indel.recal.vcf.gz GERM.vqsr.vcf.gz GERM.vf.vcf.gz GERM.sites.vcf.gz"
VARIANT_METRICS="GERM.md.metrics TUMOR.md.metrics NORMAL.md.metrics"
VARIANT_FILES="GERM.recal.table TUMOR.recal.table NORMAL.recal.table somatic.vcf.gz.stats TUMOR.pileups.table
  NORMAL.pileups.table contamination.table segments.table somatic.filtered.vcf.gz.filteringStats.tsv GERM.snp.tranches
  GERM.indel.tranches"
