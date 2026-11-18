#!/usr/bin/env bash
# Regression checks for the adversarial review's fixes (REVIEW_FINDINGS.md) that need no GATK run: each compares
# two runs of variant-rs that GATK/Picard treat as equivalent.   review_fixes.sh [OUTDIR]   (VARIANT_RS: the binary)
REPO=$(cd "$(dirname "$0")/../.." && pwd); G=${VARCALL_DATA:?set VARCALL_DATA to the test data directory}; X=${VARIANT_RS:-$REPO/target/release/variant-rs}; O=${1:-/tmp/variant-rs-review}; rm -rf $O; mkdir -p $O
data() { zcat $1 | grep -v '^##variant-rsCommandLine\|^##.*CommandLine='; }

# 1. block-compressed output names (htsjdk: .gz .gzip .bgz .bgzf, any case): BGZF + .tbi, same data as .vcf.gz
R=$G/hg002/ref/GRCh38_no_alt_analysis_set.fa; V=$G/hg002/group2/truth/warp_gvcf.g.vcf.gz
for e in vcf.gz vcf.bgz vcf.GZIP vcf.bgzf; do $X GenotypeGVCFs -R $R -V $V -O $O/gg.$e > $O/gg.$e.log 2>&1 || echo "GenotypeGVCFs -O .$e: exit $?"; done
for e in vcf.bgz vcf.GZIP vcf.bgzf; do
  ok=$( [ "$(head -c2 $O/gg.$e | od -An -tx1 | tr -d ' ')" = 1f8b ] && [ -s $O/gg.$e.tbi ] && cmp -s <(data $O/gg.vcf.gz) <(data $O/gg.$e) && echo OK || echo FAILED)
  echo "output name .$e: BGZF + .tbi, data as .vcf.gz: $ok"
done

# 2. MergeVcfs with a list file (Picard unrolls inputs without a VCF extension): same output as the inputs given directly
cd $G/mergevcfs
printf '%s\n' $PWD/in/v3.vcf.gz '' "  $PWD/in/v1.vcf.gz  " > $O/vcfs.list
$X MergeVcfs -I in/v3.vcf.gz -I in/v1.vcf.gz -O $O/mv_direct.vcf.gz > $O/mv_direct.log 2>&1
$X MergeVcfs -I $O/vcfs.list -O $O/mv_list.vcf.gz > $O/mv_list.log 2>&1 || echo "MergeVcfs with a list: exit $?"
echo "MergeVcfs list file: $(cmp -s <(data $O/mv_direct.vcf.gz) <(data $O/mv_list.vcf.gz) && echo 'same output as direct inputs: OK' || echo FAILED)"

# 3. Mutect2 germline resource with lowercase alleles (htsjdk upper-cases them): same calls as the original
H=$G/somatic/hcc1143/data; RG=$G/ref_hg38/Homo_sapiens_assembly38.fasta; REG=chr17:7000000-7700000
bcftools view -r $REG -Oz -o $O/gr.vcf.gz $H/chr17_af-only-gnomad_grch38.vcf.gz && tabix -f $O/gr.vcf.gz
zcat $O/gr.vcf.gz | awk -F'\t' 'BEGIN{OFS="\t"} /^#/ {print; next} {$4=tolower($4); $5=tolower($5); print}' | bgzip > $O/gr_lc.vcf.gz && tabix -f $O/gr_lc.vcf.gz
for v in gr gr_lc; do
  $X mutect2 -R $RG -I $H/tumor.bam -I $H/normal.bam -normal HCC1143_normal --germline-resource $O/$v.vcf.gz -L $REG \
    --af-of-alleles-not-in-resource 0.0001 -O $O/m2_$v.vcf.gz > $O/m2_$v.log 2>&1 || echo "Mutect2 with $v: exit $?"
done
n=$(zcat $O/m2_gr.vcf.gz 2>/dev/null | grep -vc "^#")
echo "germline resource, lowercase alleles: $n calls: $(cmp -s <(data $O/m2_gr.vcf.gz) <(data $O/m2_gr_lc.vcf.gz) && echo 'same as uppercase: OK' || echo FAILED)"

# 4. BaseRecalibrator --known-sites given as a plain VCF (the hg38 bundle's dbsnp138.vcf, as WARP passes it) or as
#    plain gzip: the same table as with the bgzipped files
B=$G/bqsr; KR=$B/resources
zcat $KR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz > $O/mills.vcf
zcat $KR/Homo_sapiens_assembly38.known_indels.vcf.gz | gzip -1 > $O/known_indels.vcf.gz
bq() { $X BaseRecalibrator -R $B/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa -I $B/truth/warp/na.oq.bam --use-original-qualities \
  -L chr20:10000000-20000000 --known-sites $KR/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites $1 --known-sites $2 -O $3 > $3.log 2>&1 || echo "BaseRecalibrator $3: exit $?"; }
bq $KR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz $KR/Homo_sapiens_assembly38.known_indels.vcf.gz $O/bq_bgzf.table
bq $O/mills.vcf $O/known_indels.vcf.gz $O/bq_plain.table
echo "known sites as plain VCF / plain gzip: $(cmp -s $O/bq_bgzf.table $O/bq_plain.table && echo 'same table as bgzipped: OK' || echo FAILED)"
