#!/usr/bin/env bash
# GATK 4.7 truth for M3's joint steps (truth.sh runs it), each from the pipeline's own input to that step:
#   GenomicsDBImport of the pipeline's GVCFs with sarek's arguments (intervals imported side by side: the workspace's
#   content does not depend on it), then GenotypeGVCFs once per calling interval (the pipeline's scatter_intervals,
#   sarek's per-interval tasks), records in interval order under the first run's header;
#   VariantRecalibrator SNP and INDEL on the pipeline's joint VCF; ApplyVQSR SNP then INDEL with its recal files.
#   usage: truth_joint.sh RUN_DIR        (JOBS: GATK runs at once, default 96; OUT: the pipeline output folder in RUN_DIR, default out)
set -uo pipefail
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
: "${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}"
O=$(realpath "${1:?usage: truth_joint.sh RUN_DIR}"); P=$O/${OUT:-out}; T=$O/truth
mkdir -p $T/log $T/tmp $T/gg
JOBS=${JOBS:-96}
D=$VARCALL_DATA; R=$D/ref_hg38/Homo_sapiens_assembly38.fasta
DBSNP=$D/bqsr/resources/Homo_sapiens_assembly38.dbsnp138.vcf.gz
IV=$D/somatic/wgs_ivs/wgs_calling_regions.hg38.interval_list
B=$D/bqsr/resources; V=$D/vqsr/resources; J=$P/germline/joint
SAMPLES="NA12878 NA12891 NA12892"
GB="java -Xmx48g -Djava.io.tmpdir=$T/tmp -jar $GATK_JAR"
export G="java -Xmx4g -XX:ParallelGCThreads=2 -Djava.io.tmpdir=$T/tmp -jar $GATK_JAR"
export T R DBSNP
t() {
  local name=$1; shift
  [ -e $T/log/$name.ok ] && return 0
  local s=$(date +%s)
  "$@" > $T/log/$name.log 2>&1 || { echo "FAILED $name"; return 1; }
  touch $T/log/$name.ok; echo "$name $(( $(date +%s) - s ))s"
}
export -f t
grep -v '^@' $IV | awk '{print $1":"$2"-"$3}' > $T/intervals.txt

# VQSR, from the pipeline's joint VCF and recal files, beside the import and GenotypeGVCFs
vqsr() {
  local snp="--resource:hapmap,known=false,training=true,truth=true,prior=15.0 $V/hapmap_3.3.hg38.vcf.gz --resource:omni,known=false,training=true,truth=true,prior=12.0 $V/1000G_omni2.5.hg38.vcf.gz --resource:1000G,known=false,training=true,truth=false,prior=10.0 $V/1000G_phase1.snps.high_confidence.hg38.vcf.gz --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP"
  local indel="--resource:mills,known=false,training=true,truth=true,prior=12.0 $B/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP"
  t vr_snp $GB VariantRecalibrator -V $J/joint_germline.vcf.gz -O $T/joint_germline_SNP.recal.vcf.gz \
    --tranches-file $T/joint_germline_SNP.tranches -R $R --tmp-dir $T/tmp $snp \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP &
  t vr_indel $GB VariantRecalibrator -V $J/joint_germline.vcf.gz -O $T/joint_germline_INDEL.recal.vcf.gz \
    --tranches-file $T/joint_germline_INDEL.tranches -R $R --tmp-dir $T/tmp $indel \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -mode INDEL &
  t applyvqsr_snp $GB ApplyVQSR -V $J/joint_germline.vcf.gz -O $T/joint_germline_SNP.vcf.gz -R $R \
    --tranches-file $J/joint_germline_SNP.tranches --recal-file $J/joint_germline_SNP.recal.vcf.gz --tmp-dir $T/tmp \
    --truth-sensitivity-filter-level 99.9 -mode SNP &&
  t applyvqsr_indel $GB ApplyVQSR -V $T/joint_germline_SNP.vcf.gz -O $T/joint_germline_recalibrated.vcf.gz -R $R \
    --tranches-file $J/joint_germline_INDEL.tranches --recal-file $J/joint_germline_INDEL.recal.vcf.gz \
    --tmp-dir $T/tmp --truth-sensitivity-filter-level 99.9 -mode INDEL
  wait
}
vqsr &

vs=$(for s in $SAMPLES; do echo -n "-V $P/germline/$s/$s.g.vcf.gz "; done)
[ -e $T/log/gdbimport.ok ] || rm -rf "${T:?}/joint_germline.gdb"
t gdbimport java -Xmx64g -Djava.io.tmpdir=$T/tmp -jar $GATK_JAR GenomicsDBImport $vs --genomicsdb-workspace-path $T/joint_germline.gdb \
  -L $IV --tmp-dir $T/tmp --genomicsdb-shared-posixfs-optimizations true --bypass-feature-reader \
  --max-num-intervals-to-import-in-parallel 32 || exit 1
i=0
while read -r iv; do
  i=$((i + 1))
  echo "t gg_$i \$G GenotypeGVCFs -R $R -V gendb://$T/joint_germline.gdb -L $iv --dbsnp $DBSNP -O $T/gg/$i.vcf.gz"
done < $T/intervals.txt | xargs -P $JOBS -I{} bash -c '{}'
n=$(wc -l < $T/intervals.txt)
missing=$(for i in $(seq 1 $n); do [ -e $T/log/gg_$i.ok ] || echo $i; done | wc -l)
[ $missing = 0 ] || { echo "FAILED genotypegvcfs: $missing intervals"; exit 1; }
{ zcat $T/gg/1.vcf.gz | grep '^#'; for i in $(seq 1 $n); do zcat $T/gg/$i.vcf.gz | grep -v '^#'; done; } | bgzip -@ 8 > $T/joint_germline.vcf.gz &&
  tabix -p vcf $T/joint_germline.vcf.gz && echo "joint_germline.vcf.gz: $n GATK runs, $(zcat $T/joint_germline.vcf.gz | grep -vc '^#') records"
wait
