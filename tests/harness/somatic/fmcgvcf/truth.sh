#!/usr/bin/env bash
# FilterMutectCalls on a Mutect2 -ERC GVCF (review finding 21: <NON_REF> is symbolic; GATK drops its filter data,
# writes SITE for it in AS_FilterStatus, leaves it out of the clustering model, threshold and filtering stats).
# Input: GATK Mutect2 -ERC GVCF, HCC1143 tumor-only (GVCF mode takes one sample) chr17:1-20 Mb (no PoN / germline resource: GATK's GVCF
# mode fails on them, SomaticGenotypingEngine.getGermlineAltAlleleFrequencies index -1), plus f1r2 -> GATK
# LearnReadOrientationModel for the --ob-priors case.
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/fmcgvcf" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
H=../hcc1143/data; REF=$VARCALL_DATA/ref_hg38/Homo_sapiens_assembly38.fasta
mkdir -p truth
if [ ! -s in.g.vcf.stats ]; then
  java -Xmx8g -jar $J Mutect2 -R $REF -I $H/tumor.bam -L chr17:1-20000000 -ERC GVCF --f1r2-tar-gz f1r2.tar.gz \
    -O in.g.vcf --native-pair-hmm-threads 4 > m2.log 2>&1; echo "Mutect2 -ERC GVCF rc=$?"
  java -Xmx4g -jar $J LearnReadOrientationModel -I f1r2.tar.gz -O priors.tar.gz > lrom.log 2>&1; echo "LearnReadOrientationModel rc=$?"
fi
java -Xmx8g -jar $J FilterMutectCalls -R $REF -V in.g.vcf -O truth/plain.vcf > truth/plain.log 2>&1; echo "FilterMutectCalls plain rc=$?"
java -Xmx8g -jar $J FilterMutectCalls -R $REF -V in.g.vcf --ob-priors priors.tar.gz -O truth/ob.vcf > truth/ob.log 2>&1; echo "FilterMutectCalls ob rc=$?"
