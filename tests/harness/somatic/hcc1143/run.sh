#!/usr/bin/env bash
# HCC1143 (GATK somatic workshop, hg38 chr17plus): the full Mutect2 best-practice chain with GATK 4.7 (truth/) or
# variant-rs (rust/). usage: run.sh gatk|rust
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/hcc1143" || exit 1

R=$VARCALL_DATA/ref_hg38/Homo_sapiens_assembly38.fasta
D=data; O=$([ "$1" = gatk ] && echo truth || echo rust)
if [ "$1" = gatk ]; then X="java -Xmx16g -jar ${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}"; M2=Mutect2; GPS=GetPileupSummaries; CC=CalculateContamination; LROM=LearnReadOrientationModel; FMC=FilterMutectCalls
else X=${VARIANT_RS:-$REPO/target/release/variant-rs}; M2=mutect2; GPS=getpileupsummaries; CC=calculatecontamination; LROM=learnreadorientationmodel; FMC=filtermutectcalls; fi
t() { local t0=$(date +%s); "$@"; echo "  rc=$? $(( $(date +%s) - t0 ))s"; }
echo "Mutect2"; t $X $M2 -R $R -I $D/tumor.bam -I $D/normal.bam -normal HCC1143_normal -pon $D/chr17_m2pon.vcf.gz \
  --germline-resource $D/chr17_af-only-gnomad_grch38.vcf.gz -L $D/chr17plus.interval_list --f1r2-tar-gz $O/f1r2.tar.gz -O $O/m2.vcf > $O/m2.log 2>&1
for s in tumor normal; do echo "GetPileupSummaries $s"
  t $X $GPS -I $D/$s.bam -V $D/chr17_small_exac_common_3_grch38.vcf.gz -L $D/chr17_small_exac_common_3_grch38.vcf.gz -O $O/$s.pileups.table > $O/gps_$s.log 2>&1; done
echo "CalculateContamination"; t $X $CC -I $O/tumor.pileups.table -matched $O/normal.pileups.table -O $O/contamination.table --tumor-segmentation $O/segments.table > $O/cc.log 2>&1
echo "LearnReadOrientationModel"; t $X $LROM -I $O/f1r2.tar.gz -O $O/priors.tar.gz > $O/lrom.log 2>&1
echo "FilterMutectCalls"; t $X $FMC -R $R -V $O/m2.vcf --contamination-table $O/contamination.table --tumor-segmentation $O/segments.table \
  --ob-priors $O/priors.tar.gz -O $O/filtered.vcf > $O/fmc.log 2>&1
echo "done: $(grep -vc '^#' $O/filtered.vcf) records, $(grep -v '^#' $O/filtered.vcf | awk '$7=="PASS"' | wc -l) PASS"
