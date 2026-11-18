#!/usr/bin/env bash
# GATK truth for FilterMutectCalls --mitochondria-mode (and --min-allele-fraction, --max-alt-allele-count,
# --f-score-beta) on the Mutect2 truth outputs of truth.sh. warp1/warp2: AlignAndCall.wdl's InitialFilter and
# FilterContamination (warp2 re-filters warp1's output, as WARP does; WARP's VariantFiltration step in between is
# not run). Outputs fmc/<case>.vcf (+ .filteringStats.tsv); each case's arguments in fmc/<case>.cmd.
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/group3" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
HR=$VARCALL_DATA/ref_hg38/Homo_sapiens_assembly38.fasta
WR=$VARCALL_DATA/bqsr/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
mkdir -p fmc
declare -A C
C[warp1]="-V truth/mito_warp.vcf -R $WR --stats truth/mito_warp.vcf.stats --max-alt-allele-count 4 --mitochondria-mode --min-allele-fraction 0 --contamination-estimate 0.0"
C[warp2]="-V fmc/warp1.vcf -R $WR --stats truth/mito_warp.vcf.stats --max-alt-allele-count 4 --mitochondria-mode --min-allele-fraction 0.03 --f-score-beta 1.5 --contamination-estimate 0.02"
C[mito]="-V truth/mito.vcf -R $WR --mitochondria-mode"
C[to_mito]="-V truth/to.vcf -R $HR --mitochondria-mode --f-score-beta 0.5 --min-allele-fraction 0.05"
C[to_opts]="-V truth/to.vcf -R $HR --max-alt-allele-count 2 --min-allele-fraction 0.05 --f-score-beta 2.0"
C[tn_opts]="-V truth/gg.vcf -R $HR --max-alt-allele-count 2 --min-allele-fraction 0.02 --f-score-beta 0.7"
run() {
  local c=$1 args=$2
  [ -f fmc/$c.vcf.filteringStats.tsv ] && return 0
  echo "$args" > fmc/$c.cmd
  java -Xmx4g -XX:ParallelGCThreads=2 -jar $J FilterMutectCalls $args -O fmc/$c.vcf > fmc/$c.log 2>&1
  echo "$c rc=$?"
}
for c in mito to_mito to_opts tn_opts warp1; do run $c "${C[$c]}" & done
wait
run warp2 "${C[warp2]}"
