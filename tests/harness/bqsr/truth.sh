#!/usr/bin/env bash
# GATK 4.7 BQSR reference outputs for the Rust port (BaseRecalibrator table + ApplyBQSR BAM), default arguments,
# known sites = Broad hg38 bundle (dbSNP 138, Mills + 1000G gold standard indels, known_indels).
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/bqsr" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
REF=ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
R=resources
KS="--known-sites $R/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites $R/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --known-sites $R/Homo_sapiens_assembly38.known_indels.vcf.gz"
declare -A BAM=([pcrplus]=$VARCALL_DATA/hg002_pcrplus/bam/HG002.pcrplus.bam [na12878]=$VARCALL_DATA/wgs/bam/NA12878.bam)
log() { echo "$(date +%H:%M:%S) $*" | tee -a truth/truth.log; }
run() { # name args...
  local n=$1; shift; local t0=$(date +%s)
  [ -f truth/$n.ok ] && return 0
  /usr/bin/time -v java -Xmx8g -XX:ParallelGCThreads=2 -jar $J "$@" > truth/$n.log 2>&1 && touch truth/$n.ok
  log "$n rc=$? $(( $(date +%s) - t0 ))s"
}
for s in pcrplus na12878; do
  run $s.chr20.table BaseRecalibrator -R $REF -I ${BAM[$s]} -L chr20 $KS -O truth/$s.chr20.recal.table
  run $s.chr20.apply ApplyBQSR -R $REF -I ${BAM[$s]} -L chr20 --bqsr-recal-file truth/$s.chr20.recal.table -O truth/$s.chr20.bqsr.bam
done
run pcrplus.wgs.table BaseRecalibrator -R $REF -I ${BAM[pcrplus]} $KS -O truth/pcrplus.wgs.recal.table
run pcrplus.wgs.apply ApplyBQSR -R $REF -I ${BAM[pcrplus]} --bqsr-recal-file truth/pcrplus.wgs.recal.table -O truth/pcrplus.wgs.bqsr.bam
log "all done"
