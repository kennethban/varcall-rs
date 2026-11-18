#!/usr/bin/env bash
# GATK part of run_all.sh only (step 4), so it can run while the Rust side is being fixed
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/hg002" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
REF=ref/GRCh38_no_alt_analysis_set.fa; BAM=bam/HG002.bam; P=${P:-5}
log() { echo "$(date +%H:%M:%S) $*" | tee -a run_all.log; }
gatk_run() {
  local m=$1 s=$2 e="" ext=vcf.gz
  [ $m = gvcf ] && { e="-ERC GVCF"; ext=g.vcf.gz; }
  local out=gatk/$m/$s.$ext
  [ -f $out.tbi ] && return 0
  local t0=$(date +%s)
  java -Xmx3g -XX:ParallelGCThreads=2 -jar $J HaplotypeCaller -R $REF -I $BAM -L plan/$s.intervals -O $out $e \
       --native-pair-hmm-threads 2 > gatk/logs/$m.$s.log 2>&1
  echo "$m $s rc=$? $(( $(date +%s) - t0 ))s" >> gatk/times.txt
}
export -f gatk_run; export J REF BAM
for m in vcf gvcf; do
  mkdir -p gatk/$m gatk/logs
  log "gatk $m start (P=$P)"
  tail -n +2 plan/plan.tsv | sort -t$'\t' -k6,6gr | cut -f1 | xargs -P $P -I{} bash -c "gatk_run $m {}"
  log "gatk $m done; failures: $(grep "^$m " gatk/times.txt | grep -vc 'rc=0')"
done
