#!/usr/bin/env bash
# GATK 4.7 truth for M3 (run.sh), step by step: each GATK step reads the pipeline's own input to that step, so all
# steps run side by side (one chain from the lane BAMs would take a day):
#   per sample: MarkDuplicates on the lane BAMs; BaseRecalibrator on the pipeline's MarkDuplicates BAM (whole genome,
#   as the pipeline); the GVCF from the pipeline's BAM and table: ApplyBQSR per contig, then one HaplotypeCaller
#   run per calling interval (the pipeline's scatter), records in interval order under the first run's header;
#   joint (truth_joint.sh): from the pipeline's GVCFs, GenomicsDBImport (sarek's arguments) -> GenotypeGVCFs per
#   interval; VariantRecalibrator SNP and INDEL on the pipeline's joint VCF; ApplyVQSR with the pipeline's recal files.
#   usage: truth.sh RUN_DIR        (JOBS: GATK runs at once for the per-interval steps, default 96; OUT: the
#   pipeline output folder in RUN_DIR, default out)
# Resumable: finished steps leave log/<step>.ok. Every Java temp file goes to RUN_DIR/truth/tmp.
set -uo pipefail
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
: "${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}"
O=$(realpath "${1:?usage: truth.sh RUN_DIR}"); P=$O/${OUT:-out}; T=$O/truth
mkdir -p $T/log $T/tmp $T/bqsr $T/hc
JOBS=${JOBS:-96}
D=$VARCALL_DATA; R=$D/ref_hg38/Homo_sapiens_assembly38.fasta; B=$D/bqsr/resources
DBSNP=$B/Homo_sapiens_assembly38.dbsnp138.vcf.gz
IV=$D/somatic/wgs_ivs/wgs_calling_regions.hg38.interval_list
SAMPLES="NA12878 NA12891 NA12892"
# small heaps for the many per-interval runs, large ones for the whole-genome steps
export G="java -Xmx4g -XX:ParallelGCThreads=2 -Djava.io.tmpdir=$T/tmp -jar $GATK_JAR"
export GB="java -Xmx48g -Djava.io.tmpdir=$T/tmp -jar $GATK_JAR"
export T R DBSNP
# t NAME CMD...: runs CMD once (log/NAME.ok when it succeeded), its output in log/NAME.log
t() {
  local name=$1; shift
  [ -e $T/log/$name.ok ] && return 0
  local s=$(date +%s)
  "$@" > $T/log/$name.log 2>&1 || { echo "FAILED $name"; return 1; }
  touch $T/log/$name.ok; echo "$name $(( $(date +%s) - s ))s"
}
export -f t
# the calling intervals as chr:start-end (interval_list: 1-based, inclusive), in file order
grep -v '^@' $IV | awk '{print $1":"$2"-"$3}' > $T/intervals.txt
cut -d: -f1 $T/intervals.txt | uniq > $T/contigs.txt

date "+start %T"
# the joint steps (truth_joint.sh), in the background
"$(dirname "$0")/truth_joint.sh" "$O" &
for s in $SAMPLES; do
  ins=$(ls $P/lanes/$s.*.bam | sort | sed 's/^/-I /' | tr '\n' ' ')
  t md_$s $GB MarkDuplicates $ins -O $T/$s.md.bam -M $T/$s.md.metrics --CREATE_INDEX true --TMP_DIR $T/tmp &
  t bqsr_$s $GB BaseRecalibrator -R $R -I $P/preprocessing/$s/$s.md.bam \
    --known-sites $DBSNP --known-sites $B/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --known-sites $B/Homo_sapiens_assembly38.known_indels.vcf.gz -O $T/$s.recal.table &
done
# the GVCFs: ApplyBQSR per (sample, contig), then HaplotypeCaller per (sample, interval)
for s in $SAMPLES; do
  while read -r c; do
    echo "t abqsr_${s}_$c \$G ApplyBQSR -R $R -I $P/preprocessing/$s/$s.md.bam --bqsr-recal-file $P/preprocessing/$s/$s.recal.table -L $c -O $T/bqsr/$s.$c.bam"
  done < $T/contigs.txt
done | xargs -P $JOBS -I{} bash -c '{}'
for s in $SAMPLES; do
  i=0
  while read -r iv; do
    i=$((i + 1))
    echo "t hc_${s}_$i \$G HaplotypeCaller -R $R -I $T/bqsr/$s.${iv%%:*}.bam -L $iv --dbsnp $DBSNP -ERC GVCF -O $T/hc/$s.$i.g.vcf.gz"
  done < $T/intervals.txt
done | xargs -P $JOBS -I{} bash -c '{}'
n=$(wc -l < $T/intervals.txt)
for s in $SAMPLES; do
  [ -e $T/log/concat_$s.ok ] && continue
  parts=$(for i in $(seq 1 $n); do echo $T/hc/$s.$i.g.vcf.gz; done)
  missing=$(for p in $parts; do [ -s $p ] || echo $p; done | wc -l)
  [ $missing = 0 ] || { echo "FAILED concat_$s: $missing GVCF parts missing"; continue; }
  { zcat $(echo "$parts" | head -1) | grep '^#'; for p in $parts; do zcat $p | grep -v '^#'; done; } | bgzip -@ 8 > $T/$s.g.vcf.gz &&
    tabix -p vcf $T/$s.g.vcf.gz && touch $T/log/concat_$s.ok && echo "concat_$s done"
done
wait
date "+all done %T"
for l in $T/log/*.log; do [ -e ${l%.log}.ok ] || echo "not finished: $l"; done
