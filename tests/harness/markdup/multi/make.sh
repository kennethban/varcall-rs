#!/usr/bin/env bash
# MarkDuplicates over several coordinate-sorted inputs, as nf-core/sarek runs it (one BAM per lane / FASTQ chunk of a
# sample, `--REMOVE_DUPLICATES false --VALIDATION_STRINGENCY LENIENT --REFERENCE_SEQUENCE`), and the Picard truth.
# NA12878 chr20:10-15M (12 lanes, one library): lanes = one BAM per read group; lanes_rev = the same in reverse
# order; chunks = one lane split in 3 by read name (same read group in every file, as sarek's split_fastq makes);
# mixed = the chunks plus the other lanes; chr1_lanes = all of chr1, one BAM per lane (scale).
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/markdup/multi" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
REF=$VARCALL_DATA/bqsr/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
mkdir -p in truth
if [ ! -f in/done ]; then
  # strip the duplicate flags of the source (already marked) so the inputs look like fresh alignments
  samtools view -h $VARCALL_DATA/wgs/bam/NA12878.bam chr20:10000000-15000000 | \
    awk 'BEGIN{OFS="\t"} /^@/ {print; next} { $2 = and($2, compl(1024)); print }' | samtools view -b -o in/all.bam -
  samtools split -f 'in/lane_%!.bam' in/all.bam
  for f in in/lane_*.bam; do samtools index $f; done
  first=$(ls in/lane_*.bam | head -1)
  for k in 0 1 2; do
    samtools view -h $first | awk -v k=$k '/^@/ {print; next} { n = 0; for (i = 1; i <= length($1); i++) n += index("0123456789", substr($1, i, 1)); if (n % 3 == k) print }' | \
      samtools view -b -o in/chunk$k.bam - && samtools index in/chunk$k.bam
  done
  touch in/done
fi
if [ ! -f in/chr1/done ]; then
  mkdir -p in/chr1
  samtools view -h -@4 $VARCALL_DATA/wgs/bam/NA12878.bam chr1 | \
    awk 'BEGIN{OFS="\t"} /^@/ {print; next} { $2 = and($2, compl(1024)); print }' | samtools split -@4 -f 'in/chr1/lane_%!.bam' -
  for f in in/chr1/lane_*.bam; do samtools index $f; done
  touch in/chr1/done
fi
lanes=$(ls in/lane_*.bam); first=$(echo "$lanes" | head -1)
declare -A C
C[lanes]=$(for f in $lanes; do echo -n "--INPUT $f "; done)
C[lanes_rev]=$(for f in $(echo "$lanes" | sort -r); do echo -n "--INPUT $f "; done)
C[chunks]="--INPUT in/chunk0.bam --INPUT in/chunk1.bam --INPUT in/chunk2.bam"
C[chr1_lanes]=$(for f in in/chr1/lane_*.bam; do echo -n "--INPUT $f "; done)
C[mixed]="${C[chunks]} $(for f in $(echo "$lanes" | tail -n +2); do echo -n "--INPUT $f "; done)"
S="--TMP_DIR . --REFERENCE_SEQUENCE $REF --REMOVE_DUPLICATES false --VALIDATION_STRINGENCY LENIENT"
for c in "${!C[@]}"; do
  echo "${C[$c]} $S" | $HARNESS/cmdargs.sh -w > truth/$c.cmd
  [ -f truth/$c.metrics ] || java -Xmx8g -jar $J MarkDuplicates ${C[$c]} $S --OUTPUT truth/$c.bam --METRICS_FILE truth/$c.metrics > truth/$c.log 2>&1 &
done
wait
for c in "${!C[@]}"; do echo "$c $(grep -c . truth/$c.metrics) $(grep 'Elapsed time' truth/$c.log | sed 's/.*Elapsed time: //')"; done
