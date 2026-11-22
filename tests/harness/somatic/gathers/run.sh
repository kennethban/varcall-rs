#!/usr/bin/env bash
# MergeMutectStats and GatherPileupSummaries as nf-core/sarek runs them after a scattered Mutect2 (gatk4/mergemutectstats,
# gatk4/gatherpileupsummaries). Pileup tables: GATK's per-piece tables of ../contam/sarek.sh (chr20 in three BED
# pieces, a multi-region BED) plus an empty one; stats: Mutect2 .stats of ../group3 and synthetic ones (fractions,
# 1e16) for JDK 17's compensated sum and Double.toString.  usage: run.sh truth|cmp [OUTDIR]
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/gathers" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}; G=${VARIANT_RS:-$REPO/target/release/variant-rs}
P=$VARCALL_DATA; DICT=$P/bqsr/ref/GRCh38_full_analysis_set_plus_decoy_hla.dict
T=../contam/sarek/truth; S=../group3/truth
mkdir -p in truth
if [ ! -f in/done ]; then
  # an empty table: a piece without resource sites
  java -jar $J GetPileupSummaries -I $P/wgs/bam/NA12878.bam -V ../contam/res/exac_chr20.hg38.vcf.gz -L chr20:1-1000 \
    -R ${DICT%.dict}.fa -O in/empty.pileups.table > in/empty.log 2>&1
  printf 'statistic\tvalue\ncallable\t0.1\n' > in/f1.stats; printf 'statistic\tvalue\ncallable\t0.2\n' > in/f2.stats
  printf 'statistic\tvalue\ncallable\t0.3\n' > in/f3.stats; printf 'statistic\tvalue\ncallable\t1.0E16\n' > in/big.stats
  printf 'statistic\tvalue\ncallable\t1.0\n' > in/one.stats; printf 'statistic\tvalue\ncallable\t-1.0E16\n' > in/negbig.stats
  touch in/done
fi
declare -A M G2
M[three]="--stats $S/to.vcf.stats --stats $S/af.vcf.stats --stats $S/gg.vcf.stats"
M[same_twice]="--stats $S/to.vcf.stats --stats $S/to.vcf.stats --stats $S/af.vcf.stats"
M[one]="--stats $S/mito.vcf.stats"
M[fractions]="--stats in/f1.stats --stats in/f2.stats --stats in/f3.stats"
M[cancel]="--stats in/big.stats --stats in/one.stats --stats in/negbig.stats"
G2[ordered]="--I $T/bam_p1.pileups.table --I $T/bam_p2.pileups.table --I $T/bam_p3.pileups.table"
G2[shuffled]="--I $T/bam_p3.pileups.table --I $T/bam_multi.pileups.table --I $T/bam_p1.pileups.table --I $T/bam_p2.pileups.table"
G2[with_empty]="--I in/empty.pileups.table --I $T/bam_p2.pileups.table --I $T/bam_p1.pileups.table"
G2[all_empty]="--I in/empty.pileups.table --I in/empty.pileups.table"
G2[tie]="--I $T/bam_multi_il.pileups.table --I $T/bam_string.pileups.table --I $T/bam_multi.pileups.table"
mode=${1:-cmp}; O=${2:-rust}; [ $mode = truth ] && O=truth; mkdir -p $O
for c in "${!M[@]}"; do
  [ $mode = truth ] && echo "${M[$c]}" | $HARNESS/cmdargs.sh -w > truth/mms_$c.cmd
  if [ $mode = truth ]; then java -jar $J MergeMutectStats ${M[$c]} -O $O/mms_$c.stats > $O/mms_$c.log 2>&1
  else $G --java-options "-Xmx2g" MergeMutectStats ${M[$c]} --output $O/mms_$c.stats --tmp-dir . > $O/mms_$c.log 2>&1; fi
  rc=$?; [ $mode = cmp ] && echo "  MergeMutectStats $c: rc=$rc $(cmp -s truth/mms_$c.stats $O/mms_$c.stats && echo identical || echo DIFFERENT) ($(tail -1 truth/mms_$c.stats | tr '\t' ' '))"
done
for c in "${!G2[@]}"; do
  if [ $mode = truth ]; then java -jar $J GatherPileupSummaries ${G2[$c]} --O $O/gps_$c.table --sequence-dictionary $DICT > $O/gps_$c.log 2>&1
  else $G --java-options "-Xmx2g" GatherPileupSummaries ${G2[$c]} --O $O/gps_$c.table --sequence-dictionary $DICT --tmp-dir . > $O/gps_$c.log 2>&1; fi
  rc=$?; [ $mode = cmp ] && echo "  GatherPileupSummaries $c: rc=$rc $(cmp -s truth/gps_$c.table $O/gps_$c.table && echo identical || echo DIFFERENT) ($(grep -vc '^#\|^contig' truth/gps_$c.table) sites)"
done
exit 0
