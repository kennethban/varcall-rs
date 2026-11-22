#!/usr/bin/env bash
# GetPileupSummaries as nf-core/sarek scatters it: `--intervals <piece>.bed` over a resource VCF (sarek's
# gatk4/getpileupsummaries), NA12878 chr20 (BAM; GATK cannot read this CRAM: htsjdk "unmapped slice" crai error), resource res/exac_chr20.hg38.vcf.gz. Pieces: chr20 cut in
# three BED files at resource sites (one piece ends on a site, the next starts on the following one; a gap between
# the 2nd and 3rd), a BED of several small regions, the same territory as an .interval_list and as contig:start-end.
# Truth: GATK 4.7; variant-rs compared table for table.  usage: sarek.sh truth|cmp [OUTDIR]
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/contam" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}; G=${VARIANT_RS:-$REPO/target/release/variant-rs}
P=$VARCALL_DATA; REF=$P/bqsr/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
V=res/exac_chr20.hg38.vcf.gz
mkdir -p sarek/in sarek/truth
if [ ! -f sarek/in/done ]; then
  # site positions to cut at: the 1/3 and 2/3 sites of the resource on chr20
  zcat $V | grep -v '^#' | awk '$1=="chr20" {print $2}' > sarek/in/pos.txt
  n=$(wc -l < sarek/in/pos.txt); a=$(sed -n "$((n / 3))p" sarek/in/pos.txt); b=$(sed -n "$((n / 3 + 1))p" sarek/in/pos.txt)
  c=$(sed -n "$((2 * n / 3))p" sarek/in/pos.txt); d=$(sed -n "$((2 * n / 3 + 5))p" sarek/in/pos.txt)
  printf 'chr20\t0\t%d\n' $a > sarek/in/p1.bed                       # ends on a site (BED end = 1-based last)
  printf 'chr20\t%d\t%d\n' $((b - 1)) $c > sarek/in/p2.bed           # starts on the next site
  printf 'chr20\t%d\t64444167\n' $((d - 1)) > sarek/in/p3.bed        # gap of a few sites before it
  printf 'chr20\t1000000\t1500000\nchr20\t2000000\t2000100\nchr20\t5000000\t5600000\n' > sarek/in/multi.bed
  { grep '^@SQ' ${REF%.fa}.dict; awk 'BEGIN{OFS="\t"} {print $1, $2 + 1, $3, "+", "."}' sarek/in/multi.bed; } > sarek/in/multi.interval_list
  touch sarek/in/done
fi
declare -A C
for p in p1 p2 p3 multi; do C[bam_$p]="-I $P/wgs/bam/NA12878.bam --intervals sarek/in/$p.bed"; done
C[bam_multi_il]="-I $P/wgs/bam/NA12878.bam --intervals sarek/in/multi.interval_list"
C[bam_string]="-I $P/wgs/bam/NA12878.bam --intervals chr20:1000000-1500000"
case ${1:-cmp} in
truth)
  for c in "${!C[@]}"; do
    echo "${C[$c]}" | $HARNESS/cmdargs.sh -w > sarek/truth/$c.cmd
    [ -f sarek/truth/$c.pileups.table ] || java -Xmx4g -jar $J GetPileupSummaries ${C[$c]} --variant $V --reference $REF \
      --output sarek/truth/$c.pileups.table --tmp-dir . > sarek/truth/$c.log 2>&1 &
  done
  wait
  for c in "${!C[@]}"; do echo "$c $(grep -vc '^#\|^contig' sarek/truth/$c.pileups.table) sites"; done ;;
cmp)
  O=${2:-sarek/rust}; mkdir -p $O
  for c in $(ls sarek/truth/*.cmd | xargs -n1 basename | sed 's/\.cmd$//' | sort); do
    $G --java-options "-Xmx4g -XX:-UsePerfData" GetPileupSummaries $($HARNESS/cmdargs.sh sarek/truth/$c.cmd) --variant $V --reference $REF \
      --output $O/$c.pileups.table --tmp-dir . > $O/$c.log 2>&1; rc=$?
    echo "  $c: rc=$rc $(cmp -s sarek/truth/$c.pileups.table $O/$c.pileups.table && echo identical || echo DIFFERENT) ($(grep -vc '^#\|^contig' sarek/truth/$c.pileups.table) sites)"
  done ;;
esac
