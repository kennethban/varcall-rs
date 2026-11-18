#!/usr/bin/env bash
# variant-rs HaplotypeCaller with each truth case's GATK arguments (truth/<case>.cmd, see truth.sh) vs GATK.
# usage: cmp.sh [OUTDIR] [case...]
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/hg002/group2" || exit 1; D=$PWD

G=${VARIANT_RS:-$REPO/target/release/variant-rs}
O=${1:-$D/rust}; shift; mkdir -p $O
cases=${@:-$(ls truth/*.cmd | xargs -n1 basename | sed 's/\.cmd$//' | sort)}
for c in $cases; do
  args=$(cat truth/$c.cmd); ext=vcf.gz; [[ $args == *"-ERC GVCF"* ]] && ext=g.vcf.gz
  [ -f truth/$c.$ext.tbi ] || { echo "  $c: no GATK output (GATK refused it: $(grep -m1 'USER ERROR' truth/$c.log))"; continue; }
  (cd .. && /usr/bin/time -f "%e s" $G HaplotypeCaller -R ref/GRCh38_no_alt_analysis_set.fa -I bam/HG002.bam -O $O/$c.$ext $args \
      --native-pair-hmm-threads 2 > $O/$c.log 2>&1)
  echo "  $c: $($HARNESS/hg002/cmp1.sh $O/$c.$ext truth $ext) ($(tail -1 $O/$c.log))"
done
