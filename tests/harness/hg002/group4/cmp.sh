#!/usr/bin/env bash
# variant-rs HaplotypeCaller with each truth case's GATK arguments (truth/<case>.cmd, see truth.sh) vs GATK; then the
# refusal of the same reads in two inputs.  usage: cmp.sh [OUTDIR] [case...]
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/hg002/group4" || exit 1; D=$PWD

G=${VARIANT_RS:-$REPO/target/release/variant-rs}
O=${1:-$D/rust}; shift; mkdir -p $O; O=$(realpath $O)
cases=${@:-$(ls truth/*.cmd | xargs -n1 basename | sed 's/\.cmd$//' | sort)}
for c in $cases; do
  args=$($HARNESS/cmdargs.sh truth/$c.cmd); ext=vcf.gz; [[ $args == *"-ERC GVCF"* ]] && ext=g.vcf.gz
  /usr/bin/time -f "%e s" $G HaplotypeCaller -R ../ref/GRCh38_no_alt_analysis_set.fa -O $O/$c.$ext $args \
      --native-pair-hmm-threads 2 > $O/$c.log 2>&1
  echo "  $c: $(cd .. && $HARNESS/hg002/cmp1.sh $O/$c.$ext group4/truth $ext) ($(tail -1 $O/$c.log))"
done
$G HaplotypeCaller -R ../ref/GRCh38_no_alt_analysis_set.fa -I in/L1.bam -I in/L1.bam -L chr8:1-15000000 -O $O/dup.vcf.gz > $O/dup.log 2>&1
echo "  one BAM twice: rc=$? $(grep -m1 '^variant-rs:' $O/dup.log)"
