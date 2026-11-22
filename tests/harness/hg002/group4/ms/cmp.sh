#!/usr/bin/env bash
# variant-rs multi-sample HaplotypeCaller with each truth case's GATK arguments (truth/<case>.cmd, see make.sh) vs GATK;
# then the refusal of -ERC GVCF with several samples.  usage: cmp.sh [OUTDIR] [case...]
HARNESS=$(cd "$(dirname "$0")/../../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/hg002/group4/ms" || exit 1; D=$PWD

G=${VARIANT_RS:-$REPO/target/release/variant-rs}
REF=$VARCALL_DATA/bqsr/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
O=${1:-$D/rust}; shift; mkdir -p $O; O=$(realpath $O)
cases=${@:-$(ls truth/*.cmd | xargs -n1 basename | sed 's/\.cmd$//' | sort)}
for c in $cases; do
  t0=$(date +%s.%N)
  $G HaplotypeCaller -R $REF -O $O/$c.vcf.gz $($HARNESS/cmdargs.sh truth/$c.cmd) > $O/$c.log 2>&1
  echo "  $c: $($HARNESS/hg002/cmp1.sh $O/$c.vcf.gz truth vcf.gz) ($(echo "$(date +%s.%N) - $t0" | bc | xargs printf '%.1f') s)"
done
$G HaplotypeCaller -R $REF -I in/HG002.bam -I in/NA12878.bam -L chr8:1-1000000 -ERC GVCF -O $O/gvcf.g.vcf.gz > $O/gvcf.log 2>&1
echo "  -ERC GVCF, two samples: rc=$? $(grep -m1 '^variant-rs:' $O/gvcf.log)"
