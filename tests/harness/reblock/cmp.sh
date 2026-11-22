#!/usr/bin/env bash
# variant-rs ReblockGVCF with each truth case's arguments (truth/<case>.cmd) vs GATK: data lines and header (all but
# ReblockGVCF's own command line) byte for byte.  usage: cmp.sh [OUTDIR] [case...]
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/reblock" || exit 1; D=$PWD

G=${VARIANT_RS:-$REPO/target/release/variant-rs}
REF=$VARCALL_DATA/hg002/ref/GRCh38_no_alt_analysis_set.fa
O=${1:-$D/rust}; shift; mkdir -p $O; O=$(realpath $O)
for c in ${@:-$(ls truth/*.cmd | xargs -n1 basename | sed 's/\.cmd$//' | sort)}; do
  $G ReblockGVCF -R $REF $($HARNESS/cmdargs.sh truth/$c.cmd) -O $O/$c.rb.g.vcf.gz > $O/$c.log 2>&1; rc=$?
  t=truth/$c.rb.g.vcf.gz; o=$O/$c.rb.g.vcf.gz
  n=$(zcat $t | grep -vc '^#')
  d=$(cmp -s <(zcat $t | grep -v '^#') <(zcat $o 2>/dev/null | grep -v '^#') && echo IDENTICAL || echo "DIFFERENT($(diff <(zcat $t | grep -v '^#') <(zcat $o 2>/dev/null | grep -v '^#') | grep -c '^[<>]'))")
  h=$(diff <(zcat $t | grep '^#' | grep -v '^##GATKCommandLine=<ID=ReblockGVCF') <(zcat $o 2>/dev/null | grep '^#' | grep -v '^##variant-rsCommandLine') | grep -c '^[<>]')
  echo "  $c: rc=$rc $d $n hdr=$h ($(grep -o '([0-9.]*s)' $O/$c.log))"
done
