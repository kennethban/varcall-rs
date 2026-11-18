#!/usr/bin/env bash
# variant-rs GenotypeGVCFs with each truth case's arguments (truth/<case>.cmd; multi/truth for the multi-sample cases)
# vs GATK: data lines and header (all but GenotypeGVCFs' own command line) byte for byte.  usage: cmp.sh [OUTDIR] [case...]
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/genotypegvcfs" || exit 1; D=$PWD

G=${VARIANT_RS:-$REPO/target/release/variant-rs}
O=${1:-$D/rust}; shift; mkdir -p $O; O=$(realpath $O)
cases=${@:-$(ls truth/*.cmd multi/truth/*.cmd 2>/dev/null | xargs -n1 basename | sed 's/\.cmd$//' | sort)}
for c in $cases; do
  if [ -f truth/$c.cmd ]; then T=truth; W=$D; REF=$VARCALL_DATA/hg002/ref/GRCh38_no_alt_analysis_set.fa
  else T=multi/truth; W=$D/multi; REF=$VARCALL_DATA/bqsr/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa; fi
  (cd $W && $G GenotypeGVCFs -R $REF $(cat $D/$T/$c.cmd) -O $O/$c.vcf.gz > $O/$c.log 2>&1); rc=$?
  t=$T/$c.vcf.gz; o=$O/$c.vcf.gz
  n=$(zcat $t | grep -vc '^#')
  d=$(cmp -s <(zcat $t | grep -v '^#') <(zcat $o 2>/dev/null | grep -v '^#') && echo IDENTICAL || echo "DIFFERENT($(diff <(zcat $t | grep -v '^#') <(zcat $o 2>/dev/null | grep -v '^#') | grep -c '^[<>]'))")
  h=$(diff <(zcat $t | grep '^#' | grep -v '^##GATKCommandLine=<ID=GenotypeGVCFs') <(zcat $o 2>/dev/null | grep '^#' | grep -v '^##variant-rsCommandLine') | grep -c '^[<>]')
  echo "  $c: rc=$rc $d $n hdr=$h ($(grep -o '([0-9.]*s)' $O/$c.log))"
done
