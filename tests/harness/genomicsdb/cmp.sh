#!/usr/bin/env bash
# variant-rs GenotypeGVCFs -V gendb://<ws> with each truth case's arguments (truth/<case>.cmd, truth.sh) vs GATK: data
# lines and header (all but GenotypeGVCFs' own command line) byte for byte; with RECEIVED=1 also the records the
# walker receives from the workspace vs received/<ws>.vcf (the harness dump).  usage: cmp.sh [OUTDIR] [case...]
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/genomicsdb" || exit 1; D=$PWD

G=${VARIANT_RS:-$REPO/target/release/variant-rs}
REF=$VARCALL_DATA/bqsr/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
O=${1:-$D/rust}; shift; mkdir -p $O; O=$(realpath $O)
for c in ${@:-$(ls truth/*.cmd | xargs -n1 basename | sed 's/\.cmd$//' | sort)}; do
  w=${c%.default}
  case $c in ws_wg*|ws_mid*) REF=$VARCALL_DATA/hg002/ref/GRCh38_no_alt_analysis_set.fa;; *) REF=$VARCALL_DATA/bqsr/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa;; esac
  VARIANT_RS_GDB_DUMP=$( [ -n "$RECEIVED" ] && echo $O/$c.received ) $G GenotypeGVCFs -R $REF $($HARNESS/cmdargs.sh truth/$c.cmd) -O $O/$c.vcf.gz > $O/$c.log 2>&1; rc=$?
  t=truth/$c.vcf.gz; o=$O/$c.vcf.gz
  n=$(zcat $t | grep -vc '^#')
  d=$(cmp -s <(zcat $t | grep -v '^#') <(zcat $o 2>/dev/null | grep -v '^#') && echo IDENTICAL || echo "DIFFERENT($(diff <(zcat $t | grep -v '^#') <(zcat $o 2>/dev/null | grep -v '^#') | grep -c '^[<>]'))")
  h=$(diff <(zcat $t | grep '^#' | grep -v '^##GATKCommandLine=<ID=GenotypeGVCFs') <(zcat $o 2>/dev/null | grep '^#' | grep -v '^##variant-rsCommandLine') | grep -c '^[<>]')
  r=""; [ -n "$RECEIVED" ] && [ -f received/$w.vcf ] && r="; received $(cmp -s <(grep -v '^#' received/$w.vcf) $O/$c.received && echo identical || echo DIFFERENT)"
  echo "  $c: rc=$rc $d $n hdr=$h$r ($(grep -o '([0-9.]*s)' $O/$c.log))"
done
