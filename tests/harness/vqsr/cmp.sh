#!/usr/bin/env bash
# variant-rs VariantRecalibrator / ApplyVQSR with each truth case's arguments (truth.sh), in truth/cases.txt order, vs
# GATK: tranches files byte for byte, VCFs (decompressed) on every line but the command-line header line.
#   usage: cmp.sh [OUTDIR] [case...]
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/vqsr" || exit 1; D=$PWD

G=${VARIANT_RS:-$REPO/target/release/variant-rs}
O=${1:-$D/rust}; shift; mkdir -p $O; O=$(realpath $O)
same() { [ -s $1 ] && [ -s $2 ] && cmp -s <(zcat $1 | grep -v CommandLine) <(zcat $2 | grep -v CommandLine); }
n() { [ -s $1 ] && zcat $1 | grep -vc '^#' || echo 0; }
for c in ${@:-$(cat truth/cases.txt)}; do
  args=$($HARNESS/cmdargs.sh truth/$c.cmd)
  t=$(date +%s); $G ${args//@O@/$O} > $O/$c.log 2>&1; rc=$?; s=$(( $(date +%s) - t ))
  tg=$(grep -o 'rc=[0-9]* [0-9]*s' truth/$c.log | sed 's/.* //')
  grep -q 'rc=0 ' truth/$c.log || { echo "  $c: GATK FAILED (truth/$c.log)"; continue; }
  r=""
  if [ -f truth/$c.tranches ]; then
    r="tranches $(cmp -s truth/$c.tranches $O/$c.tranches && echo identical || echo DIFFERENT), "
    r+="recal $(same truth/$c.recal.vcf.gz $O/$c.recal.vcf.gz && echo identical || echo DIFFERENT) ($(n truth/$c.recal.vcf.gz) records)"
  else
    r="vcf $(same truth/$c.vcf.gz $O/$c.vcf.gz && echo identical || echo DIFFERENT) ($(n truth/$c.vcf.gz) records)"
    [ -f truth/$c.vcf.gz.tbi ] && { [ -f $O/$c.vcf.gz.tbi ] && r+=", index written" || r+=", index MISSING"; }
  fi
  echo "  $c: rc=$rc $r; ${s}s (GATK $tg)"
done
