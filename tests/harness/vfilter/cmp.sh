#!/usr/bin/env bash
# variant-rs VariantFiltration / MakeSitesOnlyVcf with each truth case's arguments (truth.sh), in truth/cases.txt
# order, vs GATK: VCFs (decompressed) on every line but the command-line header line; where GATK failed, variant-rs
# must fail too.   usage: cmp.sh [OUTDIR] [case...]
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/vfilter" || exit 1; D=$PWD

G=${VARIANT_RS:-$REPO/target/release/variant-rs}
O=${1:-$D/rust}; shift; mkdir -p $O; O=$(realpath $O)
txt() { case $1 in *.gz) zcat $1;; *) cat $1;; esac; }
for c in ${@:-$(cat truth/cases.txt)}; do
  args=$($HARNESS/cmdargs.sh truth/$c.cmd)
  t=$(date +%s); eval "$G ${args//@O@/$O}" > $O/$c.log 2>&1; rc=$?; s=$(( $(date +%s) - t ))
  trc=$(grep -o '^rc=[0-9]*' truth/$c.log | cut -d= -f2); tg=$(grep -o '^rc=[0-9]* [0-9]*s' truth/$c.log | sed 's/.* //')
  out=$(echo "$args" | grep -oE '(-O|--OUTPUT) @O@/[^ ]*' | sed 's|.*@O@/||')
  if [ "$trc" != 0 ]; then
    [ $rc != 0 ] && r="both fail (GATK rc=$trc, variant-rs rc=$rc)" || r="DIFFERENT: GATK failed (rc=$trc), variant-rs did not"
  elif [ $rc != 0 ]; then
    r="DIFFERENT: variant-rs failed (rc=$rc, $O/$c.log)"
  else
    n=$(txt truth/$out | grep -vc '^#')
    if cmp -s <(txt truth/$out | grep -v CommandLine) <(txt $O/$out | grep -v CommandLine); then r="identical ($n records)"
    else r="DIFFERENT ($(diff <(txt truth/$out | grep -v CommandLine) <(txt $O/$out | grep -v CommandLine) | grep -c '^[<>]') lines differ, $n records)"; fi
  fi
  echo "  $c: $r; ${s}s (GATK $tg)"
done
