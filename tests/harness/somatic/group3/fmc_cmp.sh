#!/usr/bin/env bash
# variant-rs FilterMutectCalls with each fmc_truth.sh case's GATK arguments vs GATK: VCF (all but the command line)
# and filtering stats. warp2 reads variant-rs's own warp1 output.  usage: fmc_cmp.sh [OUTDIR] [case...]
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/group3" || exit 1; D=$PWD
G=${VARIANT_RS:-$REPO/target/release/variant-rs}
O=${1:-$D/rust_fmc}; shift; mkdir -p $O
cases=${@:-mito to_mito to_opts tn_opts warp1 warp2}
for c in $cases; do
  args=$(sed "s|fmc/warp1.vcf|$O/warp1.vcf|" fmc/$c.cmd)
  /usr/bin/time -f "%e s" $G FilterMutectCalls $args -O $O/$c.vcf > $O/$c.log 2>&1
  rc=$?
  n=$(grep -vc '^#' fmc/$c.vcf); p=$(grep -v '^#' fmc/$c.vcf | awk '$7=="PASS"' | wc -l)
  v=$(diff <(grep -v CommandLine fmc/$c.vcf) <(grep -v CommandLine $O/$c.vcf 2>/dev/null) > $O/$c.diff && echo identical || echo "DIFFERENT($(grep -c '^[<>]' $O/$c.diff))")
  st=$(cmp -s fmc/$c.vcf.filteringStats.tsv $O/$c.vcf.filteringStats.tsv && echo identical || echo DIFFERENT)
  echo "  $c: rc=$rc VCF $v ($n records, $p PASS), filtering stats $st ($(tail -1 $O/$c.log))"
done
