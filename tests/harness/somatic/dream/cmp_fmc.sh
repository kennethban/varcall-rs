#!/usr/bin/env bash
# FilterMutectCalls check: variant-rs filtermutectcalls on GATK's Mutect2 VCFs vs GATK FilterMutectCalls (filter_truth.sh):
# VCF (header without command lines) and filtering stats.
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/dream" || exit 1
B=${VARIANT_RS:-$REPO/target/release/variant-rs}
mkdir -p rust
for i in ${@:-1 2 3 4}; do
  $B filtermutectcalls -R data/human_g1k_v37.20.21.fasta -V truth/m2_$i.vcf -O rust/f_$i.vcf 2> rust/f_$i.log
  v=$(diff <(grep -v CommandLine truth/f_$i.vcf) <(grep -v CommandLine rust/f_$i.vcf) >/dev/null && echo identical || echo DIFFER)
  s=$(cmp -s rust/f_$i.vcf.filteringStats.tsv truth/f_$i.vcf.filteringStats.tsv && echo identical || echo DIFFER)
  echo "pair $i: VCF $v ($(grep -v '^#' rust/f_$i.vcf | awk '$7=="PASS"' | wc -l) PASS / $(grep -vc '^#' rust/f_$i.vcf)), filtering stats $s; $(tail -1 rust/f_$i.log)"
  [ "$v" = identical ] || python3 $HARNESS/somatic/tools/cmp_vcf.py truth/f_$i.vcf rust/f_$i.vcf 10
  [ "$s" = identical ] || diff truth/f_$i.vcf.filteringStats.tsv rust/f_$i.vcf.filteringStats.tsv | head -20
done
