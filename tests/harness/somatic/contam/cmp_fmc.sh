#!/usr/bin/env bash
# variant-rs filtermutectcalls with contamination / segmentation tables vs GATK (truth.sh and the fc_deep_* runs).
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/contam" || exit 1
A=${VARIANT_RS:-$REPO/target/release/variant-rs}
mkdir -p rust
chk() { # name input.vcf args...
  local n=$1 in=$2; shift 2
  $A filtermutectcalls -V $in -O rust/$n.vcf "$@" 2> rust/$n.log || tail -3 rust/$n.log
  v=$(diff <(grep -v CommandLine truth/$n.vcf) <(grep -v CommandLine rust/$n.vcf) >/dev/null && echo identical || echo DIFFER)
  s=$(cmp -s rust/$n.vcf.filteringStats.tsv truth/$n.vcf.filteringStats.tsv && echo identical || echo DIFFER)
  echo "$n: VCF $v ($(grep -v '^#' rust/$n.vcf | awk '$7=="PASS"' | wc -l) PASS), filtering stats $s"
  [ "$v" = identical ] || python3 $HARNESS/somatic/tools/cmp_vcf.py truth/$n.vcf rust/$n.vcf 10
  [ "$s" = identical ] || diff truth/$n.vcf.filteringStats.tsv rust/$n.vcf.filteringStats.tsv | head -10
}
for i in 1 2 3 4; do chk fc_$i ../dream/truth/m2_$i.vcf --contamination-table truth/dream_$i.contamination.table --tumor-segmentation truth/dream_$i.segments.table; done
chk fc_deep ../dsmerge/gatk_ds.vcf --contamination-table truth/deep.contamination.table --tumor-segmentation truth/deep.segments.table
chk fc_deep_wgs ../dsmerge/gatk_ds.vcf --contamination-table truth/wgs_to.contamination.table --tumor-segmentation truth/wgs_to.segments.table
chk fc_deep_est ../dsmerge/gatk_ds.vcf --contamination-estimate 0.05
