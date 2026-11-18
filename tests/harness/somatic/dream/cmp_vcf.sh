#!/usr/bin/env bash
# S3 check: variant-rs mutect2 unfiltered VCF + .stats vs GATK (truth.sh) on the DREAM pairs (header compared without the command line).
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/dream" || exit 1

B=${VARIANT_RS:-$REPO/target/release/variant-rs}
for i in ${@:-1 2 3 4}; do
  N=$(grep -m1 '^##normal_sample=' truth/m2_$i.vcf | cut -d= -f2)
  /usr/bin/time -f "%e s" $B mutect2 -R data/human_g1k_v37.20.21.fasta -I data/tumor_$i.bam -I data/normal_$i.bam -normal "$N" -L 20 \
    --germline-resource data/very-small-gnomad.vcf -O rust/m2_$i.vcf 2> rust/m2_$i.log
  v=$(diff <(grep -v CommandLine truth/m2_$i.vcf) <(grep -v CommandLine rust/m2_$i.vcf) >/dev/null && echo identical || echo DIFFER)
  s=$(cmp -s rust/m2_$i.vcf.stats truth/m2_$i.vcf.stats && echo identical || echo DIFFER)
  echo "pair $i: VCF $v ($(grep -vc '^#' rust/m2_$i.vcf) records), stats $s, $(tail -1 rust/m2_$i.log)"
  [ "$v" = identical ] || python3 $HARNESS/somatic/tools/cmp_vcf.py truth/m2_$i.vcf rust/m2_$i.vcf 10
done
