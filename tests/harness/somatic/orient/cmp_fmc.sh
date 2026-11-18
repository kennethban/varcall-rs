#!/usr/bin/env bash
# variant-rs filtermutectcalls --ob-priors vs GATK (truth.sh): on GATK's VCF + priors, and end to end (variant-rs M2 + priors; cmp.sh first).
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/orient" || exit 1
A=${VARIANT_RS:-$REPO/target/release/variant-rs}
for x in ${@:-1 2 3 4 deep}; do
  if [ $x = deep ]; then R=../dsmerge/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa; else R=../dream/data/human_g1k_v37.20.21.fasta; fi
  for mode in g e; do
    if [ $mode = g ]; then V=truth/m2_$x.vcf; O=truth/priors_$x.tar.gz; else V=rust/m2_$x.vcf; O=rust/priors_$x.tar.gz; fi
    n=fo_${mode}_$x
    $A filtermutectcalls -R $R -V $V --ob-priors $O -O rust/$n.vcf 2> rust/$n.log || tail -3 rust/$n.log
    v=$(diff <(grep -v CommandLine truth/fo_$x.vcf) <(grep -v CommandLine rust/$n.vcf) >/dev/null && echo identical || echo DIFFER)
    s=$(cmp -s rust/$n.vcf.filteringStats.tsv truth/fo_$x.vcf.filteringStats.tsv && echo identical || echo DIFFER)
    echo "$x ($([ $mode = g ] && echo "GATK M2 + priors" || echo "end to end")): VCF $v ($(grep -v '^#' rust/$n.vcf | awk '$7=="PASS"' | wc -l) PASS, $(grep -v '^#' rust/$n.vcf | grep -c orientation) orientation), filtering stats $s"
    [ "$v" = identical ] || python3 $HARNESS/somatic/tools/cmp_vcf.py truth/fo_$x.vcf rust/$n.vcf 10
    [ "$s" = identical ] || diff truth/fo_$x.vcf.filteringStats.tsv rust/$n.vcf.filteringStats.tsv | head -10
  done
done
