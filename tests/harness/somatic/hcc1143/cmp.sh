#!/usr/bin/env bash
# Compares the variant-rs chain (rust/) with GATK (truth/), both from run.sh.
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/hcc1143" || exit 1
same() { if cmp -s "$1" "$2"; then echo "  $3 identical"; else echo "  $3 DIFFERS"; diff "$1" "$2" | head -6; fi; }
v() { diff <(grep -v CommandLine $1) <(grep -v CommandLine $2) > /dev/null && echo "  $3 identical ($(grep -vc '^#' $2) records)" || { echo "  $3 DIFFERS"; python3 $HARNESS/somatic/tools/cmp_vcf.py $1 $2 10; }; }
v truth/m2.vcf rust/m2.vcf "Mutect2 VCF"
same truth/m2.vcf.stats rust/m2.vcf.stats "Mutect2 stats"
rm -rf rust/xt rust/xr; mkdir -p rust/xt rust/xr; tar -xzf truth/f1r2.tar.gz -C rust/xt; tar -xzf rust/f1r2.tar.gz -C rust/xr
f=$(ls rust/xt/*.alt_table); same "$f" rust/xr/$(basename "$f") "F1R2 alt table ($(($(wc -l < "$f") - 2)) sites)"
for e in ref_histogram alt_histogram; do f=$(ls rust/xt/*.$e); echo "  F1R2 $e: $(python3 $HARNESS/somatic/orient/cmp_hist.py "$f" "rust/xr/$(basename "$f")")"; done
same truth/tumor.pileups.table rust/tumor.pileups.table "tumor pileups ($(grep -vc '^#\|^contig' truth/tumor.pileups.table) sites)"
same truth/normal.pileups.table rust/normal.pileups.table "normal pileups"
same truth/contamination.table rust/contamination.table "contamination ($(tail -1 truth/contamination.table | cut -f2,3))"
same truth/segments.table rust/segments.table "segments ($(grep -vc '^#\|^contig' truth/segments.table))"
rm -rf rust/pt rust/pr; mkdir -p rust/pt rust/pr; tar -xzf truth/priors.tar.gz -C rust/pt; tar -xzf rust/priors.tar.gz -C rust/pr
f=$(ls rust/pt/*); same "$f" rust/pr/$(basename "$f") "orientation priors"
v truth/filtered.vcf rust/filtered.vcf "FilterMutectCalls VCF"
same truth/filtered.vcf.filteringStats.tsv rust/filtered.vcf.filteringStats.tsv "filtering stats"
