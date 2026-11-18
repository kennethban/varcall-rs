#!/usr/bin/env bash
# variant-rs vs GATK truth (truth.sh): records in order, incl. the downsampled deep chrM and QD jitter on later contigs
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/hcorder" || exit 1

X=${VARIANT_RS:-$REPO/target/release/variant-rs}; O=${1:-out}; mkdir -p $O
REF=../hg002/ref/GRCh38_no_alt_analysis_set.fa
for m in vcf gvcf; do
  e=""; [ $m = gvcf ] && e="-ERC GVCF"
  $X HaplotypeCaller -R $REF -I reordered.bam -O $O/out.$m.gz $e > $O/$m.log 2>&1; rc=$?
  a=$(zcat truth/out.$m.gz | grep -v '^##' | md5sum | cut -c1-32); b=$(zcat $O/out.$m.gz 2>/dev/null | grep -v '^##' | md5sum | cut -c1-32)
  n=$(zcat truth/out.$m.gz | grep -vc '^#')
  echo "  $m: rc=$rc records $( [ $a = $b ] && echo "identical ($n)" || echo "DIFFERENT (first contigs GATK: $(zcat truth/out.$m.gz | grep -v '^#' | cut -f1 | uniq | head -3 | tr '\n' ' ') variant-rs: $(zcat $O/out.$m.gz | grep -v '^#' | cut -f1 | uniq | head -3 | tr '\n' ' '))")"
done
