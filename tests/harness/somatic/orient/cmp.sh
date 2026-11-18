#!/usr/bin/env bash
# variant-rs mutect2 --f1r2-tar-gz + learnreadorientationmodel (+ filtermutectcalls --ob-priors) vs GATK (truth.sh).
# usage: cmp.sh [1 2 3 4 deep]
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/orient" || exit 1
A=${VARIANT_RS:-$REPO/target/release/variant-rs}
D=../dream/data; P=$VARCALL_DATA
mkdir -p rust
for x in ${@:-1 2 3 4 deep}; do
  N=$(grep -m1 '^##normal_sample=' truth/m2_$x.vcf | cut -d= -f2)
  if [ $x = deep ]; then
    R=../dsmerge/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
    /usr/bin/time -f "%e s" $A mutect2 -R $R -I $P/markdup/input/na12878.chr20.bam -I $P/hg002_pcrplus/bam/HG002.pcrplus.bam -normal "$N" \
      -L ../dsmerge/ds.intervals --f1r2-tar-gz rust/f1r2_$x.tar.gz -O rust/m2_$x.vcf 2> rust/m2_$x.log
  else
    R=$D/human_g1k_v37.20.21.fasta
    /usr/bin/time -f "%e s" $A mutect2 -R $R -I $D/tumor_$x.bam -I $D/normal_$x.bam -normal "$N" -L 20 \
      --germline-resource $D/very-small-gnomad.vcf --f1r2-tar-gz rust/f1r2_$x.tar.gz -O rust/m2_$x.vcf 2> rust/m2_$x.log
  fi
  v=$(diff <(grep -v CommandLine truth/m2_$x.vcf) <(grep -v CommandLine rust/m2_$x.vcf) >/dev/null && echo identical || echo DIFFER)
  echo "== $x: M2 VCF $v, $(tail -1 rust/m2_$x.log)"
  rm -rf rust/t_$x rust/r_$x; mkdir -p rust/t_$x rust/r_$x
  tar -xzf truth/f1r2_$x.tar.gz -C rust/t_$x; tar -xzf rust/f1r2_$x.tar.gz -C rust/r_$x
  for f in rust/t_$x/*.alt_table; do b=$(basename "$f"); cmp -s "$f" "rust/r_$x/$b" && echo "  alt_table identical ($(($(wc -l < "$f") - 2)) sites)" || { echo "  alt_table DIFFERS"; diff "$f" "rust/r_$x/$b" | head -6; }; done
  for e in ref_histogram alt_histogram; do f=$(ls rust/t_$x/*.$e); echo "  $e: $(python3 $HARNESS/somatic/orient/cmp_hist.py "$f" "rust/r_$x/$(basename "$f")")"; done
  # LearnReadOrientationModel on GATK's counts, and on ours
  $A learnreadorientationmodel -I truth/f1r2_$x.tar.gz -O rust/priors_g_$x.tar.gz 2> rust/lrom_g_$x.log
  $A learnreadorientationmodel -I rust/f1r2_$x.tar.gz -O rust/priors_$x.tar.gz 2> rust/lrom_$x.log
  rm -rf rust/pt_$x rust/pg_$x rust/pr_$x; mkdir -p rust/pt_$x rust/pg_$x rust/pr_$x
  tar -xzf truth/priors_$x.tar.gz -C rust/pt_$x; tar -xzf rust/priors_g_$x.tar.gz -C rust/pg_$x; tar -xzf rust/priors_$x.tar.gz -C rust/pr_$x
  f=$(ls rust/pt_$x/*); b=$(basename "$f")
  cmp -s "$f" "rust/pg_$x/$b" && echo "  priors (from GATK counts) identical" || { echo "  priors (from GATK counts) DIFFER"; diff "$f" "rust/pg_$x/$b" | head -6; }
  cmp -s "$f" "rust/pr_$x/$b" && echo "  priors (end to end) identical" || echo "  priors (end to end) DIFFER"
done
