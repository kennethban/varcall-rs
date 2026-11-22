#!/usr/bin/env bash
# GenotypeGVCFs on a GenomicsDB workspace without -L, as nf-core/sarek runs it (gatk4/genotypegvcfs: --variant
# gendb://ws --reference --dbsnp, no intervals): GATK walks every array over its whole range. GATK truth and the
# variant-rs comparison.  usage: truth_nol.sh truth|cmp [OUTDIR]
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/genomicsdb" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}; G=${VARIANT_RS:-$REPO/target/release/variant-rs}
FULL=$VARCALL_DATA/bqsr/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
NOALT=$VARCALL_DATA/hg002/ref/GRCh38_no_alt_analysis_set.fa
DBSNP=$VARCALL_DATA/bqsr/resources/Homo_sapiens_assembly38.dbsnp138.vcf.gz
declare -A C
C[nol_rb]="-R $FULL -V gendb://ws_rb"
C[nol_rb_sarek]="--reference $FULL --variant gendb://ws_rb --dbsnp $DBSNP --tmp-dir ."
C[nol_raw]="-R $FULL -V gendb://ws_raw"
C[nol_11]="-R $FULL -V gendb://ws11"
C[nol_wg]="-R $NOALT -V gendb://ws_wg"
C[nol_mid]="-R $NOALT -V gendb://ws_mid"
mkdir -p truth_nol
case ${1:-cmp} in
truth)
  for c in "${!C[@]}"; do
    echo "${C[$c]}" | $HARNESS/cmdargs.sh -w > truth_nol/$c.cmd
    [ -f truth_nol/$c.vcf.gz.tbi ] || java -Xmx6g -jar $J GenotypeGVCFs ${C[$c]} -O truth_nol/$c.vcf.gz > truth_nol/$c.log 2>&1 &
  done; wait
  for c in "${!C[@]}"; do echo "$c $(zcat truth_nol/$c.vcf.gz | grep -vc '^#') records $(grep -o 'Elapsed time: .*' truth_nol/$c.log)"; done ;;
cmp)
  O=${2:-rust_nol}; mkdir -p $O
  for c in $(printf '%s\n' "${!C[@]}" | sort); do
    t0=$(date +%s.%N)
    $G GenotypeGVCFs ${C[$c]} -O $O/$c.vcf.gz > $O/$c.log 2>&1; rc=$?
    r=$(cmp -s <(zcat truth_nol/$c.vcf.gz | grep -v '^#') <(zcat $O/$c.vcf.gz | grep -v '^#') && echo IDENTICAL || echo DIFFERENT)
    echo "  $c: rc=$rc $r ($(zcat truth_nol/$c.vcf.gz | grep -vc '^#') records, $(echo "$(date +%s.%N) - $t0" | bc | xargs printf '%.1f') s)$(grep -m1 'running GATK' $O/$c.log)"
  done ;;
esac
