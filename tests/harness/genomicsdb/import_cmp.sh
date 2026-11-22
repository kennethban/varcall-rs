#!/usr/bin/env bash
# variant-rs GenomicsDBImport with each GATK workspace's import arguments (read from its vcfheader.vcf) vs that
# workspace (truth.sh, truth2.sh, import_truth.sh): every cell (gdbdump), callset.json and vidmap.json byte for byte,
# vcfheader.vcf but the command line, the arrays, the meta files' content, the number of fragments. Then
# GenotypeGVCFs, GATK's and variant-rs's, on two of the new workspaces vs the GATK truth.
#   usage: import_cmp.sh [OUTDIR] [workspace...]
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/genomicsdb" || exit 1
X=${VARIANT_RS:-$REPO/target/release/variant-rs}
O=${1:-$PWD/rust_import}; shift; mkdir -p $O; O=$(realpath $O)
# the import arguments GATK's command line shows (explicit ones first; the defaults are dropped here)
args() {
  grep -o 'CommandLine="GenomicsDBImport [^"]*' $1/vcfheader.vcf | sed 's/CommandLine="GenomicsDBImport //' \
    | grep -oE -- '--(variant|sample-name-map|intervals|batch-size|consolidate|merge-input-intervals|bypass-feature-reader|reader-threads|genomicsdb-shared-posixfs-optimizations) [^ ]+' | tr '\n' ' '
}
metas() { for d in $1/*/; do cat $d/genomicsdb_meta_dir/genomicsdb_column_bounds.json; for f in $d/genomicsdb_meta_dir/genomicsdb_meta_*.json; do cat $f; echo; done; done | sort; }
for w in ${@:-ws_rb ws_raw ws_rb_b1 ws11 ws_wg ws_mid ws_rb_bypass ws_raw_bypass ws_sarek ws_sarek_wes}; do
  rm -rf $O/$w
  t=$(date +%s); $X GenomicsDBImport --genomicsdb-workspace-path $O/$w $(args $w) > $O/$w.log 2>&1; rc=$?; s=$(( $(date +%s) - t ))
  [ $rc = 0 ] || { echo "  $w: variant-rs FAILED (rc=$rc, $O/$w.log)"; continue; }
  r="cells $([ "$($X gdbdump $w | md5sum)" = "$($X gdbdump $O/$w | md5sum)" ] && echo identical || echo DIFFERENT)"
  for f in callset.json vidmap.json; do r+=", $f $(cmp -s $w/$f $O/$w/$f && echo identical || echo DIFFERENT)"; done
  r+=", vcfheader.vcf $(cmp -s <(grep -v CommandLine $w/vcfheader.vcf) <(grep -v CommandLine $O/$w/vcfheader.vcf) && echo identical || echo DIFFERENT)"
  r+=", arrays $(diff <(cd $w; ls -d */) <(cd $O/$w; ls -d */) > /dev/null && echo same || echo DIFFERENT)"
  r+=", meta files $(diff <(metas $w) <(metas $O/$w) > /dev/null && echo same || echo DIFFERENT)"
  r+=", fragments $([ "$(ls -d $w/*/__*/ | wc -l)" = "$(ls -d $O/$w/*/__*/ | wc -l)" ] && echo same || echo DIFFERENT)"
  echo "  $w: $r ($(grep -o 'Elapsed time: [0-9.]* minutes' $w.import.log | sed 's/Elapsed time: //') GATK, ${s}s variant-rs)"
done
# GATK and variant-rs GenotypeGVCFs on the new workspaces
J=${GATK_JAR:-}
gg() { # workspace reference truth
  [ -d $O/$1 ] || return
  local a=$($HARNESS/cmdargs.sh truth/$1.cmd | sed "s#gendb://$1#gendb://$O/$1#")
  local n=$(zcat truth/$1.vcf.gz | grep -vc '^#')
  $X GenotypeGVCFs -R $2 $a -O $O/$1.rs.vcf.gz > $O/$1.rs.log 2>&1
  echo "  GenotypeGVCFs variant-rs on $1: $(cmp -s <(zcat truth/$1.vcf.gz | grep -v CommandLine) <(zcat $O/$1.rs.vcf.gz | grep -v CommandLine) && echo identical || echo DIFFERENT) ($n records)"
  if [ -n "$J" ]; then
    java -Xmx8g -jar $J GenotypeGVCFs -R $2 $a -O $O/$1.gatk.vcf.gz > $O/$1.gatk.log 2>&1
    echo "  GenotypeGVCFs GATK on $1: $(cmp -s <(zcat truth/$1.vcf.gz | grep -v CommandLine) <(zcat $O/$1.gatk.vcf.gz | grep -v CommandLine) && echo identical || echo DIFFERENT) ($n records)"
  fi
}
gg ws_wg $VARCALL_DATA/hg002/ref/GRCh38_no_alt_analysis_set.fa
gg ws_sarek $VARCALL_DATA/bqsr/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
