#!/usr/bin/env bash
# The panel-of-normals chain in variant-rs vs GATK (truth.sh): Mutect2 on the four normals (VCF data lines), then
# GenomicsDBImport on GATK's normal VCFs (every cell, vidmap.json, callset.json; one and three arrays), the records a
# walker receives from the workspace (alleles and every sample's AD, vs GATK SelectVariants), and
# CreateSomaticPanelOfNormals on GATK's workspaces with each truth case's arguments.   usage: cmp.sh [OUTDIR]
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/pon" || exit 1
X=${VARIANT_RS:-$REPO/target/release/variant-rs}
O=${1:-$PWD/rust}; mkdir -p $O; O=$(realpath $O)
D=../dream/data; REF=$D/human_g1k_v37.20.21.fasta
for i in 1 2 3 4; do
  $X Mutect2 -R $REF -I $D/normal_$i.bam -L 20 --max-mnp-distance 0 -O $O/normal_$i.vcf.gz > $O/normal_$i.log 2>&1
  echo "  Mutect2 normal_$i: $(cmp -s <(zcat in/normal_$i.vcf.gz | grep -v '^#') <(zcat $O/normal_$i.vcf.gz | grep -v '^#') && echo identical || echo DIFFERENT) ($(zcat in/normal_$i.vcf.gz | grep -vc '^#') records)"
done
V="-V in/normal_1.vcf.gz -V in/normal_2.vcf.gz -V in/normal_3.vcf.gz -V in/normal_4.vcf.gz"
# alleles and each sample's AD of the received records ('.' for a sample without AD)
ad() { awk 'BEGIN{FS=OFS="\t"} !/^#/ {n=split($9,f,":"); k=0; for(i=1;i<=n;i++) if(f[i]=="AD") k=i; s=$1 OFS $2 OFS $4 OFS $5; for(j=10;j<=NF;j++){split($j,v,":"); x=k?v[k]:"."; s=s OFS (x==""?".":x)} print s}'; }
for c in "ws_pon:20" "ws_pon3:in/three.intervals"; do
  w=${c%%:*}; l=${c#*:}
  rm -rf $O/$w
  $X GenomicsDBImport --genomicsdb-workspace-path $O/$w -R $REF $V -L $l > $O/$w.log 2>&1
  r="cells $([ "$($X gdbdump $w | md5sum)" = "$($X gdbdump $O/$w | md5sum)" ] && echo identical || echo DIFFERENT)"
  for f in callset.json vidmap.json; do r+=", $f $(cmp -s $w/$f $O/$w/$f && echo identical || echo DIFFERENT)"; done
  echo "  GenomicsDBImport $w: $r"
  echo "  received records $w: $(cmp -s <(ad < truth/$w.received.vcf) <($X gdbreceive $w $REF | ad) && echo identical || echo DIFFERENT) ($(grep -vc '^#' truth/$w.received.vcf) records, alleles and AD)"
done
for c in $(ls truth/*.cmd | xargs -n1 basename | sed 's/\.cmd$//' | grep -v received); do
  $X CreateSomaticPanelOfNormals -R $REF $($HARNESS/cmdargs.sh truth/$c.cmd) -O $O/$c.vcf > $O/$c.log 2>&1; rc=$?
  echo "  CreateSomaticPanelOfNormals $c: rc=$rc $(cmp -s <(grep -v CommandLine truth/$c.vcf) <(grep -v CommandLine $O/$c.vcf) && echo identical || echo DIFFERENT) ($(grep -vc '^#' truth/$c.vcf) records)"
done
