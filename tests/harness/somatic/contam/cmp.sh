#!/usr/bin/env bash
# variant-rs getpileupsummaries + calculatecontamination vs the GATK truth tables (truth.sh). usage: cmp.sh [dream] [deep] [wgs]
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/contam" || exit 1
A=${VARIANT_RS:-$REPO/target/release/variant-rs}
D=../dream/data; P=$VARCALL_DATA
mkdir -p rust
same() { if cmp -s $1 $2; then echo "  $(basename $1) identical"; else echo "  $(basename $1) DIFFERS"; diff $1 $2 | head -6; fi; }
gps() { $A getpileupsummaries -I $2 -V $3 -L $3 -O rust/$1.pileups.table 2> rust/$1.gps.log || tail -3 rust/$1.gps.log; same truth/$1.pileups.table rust/$1.pileups.table; }
cc() { local m=""; [ -n "$3" ] && m="-matched truth/$3.pileups.table"
  $A calculatecontamination -I truth/$2.pileups.table $m -O rust/$1.contamination.table --tumor-segmentation rust/$1.segments.table 2> rust/$1.cc.log || tail -3 rust/$1.cc.log
  same truth/$1.contamination.table rust/$1.contamination.table; same truth/$1.segments.table rust/$1.segments.table; }
for what in ${@:-dream deep wgs}; do case $what in
dream) for i in 1 2 3 4; do echo "DREAM $i"
  gps tumor_$i $D/tumor_$i.bam res/exac_20.b37.vcf.gz; gps normal_$i $D/normal_$i.bam res/exac_20.b37.vcf.gz
  cc dream_$i tumor_$i normal_$i; cc dream_${i}_to tumor_$i; done ;;
deep) echo "deep chr20"
  gps na12878_chr20 $P/markdup/input/na12878.chr20.bam res/exac_chr20.hg38.vcf.gz
  gps hg002_chr20 $P/hg002_pcrplus/bam/HG002.pcrplus.bam res/exac_chr20.hg38.vcf.gz
  cc deep na12878_chr20 hg002_chr20; cc deep_to na12878_chr20 ;;
wgs) echo "WGS"
  gps na12878_wgs $P/wgs/bam/NA12878.bam res/small_exac_common_3.hg38.vcf.gz
  gps hg002_wgs $P/hg002_pcrplus/bam/HG002.pcrplus.bam res/small_exac_common_3.hg38.vcf.gz
  cc wgs na12878_wgs hg002_wgs; cc wgs_to na12878_wgs; cc wgs_hg002_to hg002_wgs ;;
esac; done
