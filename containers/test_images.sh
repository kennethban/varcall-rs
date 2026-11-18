#!/usr/bin/env bash
# Checks the two images:
#   1. tests/selftest.sh with align-rs and variant-rs run in the images (GPU and CPU paths when the host has a GPU):
#      align-rs vs bwa-mem2, variant-rs vs the GATK truth in tests/truth
#   2. the variant-rs image as sarek runs it: through the image's entrypoint, `gatk --version` as nf-core's modules
#      read it, the libraries, and a command handed to the image's GATK
#   3. with VARCALL_DATA (tests/harness/README.md): sarek's command forms on real data vs GATK truth
#
#   containers/test_images.sh [OUTDIR]       (default target/test_images; ALIGN_RS_IMAGE / VARIANT_RS_IMAGE: the tags)
# Needs what tests/selftest.sh needs on the host (python3, bwa-mem2, samtools, bgzip, tabix). Exit status: the number
# of differences.
set -u
REPO=$(cd "$(dirname "$0")/.." && pwd); O=$(mkdir -p "${1:-$REPO/target/test_images}" && cd "${1:-$REPO/target/test_images}" && pwd)
AI=${ALIGN_RS_IMAGE:-align-rs:0.1.0}; VI=${VARIANT_RS_IMAGE:-variant-rs:0.1.0-gatk4.7.0.0}
GPU=; nvidia-smi -L > /dev/null 2>&1 && GPU="--gpus all"
bad=0; say() { printf '  %-44s %s\n' "$1" "$2"; [[ "$2" == identical* || "$2" == ok* ]] || bad=$((bad + 1)); }
# docker run as the caller, with OUTDIR (and the data) at the same paths, in the caller's directory
d() { docker run --rm -u $(id -u):$(id -g) -e HOME=/tmp -v $O:$O ${VARCALL_DATA:+-v $VARCALL_DATA:$VARCALL_DATA} -w "$PWD" "$@"; }
mkdir -p $O/bin
printf '#!/bin/sh\nexec docker run --rm %s -u %s -e HOME=/tmp -v %s:%s -w "$PWD" %s /opt/align-rs/bin/align-rs "$@"\n' "$GPU" $(id -u):$(id -g) $O $O $AI > $O/bin/align-rs
printf '#!/bin/sh\nexec docker run --rm %s -u %s -e HOME=/tmp -v %s:%s -w "$PWD" %s /opt/variant-rs/bin/variant-rs "$@"\n' "$GPU" $(id -u):$(id -g) $O $O $VI > $O/bin/variant-rs
chmod +x $O/bin/*

echo "#### 1. tests/selftest.sh in the images ($AI, $VI${GPU:+, GPU})"
ALIGN_RS=$O/bin/align-rs VARIANT_RS=$O/bin/variant-rs SELFTEST_PAIRHMM="cpu${GPU:+ cuda}" $REPO/tests/selftest.sh $O/selftest; bad=$?
if [ -n "$GPU" ]; then
  say "align-rs in the image used the GPU" "$(grep -q 'GPU: extension' $O/selftest/align/align.log && echo ok || echo "NO ($O/selftest/align/align.log)")"
  say "variant-rs in the image used the GPU" "$(grep -q 'PairHMM CUDA' $O/selftest/variant_cuda/hc.log && echo ok || echo "NO ($O/selftest/variant_cuda/hc.log)")"
fi

echo "#### 2. the variant-rs image as sarek runs it"
say "libraries" "$(n=$(d $VI sh -c 'ldd /opt/variant-rs/bin/variant-rs | grep -c "not found"'); [ "$n" = 0 ] && echo ok || echo "$n not found")"
say "entrypoint: gatk is the wrapper" "$(docker run --rm $VI sh -c 'VARIANT_RS_DRY_RUN=1 gatk MergeVcfs -I a.vcf.gz -O b.vcf.gz' | grep -q '^RUN variant-rs' && echo ok || echo NO)"
say "gatk --version, nf-core's sed" "$(v=$(docker run --rm $VI sh -c "gatk --version | sed -n '/GATK.*v/s/.*v//p'"); [ "$v" = 4.7.0.0 ] && echo "ok ($v)" || echo "WRONG ($v)")"
cd $O/selftest/data
d $VI gatk IntervalListToBed --INPUT /dev/null --OUTPUT $O/x.bed > $O/handover.log 2>&1
say "IntervalListToBed handed to GATK" "$(grep -q 'running GATK instead' $O/handover.log && echo ok || echo "NO ($O/handover.log)")"

if [ -n "${VARCALL_DATA:-}" ]; then
  echo "#### 3. sarek's command forms on real data (VARCALL_DATA=$VARCALL_DATA)"
  P=$VARCALL_DATA; R=$O/real; mkdir -p $R; cd $R
  eng() { grep -q "Runtime.totalMemory\|Executing as" $1 && echo "GATK (Java) ran it" || echo "variant-rs ran it"; }
  H=$P/hg002; HC="gatk --java-options -Xmx4g HaplotypeCaller -R $H/ref/GRCh38_no_alt_analysis_set.fa -I $H/bam/HG002.bam -L chr8:1-15000000 --native-pair-hmm-threads 4 --tmp-dir ."
  for x in ${GPU:+gpu} cpu; do
    d $( [ $x = gpu ] && echo $GPU ) $VI $HC -O $R/hc_$x.vcf.gz > $R/hc_$x.log 2>&1
    say "HaplotypeCaller HG002 chr8:1-15M ($x)" "$(cmp -s <(zcat $P/containers/variant-rs/truth/hc_hg002_chr8_1-15M.gatk4.7.0.0.vcf.gz | grep -v '^#') <(zcat $R/hc_$x.vcf.gz | grep -v '^#') && echo identical || echo DIFFERENT) ($(zcat $R/hc_$x.vcf.gz | grep -vc '^#') records, $(grep -o 'PairHMM [A-Z0-9-]*' $R/hc_$x.log | head -1), $(grep -o 'total=[0-9.]*s' $R/hc_$x.log | tail -1); $(eng $R/hc_$x.log))"
  done
  M=$P/markdup/multi
  d $VI gatk --java-options "-Xmx4g -XX:-UsePerfData" MarkDuplicates $(sed "s#in/#$M/in/#g" $M/truth/lanes.cmd) --OUTPUT $R/md.bam --METRICS_FILE $R/md.metrics > $R/md.log 2>&1
  say "MarkDuplicates, 12 lane BAMs" "$(python3 $REPO/tests/harness/bqsr/cmpbam.py $M/truth/lanes.bam $R/md.bam 2>&1 | tail -1 | grep -q ', 0 differ$' && diff -q <(grep -v '^#' $M/truth/lanes.metrics) <(grep -v '^#' $R/md.metrics) > /dev/null && echo identical || echo DIFFERENT) (records, metrics; $(eng $R/md.log))"
  V=$P/mergevcfs
  d $VI gatk --java-options "-Xmx2g" MergeVcfs --INPUT $V/in/v1.vcf.gz --INPUT $V/in/v2.vcf.gz --INPUT $V/in/v3.vcf.gz --OUTPUT $R/mv.vcf.gz --SEQUENCE_DICTIONARY $H/ref/GRCh38_no_alt_analysis_set.dict --TMP_DIR . > $R/mv.log 2>&1
  say "MergeVcfs (sarek's form)" "$(cmp -s <(zcat $V/truth/sarek.vcf.gz) <(zcat $R/mv.vcf.gz) && echo identical || echo DIFFERENT) ($(eng $R/mv.log))"
  C=$P/somatic/contam
  d $VI gatk --java-options "-Xmx2g" GetPileupSummaries --input $P/wgs/bam/NA12878.bam --variant $C/res/exac_chr20.hg38.vcf.gz --output $R/gps.table --reference $P/bqsr/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa --intervals $C/sarek/in/p2.bed --tmp-dir . > $R/gps.log 2>&1
  say "GetPileupSummaries (BED piece)" "$(cmp -s $C/sarek/truth/bam_p2.pileups.table $R/gps.table && echo identical || echo DIFFERENT) ($(eng $R/gps.log))"
  DR=$P/somatic/dream/data
  d $GPU $VI gatk --java-options "-Xmx4g" Mutect2 --input $DR/tumor_1.bam --input $DR/normal_1.bam --normal-sample $(samtools view -H $DR/normal_1.bam | grep -m1 '^@RG' | grep -o 'SM:[^[:space:]]*' | cut -d: -f2) --intervals 20 --germline-resource $DR/very-small-gnomad.vcf --f1r2-tar-gz $R/m2.f1r2.tar.gz --output $R/m2.vcf --reference $DR/human_g1k_v37.20.21.fasta --tmp-dir . > $R/m2.log 2>&1
  say "Mutect2 (DREAM 1)" "$(cmp -s <(grep -v '^##' $P/somatic/dream/truth/m2_1.vcf) <(grep -v '^##' $R/m2.vcf) && echo identical || echo DIFFERENT) ($(grep -vc '^#' $R/m2.vcf) records; $(eng $R/m2.log))"
fi
echo; [ $bad = 0 ] && echo "images: everything identical" || echo "images: $bad difference(s)"
exit $bad
