#!/usr/bin/env bash
# Self-contained test: synthetic data (tests/synth/simulate.py) through align-rs and variant-rs.
#
#   tests/selftest.sh [OUTDIR]          (default: target/selftest; the data are made once and kept there)
#
# Needs the release build (cargo build --release; or ALIGN_RS / VARIANT_RS = the binaries, e.g. containers/test_images.sh's
# wrappers that run them in the images), python3, bwa-mem2, samtools,
# bgzip and tabix. No downloads.
#   align-rs     records, header and .bai vs bwa-mem2 + samtools sort / index, run here; with a GPU, also the CPU path
#   variant-rs   GATK's command lines (tests/variant_steps.sh) vs GATK 4.7 / Picard truth in tests/truth, made by
#                tests/make_truth.sh from the same inputs: BAM records, VCF data lines, metrics, tables; with the
#                CPU PairHMM, and with the GPU one when there is a GPU
# Exit status: the number of differences (0: everything identical).
set -u
REPO=$(cd "$(dirname "$0")/.." && pwd); O=$(mkdir -p "${1:-$REPO/target/selftest}" && cd "${1:-$REPO/target/selftest}" && pwd)
D=$O/data; T=$REPO/tests/truth; A=${ALIGN_RS:-$REPO/target/release/align-rs}; V=${VARIANT_RS:-$REPO/target/release/variant-rs}
TH=$(nproc); bad=0
for c in python3 bwa-mem2 samtools bgzip tabix; do command -v $c > /dev/null || { echo "missing: $c"; exit 1; }; done
[ -x $A ] && [ -x $V ] || { echo "build first: cargo build --release"; exit 1; }
say() { printf '  %-44s %s\n' "$1" "$2"; [[ "$2" == identical* || "$2" == ok* ]] || bad=$((bad + 1)); }
gpu=0; command -v nvidia-smi > /dev/null && nvidia-smi -L > /dev/null 2>&1 && gpu=1

echo "== data ($D)"
if [ ! -s $D/tumor_R2.fq.gz ]; then
  t=$(date +%s); python3 $REPO/tests/synth/simulate.py $D || exit 1; echo "  simulated in $(( $(date +%s) - t )) s"
fi
(cd $D && md5sum -c --quiet $REPO/tests/synth/data.md5) && say "synthetic data vs tests/synth/data.md5" identical || say "synthetic data vs tests/synth/data.md5" "DIFFERENT (simulate.py or Python changed)"
if [ ! -s $D/ref.fa.bwt.2bit.64 ]; then bwa-mem2 index $D/ref.fa > $D/index.log 2>&1 || { echo "bwa-mem2 index failed"; exit 1; }; fi
for v in known_sites germline_resource; do [ -s $D/$v.vcf.gz.tbi ] || { bgzip -c $D/$v.vcf > $D/$v.vcf.gz && tabix -f -p vcf $D/$v.vcf.gz; }; done

echo "== align-rs vs bwa-mem2"
rg() { printf '@RG\\tID:%s\\tSM:%s\\tPL:ILLUMINA\\tLB:%s\\tPU:%s' "$1" "$2" "$2" "$3"; }
LANES="GERM.L1:GERM:FCGERM1.1:germline_L1 GERM.L2:GERM:FCGERM1.2:germline_L2 TUMOR:TUMOR:FCTUMR1.1:tumor NORMAL:NORMAL:FCNORM1.1:normal"
mkdir -p $O/bwa $O/align $O/align_cpu
: > $O/lanes.tsv; : > $O/lanes_cpu.tsv
for l in $LANES; do
  IFS=: read id sm pu fq <<< "$l"
  if [ ! -s $O/bwa/$id.bam ]; then
    bwa-mem2 mem -t $TH -K 100000000 -Y -R "$(rg $id $sm $pu)" $D/ref.fa $D/${fq}_R1.fq.gz $D/${fq}_R2.fq.gz 2> $O/bwa/$id.log \
      | samtools sort -@ 4 -l 1 -o $O/bwa/$id.bam - 2>> $O/bwa/$id.log && samtools index $O/bwa/$id.bam
  fi
  printf '%s\t%s\t%s\t%s\n' $O/align/$id.bam "$(rg $id $sm $pu)" $D/${fq}_R1.fq.gz $D/${fq}_R2.fq.gz >> $O/lanes.tsv
  printf '%s\t%s\t%s\t%s\n' $O/align_cpu/$id.bam "$(rg $id $sm $pu)" $D/${fq}_R1.fq.gz $D/${fq}_R2.fq.gz >> $O/lanes_cpu.tsv
done
t=$(date +%s); $A mem -t $TH -K 100000000 -Y --gpu -l 1 --lanes $O/lanes.tsv $D/ref.fa 2> $O/align/align.log; rc=$?
echo "  align-rs --lanes (4 lanes, $( [ $gpu = 1 ] && echo GPU || echo CPU)): exit $rc, $(( $(date +%s) - t )) s"
recs() { samtools view $1 | md5sum | cut -c1-32; }
hdr() { samtools view -H --no-PG $1 | grep -v '^@PG' | md5sum | cut -c1-32; }
for l in $LANES; do
  id=${l%%:*}; a=$O/align/$id.bam; b=$O/bwa/$id.bam
  [ -s $a ] || { say "$id" "MISSING"; continue; }
  r=$( [ "$(recs $a)" = "$(recs $b)" ] && echo identical || echo DIFFERENT)
  h=$( [ "$(hdr $a)" = "$(hdr $b)" ] && echo identical || echo DIFFERENT)
  cp $a.bai $O/align/$id.own.bai; samtools index $a; i=$(cmp -s $a.bai $O/align/$id.own.bai && echo identical || echo DIFFERENT)
  say "$id ($(samtools view -c $a) records)" "$( [ $r$h$i = identicalidenticalidentical ] && echo identical || echo "records $r, header $h, .bai $i") (records, header, .bai)"
done
if [ $gpu = 1 ]; then
  $A mem -t $TH -K 100000000 -Y -l 1 --lanes $O/lanes_cpu.tsv $D/ref.fa 2> $O/align_cpu/align.log
  for l in $LANES; do id=${l%%:*}
    say "$id, CPU path" "$( [ "$(recs $O/align_cpu/$id.bam)" = "$(recs $O/align/$id.bam)" ] && echo "identical (records, vs GPU)" || echo DIFFERENT)"
  done
fi

# variant-rs: the steps in a directory of their own, inputs linked
source $REPO/tests/variant_steps.sh
workdir() {
  rm -rf $1 && mkdir -p $1 && cd $1 || exit 1
  for f in ref.fa ref.fa.fai ref.dict known_sites.vcf.gz known_sites.vcf.gz.tbi germline_resource.vcf.gz germline_resource.vcf.gz.tbi; do ln -s $D/$f .; done
  for l in $LANES; do id=${l%%:*}; ln -s $O/align/$id.bam .; ln -s $O/align/$id.bam.bai .; done
}
run() { local n=$1; shift; "$@" > $n.log 2>&1 || echo "  $n: exit $? (see $PWD/$n.log)"; }
data() { zcat $1 | grep -v '^#' | md5sum | cut -c1-32; }
compare() {
  for f in $VARIANT_BAMS; do
    say "$f" "$( [ -s $f ] && [ "$(python3 $REPO/tests/bamdigest.py $f)" = "$(awk -v f=$f '$1 == f {print $2, $3}' $T/bam_digests.tsv)" ] && echo "identical (records)" || echo DIFFERENT)"
  done
  for f in $VARIANT_VCFS; do
    say "$f" "$( [ -s $f ] && [ "$(data $f)" = "$(zcat $T/$f.data.gz | md5sum | cut -c1-32)" ] && echo "identical ($(zcat $T/$f.data.gz | wc -l) data lines)" || echo DIFFERENT)"
  done
  for f in $VARIANT_METRICS; do
    say "$f" "$(cmp -s <(grep -v '^#' $f 2>/dev/null) $T/$f && echo "identical (without # lines)" || echo DIFFERENT)"
  done
  for f in $VARIANT_FILES; do say "$f" "$(cmp -s $f <(zcat $T/$f.gz) && echo identical || echo DIFFERENT)"; done
}
# the GPU PairHMM too when there is a GPU and variant-rs was built with it (make -C variant-rs/cuda-pairhmm);
# SELFTEST_PAIRHMM sets the list instead (e.g. "cpu cuda" for a variant-rs in a container)
vgpu=$gpu; ldd $V 2>/dev/null | grep -q libcudapairhmm || vgpu=0
for hmm in ${SELFTEST_PAIRHMM:-cpu $( [ $vgpu = 1 ] && echo cuda )}; do
  echo "== variant-rs (PairHMM: $hmm) vs GATK 4.7 / Picard"
  workdir $O/variant_$hmm; X=$V; HMM="--pairhmm $hmm"
  t=$(date +%s); variant_steps; echo "  all steps: $(( $(date +%s) - t )) s"
  [ -d $T ] && compare || echo "  no truth yet: tests/make_truth.sh"
  cd $O
done

echo; [ $bad = 0 ] && echo "selftest: everything identical" || echo "selftest: $bad difference(s)"
exit $bad
