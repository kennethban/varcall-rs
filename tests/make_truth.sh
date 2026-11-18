#!/usr/bin/env bash
# GATK 4.7 / Picard truth for tests/selftest.sh: the steps of tests/variant_steps.sh run with GATK on the BAMs the
# self-test aligned (align-rs, identical to bwa-mem2), reduced to what the self-test compares and written to
# tests/truth: BAM record digests, VCF data lines (gzip), metrics without '#' lines, tables (gzip).
#
#   GATK_JAR=/path/to/gatk-package-4.7.0.0-local.jar tests/make_truth.sh [SELFTEST_OUTDIR]
#
# Run tests/selftest.sh first (it makes the data and the alignments). Needs Java 17. Takes ~15 min.
set -u
REPO=$(cd "$(dirname "$0")/.." && pwd); O=$(cd "${1:-$REPO/target/selftest}" && pwd); D=$O/data; T=$REPO/tests/truth
: "${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}"
[ -s $O/align/NORMAL.bam ] || { echo "run tests/selftest.sh first"; exit 1; }
source $REPO/tests/variant_steps.sh
W=$O/gatk; rm -rf $W && mkdir -p $W && cd $W || exit 1
for f in ref.fa ref.fa.fai ref.dict known_sites.vcf.gz known_sites.vcf.gz.tbi germline_resource.vcf.gz germline_resource.vcf.gz.tbi; do ln -s $D/$f .; done
for b in GERM.L1 GERM.L2 TUMOR NORMAL; do ln -s $O/align/$b.bam .; ln -s $O/align/$b.bam.bai .; done
run() { local n=$1; shift; local t=$(date +%s); "$@" > $n.log 2>&1; echo "  $n: exit $?, $(( $(date +%s) - t )) s"; }
X="java -Xmx8g -jar $GATK_JAR"; HMM=""
java -jar $GATK_JAR --version 2>/dev/null | head -1
variant_steps
mkdir -p $T && rm -f $T/*
for f in $VARIANT_BAMS; do echo "$f $(python3 $REPO/tests/bamdigest.py $f)"; done > $T/bam_digests.tsv
for f in $VARIANT_VCFS; do zcat $f | grep -v '^#' | gzip -n > $T/$f.data.gz; done
for f in $VARIANT_METRICS; do grep -v '^#' $f > $T/$f; done
for f in $VARIANT_FILES; do gzip -n -c $f > $T/$f.gz; done
echo "truth written to $T: $(ls $T | wc -l) files, $(du -sh $T | cut -f1)"
