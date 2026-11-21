#!/usr/bin/env bash
# variant-rs HaplotypeCaller --scatter-intervals vs GATK run once per interval.
#  1. An interval file with a duplicated, an overlapping, a contained, a touching pair and out-of-order lines.
#     Expected: overlapping intervals merged (one GATK run), touching ones separate runs, in dictionary order; the
#     truth is GATK on each expected piece, records concatenated (made once into tests/scatter_truth/).
#  2. Many small pieces (the chr22 exome targets of the pipeline's somatic_bam test, ~6,600 runs) twice: both runs
#     finish and agree (the two-stage front end's end-of-contig race was hit once in 294,884 such runs).
# usage: scatter_intervals.sh [OUTDIR]
REPO=$(cd "$(dirname "$0")/../.." && pwd); P=${VARCALL_DATA:?set VARCALL_DATA to the test data directory}; G=${VARIANT_RS:-$REPO/target/release/variant-rs}
O=${1:-$P/scatter/out}; T=$P/scatter/truth; mkdir -p $O $T
GATK="java -Xmx8g -jar ${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}"
R=$P/ref_hg38/Homo_sapiens_assembly38.fasta; BAM=$P/pipeline/nf/test/truth/germline/NA12878.recal.bam
cat > $O/odd.intervals <<'EOF'
chr20:12000000-12300000
chr20:10100000-10200000
chr20:10500000-10700000
chr20:10100000-10200000
chr20:10650000-10800000
chr20:11000000-11100000
chr20:11100001-11200000
chr20:12100000-12150000
EOF
pieces="chr20:10100000-10200000 chr20:10500000-10800000 chr20:11000000-11100000 chr20:11100001-11200000 chr20:12000000-12300000"
if [ ! -f $T/odd.vcf.gz ]; then
  parts=(); i=0
  for iv in $pieces; do i=$((i+1)); $GATK HaplotypeCaller -R $R -I $BAM -L $iv -O $T/p$i.vcf.gz > $T/p$i.log 2>&1 || { echo "GATK FAILED on $iv"; exit 1; }; parts+=($T/p$i.vcf.gz); done
  { zcat ${parts[0]} | grep "^#"; for p in "${parts[@]}"; do zcat $p | grep -v "^#"; done; } | bgzip > $T/odd.vcf.gz
fi
$G -I $BAM -R $R -L $O/odd.intervals --scatter-intervals -O $O/odd.vcf.gz --threads 16 > $O/odd.log 2>&1; rc=$?
body() { zcat $1 | grep -v -E "^##(GATKCommandLine|source|fileDate|variant-rs)"; }
echo "  odd intervals: rc=$rc $(grep -o '[0-9]* contig job(s) in [0-9]* domain(s)' $O/odd.log) (expected 5 domains); VCF $(cmp -s <(body $O/odd.vcf.gz) <(body $T/odd.vcf.gz) && echo identical || echo DIFFERENT) ($(zcat $T/odd.vcf.gz | grep -vc '^#') records)"

W=$P/somatic/seqc2/data; D=$P/pipeline/nf/test/data/wes
for k in 1 2; do
  $G -I $D/FD_N_1.chr22.bam -R $W/GRCh38.d1.vd1.fa -L $D/targets.chr22.bed --scatter-intervals -O $O/many$k.vcf.gz --threads 32 > $O/many$k.log 2>&1
  echo "  many pieces, run $k: rc=$? $(grep -o '[0-9]* domain(s)' $O/many$k.log), $(grep -o 'time: total=[0-9.]*s' $O/many$k.log)"
done
echo "  many pieces: the two runs $(cmp -s <(body $O/many1.vcf.gz) <(body $O/many2.vcf.gz) && echo identical || echo DIFFERENT) ($(zcat $O/many1.vcf.gz | grep -vc '^#') records)"
