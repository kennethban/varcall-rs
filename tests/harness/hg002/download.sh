#!/usr/bin/env bash
# HG002 validation data (public). Resumable: re-run to continue.
set -o pipefail
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/hg002" || exit 1
get() { local url=$1 out=$2; for i in 1 2 3 4 5; do curl -sSL --retry 5 -C - -o "$out" "$url" && return 0; sleep 10; done; echo "FAILED $url"; return 1; }
G=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38
N=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
D=https://storage.googleapis.com/deepvariant/case-study-testdata
date "+start %H:%M"
for f in HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed; do get $G/$f giab/$f; done
get $N/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz ref/GRCh38_no_alt_analysis_set.fna.gz
date "+ref done %H:%M"
# Google DeepVariant case study: NovaSeq PCR-free 35x, bwa mem 0.7.17 + Picard MarkDuplicates, no BQSR
get $D/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam.bai bam/HG002.bam.bai
get $D/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam bam/HG002.bam
date "+bam done %H:%M"
ls -la giab ref bam
