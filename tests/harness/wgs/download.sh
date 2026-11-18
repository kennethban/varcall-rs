#!/usr/bin/env bash
# Phase 1 whole-genome test data (public). Resumable: re-run to continue.
set -o pipefail
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/wgs" || exit 1
S3=https://1000genomes.s3.amazonaws.com
G=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38
get() { local url=$1 out=$2; for i in 1 2 3 4 5; do curl -sSL --retry 5 -C - -o "$out" "$url" && return 0; sleep 10; done; echo "FAILED $url"; return 1; }
date "+start %H:%M"
# truth set first (small)
get $G/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz      giab/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
get $G/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi  giab/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
get $G/HG001_GRCh38_1_22_v4.2.1_benchmark.bed         giab/HG001_GRCh38_1_22_v4.2.1_benchmark.bed
# reference used by the NYGC/1000G CRAMs
get $S3/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa     ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
get $S3/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai ref/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
date "+reference done %H:%M"
# NA12878 30x WGS (NYGC, BWA-MEM + MarkDuplicates + BQSR already applied)
get $S3/1000G_2504_high_coverage/data/ERR3239334/NA12878.final.cram      cram/NA12878.final.cram
get $S3/1000G_2504_high_coverage/data/ERR3239334/NA12878.final.cram.crai cram/NA12878.final.cram.crai
date "+cram done %H:%M"
ls -la giab ref cram
