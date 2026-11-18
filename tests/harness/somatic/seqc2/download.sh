#!/usr/bin/env bash
# SEQC2 HCC1395 WES (Fudan replicate 1) + resources. Resumable.
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$(dirname "$0")/data"
S=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG
G=https://storage.googleapis.com/gatk-best-practices/somatic-hg38
get() { for i in 1 2 3 4 5; do curl -sSfL --retry 5 -C - -o "$(basename "$1")" "$1" && return 0; sleep 10; done; echo "FAILED $1"; }
date "+start %H:%M"
for f in WES_FD_T_1.bwa.dedup.bam WES_FD_T_1.bwa.dedup.bai WES_FD_N_1.bwa.dedup.bam WES_FD_N_1.bwa.dedup.bai; do get $S/data/WES/$f; done
for f in GRCh38.d1.vd1.fa GRCh38.d1.vd1.dict; do get $S/technical/reference_genome/GRCh38/$f; done
get $S/technical/reference_genome/Exome_Target_bed/S07604624_Covered_human_all_v6_plus_UTR.liftover.to.hg38.bed6.gz
for f in high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz High-Confidence_Regions_v1.2.bed; do get $S/release/latest/$f; done
for f in af-only-gnomad.hg38.vcf.gz af-only-gnomad.hg38.vcf.gz.tbi 1000g_pon.hg38.vcf.gz 1000g_pon.hg38.vcf.gz.tbi small_exac_common_3.hg38.vcf.gz small_exac_common_3.hg38.vcf.gz.tbi; do get $G/$f; done
date "+done %H:%M"
ls -la
