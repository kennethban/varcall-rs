#!/usr/bin/env bash
# Picard MergeVcfs (GATK 4.7's bundled Picard) truth. Inputs as WARP makes them: HaplotypeCaller per scatter interval
# with WARP's GVCF / VCF command (HG002 chr8 in three stretches), then GermlineVariantDiscovery.wdl's MergeVCFs.
# rev: shards out of order; tie: VCF, GVCF and -contamination VCF of one stretch; twice: one shard twice (equal positions: PriorityQueue order); fmc: a FilterMutectCalls
# VCF (two ##GATKCommandLine lines) with a GVCF-less input; plain: uncompressed output without an index; sarek_*,
# dict_other: SEQUENCE_DICTIONARY.
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/mergevcfs" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
H=$VARCALL_DATA/hg002; REF=$H/ref/GRCh38_no_alt_analysis_set.fa; BAM=$H/bam/HG002.bam
GQB="-GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90"
WG="-contamination 0 -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation $GQB -ERC GVCF"
WV="-contamination 0 -G StandardAnnotation -G StandardHCAnnotation $GQB"
mkdir -p in truth
hc() { # name interval ext args
  [ -f in/$1.$3.tbi ] || java -Xmx3g -XX:ParallelGCThreads=2 -jar $J HaplotypeCaller -R $REF -I $BAM -L $2 -O in/$1.$3 $4 --native-pair-hmm-threads 2 > in/$1.log 2>&1 || echo "$1 FAILED"
}
hc g1 chr8:1-5000000 g.vcf.gz "$WG" & hc g2 chr8:5000001-10000000 g.vcf.gz "$WG" & hc g3 chr8:10000001-15000000 g.vcf.gz "$WG" &
hc v1 chr8:1-5000000 vcf.gz "$WV" & hc v2 chr8:5000001-10000000 vcf.gz "$WV" & hc v3 chr8:10000001-15000000 vcf.gz "$WV" &
hc c1 chr8:1-5000000 vcf.gz "-contamination 0.05" &
wait
F=$VARCALL_DATA/somatic/group3/fmc
declare -A C
C[gvcf]="INPUT=in/g1.g.vcf.gz INPUT=in/g2.g.vcf.gz INPUT=in/g3.g.vcf.gz OUTPUT=gvcf.g.vcf.gz"
C[vcf]="INPUT=in/v1.vcf.gz INPUT=in/v2.vcf.gz INPUT=in/v3.vcf.gz OUTPUT=vcf.vcf.gz"
C[rev]="INPUT=in/g3.g.vcf.gz INPUT=in/g1.g.vcf.gz INPUT=in/g2.g.vcf.gz OUTPUT=rev.g.vcf.gz"
C[twice]="INPUT=in/v2.vcf.gz INPUT=in/v1.vcf.gz INPUT=in/v2.vcf.gz OUTPUT=twice.vcf.gz"
# tie: different records at equal positions in four inputs (the heap decides their order)
C[tie]="INPUT=in/v1.vcf.gz INPUT=in/g1.g.vcf.gz INPUT=in/c1.vcf.gz INPUT=in/v1.vcf.gz OUTPUT=tie.g.vcf.gz"
C[fmc]="INPUT=$F/to_opts.vcf INPUT=$F/to_mito.vcf OUTPUT=fmc.vcf.gz"
C[plain]="INPUT=in/v1.vcf.gz INPUT=in/v3.vcf.gz OUTPUT=plain.vcf CREATE_INDEX=false"
# nf-core/sarek's command (gatk4/mergevcfs): the reference's .dict as SEQUENCE_DICTIONARY; dict_other: a larger
# dictionary than the headers' contigs (only the index sees it)
C[sarek_vcf]="--INPUT in/v1.vcf.gz --INPUT in/v2.vcf.gz --INPUT in/v3.vcf.gz --OUTPUT sarek.vcf.gz --SEQUENCE_DICTIONARY $H/ref/GRCh38_no_alt_analysis_set.dict --TMP_DIR ."
C[sarek_gvcf]="--INPUT in/g1.g.vcf.gz --INPUT in/g2.g.vcf.gz --INPUT in/g3.g.vcf.gz --OUTPUT sarek.g.vcf.gz --SEQUENCE_DICTIONARY $H/ref/GRCh38_no_alt_analysis_set.dict --TMP_DIR ."
C[dict_other]="INPUT=in/v3.vcf.gz INPUT=in/v1.vcf.gz OUTPUT=dict_other.vcf.gz SEQUENCE_DICTIONARY=$VARCALL_DATA/bqsr/ref/GRCh38_full_analysis_set_plus_decoy_hla.dict"
for c in "${!C[@]}"; do
  echo "${C[$c]}" > truth/$c.cmd
  o=$(echo "${C[$c]}" | grep -oE '(OUTPUT=|--OUTPUT )[^ ]*' | sed 's/.*[= ]//')
  [ -f truth/$o ] && continue
  java -Xmx2g -jar $J MergeVcfs $(echo "${C[$c]}" | sed "s|OUTPUT=|OUTPUT=truth/|; s|--OUTPUT |--OUTPUT truth/|") > truth/$c.log 2>&1; echo "$c rc=$?"
done
