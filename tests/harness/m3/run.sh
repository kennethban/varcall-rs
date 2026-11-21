#!/usr/bin/env bash
# M3: the CEPH 1463 trio (NA12878 and its parents NA12891, NA12892), 30x WGS from FASTQ (12 lanes each, 36 in all),
# through the pipeline side by side, with joint germline calling (GVCFs, GenomicsDBImport, GenotypeGVCFs, VQSR) over
# the GATK whole-genome calling regions. Lane BAMs are kept for truth.sh. Timings: RUN/trace.tsv, timeline.html.
#   usage: run.sh RUN_DIR [nextflow arguments, e.g. -resume]   (PIPELINE: the pipeline folder, default this checkout's)
# Inputs: wgs/fastq (NA12878: wgs/fastq/extract_all.sh; the parents: wgs/trio.sh), ref_hg38, bqsr/resources,
# vqsr/resources, somatic/wgs_ivs.
set -euo pipefail
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
PIPELINE=${PIPELINE:-$REPO/pipeline}
O=${1:?usage: run.sh RUN_DIR}; shift; mkdir -p "$O"; cd "$O"
D=$VARCALL_DATA; F=$D/wgs/fastq; B=$D/bqsr/resources; V=$D/vqsr/resources

# samplesheet: one row per lane, lane name <flowcell>.<lane number> as sarek's
{
  echo patient,sample,status,lane,fastq_1,fastq_2
  for s in NA12878 NA12891 NA12892; do
    for r1 in "$F"/${s}_*_R1.fq.gz; do
      fl=$(basename "$r1" _R1.fq.gz | cut -d_ -f2-)          # HCLHLDSXX_L001
      echo "$s,$s,0,${fl%_L*}.$((10#${fl##*_L})),$r1,${r1%_R1.fq.gz}_R2.fq.gz"
    done
  done
} > trio.csv
n=$(($(wc -l < trio.csv) - 1)); [ $n = 36 ] || { echo "trio.csv: $n lanes, 36 expected"; exit 1; }

cat > trio.yaml <<EOF
input: trio.csv
outdir: out
tools: haplotypecaller
step: mapping
save_lane_bams: true
fasta: $D/ref_hg38/Homo_sapiens_assembly38.fasta
known_sites: $B/Homo_sapiens_assembly38.dbsnp138.vcf.gz,$B/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz,$B/Homo_sapiens_assembly38.known_indels.vcf.gz
intervals: $D/somatic/wgs_ivs/wgs_calling_regions.hg38.interval_list
dbsnp: $B/Homo_sapiens_assembly38.dbsnp138.vcf.gz
joint_germline: true
vqsr_snp_resources: "--resource:hapmap,known=false,training=true,truth=true,prior=15.0 $V/hapmap_3.3.hg38.vcf.gz --resource:omni,known=false,training=true,truth=true,prior=12.0 $V/1000G_omni2.5.hg38.vcf.gz --resource:1000G,known=false,training=true,truth=false,prior=10.0 $V/1000G_phase1.snps.high_confidence.hg38.vcf.gz --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $B/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
vqsr_indel_resources: "--resource:mills,known=false,training=true,truth=true,prior=12.0 $B/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $B/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
EOF

# temporary files in the run folder, not /tmp
export NXF_ANSI_LOG=false TMPDIR=$PWD/tmp
mkdir -p "$TMPDIR"
t0=$(date +%s)
"$PIPELINE/run.sh" -params-file trio.yaml -w work -with-trace trace.tsv -with-timeline timeline.html \
  -with-report report.html "$@" > run.log 2>&1 && rc=0 || rc=$?
echo "exit $rc, wall $(( $(date +%s) - t0 )) s" | tee done
