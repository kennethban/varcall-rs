#!/usr/bin/env bash
# GenomicsDB workspaces made by GATK's GenomicsDBImport (JointGenotypingTasks.wdl ImportGVCFs command) and GATK
# GenotypeGVCFs truth on them (WARP's command, -V gendb://). Inputs: genotypegvcfs/multi (3 samples, chr8:1-15M; raw
# and reblocked) and multi/in11 (11 samples, chr8:1-5M, reblocked).
# ws_rb / ws_raw: 3 samples, batch 50 + --consolidate (WARP); ws_rb_b1: batch 1, no consolidation (3 fragments);
# ws11: 11 samples, batch 50 + consolidate. Truth: <ws>.vcf.gz (WARP arguments), <ws>.default.vcf.gz (no -G / -D).
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/genomicsdb" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
REF=$VARCALL_DATA/bqsr/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
DBSNP=$VARCALL_DATA/bqsr/resources/Homo_sapiens_assembly38.dbsnp138.vcf.gz
M=$VARCALL_DATA/genotypegvcfs/multi
mkdir -p truth
printf "HG002\t$M/in/HG002.rb.g.vcf.gz\nNA12878\t$M/in/NA12878.rb.g.vcf.gz\nHG002pcr\t$M/in/HG002pcr.rb.g.vcf.gz\n" > rb.map
printf "HG002\t$M/in/HG002.g.vcf.gz\nNA12878\t$M/in/NA12878.g.vcf.gz\nHG002pcr\t$M/in/HG002pcr.g.vcf.gz\n" > raw.map
for f in $M/in11/*_?.rb.g.vcf.gz; do s=$(basename $f .rb.g.vcf.gz); printf "$s\t$f\n"; done > rb11.map
imp() { # ws map interval extra
  [ -f $1/callset.json ] && return 0
  rm -rf $1
  java -Xmx8g -jar $J GenomicsDBImport --genomicsdb-workspace-path $1 -L $3 --sample-name-map $2 --reader-threads 5 --merge-input-intervals $4 > $1.import.log 2>&1 || echo "$1 import FAILED"
}
imp ws_rb rb.map chr8:1-15000000 "--batch-size 50 --consolidate" &
imp ws_raw raw.map chr8:1-15000000 "--batch-size 50 --consolidate" &
imp ws_rb_b1 rb.map chr8:1-15000000 "--batch-size 1" &
imp ws11 rb11.map chr8:1-5000000 "--batch-size 50 --consolidate" &
wait
gg() { # ws interval
  local w="-D $DBSNP -G StandardAnnotation -G AS_StandardAnnotation --only-output-calls-starting-in-intervals -V gendb://$1 -L $2 --merge-input-intervals"
  echo "$w" | $HARNESS/cmdargs.sh -w > truth/$1.cmd
  [ -f truth/$1.vcf.gz.tbi ] || java -Xmx6g -jar $J GenotypeGVCFs -R $REF $w -O truth/$1.vcf.gz > truth/$1.log 2>&1 || echo "$1 FAILED"
  echo "-V gendb://$1 -L $2" | $HARNESS/cmdargs.sh -w > truth/$1.default.cmd
  [ -f truth/$1.default.vcf.gz.tbi ] || java -Xmx6g -jar $J GenotypeGVCFs -R $REF -V gendb://$1 -L $2 -O truth/$1.default.vcf.gz > truth/$1.default.log 2>&1 || echo "$1 default FAILED"
}
gg ws_rb chr8:1-15000000 & gg ws_raw chr8:1-15000000 & gg ws_rb_b1 chr8:1-15000000 & gg ws11 chr8:1-5000000 &
wait
for f in truth/*.vcf.gz; do echo "$f $(zcat $f | grep -vc '^#')"; done
