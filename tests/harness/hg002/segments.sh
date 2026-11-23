#!/usr/bin/env bash
# Contigs inside one RNG stream, run in parallel. Each run must stay byte-identical to one GATK process over the same
# contigs, the truth of the whole-contig shards (run_gatk.sh): shard043 = chr21 + chr22, shard046 = 139 small contigs
# (deep unplaced contigs, drawing groups close to the cuts).
# VCF (variant-rs/src/segment.rs): whole contigs are cut into pieces. Small pieces (VARIANT_RS_SEGMENT_MB) put many
# more cuts near drawing groups, so segments take over their neighbours (absorbed>0 below). A run that cuts nothing
# fails.
# GVCF (rng::Replay): each contig starts from the generator state a replay of the earlier contigs' drawing groups
# gives, and checks it at its end. A run without the replay fails.
#   segments.sh [OUTDIR]
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/hg002" || exit 1

X=${VARIANT_RS:-$REPO/target/release/variant-rs}; O=${1:-/tmp/varcall-segments}; mkdir -p $O
REF=ref/GRCh38_no_alt_analysis_set.fa; BAM=bam/HG002.bam
fail=0
one() { # name shard piece_mb (empty: the default size)
  local name=$1 s=$2 mb=$3 t=gatk/vcf/$2.vcf.gz out=$O/$1.vcf.gz t0=$(date +%s)
  env VARIANT_RS_CONTIG_TIMES=1 ${mb:+VARIANT_RS_SEGMENT_MB=$mb} \
    $X HaplotypeCaller -R $REF -I $BAM -L plan/$s.intervals -O $out > $O/$name.log 2>&1; local rc=$?
  local cut=$(grep -o '^segments: [0-9]* contig' $O/$name.log | grep -o '[0-9]*')
  local jobs=$(grep -o 'cut into [0-9]* jobs' $O/$name.log | grep -o '[0-9]*')
  local absorbed=$(grep '^contig-time' $O/$name.log | grep -c 'seg=Some(([0-9]*, [0-9]*, [1-9]')
  local n=$(zcat $t | grep -vc '^#') r
  if [ $rc != 0 ]; then r="FAILED (rc=$rc)"
  elif [ "${cut:-0}" = 0 ]; then r="DIFFERENT: no contig was cut"
  elif cmp -s <(zcat $t | grep -v '^#') <(zcat $out | grep -v '^#'); then r="identical ($n)"
  else r="DIFFERENT ($(diff <(zcat $t | grep -v '^#') <(zcat $out | grep -v '^#') | grep -c '^[<>]') lines of $n)"; fi
  [[ $r == identical* ]] || fail=$((fail + 1))
  echo "  $name: $s pieces=${mb:-default} cut=${cut:-0} jobs=${jobs:-0} absorbed=$absorbed $(( $(date +%s) - t0 ))s records $r"
}
one chr21_22         shard043 ""
one chr21_22_small   shard043 0.3
one small_contigs    shard046 0.02
gvcf() { # name shard
  local name=$1 s=$2 t=gatk/gvcf/$2.g.vcf.gz out=$O/$1.g.vcf.gz t0=$(date +%s)
  $X HaplotypeCaller -R $REF -I $BAM -L plan/$s.intervals -ERC GVCF -O $out > $O/$name.log 2>&1; local rc=$?
  local rp=$(grep -o '^replay: [0-9]* contig(s), [0-9]* drawing groups' $O/$name.log)
  local n=$(zcat $t | grep -vc '^#') r
  if [ $rc != 0 ]; then r="FAILED (rc=$rc)"
  elif [ -z "$rp" ]; then r="DIFFERENT: no replay"
  elif cmp -s <(zcat $t | grep -v '^#') <(zcat $out | grep -v '^#'); then r="identical ($n)"
  else r="DIFFERENT ($(diff <(zcat $t | grep -v '^#') <(zcat $out | grep -v '^#') | grep -c '^[<>]') lines of $n)"; fi
  [[ $r == identical* ]] || fail=$((fail + 1))
  echo "  $name: $s ${rp#replay: } $(( $(date +%s) - t0 ))s records $r"
}
gvcf chr21_22_gvcf      shard043
gvcf small_contigs_gvcf shard046
exit $fail
