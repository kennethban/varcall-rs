#!/usr/bin/env bash
# variant-rs plan on a BAM whose cost is mostly a short, very deep contig (HG002 chrM, 16,569 bp, ~2.5M reads) plus a
# little chr20 (review finding 33): chrM must not be cut into pieces shorter than a 16 kb window (the planner used to
# panic in clamp); the shards must cover every contig once, in order.
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/plan_edge" || exit 1

X=${VARIANT_RS:-$REPO/target/release/variant-rs}; REF=../hg002/ref/GRCh38_no_alt_analysis_set.fa
[ -s deep_mt.bam ] || { samtools merge -f -o deep_mt.bam <(samtools view -b ../hg002/bam/HG002.bam chr20:10000000-10100000) <(samtools view -b ../hg002/bam/HG002.bam chrM) && samtools index deep_mt.bam; }
O=${1:-/tmp/plan-edge}; rm -rf $O
for p in 8 64; do
  $X plan -I deep_mt.bam -R $REF -o $O/p$p --pieces $p > $O.p$p.log 2>&1; rc=$?
  mt=$(cat $O/p$p/shard*.intervals 2>/dev/null | grep -c '^chrM')
  cov=$(cat $O/p$p/shard*.intervals 2>/dev/null | sed -E 's/^(.*):([0-9]+)-([0-9]+)$/\2 \3/' | awk '{s += $2 - $1 + 1} END {print s + 0}')
  all=$(awk '{s += $2} END {print s}' $REF.fai)
  echo "  --pieces $p: rc=$rc, chrM in $mt interval(s), shards cover $( [ "$cov" = "$all" ] && echo "the genome once ($cov bp)" || echo "$cov of $all bp: DIFFERENT")"
done
