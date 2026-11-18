#!/usr/bin/env bash
# variant-rs MergeVcfs with each truth case's arguments vs Picard: text (decompressed) byte for byte, and the .tbi
# through region queries (its bytes differ with the BGZF blocks).  usage: cmp.sh [OUTDIR] [case...]
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/mergevcfs" || exit 1; D=$PWD

G=${VARIANT_RS:-$REPO/target/release/variant-rs}
O=${1:-$D/rust}; shift; mkdir -p $O; O=$(realpath $O)
txt() { case $1 in *.gz) zcat $1;; *) cat $1;; esac; }
for c in ${@:-$(ls truth/*.cmd | xargs -n1 basename | sed 's/\.cmd$//' | sort)}; do
  args=$(cat truth/$c.cmd); o=$(echo "$args" | grep -oE '(OUTPUT=|--OUTPUT )[^ ]*' | sed 's/.*[= ]//')
  $G MergeVcfs $(echo "$args" | sed "s|OUTPUT=|OUTPUT=$O/|; s|--OUTPUT |--OUTPUT $O/|") > $O/$c.log 2>&1; rc=$?
  n=$(txt truth/$o | grep -vc '^#')
  t=$(cmp -s <(txt truth/$o) <(txt $O/$o) && echo identical || echo "DIFFERENT ($(diff <(txt truth/$o) <(txt $O/$o) | grep -c '^[<>]') lines)")
  ix="no index"
  if [ -f truth/$o.tbi ]; then
    q() { tabix -l $1; for r in chr8:100000-200000 chr8:4999000-5001000 chr8:9999990-10000010 chr8:14000000-15000000 chr17; do tabix $1 $r | md5sum; done; }
    ix="index $([ -f $O/$o.tbi ] && (diff <(q truth/$o) <(q $O/$o) > /dev/null && echo "same answers" || echo "DIFFERENT answers") || echo MISSING)"
  elif [ -f $O/$o.tbi ]; then ix="index written but Picard wrote none"; fi
  echo "  $c: rc=$rc text $t ($n records); $ix ($(grep -o '([0-9.]*s)' $O/$c.log))"
done
