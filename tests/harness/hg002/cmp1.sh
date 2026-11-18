#!/usr/bin/env bash
# cmp1.sh RUST_FILE TRUTH_DIR EXT: data lines byte-identical? header identical except the command line?
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
f=$1; ext=$3; s=$(basename $f .$ext); t=$2/$s.$ext
[ -f $t.tbi ] || { echo "$s MISSING 0 hdr=?"; exit; }
n=$(zcat $t | grep -vc "^#")
if cmp -s <(zcat $t | grep -v "^#") <(zcat $f | grep -v "^#"); then r=IDENTICAL; else r="DIFFERENT($(diff <(zcat $t | grep -v '^#') <(zcat $f | grep -v '^#') | grep -c '^[<>]'))"; fi
h=$(diff <(zcat $t | grep "^#" | grep -v "^##GATKCommandLine") <(zcat $f | grep "^#" | grep -v "^##variant-rsCommandLine") | wc -l)
echo "$s $r $n hdr=$h"
