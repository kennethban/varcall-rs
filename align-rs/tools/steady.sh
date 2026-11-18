#!/bin/bash
# steady-state time of an align-rs log: last batch done - first batch done, and per batch.   steady.sh LOG
grep "batch: " $1 | awk '{t=$NF; sub("s","",t); if (NR==1) f=t; l=t; n=NR} END {printf "batches %d, first %.1f s, steady %.1f s (%.2f s/batch)\n", n, f, l-f, (l-f)/(n-1)}'
