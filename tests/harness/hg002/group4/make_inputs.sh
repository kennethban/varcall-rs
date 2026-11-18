#!/usr/bin/env bash
# Test inputs for several -I (S3 group 4): HG002 chr8:1-15M split by read name (last character of the name, so
# mates stay together). lanes: L1.bam / L2.bam with read groups L1 / L2 of sample HG002; same_rg: S0.bam / S1.bam
# with the original read group; three: T0..T2.bam (read groups T0..T2); other: O.bam, L2's reads as sample OTHER.
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/hg002/group4" || exit 1

mkdir -p in
BAM=../bam/HG002.bam; R=chr8:1-15000000
RG="HG002.pcr-free"
# split <out> <n> <k> <new RG id or ""> <sample>: reads whose name's last digit % n == k
split() {
  local out=$1 n=$2 k=$3 id=$4 sm=$5
  samtools view -h $BAM $R 2>/dev/null | awk -v n=$n -v k=$k -v id="$id" -v sm=$sm -v rg="$RG" 'BEGIN{OFS="\t"}
    /^@RG/ { if (id != "") { sub("ID:" rg, "ID:" id); sub("SM:HG002", "SM:" sm) } print; next }
    /^@/ { print; next }
    { if (substr($1, length($1), 1) % n != k) next
      if (id != "") for (i = 12; i <= NF; i++) if ($i == "RG:Z:" rg) $i = "RG:Z:" id
      print }' | samtools view -b -o in/$out.bam - && samtools index in/$out.bam
}
split L1 2 0 L1 HG002 & split L2 2 1 L2 HG002 &
split S0 2 0 "" HG002 & split S1 2 1 "" HG002 &
split T0 3 0 T0 HG002 & split T1 3 1 T1 HG002 & split T2 3 2 T2 HG002 &
split O 2 1 L2 OTHER &
wait
ls -la in
