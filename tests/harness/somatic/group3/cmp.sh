#!/usr/bin/env bash
# variant-rs Mutect2 with each truth case's GATK arguments (truth/<case>.cmd, see truth.sh) vs GATK: VCF (all but the
# command line), .stats, and the --bam-output BAM's records and header (its .bai cannot be byte-identical: the
# BGZF blocks compress differently).  usage: cmp.sh [OUTDIR] [case...]
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/group3" || exit 1; D=$PWD

G=${VARIANT_RS:-$REPO/target/release/variant-rs}
O=${1:-$D/rust}; shift; mkdir -p $O
cases=${@:-$(ls truth/*.cmd | xargs -n1 basename | sed 's/\.cmd$//' | sort)}
for c in $cases; do
  args=$($HARNESS/cmdargs.sh truth/$c.cmd | sed "s|truth/bamout.bam|$O/bamout.bam|")
  /usr/bin/time -f "%e s" $G Mutect2 $args -O $O/$c.vcf > $O/$c.log 2>&1
  rc=$?
  n=$(grep -vc '^#' truth/$c.vcf)
  v=$(diff <(grep -v CommandLine truth/$c.vcf) <(grep -v CommandLine $O/$c.vcf 2>/dev/null) > $O/$c.diff && echo identical || echo "DIFFERENT($(grep -c '^[<>]' $O/$c.diff))")
  st=$(cmp -s truth/$c.vcf.stats $O/$c.vcf.stats && echo identical || echo DIFFERENT)
  extra=""
  [[ $args == *--bam-output* ]] && extra="; bamout $(python3 $HARNESS/bqsr/cmpbam.py truth/bamout.bam $O/bamout.bam 2>&1 | tail -1), header $(diff -q <(samtools view -H --no-PG truth/bamout.bam) <(samtools view -H --no-PG $O/bamout.bam) > /dev/null && echo identical || echo DIFFERENT), index $([ -f $O/bamout.bai ] && echo written || echo MISSING)"
  echo "  $c: rc=$rc VCF $v ($n records), stats $st$extra ($(tail -1 $O/$c.log))"
done
