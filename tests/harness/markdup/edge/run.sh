#!/usr/bin/env bash
# MarkDuplicates edge cases, variant-rs vs Picard (GATK 4.7's MarkDuplicates): records, header, @PG chain, metrics.
#   empty         header only (sarek hands MarkDuplicates such a lane when a FASTQ chunk aligns nothing)
#   empty+tiny    an empty lane merged with a non-empty one
#   unmapped      only unplaced unmapped reads (Picard counts the first one's PG before it stops)
#   noindex       empty and without a .bai (built in memory)
#   single-end    a single-end lane (the HCC1143 workshop tumor); se+pe merged with a paired lane
#   pos0          an unmapped read keeping its contig at POS 0 (pos0.py), first on chr20: copied through (Picard with
#                 VALIDATION_STRINGENCY SILENT, as WARP runs it: STRICT, the default, rejects the record)
#   hextag        an H (hex) tag on every read: htsjdk decodes it to a byte[] and writes it back as B:c
#   hiqual        quality bytes 200 in every other read (hiqual.py): Java bytes, negative, not in the duplicate score
# Inputs are made from the pipeline test's lane BAMs (pipeline/nf/test, --save_lane_bams).  usage: run.sh
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/markdup/edge" || exit 1; mkdir -p in out

L=$VARCALL_DATA/pipeline/nf/test/out
B=$L/germline/lanes/NA12878.L1.bam
[ -f in/tiny.bam ] || {
  samtools view -H -b -o in/empty.bam $B && samtools index in/empty.bam && cp in/empty.bam in/noidx.bam
  samtools view -b -o in/tiny.bam $B chr20:10000000-10001000 && samtools index in/tiny.bam
  samtools view -h $B chr20:10000000-10100000 | awk 'BEGIN{OFS="\t"} /^@/{print;next} {f=or($2,4); if(and(f,1)) f=or(f,8); f=and(f,compl(2+16+32+256+2048)); $2=f; $3="*";$4=0;$5=0;$6="*";$7="*";$8=0;$9=0; print}' \
    | samtools sort -o in/unmapped.bam - 2>/dev/null && samtools index in/unmapped.bam
  cp $L/somatic/lanes/HCC1143_tumor.C0HV4.4.bam in/se.bam && samtools index in/se.bam
}
[ -f in/hextag.bam ] || { samtools view -h in/tiny.bam | awk 'BEGIN{OFS="\t"} /^@/ {print; next} {print $0, "XH:H:1aef0b"}' | samtools view -b -o in/hextag.bam - && samtools index in/hextag.bam; }
[ -f in/pos0.bam ] || { python3 $HARNESS/markdup/edge/pos0.py in/tiny.bam in/pos0.bam && samtools index in/pos0.bam; }
[ -f in/hiqual.bam ] || { samtools view -b -o in/hisrc.bam $B chr20:10000000-10300000 && python3 $HARNESS/markdup/edge/hiqual.py in/hisrc.bam in/hiqual.bam && samtools index in/hiqual.bam && rm in/hisrc.bam; }
G=${VARIANT_RS:-$REPO/target/release/variant-rs}; J="java -Xmx4g -jar ${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}"
pgids() { samtools view -H --no-PG $1 2>/dev/null | grep '^@PG' | grep -o -P '\t(ID|PP):[^\t]*' | tr -d '\n'; }
bad=0
run() { local name=$1; shift; local ins=""; for b in "$@"; do ins="$ins -I in/$b"; done; local o=out/$name
  $G markdup $ins -O $o.rs.bam -M $o.rs.metrics --tmp-dir out > $o.rs.log 2>&1; rr=$?
  $J MarkDuplicates $ins -O $o.pic.bam -M $o.pic.metrics $PX > $o.pic.log 2>&1; pr=$?
  # record bytes (bin field included: SAM text does not show it), as bqsr/cmpbam.py compares them
  rec=$( python3 $HARNESS/bqsr/cmpbam.py $o.pic.bam $o.rs.bam 2>/dev/null | grep -q ' 0 differ' && echo same || echo DIFF )
  hd=$( [ "$(samtools view -H --no-PG $o.rs.bam 2>/dev/null | grep -v '^@PG' | md5sum)" = "$(samtools view -H --no-PG $o.pic.bam 2>/dev/null | grep -v '^@PG' | md5sum)" ] && echo same || echo DIFF )
  pg=$( [ "$(pgids $o.rs.bam)" = "$(pgids $o.pic.bam)" ] && echo same || echo DIFF )
  met=$( diff <(grep -v '^#' $o.rs.metrics 2>/dev/null) <(grep -v '^#' $o.pic.metrics 2>/dev/null) > /dev/null && echo same || echo DIFF )
  [ "$rr$pr$rec$hd$pg$met" = "00samesamesamesame" ] || bad=$((bad + 1))
  printf "%-11s variant-rs rc=%s Picard rc=%s | records %s (%s) | header %s | @PG %s | metrics %s\n" $name $rr $pr $rec $(samtools view -c $o.rs.bam 2>/dev/null) $hd $pg $met; }
run empty empty.bam; run empty+tiny empty.bam tiny.bam; run unmapped unmapped.bam; run noindex noidx.bam
run tiny tiny.bam; run single-end se.bam; run se+pe se.bam tiny.bam; run hiqual hiqual.bam; run hextag hextag.bam; PX='--VALIDATION_STRINGENCY SILENT' run pos0 pos0.bam
echo "$bad case(s) differ"
