#!/usr/bin/env bash
# M3: the pipeline's outputs (RUN/out) vs GATK's (RUN/truth, truth.sh), by file name, all files at once.
# BAM: records and header apart from @PG; VCF: records and header apart from command / source / date lines;
# tables and metrics: lines apart from '#' comments (the command line).
#   usage: cmp.sh RUN_DIR        (OUT: the pipeline output folder in RUN_DIR, default out)
set -uo pipefail
O=$(realpath "${1:?usage: cmp.sh RUN_DIR}"); P=$O/${OUT:-out}; T=$O/truth
digest() {  # file -> "md5 records"
  case $1 in
    *.bam) echo "$(samtools view --no-PG -@ 8 $1 | md5sum | cut -c1-32) $(samtools view -c -@ 8 $1)" \
                "$(samtools view -H --no-PG $1 | grep -v '^@PG' | md5sum | cut -c1-8)" ;;
    *.vcf.gz) echo "$(zcat $1 | grep -v -E '^##(GATKCommandLine|variant-rs|source=|fileDate|filtering_status)' \
                | md5sum | cut -c1-32) $(zcat $1 | grep -vc '^#')" ;;
    *) grep -v '^#' $1 | md5sum | cut -c1-32 ;;
  esac | tr '\n' ' '
}
one() {  # pipeline file
  local a=$1 rel=${1#$P/} t=$T/$(basename $1)
  [ -e $t ] || { printf '  %-55s no GATK file\n' $rel; return; }
  local da=$(digest $a) dt=$(digest $t)
  [ "$da" = "$dt" ] && printf '  %-55s identical (%s)\n' $rel "$(echo $da | cut -d' ' -f2- | xargs)" \
                    || printf '  %-55s DIFFERENT (%s vs %s)\n' $rel "$da" "$dt"
}
export -f digest one; export P T
find $P -type f ! -path '*/lanes/*' ! -name '*.bai' ! -name '*.tbi' | sort | xargs -P 16 -I{} bash -c 'one {}' | sort > $O/cmp.out
cat $O/cmp.out
echo "$(wc -l < $O/cmp.out) files, $(grep -vc ' identical' $O/cmp.out) different or missing"
