#!/usr/bin/env bash
# M3 cohort data: NA12878's parents NA12891 and NA12892 (CEPH 1463 trio), 30x NovaSeq from the 1000 Genomes
# high-coverage set (NYGC), as FASTQ lanes beside NA12878's (fastq/extract_all.sh). Resumable: re-run to continue.
#   1. the CRAMs and their reference from EBI, in parallel byte ranges (one stream gets ~1 MB/s; 48 at once got the
#      server refused for a while, so JOBS is 12), md5 / size checked (a CRAM fetched another way, e.g. from the
#      1000 Genomes S3 bucket, can be put in cram/ first: it is checked too);
#   2. each CRAM split by read group (12 lanes), each lane collated by name and written as R1/R2, primary records
#      only, as fastq/extract_all.sh does: fastq/<sample>_<flowcell>_<lane>_R{1,2}.fq.gz.
set -euo pipefail
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/wgs" || exit 1
mkdir -p cram ref fastq
EBI=https://ftp.sra.ebi.ac.uk/vol1/run/ERR398
REF_URL=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome
REF=ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
JOBS=${JOBS:-12}             # parallel byte ranges
SAMPLES=${SAMPLES:-NA12891 NA12892}
CHUNK=$((256 * 1024 * 1024))

# get URL OUT [MD5]: download in CHUNK-sized ranges, JOBS at a time, then join; skipped when OUT is complete
get() {
  local url=$1 out=$2 md5=${3:-} size n
  if [ -s "$out" ]; then
    [ -z "$md5" ] || [ "$(md5sum < "$out" | cut -d' ' -f1)" = "$md5" ] || { echo "$out: present, md5 differs"; return 1; }
    echo "$out: present${md5:+, md5 ok}"; return
  fi
  size=$(curl -sIL "$url" | grep -i '^content-length' | tail -1 | tr -dc 0-9)
  [ -n "$size" ] || { echo "$out: no size from $url"; return 1; }
  n=$(( (size + CHUNK - 1) / CHUNK ))
  mkdir -p "$out.parts"
  seq 0 $((n - 1)) | xargs -P "$JOBS" -I{} bash -c '
    i={}; s=$((i * '$CHUNK')); e=$((s + '$CHUNK' - 1)); [ $e -ge '$size' ] && e=$(('$size' - 1))
    p='"$out"'.parts/$i; want=$((e - s + 1))
    for try in $(seq 20); do
      [ "$(stat -c %s $p 2>/dev/null || echo 0)" = $want ] && exit 0
      curl -sS -L --retry 5 -r $s-$e -o $p '"$url"' || sleep 30
    done
    echo "chunk $i FAILED"; exit 1'
  for i in $(seq 0 $((n - 1))); do cat "$out.parts/$i"; done > "$out.part"
  [ "$(stat -c %s "$out.part")" = "$size" ] || { echo "$out: size differs"; return 1; }
  if [ -n "$md5" ]; then
    [ "$(md5sum < "$out.part" | cut -d' ' -f1)" = "$md5" ] || { echo "$out: md5 differs"; return 1; }
  fi
  mv "$out.part" "$out"; rm -rf "${out:?}.parts"
  echo "$out: $size bytes${md5:+, md5 ok}"
}

date "+start %T"
get $REF_URL/GRCh38_full_analysis_set_plus_decoy_hla.fa     $REF
get $REF_URL/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai $REF.fai
declare -A RUN=([NA12891]=ERR3989341 [NA12892]=ERR3989342)
declare -A MD5=([NA12891]=965a44721c7a426350cae5baeb696738 [NA12892]=b9750c2b608e468a3e844ae74da0c562)
for s in $SAMPLES; do get $EBI/${RUN[$s]}/$s.final.cram cram/$s.final.cram ${MD5[$s]}; done
date "+downloads done %T"

# one lane: read group ID (NA12891_<barcode>_<flowcell>_<lane>) -> fastq/<sample>_<flowcell>_<lane>_R{1,2}.fq.gz
lane() {
  local id=$1 s fl
  s=$(cut -d_ -f1 <<< "$id"); fl=$(cut -d_ -f3- <<< "$id")
  [ -s fastq/${s}_${fl}_R2.fq.gz ] && { rm -f fastq/split/$id.bam; return; }
  # .part.fq.gz, not .fq.gz.part: samtools fastq gzips only names ending in .gz
  samtools collate -@ 8 -u -O -T fastq/split/$id.tmp fastq/split/$id.bam \
    | samtools fastq -@ 8 -F 0x900 -1 fastq/${s}_${fl}_R1.part.fq.gz -2 fastq/${s}_${fl}_R2.part.fq.gz \
        -0 /dev/null -s /dev/null -n - 2> fastq/${s}_${fl}.log
  mv fastq/${s}_${fl}_R1.part.fq.gz fastq/${s}_${fl}_R1.fq.gz; mv fastq/${s}_${fl}_R2.part.fq.gz fastq/${s}_${fl}_R2.fq.gz
  rm -f fastq/split/$id.bam
  echo "${s}_${fl} done $(date +%T)"
}
export -f lane
mkdir -p fastq/split
for s in $SAMPLES; do
  [ "$(ls fastq/${s}_*_R2.fq.gz 2>/dev/null | wc -l)" = 12 ] && { echo "$s: 12 lanes present"; continue; }
  [ -e fastq/split/$s.done ] || { samtools split -@ 32 --reference $REF --output-fmt BAM,level=1 \
      -f "fastq/split/%!.%." cram/$s.final.cram && touch fastq/split/$s.done; }
  echo "$s split $(date +%T)"
  samtools view -H cram/$s.final.cram | grep '^@RG' | grep -o 'ID:[^[:space:]]*' | cut -c4- | xargs -P 6 -I{} bash -c 'lane {}'
done
rm -f fastq/split/*.done; rmdir fastq/split 2>/dev/null || true
date "+all done %T"
ls fastq/NA1289[12]_*_R1.fq.gz | wc -l
