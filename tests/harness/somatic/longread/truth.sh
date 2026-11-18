#!/usr/bin/env bash
# Mutect2 (tumor-only) with long reads (review finding 19 + the ring limit): HCC1143 tumor chr17:7.5-8.5 Mb plus
# 75M20000D75M reads that cross interval gaps (variant-rs splits contigs at gaps >= 5 kb), 75M8000D75M reads, and a
# 200,000 bp all-match read; intervals of 1 kb every 10 kb.
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/longread" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
H=../hcc1143/data; REF=$VARCALL_DATA/ref_hg38/Homo_sapiens_assembly38.fasta
if [ ! -s tumor.bam ]; then
  samtools view -h $H/tumor.bam chr17:7500000-8500000 > src.sam
  rg=$(grep -m1 '^@RG' src.sam | tr '\t' '\n' | grep '^ID:' | cut -c4-)
  python3 - "$rg" "$REF" >> src.sam <<'PY'
import sys, subprocess
rg, ref = sys.argv[1:]
def seq(s, e): return ''.join(subprocess.run(['samtools', 'faidx', ref, f'chr17:{s}-{e}'], capture_output=True, text=True).stdout.split('\n')[1:]).upper()
def rec(name, flag, pos, cigar, s): print('\t'.join([name, str(flag), 'chr17', str(pos), '60', cigar, '*', '0', '0', s, 'F' * len(s), 'RG:Z:' + rg]))
for k in range(10):
    p = 7500501 + k * 100000          # inside an interval; the deletion crosses two gaps
    rec(f'cross{k}', 16 * (k % 2), p, '75M20000D75M', seq(p, p + 74) + seq(p + 20075, p + 20149))
    p = 7530401 + k * 100000          # 8 kb deletion
    rec(f'del8k{k}', 16 * (k % 2), p, '75M8000D75M', seq(p, p + 74) + seq(p + 8075, p + 8149))
p = 7700001
rec('longmatch', 0, p, '200000M', seq(p, p + 199999))
PY
  samtools sort -o tumor.bam src.sam && samtools index tumor.bam && rm src.sam
  python3 -c "
for k in range(100): print(f'chr17\t{7500000 + k * 10000}\t{7501000 + k * 10000}')" > ivs.bed
fi
mkdir -p truth
java -Xmx8g -jar $J Mutect2 -R $REF -I tumor.bam -L ivs.bed -O truth/out.vcf > truth/m2.log 2>&1; echo "GATK rc=$?"
