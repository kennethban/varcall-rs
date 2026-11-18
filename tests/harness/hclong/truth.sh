#!/usr/bin/env bash
# HaplotypeCaller with reads spanning >= 131,072 bp (review finding 17): a 75M140000D75M read and a 200,000 bp
# all-match read added to HG002 chr20:10.0-10.4 Mb. GATK has no span limit.
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/hclong" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
REF=../hg002/ref/GRCh38_no_alt_analysis_set.fa; SRC=../hg002/bam/HG002.bam
if [ ! -s long.bam ]; then
  samtools view -h $SRC chr20:10000000-10400000 > src.sam
  rg=$(grep -m1 '^@RG' src.sam | tr '\t' '\n' | grep '^ID:' | cut -c4-)
  python3 - "$rg" "$REF" >> src.sam <<'PY'
import sys, subprocess
rg, ref = sys.argv[1:]
def seq(c, s, e): return ''.join(subprocess.run(['samtools', 'faidx', ref, f'{c}:{s}-{e}'], capture_output=True, text=True).stdout.split('\n')[1:]).upper()
q = lambda n: 'F' * n
s1 = 10050001; a = seq('chr20', s1, s1 + 74); b = seq('chr20', s1 + 75 + 140000, s1 + 75 + 140000 + 74)
print('\t'.join(['longdel', '0', 'chr20', str(s1), '60', '75M140000D75M', '*', '0', '0', a + b, q(150), 'RG:Z:' + rg]))
s2 = 10120001; c = seq('chr20', s2, s2 + 199999)
print('\t'.join(['longmatch', '16', 'chr20', str(s2), '60', '200000M', '*', '0', '0', c, q(200000), 'RG:Z:' + rg]))
PY
  samtools sort -o long.bam src.sam && samtools index long.bam && rm src.sam
fi
mkdir -p truth
for m in vcf gvcf; do
  e=""; [ $m = gvcf ] && e="-ERC GVCF"
  java -Xmx8g -jar $J HaplotypeCaller -R $REF -I long.bam -L chr20:10000000-10400000 -O truth/out.$m.gz $e > truth/$m.log 2>&1; echo "GATK $m rc=$?"
done
