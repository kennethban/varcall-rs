#!/usr/bin/env bash
# FilterMutectCalls with a record whose tumor AD sums to 0 (as a forced allele with no informative reads gives;
# review finding 34): its allele fraction 0/0 is a NaN, which Java's Double.compare sorts last. Input: the
# somatic/longread Mutect2 calls (GATK) with one SNV's AD set to 0,0 (and DP 0).
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/somatic/fmcnan" || exit 1

J="java -Xmx4g -jar ${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}"; REF=$VARCALL_DATA/ref_hg38/Homo_sapiens_assembly38.fasta
S=../longread/truth/out.vcf
awk -F'\t' 'BEGIN{OFS="\t"} /^#/ {print; next} !done && length($4) == 1 && length($5) == 1 {
  n = split($9, k, ":"); m = split($10, v, ":"); for (i = 1; i <= n; i++) { if (k[i] == "AD") v[i] = "0,0"; if (k[i] == "DP") v[i] = "0" }
  s = v[1]; for (i = 2; i <= m; i++) s = s ":" v[i]; $10 = s; done = 1 } {print}' $S > in.vcf
cp $S.stats in.vcf.stats
mkdir -p truth
$J FilterMutectCalls -R $REF -V in.vcf -O truth/out.vcf > truth/fmc.log 2>&1; echo "GATK FilterMutectCalls rc=$?"
