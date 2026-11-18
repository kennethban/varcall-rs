#!/usr/bin/env bash
# 11-sample GenotypeGVCFs truth (InbreedingCoeff / AS_InbreedingCoeff need >= 10 samples): the three GVCFs of
# make.sh (and their reblocked versions) copied under 11 sample names, chr8:1-5M, CombineGVCFs, GenotypeGVCFs.
HARNESS=$(cd "$(dirname "$0")/../.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/genotypegvcfs/multi" || exit 1

J=${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar}
REF=$VARCALL_DATA/bqsr/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
DBSNP=$VARCALL_DATA/bqsr/resources/Homo_sapiens_assembly38.dbsnp138.vcf.gz
mkdir -p in11 truth
copy() { # src sample newname suffix
  local o=in11/$3$4.g.vcf.gz
  [ -f $o.tbi ] || { tabix -h in/$1$4.g.vcf.gz chr8:1-5000000 | sed "/^#CHROM/s/\t$2\$/\t$3/" | bgzip > $o && tabix -p vcf $o; }
}
for x in "" .rb; do
  for k in a b c d; do copy HG002 HG002 HG002_$k $x; copy NA12878 NA12878 NA12878_$k $x; done
  for k in a b c; do copy HG002pcr HG002pcr HG002pcr_$k $x; done
done
for x in "" .rb; do
  n=combined11${x//./_}
  [ -f in11/$n.g.vcf.gz.tbi ] || java -Xmx6g -jar $J CombineGVCFs -R $REF $(ls in11/*_[a-d]$x.g.vcf.gz | grep -v "\.rb\." | sed 's/^/-V /' | tr '\n' ' ' | ( [ -z "$x" ] && cat || true )) $( [ -n "$x" ] && ls in11/*_[a-d].rb.g.vcf.gz | sed 's/^/-V /' | tr '\n' ' ' ) -G StandardAnnotation -G AS_StandardAnnotation -L chr8:1-5000000 -O in11/$n.g.vcf.gz > in11/$n.log 2>&1 &
done
wait
W="-D $DBSNP -G StandardAnnotation -G AS_StandardAnnotation --only-output-calls-starting-in-intervals -L chr8:1-5000000 --merge-input-intervals"
declare -A C
C[m11_warp]="-V in11/combined11.g.vcf.gz $W"
C[m11_warp_rb]="-V in11/combined11_rb.g.vcf.gz $W"
C[m11_default]="-V in11/combined11.g.vcf.gz"
for c in "${!C[@]}"; do
  echo "${C[$c]}" > truth/$c.cmd
  [ -f truth/$c.vcf.gz.tbi ] || java -Xmx6g -jar $J GenotypeGVCFs -R $REF ${C[$c]} -O truth/$c.vcf.gz > truth/$c.log 2>&1 &
done
wait
for c in "${!C[@]}"; do echo "$c $(zcat truth/$c.vcf.gz | grep -vc '^#') $(zcat truth/$c.vcf.gz | grep -v '^#' | grep -c InbreedingCoeff)"; done
