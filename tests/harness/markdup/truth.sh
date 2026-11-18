#!/usr/bin/env bash
# Picard (GATK 4.7 bundle, Picard 3.5.0) MarkDuplicates reference run on the whole-genome test BAM.
HARNESS=$(cd "$(dirname "$0")/.." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA/markdup" || exit 1

until grep -q done input/make.log; do sleep 30; done
mkdir -p truth/tmp
/usr/bin/time -v java -Xmx12g -jar ${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar} MarkDuplicates -I input/pcrplus.undup.bam \
  -O truth/pcrplus.md.bam -M truth/pcrplus.metrics.txt --TMP_DIR truth/tmp --CREATE_INDEX true > truth/pcrplus.log 2>&1
echo "rc=$?" >> truth/pcrplus.log
