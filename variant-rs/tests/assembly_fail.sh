#!/usr/bin/env bash
# A region whose assembly fails must fail the run, as in GATK (AssemblyBasedCallerUtils.assembleReads rethrows):
# non-zero exit and no finalized output (no .tbi). REVIEW_FINDINGS.md: the port used to skip the region and exit 0.
# Uses the test hook VARIANT_RS_TEST_ASSEMBLY_FAIL=POS (assembly fails for regions starting at or after POS).
REPO=$(cd "$(dirname "$0")/../.." && pwd); G=${VARCALL_DATA:?set VARCALL_DATA to the test data directory}; O=${1:-/tmp/hc-assembly-fail}; rm -rf $O; mkdir -p $O
VARIANT_RS_TEST_ASSEMBLY_FAIL=10200000 ${VARIANT_RS:-$REPO/target/release/variant-rs} -I $G/data/NA12878.chr20_10-15M.bam -R $G/data/chr20.fa \
  -L chr20:10000000-10500000 -O $O/out.vcf.gz --pairhmm cpu --threads 8 > $O/log 2>&1; rc=$?
grep -m1 -E "assembly failed|test assembly failure" $O/log
if [ $rc -ne 0 ] && [ ! -e $O/out.vcf.gz.tbi ]; then echo "assembly failure: run failed (exit $rc), no index written: OK"
else echo "assembly failure: exit $rc, index $( [ -e $O/out.vcf.gz.tbi ] && echo written || echo absent) (FAILED: GATK fails the run)"; fi
