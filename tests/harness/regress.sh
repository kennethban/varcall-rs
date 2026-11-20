#!/usr/bin/env bash
# The full regression harness: align-rs vs bwa-mem2 and variant-rs vs GATK 4.7 / Picard on real data (NA12878,
# HG002, HCC1143, DREAM, SEQC2; tests/harness/README.md), ~75 min. Each check prints its own result lines; the summary
# at the end lists every line reporting a difference.
#   VARCALL_DATA=/data GATK_JAR=/path/gatk-package-4.7.0.0-local.jar tests/harness/regress.sh [OUTDIR]
# Run in the environment with samtools, bcftools, java, bwa-mem2 and nextflow (environment.yml), after
# cargo build --release. The pipeline checks rerun pipeline/main.nf (alignments from Nextflow's cache when unchanged).
HARNESS=$(cd "$(dirname "$0")/." && pwd); REPO=$(cd "$HARNESS/../.." && pwd)
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory (tests/harness/README.md)}"
cd "$VARCALL_DATA" || exit 1; R=$PWD
O=${1:-/tmp/varcall-regress}; mkdir -p $O
export VARCALL_DATA=$R

run() { local name=$1; shift; local t=$(date +%s); echo "=== $name"; ( "$@" ) > $O/$name.log 2>&1; local rc=$?
  echo "rc=$rc $(( $(date +%s) - t ))s"; grep -i -E "identical|differ|mismatch|same|rc=|records" $O/$name.log | tail -12; }

run aln_rs           bash $REPO/align-rs/tests/cmp.sh $O/aln_rs
run jmath            bash $REPO/variant-rs/tests/jmath/run.sh 10
run hc_group2        bash $HARNESS/hg002/group2/cmp.sh $O/hc_group2
run hc_group4        bash $HARNESS/hg002/group4/cmp.sh $O/hc_group4
run hc_multisample   bash $HARNESS/hg002/group4/ms/cmp.sh $O/hc_ms
run genotypegvcfs    bash $HARNESS/genotypegvcfs/cmp.sh $O/gg
run reblock          bash $HARNESS/reblock/cmp.sh $O/rb
run genomicsdb       bash $HARNESS/genomicsdb/cmp.sh $O/gdb
run genomicsdb_import bash $HARNESS/genomicsdb/import_cmp.sh $O/gdbimport
run markdup_multi    bash $HARNESS/markdup/multi/cmp.sh $O/md_multi
run markdup_warp     bash $HARNESS/markdup/cmp_warp.sh
run markdup_edge     bash $HARNESS/markdup/edge/run.sh
run bqsr_gather      bash $HARNESS/bqsr/gather/cmp.sh
run bqsr_warp        bash $HARNESS/bqsr/cmp_warp.sh
run m2_dream         bash $HARNESS/somatic/dream/cmp_vcf.sh
run fmc_dream        bash $HARNESS/somatic/dream/cmp_fmc.sh
run fmc_contam       bash $HARNESS/somatic/contam/cmp_fmc.sh
run fmc_orient       bash $HARNESS/somatic/orient/cmp_fmc.sh
run m2_group3        bash $HARNESS/somatic/group3/cmp.sh $O/g3
run fmc_group3       bash $HARNESS/somatic/group3/fmc_cmp.sh $O/g3fmc
run hcc1143          bash -c "bash $HARNESS/somatic/hcc1143/run.sh rust && bash $HARNESS/somatic/hcc1143/cmp.sh"
run gatk_cli         bash $REPO/variant-rs/tests/gatk_cli_check.sh $O/cli
run scatter_intervals bash $REPO/variant-rs/tests/scatter_intervals.sh $O/scatter
# adversarial review fixes (REVIEW_FINDINGS.md): each failed before its fix
run review_fixes     bash $REPO/variant-rs/tests/review_fixes.sh $O/review
run hc_assembly_fail bash $REPO/variant-rs/tests/assembly_fail.sh $O/asmfail
run markdup_cram     bash $HARNESS/markdup/cram/run.sh
run bqsr_qualhi      bash $HARNESS/bqsr/qualhi/cmp.sh $O/qualhi
run hc_contig_order  bash $HARNESS/hcorder/cmp.sh $O/hcorder
run hc_long_reads    bash $HARNESS/hclong/cmp.sh $O/hclong
run m2_long_reads    bash $HARNESS/somatic/longread/cmp.sh $O/m2long
run fmc_gvcf         bash $HARNESS/somatic/fmcgvcf/cmp.sh $O/fmcgvcf
run aln_badfq        bash $REPO/align-rs/tests/badfq.sh ${ALIGN_RS:-$REPO/target/release/align-rs} $O/badfq
run aln_gpufail      bash $REPO/align-rs/tests/gpufail.sh ${ALIGN_RS:-$REPO/target/release/align-rs} $O/gpufail
run pipeline_edge    bash $REPO/pipeline/test/edge.sh
run bqsr_pos0        bash $HARNESS/bqsr/pos0/cmp.sh $O/bqsr_pos0
run fmc_nan          bash $HARNESS/somatic/fmcnan/cmp.sh $O/fmcnan
run fmc_dup_tables   bash $HARNESS/somatic/fmcdup/run.sh $O/fmcdup
run mergevcfs_qual   bash $HARNESS/mergevcfs/qual/run.sh $O/mvqual
run vqsr             bash $HARNESS/vqsr/cmp.sh $O/vqsr
run vfilter          bash $HARNESS/vfilter/cmp.sh $O/vfilter
run plan_edge        bash $HARNESS/plan_edge/run.sh $O/plan_edge
run aln_alt          bash $REPO/align-rs/tests/alt.sh ${ALIGN_RS:-$REPO/target/release/align-rs} $O/aln_alt
# sarek whole-genome GatherBQSRReports (356 tables, both runs' input orders) vs GATK's own gather
gather_wgs() {
  for v in stock gatkrs; do
    d=$(grep -l GatherBQSRReports sarek/e2e/work/wgs_germline_$v/*/*/.command.sh | head -1 | xargs dirname)
    ins=$(grep -o -- '--input [^ ]*' $d/.command.sh | sed "s#--input #--input $d/#" | tr '\n' ' ')
    ${VARIANT_RS:-$REPO/target/release/variant-rs} GatherBQSRReports $ins --output $O/gather_$v.table > /dev/null 2>&1
    java -Xmx4g -jar ${GATK_JAR:?set GATK_JAR to the GATK 4.7.0.0 jar} GatherBQSRReports $ins --output $O/gather_$v.gatk.table > /dev/null 2>&1
    cmp -s $O/gather_$v.table $O/gather_$v.gatk.table && echo "  $v order: identical" || echo "  $v order: DIFFERENT"
  done
}
[ -d sarek/e2e/work/wgs_germline_stock ] && run sarek_wgs_gather gather_wgs
# the Nextflow pipeline from FASTQ vs GATK (truth from test/truth.sh). Resumes the test's last successful run, so
# unchanged steps (the alignments) come from the cache; the variant-rs steps rerun when the binary was rebuilt
# (main.nf puts its build time in their commands). REGRESS_FRESH=1 starts from scratch.
pipe() {
  cd $REPO/pipeline; export NXF_ANSI_LOG=false; bash test/render.sh || exit 1
  local sid=
  [ -z "$REGRESS_FRESH" ] && sid=$(awk -F'\t' -v p="-params-file test/$1.yaml " '$4 == "OK" && index($7, p) { s = $6 } END { print s }' .nextflow/history 2>/dev/null)
  [ -z "$sid" ] && rm -rf test/work/$1
  rm -rf test/out/$1
  ${NEXTFLOW:-nextflow} run main.nf -params-file test/$1.yaml -w test/work/$1 -with-trace test/log/$1.trace.tsv ${sid:+-resume $sid} | tail -2
  echo "steps run (not from the cache): $(grep -c COMPLETED test/log/$1.trace.tsv), cached: $(grep -c CACHED test/log/$1.trace.tsv)"
  python3 test/compare.py $1
}
run pipeline_germline pipe germline
run pipeline_germline_bam pipe germline_bam
run pipeline_somatic  pipe somatic
run pipeline_somatic_bam pipe somatic_bam

echo; echo "=== lines reporting a difference (empty: all identical)"
grep -H -E "DIFFER|DIFFERENT|FAILED|panicked|Traceback|No such file|command not found|unbound variable|parameter null or not set|mismatches +[1-9]|[1-9][0-9]* different or missing|[1-9][0-9]* case\(s\) differ|[1-9][0-9]* columns differ|records, [1-9][0-9]* differ|\| [a-z ]*DIFF" $O/*.log | head -40
