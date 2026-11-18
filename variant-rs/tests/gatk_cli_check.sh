#!/usr/bin/env bash
# variant-rs with GATK's own command lines (as in hg002/run_gatk.sh, bqsr/truth.sh, markdup/truth.sh,
# pipeline/germline.sh, somatic/orient/truth.sh, somatic/contam/truth.sh, WARP's BamProcessing.wdl via
# bqsr/cmp_warp.sh and markdup/cmp_warp.sh, WARP's GermlineVariantDiscovery.wdl via hg002/group2/cmp.sh, Mutect2 options via somatic/group3/cmp.sh) vs the GATK 4.7 outputs; then GATK_FALLBACK with a stub `gatk` that records its arguments.
# usage: gatk_cli_check.sh [OUTDIR]

REPO=$(cd "$(dirname "$0")/../.." && pwd); P=${VARCALL_DATA:?set VARCALL_DATA to the test data directory}; G=${VARIANT_RS:-$REPO/target/release/variant-rs}
O=${1:-$P/variant-rs/tests/out}; rm -rf $O; mkdir -p $O/tmp
same() { [ "$1" = "$2" ] && echo "identical" || echo "DIFFERENT"; }

echo "== HaplotypeCaller (hg002/run_gatk.sh command line, shard023)"
cd $P/hg002
for m in vcf gvcf; do
  e=""; ext=vcf.gz; [ $m = gvcf ] && { e="-ERC GVCF"; ext=g.vcf.gz; }
  /usr/bin/time -f "%e s" $G HaplotypeCaller -R ref/GRCh38_no_alt_analysis_set.fa -I bam/HG002.bam -L plan/shard023.intervals \
    -O $O/shard023.$ext $e --native-pair-hmm-threads 2 > $O/hc.$m.log 2>&1
  echo "  $m: $($REPO/tests/harness/hg002/cmp1.sh $O/shard023.$ext gatk/$m $ext) ($(tail -1 $O/hc.$m.log))"
done

echo "== HaplotypeCaller with WARP's arguments and S3 group 2 (hg002/group2/truth.sh: -G, -GQB, -contamination, --dbsnp, -ip, -XL)"
$REPO/tests/harness/hg002/group2/cmp.sh $O/hc_group2

echo "== HaplotypeCaller with several -I (S3 group 4: hg002/group4/truth.sh, inputs from make_inputs.sh)"
$REPO/tests/harness/hg002/group4/cmp.sh $O/hc_group4

echo "== HaplotypeCaller, several samples called jointly (S3 group 4: hg002/group4/ms/make.sh; HG002, NA12878, HG002 PCR+)"
$REPO/tests/harness/hg002/group4/ms/cmp.sh $O/hc_ms

echo "== BaseRecalibrator + ApplyBQSR (bqsr/truth.sh command lines, NA12878 chr20)"
cd $P/bqsr
REF=ref/GRCh38_full_analysis_set_plus_decoy_hla.fa; R=resources; BAM=$P/wgs/bam/NA12878.bam
KS="--known-sites $R/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites $R/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --known-sites $R/Homo_sapiens_assembly38.known_indels.vcf.gz"
$G BaseRecalibrator -R $REF -I $BAM -L chr20 $KS -O $O/na.recal.table > $O/br.log 2>&1
echo "  recal table: $(cmp -s $O/na.recal.table truth/na12878.chr20.recal.table && echo identical || echo DIFFERENT)"
$G ApplyBQSR -R $REF -I $BAM -L chr20 --bqsr-recal-file truth/na12878.chr20.recal.table -O $O/na.bqsr.bam > $O/ab.log 2>&1
echo "  recalibrated BAM records: $(python3 $REPO/tests/harness/bqsr/cmpbam.py truth/na12878.chr20.bqsr.bam $O/na.bqsr.bam 2>&1 | tail -1)"

echo "== MarkDuplicates (markdup/truth.sh and germline.sh command lines, smoke BAM)"
cd $P/markdup
$G MarkDuplicates -I input/smoke.bam -O $O/smoke.md.bam -M $O/smoke.metrics.txt --TMP_DIR $O/tmp --CREATE_INDEX true > $O/md.log 2>&1
echo "  records: $(python3 $REPO/tests/harness/bqsr/cmpbam.py truth/smoke.md.bam $O/smoke.md.bam 2>&1 | tail -1); metrics: $(diff -q <(grep -v '^#' truth/smoke.metrics.txt) <(grep -v '^#' $O/smoke.metrics.txt) > /dev/null && echo identical || echo DIFFERENT)"
$G MarkDuplicates I=input/smoke.bam O=$O/smoke2.md.bam M=$O/smoke2.metrics.txt TMP_DIR=$O/tmp > $O/md2.log 2>&1
echo "  Picard NAME=value syntax: records $(python3 $REPO/tests/harness/bqsr/cmpbam.py truth/smoke.md.bam $O/smoke2.md.bam 2>&1 | tail -1)"

echo "== WARP arguments (tasks/wdl/BamProcessing.wdl): bqsr/cmp_warp.sh, markdup/cmp_warp.sh"
$REPO/tests/harness/bqsr/cmp_warp.sh $O/bqsr_warp
$REPO/tests/harness/markdup/cmp_warp.sh $O/markdup_warp
echo "== MarkDuplicates over several coordinate-sorted inputs, as nf-core/sarek runs it (markdup/multi/make.sh)"
$REPO/tests/harness/markdup/multi/cmp.sh $O/markdup_multi

echo "== GatherBQSRReports (S3 group 4: bqsr/gather/truth.sh, BamProcessing.wdl's command on scattered BaseRecalibrator tables)"
$REPO/tests/harness/bqsr/gather/cmp.sh $O/gather

echo "== GatherBamFiles (S3 group 4: gatherbam/truth.sh, BamProcessing.wdl's GatherSortedBamFiles / GatherUnsortedBamFiles on scattered ApplyBQSR outputs)"
$REPO/tests/harness/gatherbam/cmp.sh $O/gatherbam

echo "== SortSam (S3 group 4: sortsam/truth.sh, BamProcessing.wdl's command on the WARP MarkDuplicates outputs)"
$REPO/tests/harness/sortsam/cmp.sh $O/sortsam

echo "== MergeVcfs (S3 group 4: mergevcfs/truth.sh, GermlineVariantDiscovery.wdl's MergeVCFs on scattered HaplotypeCaller outputs)"
$REPO/tests/harness/mergevcfs/cmp.sh $O/mergevcfs

echo "== ReblockGVCF (S3 group 4: reblock/truth.sh, GermlineVariantDiscovery.wdl's Reblock on HaplotypeCaller GVCFs; whole genome: reblock/wgs_cmp.sh)"
$REPO/tests/harness/reblock/cmp.sh $O/reblock

echo "== GenotypeGVCFs (S3 group 4: genotypegvcfs/truth.sh single sample, multi/make.sh 3 samples, multi/make11.sh 11 samples; JointGenotypingTasks.wdl's arguments on GVCF files)"
$REPO/tests/harness/genotypegvcfs/cmp.sh $O/genotypegvcfs

echo "== GenotypeGVCFs -V gendb:// (S3 group 4: genomicsdb/truth.sh + truth2.sh, workspaces from GATK GenomicsDBImport with WARP's import command)"
$REPO/tests/harness/genomicsdb/cmp.sh $O/genomicsdb
$REPO/tests/harness/genomicsdb/truth_nol.sh cmp $O/gendb_nol   # without -L, as nf-core/sarek runs it

echo "== Mutect2 family (somatic/orient/truth.sh and somatic/contam/truth.sh command lines, DREAM pairs 1-4)"
cd $P/somatic; D=dream/data; R=$D/human_g1k_v37.20.21.fasta; OT=orient/truth; CT=contam/truth
vcfsame() { diff -q <(grep -v CommandLine $1) <(grep -v CommandLine $2) > /dev/null && echo identical || echo DIFFERENT; }
for x in 1 2 3 4; do
  N=$(samtools view -H $D/normal_$x.bam | grep -m1 '^@RG' | tr '\t' '\n' | grep '^SM:' | cut -c4-)
  $G Mutect2 --java-options "-Xmx16g" -R $R -I $D/tumor_$x.bam -I $D/normal_$x.bam -normal "$N" -L 20 --germline-resource $D/very-small-gnomad.vcf \
    --f1r2-tar-gz $O/f1r2_$x.tar.gz -O $O/m2_$x.vcf > $O/m2_$x.log 2>&1
  rm -rf $O/t $O/r; mkdir -p $O/t $O/r; tar -xzf $OT/f1r2_$x.tar.gz -C $O/t; tar -xzf $O/f1r2_$x.tar.gz -C $O/r
  f1=identical; for f in $O/t/*.alt_table; do cmp -s $f $O/r/$(basename $f) || f1=DIFFERENT; done
  for e in ref_histogram alt_histogram; do f=$(ls $O/t/*.$e); python3 $REPO/tests/harness/somatic/orient/cmp_hist.py $f $O/r/$(basename $f) | grep -q '^identical' || f1=DIFFERENT; done
  $G LearnReadOrientationModel -I $OT/f1r2_$x.tar.gz -O $O/priors_$x.tar.gz > $O/lrom_$x.log 2>&1
  rm -rf $O/t $O/r; mkdir -p $O/t $O/r; tar -xzf $OT/priors_$x.tar.gz -C $O/t; tar -xzf $O/priors_$x.tar.gz -C $O/r
  pr=$(diff -rq $O/t $O/r > /dev/null && echo identical || echo DIFFERENT)
  $G FilterMutectCalls -R $R -V $OT/m2_$x.vcf --ob-priors $OT/priors_$x.tar.gz -O $O/fo_$x.vcf > $O/fo_$x.log 2>&1
  fo="$(vcfsame $OT/fo_$x.vcf $O/fo_$x.vcf) (stats $(cmp -s $OT/fo_$x.vcf.filteringStats.tsv $O/fo_$x.vcf.filteringStats.tsv && echo identical || echo DIFFERENT))"
  gps=identical
  for s in tumor normal; do
    $G GetPileupSummaries -I $D/${s}_$x.bam -V contam/res/exac_20.b37.vcf.gz -L contam/res/exac_20.b37.vcf.gz -O $O/${s}_$x.pileups.table > $O/gps_$x.log 2>&1
    cmp -s $CT/${s}_$x.pileups.table $O/${s}_$x.pileups.table || gps=DIFFERENT
  done
  $G CalculateContamination -I $CT/tumor_$x.pileups.table -matched $CT/normal_$x.pileups.table -O $O/dream_$x.contamination.table \
    --tumor-segmentation $O/dream_$x.segments.table > $O/cc_$x.log 2>&1
  cc=$(cmp -s $CT/dream_$x.contamination.table $O/dream_$x.contamination.table && cmp -s $CT/dream_$x.segments.table $O/dream_$x.segments.table && echo identical || echo DIFFERENT)
  $G FilterMutectCalls -R $R -V dream/truth/m2_$x.vcf --contamination-table $CT/dream_$x.contamination.table \
    --tumor-segmentation $CT/dream_$x.segments.table -O $O/fc_$x.vcf > $O/fc_$x.log 2>&1
  echo "  pair $x: Mutect2 VCF $(vcfsame $OT/m2_$x.vcf $O/m2_$x.vcf), stats $(cmp -s $OT/m2_$x.vcf.stats $O/m2_$x.vcf.stats && echo identical || echo DIFFERENT), F1R2 counts $f1;"\
    "LearnReadOrientationModel $pr; FilterMutectCalls --ob-priors $fo; GetPileupSummaries $gps; CalculateContamination $cc;"\
    "FilterMutectCalls --contamination-table $(vcfsame $CT/fc_$x.vcf $O/fc_$x.vcf)"
done

echo "== GetPileupSummaries over BED / interval-list / contig:start-end pieces, as nf-core/sarek scatters it (somatic/contam/sarek.sh)"
$REPO/tests/harness/somatic/contam/sarek.sh cmp $O/gps_sarek
echo "== MergeMutectStats + GatherPileupSummaries, sarek's gathers after a scattered Mutect2 (somatic/gathers/run.sh)"
$REPO/tests/harness/somatic/gathers/run.sh cmp $O/m2_gathers

echo "== Mutect2 S3 group 3 (somatic/group3/truth.sh: tumor-only, --genotype-germline-sites / --genotype-pon-sites, --af-of-alleles-not-in-resource, --mitochondria-mode incl. WARP's mito command, --bam-output)"
$REPO/tests/harness/somatic/group3/cmp.sh $O/m2_group3
echo "== FilterMutectCalls --mitochondria-mode, --min-allele-fraction, --max-alt-allele-count, --f-score-beta (somatic/group3/fmc_truth.sh, incl. AlignAndCall.wdl's two Filter calls)"
$REPO/tests/harness/somatic/group3/fmc_cmp.sh $O/fmc_group3

echo "== refusals and GATK_FALLBACK"
cat > $O/gatk <<'STUB'
#!/usr/bin/env bash
echo "stub gatk called with: $*"
STUB
chmod +x $O/gatk
cd $P/hg002
r=$($G HaplotypeCaller -R ref/GRCh38_no_alt_analysis_set.fa -I bam/HG002.bam -O $O/x.vcf -ploidy 1 2>&1); echo "  refused (rc=$?): $(echo "$r" | sed -n 2p)"
echo "  fallback, refused argument: $(GATK_FALLBACK=$O/gatk $G HaplotypeCaller -R r.fa -I x.bam -O x.vcf -ploidy 1 2>/dev/null)"
echo "  fallback, unsupported tool: $(GATK_FALLBACK=$O/gatk $G VariantRecalibrator -R r.fa -V x.vcf.gz -O y.recal 2>/dev/null)"
echo "  fallback, launcher option first (nf-core form): $(GATK_FALLBACK=$O/gatk $G --java-options "-Xmx2g -XX:-UsePerfData" VariantRecalibrator -V x.vcf.gz -O y.recal 2>/dev/null)"
r=$($G VariantRecalibrator -V x 2>&1); echo "  unsupported tool without fallback (rc=$?): $(echo "$r" | head -1 | cut -c1-80)"
