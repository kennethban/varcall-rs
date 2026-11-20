# Validation

"Identical" means byte for byte:

* **BAM:** the records. The header too, apart from the tool's own `@PG` line.
* **VCF:** every data line. The header too, apart from the command-line line.
* **Metrics, recalibration tables, pileup and contamination tables, `.stats`:** the files.
* **Indexes (`.bai`, `.tbi`):** either byte-identical or giving identical answers to region queries.

The reference outputs come from bwa-mem2 2.2.1 (+ samtools 1.25 sort / index) and GATK 4.7.0.0 (Picard 3.5.0,
htsjdk 5.0.0).

## Test tiers

| Tier | Command | Data | Time |
|---|---|---|---|
| Unit tests | `cargo test --release --workspace` | none | seconds |
| Self-contained | `tests/selftest.sh` | synthetic genome and reads made by `tests/synth/simulate.py` (seeded, byte-reproducible); GATK truth committed in `tests/truth` | ~1 min |
| Images | `containers/test_images.sh` | the self-contained test run inside the Docker images (+ real-data checks with `VARCALL_DATA`) | ~5 min |
| Full harness | `tests/harness/regress.sh` | ~300 GB of public data ([tests/harness/README.md](../tests/harness/README.md)) | ~75 min |

The unit tests and the self-contained test run in CI ([.github/workflows/ci.yml](../.github/workflows/ci.yml))
on every push and pull request: CPU only, as GitHub's runners have no GPU.

**Self-contained test.** The synthetic data is a 1.9 Mb genome with these features:

* GC drift and a repeat family;
* STRs, homopolymers, segmental duplications and N gaps;
* germline SNPs, indels, MNPs and multi-allelic sites;
* three somatic clones;
* Illumina-like quality profiles, adapters, and PCR and optical duplicates.

It is aligned by bwa-mem2 and by align-rs (4 lanes; GPU and CPU paths). Then variant-rs runs GATK's command lines
for MarkDuplicates, BQSR, HaplotypeCaller (VCF and GVCF), GenotypeGVCFs (from a GVCF and through a GenomicsDB workspace), Mutect2 and its filters, a panel of normals, VQSR
(VariantRecalibrator, ApplyVQSR) and VariantFiltration / MakeSitesOnlyVcf, with the CPU and the GPU PairHMM. 34
outputs per PairHMM are compared with GATK's. The truth in `tests/truth` comes from GATK's release jar
(`gatk-package-4.7.0.0-local.jar` from `gatk-4.7.0.0.zip` on GitHub, SHA-256 `882e1807…3861365`);
`tests/make_truth.sh` regenerates it.

**Full harness.** These are the comparisons that were run, at full scale, during development. All are in
`regress.sh` or the harness folders:

| Tool | Data | Result |
|---|---|---|
| align-rs | NA12878 30x lanes (63.8 M records each), 20 k / 1 M-pair subsets, synthetic indel and long-read sets, malformed FASTQ, `.alt` files | identical to bwa-mem2 |
| HaplotypeCaller | NA12878 30x WGS (BAM and CRAM; 23 and 97 shards); HG002 35x WGS, HG002 PCR+ 30x WGS, HG002 exome | VCF and GVCF identical (e.g. NA12878: 5,064,948 VCF and 461,802,035 GVCF lines; HG002 PCR+: 583,697,036 GVCF lines) |
| HaplotypeCaller, other arguments | WARP's command lines, `--dbsnp`, `-contamination`, interval padding, several inputs, joint calling of up to 11 samples | identical |
| BaseRecalibrator / ApplyBQSR | HG002 PCR+ 30x WGS, NA12878 chr20 (16 read groups), WARP's `-OQ` / static quantization | table and 721,883,698 records identical |
| MarkDuplicates | HG002 PCR+ 30x WGS, 12-lane inputs, query-grouped (WARP), CRAM, edge cases | records and metrics identical |
| Mutect2 + filters | DREAM synthetic tumors 1-4, HCC1143, SEQC2 WES; tumor-only, PoN, mitochondria mode, `--bam-output` | VCF, stats, filtering stats identical |
| GenotypeGVCFs, ReblockGVCF | HG002 whole-genome GVCF (47 shards), 3 and 11 samples, GenomicsDB workspaces | identical |
| GenomicsDBImport | 3 and 11 samples (raw and reblocked GVCFs), whole-genome GVCF (three arrays, an array starting mid-contig); WARP's sample-map commands (batch 1 and 50, `--consolidate`), sarek's (`-V`, BED intervals, `--merge-input-intervals`, `--bypass-feature-reader`) | every cell, `vidmap.json`, `callset.json` identical (10 workspaces); GenotypeGVCFs (GATK's and variant-rs's) on them gives GATK's output |
| CreateSomaticPanelOfNormals | the four DREAM normals (b37 chr20): Mutect2 tumor-only, GenomicsDBImport (one and three arrays), the panel with five argument sets | normal VCFs, workspace cells, received records and panels identical |
| VariantRecalibrator / ApplyVQSR | NA12878 30x VCF (5.06 M records); 1000 Genomes 30x joint calls chr20:5-17 Mb (50 of 3202 samples, 500 k records); WARP's, sarek's and best-practice (InbreedingCoeff) settings, SNP / INDEL / BOTH, ApplyVQSR chained | tranches, recal VCF and filtered VCF identical (27 outputs) |
| VariantFiltration / MakeSitesOnlyVcf | NA12878 30x VCF, the 1000 Genomes cohort slice, 11-sample joint calls; WARP's ExcessHet, exome-mask and DRAGEN filters, GATK's hard-filtering examples, the expression language's edge cases (missing values, lists, strings, inversion, previous filters, masks) | identical; where GATK fails, variant-rs fails |
| MergeVcfs, SortSam, GatherBamFiles, GatherBQSRReports | WARP's and sarek's command lines | identical |
| Pipeline | NA12878 chr20 slice (2 lanes), HCC1143 tumor / normal, SEQC2 WES; joint calling of NA12878 and HCC1143_normal (GVCFs, GenomicsDBImport, GenotypeGVCFs, VQSR); each step vs GATK | identical |

Accuracy follows from identity. On HG002 against the GIAB v4.2.1 benchmark (chr1–22), variant-rs equals GATK:
precision 0.9887, recall 0.9923, F 0.9905.

## Adversarial review

After the first real-data crash, one that the final-output tests had never reached, the code was reviewed for the
same class of bug. Eight independent reviewers (AI agents) each took one area: four in align-rs, four in variant-rs.
Each read the code against the reference source (bwa-mem2, GATK, Picard, htsjdk) and reported concrete inputs that
would diverge. A separate skeptic then tried to refute each finding.

Every fix came with a regression test that failed before the fix, and those tests are now part of the harness.

| | |
|---|---|
| Findings | 45: 36 confirmed, 3 unverified, 6 refuted |
| Fixed | 37 of the 39 confirmed or unverified, plus 1 refuted finding |
| Not fixed | 2: see below |

What the review found, by kind:

* **Crashes on rare inputs:**
  * reads of 725 bp or more;
  * an empty mate;
  * references under ~500 Mb with `--gpu`;
  * a BAM bucket over 4 GiB;
  * reads spanning over 131 kb (HaplotypeCaller) or 4.8 kb (Mutect2);
  * a short, very deep contig in the shard planner;
  * an empty input in GatherBamFiles;
  * Picard list-file inputs.
* **Silent wrong output:**
  * an assembly failure skipped the region where GATK fails the run;
  * MarkDuplicates and ApplyBQSR on CRAM dropped unplaced reads;
  * quality scores ≥ 94 indexed the wrong table cells;
  * the BAM `bin` field was copied instead of recomputed;
  * `.bgz` outputs were written as plain text;
  * a failed GPU allocation aborted the run instead of falling back to the CPU.
* **Edge-case divergence from the reference:**
  * malformed FASTQ, now read exactly as `kseq_read` does;
  * `.alt` parsing;
  * NaN ordering (`Double.compare`);
  * case-insensitive allele matching;
  * `<NON_REF>` in FilterMutectCalls;
  * contig order without `-L`;
  * htsjdk's QUAL round trip;
  * H tags;
  * records at POS 0;
  * Java's signed bytes for qualities ≥ 128;
  * duplicate sample tables.
* **Pipeline robustness:** `pipefail` on the bwa-mem2 branch; duplicate (sample, lane) rows.

## Known differences and limits

* **bwa-mem2's AVX-512 mate rescue** (`kswv`) counts rows in 16-bit integers, which overflows for rescue windows
  of 32 kb or more. align-rs gives the result of bwa's scalar `ksw_align2`, which bwa-mem2 also uses on CPUs
  without AVX-512. Such windows need very large insert sizes; none occurred in the test data.
* **MarkDuplicates** can pair reads differently from Picard when one read-group + read-name occurs in two pairs
  whose reads fall in different 2 Mb chunks. That is malformed input, and Picard's result for it depends on order
  too.
* **The same read in two HaplotypeCaller inputs** is refused. GATK collapses identical reads in some of its
  hash-keyed maps but not in others.
* Inputs or arguments outside what is listed in [reference.md](reference.md) are refused, or handed to GATK under
  `GATK_FALLBACK`. They never silently give a different result.
* The full harness has run on one GPU model (RTX 5070 Ti, sm_120) and one CPU (16-core Zen 5, with AVX-512). On a
  server with L40S GPUs (sm_89) and EPYC 9654 CPUs, the pipeline's tests, NA12878 30x from FASTQ on 1 and 6 GPUs, and
  the SEQC2 exome pair gave the same outputs as on the workstation. The CPU PairHMM without AVX-512 is checked by the self-contained test (`VARIANT_RS_NO_AVX512=1`, and CI runners
  without AVX-512).
