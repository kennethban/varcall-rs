# Full regression harness

`tests/selftest.sh` runs anywhere in a minute on synthetic data. This harness is the heavier half: align-rs vs
bwa-mem2 and variant-rs vs GATK 4.7 / Picard on public sequencing data, 50 checks in ~75 minutes on a 16-core
machine with a GPU, ~300 GB of data. Every check compares outputs byte for byte (VCF data lines, BAM records,
metrics, tables), so a change that alters any result shows up here.

```bash
export VARCALL_DATA=/data/varcall            # the data directory (layout below)
export GATK_JAR=/opt/gatk-4.7.0.0/gatk-package-4.7.0.0-local.jar
cargo build --release                        # in the repository; CUDA libraries optional
tests/harness/regress.sh /tmp/regress        # prints each check; the summary lists every difference
```

Run it in an environment with `samtools`, `bcftools`, `bgzip`/`tabix`, Java 17, `bwa-mem2`, `nextflow` and Python 3
(`environment.yml`). Optional: `NEXTFLOW` (the launcher, if not `nextflow` on `PATH`), `VARIANT_RS` / `ALIGN_RS`
(binaries other than `target/release/*`).

## How it is organised

Each folder here matches a folder of `$VARCALL_DATA` holding that area's inputs, GATK truth and outputs. A script
`cd`s into its data folder, so its relative paths (`truth/`, `in/`, `../ref/`) refer to the data, while helpers it
calls come from this folder (`$HARNESS/...`).

| Script | What it does |
|---|---|
| `*/download.sh` | fetches public data (resumable) |
| `*/make*.sh`, `*/make.py` | derives inputs from the downloads (slices, scattered pieces, edge-case BAMs) |
| `*/truth*.sh` | runs GATK / Picard with the command lines the check uses; writes `truth/` in the data folder |
| `*/cmp*.sh`, `*/run.sh` | runs variant-rs the same way and compares with `truth/` |
| `regress.sh` | runs every check (plus the crates' own data tests in `align-rs/tests`, `variant-rs/tests` and the pipeline tests) |
| `cmdargs.sh` | reads / writes a truth case's arguments (`truth/<case>.cmd`) |

A truth case's arguments are kept in `truth/<case>.cmd`, the same for GATK and variant-rs. Paths into the data
directory are written there as `$VARCALL_DATA/...` (`cmdargs.sh -w` when writing, `cmdargs.sh FILE` expands them
when reading), so the truth folders work wherever the data directory is.

To rebuild a truth folder: run the area's `truth*.sh` (needs `GATK_JAR`). Truth is GATK 4.7.0.0 output, so it only
changes when the GATK version does.

Most truth folders were made during development with GATK 4.7.0.0 built from its source tag, plus dump hooks that
are off unless their environment variables are set. On the inputs checked, it gives the same output as the
release jar: the self-contained test's 18 truth files and HaplotypeCaller on HG002 chr8:1-15M (48,066 records) are
byte-identical. Use the release jar (`gatk-package-4.7.0.0-local.jar` from `gatk-4.7.0.0.zip`) for new truth.

## Data (`$VARCALL_DATA`)

All public. Sizes are what the harness keeps, outputs included.

| Folder | Contents | Source |
|---|---|---|
| `ref_hg38/` (20 GB) | GATK GRCh38 `Homo_sapiens_assembly38` + bwa-mem2 index | `ref_hg38/download.sh` (AWS iGenomes) |
| `wgs/` (120 GB) | NA12878 30x WGS (1000 Genomes / NYGC CRAM and a BAM made from it), its reference, GIAB HG001 benchmark, FASTQ lanes | `wgs/download.sh` (1000 Genomes S3, GIAB); `wgs/bam/NA12878.bam`: the CRAM as BAM; `wgs/fastq/NA12878_<flowcell>_<lane>_R{1,2}.fq.gz`: `wgs/fastq/extract_all.sh` (split by read group, `samtools collate` + `fastq`); `wgs/fastq/sub/` = the first pairs of lane HCLHLDSXX_L001 (tiny 20k, small 1M, mid) plus `align-rs/tests/synth.py` sets, each with `*.bwamem2.sam` truth; NA12878's parents NA12891 and NA12892 (M3's trio, 12 lanes each, 111 GB of FASTQ): `wgs/trio.sh` (1000 Genomes 30x CRAMs from EBI, split as `extract_all.sh`) |
| `data/` | NA12878 chr20:10-15 Mb slice and `chr20.fa` | slice of `wgs/bam` |
| `hg002/` (55 GB) | HG002 NovaSeq PCR-free 35x BAM (DeepVariant case study), GRCh38 no-alt reference, GIAB v4.2.1 benchmark, 96 shard plans, GATK truth per shard | `hg002/download.sh`; `hg002/run_gatk.sh` (truth); `group2/`, `group4/` (`make_inputs.sh`, `ms/make.sh`) |
| `bqsr/` (7 GB) | GRCh38 + decoys + HLA reference, Broad hg38 known sites (dbSNP 138, Mills, known indels), BQSR truth | reference from `wgs/download.sh`; known sites from the Broad hg38 bundle (`gs://gcp-public-data--broad-references/hg38/v0`); `truth.sh`, `truth_warp.sh`, `gather/`, `qualhi/`, `pos0/` |
| `markdup/` (72 GB) | MarkDuplicates inputs (smoke, WARP-style, several lanes, edge cases, CRAM) and Picard truth | `truth.sh`, `truth_warp.sh`, `multi/make.sh`, `edge/run.sh` (makes its inputs), `cram/run.sh` |
| `gatherbam/`, `sortsam/`, `mergevcfs/`, `reblock/`, `genotypegvcfs/`, `genomicsdb/` (9 GB) | scattered outputs of the steps above, Picard / GATK truth | each folder's `truth*.sh` (and `genotypegvcfs/multi/make*.sh`; `genomicsdb/import_truth.sh`: GenomicsDBImport truth; `genomicsdb/scatter_truth.sh`: GenotypeGVCFs once per interval) |
| `m3` (run folder, ~460 GB; outside `$VARCALL_DATA`) | the trio's pipeline run and GATK truth | `m3/run.sh RUN_DIR` (the pipeline), `m3/truth.sh RUN_DIR` (GATK, each step from the pipeline's input; ~3.5 h, set by GATK's single-threaded BaseRecalibrator), `m3/cmp.sh RUN_DIR` |
| `vqsr/` (6 GB) | HapMap / Omni / 1000G / Axiom resources, NA12878 30x HaplotypeCaller VCF (`wgs/truth` shards joined), 1000 Genomes 30x joint calls chr20:5-17 Mb (50 samples), VQSR truth | `vqsr/download.sh` (Broad hg38 bundle, 1000 Genomes FTP), `vqsr/truth.sh` |
| `vfilter/` (2 GB) | exome target masks, edge-case VCFs, VariantFiltration truth (inputs from `vqsr/`, `genotypegvcfs/`, `somatic/seqc2/`) | `vfilter/truth.sh` (makes its inputs) |
| `somatic/dream/` (2 GB) | ICGC-TCGA DREAM synthetic tumor/normal pairs 1-4 (chr20/21, b37) as in GATK's own tests | `dream/fetch_lfs.py` (GATK repository Git LFS: `src/test/resources/large/mutect`) |
| `somatic/hcc1143/` (1 GB) | HCC1143 tumor/normal (GATK somatic workshop, hg38 chr17), PoN, gnomAD AF-only, ExAC common | `gs://gatk-best-practices/somatic-hg38` |
| `somatic/seqc2/` (14 GB) | SEQC2 HCC1395 WES tumor/normal (FD_1), GRCh38.d1.vd1, targets, high-confidence calls | `somatic/seqc2/download.sh` (NCBI SEQC2, GATK best-practices bucket) |
| `somatic/pon/` | Mutect2 calls on the four DREAM normals, their GenomicsDB workspaces, panel-of-normals truth | `somatic/pon/truth.sh` (makes its inputs) |
| `somatic/{contam,orient,group3,longread,fmcgvcf,fmcnan,fmcdup,gathers}/` | Mutect2-family truth on the data above | each folder's `truth*.sh` / `run.sh` |
| `hcorder/`, `hclong/`, `plan_edge/` | HaplotypeCaller / planner edge cases cut from HG002 | `truth.sh` / `run.sh` (make their inputs) |
| `containers/variant-rs/truth/` | GATK 4.7.0.0 HaplotypeCaller on HG002 chr8:1-15M, run in sarek's GATK image (`containers/test_images.sh`) | GATK, from `hg002/` |
| `pipeline/nf/test/` (3.5 GB) | the pipeline tests' FASTQs (NA12878 chr20 slice in two lanes, HCC1143 lanes, SEQC2 WES), GATK truth, lane digests | `pipeline/test/make_fastq.sh`, `pipeline/test/truth.sh` (in the repository) |
| `scatter/` | GATK truth for HaplotypeCaller `--scatter-intervals` (`variant-rs/tests/scatter_intervals.sh`) | made by that test |
| `sarek/` (36 GB) | nf-core/sarek 3.12.0 source (`src/`) and its end-to-end runs with variant-rs (`e2e/`); `bin/nextflow` | the nf-core/sarek repository; `e2e/` runs |
| `tools/` | GATK 4.7.0.0 release (`gatk-release/`), bwa-mem2 2.2.1 (with the dump hooks of `docs/notes/align-rs.md`), `nextflow` (link to `sarek/bin`) | the GATK and bwa-mem2 releases |
| `upstream/` (1 GB) | the source code read while porting: GATK, Picard, htsjdk, GenomicsDB, GKL, commons-math, JEXL | their repositories at the reference versions |
| `runs/` | ad-hoc pipeline runs | |
| `archive/` (14 GB) | the prototype-era folders (before the code moved into this repository, 2026-11-17), kept as they were | |

Reference genomes are large but shared: `ref_hg38` (pipeline, HCC1143, align-rs), `wgs/ref` and `bqsr/ref` (the
1000 Genomes GRCh38 analysis set: NA12878 BQSR and MarkDuplicates), `hg002/ref` (no-alt set), `somatic/dream/data`
(b37 chr20/21).
