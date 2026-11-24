# varcall-rs FASTQ → VCF pipeline (Nextflow)

This Nextflow pipeline calls germline short variants (HaplotypeCaller, per sample or jointly over a cohort) and
somatic short variants (Mutect2) from FASTQ files. It follows the order of the GATK Best Practices and uses align-rs
for alignment and variant-rs for every GATK step. The output of each step is identical to the output of GATK 4.7 for
the same input ([Tests](#tests)).

The steps of the pipeline are as follows:

```
per sample   align-rs mem (bwa-mem2 2.2.1's alignments, GPU)      a sorted BAM + index per lane, written directly
per sample   MarkDuplicates over all lanes -> BaseRecalibrator (whole genome)
germline     HaplotypeCaller, recal table applied on the fly (no recalibrated BAM), each interval its own GATK run
             (as sarek's scatter, but in one process) -> hard filters
joint        (--joint_germline) HaplotypeCaller GVCFs -> GenomicsDBImport -> GenotypeGVCFs -> VQSR (SNPs, then INDELs),
             as sarek's joint germline calling: the intervals imported side by side, each interval genotyped as
             its own GATK run (sarek's per-interval tasks, in one process)
somatic      ApplyBQSR -> Mutect2 (tumor + matched normal, or tumor only, with F1R2 counts)
             GetPileupSummaries -> CalculateContamination; LearnReadOrientationModel -> FilterMutectCalls
```

The pipeline is faster than nf-core/sarek configured with variant-rs (2.8 min compared with 15.3 min on the somatic
exome described under [Benchmarks](#benchmarks)). Five differences in design account for this:

* **No external scatter.** Each step runs once per sample on the whole genome and uses all of the CPUs given to it
  (`--cpus`), and variant-rs parallelizes the work internally (threads, contig workers, GPU PairHMM). In contrast,
  nf-core/sarek divides the genome into ~350 intervals and runs one 2-CPU task per interval, each of which reloads
  the reference and the resources, and then gathers the pieces. HaplotypeCaller still treats each interval as a
  separate GATK run (`scatter_intervals`), so its output is identical to that of nf-core/sarek.
* **No recalibrated BAM for germline calling.** HaplotypeCaller applies the recalibration table to the reads as it
  reads them.
* **BAM files throughout, written by multi-threaded writers.** The pipeline has no single-threaded samtools CRAM
  step.
* **Indexed lane BAMs as the input of MarkDuplicates,** so that no input is indexed in memory.
* **Alignment on the GPU** (`aligner = 'align-rs'`, the default). [align-rs](../align-rs) produces the records of
  bwa-mem2 2.2.1 and writes the coordinate-sorted BAM and its `.bai` directly. One task aligns all lanes of a sample
  (`--lanes`), so the index is loaded once and the final BAM pass of each lane overlaps the alignment of the next
  lane. A 30x lane takes 124 s, including sorting and indexing, compared with 572 s for `bwa-mem2 | samtools sort`
  and a further 16 s for `samtools index`. Without a GPU, align-rs runs on the CPU. With `aligner = 'bwa-mem2'`, the
  pipeline runs `bwa-mem2 | samtools sort | samtools index` on each lane instead.

## Run

```bash
# in an environment with samtools/htslib, Java for Nextflow (and bwa-mem2): environment.yml at the repo root
./run.sh -params-file my.yaml
```

`run.sh` selects the GPUs and then runs `nextflow run main.nf` with the arguments given to it. The pipeline can also
be started directly with `nextflow run main.nf`, in which case `CUDA_VISIBLE_DEVICES` lists the GPUs
([Several GPUs](#several-gpus)).

### Requirements

* **Binaries.** variant-rs and align-rs should be built for the machine on which they run
  (`RUSTFLAGS="-C target-cpu=native"`, and `make -C align-rs/cuda` from the repository root for the GPU kernels of
  align-rs). The parameters `--variant_rs` and `--align_rs` give the paths of the binaries.
* **Reference.** The reference FASTA file requires its `.fai`, `.dict` and bwa-mem2 index in the same directory.
  align-rs reads the index files of bwa-mem2.
* **Samplesheet.** The samplesheet (`--input`) has the format of nf-core/sarek, with the columns
  `patient,sample,status,lane,fastq_1,fastq_2`. A status of 0 denotes a normal sample and 1 a tumor sample, and an
  empty `fastq_2` denotes a single-end lane.

The pipeline can also start from coordinate-sorted BAM files. This requires `--step markduplicates` and a samplesheet
with the columns `patient,sample,status,bam[,bai]` and one row per BAM file. MarkDuplicates merges the BAM files of a
sample, and a BAM file without an index is indexed in memory.

HaplotypeCaller runs on the normal samples. Mutect2 runs once per tumor sample, with the normal sample of the same
patient when one is present.

### GPUs and threads

`run.sh` uses the GPUs on which no process is running (`VARCALL_GPUS=N` limits them to N). It leaves the selection
unchanged when `CUDA_VISIBLE_DEVICES` is already set, either by the user or by a scheduler such as SLURM. The number
of threads follows from the machine and from the number of GPUs in use:

| Step | Tasks | Threads per task |
|---|---|---|
| align-rs | one task per GPU, run concurrently | the threads of the machine divided by the number of GPUs, at most 64 (`--align_cpus`) |
| all later steps | one sample at a time | the threads of the machine, at most 128 (`--cpus`) |

The two limits are the thread counts beyond which no further gain was measured on a 30x genome on the L40S server.
MarkDuplicates took 75 s at 64 threads and 56 s at 128. BaseRecalibrator took 30, 19 and 16 s at 64, 128 and 192
threads, and HaplotypeCaller took 121 and 118 s at 64 and 128 threads. Alignment was measured with up to 64 threads
per task.

### Resuming

With `-resume <session id>`, Nextflow takes the results of unchanged steps from the work directory, which saves the
alignments in particular (`nextflow log` lists the session ids). A `-resume` without a session id resumes the last
run that was launched from the current directory, whichever run that was.

Each command includes the build time of its binary (`VARIANT_RS_BUILD=...`, `ALIGN_RS_BUILD=...`). After variant-rs
has been rebuilt, a resumed run therefore repeats the variant-rs steps and no others, and after align-rs has been
rebuilt it repeats the alignments. [`tests/harness/regress.sh`](../tests/harness/regress.sh) runs the pipeline tests
in this way, which takes 2.5 min with the 40 alignments cached and ~12 min without the cache.

### Parameters

| Parameter | Description |
|---|---|
| `tools` | The callers to run: `haplotypecaller`, `mutect2` (comma-separated). |
| `fasta`, `known_sites` | The reference, and the known sites for BQSR as a comma-separated list (dbSNP, Mills, known indels). |
| `intervals` | The intervals (`-L`) for HaplotypeCaller and Mutect2. BQSR always uses the whole genome. |
| `scatter_intervals` | With `true`, HaplotypeCaller and the GenotypeGVCFs step of joint calling treat each interval as a separate GATK run (`--scatter-intervals`, run concurrently), as nf-core/sarek scatters them. With `false`, they perform one GATK run over all intervals. With `auto` (the default), the value is `true` for up to 1,000 intervals (whole-genome calling regions) and `false` for more (exome targets, for which one run takes 30 s and one run per target takes 317 s). |
| `hc_workers` | The number of contig workers of HaplotypeCaller. A value of 0 selects the default of variant-rs, which is 8 at 16+ threads. |
| `dbsnp`, `gvcf` | `dbsnp` gives the `--dbsnp` file for HaplotypeCaller and for GenotypeGVCFs in joint calling. `gvcf` selects `-ERC GVCF`, in which case no hard filters are applied. |
| `joint_germline` | Joint calling over all germline samples, equivalent to `--joint_germline` in nf-core/sarek. The GVCFs are imported into one GenomicsDB workspace (one array per interval, imported concurrently; every contig of the reference when `intervals` is not given) and genotyped with GenotypeGVCFs (per interval with `scatter_intervals`, as the tasks of nf-core/sarek). VQSR follows when its resources are given. No per-sample hard filters are applied. |
| `wes` | For joint calling on exome targets, GenomicsDBImport runs with `--merge-input-intervals` (equivalent to `--wes` in nf-core/sarek). |
| `vqsr_snp_resources`, `vqsr_indel_resources` | The resources of VariantRecalibrator, given as `--resource:name,known=…,training=…,truth=…,prior=… file` arguments (in nf-core/sarek, `known_snps_vqsr` / `known_indels_vqsr` plus `dbsnp_vqsr`). When both are given, VQSR runs with the annotations used by nf-core/sarek, and ApplyVQSR runs at a truth sensitivity of 99.9 (SNPs, then INDELs). |
| `germline_resource`, `pon`, `pileup_sites` | The resources for Mutect2 and GetPileupSummaries. |
| `aligner`, `align_rs` | The aligner, `align-rs` (the default) or `bwa-mem2`, and the path of the align-rs binary. |
| `cpus`, `align_cpus`, `align_memory` | `cpus` is the number of threads for each heavy step other than alignment (by default the threads of the machine, at most 128). `align_cpus` is the number of threads for each aligner task (by default the threads of the machine divided by the number of GPUs in use, at most 64). `align_memory` is the memory for each aligner task; the ~17 GB hg38 index is memory-mapped and shared by aligner tasks that run concurrently. |
| `hc_gpus` | The number of GPUs for each HaplotypeCaller task, taken from the GPUs in use (default 2). |
| `save_lane_bams`, `save_recal_bam` | Publish the lane BAMs, and the recalibrated BAMs of the somatic workflow. |
| `publish_dir_mode` | The default, `link`, publishes hard links into the work directory, so no second copy is written. |

## Several GPUs

`run.sh` uses the idle GPUs (`VARCALL_GPUS=4 ./run.sh ...` limits them to 4). To select particular GPUs, they are
listed in `CUDA_VISIBLE_DEVICES` when the pipeline is launched. Each alignment task receives the threads of the
machine divided by the number of GPUs, so that all GPUs align concurrently:

```bash
CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7 ./run.sh -params-file my.yaml   # 8 GPUs, 128 threads: 16 per alignment task
```

* Nextflow assigns one of the listed GPUs to each GPU step (align-rs, HaplotypeCaller, Mutect2) and runs at most one
  GPU step per GPU at a time.
* The lanes of each sample are distributed over one align-rs task per GPU, so that the alignment of one sample uses
  all GPUs. The BAM of a lane is the same regardless of the task that aligns it.
* HaplotypeCaller runs once per sample on up to 2 GPUs (`hc_gpus`) and distributes its PairHMM over them. On a 30x
  genome with L40S GPUs, it took 115 s on 1 GPU, 98 s on 2 GPUs and 97 s on 4 GPUs, because the longest calling
  interval, 141 Mb of chr2, determines the time. Mutect2 runs once per tumor sample on one GPU. Several samples are
  called concurrently as far as the GPUs and threads allow.
* Steps run concurrently only when their `cpus` and memory fit the machine. Each heavy step is allotted 24 GB, and
  each aligner task is allotted `align_memory`. With 8 GPUs this requires about 200 GB, or a lower `align_memory`.

When `nextflow run main.nf` is started without `CUDA_VISIBLE_DEVICES`, it uses the first GPU for every step and
prints a warning if it finds more than one. The outputs do not depend on the number of GPUs. A change in the number
of GPUs does cause `-resume` to repeat the alignments, because the lanes are then grouped into different tasks.

The use of several GPUs was tested on 2026-11-19 on a server with 8 L40S GPUs, of which 6 were free. NA12878 30x
from FASTQ on 6 GPUs gave a VCF identical to that of the run on one GPU, and the SEQC2 exome pair gave outputs
identical to those of the workstation ([Benchmarks](#benchmarks)). On a machine with one GPU, a list of more GPUs
(`CUDA_VISIBLE_DEVICES=0,1,2,3`) does no harm, because the tasks that are assigned a GPU that does not exist run on
the CPU and give the same outputs.

## Outputs

| Folder in `outdir` | Contents |
|---|---|
| `preprocessing/<sample>/` | MarkDuplicates BAM and metrics, recalibration table |
| `germline/<sample>/` | VCF and filtered VCF; the GVCF with `joint_germline` |
| `germline/joint/` | `joint_germline.vcf.gz`, the VQSR recal files (indexed) and tranches, `joint_germline_recalibrated.vcf.gz` |
| `somatic/<tumor>_vs_<normal>/` | Mutect2 VCF and stats, F1R2 counts, priors, contamination, segments, filtered VCF |
| `somatic/pileups/` | GetPileupSummaries tables |

## Tests

The tests of the pipeline compare every output file with the corresponding output of GATK. They were first run on
2026-11-16 and consist of three scripts:

* `test/make_fastq.sh` prepares the inputs;
* `test/truth.sh germline|somatic|somatic_bam|joint` runs the same sequence of steps with GATK 4.7, starting from the
  lane BAMs written by the pipeline, because alignment is not a GATK step;
* `test/compare.py <test>` compares the two sets of outputs.

| Test | Input | Result |
|---|---|---|
| germline | NA12878 chr20:10-15M, 1.19M pairs divided into 2 lanes; HaplotypeCaller on two intervals (`test/germline.intervals`, chr20:10-12.5M and 12.5-15M; truth: one GATK run per interval, records concatenated). The same outputs are obtained from the 2 lane BAMs with `--step markduplicates`. | 5 files identical: MarkDuplicates BAM (1,192,285 records) and metrics, recal table, VCF (8,039 records), filtered VCF |
| somatic_bam | SEQC2 HCC1395 exome pair, chr22 (2.2M / 1.6M reads), from BAMs whose read-group SM (`WES_FD_*`) differs from the sample name in the samplesheet; the chr22 part of the real target BED (6,644 regions, including duplicated and overlapping regions); Mutect2 on the pair, HaplotypeCaller on the normal | 17 files identical: MarkDuplicates BAMs, metrics and tables of both samples; Mutect2 VCF (94 records) and stats, F1R2, priors, pileups, contamination, segments, filtered VCF; HaplotypeCaller VCF (2,299 records) and filtered VCF |
| somatic | HCC1143 tumor / normal (GATK somatic workshop, hg38 chr17), one FASTQ per read group: normal 16 paired lanes, tumor 18 paired + 4 single-end lanes | 17 files identical: MarkDuplicates and recalibrated BAMs, metrics and tables of both samples; Mutect2 VCF (359 records) and stats, F1R2 counts, priors, pileups, contamination, segments, filtered VCF |
| joint (2026-11-21) | NA12878 (the 2 lane BAMs of the germline test) and HCC1143_normal (the 16 lane BAMs of the somatic test), `--step markduplicates`, `joint_germline` over chr17:7-9.5M and the two chr20 intervals of the germline test, `--dbsnp`, VQSR with the HapMap / Omni / 1000G / dbSNP and Mills / dbSNP resources | 14 files identical: MarkDuplicates BAMs, metrics, recal tables and GVCFs of both samples (805,650 and 58,712 records); the joint VCF (8,400 records), both recal files and tranches, the VQSR-filtered VCF. The joint steps take 3-4 s each |

The F1R2 alt histograms are compared by column name. GATK writes their columns in an order that follows Java
identity hash codes (a HashMap keyed by enum pairs), and this order differs between GATK runs. The counts are
identical, and LearnReadOrientationModel reads the columns by name.

The somatic test takes 10.6 min of wall-clock time with bwa-mem2, most of which is spent in the 38 alignments,
because each loads the whole-genome index (15-27 s per lane). With align-rs it takes 2.9 min (2.5 s per lane, because
the memory-mapped index remains in the page cache), and the records of the lane BAMs are identical to those of
bwa-mem2 for all 40 test lanes (`test/log/lanes.bwamem2.md5`). The germline test takes 36 s. With one ALIGN task per
sample (2026-11-17), the somatic test takes 103 s and the germline test 22 s, with unchanged outputs. The GATK steps
take a few seconds in variant-rs (MarkDuplicates 6 s, BaseRecalibrator 14-25 s, Mutect2 9 s), whereas GATK took
20-40 s for each step and 109 s for Mutect2.

The preparation of these tests revealed one defect. The first version of the test data contained empty tumor FASTQ
files, because the tumor BAM of the workshop is single-end, and variant-rs crashed on the resulting empty lane BAM.
The defect has been corrected and is covered by
[`tests/harness/markdup/edge/run.sh`](../tests/harness/markdup/edge/run.sh).

## Benchmarks

Every run started with no cached results. The outputs were identical in all runs of the same data: for the germline
genome, a VCF of 4,951,704 records (4,944,053 after filtering), and for the somatic exome, all 17 outputs (Mutect2:
5,485 records, of which 1,321 PASS). Step times are in minutes unless a unit is given.

### Workstation

The workstation has a 16-core Ryzen 9 9950X3D (32 threads) and one RTX 5070 Ti (16 GB), and runs Windows 11 with
WSL2 (48 GB of RAM given to Linux; ext4 on the virtual disk of WSL2).

| Run | Wall | Steps | Reference tools | Date, commit |
|---|---:|---|---|---|
| **Germline whole genome**: NA12878 30x, 12 FASTQ lanes, Broad hg38; three runs | **34.0 min** (2030, 2055, 2036 s) | ALIGN 23.5 min, MarkDuplicates 3.4, BaseRecalibrator 1.3, HaplotypeCaller 5.6; every output identical to that of commit 7d6ffa6 | stock tools on the same machine: ~5.2 h (bwa-mem2 + samtools ~117 min, stock GATK steps of nf-core/sarek ~196 min) | 2026-11-23, 736009f |
| The same, with memory limited to 30 GB (a cgroup, 16 GB of swap; Nextflow told 30 GB) | **41.6 min** | ALIGN 31.0 min (peak 26.9 GB, swap at most 0.9 GB, because pages of the hg38 index are dropped and read again), BaseRecalibrator 1.4, the other steps as above; VCF identical | | |
| **Somatic exome**: SEQC2 HCC1395 tumor / normal BAMs | **2.8 min** (167 s) | MarkDuplicates 0.2 min, BaseRecalibrator 0.1-0.15, ApplyBQSR 0.3-0.4, HaplotypeCaller 0.5, Mutect2 0.5, the other steps seconds (peak 5.1 GB, Mutect2); all 17 outputs identical to those of 2026-11-19 | nf-core/sarek: 46.5 min with stock GATK, 15.3 min with variant-rs | 2026-11-22 |

### Server

The server has 2 x AMD EPYC 9654 (192 cores, 384 threads, 2 NUMA nodes), 2.2 TiB of RAM, 8 x NVIDIA L40S (48 GB; two
of them in use by another job) and an NVMe SSD, and runs Ubuntu 24.04.

**Germline whole genome with matched resources.** The input is the same as on the workstation. The runs used GPUs
4-7 and the 128 threads of their socket (`taskset`: 64 cores + SMT), with 128 threads for each step after alignment
and HaplotypeCaller on 2 GPUs.

| GPUs | Wall | Steps |
|---|---:|---|
| 4 (32 threads per alignment task), three runs | **9.3 min** (567, 552, 558 s) | ALIGN 5.4 min, MarkDuplicates 1.0, BaseRecalibrator 0.5, HaplotypeCaller 2.2; every output identical to that of commit 7b3df92 |
| 2 (GPUs 4-5; 64 threads per alignment task) | 11.1 min (667 s) | ALIGN 6.8 min (2 tasks of 6 lanes), MarkDuplicates 1.1, BaseRecalibrator 0.5, HaplotypeCaller 2.5 |
| 1 (GPU 4; one alignment task of 128 threads) | 17.6 min (1053 s) | ALIGN 12.7 min, MarkDuplicates 1.2, BaseRecalibrator 0.6, HaplotypeCaller 2.8 |

The three runs on 4 GPUs were made on 2026-11-23 at commit 736009f.

**NVIDIA Parabricks 4.8.0-1 `germline`** on the same genome and the same resources. Parabricks produces outputs that
differ from those of GATK, and the last column compares its VCF with the VCF that is identical to that of GATK.

| Resources | Wall | Steps | Output compared with GATK |
|---|---:|---|---|
| GPUs 4-7 and the same 128 threads as varcall-rs (`--cpuset-cpus`) | 16.8 min | BWA + sort 5.9 min (64 CPU worker threads, chosen by its auto-configuration), MarkDuplicates + BQSR 3.0, HaplotypeCaller 6.8 | 99.85% of sites with the same alleles, and 99.98% of those with the same genotype; 3.18% identical records (SNVs 0.55%, indels 14.65%); QD, PL and QUAL differ in 81-85% of shared records |

This run was made on 2026-11-22.

**Germline whole genome with the defaults of `run.sh`,** which used the 6 free GPUs and all 384 threads.

| Configuration | Wall | Steps |
|---|---:|---|
| 6 GPUs, 64 threads per alignment task, 128 per other step, HaplotypeCaller on 2 GPUs | **6.8 min** (409 s) | ALIGN 3.7-3.8 min (6 tasks of 2 lanes, run concurrently; peak 37 GB each), MarkDuplicates 0.8, BaseRecalibrator 0.3, HaplotypeCaller 1.7, hard filters 0.1 |

This run was made on 2026-11-22.

**Other runs.**

| Run | Wall | Steps | Date |
|---|---:|---|---|
| Joint germline, CEPH trio (NA12878, NA12891, NA12892; 30x each, 36 FASTQ lanes) with VQSR | **~34 min** | ALIGN 13.7 min (18 tasks of 2 lanes); GVCFs complete at 21 min (per sample: MarkDuplicates 1.8 min, BaseRecalibrator 0.5-0.9, HaplotypeCaller 4.6-5.2); GenomicsDBImport 5.1, GenotypeGVCFs 5.0, VariantRecalibrator 1.9 (SNP) / 0.5 (INDEL) run concurrently, ApplyVQSR 1.0; all 18 outputs identical to GATK (`tests/harness/m3`) | 2026-11-21 |
| Somatic exome, SEQC2 tumor / normal BAMs, with the defaults on GPUs 2-7 (128 threads per step) | **63 s** | MarkDuplicates 11 s, BaseRecalibrator 3.5, ApplyBQSR 11-12, Mutect2 9.8, HaplotypeCaller on the normal 41.9; all 17 outputs identical to those of the workstation | 2026-11-19 |
| The same, on 4 GPUs (4-7) and the 128 threads of their socket (`taskset`, `--cpus 128`) | **109 s** | MarkDuplicates 4.2-4.5 s, BaseRecalibrator 3.4-3.7, ApplyBQSR 10-11, Mutect2 10.3, HaplotypeCaller on the normal 37.2. Each step is as fast or faster, but the total is longer because with 128 CPUs Nextflow runs one 128-thread step at a time, whereas with 384 it ran the steps of the tumor and of the normal concurrently; all 17 outputs identical | 2026-11-22 |

The [measurement log](../docs/notes/pipeline-benchmarks.md) records how each number was obtained, together with the
results of earlier configurations.
