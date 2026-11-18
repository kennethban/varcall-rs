# variant-rs FASTQ -> VCF pipeline (Nextflow)

Germline (HaplotypeCaller) and somatic (Mutect2) short variants from FASTQ, in GATK Best Practices order, with
variant-rs for every GATK step. Each step's output is GATK 4.7's (`test/compare.py` against `test/truth.sh`).

```
per sample   align-rs mem (bwa-mem2 2.2.1's alignments, GPU)      a sorted BAM + index per lane, written directly
per sample   MarkDuplicates over all lanes -> BaseRecalibrator (whole genome)
germline     HaplotypeCaller, recal table applied on the fly (no recalibrated BAM), each interval its own GATK run
             (as sarek's scatter, but in one process) -> hard filters
somatic      ApplyBQSR -> Mutect2 (tumor + matched normal, or tumor only, with F1R2 counts)
             GetPileupSummaries -> CalculateContamination; LearnReadOrientationModel -> FilterMutectCalls
```

What makes it faster than nf-core/sarek with the same variant-rs (measured in end-to-end sarek runs whose notes are not
in this repository):

* **No external scatter.** Each step runs once per sample on the whole genome and takes the machine (`--cpus`,
  default 16); variant-rs parallelises inside (threads, contig workers, GPU PairHMM). Sarek splits the genome into ~350
  intervals and runs 2-CPU tasks per interval, each reloading the reference and resources, then gathers the pieces.
  HaplotypeCaller still treats each interval as its own GATK run (`scatter_intervals`), so its output is sarek's.
* **No recalibrated BAM for germline calling.** HaplotypeCaller applies the recal table while reading.
* **BAM throughout, multi-threaded writers.** No single-threaded samtools CRAM step.
* MarkDuplicates gets indexed lane BAMs (no in-memory indexing).
* **Alignment on the GPU** (`aligner = 'align-rs'`, default): [`align-rs`](../align-rs), a port of bwa-mem2 2.2.1 `mem` whose
  records are bwa-mem2's (checked record for record on a 30x lane), with seeding, extension and mate rescue on the
  GPU, writing the coordinate-sorted BAM and `.bai` itself (records as `samtools sort` gives them, index bytes as
  `samtools index`). One task per sample (`--lanes`): the index loads once and each lane's final BAM pass overlaps
  the next lane's alignment; each lane's BAM is the one aligning it alone gives. A 30x lane: 236 s + 28 s final
  pass vs 549 s for bwa-mem2 | samtools sort, + 17 s samtools index. Without a GPU it runs on the CPU;
  `aligner = 'bwa-mem2'` runs bwa-mem2 | samtools sort | samtools index per lane instead.

## Run

```bash
# in an environment with samtools/htslib, Java for Nextflow (and bwa-mem2): environment.yml at the repo root
nextflow run main.nf -params-file my.yaml
```

Rerunning with `-resume <session id>` (`nextflow log` lists them; a bare `-resume` takes the last run launched from
this directory, whichever it was) takes unchanged steps, the alignments above all, from the work directory. Each
variant-rs command carries the binary's build time (`VARIANT_RS_BUILD=...`), so a rebuilt variant-rs reruns the variant-rs steps
and only them. [`tests/harness/regress.sh`](../tests/harness/regress.sh) runs the tests this way: 2.5 min instead of ~12 with the 40 alignments cached.

The reference needs `.fai`, `.dict` and the bwa-mem2 index beside it (align-rs reads the same index files). variant-rs
and align-rs are built for the machine (`RUSTFLAGS="-C target-cpu=native"`; align-rs's kernels: `make -C
align-rs/cuda` from the repository root first); `--variant_rs` / `--align_rs` point at the binaries. align-rs's commands carry its build time
too (`ALIGN_RS_BUILD=...`), so a rebuilt align-rs realigns.

Samplesheet (`--input`), as nf-core/sarek's: `patient,sample,status,lane,fastq_1,fastq_2`. Status 0 is normal, 1 is
tumor; an empty `fastq_2` is a single-end lane. `--step markduplicates` starts from coordinate-sorted BAMs instead:
`patient,sample,status,bam[,bai]`, one row per BAM (a sample's BAMs are merged by MarkDuplicates; a BAM without an
index is indexed in memory). HaplotypeCaller runs on the normal samples, Mutect2 once per tumor
(with the patient's normal when there is one).

| parameter | |
|---|---|
| `tools` | `haplotypecaller`, `mutect2` (comma-separated) |
| `fasta`, `known_sites` | reference; BQSR known sites, comma-separated (dbSNP, Mills, known indels) |
| `intervals` | `-L` for HaplotypeCaller and Mutect2 (BQSR is whole genome) |
| `scatter_intervals` | `true`: HaplotypeCaller runs each interval as its own GATK run (`--scatter-intervals`), as sarek scatters it; `false`: one GATK run over all intervals; `auto` (default): `true` for up to 1,000 intervals (whole-genome calling regions), `false` above (exome targets: one run 30 s, per target 317 s) |
| `hc_workers` | HaplotypeCaller contig workers (0: variant-rs default, 8 at 16+ threads) |
| `dbsnp`, `gvcf` | HaplotypeCaller `--dbsnp`; `-ERC GVCF` (then no hard filters) |
| `germline_resource`, `pon`, `pileup_sites` | Mutect2 / GetPileupSummaries resources |
| `aligner`, `align_rs` | `align-rs` (default) or `bwa-mem2`; the align-rs binary |
| `cpus`, `align_memory` | threads for each heavy step; memory for each aligner task (the ~17 GB hg38 index is mapped, and shared by aligner tasks running together) |
| `save_lane_bams`, `save_recal_bam` | publish the lane BAMs / the somatic recalibrated BAMs |
| `publish_dir_mode` | `link` (default: hard links into the work directory, no second copy) |

## Several GPUs

List the GPUs in `CUDA_VISIBLE_DEVICES` when you launch, and give each heavy step a share of the cores:

```bash
# 8 GPUs, 128 cores
CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7 nextflow run main.nf -params-file my.yaml --cpus 16
```

* Nextflow gives each GPU step (align-rs, HaplotypeCaller, Mutect2) one of the listed GPUs, and runs at most one GPU
  step per GPU at a time.
* Each sample's lanes are dealt out over one align-rs task per GPU, so one sample's alignment uses all of them. A
  lane's BAM is the same whichever task aligns it.
* HaplotypeCaller and Mutect2 run once per sample (or tumor), each on one GPU, so several samples call side by side.
  One sample's calling uses one GPU.
* Steps run side by side only while their `cpus` and memory fit the machine: 24 GB per heavy step and `align_memory`
  per aligner task. With 8 GPUs, that is about 200 GB, or lower `align_memory`.

Without `CUDA_VISIBLE_DEVICES`, every step uses the first GPU, as before, and the pipeline warns if it finds more.
Outputs do not depend on the number of GPUs. Changing it does make `-resume` realign, because the lanes are grouped
into different tasks.

Tested (2026-11-18) on one GPU with `CUDA_VISIBLE_DEVICES=0,1,2,3`, `--cpus 8 --align_memory 6.GB`. Four align-rs
tasks ran at once. The tasks given GPUs 1-3, which do not exist here, fell back to the CPU. The somatic test's 17
outputs were identical, and all 38 lane BAMs had the records bwa-mem2 writes. Not yet run on a machine with several
GPUs.

Outputs: `preprocessing/<sample>/` (MarkDuplicates BAM + metrics, recal table), `germline/<sample>/` (VCF, filtered
VCF), `somatic/<tumor>_vs_<normal>/` (Mutect2 VCF + stats, F1R2, priors, contamination, segments, filtered VCF),
`somatic/pileups/`.

## Benchmarks

Final runs (2026-11-18), on one workstation (16-core Ryzen 9 9950X3D, 30 GB RAM, RTX 5070 Ti), nothing cached:

| Run | Input | Wall | Steps | Reference |
|---|---|---:|---|---|
| germline whole genome | NA12878 30x, 12 FASTQ lanes, Broad hg38 | **47.8 min** | ALIGN 35.4 min, MarkDuplicates 3.9, BaseRecalibrator 1.2, HaplotypeCaller 7.0, hard filters 0.1 | stock tools on the same machine: ~5.2 h (bwa-mem2 + samtools ~113 min, sarek's stock GATK steps ~196 min) |
| somatic exome | SEQC2 HCC1395 tumor / normal BAMs | **2.7 min** | as in the 2026-11-16 run below | nf-core/sarek: 46.5 min with stock GATK, 15.3 min with variant-rs |

Both runs' outputs are identical to the earlier runs' below (germline: VCF 4,951,704 records, 4,944,053 after
filtering, recal table, duplicate metrics; somatic: all 13 outputs). The sections that follow are the measurement
log: how each number was taken, and what changed between runs.

## Tests (2026-11-16)

`test/make_fastq.sh` makes the inputs; `test/truth.sh germline|somatic` runs the same chain with GATK 4.7 from the
pipeline's own lane BAMs (alignment is not a GATK step); `test/compare.py germline|somatic` compares.

| test | input | result |
|---|---|---|
| germline | NA12878 chr20:10-15M, 1.19M pairs cut into 2 lanes; HaplotypeCaller on two intervals (`test/germline.intervals`, chr20:10-12.5M and 12.5-15M; truth: one GATK run per interval, records concatenated) (also from the 2 lane BAMs, `--step markduplicates`: same outputs) | 5 files identical: MarkDuplicates BAM (1,192,285 records) and metrics, recal table, VCF (8,039 records), filtered VCF |
| somatic_bam | SEQC2 HCC1395 exome pair, chr22 (2.2M / 1.6M reads), from BAMs whose read-group SM (`WES_FD_*`) is not the samplesheet's sample name; the real target BED's chr22 part (6,644 regions, duplicated and overlapping ones included); Mutect2 on the pair, HaplotypeCaller on the normal | 17 files identical: both samples' MarkDuplicates BAMs, metrics and tables; Mutect2 VCF (94 records) and stats, F1R2, priors, pileups, contamination, segments, filtered VCF; HaplotypeCaller VCF (2,299 records) and filtered VCF |
| somatic | HCC1143 tumor / normal (GATK somatic workshop, hg38 chr17), one FASTQ per read group: normal 16 paired lanes, tumor 18 paired + 4 single-end lanes | 17 files identical: both samples' MarkDuplicates and recalibrated BAMs, metrics and tables; Mutect2 VCF (359 records) and stats, F1R2 counts, priors, pileups, contamination, segments, filtered VCF |

The F1R2 alt histograms are compared by column name: GATK writes their columns in an order that follows Java identity
hash codes (a HashMap keyed by enum pairs), which differs between GATK runs; the counts are identical and
LearnReadOrientationModel reads them by name.

Somatic test: 10.6 min wall with bwa-mem2, most of it the 38 alignments (each loads the whole-genome index, 15-27 s
per lane); with align-rs 2.9 min (2.5 s per lane: the mapped index stays in the page cache), lane BAMs' records
identical to bwa-mem2's (all 40 test lanes, `test/log/lanes.bwamem2.md5`); germline test 36 s. With one ALIGN
task per sample (2026-11-17): somatic 103 s, germline 22 s, outputs unchanged. The GATK steps
take seconds (MarkDuplicates 6 s, BaseRecalibrator 14-25 s, Mutect2 9 s) where GATK took 20-40 s each and 109 s for
Mutect2.

Found on the way: the pipeline's test data first had empty tumor FASTQs (the workshop's tumor BAM is single-end),
which exposed a variant-rs crash on an empty lane BAM, now fixed ([`tests/harness/markdup/edge/run.sh`](../tests/harness/markdup/edge/run.sh)).

## Whole genome (2026-11-16): NA12878 30x from its 41 GB BAM, as the sarek `wgs_germline` run

`profile/run_profiled.sh wgs/na12878.yaml NAME` (the profiling scripts in `profile/` and the run configurations and logs
in `wgs/` and `wes/` were development files, not in the repository): same reference, known sites and HaplotypeCaller intervals as sarek;
32 CPUs per step. Machine-wide samples every 2 s (`profile/sample.py`), per task (`profile/report.py`; "busy" is the
mean number of busy cores, "< half" the share of the task's time with fewer than 16 busy):

| step | run 1 | run 2 (parallel BQSR in HaplotypeCaller) | run 2 busy cores / < half | iowait | disk read / write | GPU |
|---|---|---|---|---|---|---|
| MarkDuplicates | 4.3 min | 4.5 min | 18.5 / 44% | 0.6 | 533 / 419 MB/s | - |
| BaseRecalibrator (whole genome) | 1.3 min | 1.2 min | 30.7 / 8% | 0.1 | 622 MB/s | - |
| HaplotypeCaller (BQSR on the fly, GPU) | 17.9 min | **10.9 min** | 13.0 / 73% | 0.0 | 65 MB/s | 22% |
| hard filters | 6 s | 6 s | | | | |
| **total (wall)** | **23.7 min** | **16.8 min** | | | | |

Run 2's VCF (4,951,737 records) and recal table are byte-identical to run 1's. Sarek with variant-rs, same steps: ~95 min
of wall clock (MarkDuplicates 55.9, BaseRecalibrator 17.5, ApplyBQSR 10.6, HaplotypeCaller 10.3), plus ~60 min of
sarek's own samtools / mosdepth steps. (The peak busy-core samples can exceed 32: WSL's /proc/stat counters arrive
unevenly; the means are what to read.)

No step waits on the disk (iowait < 1 core). Run 1's profile showed HaplotypeCaller's front end as the bottleneck:
7,993 s summed over contigs against ~700 s for the whole back end, because on-the-fly BQSR recalibrated every read on
the per-contig reading thread (chr20: 19.0 s without BQSR, 31.5-32.4 s with it). Reads are now recalibrated in batches
on a small pool of their own, several batches in flight, while parsing and the rest stay in file order: chr20 with
BQSR 20.1-20.8 s, the same as without; output identical ([`tests/harness/regress.sh`](../tests/harness/regress.sh) all identical).

### HaplotypeCaller per interval (`scatter_intervals`, default since 2026-11-16)

One GATK run draws from one random generator (`java.util.Random`, seed 47382911) in genome order: downsampling at
start positions with more than 50 reads (16,256 of them here) and the QD jitter. variant-rs reproduces that order, so in
one run each contig waits at its first such position until every earlier contig is finished: this, not contig
length, held run 2's HaplotypeCaller at 12.5 busy threads. Cutting contigs inside one run would only add waits.
Treating each interval as its own GATK run (each with its own generator, as sarek's 356 HaplotypeCaller tasks are)
removes the chain:

| HaplotypeCaller, same BAM and recal table | wall | busy threads (mean) | peak RSS |
|---|---|---|---|
| one GATK run over all intervals (run 2) | 672 s | 12.5 | 6.3 GB |
| each of the 356 intervals its own run, 8 workers | 407-419 s | 22.6 | 8.8 GB |
| same, 12 / 16 workers | 405 / 407 s | 23.3 | 12.1 / 15.5 GB |

* **Equals sarek's stock GATK:** with sarek's gathered recal table and `--dbsnp`, the VCF's 4,954,925 records are
  identical to sarek's stock `NA12878.haplotypecaller.vcf.gz` (356 GATK runs merged by MergeVcfs).
* Against one GATK run over all intervals: 4,090,492 records identical, 856,064 differ only in QD (the jitter), about
  5,600 high-depth sites differ otherwise or are in only one file (downsampling draws). GATK itself differs the same
  way between a scattered and an unscattered run.
* More workers do not help: the 9950X3D has 16 cores (32 threads), and ~23 busy threads is close to what SMT adds;
  every stage's summed time grows with the worker count. 8 stays the default.
* Run 3, the whole pipeline with it (`wgs/na12878_run3.*`): **12.75 min wall** (765 s; run 2: 16.8 min). MarkDuplicates
  4.5 min (19.1 busy threads), BaseRecalibrator 1.2 min (32.0), HaplotypeCaller **6.8 min** (24.4 busy, under half
  the threads 13% of its time, was 73%; GPU 32%), hard filters 6 s. VCF records identical to the 356-shard run's,
  recal table identical to run 2's, MarkDuplicates metrics identical apart from the start time.

MarkDuplicates (2026-11-16, after run 3): 247 -> 196 s on this BAM (the dense small contigs split, pass 2 pipelined;
see [docs/notes/variant-rs.md](../docs/notes/variant-rs.md)), under 16 busy threads 20% of its time instead of 45%. Whole run estimate: ~12 min.

### From FASTQ with align-rs (2026-11-17)

`wgs/na12878_fastq.yaml`: the 12 lanes as FASTQ (extracted from the BAM by [`tests/harness/wgs/fastq/extract_all.sh`](../tests/harness/wgs/fastq/extract_all.sh)), Broad
hg38 (`ref_hg38`, which has the bwa-mem2 index), same known sites and intervals. **81 min wall**: 12 ALIGN tasks
69.3 min (5.6-6.0 min each, one at a time: GPU and memory), then (`wgs/na12878_fastq_md.*`, run from the lane BAMs
after a full disk stopped the first run at MarkDuplicates) MarkDuplicates over the 12 lane BAMs 3.5 min,
BaseRecalibrator 1.2, HaplotypeCaller 6.7, hard filters 0.1. VCF 4,951,704 records (4,944,053 after filtering).
Stock tools on this machine, same steps: ~5.2 h (bwa-mem2 + samtools ~113 min, sarek's stock GATK spans ~196 min).

### From FASTQ again, after the aligner speed work (2026-11-17)

`wgs/na12878_fastq2.yaml` (same inputs, one ALIGN task for all 12 lanes): **48.8 min wall** (2,925 s), one run,
nothing cached. ALIGN 36.1 min (index load, 12 lanes at ~2.0 s per batch, each lane's final BAM pass overlapping
the next lane), MarkDuplicates 4.1, BaseRecalibrator 1.2, HaplotypeCaller 7.0, hard filters 0.1. Outputs
identical to the first FASTQ run's: VCF data lines (4,951,704; 4,944,053 after filtering) and the duplicate
metrics. Against the first run (81 min) 1.66x; against stock tools on this machine (~5.2 h) ~6.4x.

## Somatic exome (2026-11-16): SEQC2 HCC1395 tumor/normal, as the sarek `wes_somatic` run

`profile/run_profiled.sh wes/seqc2.yaml wes_seqc2`: the sarek run's inputs and resources (GRCh38.d1.vd1, ~70M reads
per sample from BAMs, 295k-region target BED, af-only gnomAD, 1000G PoN), HaplotypeCaller on the normal, Mutect2 on
the pair; 32 CPUs per step. **3.3 min wall** (200 s) against sarek 46.5 min (stock GATK) and 15.3 min (with variant-rs).
Steps: MarkDuplicates 0.2-0.3 min, BaseRecalibrator 0.1-0.2, ApplyBQSR 0.3, HaplotypeCaller 0.5 (one run over the
targets: `scatter_intervals` auto), GetPileupSummaries 0.4 each (~3 busy threads; since fixed, 23.4 -> 3.4 s each: see below), Mutect2 0.4, the rest seconds.
HaplotypeCaller VCF 89,817 records, Mutect2 5,485 (1,321 PASS); sarek: 89,885 / 5,638, from its gathered recal
tables and per-interval runs (not expected to be equal; the tools themselves match GATK in the regression suite).

Found on the way: the target BED repeats regions (now merged for `--scatter-intervals`; touching intervals stay
separate runs); with BAM input, Mutect2's `-normal` now takes the normal's SM from its BAM (here `WES_FD_N_1`, not
the samplesheet's `FD_N_1`); a rare race at the end of a contig in HaplotypeCaller's two-stage front end, hit once
in 294,884 per-target runs, fixed (`wait_until` re-checks progress after the other stage's stop flag).

GetPileupSummaries over the 295k targets spent 23 of its 26 s choosing sites: the gnomAD index sends every target's
query to the start of its 16 kb window, so 96M resource lines were parsed to keep 165k sites, on one thread. Lines
starting outside the targets are now only peeked at (POS and REF length, which still decide GATK's first-record rule)
and the contigs run in parallel: 3.4 s, tables identical. Tests added for what this run found:
`test/somatic_bam` (above) and [`variant-rs/tests/scatter_intervals.sh`](../variant-rs/tests/scatter_intervals.sh) (duplicate / overlapping / contained /
touching intervals vs GATK per interval; 6,641 exome pieces twice), both in [`tests/harness/regress.sh`](../tests/harness/regress.sh).
