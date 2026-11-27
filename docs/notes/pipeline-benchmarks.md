# Pipeline benchmarks: measurement log

This log was written during development and records how each number in
[pipeline/README.md](../../pipeline/README.md#benchmarks) was taken. It is partly historical, and later sections
supersede the numbers of earlier ones.

## Server comparisons (2026-11-23 to 2026-11-27, at 736009f)

The results are in [pipeline/README.md](../../pipeline/README.md#server); this section adds what is not there. The
runs are in `runs/platforms/` of the server's working area (one directory per run, each with a `done` file holding
the exit status and wall time, and `identity.txt` where outputs were compared). Drivers ran the steps one after
another and logged the commit, a clean tree and the binary checksums after each one.

**Single steps on S** (varcall-rs, three runs; GATK one run, from the earlier single-step measurement). MarkDuplicates,
BaseRecalibrator and ApplyBQSR on HG002 PCR-plus, HaplotypeCaller on NA12878 in 23 shards. Every output identical to
GATK's in every run.

| Step | GATK (s) | varcall-rs (s) | Mean (s) |
|---|---:|---:|---:|
| MarkDuplicates | 3,477.3 | 297.7 / 297.1 / 265.7 | 286.8 |
| BaseRecalibrator | 8,843.4 | 86.7 / 94.4 / 90.7 | 90.6 |
| ApplyBQSR | 7,242.0 | 340.2 / 404.2 / 406.2 | 383.5 |
| HaplotypeCaller VCF | 5,525.9 | 515.9 / 438.6 / 478.7 | 477.7 |
| HaplotypeCaller GVCF | 6,795.4 | 523.8 / 659.2 / 642.3 | 608.4 |

The spread between these runs (9-26%) is larger than that of the whole pipeline runs (under 4%); the first run was
made on 2026-11-24, the other two in the night of 2026-11-25. The cause was not investigated.

**HaplotypeCaller on the CPU only** (C, three runs each, shards identical to GATK's): VCF 1,836.7 / 1,811.2 /
1,837.9 s, GVCF 1,221.7 / 1,277.6 / 1,259.5 s. VCF mode is consistently about 46% slower than GVCF mode on the CPU;
not investigated.

**Mutect2 alone, inputs and checks.** The pipeline does not keep its recalibrated BAMs, so they were recreated from
the first run's duplicate-marked BAMs and recalibration tables with variant-rs ApplyBQSR, and variant-rs Mutect2 on
them gave the first run's VCF and statistics again. Each Parabricks, GATK and repeat run read the two BAMs (218 GB)
once before its timer started. A single earlier standalone variant-rs run, made straight after the BAMs were
written, took 644.8 s; the three repeats (532.5-535.1 s) and the run with the rewritten resources (536.4 s) show it
was an outlier. The rewritten resources keep 2,536,803 of the panel's 2,609,566 records and 268,169,676 of gnomAD's
268,225,276. GATK Mutect2 as one process: max RSS 8.4 GiB, 805% CPU; FilterMutectCalls 125.3 s.

**CPU baseline, memory:** max RSS of MarkDuplicates 13.5 GiB, BaseRecalibrator 1.5 GiB, ApplyBQSR 0.8 GiB; the
alignment's cgroup peak (page cache included) 49.3 GiB.

**Other load on the server.** Another user's processes were on GPUs 0-1 throughout. An unrelated job of ours used
NUMA node 0 and GPUs 2-3, and GPUs 5-7 during two periods: on 2026-11-25 until 16:50 (start not recorded; it also
used some node-1 threads from 15:07 to 15:28, during the stock sarek somatic run, adding at most about 2% to it) and
during the single-process GATK Mutect2 run, which used no GPU. GPUs 5-7 were empty for the GPU scaling, Parabricks
Mutect2 and Mutect2 repeat runs.

## Earlier configurations

These rows were superseded by the current results in [pipeline/README.md](../../pipeline/README.md#benchmarks). The
machines and inputs are those described there.

Workstation, germline whole genome (NA12878 30x, 12 FASTQ lanes, Broad hg38). The somatic exome took 166 s on
2026-11-19.

| Run | Wall | Steps |
|---|---:|---|
| The same, one run (2026-11-22, at 7d6ffa6) | 34.3 min | ALIGN 23.6 min (peak 32.5 GB), MarkDuplicates 3.5, BaseRecalibrator 1.3, HaplotypeCaller 5.7, hard filters 0.1 (earlier on 2026-11-22: 36.9 and 40.1 min; 2026-11-19: 44.6 min) |

Server, germline whole genome. In the first row, "4" denotes 4 GPUs (4-7) with the 128 threads of their socket. In
the rows that begin with "the same", the configuration is that of the defaults of `run.sh` (6 GPUs, 64 threads per
alignment task, 128 per other step, HaplotypeCaller on 2 GPUs), and these rows and the three after them date from
2026-11-19. The last row is NVIDIA Parabricks 4.8.0-1 `germline`.

| Configuration | Wall | Steps |
|---|---:|---|
| 4, one run (2026-11-22, at 7b3df92) | 9.3 min (555 s) | ALIGN 5.2 min (4 tasks of 3 lanes, side by side; peak ~36 GB each, GPUs 46% busy), MarkDuplicates 1.0, BaseRecalibrator 0.5, HaplotypeCaller 2.3, hard filters 0.1 |
| the same on 2026-11-19 | 7.8 min (465 s) | ALIGN 4.0-4.6 min, MarkDuplicates 0.9, BaseRecalibrator 0.3, HaplotypeCaller 1.7 |
| the same, HaplotypeCaller on 1 GPU | 8.1 min (484 s) | HaplotypeCaller 1.9 |
| 6 GPUs, 64 threads every step | 8.4 min (505 s) | MarkDuplicates 1.2, BaseRecalibrator 0.5, HaplotypeCaller 2.0 |
| 6 GPUs, 64 threads every step, before the HaplotypeCaller GPU engine pool | 9.9 min (593 s) | HaplotypeCaller 3.5 |
| 1 GPU, 64 threads every step, before the engine pool | 20.9 min | ALIGN 15.3 min (12 lanes, one task), MarkDuplicates 1.3, BaseRecalibrator 0.5, HaplotypeCaller 3.5 |
| 4 GPUs, 64 threads (2026-11-18) | 16.3 min | 99.8% of its calls have the same site, alleles and genotype as GATK's, 3% are identical records |

## Server, from FASTQ (2026-11-19)

`~/projects/wgs*/` on the server: the workstation's inputs (12 FASTQ lanes, Broad hg38, same known sites and
calling intervals), GPUs 2-7. Each change was measured on the whole run, VCF data lines identical every time
(md5 741f6ebc50be; filtered 17d3db3ded77):

* 593 s: 6 alignment tasks of 2 lanes (3.8-4.4 min each, ~50 of their 64 threads busy: the 384 threads are shared),
  MarkDuplicates 77 s, BaseRecalibrator 30 s, HaplotypeCaller 211 s on one GPU at ~40% (its 8 contig workers queued
  on one PairHMM engine; more threads made it slower: 212 s at 64, 250 s at 384).
* 505 s: HaplotypeCaller's GPU PairHMM engine pool and 16 contig workers: 211 -> 120 s.
* 484 s: `run.sh`'s defaults, 128 threads for the steps after alignment: MarkDuplicates 74 -> 57 s, BaseRecalibrator
  30 -> 19 s, HaplotypeCaller 120 -> 111 s. The alignment tasks' times vary by ~±10 s between runs (4m06-4m43 here).
* 465 s: HaplotypeCaller's PairHMM on 2 GPUs (`hc_gpus`): 111 -> 103 s.
* 460 s (2026-11-22, `runs/wgs5`, at 98e7200): the PairHMM tile kernel (PLAN.md stage 3) and the align-rs kernels
  of stage 7 (mate SW wavefront, ungapped extension, certified CIGAR band). Alignment 3m51-4m43, HaplotypeCaller
  98 s: the run is CPU-bound (GPUs 19-32% busy), so the faster kernels hardly show.
* 409 s (2026-11-22, `runs/wgs6`, at 401e274): align-rs's bounded dropper queue (`docs/notes/align-rs.md` section
  15). Each alignment task's memory had grown with the batches its one freeing thread was behind on (peak 54-169
  GB a task, the largest the slowest); now 37 GB each, and the six tasks finish together at 3m43-3m47 (before
  3m51-4m43).
* 555 s with 4 GPUs and 128 threads (2026-11-22, `runs/wgs7`, at 7b3df92, stage 8 done): GPUs 4-7 and the 128
  threads of their socket (NUMA node 1, CPUs 96-159 and 288-351 by `taskset`; GPUs 0-1 were another user's), 32
  per alignment task (the workstation's threads per GPU), 128 per other step. 4 alignment tasks of 3 lanes,
  5m13-5m15, ~36 GB each, the GPUs 46% busy; MarkDuplicates 59 s, BaseRecalibrator 29 s, HaplotypeCaller 138 s,
  hard filters 8 s. VCF data lines, recalibration table and duplicate metrics identical to the workstation's
  (`wgs/stage9`). Parabricks `germline` took 16.3 min on 4 of these GPUs with 64 threads, its own outputs; 16.8 min on
  the same GPUs 4-7 and 128 threads (2026-11-22).

One GPU, 64 threads (`~/projects/wgs1`, before the engine pool): 20.9 min, alignment 73% of it (GPU ~75% busy), peak
151 GB for the 12-lane align-rs task.

Somatic exome (`~/projects/seqc2`): 63 s, against 166 s on the workstation (2026-11-22: 109 s held to 4 GPUs and 128
threads, 167 s on the workstation, outputs unchanged). HaplotypeCaller on the normal is the
longest step (41.9 s, 30 s on the workstation): over an exome's targets it runs once (`scatter_intervals` auto), and
one GATK run draws from one random generator in genome order, so its front end reads the BAM in one pass (288 of its
summed stage seconds; the small chromosomes at the end of the file finish last). Its time follows one core's speed:
the same 39.6-41.0 s, VCF identical, with 32 to 128 threads, 8 to 28 contig workers, 8 or 16 GPU engines (GPU 2%
busy). `scatter_intervals true` calls each target on its own (as sarek, and GATK per interval), in parallel.

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
BQSR 20.1-20.8 s, the same as without; output identical ([`tests/harness/regress.sh`](../../tests/harness/regress.sh) all identical).

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
see [docs/notes/variant-rs.md](variant-rs.md)), under 16 busy threads 20% of its time instead of 45%. Whole run estimate: ~12 min.

### From FASTQ with align-rs (2026-11-17)

`wgs/na12878_fastq.yaml`: the 12 lanes as FASTQ (extracted from the BAM by [`tests/harness/wgs/fastq/extract_all.sh`](../../tests/harness/wgs/fastq/extract_all.sh)), Broad
hg38 (`ref_hg38`, which has the bwa-mem2 index), same known sites and intervals. **81 min wall**: 12 ALIGN tasks
69.3 min (5.6-6.0 min each, one at a time: GPU and memory), then (`wgs/na12878_fastq_md.*`, run from the lane BAMs
after a full disk stopped the first run at MarkDuplicates) MarkDuplicates over the 12 lane BAMs 3.5 min,
BaseRecalibrator 1.2, HaplotypeCaller 6.7, hard filters 0.1. VCF 4,951,704 records (4,944,053 after filtering).
Stock tools on this machine, same steps: ~5.2 h (bwa-mem2 + samtools ~117 min, sarek's stock GATK spans ~196 min).

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
`test/somatic_bam` ([pipeline tests](../../pipeline/README.md#tests)) and [`variant-rs/tests/scatter_intervals.sh`](../../variant-rs/tests/scatter_intervals.sh) (duplicate / overlapping / contained /
touching intervals vs GATK per interval; 6,641 exome pieces twice), both in [`tests/harness/regress.sh`](../../tests/harness/regress.sh).
