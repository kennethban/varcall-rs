# Plan: speed (one GPU first, several as an option), then GATK coverage

Most users have one GPU (a workstation, or one cloud GPU), so the stages below speed up that case first. Several
GPUs (large servers, consortium runs) come after, built on the same pieces. Every stage keeps the outputs identical
(the rule in [CONTRIBUTING.md](CONTRIBUTING.md)) and reports its timing and that check before the next one starts.

## Status (2026-11-22)

**Alignment (align-rs): done for now.** Stage 7 (2026-11-22): the mate-SW wavefront, the ungapped extension
shortcut and the certified CIGAR band are in (same records), idle-lane refill measured and dropped; the bounded
dropper queue took the server's germline run 460 -> 409 s. Earlier work in: the sorted BAM's final pass in parallel (a 30x lane 25 -> 5 s),
extension's empty redo scan skipped and its own thread pool, and the pipeline's thread split per GPU. What was tried
and not pursued, with the measurements in `docs/notes/align-rs.md` sections 8-11: other extension schedules, GPU
stream priorities, pinned copies, chaining and mate-rescue CPU cuts (stage 2), a one-thread-per-alignment mate-SW
kernel (stage 4). Each machine is now held by something no small change moves:

* one alignment on one L40S: the GPU (its kernels alone ~720 of ~850 ms a batch);
* several alignments on one server: the CPU (6 x 64 threads fill the 384);
* the workstation: extension's CPU work (setup, applying each round's results, replay), on which its GPU rounds
  already overlap. Pinned buffers filled in place made extension's copies 3-4x faster and left the wall clock as
  it was (stage 8).

Going further: the candidates measured on 2026-11-22 are in [stage 8](#stage-8-candidates-measured-2026-11-22--not-started);
none moves a machine's pace on its own.

**Calling (variant-rs):** HaplotypeCaller's GPU engine pool (stage 1), several GPUs (M1) and a cohort (M3: a 30x
trio, FASTQ -> joint VCF in ~34 min on the server) are in; MarkDuplicates on the workstation (stage 5) 281 -> 246 s;
HaplotypeCaller on the workstation (stage 6) 456 -> 400 s; the PairHMM kernel (stage 3) 2.4-2.8x faster, HaplotypeCaller
on 1 L40S 106 -> 94 s (as fast as 2 GPUs) and on the workstation 400 -> 355 s.

**Coverage:** the GATK tools to add, researched 2026-11-20: [Coverage](#coverage-gatk-tools-to-add-next). C1
VQSR, C2 VariantFiltration, C3 GenomicsDBImport and C4 CreateSomaticPanelOfNormals are done, and the pipeline runs
joint germline calling with them (2026-11-21). No further tool is planned.

## Where the time goes (2026-11-22)

NA12878 30x, 12 lanes, FASTQ -> VCF (`pipeline/`), `main` after stages 3, 5, 6 and 7 (workstation: 98e7200 plus the
dropper fix; server: 401e274); VCF data lines identical in every run (md5 741f6ebc50be, filtered 17d3db3ded77):

| Step | Workstation, 1 RTX 5070 Ti, 32 threads (2026-11-19 stage 0 -> now) | Server, 6 L40S, `run.sh` defaults (2026-11-19 -> now) |
|---|---:|---:|
| Alignment | 32.1 -> 29.7 min (12 lanes in one task; peak RSS 34.6 GB) | 4.0-4.6 -> 3.7-3.8 min (6 tasks of 2 lanes, side by side; 37 GB each) |
| MarkDuplicates | 3.9 -> 3.3 min | 0.9 -> 0.8 min |
| BaseRecalibrator | 1.4 -> 1.25 min | 0.3 -> 0.3 min |
| HaplotypeCaller | 6.9 -> 5.5 min | 1.7 -> 1.7 min |
| Total | 44.6 -> **40.1 min** | 465 -> **409 s** |

The server's run is paced by the CPU (6 alignments on 384 threads; GPUs 23-32% busy on average), so the kernel work
of stages 3 and 7 hardly shows there (465 -> 460 s, `runs/wgs5`); the dropper fix took it to 409 s (`runs/wgs6`:
the alignment tasks no longer grow to 54-169 GB, and all six finish at ~3m45). The workstation gains 10%
(`wgs/stage7`). With its memory held to 30 GB (a cgroup, `MemoryMax=30G`, 16 GB of swap; Nextflow told 30 GB)
the same run took 41.6 min, VCF identical (`wgs/mem30`): alignment 31.0 min (peak 26.9 GB; the hg38 index's pages
are dropped and re-read rather than swapped, swap at most 0.9 GB), BaseRecalibrator 1.4, the rest as before.

What limits each step, measured on the server (2026-11-19; HaplotypeCaller has had its engine pool since):

* **HaplotypeCaller is limited by how the GPU is fed.** Its 8 contig workers share one CUDA PairHMM engine behind
  one lock (`crates/cudahmm/src/lib.rs`): one stream, one pinned buffer, and the packing of a batch on the CPU
  (74 s in all) inside the lock, so it never overlaps another batch's GPU time (109 s; kernel 89 s at 245 GCUPS).
  More threads do not help: 212 s at 64 threads, 218 at 128, 225 at 192, 250 at 384; GPU busy ~40%.
* **align-rs with one GPU is GPU-bound.** Alone, its kernels take ~720 ms a batch (mate SW 289, extension 241,
  seeding 162, SA 30) against ~850 ms per batch in the pipeline (one 30x lane in 87.6 s at 64 threads).
  Scheduling the kernels differently only moves time between them (`docs/notes/align-rs.md` section 8).
* **align-rs on a large server is CPU-bound**: 6 alignments fill the 384 hardware threads (a lane ~125 s against
  88 s alone), so GPU-side gains do not show there.
* MarkDuplicates and BaseRecalibrator scale with threads on a large server (outputs identical at every count):
  MarkDuplicates 75 s at 64 threads, 56 at 128, 55 at 192; BaseRecalibrator 30, 19, 16.

## Stages

### Stage 0: single-GPU baseline (measurement only) — done 2026-11-19

The full run with one GPU (same inputs; VCF data lines identical on both machines and to the 6-GPU run: 4,951,704,
md5 741f6ebc50be):

| Step | Workstation: RTX 5070 Ti, 32 threads, 30 GB | Server: one L40S, 64 threads |
|---|---:|---:|
| Alignment (12 lanes, one task) | 32.1 min (72%): CPU 2,667% of 3,200, GPU 48% | 15.3 min (73%): GPU 75% |
| MarkDuplicates | 3.9 min | 1.3 min |
| BaseRecalibrator | 1.4 min | 0.5 min |
| HaplotypeCaller | 6.9 min (15%) | 3.5 min (17%) |
| Total | 44.6 min | 20.9 min |
| Alignment peak memory | 24.4 GB | 151 GB |

* **On a workstation, alignment is CPU-bound** (83% of the threads busy, the GPU half idle), so the typical
  one-GPU user gains from less CPU work per lane (stage 2), not from faster alignment kernels (stage 4).
* **HaplotypeCaller's PairHMM path behind one lock costs on both:** on the workstation packing + GPU took 250 of its
  453 s (kernel 169 GCUPS; pack 71 s, GPU 167 s), next to heavy CPU stages (front end, assembly), so stage 1 gains
  less there than on the server.
* Workstation folder: `gatk/wgs/stage0/` (run.sh, trace.tsv, gpu.csv); server: `~/projects/wgs1/`.

### Stage 1: HaplotypeCaller, packing overlapped with the GPU (variant-rs `cudahmm`) — done 2026-11-19

A small pool of PairHMM engines, each with its own CUDA stream and pinned buffers, in place of one engine behind one
lock: workers pack their batches in parallel and their kernels overlap. Each (read, haplotype) pair is computed on
its own by the same kernel, so the likelihoods cannot change; checked by `tests/selftest.sh` and the 30x VCF.
Expected on the L40S: 211 s -> ~120 s (the GPU's own 109 s becomes the floor).

Result: a pool of 8 engines (`VARIANT_RS_GPU_STREAMS`), and more contig workers beyond 32 threads (a quarter of the
threads, up to 32). L40S, 64 threads: 190 s -> 121 s; 128 threads: 114 s. RTX 5070 Ti, 32 threads: unchanged
(397-403 s whatever the engines or workers: its CPU stages set the pace), so its defaults are as before. VCF
identical in every run; selftest identical. Details: `docs/notes/variant-rs.md` (Performance).

### Stage 2: align-rs, less CPU work per lane — tried and dropped 2026-11-19

Chaining and mate rescue's dedup were cut (-10.5% user CPU per lane) without changing the workstation's wall
clock: its lane is paced by GPU stages whose calls are mostly copies through pageable memory, slow when the CPU is
busy. The changes were reverted as not worth their code (`docs/notes/align-rs.md` section 10). The structural fix
for such machines, pinned buffers filled and read in place by every GPU stage, is not planned for now.

### Stage 3: PairHMM kernel — done 2026-11-21

Was 245 GCUPS on the L40S, 169 on the RTX 5070 Ti (one warp per pair); 8.9 M of 1.15 G pairs go to the
double-precision pass.

* **Float pass, tiled:** G threads per pair (4-32 lanes), thread t holding read rows t*K+1 .. t*K+K in registers and
  sweeping the haplotype one step behind thread t-1 (gpuPairHMM's layout, Schmidt et al., arXiv 2411.11547), with
  GKL's operations per cell as before. Pairs are binned by read length into 12 tiles (16 .. 256 rows, (G, K) tuned per
  tile on both GPUs); longer reads keep the warp kernel. The 5M-pair dump (1,087,302 pairs, 400 regions a launch),
  every value bit-identical to GKL AVX-512: float bins up to 1.33 TCUPS on the L40S.
* **Double pass:** the redo list is spread one warp per SM (it had 4 warps on each of the first SMs, queued for the
  SM's 2 FP64 units): 15 -> 7 ms on the L40S for the dump's 2,657 pairs. A tiled double kernel was slower on both GPUs
  (13 / 23 ms against 8 / 15) and was dropped.
* **Kernel, the dump:** L40S 282 -> 663 GCUPS, RTX 5070 Ti 191 -> ~420 (summed kernel time).
* **HaplotypeCaller, 30x genome, VCF identical in every run:** server (m1's NA12878, 128 threads, 2 runs each,
  alternating): 1 GPU 103.6 / 108.5 -> 94.8 / 92.6 s (summed kernel time 236 -> 62 s, callers' wait for an engine
  ~1,370 -> ~120 s; one GPU now about as fast as two); 2 GPUs 95.3 / 92.8 -> 91.4 / 98.4 s (the 141 Mb chr2 interval
  sets the time). Workstation (32 threads, 1 run): 400 -> 355 s, kernel 142 -> 61 s, GPU busy 37 -> 20%.
  `hmmcmp` (synthetic reads 30-300 bp, 3 seeds) bit-identical to the CPU port; unit tests and the 5 pipeline tests
  (run fresh, not resumed) identical.

### Stage 4: align-rs kernels — prototyped, not pursued (2026-11-19)

Mate SW first (289 ms a batch on the L40S, ~490 G cells/s), then extension (241 ms). Would help only where alignment
is GPU-bound (one GPU with CPU to spare); not the workstation (copy-bound) or several alignments per server
(CPU-bound). A prototype of the mate-SW forward pass as one thread per alignment over 32-column strips of the query
(exact: all 65,536 jobs of a dumped chunk) was 1.25-1.3x the current kernel (L40S 12.9 vs ~16.5 ms; RTX 5070 Ti
22.8 vs ~26.9 ms): ~60-70 ms of ~720 ms GPU time a batch, perhaps 5-8% of a one-GPU alignment, for a new kernel with
its reverse pass and a fallback. Not pursued; details in `docs/notes/align-rs.md` section 11.

### Stage 5: MarkDuplicates and BaseRecalibrator — done 2026-11-21

5.3 of the workstation's 44.6 min. Profiled there (2026-11-21, 30x, 32 threads, page cache dropped before each
run; outputs identical, `VARIANT_RS_TRACE_CPU=1` prints the CPU by step):

| MarkDuplicates | Total | Pass 1 | Duplicates | Pass 2 | CPU |
|---|---:|---:|---:|---:|---:|
| Before | 281 s | 81 s | 43 s | 157 s | 6,347 s |
| Smaller bucket ends, direct bucket reads, cheaper pass 1 | 255 s | 66 s | 40 s | 148 s | 6,090 s |
| And BGZF level 2 in the pipeline, as sarek's | 246 s | 67 s | 40 s | 139 s | 5,677 s |

The same on the server (L40S machine, 2 x EPYC 9654; m1's 12 lanes, inputs in the page cache, tracing off: on, it
made the 128-thread run 53.5 -> 69.3 s): records and metrics identical to m1's MarkDuplicates BAM in every run.

| MarkDuplicates, server | 64 threads | 128 threads |
|---|---:|---:|
| Before | 77.5 s | 58.5 s, 56.4 s |
| After, level 5 | 71.4 s | 53.5 s |
| After, level 2 (the pipeline) | 60.7 s | 48.1 s |

* **Done:** the on-disk ends packed (35 bytes a fragment, 50 a pair; the buckets 60 -> 43 GB) and read straight into
  their arrays; pass 1's mate lookup without a per-record allocation, the @PG / RG tags found in one scan. The
  pipeline writes the BAM at level 2, as sarek does through the gatk launcher (`-Dsamjdk.compression_level=2`;
  Picard's own default is 5): the same records, the BAM ~7% larger.
* **Left:** deflate is still a third of MarkDuplicates' CPU, file reads and inflate a sixth. mimalloc without
  decommit (`MIMALLOC_PURGE_DECOMMITS=0`) gave another ~4% (237 s), not taken: memory returned later.
  BaseRecalibrator (107 s) is 58% `process_read` (the covariates and the table updates); not changed.

### Stage 6: HaplotypeCaller's CPU stages on the workstation — done 2026-11-21

On the RTX 5070 Ti (32 threads) HaplotypeCaller is CPU-bound, not GPU-bound: GPU busy 32%, ~10,000 CPU-s in a
445 s run. Thread CPU by stage (`VARIANT_RS_TRACE_CPU=1`): assembly 2,674 s, BAM decode + filters 774, prepare 663,
the PairHMM call 650 (mostly WSL GPU driver time; packing 91), genotyping 568, read scatter 552, realignment 389,
finalizeRegion 384. Only ~23 of the 32 threads were busy, but more parallelism made it slower, not faster (more batches
in flight 449 s, deeper front-to-back queue 451, 12 / 16 contig workers 454 / 460, a 48-thread pool 476): the pool
threads woke and slept ~38,000 times a second, and every stage was a parallel pass over the batch with a barrier
after it.

* **Done:** each region's work before the PairHMM (finalizeRegion, assembly, PairHMM inputs) is one task, and so is
  its work after it (likelihoods, realignment): two parallel passes over a batch instead of five. Same calls per
  region, so the same output. 30x genome, 2 runs each, alternating: 456 -> 400 s (-12%), pool context switches
  ~38,000 -> ~15,000/s, ~26 threads busy instead of ~23; VCF identical in every run, unit and pipeline tests
  identical.
* **Server (L40S, m1's NA12878, 2 runs each alternating, VCF identical in every run):** smaller gains where the GPU
  sets the pace: 1 GPU at 64 threads 121.6 -> 112.4 s (-8%), at 128 threads 112.1 -> 108.1 s (-4%); 2 GPUs at 128
  threads (the pipeline's) 98.4 -> 94.5 s (-4%; runs vary by ~7 s).
* **What a GPU could still take:** assembly on chr20 (`VARIANT_RS_ASM_PROFILE=1`) is 60% building the graphs (reads
  threaded k-mer by k-mer into hash tables: CPU work) and 39% haplotypes and their alignment to the reference, nearly
  all the haplotype Smith-Waterman (43 G cells on chr20 with prefix reuse, 116 G without). The haplotype and read
  realignment Smith-Waterman together are ~12-15% of HaplotypeCaller's CPU: integer DP that ports exactly, worth
  ~50 s of the workstation's 400, at the cost of splitting each region's task at a batched GPU call. Not started.

### Stage 7: align-rs kernels from the 2026-11-21 research — done 2026-11-22

A survey of GPU work on BWA-MEM-style alignment (2022-2026) found no GPU port that is both faster than align-rs and
byte-identical to bwa-mem2 (Pham et al. ICS 2023 and G3SA ICS 2025 extend with a static band and no per-row
z-drop; Parabricks documents small mismatches; minibwa and bwa-mem3's faster options are inexact by design). Exact
ideas, in order:

1. **Mate SW on a register-tiled wavefront** (CUDASW++4.0, BMC Bioinformatics 2024, 5 TCUPS on an L40S; the layout
   stage 3 used): G threads per job, K query columns each in registers, thread t one target row behind thread t-1,
   the boundary H / F and the running row maximum passed by shuffle: no lazy-F loop, no per-row vote or reduction.
   bwa's striped semantics (row maximum and E from the values before the lazy-F pass, F restarting at each 16-lane
   segment) are carried as a second F chain. **Done 2026-11-21, the default:** kernel 21 -> 12 ms per 65,536 jobs on
   the L40S (31 -> 18 on the RTX 5070 Ti), every field identical to the striped kernel (a lane's chunk, 600,000
   synthetic jobs); one 30x lane on one L40S 86.5 -> 82.4 s (2 runs each), workstation `mid` 1.60 -> 1.47 s a
   batch; records identical to bwa-mem2, `tests/cmp.sh` and the pipeline tests identical
   (`docs/notes/align-rs.md` section 12).
2. **Ungapped shortcut for extension** (bwa-mem3 PRs #544 / #540, output identical at default parameters): an
   extension whose diagonal has at most one mismatch takes the ungapped walk's result instead of banded SW; must
   keep bwa's strict new-maximum rule, fall back when the mismatches cost more than `-d`, and not cut the w / 2w
   retries short (#540 changed XS). **Done 2026-11-22, the default:** a quarter of the jobs, every result identical
   to the DP (`ALIGN_RS_UNGAPPED=check` on a 30x lane: 229.8 M of 920.8 M jobs); the walk runs on the GPU (on the
   CPU it cost more than it saved on the workstation). L40S, one 30x lane: extension kernel 29.5 -> 25.7 s (-12%),
   the lane ~85 -> ~76 s in the last three of five pairs (noisy: ~10 s between runs). Records identical to bwa-mem2
   (`docs/notes/align-rs.md` section 13).
3. **Certified band for ksw_global2** (bwa-mem3 #575, score / CIGAR / NM / MD unchanged): about half the band
   cells in the CIGAR stage, ~5% of a lane's CPU (helps where alignment is CPU-bound). **Done 2026-11-22, the
   default:** equal-length regions run in the band the ungapped score certifies; workstation `mid`: band cells
   22.4 G -> 10.9 G, CIGAR / NM / MD 140 -> 96 thread-s (user CPU 674 -> 651 s; wall unchanged, the GPU paces it);
   `ALIGN_RS_CIGAR_BAND=check` found no difference; records identical to bwa-mem2, `tests/cmp.sh` and the pipeline
   tests identical (`docs/notes/align-rs.md` section 14). bwa-mem3 #580 (unequal lengths, 1.2 M of the 3.9 M DP
   calls here) would take another ~25% of the cells; not done.
4. **Lanes that finish early take the next job** in the extension kernel (AGAThA, PPoPP 2024): only if its warps
   idle; measure warp efficiency first. **Measured 2026-11-22, not pursued:** the DP's lanes are 53% busy (a 30x
   lane on the L40S), but the kernel's time follows its working lanes, not its warp steps: one batch's jobs with
   only k lanes of each warp given work take x1.03 / x1.10 / x1.31 of the full-warp time at k = 24 / 16 / 8 (x1.33
   / x2 / x4 the warp steps). Packing every lane would save ~9% of the kernel at best, ~2 s a lane on a GPU the
   server does not wait for (`docs/notes/align-rs.md` section 16).
5. Seeding k-mer interval table and read reordering (Pham et al.): a few % of seeding at best. Not planned.

### Stage 8: candidates measured 2026-11-22 — not started

Measured on the workstation, outputs unchanged in every case; details in `docs/notes/align-rs.md` section 17 and
`runs/stage7`:

1. **MarkDuplicates' mate map with foldhash** (or rapidhash) instead of its FNV-1a: 65 -> 35 ns a record on 1.45 M
   chr2 records, one thread (std's SipHash 42 ns: FNV's bits spread badly for the table). ~1% of pass 1; three lines.
2. **ISA-L igzip for BGZF compression** in MarkDuplicates' output (deflate is a third of its CPU): level 1-2
   ~1,200 MB/s against libdeflate level 2's 472 MB/s on BAM records, one thread; the BAM ~11% larger (15.1 against
   13.6% of the raw records). Perhaps 20-30 s of the workstation's 3.3 min, for a C dependency. libdeflate stays
   faster at inflating (2,410 against 1,946 MB/s); igzip inflates FASTQ.gz 1.8x faster than zlib-rs (1,491 against
   841 MB/s), but that runs on its own thread per file, off the critical path.
3. **Extension's CPU side** (setup, apply, replay + dedup: ~16 s of the stage's ~20 s on the `mid` set), which paces
   the workstation's lanes: not yet profiled.
4. **GPU work the machines do not wait on:** bwa-mem3's mate-rescue pruning (a third of the mate-SW cells on a
   dumped chunk) and a warp-uniform column loop in `k_extend` (its lanes read 3.44 local-memory lines a step).
   Only for a GPU-bound setup.
5. **Dropped:** bwa-mem3's certified extension band (x1.08 the kernel time here: the ungapped shortcut has taken
   the jobs it certifies cheaply) and pinned extension buffers (wall clock unchanged).

## Several GPUs (an option, after stage 1)

* **M1, HaplotypeCaller on several GPUs — done 2026-11-19:** the engine pool spread over every visible GPU, in one
  process and one VCF. 30x genome, L40S, 128 threads: 1 GPU 115 s, 2 GPUs 98, 4 GPUs 97 (the 141 Mb chr2 calling
  interval, one sequential GATK run, sets the time); the pipeline gives HaplotypeCaller up to 2 GPUs (`hc_gpus`).
* **M2, alignment on large servers — done for now, with alignment (see Status):** CPU-bound (6 alignments fill the
  384 threads); the alignment threads are split evenly per GPU (`align_cpus`). Stage 2's CPU cuts were not measured
  there before being dropped.
* **M4, threads per step on large servers — done 2026-11-19:** the steps after alignment run one sample at a time,
  so they take the machine's threads (at most 128, `cpus`) and each alignment task its share (`align_cpus`):
  MarkDuplicates 75 -> 56 s, BaseRecalibrator 30 -> 19 s; 6 GPUs 505 -> 484 s.
* **M3, many samples — done 2026-11-21:** the CEPH trio (NA12878, NA12891, NA12892; 30x, 36 FASTQ lanes) with joint
  calling and VQSR on the server (6 L40S): ~34 min FASTQ -> filtered joint VCF (alignment 13.7 min, the GVCFs at
  21 min, the joint steps 13 min); all 18 outputs identical to GATK 4.7 (`tests/harness/m3`). Found on the way: the
  joint steps used ~4 cores (GenomicsDBImport ~2.5 h): the intervals are now imported side by side, and
  GenotypeGVCFs genotypes each interval as its own GATK run (sarek's scatter), side by side; a HaplotypeCaller GVCF
  rule (a call with `<NON_REF>` written 0/0, PLs 0) seen once in the three genomes. Larger cohorts not measured.

Order: stages 1, 3, 5 and 6, M1, M3 and M4 done; stages 2 and 4 tried and not pursued; stage 7 items 1-3 done, item 4
measured and not pursued, item 5 not planned; stage 8 measured, not started.

## Coverage: GATK tools to add next

Researched 2026-11-20 against nf-core/sarek 3.12.0 (`$VARCALL_DATA/sarek/src`), Broad's WARP (develop: JointGenotyping
1.7.4, Whole-Genome / Exome Germline Single Sample 3.3.9 / 3.2.9) and the GATK 4.7.0.0 source. What the pipelines run
by default, and whether GATK's output can be reproduced byte for byte:

| Tool | Default in | Exact port | Effort |
|---|---|---|---|
| VariantRecalibrator + ApplyVQSR | joint germline in sarek (`--joint_germline`, any cohort size) and WARP (`run_vets = false`) | yes | medium-large |
| VariantFiltration, general expressions | WARP joint calling (`ExcessHet > 54.69`); sarek only the hard filters (`variant-rs filter`) | yes, for a whitelisted expression subset | small-medium |
| GenomicsDBImport | joint germline in sarek and WARP; GATK's `mutect2_pon.wdl` | same content (GATK's own bytes vary run to run) | medium |
| CreateSomaticPanelOfNormals | GATK's `mutect2_pon.wdl` only (beta); sarek takes a PoN, Broad's 1000G one by default | yes | small |
| FilterVariantTranches | sarek single-sample, after CNNScoreVariants | yes | small, but little value alone |
| CNNScoreVariants / NVScoreVariants | sarek single-sample (pinned to a GATK 4.5.0.0 image: removed in 4.6.1) | no: PyTorch, GATK's own output varies by hardware, VCF written by pysam | leave to GATK |
| VETS (Extract / Train / ScoreVariantAnnotations) | opt-in in WARP; beta | no: Python scikit-learn model, pickle and HDF5 outputs | leave to GATK |
| CombineGVCFs, SelectVariants, Funcotator, gCNV | not a default step in sarek or WARP | mostly | skip unless asked |

Steps, in order, each validated against GATK 4.7.0.0 before the next:

* **C1, VQSR: done (2026-11-20).** Identical to GATK 4.7.0.0 on 27 outputs (`tests/harness/vqsr`): NA12878 30x
  (5.06M records) and a 1000 Genomes cohort slice (50 of 3202 samples, chr20:5-17 Mb, 500k records, InbreedingCoeff),
  WARP / sarek / best-practice settings, SNP, INDEL and BOTH, ApplyVQSR chained as the pipelines run it. A 30x
  genome's SNP model ~1 min (GATK ~11). Found on the way: htsjdk repairs standard INFO / FORMAT header lines on
  reading (`AD` `Number=.` -> `R`, `src/vcfheader.rs`), now done in ApplyVQSR and MergeVcfs (it only shows on VCFs GATK
  did not write, such as the 1000 Genomes calls). The pipeline runs it in joint calling (`joint_germline`,
  2026-11-21). Not done: `-AS`, `--aggregate`, model in / out, scattered tranches. The plan as written:
  VariantRecalibrator (~3.4k Java lines: VariantDataManager, GaussianMixtureModel,
  MultivariateGaussian, VariantRecalibratorEngine, TrancheManager, tranches; plus Jama's LU / inverse / det) and
  ApplyVQSR (small). Deterministic: sequential, one `java.util.Random(47382911)` (`Utils.java:52`) for the whole run.
  The risk is that stream: the datums' number and order (multi-input merging, same-start queueing, AS mode, training
  overlaps), the decode-time jitter, the k-means initialization and the marginalization draws must all happen in
  GATK's order, or every later draw differs. Exact pieces needed: Jama LU, commons-math `digamma`, jmath,
  `DoubleStream.sum` (Kahan), `%.4f` VQSLOD and tranche formatting, stable sorts. Refuse `--rscript-file` (its PDF
  is not reproducible). Test data: a cohort large enough for VQSR (WGS, or tens of exomes, e.g. 1000 Genomes /
  GIAB trios) with the hg38 HapMap, Omni, 1000G and dbSNP resources, SNP and INDEL modes as sarek and WARP run them.
* **C2, VariantFiltration: done (2026-11-20).** `src/vfilter/` (JEXL 2.1.1's semantics in `jexl.rs`), plus Picard
  MakeSitesOnlyVcf, which WARP runs next. Identical to GATK on `tests/harness/vfilter` (NA12878 30x, the 1000 Genomes
  cohort slice, 11-sample joint calls; WARP's ExcessHet / exome mask / DRAGEN filters, GATK's examples, the
  language's edge cases; GATK's own failures reproduced). Found on the way: htsjdk re-encodes every genotype of a VCF
  whose samples are not sorted (refused now, also in ApplyVQSR and MergeVcfs); `KEY=` is written `KEY=.`. Not done:
  genotype filters, SNP clusters, `-L`; `variant-rs filter` stays as it is. The plan as written:
  Extend `filter.rs`: identifiers, numeric literals, comparisons,
  `&& || !`, parentheses, JEXL's string-to-number coercion, and a missing attribute aborting the whole expression as
  htsjdk's `JEXLMap` does; genotype filters (`GQ`, `DP`) through the genotype context. Refuse method calls
  (`vc.isSNP()`) and multi-valued attributes. Covers WARP's ExcessHet and DRAGEN QUAL filters.
* **C3, GenomicsDBImport: done (2026-11-21).** `src/gdb_import.rs` (cells, vidmap / callset JSON, arrays, batches),
  `src/gdb_write.rs` (TileDB schema and fragments). On the 10 GATK workspaces of `tests/harness/genomicsdb` (WARP's
  sample-map imports, batch 1 and 50, three arrays, an array starting mid-contig, sarek's `-V` + BED form with and
  without `--merge-input-intervals`, `--bypass-feature-reader`): every cell identical, `vidmap.json` and
  `callset.json` byte-identical, the header identical but for the command line; GATK's and variant-rs's
  GenotypeGVCFs on them give GATK's output. Specs: `docs/specs/genomicsdb_import.md`. The pipeline runs it in joint
  calling (`joint_germline`, 2026-11-21). Not done: incremental import, `--merge-contigs-into-num-partitions`,
  filtered or multi-sample GVCFs. The plan as written: A workspace whose GenotypeGVCFs output (GATK's and variant-rs's)
  equals that from GATK's workspace: htsjdk parse, BCF2 round trip (floats to f32), the vidmap combine operations,
  callset order, fragments per batch. Or skip the workspace inside the pipeline and feed `gdb_merge` directly. Check
  the clipping of records that start before an interval empirically.
* **C4, CreateSomaticPanelOfNormals: done (2026-11-21).** `src/mutect/pon.rs`. GATK's `mutect2_pon.wdl` chain on the
  four DREAM normals (`tests/harness/somatic/pon`): Mutect2 `--max-mnp-distance 0`, GenomicsDBImport (Mutect2 VCFs,
  flag INFO fields: one char, 1 or 0x7F; one and three arrays: every cell identical), the records a walker receives
  from the workspace (alleles and AD vs GATK SelectVariants, 2170), and the panel with five argument sets: all
  identical. The whole-workspace read (no -L, arrays in `ChrArrayFolderComparator` order) is shared with
  GenotypeGVCFs. Not done: `-L`, a VCF input. The plan as written: Small: the commons-math binomial /
  beta-binomial and Brent fit already ported for FilterMutectCalls (`src/mutect/fmc_math.rs`). Validate the GenomicsDB merge on Mutect2 VCFs (no `<NON_REF>`,
  no-call samples without AD) first.

With C1-C3 the pipeline runs joint germline calling as sarek does (`joint_germline`: GVCFs, GenomicsDBImport,
GenotypeGVCFs, VQSR; identical to GATK in `pipeline/test`, test `joint`); C4 completes GATK's PoN workflow.
Also seen, not planned: WARP's default Picard QC metrics (CollectWgsMetrics, CollectHsMetrics,
CollectMultipleMetrics), medium; sarek uses mosdepth and samtools instead.
