# Speed plan: one GPU first, several as an option

Most users have one GPU (a workstation, or one cloud GPU), so the stages below speed up that case first. Several
GPUs (large servers, consortium runs) come after, built on the same pieces. Every stage keeps the outputs identical
(the rule in [CONTRIBUTING.md](CONTRIBUTING.md)) and reports its timing and that check before the next one starts.

## Status (2026-11-19)

**Alignment (align-rs): done for now.** What is in: the sorted BAM's final pass in parallel (a 30x lane 25 -> 5 s),
extension's empty redo scan skipped and its own thread pool, and the pipeline's thread split per GPU. What was tried
and not pursued, with the measurements in `docs/notes/align-rs.md` sections 8-11: other extension schedules, GPU
stream priorities, pinned copies, chaining and mate-rescue CPU cuts (stage 2), a one-thread-per-alignment mate-SW
kernel (stage 4). Each machine is now held by something no small change moves:

* one alignment on one L40S: the GPU (its kernels alone ~720 of ~850 ms a batch);
* several alignments on one server: the CPU (6 x 64 threads fill the 384);
* the workstation: the copies between host and GPU through pageable memory, slow while the CPU is busy.

Going further means a larger project: pinned buffers filled and read in place by every GPU stage (the workstation),
or new alignment kernels (one GPU with CPU to spare), each for an uncertain 5-25%.

**Calling (variant-rs):** HaplotypeCaller's GPU engine pool (stage 1) and several GPUs (M1) are in; the PairHMM
kernel (stage 3) and MarkDuplicates / BaseRecalibrator on the workstation (stage 5) are open.

## Where the time goes (2026-11-19)

NA12878 30x, 12 lanes, FASTQ -> VCF (`pipeline/`), `main` at 2e455f7:

| Step | Workstation, 1 RTX 5070 Ti, 32 threads (2026-11-17, before the align-rs work below) | Server, 6 L40S, 64 threads per step |
|---|---:|---:|
| Alignment | 36.1 min (12 lanes in one task) | 3.8-4.4 min (6 tasks of 2 lanes, side by side) |
| MarkDuplicates | 4.1 min | 1.3 min |
| BaseRecalibrator | 1.2 min | 0.5 min |
| HaplotypeCaller | 7.0 min | 3.5 min |
| Total | 48.8 min | 9.9 min |

What limits each step, measured on the server:

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

### Stage 3: PairHMM kernel (if stage 1 leaves HaplotypeCaller GPU-bound)

245 GCUPS on the L40S, 169 on the RTX 5070 Ti; 8.9 M of 1.15 G pairs go to the double-precision pass.

### Stage 4: align-rs kernels — prototyped, not pursued (2026-11-19)

Mate SW first (289 ms a batch on the L40S, ~490 G cells/s), then extension (241 ms). Would help only where alignment
is GPU-bound (one GPU with CPU to spare); not the workstation (copy-bound) or several alignments per server
(CPU-bound). A prototype of the mate-SW forward pass as one thread per alignment over 32-column strips of the query
(exact: all 65,536 jobs of a dumped chunk) was 1.25-1.3x the current kernel (L40S 12.9 vs ~16.5 ms; RTX 5070 Ti
22.8 vs ~26.9 ms): ~60-70 ms of ~720 ms GPU time a batch, perhaps 5-8% of a one-GPU alignment, for a new kernel with
its reverse pass and a fallback. Not pursued; details in `docs/notes/align-rs.md` section 11.

### Stage 5: MarkDuplicates and BaseRecalibrator

5.3 of the workstation's 44.6 min; profile there if they matter after stages 1-2.

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
* **M3, many samples:** the pipeline already runs samples side by side, one GPU per GPU step; to check at cohort
  scale.

Order: alignment done for now (stages 2 and 4 tried and not pursued; M2); stage 1, M1 and M4 done; stages 3 and 5
and M3 open, as measured.
