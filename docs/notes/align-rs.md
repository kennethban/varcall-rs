# align-rs development notes

How align-rs was built up stage by stage against bwa-mem2 2.2.1, what was measured, and what was kept or dropped,
in the order it happened (November 2026). Paths are relative to `align-rs/` unless they start with
`tests/harness/` or `$VARCALL_DATA` (the harness data, see `tests/harness/README.md`). For the current picture see
[architecture.md](../architecture.md) and [optimizations.md](../optimizations.md).

Goal: `bwa-mem2 mem -K 100000000 -Y -R RG ref.fa R1.fq.gz R2.fq.gz` output, record for record, as fast as the GPU
allows (the way fq2bam accelerates BWA-MEM, but with bwa-mem2's exact output), writing BAM directly.

## Baseline (2026-11-16, 9950X3D 16 cores / 32 threads, 30 GB RAM, RTX 5070 Ti 16 GB)

NA12878 30x, one lane (31.5M pairs, 1/12 of the genome) against `ref_hg38/Homo_sapiens_assembly38.fasta`:
bwa-mem2 + samtools sort 549 s (13.4 min before the sort memory fix: 8 GB of sort buffers swapped bwa-mem2's index
out), so ~1.8 h for the genome. bwa-mem2's own split, per thread: SMEM 214 s, SA lookup 63 s, banded SW 118 s,
SAM processing (pairing, CIGAR, MAPQ, text) 149 s, peak RSS 21 GB.

## Structure of bwa-mem2 (src/) and the port

Per batch of 100M bases (`-K`, fixed regardless of threads, so results do not depend on the thread count):

| bwa-mem2 | what | where it should run |
|---|---|---|
| `worker_bwt` -> `mem_kernel1_core` | SMEMs (`mem_collect_smem`: forward/backward extension, re-seeding, 3rd round), SA lookups, `mem_chain_seeds`, `mem_chain_flt`, `mem_flt_chained_seeds` | SMEM + SA: GPU; chaining: CPU |
| `worker_aln` -> `mem_kernel2_core` | seed extension, banded SW (`bandedSWA`, = `ksw_extend2`), `mem_sort_dedup_patch` | SW: GPU |
| `mem_pestat` | insert-size distribution of the batch | CPU |
| `worker_sam` | `mem_sam_pe`: mate rescue (`mem_matesw`, `ksw_align2`), pairing, primary marking, `mem_reg2aln` (global alignment for CIGAR, NM/MD), MAPQ, XA/SA, SAM text | CPU (mate rescue SW: GPU later) |

Index (`ref.fa.*`): `.bwt.2bit.64` (BWT checkpoints, 64 bytes per 64 positions, + compressed SA, sampled every 8),
`.0123` (reference, one byte per base), `.pac`, `.ann`, `.amb`, `.alt`. Memory-mapped, so repeated test runs reuse
the page cache.

## GPUs and settings

`make -C cuda` builds kernels for sm_75 (Turing: T4, RTX 20xx) through sm_120 (Blackwell: RTX 50xx), plus
sm_120 PTX that the driver compiles for newer GPUs (about 2.5 min; one architecture only, faster:
`make -C cuda ARCH="-gencode arch=compute_120,code=sm_120"`). A GPU the build has no kernels for, or none at all,
runs everything on the CPU: same output, much slower; align-rs says which at start. The warp reductions use
`redux.sync` from sm_80 and shuffles before (checked on sm_120 with `-DALN_NO_REDUX`: mate-SW results identical).

The settings below spread the same work differently; **no setting changes any output**. Defaults come from the
device's properties and are printed at start, e.g.
`[align-rs-gpu] NVIDIA GeForce RTX 5070 Ti (sm_120, 70 SMs, ~128 KB L1/SM, 1536 threads/SM): seeding 4 lanes per read,
extension 6 blocks/SM, mate SW 56 registers`.

| environment | what | default | measured on RTX 5070 Ti (sm_120) | elsewhere |
|---|---|---|---|---|
| `ALIGN_RS_SMEM_G` | seeding: lanes per read (1, 2, 4, 8, 16, 32; 1 = a thread per read) | 4 | 1M pairs: 1: 1.22 s, 2: 0.79, 4: 0.68, 8: 0.68, 16: 0.99, 32: 1.60 | depends on the reads more than the GPU: 4 or 8 |
| `ALIGN_RS_EXT_BPSM` | extension: blocks of 128 per SM (0: a block per 128 jobs) | (L1 KB) / 21, L1 ~ max shared per SM + 28 KB | 1/6 lane: 4: 5.45 s, 6: 5.06, 8: 5.07, 10: 5.68, 0: 6.14 | sized to L1: ~6 on 128 KB (Ada, consumer Blackwell), ~9 on 192 KB (A100), ~12 on 256 KB (H100, B200) |
| `ALIGN_RS_MSW_MINB` | mate SW, striped kernel (`ALIGN_RS_MSW=striped`): blocks of 64 per SM the registers are capped for (12, 18, 24); the wavefront kernel is fixed at 12 | 3/4 of the SM's threads / 64 | lane chunk: 18 (56 registers) 306 GCUPS, uncapped (161) 164 | 1536 threads/SM -> 18; 2048 (A100, H100) -> 24; 1024 (Turing) -> 12 |

Only the sm_120 values are measured; the rules for other GPUs follow what limited each kernel here (extension:
its DP rows in local memory, so L1 size; mate SW: resident warps, so the register file and thread slots).

**`align-rs tune`** measures them on the GPU at hand (~35 s, plus the index load):

```
align-rs tune -t 32 -K 100000000 ref.fa R1.fq.gz R2.fq.gz
```

It aligns the first batch of the reads (output discarded) while capturing what the GPU is given (the seeding
batch, every extension round's jobs, the mate-SW jobs), then replays those inputs under each candidate value,
one call at a time, timing kernels only (3 runs each, the fastest kept): seeding 1 / 2 / 4 / 8 / 16 lanes,
extension 0 / 2-8 / 10 / 12 / 16 blocks per SM, mate SW 12 / 18 / 24. Every candidate's results are compared
with the first's, and the tune fails if any differ. A value replaces the device default only if it is 3% faster
(run-to-run spread here is 1-3%). The result goes to `~/.cache/align-rs/gpu-settings.tsv` (or
`$XDG_CACHE_HOME/align-rs/`, or `$ALIGN_RS_TUNE_CACHE`), one line per GPU model and architecture; `align-rs mem` applies
it at start and says so (with a note if the driver changed since). The environment overrides it. Re-run after
updating align-rs or the driver. On the RTX 5070 Ti the tune keeps the defaults (4 / 6 / 18): extension 5-8 blocks
per SM are within 3% of each other.

By hand, kernel times one kernel at a time (each setting's runs give identical output):

```
ALIGN_RS_PROFILE=1 ALIGN_RS_GPU_SERIAL=1 ALIGN_RS_EXT_BPSM=8 align-rs mem ...   # "[align-rs-gpu] kernels only: extension ... s"
tools/timerun.sh OUT NAME DIGEST R1 R2                          # end to end: s/batch, records vs digests.txt
```

and the mate-SW kernel alone on a dumped chunk (`ALIGN_RS_DUMP_MSW=FILE align-rs mem ...`, then `cuda/bench_msw.cu`,
which also checks every result field). Measured under WSL2, where pinned staging buffers were slower than pageable
copies and destroying CUDA events blocked; native Linux may differ.

When the device fails (typically out of memory: another process took it, or a very large `-K`), the work of that
call is redone on the CPU, with identical results; align-rs says so once per stage. `tests/gpufail.sh` forces it
(`ALIGN_RS_GPU_FAIL_GROW=N`: device allocations after the first N fail).

## Verification

* Truth: bwa-mem2 2.2.1 (conda build) on `$VARCALL_DATA/wgs/fastq/sub/{tiny,small}` (20k / 1M pairs, first pairs of lane
  HCLHLDSXX_L001; truth SAMs beside them), then the lane. bwa-mem2 is deterministic (rerun identical).
* Edge inputs vs the bwa-mem2 binary (adversarial review, [validation.md](../validation.md)): `tests/badfq.sh` (FASTQ read
  as kseq does, byte for byte: malformed records dropped and the batch ended there, the input ended at the second
  empty batch as bwa-mem2's two pipeline workers do, BOM / junk before '@', spaces in sequence lines, CRLF, Latin-1
  names and qualities copied verbatim, gzip with trailing bytes), `tests/alt.sh` (.alt parsing), `tests/gpufail.sh`.
  `--lanes` refuses two lanes with one output file before doing any work.
* Stage by stage: bwa-mem2 built from `$VARCALL_DATA/tools/bwa-mem2-2.2.1` with dump hooks (SMEMs, chains, regions after
  extension, pestat, final records), compared with the port's dumps.

## Milestones

1. Index loading, FASTQ batching (kseq semantics), SMEMs + SA positions identical. **Done 2026-11-16**: tiny
   2,672,071 and small (3 batches) 133,831,760 dump lines identical to bwa-mem2's (`--dump` vs `BWA_DUMP`). Notes:
   the SA lookup that reaches the sentinel returns 0 (`call_one_step`), and a 512-read chunk with a single seed gets
   no SA lookups or chains (`mem_chain_seeds` loops while pos < num_smem - 1); both reproduced.
2. Chaining and chain filtering identical. **Done 2026-11-16**: small (1M pairs) 20,896,857 chains and 26,970,816
   chained seeds identical (rid, pos, ALT, weight, kept, first, frac_rep bits, seeds), tiny likewise. klib's B-tree
   (`kbtree.h`, t = 5 from sizeof(mem_chain_t) = 48) and `ks_introsort` are ported as such: equal chain positions
   and equal weights come out in their order. `mem_flt_chained_seeds` is a no-op below ~700 bp reads (min_l > 0.05
   l_query); its seed SW (`mem_seed_sw`) is not ported yet and stops the run if needed.
3. Extension (scalar `ksw_extend2` first), regions identical. **Done 2026-11-16**: small 20,846,330 and tiny 417,684
   regions identical after the skip replay, the purge and `mem_sort_dedup_patch` (rid, rb, re, qb, qe, score,
   truesc, w, seedcov, seedlen0, n_comp, ALT, frac_rep); synthetic indel pairs (`tests/harness/wgs/fastq/sub/indel.*`, 2,822
   pairs, deletions 20-400 bp, insertions up to 60) likewise. bwa-mem2 extends every seed (SIMD, 8- and 16-bit
   kernels) and then replays bwa's skip test; the scalar `ksw_extend2` port gives the same results here.
   Not exercised: `mem_patch_reg` merges (none in these sets: the relative-bandwidth limit, 5% of the span, keeps
   150 bp reads from ever merging; long reads would) and `mem_seed_sw`.
4. Insert-size stats, mate rescue, pairing, SAM records identical (tiny, small, lane). **Done 2026-11-16**: SAM
   byte-identical to bwa-mem2 (header apart from @PG) on tiny (40,523 records), small (2,027,515, 3 batches),
   the synthetic indel set (6,824, 1,180 supplementary) and **the full lane: 63,787,017 records identical** after the
   pipeline's `samtools sort`. Notes: bwa-mem2's batched AVX-512 pairing (`mem_sam_pe_batch_*`, SIMD `kswv` mate SW)
   gives the same results as bwa's sequential `mem_sam_pe` + `ksw_align2`, which is what is ported (the striped
   u8 / i16 SW emulated lane by lane); `mapQ_coef_fac` is an `int` in mem_opt_t, so log(50) = 3.91 is used as 3;
   `mem_pair`'s pair id is an `int` and `id<<8` wraps at 32 bits before the 64-bit hash.
   Lane time: align-rs 20.1 min (28 busy threads) vs bwa-mem2 9.2 min, same sort; nothing tuned yet.
5. **Sorted BAM + index output. Done 2026-11-17** (`-o x.bam [-l level]`, `src/bam.rs`, `hts/alnhts.c` against the
   conda env's htslib 1.25): each SAM line is parsed by htslib's `sam_parse1` (samtools' own parser); per batch
   the records go to 1024 buckets of the linear genome (unplaced last), deflated into one spill file
   (`x.bam.spill.tmp`); at the end each bucket is stably sorted by (tid, pos, strand), samtools sort's order, and
   written through htslib's BGZF writer with `sam_idx_init` (as `samtools sort --write-index`). Lane: **records
   identical to bwa-mem2 | samtools sort, `.bai` byte-identical to `samtools index`**; final pass 27 s.
6. GPU kernels: SMEM/SA, banded SW; then mate rescue. **In progress 2026-11-16**:
   * Profile (1M pairs, thread-seconds, CPU port): seeding 159, SA + chaining 145, extension 495, pairing + SAM
     390 (mate SW 279 of it, CIGAR 24, text 5); bwa-mem2's equivalents ~217 / ~64 / ~120 / ~151.
   * **Seed extension on the GPU** (`cuda/aln_ext.cu`, `--gpu`): extension restructured into batch-wide rounds as
     bwa-mem2 runs them (left at w, left retried at 2w, right with the left score as h0, retried), jobs run one
     per CUDA thread with the scalarBandedSWA recurrence (H and E packed in one word, jobs sorted by size), the
     2-bit `.pac` resident on the device (checked against `.0123` at start), reads uploaded once per batch.
     Kernel 1.75-1.86 s for 29.3M jobs (one-thread-per-job first try: 18.6 s); extension phase 18.8 s on 32 CPU
     threads -> 6.6-7.0 s. SAM identical (tiny, small, indel).
   * Mate SW (`ksw_u8` / `ksw_i16`) on the CPU with the same SSE2 intrinsics as ksw.cpp (the lane-by-lane
     emulation stays as `ALIGN_RS_KSW_EMULATE=1`): 279 -> 127 thread-s, SAM identical.
   * 1M pairs end to end: 28.7 s (bwa-mem2 26 s).
   * **Seeding and SA lookups on the GPU** (`k_smem`: one thread per read runs all three seeding rounds on the
     resident FM-index checkpoints; `k_sa`: one thread per SA position walks LF to a sampled entry, the sentinel
     giving 0 as call_one_step does). Reads over 256 bases or with more than 64 seeds are redone on the CPU.
     SAM identical (tiny, small, indel; CPU-only path likewise). 1M pairs: seeding 1.33 s and SA 0.36 s on the
     device; **25.5 s end to end** (CPU only 39.8 s, bwa-mem2 26 s); the GPU is idle most of that time.
   * Device memory tiers (`gpu.rs`, chosen at start from free memory, 3 GB kept for per-batch buffers): .pac
     (0.8 GB): extension; + FM checkpoints (6.4 GB): seeding; + SA samples (4.0 GB): SA lookups. GRCh38 on the
     16 GB RTX 5070 Ti: all three, 2.8 GB left. A larger card would also hold `.0123` (6.4 GB, one byte per base)
     and room for several batches in flight; the batch size stays `-K` (bwa-mem2's output depends on it).
   * **Pipeline** (`main.rs`): five stages, one batch each, handed over without buffering (sync_channel(0)):
     read + gunzip; seeds + SA + chains; extension; insert sizes + mate rescue (GPU); pairing + SAM (+ BAM
     buckets). Pairing/SAM and rescue share their own rayon pool (in the global pool one stage starved the
     other); large structures are freed on a dropper thread; per-stage timeline with `ALIGN_RS_PROFILE=1`.
   * **Mate SW on the GPU** (`k_matesw`: a half-warp per job emulating the SSE lanes, u8 jobs; results looked up
     by sam_pe as mem_matesw_batch_pre does).
   * **Memory**: the lane first ran 741 s, swapping (align-rs 26 GB + samtools sort 5 GB on 30 GB). mimalloc,
     unbuffered stage hand-off and shrinking regions after dedup took the 5M-pair anonymous peak 17.9 -> 12 GB.
   * **Lane (2026-11-17)**: `| samtools sort` 366 s; **`-o lane.bam` (sorted BAM + .bai) 357 s** vs bwa-mem2 +
     samtools sort 549 s + samtools index 17 s; 19 busy threads, GPU ~58% busy. Per batch: seeds + chains 2.9 s,
     extension 3.2 s, rescue 3.4 s (2.1 s GPU), SAM 1.3 s.
7. **Speed, measured first (2026-11-17)**. `ALIGN_RS_PROFILE=1` now also reports thread CPU time per piece of work
   (`prof.rs`: per-thread counters of CLOCK_THREAD_CPUTIME, which other stages' threads do not inflate), counters,
   and kernel-only GPU times (`ALIGN_RS_GPU_SERIAL=1`: one GPU call at a time, each kernel's own cost). Stack samples
   with gdb (`ALIGN_RS_PTRACER=1` lets a sampler attach under yama's ptrace_scope 1) located the hot lines.
   * Measured, then dropped: running only the extensions bwa's replay keeps (it keeps 78% of regions, 73% of
     jobs: little to gain); de-duplicating SA lookups (35% repeats, but the SA kernel is 2 s per lane); pinned
     staging buffers and a high-priority seeding stream (no gain); concurrent vs serialized kernels (equal).
   * Kept (all outputs identical: tests/cmp.sh, a 30x lane's records vs bwa-mem2):
     - `--lanes FILE`: a sample's lanes in one process (index loaded once, each lane's final BAM pass overlaps the
       next lane); each lane's BAM and .bai byte-identical to its single-lane run (in tests/cmp.sh). The
       pipeline's ALIGN is one task per sample (germline test 36 -> 22 s, somatic 177 -> 103 s).
     - region sorts (`a_idx_sort`): a fast sort, with ks_introsort only when regions with equal keys differ
       (0.4% of sorts); regions are `Copy`. Mate rescue's CPU part (no SW: its sort_dedup per rescue hit) 218 ->
       125 thread-s on 1/6 of a lane.
     - chaining: chain_flt's pairwise test on flat arrays, the B-tree with keys' positions inline (188 -> 141).
     - extension jobs ordered by size on the GPU (CUB radix sort; the host sort and scatter of 10M jobs per round
       took ~8% of all CPU); capacities reserved in setup; CUDA blocking sync (threads spun in cudaMemcpy).
     - mate SW kernel: slen as a template parameter (registers, not local memory), ints instead of packed bytes,
       the score from bwa's matrix form, no u8 saturation where H cannot reach 255, an early exit from the
       lazy-F loop: 19.7 -> 11.9 s on 1/6 of a lane (bwa's striped semantics kept: the row maximum before the
       lazy-F pass, E not updated by it).
     - seeding output packed on the device (CUB scan) instead of 64 fixed slots per read, cap 256: reads redone
       on the CPU 1.1% -> 0.1%.
     - rescue its own rayon pool (its job building waited for the previous batch's SAM text); seeding and SA +
       chaining two pipeline stages.
   * **Lane: 2.33 s per batch (was 3.35), alignment done at 236 s (330 s), + final BAM pass 28 s**; 4,420 s user
     CPU (5,717), peak RSS 24.3 GB. Now paced by extension (2.8 s per batch: host-side setup and replay).
   * Extension: the replay's skip test looks up kept regions by position (whether any contains the seed is what
     bwa's ordered scan decides: 37.9 -> 15.1 thread-s on 1/6 of a lane); fetch_window takes the chain's rid
     (pos2rid on packed offsets); the batch goes through in two halves, the CPU setting up / applying /
     replaying one while the device runs the other's round.
   * Memory reuse: mimalloc purges freed memory with MADV_FREE instead of decommitting (`ALIGN_RS_MI_DECOMMIT=1`
     restores it): the next batch reuses the pages unless the kernel reclaims them. 1/6 of a lane: minor faults
     10.6M -> 2.0M, 2.25 -> 2.13 s per batch; anonymous RSS peak 12.2 -> 14.4 GB, the difference reclaimable.
   * **Lane: 2.27 s per batch, alignment done at 221 s** (first batch at 8 s, index in the page cache), system
     time 624 -> 337 s, user 4,348 s; records identical to bwa-mem2. 2,060% CPU: the stages, not the core
     count, now set the pace.
   * CIGAR: ksw_global2 on per-thread scratch with a bounds-check-free band loop, bns_get_seq sized up front,
     numbers written into the SAM / MD / XA text without a String each: CIGAR 115 -> 93 thread-s, SAM pool
     379 -> 343 on 1/6 of a lane; wall unchanged (2.11-2.15 s per batch): the SAM stage does not set the pace.
   * Then GPU-bound: kernels 1.65 s per 2.1 s batch (mate SW 0.80, seeding 0.46, extension 0.36).
   * Mate SW kernel, tuned on a dumped lane chunk (`ALIGN_RS_DUMP_MSW=FILE`; `cuda/bench_msw.cu` runs the kernel on it
     and checks every result field): it waits on its per-row shuffles and votes, so occupancy decides. 64-thread
     blocks capped at 56 registers (some spilling in long-query variants), and qe chosen when the row max
     improves instead of saving the row: 164 -> 306 GCUPS, results identical. Nsight Compute cannot read the
     counters under WSL without the Windows-side permission (ERR_NVGPUCTRPERM).
   * Seeding kernel: register caps / block sizes change nothing (1.19-1.33 s on 1M pairs): it waits on dependent
     FM-index reads; not on the critical path now (seed stage 1.7 s, extension 2.15 s).
   * **Lane: 1.93 s per batch, alignment done at 189 s** (from 357 s this morning, records identical to
     bwa-mem2); 2,367% CPU. Paced by extension (2.15 s) and SA + chaining (2.0 s).
   * Timing tools in `tools/` (`timerun.sh`: steady batch rate, RSS, faults, records vs `digests.txt`), since
     /tmp does not survive a WSL restart. `ALIGN_RS_PROFILE=1` is now cheap (timeline, counters); the per-item CPU
     timers (a system call each) are `ALIGN_RS_PROFILE=2`.
   * Each batch's reads go to the device once (keyed buffers shared by seeding, extension and mate SW, instead
     of three uploads). Stage slack (`ALIGN_RS_SLACK*`, batches a stage may run ahead): no gain on the lane (1.91 /
     1.92 / 1.96 s per batch for none / read + seed / all) and up to 4 GB more memory: default none.
   * Extension kernel (Nsight Compute: L2-bound at 94%, the per-thread DP row in local memory): the query base
     packed into the H/E word, one local load and store per cell instead of three accesses: 5.07 -> 4.67 s on
     1/6 of a lane. Jobs whose scores could reach 2^14 go to the CPU (none on the lane).
   * The device is now the shared limit: seeding, SA lookup, extension and rescue keep ~2.3 GPU calls in flight;
     the SA + chaining stage waits on the device, not on its own work. Largest GPU cost, seeding: 3.2 of 32
     threads active per warp (one thread per read; a few repetitive reads per warp run ~10x longer than the rest),
     60% of stalls in the occurrence lookups; 69% of backward_ext calls are in the backward phase, independent
     across candidates (forward walk 17%, round 3 13%).
   * **Lane: 1.91 s per batch** (alignment ~185 s), records identical to bwa-mem2.
   * **Seeding with a group of lanes per read** (`k_smem_g<G>`, `ALIGN_RS_SMEM_G`, default 4; 1 = the former one
     thread per read): forward walks and rounds 2 / 3 in lockstep on the group's lanes; the backward phase's
     candidates extended G at a time, then every lane replays bwa's selection in order on shuffled values (first
     loop's break, curr_s with its int32 truncation, in-place compaction), so the seeds are one_pos's. Persistent
     groups take reads from a counter (a slow repetitive read no longer holds 31 other lanes). 1M pairs: seeding
     kernel 1.22 -> 0.68 s (G = 2: 0.79, 8: 0.68, 16: 0.99, 32: 1.60); SAM identical at every G.
   * FASTQ decompression on a thread per file (4 MB blocks to the parser): read stage 1.62 -> 1.05 s per batch.
   * **Lane: 1.80 s per batch, alignment done at ~183 s**, records identical to bwa-mem2. Now paced by extension
     (1.7 s: ~0.9 s in GPU calls, half of it queueing, ~1.2 s CPU in two overlapping halves).
   * Extension round: 24-byte jobs (band, end bonus and direction per round, as kernel arguments; upload
     2.41 -> 1.81 s on 1/6 of a lane); the kernel grid-strided over 6 blocks per SM (`ALIGN_RS_EXT_BPSM`: with fewer
     registers more threads were resident and their DP rows no longer fit L1: 6.92 -> 5.06 s); timing events
     created once per thread (`cudaEventDestroy` blocked ~8 ms per call under WSL, 13.7 s of stage time).
     Pinned staging measured again with a timeline: slower (1.94-2.03 vs 1.73 s per batch), removed.
     Nsight Systems: the GPU executes kernels 41% of the time; pageable copies run at ~1.7 GB/s H2D in the
     pipeline. **Lane: 1.80-1.82 s per batch, unchanged**: the stages are within ~10% of each other and the
     remaining GPU waits are copies and driver calls under WSL.
   * **Bug found by the 12-lane run (fixed)**: the replay's position lookup panicked when a seed was longer than
     every kept region (its window [seed end - longest region, seed start] empty, lo past hi; no region can then
     contain the seed). It never occurred in tiny / small / mid / lane L001, first in lane L002. Checked on L002
     with `ALIGN_RS_REPLAY_CHECK=1` (the lookup and bwa's linear scan run side by side and must agree on every seed).
     The crash also showed that a panicking thread let the other stages finalize BAMs (lane L002 cut short at 18M
     records, the lanes not reached written empty) before the error exit: release builds now abort on panic.

8. **Extension on a many-core server (2026-11-19, L40S, 64 of 384 threads, one 30x lane)**. A batch every 900 ms,
   extension busy 880 of them, mate rescue 836. `ALIGN_RS_TRACE` (one line per device call, and then per extension CPU step: those events were removed once
   used) showed:
   * the device thread's scan for jobs the kernel left to the CPU cost 80 ms a batch (a parallel loop waiting for
     the thread pool) and found none: the kernel now counts them and the scan runs only when there are some;
   * the first half's setup took 240 ms with the device idle, the same work 67 ms for the second half: it starts
     when the other stages start their next batch, behind SA positions' and chaining's loops in the global pool.
     Extension's CPU steps now have their own pool, as SAM and mate rescue already did (the device calls stay
     outside it): first setup 240 -> 203 ms (the cores are still shared with the other stages);
   * extension's device calls overlap mate rescue's for 470 of their 503 ms: extension kernel 36.2 s on the lane
     against 22.9 s with `ALIGN_RS_GPU_SERIAL=1`, mate SW 35.8 against 27.4 s.

   Three runs each, alternating, records identical: 64 threads 113.6 -> **108.5 s** (-4.5%), 16 threads 197.2 ->
   194.1 s. Now extension (841 ms) and rescue (836 ms) pace the batches equally; the final BAM pass is 25 s of it.
   Tried and dropped: extension in K parts of a batch on their own threads, one device thread serving their rounds
   as they come (no lockstep wait): 2 parts as fast as the halves, 4 and 8 slower (113.1 / 117.5 / 120.6 s against
   114.1 s): smaller device calls, and more CPU loops competing with the other stages.

   Mate rescue's device time, traced per chunk (64K jobs; ~898K jobs, 143 G cells a batch): 555 ms of calls a
   batch, kernel 370 (293 alone), result download 175 (9 alone: median 0.4 ms, one chunk in ten 55-157 ms, during
   extension's and SA lookups' large pageable copies). Tried and dropped, records identical:
   * results through a pinned buffer: download 175 -> 1 ms, rescue 832 -> 599 ms a batch, but its kernels then
     take the device from extension (808 -> 910 ms): 107.9 -> 114.2 s;
   * extension on a highest-priority stream: its kernel slower (320 -> 423 ms a batch), seeding too: 115.5 s
     (with pinned results 117.1 s).

   The device is full: alone, the kernels take ~720 ms a batch (mate SW 289, extension 241, seeding 162, SA 30)
   against 848 ms per batch in the pipeline, where they overlap and slow each other. Scheduling only moves time
   between the stages; what is left is less device work (the mate-SW kernel: ~490 G cells/s on the L40S).

9. **Final pass of the sorted BAM, in parallel (2026-11-19)**. It took 24-25 s a 30x lane (22% of the run), one
   thread busy 90% of it handing records to htslib's writer: gdb samples of it (`perf` is closed to users on the
   server) were the copy of each 64 KB block into a compression job (38%), our copy of each record into a bam1_t
   (26%), the index push waiting for the writer's lock (23%). Now:
   * parsing stores each record as `bam_write1` writes it (`ondisk` in `hts/alnhts.c`; byte-identical to
     `bam_write1` on every record of the small set and on a record of 70,000 CIGAR operations, stored with a CG
     tag), with its end (`bam_endpos`) and mapped flag in the spill frame;
   * the final pass cuts the sorted records into BGZF blocks as htslib's writer does (`bgzf_flush_try`: a record
     that does not fit starts a new block), compresses them in parallel with htslib's `bgzf_compress`, and writes
     them after the header (`bgzf_raw_write`); the `.bai` comes from htslib's `hts_idx_push` with each record's
     virtual offset as a reader's `bgzf_tell` gives it (a record that ends a block ends at the next block's
     start), as `samtools index` builds it;
   * per window of 64 buckets, the sort, the cutting and compression, the write and the index run on their own
     threads, overlapping.

   One 30x lane, L40S, three runs each against main, alternating: final pass 24.6-25.2 -> **4.7-5.0 s** at 64
   threads, lane **107.6 -> 87.6 s**; 23.0-24.5 -> 9.8-10.0 s at 16 threads, lane 193.9 -> 180.0 s; records
   identical, `.bai` identical to `samtools index`. The BAM is 23 KB larger (each window ends its last block). With
   several lanes a run, only the last lane's final pass is not overlapped by alignment. (`ALIGN_RS_TRACE` is not
   free: a traced lane took 99.6 s.)

10. **Less CPU work per lane: tried and dropped (2026-11-19; PLAN.md stage 2).** On a workstation (RTX 5070 Ti, Ryzen 9
    9950X3D, 32 threads) alignment looked CPU-bound (25 of 32 threads busy, GPU half busy). Thread-seconds per 30x
    lane (`ALIGN_RS_PROFILE=2`, inflated by its per-item timers): chaining 910, mate rescue's CPU part 828 (mostly
    its region dedup), CIGAR/NM/MD 587, extension setup 519. Two exact changes cut CPU but not wall clock, and were
    dropped as not worth their code:
    * `pos2rid` by 64 kb block (a table of the sequences per block, against bwa's binary search over ~3,300
      sequences): `chain_seeds` 448 -> 369 thread-s, ~1.5% of the lane's CPU;
    * mate rescue's dedup recomputing only what the region added since the last dedup changes (bwa sorts the
      mate's whole list twice after every rescue direction, ~130 regions, 2.7 M calls per 1 M pairs): 828 -> 353
      thread-s, user CPU -10.5%; ~150 lines with a subtle exactness argument (its check switch found one missed
      case: bwa leaves a list of one region as it is, `n_comp` included).

    Wall clock unchanged within the machine's noise (±8 s on ~2:45 a lane): the lane is paced by its four GPU
    stages, each ~1.4-1.5 s a batch, and their GPU calls are mostly copies through pageable memory, which the
    driver stages on the calling CPU thread: slow when the CPU is busy (in the pipeline ~1 GB/s against 16.7 GB/s
    idle; pinned 51 GB/s on this PCIe 5 card). Pinned SA lookups (a test): their device time 54.5 -> 15.9 s a lane,
    the lane unchanged (the other stages set the pace, and the copies into and out of the pinned buffer cost the
    CPU). Doing that for every GPU stage, filling and reading pinned buffers in place, is the structural fix for
    such machines; not done. Kept from this round: the profile-only size sums (`mbsz`) are evaluated only when
    printed.

11. **Mate SW as one thread per alignment: prototyped, not pursued (2026-11-19; PLAN.md stage 4).** `k_matesw` runs
    bwa's striped SSE ksw_u8 with 16 threads per alignment: per target row a shuffle, a 16-lane maximum and a vote
    (lazy F). The alternative tried: one thread per alignment, the plain recurrence over strips of the query held in
    registers, each row's boundary values (H, F, the row's maximum and its first position, all u8) passed between
    strips through scratch. It reproduces ksw_u8 while nothing saturates (16 x slen <= 250), with the query padded to
    16 x slen positions scoring 0 as `ksw_qinit` pads it (those cells count in the row maxima, which feed score2); the
    forward pass's score, te, qe, score2, te2 were identical on all 65,536 jobs of a dumped chunk (1.04e10 cells).
    Forward pass only, kernel time:

    | | RTX 5070 Ti | L40S |
    |---|---:|---:|
    | `k_matesw` (incl. transfers) | ~26.9 ms | ~16.5 ms |
    | one thread, 16-column strips (96 registers) | 38.5 ms | 23.7 ms |
    | one thread, 32-column strips (188 registers) | 22.8 ms | 12.9 ms |

    The per-row work of a strip (the target base, the boundary load and store) dominates, so wide strips win, and 188
    registers already limit the resident threads. 1.25-1.3x on mate SW is ~60-70 ms of the ~720 ms of kernel time a
    batch on the L40S, perhaps 5-8% of an alignment on one GPU (and nothing where alignment is CPU- or copy-bound),
    for a second kernel with its reverse pass and a fallback for long reads: not pursued.

12. **Mate SW as a wavefront (2026-11-21; PLAN.md stage 7): the default.** `k_matesw_wave` keeps 16 threads per job and
    the striped layout (thread t holds query segment t in registers) but runs thread t one target row behind thread
    t-1 (CUDASW++4.0 / gpuPairHMM), so the row's boundary comes by shuffle: the final H at the segment's end, the F
    leaving it, the row maximum so far, and (only when the row may beat the best so far) the row's largest final H
    with its first position. No lazy-F loop, no per-row vote or reduction. ksw_u8's striped semantics are carried
    exactly: the first sweep's F restarts at each segment and feeds the row maximum and E (hf), while H takes the full
    F (hF), so two F chains; E and the F chains need no clamp at 0 (H >= 0 absorbs a negative value). The last thread
    keeps lane 0's state (b, gmax, te, qe, score2). ALU-bound, so per-cell work decides: computing the position of the
    row's largest H only for segments that could hold it (a bound from the row maximum and the incoming F) saved
    ~1 ms of 19. Register cap 12 blocks of 64 per SM (80 registers), best on both GPUs; `ALIGN_RS_MSW=striped` runs
    `k_matesw` (and `ALIGN_RS_MSW_MINB` / `tune` apply to it only).

    | Kernel time, a lane's chunk (65,536 jobs) | RTX 5070 Ti | L40S |
    |---|---:|---:|
    | `k_matesw` (striped, cap 18) | 31.0 ms | 21.2 ms |
    | `k_matesw_wave` (cap 12) | 18.3 ms | 12.2 ms |

    Checked: every result field identical to `k_matesw` on the dumped chunk and on 600,000 synthetic jobs (queries of
    1-256 bases, every segment count including the saturating 16-segment case, Ns, both strands, every start / minsc /
    endsc flag combination); `tests/cmp.sh` and the pipeline tests identical. Lanes, records identical to bwa-mem2 in
    every run: workstation `mid` 1.60 -> 1.47 s per batch; one 30x lane on one L40S at 64 threads 86.2 / 86.9 ->
    80.7 / 84.1 s (two runs each, alternating). Runs: `$VARCALL_DATA/runs/stage7`, server `runs/stage7`.

13. **Ungapped shortcut for extension (2026-11-22; PLAN.md stage 7): the default.** From bwa-mem3 (PRs #540, #544),
    with the conditions derived for ksw_extend2 (`ext.rs` `ungapped_extend`): when the main diagonal's total loss
    (a - score per cell: a mismatch a + b, an N a + 1) stays below min(o_del + e_del, o_ins + e_ins), every path
    with a gap scores at least 2 below the diagonal in its row and in its column, so the row maxima, the new maxima
    (max_off 0) and gscore at row qlen - 1 are the diagonal's; with tlen >= qlen, every diagonal value > 0 and no
    z-drop on the diagonal, the DP reaches row qlen - 1 as the walk does, and the rows past the query's end are
    gapped. Default scoring: at most one mismatch (or up to three Ns).
    The walk first ran on the CPU while building each round's jobs: exact, but on the workstation the lane slowed
    (1.59 -> 1.74 s a batch: job build 5 -> 7.6 s of CPU against 1.3 s of kernel saved). It now runs on the GPU in
    `k_extkey` (one thread per job, before the size sort): a job it takes gets its result there and the largest key,
    and `k_extend` skips the last *nug jobs; no CPU cost, uploads unchanged. The CPU extender (no GPU) runs the same
    walk. `ALIGN_RS_UNGAPPED=0` turns it off; `=check` runs every job through the DP as well and aborts on a
    difference.
    * Checked: `=check` on the `mid` set (36.6 M of 146.6 M jobs taken, 25%) and on a 30x lane on the L40S (229.8 M
      of 920.8 M), every result identical to the DP; records identical to bwa-mem2 in every run; the CPU extender's
      records identical to bwa-mem2 on tiny / indel / longdel with the shortcut on and off.
    * L40S, one 30x lane at 64 threads, five pairs alternating: extension kernel 29.4-29.5 -> 25.5-26.0 s a lane
      (-12%); the lane 77.4 / 75.6 / 85.4 / 87.4 / 81.5 -> 92.8 / 76.9 / 75.7 / 79.5 / 72.0 s (the last three pairs
      -8 to -10 s; runs vary by ~10 s). Workstation `mid`: kernel 5.35 -> 4.1-4.5 s, 1.53-1.54 -> 1.55-1.57 s a batch
      (not paced by extension).
14. **Certified band for the CIGAR's global alignment (2026-11-22; PLAN.md stage 7): the default.** From bwa-mem3
    (PR #575; `ext.rs` `cert_band`). gen_cigar2 runs ksw_global2 in a band of half-width w around the main diagonal.
    When the query and the reference span have equal lengths, a path with gaps has I >= 1 insertions and as many
    deletions, so it scores at most amax * len - o_ins - o_del - I * (amax + e_ins + e_del), and a path through a
    cell k off the diagonal has I >= k. With su the ungapped score, kc = (amax * len - o_ins - o_del - su) /
    (amax + e_ins + e_del) bounds the cells of every path scoring at least su, so every path to an optimal cell, and
    every predecessor that ties or wins a comparison the backtrack reads, lies within kc: band min(w, kc) gives the
    same score and direction bits on the backtrack, so the same CIGAR (and NM / MD). kc = 0 (the ungapped alignment is
    the only optimum) skips the DP. Unequal lengths keep the usual band. `ALIGN_RS_CIGAR_BAND=0` turns it off;
    `=check` runs the usual band as well and aborts on a difference (score and CIGAR; the score-only variant used
    when merging regions does the same).
    * Workstation `mid` (5 M pairs, 32 threads): of 18.1 M gen_cigar2 calls, 3.9 M run the DP; 2.65 M of those now
      run in a narrower band and 52 k skip it; 1.17 M have unequal lengths. Band cells 22.4 G -> 10.9 G (bwa-mem3
      reports 22.2 -> 10.8 G on its WGS set); CIGAR / NM / MD 140 -> 96 thread-s, user CPU 674 -> 651 s; wall
      unchanged (28-31 s, the GPU paces the lane). `=check`: no difference; records identical in all three modes and
      to bwa-mem2; `tests/cmp.sh`, the germline and somatic pipeline tests and their 40 lane BAMs identical.
    * Not done: bwa-mem3 #580 extends the bound to unequal lengths (a single-gap lower bound; band |d| + x, at least
      |d| + 3); there it took the CIGAR pass's cells down another ~26%, under 1% of bwa-mem3's wall.

