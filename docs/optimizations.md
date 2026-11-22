# Optimizations

Every optimization here keeps the output identical to the reference tool. The pattern is always the same:

1. find the invariant the reference's result actually depends on (usually an order, or a rounding);
2. compute the same thing faster while keeping that invariant;
3. check it. Either keep the plain version beside the fast one behind a check switch, or compare full outputs on
   data that exercises the path.

Measured first, every time: `ALIGN_RS_PROFILE=1` / `VARIANT_RS_TRACE_STAGES=1` give per-stage times, thread CPU
time and counters; `ALIGN_RS_TRACE=file` writes a per-batch timeline of align-rs's stages and GPU calls. The
numbers below are from one 16-core machine with an RTX 5070 Ti unless they name an L40S (a 2 × EPYC 9654 server;
details in [notes/](notes/)).

## align-rs

### Exact fast paths, with a fallback for ties

bwa-mem2 sorts regions with klib's `ks_introsort`, which is not stable: regions with equal keys come out in an
order only that algorithm produces. When no two regions with equal keys differ, any sort gives the same result.
That holds for 99.6% of sorts on real reads, so a fast sort runs first and the ported introsort runs only for the
rest (`align-rs/src/ext.rs`):

```rust
pub fn a_idx_sort<K: Copy>(a: &mut [Reg], key: impl Fn(&Reg) -> K, lt: impl Fn(&K, &K) -> bool) {
    let mut ki: Vec<(K, u32)> = a.iter().enumerate().map(|(i, r)| (key(r), i as u32)).collect();
    let mut fast = ki.clone();
    fast.sort_unstable_by(|x, y| ord(&x.0, &y.0));
    // regions of equal keys must all be the same region
    let ok = fast.windows(2).all(|w| lt(&w[0].0, &w[1].0) || a[w[0].1 as usize] == a[w[1].1 as usize]);
    if ok { ki = fast; } else { crate::ksort::introsort(&mut ki, |x, y| lt(&x.0, &y.0)); }
    // ... permute `a` by ki
}
```

### An incremental dedup, exact by an invariant

bwa-mem2's mate rescue runs `mem_sort_dedup_patch` on the mate's regions after each of its four directions: two
sorts and an overlap loop over every region, for at most one region added. Reads with hundreds of regions made it
18% of the workstation's CPU. align-rs keeps a flag per mate: its regions are *clean* when their ends differ, the
overlap loop (by end) would drop none of them, and their final keys differ. A region added to a clean set can only
end the loop's windows earlier, so only the tests with the new region can drop anything: its own turn (the regions
ending before it, back to the first that ends its window) and the turns of the later regions whose windows reach
it, in end order. `dedup_add` runs just those and puts the region in its place in the score order; nothing of the
old set dropped, the result is clean again. A direction that rescued nothing leaves a clean set unchanged. Ties,
equal ends or keys, and the first call per mate run the full function, which reports whether its result is clean
(the loop rerun as a check). `ALIGN_RS_DEDUP_CHECK=1` runs the full function on a copy every time: equal in 84 M
calls on a 30x lane, of which 74% were incremental and 17% skipped. With the region sorts on one packed integer
key, permuted in place: the lane 149.5 -> 138.8 s, its CPU time -21%.

### Replacing a linear scan, with a check mode

After extension, bwa-mem2 replays bwa's skip test: for every seed, it scans the kept regions in order and stops at
the first that contains the seed. Only *whether* one contains it matters, so align-rs keeps the regions sorted by
start and tests only those that could contain the seed. That is a binary-search window `[seed end − longest
region, seed start]`. `ALIGN_RS_REPLAY_CHECK=1` runs bwa's scan beside the lookup and stops on any disagreement.
This check found a seed longer than every kept region (an empty window) on the second lane of real data.

```rust
let lo = kept.partition_point(|x| x.0 < s.rbeg + s.len as i64 - max_len);
let hi = kept.partition_point(|x| x.0 <= s.rbeg);
let contained = lo < hi && kept[lo..hi].iter().any(|&(_, pi)| test(&av[pi as usize]));
if replay_check() {
    // bwa's scan: the first `lim` regions not purged, in order
    /* ... */
    assert_eq!(contained, scan, "replay lookup differs from bwa's scan");
}
```

Effect: 37.9 → 15.1 thread-seconds on 1/6 of a lane.

### Seed extension on the GPU: one memory word per DP cell

Banded Smith-Waterman (`ksw_extend2`) runs one job per thread. Each thread's DP row lives in local memory, and
Nsight Compute showed the kernel bound by L2 traffic. H, E and the query base are packed into one 32-bit word, so
a cell costs one load and one store (`align-rs/cuda/aln_ext.cu`, `k_extend`):

```cuda
// H (bits 0-13), E (14-27) and the query base (28-30): one local load and one store per cell.
// H and E stay below h0 + qlen * max(mat); jobs where that could reach 2^14 go to the CPU.
#define EH(h, e) ((uint32_t)(h) | ((uint32_t)(e) << 14))
...
uint32_t x = eh[j];
int M = EH_H(x), e = EH_E(x), h, tt;
M = M ? M + mrow[EH_Q(x)] : 0;
...
eh[j] = EH(h1, e) | (x & EH_QM);   // H(i,j-1) for the next row, E(i+1,j)
```

Along with that:

* jobs sorted by size on the device (CUB radix sort);
* 24-byte jobs, with per-round constants passed as kernel arguments;
* a grid-stride loop over a fixed number of blocks per SM, so that the rows stay in L1 (`ALIGN_RS_EXT_BPSM`).

Kernel 6.9 → 4.7 s on 1/6 of a lane.

### Seeding: a group of lanes per read

One thread per read left 3 of 32 lanes busy per warp: a few repetitive reads per warp run ~10× longer than the
others. `k_smem_g<G>` gives each read G lanes. The forward walks run in lockstep, and the backward phase's
candidates are extended G at a time. Then every lane replays bwa's selection in order on shuffled values, so the
seeds are exactly `one_pos`'s, including an `int32` truncation. Groups are persistent and take reads from a
counter, so a slow read no longer holds 31 idle lanes:

```cuda
for (;;) {
    int r = 0;
    if (lane == 0) r = atomicAdd(next_read, 1);
    r = __shfl_sync(mask, r, 0, G);
    if (r >= n_reads) return;
    ...
    while (x < len && !out.over) x = one_pos_g<G>(cp, q, len, x, 1, min_seed_len, out, prev, lane, mask);
```

Seeding kernel 1.22 → 0.68 s on 1M pairs (G = 4).

### Mate rescue kernel: occupancy

Mate SW emulates bwa's striped SSE `ksw_u8` with a half-warp per job, and it waits on shuffles and votes, so
resident warps decide its speed. These changes took it from 164 to 306 GCUPS, with every result field identical
(`cuda/bench_msw.cu` checks them on a dumped chunk):

* the query length as a template parameter (rows in registers, not local memory);
* no u8 saturation where H provably cannot reach 255;
* registers capped at 56 for 64-thread blocks.

### Mate rescue as a wavefront

The striped kernel still waited on its lazy-F loop and on a vote and reduction per target row. `k_matesw_wave`
keeps 16 threads per job but gives thread t the query segment t (the columns SSE lane t holds) in registers, one
target row behind thread t-1 (the register-tiled wavefront of CUDASW++4.0). Thread t-1 hands over, by shuffle, the
row's boundary H, the F leaving its segment and the running row maximum, so no row needs a vote or a second pass.
bwa's striped semantics are kept exactly: the row maximum and E use the values before the lazy-F pass (F within
the segment), H the full F, so the kernel carries two F chains. The argmax over a row is computed only when the
row could beat the maximum seen a few rows back (a bound the shuffles carry for free).

Every result field identical to the striped kernel (a lane's chunk, 600,000 synthetic jobs); 21 -> 12 ms per
65,536 jobs on an L40S, 31 -> 18 on the RTX 5070 Ti. `ALIGN_RS_MSW=striped` runs the old kernel.

### Exact shortcuts from bwa-mem3

bwa-mem3 found two places where bwa's DP provably gives what a cheaper computation gives:

* **Ungapped extension.** When the seed diagonal's total loss (a - score per cell) stays below the cheaper gap
  open + extend, every gapped path scores at least 2 below the diagonal in its row and column, so ksw_extend2's
  row maxima, new maxima and gscore are the diagonal's. Conditions as derived for this code: tlen >= qlen, every
  diagonal value > 0, no z-drop on the diagonal. The walk runs on the GPU in the key kernel before the size sort;
  a quarter of the jobs take it and skip the DP (extension kernel -12% on the L40S). `ALIGN_RS_UNGAPPED=check`
  runs the DP as well and aborts on a difference: none in 229.8 M shortcut jobs of a 30x lane.
* **Certified CIGAR band.** For equal-length regions, any path more than kc cells off the diagonal scores below the
  ungapped score, kc = (amax·len − o_ins − o_del − su) / (amax + e_ins + e_del); band min(w, kc) gives ksw_global2
  the same score and the same direction bits along the backtrack. Half the CIGAR stage's band cells (22.4 G ->
  10.9 G on 5 M pairs), its CPU 140 -> 96 thread-s; kc = 0 skips the DP. `ALIGN_RS_CIGAR_BAND=check` compares.

```rust
pub fn cert_band(o: &ExtOpt, q: &[u8], r: &[u8]) -> i64 {
    // ... equal lengths, non-negative gap opens, a positive slope
    let su: i64 = q.iter().zip(r).map(|(&a, &b)| o.mat[b as usize * 5 + a as usize] as i64).sum();
    let num = amax * q.len() as i64 - o.o_ins as i64 - o.o_del as i64 - su;
    if num < 0 { 0 } else { num / slope }
}
```

### Pipeline and memory

* **Stages:** five stages connected by unbuffered channels; pairing / SAM and rescue in separate rayon pools.
  One lane's batch went from 3.35 s to 1.80 s, and the stages now take within 10% of each other.
* **Gunzip:** one thread per FASTQ file. Read stage 1.62 → 1.05 s per batch.
* **Allocator:** mimalloc, set to purge freed memory with `MADV_FREE` instead of decommitting it, so the next batch
  reuses the pages. Minor page faults fell from 10.6 M to 2.0 M.
* **BAM output** (bam.rs): instead of SAM → `samtools sort`, records go per batch into 1024 genome buckets and are
  spilled compressed. At the end, each bucket is sorted stably by (tid, pos, strand), which is samtools sort's
  order, and written with htslib.
* **Final pass in parallel** (bam.rs, hts/alnhts.c): one thread handing records to htslib's writer took 25 s of a
  30x lane. Now records are stored as `bam_write1` writes them, cut into BGZF blocks by htslib's own rule
  (`bgzf_flush_try`), compressed in parallel with `bgzf_compress`, and indexed with `hts_idx_push` at the virtual
  offsets a reader sees; sort, compression, write and index overlap window by window. The records are unchanged,
  and the `.bai` is byte for byte what `samtools index` builds for the BAM. L40S, 64 threads: final pass 25 → 5 s, lane 108 → 88 s.
* **Extension on many cores:** extension's CPU steps get their own rayon pool, like SAM and rescue, so they no
  longer queue behind the other stages' loops in the global pool (`ALIGN_RS_EXT_POOL=0` to compare); the kernel
  counts the jobs it leaves to the CPU, so the scan for them runs only when there are some. L40S, 64 threads:
  lane 113.6 → 108.5 s.
* **Freeing off the stage threads, bounded:** each batch's large structures are freed by dropper threads, so a
  stage does not pay for it. One dropper thread fell behind once batches got faster: its queue held ~1 GB a
  batch, a workstation's 12 lanes filled 48 GB in ~6 minutes, and the server's alignment tasks grew to 54-169 GB.
  Four threads now share a queue of at most 16 structures (~4 batches), so a stage waits instead: 37 GB a task,
  the six tasks of the server's germline run finish together, 460 → 409 s for the whole pipeline.

**Measured and dropped:**

* running only the extensions bwa's replay keeps (73% of jobs);
* de-duplicating SA lookups;
* pinned staging buffers (slower under WSL2);
* stage slack (more memory, no gain);
* less CPU work per lane (`pos2rid` by block, incremental rescue dedup): CPU down, wall clock unchanged, since the
  lane is paced by its GPU stages;
* mate SW as one thread per alignment: 1.25-1.3× on that kernel, ~5% of a lane, for a second kernel;
* idle lanes in the extension kernel taking the next job (AGAThA): its lanes are 53% busy, but with only k of 32
  lanes given work the kernel takes ×1.03 / ×1.10 / ×1.31 the time at k = 24 / 16 / 8. It is paced by its
  working lanes' local-memory traffic, so filling the idle ones would save ~9% of it at best;
* bwa-mem3's certified adaptive extension band (probe at band 20, keep the results a score bound proves final,
  rerun the rest): a quarter of the DP jobs certify, since the ungapped shortcut has taken the easy ones, and probe
  plus rerun take ×1.08 the kernel time;
* pinned extension buffers filled and read in place: the copies 3-4× faster, the wall clock unchanged, since
  extension's GPU rounds already overlap its CPU work, which is the longer;
* merging per-read vectors into fewer, larger ones (extension's left and right templates as pairs, `chain_flt`'s
  per-chain arrays as tuples): a repetitive read's allocation then lands in a larger mimalloc class, whose pages
  are purged and faulted in again (page faults ×3-5, a lane no faster). Allocation count is not the only cost:
  each one's size class matters too.

Details of each: [notes/align-rs.md](notes/align-rs.md).

## variant-rs

### Integers where GATK sums floating point

BaseRecalibrator without BAQ counts errors as 0 or 1, so each table cell is a pair of integers. Integer sums don't
depend on order, so per-chunk tables can be summed in any order and still give GATK's sequential result. Only
GATK's floating-point finalization (`RecalDatum.combine`, rounding, the Bayesian quality) runs once, in GATK's
order (`variant-rs/src/bqsr/recal.rs`):

```rust
/// Without BAQ (GATK's default) every error is 0 or 1, so a table cell is just (observations,
/// errors) integer counts: per-thread tables can be summed in any order and give GATK's values.
pub struct Counts {
    pub qual: Vec<[u64; 2]>,      // [rg][qual]
    pub context: Vec<[u64; 2]>,   // [rg][qual][context]
    pub cycle: Vec<[u64; 2]>,     // [rg][qual][cycle]
    ...
}
```

Whole genome: 95 s vs GATK's 99 min.

### No pileup objects

GATK builds a pileup object per locus. variant-rs walks each read's CIGAR once and adds its bases straight into a
ring buffer of per-locus accumulators (`frontend.rs`, `Ring`). Reads are added in file order, so the
floating-point sums are added in the same order as GATK's and come out bit-identical. One thread adds reads to
loci ahead of the frontier while another evaluates the loci behind it. A read spanning more than the ring grows
it first.

### Parallel contigs, sequential random numbers

GATK draws downsampling and QD-jitter random numbers from one global `java.util.Random`, in traversal order. To
run contigs in parallel and still draw in that order, the generator is a baton (`variant-rs/src/rng.rs`):

```rust
/// Blocks until contig `k`'s turn (every earlier contig of the domain is done with the generator).
pub fn take(&self, k: usize) -> T { /* wait on a Condvar until cursor == k */ }

/// Contig `k` is done with the generator.
pub fn finish(&self, k: usize, rng: T) {
    let mut s = self.m.lock().unwrap();
    assert!(s.cursor == k && s.rng.is_none(), "RNG baton out of order");
    s.rng = Some(rng);
    s.cursor = k + 1;
    Self::skip(&mut s);   // contigs that draw nothing are passed over
    ...
}
```

Two more steps keep the draws in order without waiting:

* Genotyping writes a placeholder for each QD jitter and fills it in when the baton arrives.
* The downsampler, which needs the generator's true state, waits only at its first drawing group.

So contigs without drawing groups (most of them) never wait.

### A bit-exact GPU PairHMM

GATK's PairHMM is Intel GKL's AVX-512 kernel. Its rounding decides PLs and QUAL, so an approximately-equal GPU
kernel gives different VCF lines: an earlier CUDA version differed from GATK in 15 records of the first 3.8 M. `cuda_pairhmm.cu` is
compiled with `-fmad=false -ftz=true`, and every operation that GKL performs is spelled out with an explicit
rounding intrinsic. Fused multiply-adds appear only where GKL fuses, and doubles are flushed to zero by hand,
because GKL's MXCSR flushes doubles too:

```cuda
// Flush-to-zero. GKL sets MXCSR.FTZ, which flushes denormal results of float AND double operations.
// nvcc -ftz=true covers float only, so every double result is flushed explicitly.
// |x| < DBL_MIN exactly when the exponent field is 0: tested on the integer bits, off the slow FP64 units.
__device__ __forceinline__ double FZ(double x) {
  const long long u = __double_as_longlong(x);
  return (u & 0x7ff0000000000000LL) == 0 ? __longlong_as_double(u & (long long)0x8000000000000000ULL) : x;
}
// Explicit fused multiply-add (the file is compiled with -fmad=false, so nothing else is contracted).
__device__ __forceinline__ float  FMA(float a, float b, float c) { return __fmaf_rn(a, b, c); }
__device__ __forceinline__ double MUL(double a, double b) { return FZ(__dmul_rn(a, b)); }
```

Pairs are batched across regions into one launch. At first each warp computed one (read, haplotype) pair (whole
genome: 23/23 shards identical, GPU busy ~220 s of a ~12 min run). The float pass is now tiled, as in gpuPairHMM
(Schmidt et al., arXiv 2411.11547): a few threads per pair, each holding K read rows in registers and sweeping the
haplotype one step behind its neighbour, with pairs binned by read length. Same operations per cell, so the same
bits, at 2.4-2.8x the speed (L40S 282 -> 663 GCUPS).

The contig workers share a pool of engines (a CUDA stream and buffers each, 8 per GPU by default) instead of one
engine behind a lock, so one worker's copies overlap another's kernels; the engines spread over every visible GPU.
30x genome on an L40S at 64 threads: 190 → 122 s; on 2 GPUs 115 → 98 s (the longest calling interval, 141 Mb of
chr2, sets the floor).

### Java semantics as small, shared functions

A single wrong tie-break or rounding changes a VCF line, so Java's behaviour is reproduced in shared helpers, each
with a unit test. For example, `f64::total_cmp` sorts the NaN that `0.0/0.0` produces on x86-64 *first*, while
Java's `Double.compare` puts every NaN *last* (`crates/assembly/src/jmath.rs`):

```rust
/// `Double.compare`: every NaN is equal to itself and above +Infinity, -0.0 < 0.0.
pub fn java_cmp(a: f64, b: f64) -> std::cmp::Ordering {
    match (a.is_nan(), b.is_nan()) {
        (true, true) => std::cmp::Ordering::Equal,
        (true, false) => std::cmp::Ordering::Greater,
        (false, true) => std::cmp::Ordering::Less,
        _ => a.total_cmp(&b),
    }
}
```

Others of the same kind:

* `Math.log10` / `exp` / `pow` as the JDK computes them;
* `String.format` rounding;
* htsjdk's BAM tag re-encoding: tags sorted by binary key, integers narrowed to the smallest type, H tags written
  back as `B:c` (`bamwrite.rs`, `normalize_tags`);
* Picard's `HashSet` iteration order.

### Streaming instead of re-reading

* **MarkDuplicates:** pass 1 runs over 2 Mb chunks in parallel. Read ends go to on-disk buckets keyed by
  (contig, position), and buckets are grouped in parallel with Picard's comparator. Pass 2 rewrites records in
  parallel and appends them in order. Whole genome: 6 min vs Picard's 42 min. The ends on disk are packed
  (35 bytes a fragment, 50 a pair) and read straight into their arrays.
* **ApplyBQSR:** rewrites each record's quality bytes in place and re-encodes the tags as htsjdk does.
* **BGZF:** compressed and decompressed in parallel with libdeflate everywhere. The `.bai` / `.tbi` are built while
  writing.

## Proposing an optimization

Show that the output is unchanged on data that exercises the path (see [CONTRIBUTING.md](../CONTRIBUTING.md)).
Ideally, add a check mode like `ALIGN_RS_REPLAY_CHECK` that compares the fast path with the plain one, and include
before / after timings from the profile switches above.
