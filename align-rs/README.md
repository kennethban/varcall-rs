# align-rs

`bwa-mem2 mem` 2.2.1 in Rust and CUDA: the same records as bwa-mem2, with seeding, SA lookups, seed extension and
mate rescue on the GPU, and the coordinate-sorted BAM + `.bai` written directly (records as `samtools sort`,
index bytes as `samtools index`). One NA12878 30x lane: 124 s vs 587 s for bwa-mem2 + samtools sort + index.

```bash
make -C cuda                      # GPU kernels (optional; without them align-rs runs on the CPU, same output)
cargo build --release             # from the repository root
align-rs mem -t 32 -K 100000000 -Y --gpu -R '@RG\tID:L1\tSM:S1' -o S1.L1.bam ref.fa R1.fq.gz R2.fq.gz
align-rs mem -t 32 -K 100000000 -Y --gpu --lanes lanes.tsv ref.fa        # several lanes, index loaded once
align-rs tune -t 32 -K 100000000 ref.fa R1.fq.gz R2.fq.gz                # measure GPU settings for this GPU
```

| | |
|---|---|
| `src/main.rs` | options, the stage pipeline, lanes |
| `src/fastq.rs` | FASTQ / gzip reading, `kseq_read` byte for byte |
| `src/index.rs` | the bwa-mem2 index (`.bwt.2bit.64`, `.0123`, `.pac`, `.ann`, `.amb`, `.alt`), memory-mapped |
| `src/smem.rs`, `src/chain.rs`, `src/ksort.rs` | seeding, chaining (klib B-tree and introsort ported) |
| `src/ext.rs`, `src/ksw.rs` | seed extension, bwa's skip replay, dedup; banded and global SW |
| `src/sam.rs` | insert sizes, mate rescue, pairing, MAPQ, CIGAR / MD / XA / SA, SAM text |
| `src/bam.rs`, `hts/alnhts.c` | sorted BAM + index: htslib's record encoding, BGZF blocks and index, the final pass in parallel |
| `src/gpu.rs`, `cuda/aln_ext.cu` | device memory, kernels (`k_smem_g`, `k_sa`, `k_extend`, `k_matesw`), CPU fallback |
| `src/tune.rs`, `src/prof.rs` | `align-rs tune`, profiling counters, the `ALIGN_RS_TRACE` timeline |
| `tests/` | edge cases vs bwa-mem2 (`badfq.sh`, `alt.sh`, `gpufail.sh`, `big_bucket.sh`, `cmp.sh`) |
| `tools/` | timing scripts, lane digests |

See [docs/reference.md](../docs/reference.md) (options, environment variables), [docs/architecture.md](../docs/architecture.md),
[docs/optimizations.md](../docs/optimizations.md) and the [development notes](../docs/notes/align-rs.md).
