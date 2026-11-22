# Building

## Requirements

| | version used | needed for |
|---|---|---|
| Rust | 1.99 | everything |
| a C compiler, make | gcc 13 | the htslib shims (`align-rs/hts`, `variant-rs/csrc`) |
| htslib | 1.25 | align-rs BAM output, variant-rs CRAM input |
| CUDA toolkit (nvcc) | 12.9 | the GPU code (optional) |
| NVIDIA driver | one that supports CUDA 12.9 | running on the GPU |

`environment.yml` installs all of them, plus the tools the tests compare against: bwa-mem2, samtools, bcftools,
Python, Java 17 and Nextflow.

```bash
conda env create -f environment.yml
conda activate varcall-rs
```

Without conda, set `HTSLIB_PREFIX` to an htslib installation (the directory that holds `include/` and `lib/`).

## Hardware

What a 30x human genome (NA12878, 12 lanes, hg38) needed, measured by Nextflow (`peak_rss`):

| Step | Workstation: 32 threads, one GPU | Server: 128 threads, 4 GPUs |
|---|---:|---:|
| align-rs, all lanes | 32.5 GB (one task, 12 lanes) | ~36 GB per task (4 tasks of 3 lanes, 32 threads each) |
| MarkDuplicates | 16.5 GB | 60.7 GB |
| BaseRecalibrator | 14.5 GB | 45 GB |
| HaplotypeCaller | 7.3 GB | 20.6 GB |
| hard filters | 0.1 GB | 0.2 GB |
| wall clock | 34.3 min | 9.3 min |

The somatic exome pair (SEQC2, from BAMs) needs much less: at most 5.1 GB (Mutect2) on the workstation, 13.5 GB with
128 threads.

* **Memory:** align-rs needs the most: the hg38 index is 17 GB and is memory-mapped, plus ~15 GB of batches in
  flight at 32 threads. Each alignment task needs its own, so with one task per GPU, plan ~36 GB per GPU. The other
  steps run one at a time and use more memory with more threads. 48 GB is enough for 32 threads. With 30 GB (and
  16 GB of swap) the genome still ran, in 41.6 min instead of 34.3, with the same output: the index's pages were
  dropped and re-read. For comparison, bwa-mem2 + samtools sort on one lane peaked at 21 GB.
* **GPU:** an NVIDIA GPU of compute capability 7.5 or newer (7.0 for variant-rs). align-rs loads as much of the
  index as fits, checked at start, and keeps 2 GB for batches: the 2-bit reference (0.8 GB) for extension and
  mate rescue, + the FM-index checkpoints (6.4 GB) for seeding, + the sampled suffix array (4.0 GB) for SA
  lookups. A 16 GB card holds all three for hg38. With less, the stages that do not fit run on the CPU, slower, with
  the same output. The PairHMM needs little device memory. Without a GPU, everything runs on the CPU.
* **CPU:** x86-64 with AVX2 for full speed ([target CPU](#build)). Both machines above are AMD (Ryzen 9 9950X3D,
  2 x EPYC 9654).
* **Disk:** the inputs are 60 GB of FASTQ and the 20 GB hg38 bundle (index, known sites). The pipeline's work
  directory held 96 GB at the end (the lane BAMs, the duplicate-marked BAM, the VCFs); MarkDuplicates' temporary
  files add to it while it runs.

## Build

```bash
make -C align-rs/cuda                # align-rs/cuda/build/libalngpu.so        (~2.5 min, every architecture)
make -C variant-rs/cuda-pairhmm      # variant-rs/cuda-pairhmm/build/libcudapairhmm.so
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

This gives `target/release/align-rs` and `target/release/variant-rs`. Their run paths point at the CUDA libraries
and htslib they were linked against, so they run in place.

* **CPU target.** `-C target-cpu=native` is for the machine that builds. Leave it out and the PairHMM and other
  hot loops run several times slower: the CPU PairHMM is about 20× slower without AVX-512 / AVX2 code generation.
  Builds for other machines (the Docker images, binary releases) use `-C target-cpu=x86-64-v3`, which needs AVX2.
  The AVX-512 PairHMM is then still chosen at run time on CPUs that have it. Outputs are the same for every target.
* **GPU architectures.** Both Makefiles build a fat binary: sm_75 … sm_120 for align-rs, sm_70 … sm_120 for
  variant-rs, plus sm_120 PTX that the driver compiles for newer GPUs. For a faster build, choose one
  architecture:
  `make -C align-rs/cuda ARCH="-gencode arch=compute_120,code=sm_120"` (and the same for `variant-rs/cuda-pairhmm`).
  At start, align-rs reports whether the GPU has kernels. When the CUDA PairHMM is not usable, variant-rs falls back
  to the CPU PairHMM (`--pairhmm auto`) or stops (`--pairhmm cuda`).
* **CPU only.** Skip the two `make` steps. Cargo prints a warning that the CUDA libraries were not found and builds
  the CPU paths only (align-rs without `cfg(gpu)`, variant-rs with the CPU PairHMM). The outputs are the same.
* **Host compiler.** nvcc uses the conda toolchain's `x86_64-conda-linux-gnu-g++` when there is one, so that the
  CUDA libraries need the same libstdc++ as the rest. Otherwise it uses `g++`. Override it with `CXX_HOST=...`.
* **Release profile.** One profile for the workspace: LTO, one codegen unit, line-level debug info (for profilers
  and backtraces). align-rs aborts on the first panic in any thread, so that no half-written BAM is finalized.
  variant-rs unwinds, because it catches panics in a few places, and a panic still ends the run with an error
  before any output is finalized.

The other workspace crates (`variant-rs/crates/*`) also have development binaries: validators against GATK dumps,
and fuzzers. `cargo build --release` does not build them. Build one with `cargo build --release -p hc-<crate>`.

## Tests

```bash
cargo test --release --workspace           # unit tests (Java semantics, RNGs, formatting, CLI translation, ...)
tests/selftest.sh                          # synthetic data: align-rs vs bwa-mem2, variant-rs vs GATK truth (~1 min)
containers/test_images.sh                  # the same through the Docker images
VARCALL_DATA=... GATK_JAR=... tests/harness/regress.sh OUT    # the full harness on public data (~300 GB, ~75 min)
```

See [validation.md](validation.md) for what each one checks. CI ([.github/workflows/ci.yml](../.github/workflows/ci.yml))
runs the first two, CPU only. `VARIANT_RS_NO_AVX512=1` makes the CPU PairHMM take the path for CPUs without AVX-512.

## Docker

```bash
docker build -f containers/align-rs/Dockerfile   -t align-rs:0.1.0 .
docker build -f containers/variant-rs/Dockerfile -t variant-rs:0.1.0-gatk4.7.0.0 .
```

Run both from the repository root. They build for x86-64-v3 and every GPU architecture listed above, and pin
Rust, htslib and the CUDA image. Build arguments: `CUDA_ARCH`, `RUST_VERSION`, `HTSLIB_VERSION`,
`BWA_MEM2_VERSION` (align-rs), and `GATK_IMAGE` (variant-rs: the GATK image it extends). See
[containers/README.md](../containers/README.md).
