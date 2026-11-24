# Building

## Requirements

| | version used | needed for |
|---|---|---|
| Rust | 1.99 | everything |
| a C compiler, make | gcc 13 | the htslib shims (`align-rs/hts`, `variant-rs/csrc`) |
| htslib | 1.25 | align-rs BAM output, variant-rs CRAM input |
| CUDA toolkit (nvcc) | 12.9 | the GPU code (optional) |
| NVIDIA driver | one that supports CUDA 12.9 | running on the GPU |

The file `environment.yml` installs all of these, together with the tools that the tests require: bwa-mem2,
samtools, bcftools, Python, Java 17 and Nextflow.

```bash
conda env create -f environment.yml
conda activate varcall-rs
```

Without conda, `HTSLIB_PREFIX` must be set to an htslib installation (the directory that contains `include/` and
`lib/`).

## Hardware

The following table gives the peak memory of each step for a 30x human genome (NA12878, 12 lanes, hg38), as
measured by Nextflow (`peak_rss`), and the elapsed time of the run:

| Step | Workstation: 32 threads, one GPU | Server: 128 threads, 4 GPUs |
|---|---:|---:|
| align-rs, all lanes | 32.5 GB (one task, 12 lanes) | ~36 GB per task (4 tasks of 3 lanes, 32 threads each) |
| MarkDuplicates | 16.5 GB | 60.7 GB |
| BaseRecalibrator | 14.5 GB | 45 GB |
| HaplotypeCaller | 7.3 GB | 20.6 GB |
| hard filters | 0.1 GB | 0.2 GB |
| wall clock | 34.0 min | 9.3 min |

The somatic exome pair (SEQC2, starting from BAM files) requires much less memory: at most 5.1 GB (Mutect2) on the
workstation, and 13.5 GB with 128 threads.

* **Memory.** Alignment has the largest requirement. The hg38 index occupies 17 GB and is memory-mapped, and at 32
  threads the batches in progress occupy a further ~15 GB. Each alignment task requires this amount separately, so
  about 36 GB should be allowed per GPU when one task runs per GPU. The other steps run one at a time and use more
  memory with more threads. 48 GB is sufficient for 32 threads. With 30 GB and 16 GB of swap, the genome was
  processed with the same output in 41.6 min instead of 34 min, because pages of the index were dropped and read
  again. For comparison, bwa-mem2 with `samtools sort` on one lane had a peak of 21 GB.
* **GPU.** An NVIDIA GPU of compute capability 7.5 or newer is required (7.0 for variant-rs). At start-up, align-rs
  checks the free device memory, reserves 2 GB for batches, and loads as many of the three parts of the index as
  fit. For hg38 these are:

  | Part of the index | Device memory | Stages it puts on the GPU |
  |---|---:|---|
  | 2-bit reference | 0.8 GB | seed extension, mate rescue |
  | FM-index checkpoints | 6.4 GB | seeding |
  | sampled suffix array | 4.0 GB | SA lookups |

  A GPU with 16 GB holds all three. With less memory, the stages whose part of the index does not fit run on the
  CPU, which is slower and gives the same output. The PairHMM requires little device memory. Without a GPU, all
  steps run on the CPU.
* **CPU.** An x86-64 processor with AVX2 is required for full speed ([target CPU](#build)). Both machines in the
  table have AMD processors (Ryzen 9 9950X3D and 2 x EPYC 9654).
* **Disk.** The inputs are 60 GB of FASTQ files and the 20 GB hg38 bundle (index and known sites). The work
  directory of the pipeline contained 96 GB at the end of the run (the lane BAMs, the duplicate-marked BAM and the
  VCF files), and the temporary files of MarkDuplicates add to this while it runs.

## Build

```bash
make -C align-rs/cuda                # align-rs/cuda/build/libalngpu.so        (~2.5 min, every architecture)
make -C variant-rs/cuda-pairhmm      # variant-rs/cuda-pairhmm/build/libcudapairhmm.so
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

The build produces `target/release/align-rs` and `target/release/variant-rs`. Their run paths refer to the CUDA
libraries and the htslib against which they were linked, so the binaries can be run in place.

* **CPU target.** `-C target-cpu=native` generates code for the machine on which the build runs. Without it, the
  PairHMM and other intensive loops are several times slower, and the CPU PairHMM is about 20× slower without
  AVX-512 or AVX2 code generation. Builds intended for other machines (the Docker images and binary releases) use
  `-C target-cpu=x86-64-v3`, which requires AVX2. In such builds the AVX-512 PairHMM is still selected at run time
  on CPUs that support it. The outputs are the same for every target.
* **GPU architectures.** Both Makefiles build a binary for several architectures: sm_75 to sm_120 for align-rs and
  sm_70 to sm_120 for variant-rs, together with sm_120 PTX that the driver compiles for newer GPUs. A single
  architecture can be selected for a faster build:
  `make -C align-rs/cuda ARCH="-gencode arch=compute_120,code=sm_120"` (and likewise for `variant-rs/cuda-pairhmm`).
  At start-up, align-rs reports whether kernels are available for the GPU. When the CUDA PairHMM cannot be used,
  variant-rs falls back to the CPU PairHMM (`--pairhmm auto`) or stops (`--pairhmm cuda`).
* **CPU-only build.** The two `make` steps are omitted. Cargo then prints a warning that the CUDA libraries were
  not found and builds only the CPU paths (align-rs without `cfg(gpu)`, variant-rs with the CPU PairHMM). The
  outputs are the same.
* **Host compiler.** nvcc uses `x86_64-conda-linux-gnu-g++` of the conda toolchain when it is present, so that the
  CUDA libraries require the same libstdc++ as the other components. Otherwise it uses `g++`. The compiler can be
  overridden with `CXX_HOST=...`.
* **Release profile.** The workspace has a single profile, with LTO, one codegen unit and line-level debug
  information for profilers and backtraces. align-rs aborts on the first panic in any thread, so that no partially
  written BAM file is finalized. variant-rs unwinds, because it catches panics in a few places, and a panic still
  ends the run with an error before any output is finalized.

The other crates of the workspace (`variant-rs/crates/*`) also contain development binaries, namely validators
against GATK dumps and fuzzers. `cargo build --release` does not build them; an individual binary is built with
`cargo build --release -p hc-<crate>`.

## Tests

```bash
cargo test --release --workspace           # unit tests (Java semantics, RNGs, formatting, CLI translation, ...)
tests/selftest.sh                          # synthetic data: align-rs vs bwa-mem2, variant-rs vs GATK truth (~1 min)
containers/test_images.sh                  # the same through the Docker images
VARCALL_DATA=... GATK_JAR=... tests/harness/regress.sh OUT    # the full harness on public data (~300 GB, ~75 min)
```

[validation.md](validation.md) describes what each test checks. Continuous integration
([.github/workflows/ci.yml](../.github/workflows/ci.yml)) runs the first two on the CPU only. Setting
`VARIANT_RS_NO_AVX512=1` makes the CPU PairHMM use the code path for CPUs without AVX-512.

## Docker

```bash
docker build -f containers/align-rs/Dockerfile   -t align-rs:0.1.0 .
docker build -f containers/variant-rs/Dockerfile -t variant-rs:0.1.0-gatk4.7.0.0 .
```

Both commands are run from the repository root. The images are built for x86-64-v3 and for every GPU architecture
listed above, and they pin the versions of Rust, htslib and the CUDA image. The build arguments are `CUDA_ARCH`,
`RUST_VERSION`, `HTSLIB_VERSION`, `BWA_MEM2_VERSION` (align-rs) and `GATK_IMAGE` (variant-rs; the GATK image that it
extends). Further details are given in [containers/README.md](../containers/README.md).
