# Building

## Requirements

| | version used | needed for |
|---|---|---|
| Rust | 1.99 | everything |
| a C compiler, make | gcc 13 | the htslib shims (`align-rs/hts`, `variant-rs/csrc`) |
| htslib | 1.25 (1.21 in the images) | align-rs BAM output, variant-rs CRAM input |
| CUDA toolkit (nvcc) | 12.9 | the GPU code (optional) |
| NVIDIA driver | one that supports CUDA 12.9 | running on the GPU |

`environment.yml` installs all of them, plus the tools the tests compare against: bwa-mem2, samtools, bcftools,
Python, Java 17 and Nextflow.

```bash
conda env create -f environment.yml
conda activate varcall-rs
```

Without conda, set `HTSLIB_PREFIX` to an htslib installation (the directory that holds `include/` and `lib/`).

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
