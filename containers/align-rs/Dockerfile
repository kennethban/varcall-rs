# align-rs image: `align-rs mem` (bwa-mem2 2.2.1's records; seeding, extension and mate rescue on the GPU) with its
# CUDA kernels and htslib. The index is bwa-mem2's (`bwa-mem2 index`, also in the image).
#
#   docker build -f containers/align-rs/Dockerfile -t align-rs:0.1.0 .        # from the repository root
#   docker run --rm --gpus all -v /data:/data align-rs:0.1.0 align-rs mem -t 16 -Y --gpu ... ref.fa r1.fq.gz r2.fq.gz
#
# CPU code: x86-64-v3 (AVX2). GPU code: a fat binary for sm_75 ... sm_120 plus PTX. Without a usable GPU (or without
# --gpus) align-rs aligns on the CPU, with the same output.

# --- build stage (the first steps are the same in containers/variant-rs/Dockerfile, so the two builds share them)
FROM nvidia/cuda:12.9.1-devel-ubuntu22.04 AS build
ARG HTSLIB_VERSION=1.21
ARG RUST_VERSION=1.99.1
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        curl ca-certificates build-essential pkg-config zlib1g-dev libbz2-dev liblzma-dev libdeflate-dev bzip2 \
    && rm -rf /var/lib/apt/lists/*
RUN curl -fsSL https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 | tar xj -C /tmp \
    && cd /tmp/htslib-${HTSLIB_VERSION} && ./configure --prefix=/opt/htslib --disable-libcurl --disable-gcs --disable-s3 \
    && make -j"$(nproc)" && make install && rm -rf /tmp/htslib-${HTSLIB_VERSION}
RUN curl -fsSL https://sh.rustup.rs | sh -s -- -y --profile minimal --default-toolchain ${RUST_VERSION}
ENV PATH=/root/.cargo/bin:$PATH
WORKDIR /src
COPY . .
# --- align-rs
ARG CUDA_ARCH="-gencode arch=compute_75,code=sm_75 -gencode arch=compute_80,code=sm_80 -gencode arch=compute_86,code=sm_86 -gencode arch=compute_89,code=sm_89 -gencode arch=compute_90,code=sm_90 -gencode arch=compute_100,code=sm_100 -gencode arch=compute_120,code=sm_120 -gencode arch=compute_120,code=compute_120"
ARG BWA_MEM2_VERSION=2.2.1
RUN make -C align-rs/cuda ARCH="${CUDA_ARCH}" CXX_HOST=g++
RUN HTSLIB_PREFIX=/opt/htslib \
        RUSTFLAGS='-C target-cpu=x86-64-v3 -C link-arg=-Wl,--disable-new-dtags -C link-arg=-Wl,-rpath,$ORIGIN/../lib' \
        cargo build --release --locked -p align-rs \
    && mkdir -p /opt/align-rs/bin /opt/align-rs/lib \
    && cp target/release/align-rs /opt/align-rs/bin/ \
    && cp align-rs/cuda/build/libalngpu.so /opt/align-rs/lib/ \
    && cp -P /opt/htslib/lib/libhts.so* /opt/align-rs/lib/
# bwa-mem2 for `bwa-mem2 index`, built from source with gcc as bioconda does (the release's Intel-compiled AVX-512
# binary refuses non-Intel CPUs), with its own __rdtsc() replaced by the compiler's (GCC >= 11 defines it; profiling
# timers only); `make multi`: one build per instruction set, the launcher picks the CPU's
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends git \
    && rm -rf /var/lib/apt/lists/* \
    && git clone --recursive --branch v${BWA_MEM2_VERSION} https://github.com/bwa-mem2/bwa-mem2 /tmp/bwa-mem2 \
    && perl -0pi -e 's/#if defined\(__GNUC__\) && !defined\(__clang__\)\n#if defined\(__i386__\).*?#endif\n#endif\n/#include <x86intrin.h>\n/s' /tmp/bwa-mem2/src/utils.h \
    && grep -q '^#include <x86intrin.h>' /tmp/bwa-mem2/src/utils.h \
    && make -C /tmp/bwa-mem2 multi -j"$(nproc)" \
    && cp /tmp/bwa-mem2/bwa-mem2 /tmp/bwa-mem2/bwa-mem2.* /opt/align-rs/bin/

FROM ubuntu:22.04
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        zlib1g libbz2-1.0 liblzma5 libdeflate0 \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /opt/align-rs /opt/align-rs
ENV PATH=/opt/align-rs/bin:$PATH
RUN ! ldd /opt/align-rs/bin/align-rs | grep 'not found' && align-rs 2>&1 | grep -q mem
