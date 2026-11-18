# variant-rs drop-in image for nf-core/sarek: sarek's own GATK 4.7.0.0 + samtools image (the fallback for everything
# variant-rs does not run) with variant-rs, its CUDA PairHMM and a `gatk` wrapper in place of the image's `gatk`.
#
#   docker build -f containers/variant-rs/Dockerfile -t variant-rs:0.1.0-gatk4.7.0.0 .     # from the repository root
#
# DT_RPATH (not RUNPATH) so that libcudapairhmm.so also finds the bundled libstdc++ in /opt/variant-rs/lib.
# CPU code: x86-64-v3 (AVX2; the AVX-512 PairHMM is chosen at run time where the CPU has it). GPU code: a fat binary
# for sm_70 ... sm_120 plus PTX; without a usable GPU the CPU PairHMM runs (same output).

# sarek's GATK 4.7.0.0 + samtools image (its MarkDuplicates / ApplyBQSR modules convert to CRAM in the same task; the
# other GATK modules use the same GATK without samtools, gatk4-lite:4.7.0.0--79918b8b2632f4b1)
ARG GATK_IMAGE=community.wave.seqera.io/library/gatk4-lite_samtools:9b6806e26dd32afe

# --- build stage (the first steps are the same in containers/align-rs/Dockerfile, so the two builds share them)
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
# --- variant-rs
ARG CUDA_ARCH="-gencode arch=compute_70,code=sm_70 -gencode arch=compute_75,code=sm_75 -gencode arch=compute_80,code=sm_80 -gencode arch=compute_86,code=sm_86 -gencode arch=compute_89,code=sm_89 -gencode arch=compute_90,code=sm_90 -gencode arch=compute_100,code=sm_100 -gencode arch=compute_120,code=sm_120 -gencode arch=compute_120,code=compute_120"
RUN make -C variant-rs/cuda-pairhmm ARCH="${CUDA_ARCH}" CXX_HOST=g++
RUN HTSLIB_PREFIX=/opt/htslib \
        RUSTFLAGS='-C target-cpu=x86-64-v3 -C link-arg=-Wl,--disable-new-dtags -C link-arg=-Wl,-rpath,$ORIGIN/../lib' \
        cargo build --release --locked -p variant-rs \
    && mkdir -p /opt/variant-rs/bin /opt/variant-rs/lib \
    && cp target/release/variant-rs containers/variant-rs/gatk /opt/variant-rs/bin/ \
    && cp variant-rs/cuda-pairhmm/build/libcudapairhmm.so /opt/variant-rs/lib/ \
    && cp -P /opt/htslib/lib/libhts.so* /opt/variant-rs/lib/ \
    && for l in libdeflate.so.0 libbz2.so.1.0 liblzma.so.5 libz.so.1 libstdc++.so.6 libgcc_s.so.1; do \
           cp -L "$(ldconfig -p | awk -v l=$l '$1 == l {print $NF; exit}')" /opt/variant-rs/lib/; done

FROM ${GATK_IMAGE}
COPY --from=build /opt/variant-rs /opt/variant-rs
# The wrapper takes the place of the image's `gatk` (the entrypoint /shell-hook.sh resets PATH, so a PATH entry
# would be skipped); the real launcher moves to gatk.real in the same directory (it finds its jar from there).
RUN real="$(command -v gatk)" && test -n "$real" && mv "$real" "$real.real" \
    && sed -i "s#@GATK_REAL@#${real}.real#" /opt/variant-rs/bin/gatk && cp /opt/variant-rs/bin/gatk "$real" && chmod 755 "$real" \
    && ! ldd /opt/variant-rs/bin/variant-rs | grep 'not found' \
    && gatk --version && VARIANT_RS_DRY_RUN=1 gatk MergeVcfs -I a.vcf.gz -O b.vcf.gz | grep -q '^RUN variant-rs'
ENV PATH=/opt/variant-rs/bin:$PATH
