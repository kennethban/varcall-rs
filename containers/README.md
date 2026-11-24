# Containers

The repository provides two container images. Both are built from the repository root, where `.dockerignore`
limits the build context to the sources:

| Image | Contents | Build |
|---|---|---|
| `align-rs:0.1.0` | `align-rs` with its CUDA kernels and htslib, `bwa-mem2` 2.2.1 (for `bwa-mem2 index`), on Ubuntu 22.04 | `docker build -f containers/align-rs/Dockerfile -t align-rs:0.1.0 .` |
| `variant-rs:0.1.0-gatk4.7.0.0` | `variant-rs` with its CUDA PairHMM, added to nf-core/sarek's GATK 4.7.0.0 + samtools image as a replacement for `gatk` ([variant-rs/README.md](variant-rs/README.md)) | `docker build -f containers/variant-rs/Dockerfile -t variant-rs:0.1.0-gatk4.7.0.0 .` |

Both images are built for any x86-64 CPU with AVX2 (`-C target-cpu=x86-64-v3`) and contain GPU code for every
NVIDIA architecture from Volta or Turing to Blackwell, together with PTX (build argument `CUDA_ARCH`). The option
`--gpus all` (NVIDIA Container Toolkit) gives a container access to the GPU. Without this option, or without a
GPU, the CPU paths run and give the same output.

```bash
docker run --rm --gpus all -u $(id -u):$(id -g) -v /data:/data align-rs:0.1.0 \
    align-rs mem -t 16 -K 100000000 -Y --gpu -R '@RG\tID:L1\tSM:S1' -o /data/L1.bam /data/ref.fa /data/R1.fq.gz /data/R2.fq.gz
docker run --rm --gpus all -u $(id -u):$(id -g) -v /data:/data variant-rs:0.1.0-gatk4.7.0.0 \
    gatk HaplotypeCaller -R /data/ref.fa -I /data/S1.bam -O /data/S1.vcf.gz
```

## Testing the images

```bash
containers/test_images.sh            # outputs in target/test_images
```

This script runs `tests/selftest.sh` with both binaries inside the images, on the GPU and CPU paths when the host
has a GPU. It then checks the variant-rs image as sarek uses it: the entrypoint, `gatk --version`, and the transfer
of unsupported commands to GATK. When `VARCALL_DATA` is set ([tests/harness/README.md](../tests/harness/README.md)),
it also runs the command forms of sarek on real data and compares the results with GATK outputs. The host requires
Docker and the dependencies of `tests/selftest.sh` (python3, bwa-mem2, samtools, bgzip and tabix).
