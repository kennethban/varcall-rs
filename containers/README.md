# Containers

Two images, both built from the repository root (`.dockerignore` limits the context to the sources):

| Image | Contents | Build |
|---|---|---|
| `align-rs:0.1.0` | `align-rs` with its CUDA kernels and htslib, `bwa-mem2` 2.2.1 (for `bwa-mem2 index`), on Ubuntu 22.04 | `docker build -f containers/align-rs/Dockerfile -t align-rs:0.1.0 .` |
| `variant-rs:0.1.0-gatk4.7.0.0` | `variant-rs` with its CUDA PairHMM, added to nf-core/sarek's GATK 4.7.0.0 + samtools image as a drop-in `gatk` ([variant-rs/README.md](variant-rs/README.md)) | `docker build -f containers/variant-rs/Dockerfile -t variant-rs:0.1.0-gatk4.7.0.0 .` |

Both are built for any x86-64 CPU with AVX2 (`-C target-cpu=x86-64-v3`) and carry GPU code for every NVIDIA
architecture from Volta / Turing to Blackwell plus PTX (`CUDA_ARCH` build argument). Run them with `--gpus all`
(NVIDIA Container Toolkit) to use the GPU; without it, or without a GPU, the CPU paths run and give the same output.

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

runs `tests/selftest.sh` with both binaries inside the images (the GPU and CPU paths when the host has a GPU), then
checks the variant-rs image the way sarek uses it (entrypoint, `gatk --version`, hand-over to GATK). With
`VARCALL_DATA` set (see [tests/harness/README.md](../tests/harness/README.md)) it also runs sarek's command forms on
real data against GATK truth. The host needs what `tests/selftest.sh` needs (python3, bwa-mem2, samtools, bgzip,
tabix) and Docker.
