# HaplotypeCaller in Rust (+ CUDA PairHMM)

A from-scratch port of GATK 4.7 HaplotypeCaller (default arguments; one sample, or several in VCF mode). It
produces a VCF whose data lines are **byte-identical to GATK's** on NA12878 30x, chr20:10–15 Mb
(8,064 / 8,064 records, same order, all 10 columns).

| Stage | Crate | Ported from (GATK 4.7) |
|---|---|---|
| BAM decode, read filters, isActive, activity profile, regions, region read filling, pipeline driver | `variant-rs` (binary, [`..`](..)) | `AssemblyRegionIterator`, `HaplotypeCallerEngine.isActive`, `BandPassActivityProfile` |
| Region read preprocessing | `preproc` | `finalizeRegion`, `ReadClipper`, overlapping-pair quality fix |
| Local assembly | `assembly` | `ReadThreadingAssembler`, read-threading/seq graphs, k-best, Smith-Waterman, CIGAR utils |
| Trimming, genotyping reads, likelihood engine | `likelihood` | `AssemblyRegionTrimmer`, `trimTo`, `PairHMMLikelihoodCalculationEngine` (PCR model, caps, normalisation, poorly-modelled filter) |
| PairHMM on GPU | `cudahmm` + [`../cuda-pairhmm`](../cuda-pairhmm/native/cuda_pairhmm.cu) | GKL AVX-512 PairHMM, bit-exact (FMA placement, separate M/X sums, IEEE division, double FTZ, glibc log10f on the host) |
| Realignment, genotyping, annotations, VCF text | `genotype` | `realignReadsToTheirBestHaplotype`, `HaplotypeCallerGenotypingEngine`, 13 default annotations, htsjdk VCF encoding; ERC genotyping, ERC annotations and physical phasing (`erc.rs`) |
| Reference confidence model, GVCF banding (`-ERC GVCF`) | `refconf` | `ReferenceConfidenceModel`, `GVCFWriter` / `GVCFBlockCombiner` / `HomRefBlock` |
| Positional downsampling (`--max-reads-per-alignment-start 50`) and its shared RNG schedule | `downsample` | `PositionalDownsampler`, `ReservoirDownsampler`, `Utils.getRandomGenerator()` |
| Whole-genome driver: BAI queries, interval/shard semantics, contigs in parallel, bgzip + tabix, generated header | `variant-rs` (binary) | `AssemblyRegionWalker.makeReadShards`, `MultiIntervalLocalReadShard`, htsjdk query/VCF writer |

Each crate has a README with the GATK behaviours it reproduces. Each was validated against GATK dumps at its
stage boundary and fuzzed against GATK's own classes called from gatk.jar.

## End-to-end result (BAM → VCF)

```bash
target/release/variant-rs -I $VARCALL_DATA/data/NA12878.chr20_10-15M.bam -R $VARCALL_DATA/data/chr20.fa -L chr20:10000000-15000000 \
    --assembly-region-out regions.igv -O out.vcf.gz --pairhmm cuda --threads 16
```

Whole genome (one GATK JVM per shard file; see [docs/notes/variant-rs.md](../../docs/notes/variant-rs.md) for modes, results and timings):

```bash
target/release/variant-rs -I $VARCALL_DATA/wgs/bam/NA12878.bam -R $VARCALL_DATA/wgs/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    $(for i in $(seq -w 0 22); do echo --shard $VARCALL_DATA/wgs/shards/shard$i.intervals; done) \
    --shard-output-dir per_shard -O all.vcf.gz [-ERC GVCF] --threads 16
```

| Run | Wall (3 runs) | vs GATK default | VCF data vs GATK default |
|---|---:|---:|---|
| GATK HaplotypeCaller default (GKL AVX-512, 4 OpenMP threads) | 33.0–36.2 s (median 33.8) | 1.0x | — |
| GATK + CUDA PairHMM (per-region JNI calls) | 37.3–38.2 s | 0.9x | identical |
| **Rust, CUDA PairHMM, 1 thread** | 8.3–9.6 s | **3.8x** | **identical** |
| Rust, CPU PairHMM (exact GKL AVX-512 port), 1 thread | 15.0–15.5 s | 2.2x | identical |
| **Rust, CUDA PairHMM, 4 threads** | 3.9–4.2 s | **8.2x** | identical |
| **Rust, CUDA PairHMM, 16 threads** | 3.1 s | **10.8x** | identical |
| Rust, CPU PairHMM, 16 threads | 3.9 s | 8.7x | identical |

Stage breakdown, Rust with 16 threads and CUDA (2.9 s internal total):
- BAM decode, pileup scan and region filling: about 1.6 s, single-threaded, and now the bottleneck.
- Assembly: 0.62 s.
- PairHMM on the GPU: 0.15 s for 1.09 M pairs in 26 launches of 256 regions.
- Realignment: 0.11 s.
- Genotyping and VCF formatting: 0.10 s.
- Preprocessing: 0.11 s.
- Likelihood prep and finish: 0.11 s.

Why the GPU helps here and not inside GATK: the Rust driver batches about 256 active regions (about 40k
read×haplotype pairs) per launch, with no JVM marshaling. GATK calls the PairHMM once per region (about 170 pairs).

## GVCF mode (`-ERC GVCF`)

```bash
target/release/variant-rs -I $VARCALL_DATA/data/NA12878.chr20_10-15M.bam -R $VARCALL_DATA/data/chr20.fa -L chr20:10000000-15000000 \
    -ERC GVCF -O out.g.vcf.gz [--gvcf-raw raw.vcf] --pairhmm cuda --threads 16
```

The banded GVCF is byte-identical to GATK's `$VARCALL_DATA/gvcf/gatk.g.vcf.gz` (798,047 / 798,047 data lines; the generated
header differs only in the command-line line).
The unbanded per-position records are byte-identical to `$VARCALL_DATA/gvcf/gatk_raw_calls.vcf.gz` (5,000,001 / 5,000,001).
Results, timings and the GATK quirks are in [`refconf/README.md`](refconf/README.md).

## Whole genome (NA12878 30x WGS, GRCh38 + decoys + HLA; see [docs/notes/variant-rs.md](../../docs/notes/variant-rs.md))

Per-shard mode vs GATK's per-shard truth (23 shards, one GATK JVM each), raw text compared with zcat:

| Output | `--pairhmm cuda` (bit-exact kernel, 2026-11-12) | `--pairhmm cpu` |
|---|---|---|
| VCF, 5,064,948 lines | **23/23 shards byte-identical** | identical where re-checked (shard00, 18, 19) |
| GVCF, 461,802,035 lines | **23/23 shards byte-identical** | identical where re-checked (shard17, 21) |

(Before the kernel fix, CUDA differed in 20 VCF records and 1,030 GVCF lines, all ±1 in QUAL/PL/GQ; the kernel's
comments in [`../cuda-pairhmm/native/cuda_pairhmm.cu`](../cuda-pairhmm/native/cuda_pairhmm.cu) describe the exact arithmetic it now uses.)

Headers are identical except `##GATKCommandLine` → `##variant-rsCommandLine`.

GIAB HG001 v4.2.1 (`rtg vcfeval`, chr1–22, high-confidence BED, unfiltered calls), whole-genome VCF:

| Caller | TP (baseline) | FN | FP | Recall | Precision | F |
|---|---:|---:|---:|---:|---:|---:|
| GATK 4.7 | 3,698,889 | 23,122 | 32,245 | 0.9938 | 0.9914 | 0.9926 |
| Rust + CUDA | 3,698,889 | 23,122 | 32,245 | 0.9938 | 0.9914 | 0.9926 |

SNPs: recall 99.40%, precision 99.12%. Indels: recall 99.25%, precision 99.24%. The counts are identical, per type.

Whole-genome timings, same machine (Ryzen 9 9950X3D 16C/32T, RTX 5070 Ti):

| Run | VCF | GVCF |
|---|---:|---:|
| GATK 4.7, 23 shard JVMs, 5 at a time (a development script, not in the repository) | 2 h 11 min wall | 2 h 23 min wall |
| GATK 4.7, sum of the 23 per-shard JVM wall times | 9.7 h | 10.7 h |
| **Rust + CUDA, 16 threads, idle machine** | **12 min 51 s** (11,044 s CPU, 9.3 GB RSS) | **17 min 4 s** (14,689 s CPU, 11.2 GB RSS) |

chr1 alone takes most of the Rust wall time (one contig = one sequential pipeline); splitting large contigs is
the next speed-up.

## Limits

Defaults only: diploid, VCF mode or `-ERC GVCF` (several samples: VCF mode, driven by `variant-rs`). DRAGEN, flow, PDHMM, force-calling,
are not ported (contamination and multi-sample are, through `variant-rs`). Input is BAM + BAI or CRAM + CRAI (htslib). Both the CUDA PairHMM and `--pairhmm cpu` are bit-exact
to GKL AVX-512. Validated on one
sample (NA12878, NovaSeq, PCR-free, 30x); exome, PCR+ and non-Illumina data are untested.
