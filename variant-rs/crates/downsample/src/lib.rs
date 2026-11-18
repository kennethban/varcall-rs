//! # hc_downsample
//!
//! Exact port of GATK 4.7 HaplotypeCaller's read downsampling
//! (`--max-reads-per-alignment-start`, default 50): `AssemblyRegionWalker.createDownsampler()` =
//! `PositionalDownsampler` over a `ReservoirDownsampler`, fed by `Utils.getRandomGenerator()`,
//! the same global `java.util.Random(47382911)` that `QualByDepth` jitters QD with.
//!
//! * [`javarandom::JavaRandom`]: bit-exact `java.util.Random` (`nextInt(bound)`, `nextGaussian`).
//!   ONE instance must serve both the downsampler and the genotyping engine's QD jitter, for the
//!   whole run (all contigs, in GATK's traversal order).
//! * [`downsampler::PositionalDownsampler`]: per-shard (contig) downsampler; `push` buffers,
//!   finalizing a group runs the reservoir and makes its draws.
//! * [`schedule`]: the rule fixing where each group's draws fall relative to `apply()` of each
//!   assembly region, plus the helpers to evaluate it ([`schedule::ReadyTracker`],
//!   [`schedule::lookahead_locus`], [`schedule::LibsCoverage`], `PositionalDownsampler::horizon`).
//!
//! See README.md for semantics, the integration recipe and the evidence.

pub mod downsampler;
pub mod javarandom;
pub mod schedule;

pub use downsampler::{reservoir_downsample, DownsampleStats, Horizon, PositionalDownsampler};
pub use javarandom::{JavaRandom, JavaRandomState, GATK_RANDOM_SEED};
pub use schedule::{lookahead_locus, pileup_span, LibsCoverage, ReadyAt, ReadyTracker};
