//! Port of the last stage of GATK 4.7 HaplotypeCaller (`HaplotypeCallerEngine.callRegion` after
//! `computeReadLikelihoods`): realignment of reads to their best haplotype, genotyping,
//! annotation and VCF record output, with default arguments. See README.md.

pub mod annotate;
pub mod engine;
pub mod contamination;
pub mod erc;
pub mod events;
pub mod genotyper;
pub mod javamath;
pub mod likelihoods;
pub mod mathutils;
pub mod realign;
pub mod stats;
pub mod vcf;
pub mod well;

pub use engine::{resolve_qd_jitter, EngineConfig, GenotypingEngine, RegionInput, QD_JITTER_PLACEHOLDER};
pub use javamath::{JavaRandom, GATK_RANDOM_SEED};
pub use hc_assembly::{Haplotype, Read};
pub use vcf::{Phase, VcfRecord};
