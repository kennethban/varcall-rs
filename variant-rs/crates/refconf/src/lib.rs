//! GATK 4.7 HaplotypeCaller `-ERC GVCF`: the reference confidence model (`ReferenceConfidenceModel`) and
//! the GVCF writer (`GVCFWriter`, `GVCFBlockCombiner`, `HomRefBlock`). See README.md.

pub mod gvcf;
pub mod model;
pub mod stats;

pub use gvcf::{default_gq_bands, CallRecord, GvcfWriter, Item};
pub use model::{calculate_ref_confidence, passes_hc_read_filter, CallSpan, Entry, SiteRecord};

/// GATK's GVCF header for this run (taken verbatim from `gvcf/gatk.g.vcf.gz`).
pub const GVCF_HEADER: &str = include_str!("../tables/gvcf_header.txt");
/// Header of the unbanded records (`GATK_RAW_CALLS_VCF`).
pub const RAW_HEADER: &str = include_str!("../tables/raw_header.txt");

/// `referenceModelForNoVariation` after the (optional) `finalizeRegion`: `filterNonPassingReads`, the
/// reference haplotype over the padded span, and `calculateRefConfidence` with the dummy read map.
/// `reads` are the region's (finalized) reads, `ref_padded_span` the reference bases of `padded_start..`.
pub fn ref_model_for_no_variation(reads: &[hc_assembly::Read], span_start: i64, span_end: i64, padded_start: i64, ref_padded_span: &[u8]) -> Vec<SiteRecord> {
    let kept: Vec<&hc_assembly::Read> = reads.iter().filter(|r| passes_hc_read_filter(r)).collect();
    calculate_ref_confidence(ref_padded_span, padded_start, span_start, span_end, &kept, &[])
        .into_iter()
        .map(|e| match e {
            Entry::Site(s) => s,
            Entry::Call(_) => unreachable!(),
        })
        .collect()
}
