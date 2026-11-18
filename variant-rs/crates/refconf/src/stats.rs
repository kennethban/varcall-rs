//! Code-path hit counters (printed by variant-rs in GVCF mode). Relaxed atomics; negligible cost.
use std::sync::atomic::{AtomicU64, Ordering};

pub struct Counter(pub &'static str, pub AtomicU64);

macro_rules! counters {
    ($($name:ident),*) => {
        $(pub static $name: Counter = Counter(stringify!($name), AtomicU64::new(0));)*
        pub fn all() -> Vec<&'static Counter> { vec![$(&$name),*] }
    };
}

counters!(
    SITES,
    PILEUP_ELEMENTS,
    LOW_QUAL_BASES_SKIPPED,
    DELETION_ELEMENTS,
    ADAPTOR_EXCLUDED,
    IDL_CACHE_BUILT,
    IDL_CACHE_EMPTY,
    IDL_REF_SHORTER,
    INDEL_INFORMATIVE_CAPPED,
    INDEL_GLS_CHOSEN,
    CALL_STARTS,
    SITES_UNDER_CALL,
    GVCF_SITES_DROPPED_UNDER_VARIANT,
    GVCF_BLOCKS,
    GVCF_VARIANT_RECORDS,
    GVCF_CALL_AS_HOM_REF_SITE,
    RAW_PL_ALIASED
);

#[inline]
pub fn inc(c: &Counter) {
    c.1.fetch_add(1, Ordering::Relaxed);
}

pub fn report() -> String {
    all().iter().map(|c| format!("{}={}", c.0, c.1.load(Ordering::Relaxed))).collect::<Vec<_>>().join(" ")
}
