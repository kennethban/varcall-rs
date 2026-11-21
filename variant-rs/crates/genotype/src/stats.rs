//! Code-path hit counters (for the validator's coverage report). Relaxed atomics; negligible cost.
use std::sync::atomic::{AtomicU64, Ordering};

pub struct Counter(pub &'static str, pub AtomicU64);

macro_rules! counters {
    ($($name:ident),*) => {
        $(pub static $name: Counter = Counter(stringify!($name), AtomicU64::new(0));)*
        pub fn all() -> Vec<&'static Counter> { vec![$(&$name),*] }
    };
}

counters!(
    LOCI_EVALUATED,
    CALLS,
    MULTIALLELIC_LOCI,
    SPAN_DEL_LOCI,
    SPAN_DEL_SPURIOUS,
    SPAN_DEL_EMITTED,
    ONLY_SPAN_DEL_DROPPED,
    MAX_ALT_ALLELES_REDUCTION,
    TOO_MANY_GENOTYPES_REDUCTION,
    REVERSE_TRIMMED,
    ALLELES_DROPPED,
    MWU_EXACT,
    MWU_NORMAL,
    FS_TABLE_NORMALIZED,
    NO_CALL_GENOTYPE,
    QD_JITTER,
    SW_FAILED_KEEP_ORIGINAL,
    REALIGN_EXACT_MATCH,
    REALIGN_FULL_SW,
    FILTERED_READS_IN_ANNOTATION,
    COMPOUND_EVENTS,
    ERC_NONREF_MEDIAN_NAN,
    ERC_NONREF_RECOMPUTED_FOR_CALL_ALLELES,
    ERC_MONOMORPHIC_EMITTED_NONREF_FIRST,
    ERC_NON_REF_CALL_TO_HOM_REF,
    ERC_HOM_REF_CALL,
    ERC_NO_CALL,
    ERC_PHASED_CALLS,
    ERC_PHASING_ABORTED,
    ERC_PHASED_GT_REVERSED
);

#[inline]
pub fn inc(c: &Counter) {
    c.1.fetch_add(1, Ordering::Relaxed);
}

pub fn report() -> String {
    all().iter().map(|c| format!("{}={}", c.0, c.1.load(Ordering::Relaxed))).collect::<Vec<_>>().join(" ")
}
