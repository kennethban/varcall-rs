//! BaseRecalibrator: per-read error/skip arrays and covariate counts (`BaseRecalibrationEngine`).
//!
//! Without BAQ (GATK's default) every error is 0 or 1, so a table cell is just (observations,
//! errors) integer counts: per-thread tables can be summed in any order and give GATK's values.

use super::covariates::{self, Scratch, CONTEXT_DIM, CYCLE_DIM, PRESERVE_QSCORES_LESS_THAN, QUAL_DIM};
use super::known::KnownSites;
use hc_preproc::{read_index_for_ref, Read, READ_INDEX_NOT_FOUND};

/// Count tables for base substitutions: `[rg][qual]`, `[rg][qual][context]`, `[rg][qual][cycle]`;
/// each cell (observations, errors).
pub struct Counts {
    pub n_rg: usize,
    pub qual: Vec<[u64; 2]>,
    pub context: Vec<[u64; 2]>,
    pub cycle: Vec<[u64; 2]>,
    pub reads_processed: u64,
}

impl Counts {
    pub fn new(n_rg: usize) -> Counts {
        Counts {
            n_rg,
            qual: vec![[0; 2]; n_rg * QUAL_DIM],
            context: vec![[0; 2]; n_rg * QUAL_DIM * CONTEXT_DIM],
            cycle: vec![[0; 2]; n_rg * QUAL_DIM * CYCLE_DIM],
            reads_processed: 0,
        }
    }
    pub fn add(&mut self, o: &Counts) {
        for (a, b) in [(&mut self.qual, &o.qual), (&mut self.context, &o.context), (&mut self.cycle, &o.cycle)] {
            for (x, y) in a.iter_mut().zip(b.iter()) {
                x[0] += y[0];
                x[1] += y[1];
            }
        }
        self.reads_processed += o.reads_processed;
    }
}

/// Per-thread scratch.
#[derive(Default)]
pub struct Work {
    cov: Scratch,
    snp: Vec<u8>,
    known: Vec<bool>,
    pub cigar_scratch: Vec<(u8, u32)>,
}

/// BaseRecalibrator's read filters (`getStandardBQSRReadFilterList`) on a decoded record, in GATK's
/// order: MappingQualityNotZero, MappingQualityAvailable, Mapped, NotSecondaryAlignment,
/// NotDuplicate, PassesVendorQualityCheck, WellformedRead.
pub fn passes_filters(r: &Read, has_rg: bool, qual_missing: bool, ref_len: i64) -> bool {
    if r.mapq == 0 || r.mapq == 255 {
        return false;
    }
    if r.flag & 0x4 != 0 || r.ref_id < 0 || r.start < 1 {
        return false;
    }
    if r.flag & 0x100 != 0 || r.flag & 0x400 != 0 || r.flag & 0x200 != 0 {
        return false;
    }
    // WellformedReadFilter
    let cigar_read_len: usize = r.cigar.iter().filter(|e| matches!(e.0, b'M' | b'I' | b'S' | b'=' | b'X')).map(|e| e.1 as usize).sum();
    has_rg
        && !r.seq.is_empty()
        && !qual_missing
        && r.qual.len() == r.seq.len()
        // ReadLengthEqualsCigarLengthReadFilter: a mapped read with no CIGAR ('*') has CIGAR length 0 and is dropped
        && cigar_read_len == r.seq.len()
        && !r.cigar.iter().any(|e| e.0 == b'N')
        && r.start <= ref_len
        && r.end - r.start + 1 >= 0
}

/// `BaseRecalibrationEngine.processRead` after the read transform (`r` is already clipped).
/// `orig_span` is the alignment span of the read before clipping (the FeatureContext interval);
/// `refseq`/`ref_from` hold the reference bases from 1-based position `ref_from`.
pub fn process_read(
    r: &Read,
    rg: usize,
    tid: usize,
    orig_span: (i64, i64),
    known: &KnownSites,
    refseq: &[u8],
    ref_from: i64,
    c: &mut Counts,
    w: &mut Work,
) {
    let len = r.seq.len();
    // calculateIsSNPOrIndel (only the SNP array is used without indel tables)
    w.snp.clear();
    w.snp.resize(len, 0);
    let mut read_pos = 0usize;
    let mut ref_pos = (r.start - ref_from) as usize;
    for &(op, l) in &r.cigar {
        let l = l as usize;
        match op {
            b'M' | b'=' | b'X' => {
                for _ in 0..l {
                    let rb = refseq[ref_pos];
                    w.snp[read_pos] = (covariates::base_index(r.seq[read_pos]) != covariates::base_index(rb)) as u8;
                    read_pos += 1;
                    ref_pos += 1;
                }
            }
            b'D' | b'N' => ref_pos += l,
            b'I' | b'S' => read_pos += l,
            _ => {}
        }
    }
    // calculateKnownSites over the clipped read (soft clips are gone, so softStart = start)
    w.known.clear();
    w.known.resize(len, false);
    let soft_start = r.soft_start();
    let soft_end = {
        let mut e = r.end;
        for &(op, l) in r.cigar.iter().rev() {
            if op == b'S' {
                e += l as i64;
            } else if op != b'H' {
                break;
            }
        }
        e
    };
    let cigar = &r.cigar;
    let known_arr = &mut w.known;
    known.for_each_overlapping(tid, orig_span.0, orig_span.1, |ks, ke| {
        if ke < soft_start || ks > soft_end {
            return;
        }
        let (i0, op0) = read_index_for_ref(soft_start, cigar, ks);
        let mut fs = if i0 == READ_INDEX_NOT_FOUND { 0 } else { i0 };
        if op0 == Some(b'D') {
            fs -= 1;
        }
        let (i1, _) = read_index_for_ref(soft_start, cigar, ke);
        let mut fe = if i1 == READ_INDEX_NOT_FOUND { len as i64 } else { i1 };
        if fs > len as i64 {
            fs = len as i64;
            fe = len as i64;
        }
        let a = fs.max(0) as usize;
        let b = (fe + 1).min(len as i64).max(0) as usize;
        for k in a..b.max(a) {
            known_arr[k] = true;
        }
    });
    // covariates
    covariates::compute(&r.seq, &r.qual, r.flag, &mut w.cov);
    // skip array and counts
    let qbase = rg * QUAL_DIM;
    for i in 0..len {
        let q = r.qual[i];
        // signed, as Java compares bytes: >= 128 is negative and skipped
        if covariates::base_index(r.seq[i]) == -1 || (q as i8) < PRESERVE_QSCORES_LESS_THAN as i8 || w.known[i] {
            continue;
        }
        let q = q as usize;
        // GATK's NestedIntegerArray refuses the key and the run fails
        assert!(q < QUAL_DIM, "Key {q} is too large for dimension 1 (max is {})", QUAL_DIM - 1);
        let e = w.snp[i] as u64;
        let cell = &mut c.qual[qbase + q];
        cell[0] += 1;
        cell[1] += e;
        let ctx = w.cov.context[i];
        if ctx >= 0 {
            let cell = &mut c.context[(qbase + q) * CONTEXT_DIM + ctx as usize];
            cell[0] += 1;
            cell[1] += e;
        }
        let cyc = w.cov.cycle[i] as usize;
        let cell = &mut c.cycle[(qbase + q) * CYCLE_DIM + cyc];
        cell[0] += 1;
        cell[1] += e;
    }
    c.reads_processed += 1;
}
