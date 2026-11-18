//! `AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype` and
//! `AlignmentUtils.createReadAlignedToRef`.

use hc_assembly::cigar::{self, Cigar, CigarBuilder};
use hc_assembly::sw::{Overhang, SwAligner, SwParameters};
use hc_assembly::{Haplotype, Read};


/// `--smith-waterman-read-to-haplotype-*` defaults (match 10, mismatch -15, open -30, extend -5).
pub const READ_TO_HAPLOTYPE_SW: SwParameters =
    SwParameters { match_value: 10, mismatch_penalty: -15, gap_open_penalty: -30, gap_extend_penalty: -5 };

/// `AssemblyBasedCallerUtils.HAPLOTYPE_ALIGNMENT_TIEBREAKING_PRIORITY`: isRef + (1 - #cigar elements).
pub fn haplotype_tiebreak_priority(h: &Haplotype) -> f64 {
    (if h.is_ref { 1.0 } else { 0.0 }) + (1.0 - h.cigar.len() as f64)
}

/// Bases left after `ReadClipper.hardClipSoftClippedBases` (soft-clipped bases removed).
pub fn bases_without_soft_clips(r: &Read) -> &[u8] {
    if r.seq.is_empty() {
        return &r.seq;
    }
    let mut read_index: i64 = 0;
    let mut cut_left: i64 = -1;
    let mut cut_right: i64 = -1;
    let mut right_tail = false;
    for &(op, len) in &r.cigar {
        if op == b'S' {
            if right_tail {
                cut_right = read_index;
            } else {
                cut_left = read_index + len as i64 - 1;
            }
        } else if op != b'H' {
            right_tail = true;
        }
        if cigar::consumes_read(op) {
            read_index += len as i64;
        }
    }
    let mut end = r.seq.len();
    let mut start = 0usize;
    if cut_right >= 0 {
        end = cut_right as usize;
    }
    if cut_left >= 0 {
        start = (cut_left + 1) as usize;
    }
    if start >= end {
        return &r.seq[0..0];
    }
    &r.seq[start..end]
}

/// `AlignmentUtils.readStartOnReferenceHaplotype`
fn read_start_on_reference_haplotype(hap_vs_ref: &[(u8, u32)], read_start_on_hap: i64) -> i64 {
    if read_start_on_hap == 0 {
        return 0;
    }
    let mut ref_consumed = 0i64;
    let mut hap_consumed = 0i64;
    for &(op, len) in hap_vs_ref {
        if cigar::consumes_ref(op) {
            ref_consumed += len as i64;
        }
        if cigar::consumes_read(op) {
            hap_consumed += len as i64;
        }
        if hap_consumed >= read_start_on_hap {
            let excess = if cigar::consumes_ref(op) { hap_consumed - read_start_on_hap } else { 0 };
            return ref_consumed - excess;
        }
    }
    panic!("Cigar doesn't reach the read start");
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum OpClass {
    M,
    I,
    D,
}
fn op_class(op: u8) -> OpClass {
    match op {
        b'M' | b'=' | b'X' => OpClass::M,
        b'I' | b'S' => OpClass::I,
        b'D' => OpClass::D,
        _ => panic!("Unexpected state {}", op as char),
    }
}

/// `AlignmentUtils.applyCigarToCigar` (result through `new CigarBuilder()`).
pub fn apply_cigar_to_cigar(first_to_second: &[(u8, u32)], second_to_third: &[(u8, u32)]) -> Cigar {
    let mut nb = CigarBuilder::new(true);
    let (n12, n23) = (first_to_second.len(), second_to_third.len());
    let (mut c12, mut c23, mut e12, mut e23) = (0usize, 0usize, 0u32, 0u32);
    while c12 < n12 && c23 < n23 {
        let elt12 = first_to_second[c12];
        let elt23 = second_to_third[c23];
        // (op13, advance12, advance23) in the order of GATK's cigarPairTransformers
        let (op13, a12, a23): (Option<u8>, u32, u32) = match (op_class(elt12.0), op_class(elt23.0)) {
            (OpClass::M, OpClass::M) => (Some(b'M'), 1, 1),
            (OpClass::M, OpClass::I) => (Some(b'I'), 1, 1),
            (OpClass::M, OpClass::D) => (Some(b'D'), 0, 1),
            (OpClass::D, OpClass::M) => (Some(b'D'), 1, 1),
            (OpClass::D, OpClass::D) => (Some(b'D'), 0, 1),
            (OpClass::D, OpClass::I) => (None, 1, 1),
            (OpClass::I, OpClass::M) => (Some(b'I'), 1, 0),
            (OpClass::I, OpClass::D) => (Some(b'I'), 1, 0),
            (OpClass::I, OpClass::I) => (Some(b'I'), 1, 0),
        };
        if let Some(op) = op13 {
            nb.add((op, 1));
        }
        e12 += a12;
        e23 += a23;
        if e12 == elt12.1 {
            c12 += 1;
            e12 = 0;
        }
        if e23 == elt23.1 {
            c23 += 1;
            e23 = 0;
        }
    }
    nb.make(false)
}

/// `AlignmentUtils.appendClippedElementsFromCigarToCigar`
fn append_clipped_elements(new_cigar: &[(u8, u32)], original: &[(u8, u32)]) -> Cigar {
    let mut first_index = 0usize;
    let mut last_index = original.len() - 1;
    let mut first = original[0];
    let mut last = original[last_index];
    let mut out: Cigar = Vec::with_capacity(new_cigar.len() + 4);
    while cigar::is_clipping(first.0) && first_index != last_index {
        out.push(first);
        first_index += 1;
        first = original[first_index];
    }
    out.extend_from_slice(new_cigar);
    let mut tail: Vec<(u8, u32)> = Vec::new();
    while cigar::is_clipping(last.0) && first_index != last_index {
        tail.push(last);
        last_index -= 1;
        last = original[last_index];
    }
    out.extend(tail.iter().rev());
    out
}

/// Result of realigning one read: new alignment start and cigar (everything else is unchanged).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Realignment {
    pub start: i64,
    pub cigar: Cigar,
}

/// `AlignmentUtils.createReadAlignedToRef`. Returns `None` when GATK returns the original read
/// (Smith-Waterman offset of -1).
pub fn create_read_aligned_to_ref(
    sw: &mut SwAligner,
    read: &Read,
    hap: &Haplotype,
    ref_hap: &Haplotype,
    reference_start: i64,
    params: &SwParameters,
) -> Option<Realignment> {
    let bases = bases_without_soft_clips(read);
    let soft_clipped = read.seq.len() - bases.len();
    // SWNativeAlignerWrapper: exact substring shortcut (SOFTCLIP), else full Smith-Waterman
    let idx = hc_assembly::sw::last_index_of(&hap.bases, bases);
    let aln = if idx != -1 {
        crate::stats::inc(&crate::stats::REALIGN_EXACT_MATCH);
        hc_assembly::sw::SwAlignment { cigar: vec![(b'M', bases.len() as u32)], offset: idx }
    } else {
        crate::stats::inc(&crate::stats::REALIGN_FULL_SW);
        sw.align(&hap.bases, bases, params, Overhang::SoftClip)
    };
    if aln.offset == -1 {
        crate::stats::inc(&crate::stats::SW_FAILED_KEEP_ORIGINAL);
        return None;
    }
    let sw_cigar = {
        let mut b = CigarBuilder::new(true);
        b.add_all(&aln.cigar);
        b.make(false)
    };
    // haplotype.getConsolidatedPaddedCigar(1000)
    let padded_hap_cigar = {
        let mut b = CigarBuilder::new(true);
        b.add_all(&hap.cigar);
        b.add((b'M', 1000));
        b.make(false)
    };
    let offset = aln.offset as i64;
    let read_start_on_ref_hap = read_start_on_reference_haplotype(&padded_hap_cigar, offset);
    let read_start_on_ref = reference_start + hap.align_start_wrt_ref + read_start_on_ref_hap;
    let hap_to_ref =
        cigar::trim_cigar_by_bases(&padded_hap_cigar, offset, cigar::read_length(&padded_hap_cigar) - 1).cigar;
    let read_to_ref = apply_cigar_to_cigar(&sw_cigar, &hap_to_ref);
    let la = cigar::left_align_indels(&read_to_ref, &ref_hap.bases, bases, read_start_on_ref_hap);
    let start = read_start_on_ref + la.leading_deletion_bases_removed as i64;
    let new_cigar = append_clipped_elements(&la.cigar, &read.cigar);
    let la_len = cigar::read_length(&la.cigar) as usize;
    assert!(
        la_len + soft_clipped == read.seq.len(),
        "realigned cigar {} has read length {} != {} for read {}",
        cigar::to_string(&la.cigar),
        la_len,
        read.seq.len() - soft_clipped,
        read.name
    );
    Some(Realignment { start, cigar: new_cigar })
}

/// `realignReadsToTheirBestHaplotype` followed by `changeEvidence`: realigns every read in place.
/// `matrix[h][r]` are the log10 read likelihoods; `haplotypes` is the likelihood allele list.
pub fn realign_reads_to_best_haplotype(
    sw: &mut SwAligner,
    reads: &mut [Read],
    haplotypes: &[Haplotype],
    matrix: &[Vec<f64>],
    ref_padded_start: i64,
) {
    realign_reads_to_best_haplotype_t(sw, reads, haplotypes, matrix, ref_padded_start, crate::likelihoods::LOG_10_INFORMATIVE_THRESHOLD);
}

/// Per read: its best haplotype (`HAPLOTYPE_ALIGNMENT_TIEBREAKING_PRIORITY`), whether that is informative
/// (`BestAllele.isInformative`: confidence > 0.2 even for natural-log likelihoods), and
/// whether the read was realigned (`createReadAlignedToRef` did not return the original read).
#[derive(Clone, Copy, Debug)]
pub struct Realigned {
    pub best: usize,
    pub informative: bool,
    pub realigned: bool,
}

/// [`realign_reads_to_best_haplotype`] with the tie threshold of the matrix's log base (Mutect2 realigns
/// after `switchToNaturalLog`: 0.2 * ln 10).
pub fn realign_reads_to_best_haplotype_t(
    sw: &mut SwAligner,
    reads: &mut [Read],
    haplotypes: &[Haplotype],
    matrix: &[Vec<f64>],
    ref_padded_start: i64,
    tie_threshold: f64,
) -> Vec<Realigned> {
    let mut out = Vec::with_capacity(reads.len());
    let ref_hap = haplotypes.iter().find(|h| h.is_ref).expect("no reference haplotype");
    let priorities: Vec<f64> = haplotypes.iter().map(haplotype_tiebreak_priority).collect();
    let mut column = vec![0f64; haplotypes.len()];
    for (r, read) in reads.iter_mut().enumerate() {
        for (h, row) in matrix.iter().enumerate() {
            column[h] = row[r];
        }
        let ba = crate::likelihoods::best_allele_t(&column, Some(&priorities), tie_threshold);
        let best = ba.allele.expect("no best haplotype");
        let mut realigned = false;
        if let Some(re) =
            create_read_aligned_to_ref(sw, read, &haplotypes[best], ref_hap, ref_padded_start, &READ_TO_HAPLOTYPE_SW)
        {
            read.start = re.start;
            read.end = re.start + cigar::reference_length(&re.cigar) - 1;
            read.cigar = re.cigar;
            realigned = true;
        }
        // BestAllele.isInformative() compares with LOG_10_INFORMATIVE_THRESHOLD whatever the log base
        out.push(Realigned { best, informative: ba.is_informative(), realigned });
    }
    out
}
