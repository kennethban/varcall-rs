//! Port of CigarUtils.calculateCigar (haplotype -> reference alignment with padding,
//! trimming and indel left-alignment).

use crate::cigar::{self, Cigar};
use crate::sw::{Overhang, FixedRefAligner, SwAligner, SwAlignment, SwParameters};

const SW_PAD: &[u8] = b"NNNNNNNNNN";

/// Returns None if the SW alignment "failed" (offset > 0 or contains soft clips).
pub fn calculate_cigar(
    aligner: &mut SwAligner,
    ref_seq: &[u8],
    alt_seq: &[u8],
    params: &SwParameters,
    strategy: Overhang,
    pad_buf: &mut (Vec<u8>, Vec<u8>),
) -> Option<Cigar> {
    calculate_cigar_impl(ref_seq, alt_seq, pad_buf, |r, a| aligner.align(r, a, params, strategy))
}

pub fn pad(seq: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.extend_from_slice(SW_PAD);
    out.extend_from_slice(seq);
    out.extend_from_slice(SW_PAD);
}

/// calculateCigar with the SOFTCLIP strategy for each of `alt_seqs`, using an aligner whose reference has been
/// set to the padded `ref_seq` (sorted alternates share the most Smith-Waterman columns).
pub fn calculate_cigars_softclip(aligner: &mut FixedRefAligner, ref_seq: &[u8], alt_seqs: &[&[u8]]) -> Vec<Option<Cigar>> {
    let mut out: Vec<Option<Cigar>> = vec![None; alt_seqs.len()];
    let mut todo: Vec<(usize, Vec<u8>)> = Vec::new();
    for (k, alt) in alt_seqs.iter().enumerate() {
        match early_cigar(ref_seq, alt) {
            Some(c) => out[k] = c,
            None => {
                let mut padded = Vec::with_capacity(alt.len() + 2 * SW_PAD.len());
                pad(alt, &mut padded);
                todo.push((k, padded));
            }
        }
    }
    let padded: Vec<&[u8]> = todo.iter().map(|(_, p)| &p[..]).collect();
    let alns = aligner.align_batch(&padded);
    for ((k, p), aln) in todo.iter().zip(alns) {
        out[*k] = cigar_from_alignment(ref_seq, alt_seqs[*k], p.len(), aln);
    }
    out
}

/// The cases calculateCigar answers without aligning.
fn early_cigar(ref_seq: &[u8], alt_seq: &[u8]) -> Option<Option<Cigar>> {
    if alt_seq.is_empty() {
        return Some(Some(vec![(b'D', ref_seq.len() as u32)]));
    }
    if ref_seq == alt_seq {
        return Some(Some(vec![(b'M', ref_seq.len() as u32)]));
    }
    None
}

fn calculate_cigar_impl(
    ref_seq: &[u8],
    alt_seq: &[u8],
    pad_buf: &mut (Vec<u8>, Vec<u8>),
    mut align: impl FnMut(&[u8], &[u8]) -> SwAlignment,
) -> Option<Cigar> {
    if let Some(c) = early_cigar(ref_seq, alt_seq) {
        return c;
    }
    let (padded_ref, padded_path) = pad_buf;
    pad(ref_seq, padded_ref);
    pad(alt_seq, padded_path);
    let alignment = align(padded_ref, padded_path);
    cigar_from_alignment(ref_seq, alt_seq, padded_path.len(), alignment)
}

/// The rest of calculateCigar, from the alignment of the padded sequences (`padded_len`: padded alternate length).
fn cigar_from_alignment(ref_seq: &[u8], alt_seq: &[u8], padded_len: usize, alignment: SwAlignment) -> Option<Cigar> {
    if alignment.offset > 0 || alignment.cigar.iter().any(|e| e.0 == b'S') {
        return None;
    }
    let base_start = SW_PAD.len() as i64;
    let base_end = padded_len as i64 - SW_PAD.len() as i64 - 1;
    let trimmed = cigar::trim_cigar_by_bases(&alignment.cigar, base_start, base_end);
    let mut non_standard = trimmed.cigar;
    let trimmed_leading = trimmed.leading_deletion_bases_removed;
    let trimmed_trailing = trimmed.trailing_deletion_bases_removed;
    if trimmed_trailing > 0 {
        non_standard.push((b'D', trimmed_trailing));
    }
    let la = cigar::left_align_indels(&non_standard, ref_seq, alt_seq, trimmed_leading as i64);
    let total_leading = trimmed_leading + la.leading_deletion_bases_removed;
    let total_trailing = la.trailing_deletion_bases_removed;
    if total_leading == 0 && total_trailing == 0 {
        Some(la.cigar)
    } else {
        let mut res = Vec::with_capacity(la.cigar.len() + 2);
        if total_leading > 0 {
            res.push((b'D', total_leading));
        }
        res.extend_from_slice(&la.cigar);
        if total_trailing > 0 {
            res.push((b'D', total_trailing));
        }
        Some(res)
    }
}
