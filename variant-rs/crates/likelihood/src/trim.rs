//! `AssemblyRegionTrimmer.trim` (non-legacy), the tandem-repeat helpers it uses
//! (`TandemRepeat.getNumTandemRepeatUnits`, `GATKVariantContextUtils.getNumTandemRepeatUnits`,
//! `findRepeatedSubstring`, `findNumberOfRepetitions`), `Haplotype.trim` and
//! `AssemblyResultSet.trimTo`.

use crate::events::Event;
use crate::{Haplotype, Params, Span};
use hc_assembly::cigar::{self, CigarBuilder, CigarElem};

// ------------------------------------------------------------------------------------------------
// Tandem repeats
// ------------------------------------------------------------------------------------------------

/// `GATKVariantContextUtils.findNumberOfRepetitions(repeatUnit, testString, leadingRepeats)`
/// on sub-slices (the offset version: empty test string gives 0).
pub fn find_number_of_repetitions(unit: &[u8], test: &[u8], leading: bool) -> usize {
    if test.is_empty() {
        return 0;
    }
    let ul = unit.len();
    let diff = test.len() as i64 - ul as i64;
    let mut n = 0;
    if leading {
        let mut start: i64 = 0;
        while start <= diff {
            let s = start as usize;
            if &test[s..s + ul] == unit {
                n += 1;
            } else {
                return n;
            }
            start += ul as i64;
        }
    } else {
        let mut start = diff;
        while start >= 0 {
            let s = start as usize;
            if &test[s..s + ul] == unit {
                n += 1;
            } else {
                return n;
            }
            start -= ul as i64;
        }
    }
    n
}

/// `GATKVariantContextUtils.findRepeatedSubstring(bases)`.  Pieces past the end are compared
/// zero-padded (`Arrays.copyOfRange`), so a partial last piece never matches.
pub fn find_repeated_substring(bases: &[u8]) -> usize {
    let n = bases.len();
    let mut rep = 1;
    while rep <= n {
        let unit = &bases[..rep];
        let mut all = true;
        let mut start = rep;
        while start < n {
            if start + rep > n || &bases[start..start + rep] != unit {
                all = false;
                break;
            }
            start += rep;
        }
        if all {
            return rep;
        }
        rep += 1;
    }
    rep
}

/// `GATKVariantContextUtils.getNumTandemRepeatUnits(byte[] refBases, byte[] altBases, byte[] remainingRefContext)`
/// -> (repetition counts [ref, alt], repeat unit).
fn tandem_repeat_units_bases(ref_bases: &[u8], alt_bases: &[u8], remaining: &[u8]) -> ([i64; 2], Vec<u8>) {
    let long_b = if alt_bases.len() > ref_bases.len() { alt_bases } else { ref_bases };
    let ul = find_repeated_substring(long_b);
    let mut unit = long_b[..ul.min(long_b.len())].to_vec();
    unit.resize(ul, 0); // Arrays.copyOf pads with zeros (unreachable for indels)
    let reps_in_ref = find_number_of_repetitions(&unit, ref_bases, true) as i64;
    let mut buf = Vec::with_capacity(ref_bases.len().max(alt_bases.len()) + remaining.len());
    buf.extend_from_slice(ref_bases);
    buf.extend_from_slice(remaining);
    let r0 = find_number_of_repetitions(&unit, &buf, true) as i64 - reps_in_ref;
    buf.clear();
    buf.extend_from_slice(alt_bases);
    buf.extend_from_slice(remaining);
    // can be negative (e.g. a deletion of a single unit not followed by another copy); Java keeps it
    let r1 = find_number_of_repetitions(&unit, &buf, true) as i64 - reps_in_ref;
    ([r0, r1], unit)
}

/// `TandemRepeat.getNumTandemRepeatUnits(ReferenceContext, Event)` followed by the trimmer's use
/// of it: returns `Some((max repeat count over ref/alt, unit length))` when GATK returns a non-null
/// result, `None` otherwise.  `ctx` is the reference context window (the region's padded span)
/// starting at `ctx_start`.
pub fn event_tandem_repeat(ctx: &[u8], ctx_start: i64, e: &Event) -> Option<(i64, usize)> {
    if e.alt_bases.len() == e.ref_bases.len() {
        return None; // only indels are tandem repeats
    }
    let start_index = e.start + 1 - ctx_start;
    assert!(start_index >= 0 && start_index as usize <= ctx.len(), "event outside reference context");
    let remaining = &ctx[start_index as usize..];
    let (counts, unit) = tandem_repeat_units_bases(&e.ref_bases[1..], &e.alt_bases[1..], remaining);
    let (c0, c1) = (counts[0], counts[1]);
    if c0 == 0 || c1 == 0 {
        return None;
    }
    Some((c0.max(c1), unit.len()))
}

// ------------------------------------------------------------------------------------------------
// AssemblyRegionTrimmer.trim
// ------------------------------------------------------------------------------------------------

/// Result of `AssemblyRegionTrimmer.trim`: `None` = no variation (`isVariationPresent() == false`),
/// else (variant span, padded variant span).
pub fn trim_region(region: Span, padded: Span, events: &[Event], ctx: &[u8], ctx_start: i64, p: &Params) -> Option<(Span, Span)> {
    let mut min_start = i64::MAX;
    let mut max_end = i64::MIN;
    let mut any = false;
    for e in events {
        if overlaps(region.start, region.end, e.start, e.end()) {
            any = true;
            min_start = min_start.min(e.start);
            max_end = max_end.max(e.end());
        }
    }
    if !any {
        return None;
    }
    let variant = intersect(Span { start: min_start, end: max_end }, region);
    for e in events {
        if !overlaps(region.start, region.end, e.start, e.end()) {
            continue;
        }
        let mut padding = p.snp_padding;
        if e.is_indel() {
            padding = p.indel_padding;
            if let Some((most_repeats, unit_len)) = event_tandem_repeat(ctx, ctx_start, e) {
                padding = p.str_padding + most_repeats * unit_len as i64;
            }
        }
        min_start = min_start.min((e.start - padding).max(1));
        max_end = max_end.max(e.end() + padding);
    }
    let padded_variant = intersect(Span { start: min_start, end: max_end }, padded);
    Some((variant, padded_variant))
}

/// htsjdk `CoordMath.overlaps` (Locatable.overlaps on one contig).
#[inline]
pub fn overlaps(start: i64, end: i64, start2: i64, end2: i64) -> bool {
    (start2 >= start && start2 <= end) || (end2 >= start && end2 <= end) || (start >= start2 && end <= end2)
}

/// `SimpleInterval.intersect` (caller guarantees overlap).
#[inline]
pub fn intersect(a: Span, b: Span) -> Span {
    assert!(overlaps(a.start, a.end, b.start, b.end), "the two intervals need to overlap");
    Span { start: a.start.max(b.start), end: a.end.min(b.end) }
}

// ------------------------------------------------------------------------------------------------
// Haplotype.trim / AssemblyResultSet.trimTo
// ------------------------------------------------------------------------------------------------

/// `AlignmentUtils.getBasesCoveringRefInterval(refStart, refEnd, bases, 0, cigar)`;
/// `None` when either end falls in a deletion.
fn bases_covering_ref_interval(ref_start: i64, ref_end: i64, bases: &[u8], cig: &[CigarElem]) -> Option<(usize, usize)> {
    let mut ref_pos: i64 = 0;
    let mut bases_pos: i64 = 0;
    let mut bstart: i64 = -1;
    let mut bstop: i64 = -1;
    'outer: for &(op, len) in cig {
        match op {
            b'I' => bases_pos += len as i64,
            b'M' | b'X' | b'=' => {
                for _ in 0..len {
                    if ref_pos == ref_start {
                        bstart = bases_pos;
                    }
                    if ref_pos == ref_end {
                        bstop = bases_pos;
                        break 'outer;
                    }
                    ref_pos += 1;
                    bases_pos += 1;
                }
            }
            b'D' => {
                for _ in 0..len {
                    if ref_pos == ref_end || ref_pos == ref_start {
                        return None;
                    }
                    ref_pos += 1;
                }
            }
            _ => panic!("Unsupported operator {}", op as char),
        }
    }
    assert!(bstart != -1 && bstop != -1, "Never found start or stop");
    let _ = bases;
    Some((bstart as usize, bstop as usize + 1))
}

/// `AlignmentUtils.trimCigarByReference(cigar, start, end).getCigar()` (CigarBuilder removes
/// leading/trailing deletions).
pub fn trim_cigar_by_reference(cig: &[CigarElem], start: i64, end: i64) -> Vec<CigarElem> {
    assert!(start >= 0 && end >= start);
    let mut nb = CigarBuilder::new(true);
    let mut element_end: i64 = 0;
    for &(op, len) in cig {
        let element_start = element_end;
        element_end = element_start + if cigar::consumes_ref(op) { len as i64 } else { 0 };
        if element_end < start || (element_end == start && element_start < start) {
            continue;
        } else if element_start > end && element_end > end + 1 {
            break;
        }
        let overlap = if element_end == element_start { len as i64 } else { (end + 1).min(element_end) - start.max(element_start) };
        nb.add((op, overlap as u32));
    }
    assert!(element_end > end, "cigar elements don't reach end position");
    nb.make(false)
}

/// `Haplotype.trim(loc, ignoreRefState = true)`: `None` when the haplotype cannot be trimmed
/// (an end falls in a deletion, or nothing but insertion remains).
pub fn trim_haplotype(h: &Haplotype, loc: Span) -> Option<Haplotype> {
    assert!(h.loc_start <= loc.start && loc.end <= h.loc_end, "Can only trim a Haplotype to a containing span");
    let new_start = loc.start - h.loc_start;
    let new_stop = new_start + loc.end - loc.start;
    let (b0, b1) = bases_covering_ref_interval(new_start, new_stop, &h.bases, &h.cigar)?;
    if b1 <= b0 {
        return None;
    }
    let new_cigar = trim_cigar_by_reference(&h.cigar, new_start, new_stop);
    let leading_ins = !cigar::consumes_ref(new_cigar[0].0);
    let trailing_ins = !cigar::consumes_ref(new_cigar[new_cigar.len() - 1].0);
    let first = if leading_ins { 1 } else { 0 };
    let last = new_cigar.len() - if trailing_ins { 1 } else { 0 };
    if last <= first {
        return None;
    }
    let c = if !(leading_ins || trailing_ins) {
        new_cigar
    } else {
        let mut b = CigarBuilder::new(false);
        b.add_all(&new_cigar[first..last]);
        b.make(false)
    };
    let bases = h.bases[b0..b1].to_vec();
    // Haplotype.setCigar: CigarBuilder(false) normalisation + read-length check
    let c = cigar::haplotype_normalize(&c);
    assert_eq!(cigar::read_length(&c), bases.len() as i64, "Read length not equal to the read length of the cigar");
    Some(Haplotype {
        bases,
        cigar: c,
        align_start_wrt_ref: new_start + h.align_start_wrt_ref,
        loc_start: loc.start,
        loc_end: loc.end,
        is_ref: false,
        score: h.score,
        kmer: h.kmer,
    })
}

/// `AssemblyResultSet.trimTo(region).getHaplotypeList()`: trims every haplotype to the padded
/// span, merges haplotypes that become equal in bases (the reference wins, otherwise the first one
/// is kept), restores the reference flag (the restored reference haplotype has kmer size 0), and
/// sorts by `Haplotype.SIZE_AND_BASE_ORDER` (length, then bases).
pub fn trim_haplotypes(haplotypes: &[Haplotype], padded: Span) -> Vec<Haplotype> {
    // keyed by bases; value = (trimmed, original-is-ref)
    let mut out: Vec<(Haplotype, bool)> = Vec::with_capacity(haplotypes.len());
    let mut index: std::collections::HashMap<Vec<u8>, usize> = std::collections::HashMap::with_capacity(haplotypes.len());
    for h in haplotypes {
        match trim_haplotype(h, padded) {
            Some(t) => {
                if let Some(&i) = index.get(&t.bases) {
                    if h.is_ref {
                        // remove(trimmed) + put(trimmed, h): the key becomes this trimmed haplotype
                        out[i] = (t, true);
                    }
                } else {
                    index.insert(t.bases.clone(), out.len());
                    out.push((t, h.is_ref));
                }
            }
            None => {
                assert!(!h.is_ref, "trimming eliminates the reference haplotype");
            }
        }
    }
    let mut res: Vec<Haplotype> = out
        .into_iter()
        .map(|(mut t, was_ref)| {
            if was_ref {
                // new Haplotype(bases, true) + cigar, location, score, alignment start (kmer size not copied)
                t.is_ref = true;
                t.kmer = 0;
            }
            t
        })
        .collect();
    res.sort_by(|a, b| a.bases.len().cmp(&b.bases.len()).then_with(|| a.bases.cmp(&b.bases)));
    assert!(res.iter().any(|h| h.is_ref), "missing reference haplotype in the trimmed set");
    res
}
