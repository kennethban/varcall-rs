//! Variation events of assembled haplotypes: `EventMap.getEvents` / `EventMap.addEvent` /
//! `Event` (minimal representation) and `AssemblyResultSet.getVariationEvents(maxMnpDistance)`.

use crate::Haplotype;
use std::cmp::Ordering;

/// A biallelic variation event (GATK `Event`), 1-based.  `ref_bases`/`alt_bases` are already in
/// GATK's minimal representation (common suffix removed when both alleles are longer than 1 base).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Event {
    pub start: i64,
    pub ref_bases: Vec<u8>,
    pub alt_bases: Vec<u8>,
}

impl Event {
    /// `new Event(contig, start, ref, alt)`: applies `makeMinimalRepresentation`.
    pub fn new(start: i64, mut r: Vec<u8>, mut a: Vec<u8>) -> Event {
        let different_last = r.is_empty() || a.is_empty() || r[r.len() - 1] != a[a.len() - 1];
        if !(r.len() == 1 || a.len() == 1 || different_last) {
            assert!(r != a, "ref and alt alleles are identical");
            let min_len = r.len().min(a.len());
            let mut overlap = 0;
            while overlap < min_len && r[r.len() - 1 - overlap] == a[a.len() - 1 - overlap] {
                overlap += 1;
            }
            r.truncate(r.len() - overlap);
            a.truncate(a.len() - overlap);
        }
        Event { start, ref_bases: r, alt_bases: a }
    }
    #[inline]
    pub fn end(&self) -> i64 {
        self.start + self.ref_bases.len() as i64 - 1
    }
    #[inline]
    pub fn is_snp(&self) -> bool {
        self.ref_bases.len() == 1 && self.alt_bases.len() == 1
    }
    /// alt alleles are never symbolic here
    #[inline]
    pub fn is_indel(&self) -> bool {
        self.ref_bases.len() != self.alt_bases.len()
    }
    #[inline]
    pub fn is_simple_insertion(&self) -> bool {
        self.ref_bases.len() == 1 && self.alt_bases.len() > 1
    }
    #[inline]
    pub fn is_simple_deletion(&self) -> bool {
        self.ref_bases.len() > 1 && self.alt_bases.len() == 1
    }
    #[inline]
    pub fn is_mnp(&self) -> bool {
        self.ref_bases.len() > 1 && self.ref_bases.len() == self.alt_bases.len()
    }
}

/// `AssemblyResultSet.HAPLOTYPE_EVENT_COMPARATOR`: start, ref length, then alt allele
/// (`Allele.compareTo` = base-string order for two non-reference alleles).
pub fn haplotype_event_cmp(a: &Event, b: &Event) -> Ordering {
    a.start
        .cmp(&b.start)
        .then(a.ref_bases.len().cmp(&b.ref_bases.len()))
        .then_with(|| a.alt_bases.cmp(&b.alt_bases))
}

/// `BaseUtils.isRegularBase`: ACGT (either case) and `*`.
#[inline]
pub fn is_regular_base(b: u8) -> bool {
    matches!(b, b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't' | b'*')
}
#[inline]
fn all_regular(bs: &[u8]) -> bool {
    bs.iter().all(|&b| is_regular_base(b))
}

/// `EventMap.makeCompoundEvents(e1, e2)` (same start).
fn make_compound(e1: Event, e2: Event) -> Event {
    if e1.is_snp() || e2.is_snp() {
        assert!(!(e1.is_snp() && e2.is_snp()), "Trying to put two overlapping SNPs in one EventMap");
        let (snp, indel) = if e1.is_snp() { (e1, e2) } else { (e2, e1) };
        if snp.ref_bases == indel.ref_bases {
            // SNP + insertion: A->G + A->CT = A->GT
            let mut alt = snp.alt_bases.clone();
            alt.extend_from_slice(&indel.alt_bases[1..]);
            Event::new(snp.start, snp.ref_bases, alt)
        } else {
            // SNP + deletion: A->T + AC->A = AC->T
            Event::new(snp.start, indel.ref_bases, snp.alt_bases)
        }
    } else {
        assert!(
            (e1.is_simple_deletion() && e2.is_simple_insertion()) || (e1.is_simple_insertion() && e2.is_simple_deletion()),
            "Can only merge single insertion with deletion (or vice versa)"
        );
        let start = e1.start;
        let (ins, del) = if e1.is_simple_insertion() { (e1, e2) } else { (e2, e1) };
        Event::new(start, del.ref_bases, ins.alt_bases)
    }
}

/// `EventMap.fromHaplotype(haplotype, ref, refLoc, maxMnpDistance).getEvents()`: the haplotype's
/// events in start order (one per start; events at the same start are compounded).
/// `reference` is the full padded reference (`getFullReferenceWithPadding`) starting at `ref_start`.
pub fn haplotype_events(h: &Haplotype, reference: &[u8], ref_start: i64, max_mnp_distance: usize) -> Vec<Event> {
    let mut proposed: Vec<Event> = Vec::new();
    let mut ref_pos = h.align_start_wrt_ref;
    if ref_pos < 0 {
        return proposed;
    }
    let alignment = &h.bases;
    let mut align_pos: usize = 0;
    let n = h.cigar.len();
    let mut mismatches: Vec<usize> = Vec::new();
    for (ci, &(op, len)) in h.cigar.iter().enumerate() {
        let len_u = len as usize;
        match op {
            b'I' => {
                if ref_pos > 0 && ci > 0 && ci < n - 1 {
                    let ref_byte = reference[(ref_pos - 1) as usize];
                    let mut ins = Vec::with_capacity(len_u + 1);
                    ins.push(ref_byte);
                    ins.extend_from_slice(&alignment[align_pos..align_pos + len_u]);
                    if all_regular(&ins) {
                        proposed.push(Event::new(ref_start + ref_pos - 1, vec![ref_byte], ins));
                    }
                }
                align_pos += len_u;
            }
            b'S' => align_pos += len_u,
            b'D' => {
                if ref_pos > 0 {
                    let rp = ref_pos as usize;
                    let del = reference[rp - 1..rp + len_u].to_vec();
                    let ref_byte = reference[rp - 1];
                    if is_regular_base(ref_byte) && all_regular(&del) {
                        proposed.push(Event::new(ref_start + ref_pos - 1, del, vec![ref_byte]));
                    }
                }
                ref_pos += len as i64;
            }
            b'M' | b'=' | b'X' => {
                let rp = ref_pos as usize;
                mismatches.clear();
                for off in 0..len_u {
                    let rb = reference[rp + off];
                    let ab = alignment[align_pos + off];
                    if rb != ab && is_regular_base(rb) && is_regular_base(ab) {
                        mismatches.push(off);
                    }
                }
                let mut k = 0;
                while k < mismatches.len() {
                    let start = mismatches[k];
                    let mut end = start;
                    k += 1;
                    while k < mismatches.len() && mismatches[k] - end <= max_mnp_distance {
                        end = mismatches[k];
                        k += 1;
                    }
                    proposed.push(Event::new(
                        ref_start + ref_pos + start as i64,
                        reference[rp + start..rp + end + 1].to_vec(),
                        alignment[align_pos + start..align_pos + end + 1].to_vec(),
                    ));
                }
                ref_pos += len as i64;
                align_pos += len_u;
            }
            _ => panic!("Unsupported cigar operator created during SW alignment: {}", op as char),
        }
    }
    // EventMap (TreeMap by start): addEvent compounds events with the same start
    let mut map: Vec<Event> = Vec::with_capacity(proposed.len());
    for e in proposed {
        match map.binary_search_by(|x| x.start.cmp(&e.start)) {
            Ok(i) => {
                let old = std::mem::replace(&mut map[i], Event { start: 0, ref_bases: vec![], alt_bases: vec![] });
                map[i] = make_compound(old, e);
            }
            Err(i) => map.insert(i, e),
        }
    }
    map
}

/// `AssemblyResultSet.getVariationEvents(maxMnpDistance)`: union of all haplotypes' events, as a
/// set sorted by `HAPLOTYPE_EVENT_COMPARATOR` (duplicates under the comparator removed).
pub fn variation_events(haplotypes: &[Haplotype], reference: &[u8], ref_start: i64, max_mnp_distance: usize) -> Vec<Event> {
    let mut all: Vec<Event> = Vec::new();
    for h in haplotypes {
        all.extend(haplotype_events(h, reference, ref_start, max_mnp_distance));
    }
    // TreeSet insertion keeps the first of equal elements; equal elements are identical here
    all.sort_by(haplotype_event_cmp);
    all.dedup_by(|b, a| haplotype_event_cmp(a, b) == Ordering::Equal);
    all
}
