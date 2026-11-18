//! `EventMap`, `Event`, the merged per-locus VariantContext (`getVariantsFromActiveHaplotypes`,
//! `replaceSpanDels`, `makeMergedVariantContext` / `simpleMerge`) and `createAlleleMapper`.

use std::collections::{BTreeMap, BTreeSet, HashSet};

use hc_assembly::Haplotype;

/// htsjdk `Allele`: equality is bases + reference flag. The spanning deletion is `*` (non-ref).
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Allele {
    pub bases: Vec<u8>,
    pub is_ref: bool,
}

impl Allele {
    pub fn new(bases: &[u8], is_ref: bool) -> Self {
        Allele { bases: bases.to_vec(), is_ref }
    }
    pub fn span_del() -> Self {
        Allele { bases: b"*".to_vec(), is_ref: false }
    }
    #[inline]
    pub fn is_span_del(&self) -> bool {
        !self.is_ref && self.bases == b"*"
    }
    #[inline]
    pub fn len(&self) -> usize {
        self.bases.len()
    }
    pub fn display(&self) -> String {
        String::from_utf8_lossy(&self.bases).into_owned()
    }
}

/// `BaseUtils.isRegularBase`: AaCcGgTt and the wildcard `*`.
#[inline]
pub fn is_regular_base(b: u8) -> bool {
    matches!(b, b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't' | b'*')
}
fn all_regular(bs: &[u8]) -> bool {
    bs.iter().all(|&b| is_regular_base(b))
}

/// `utils.haplotype.Event` (minimal representation applied on construction).
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Event {
    pub start: i64,
    pub ref_bases: Vec<u8>,
    pub alt_bases: Vec<u8>,
}

impl Event {
    pub fn new(start: i64, r: &[u8], a: &[u8]) -> Self {
        // makeMinimalRepresentation
        if r.len() == 1 || a.len() == 1 || r.is_empty() || a.is_empty() || r[r.len() - 1] != a[a.len() - 1] {
            return Event { start, ref_bases: r.to_vec(), alt_bases: a.to_vec() };
        }
        assert!(r != a, "ref and alt alleles are identical");
        let min_len = r.len().min(a.len());
        let mut overlap = 0;
        while overlap < min_len && r[r.len() - 1 - overlap] == a[a.len() - 1 - overlap] {
            overlap += 1;
        }
        Event { start, ref_bases: r[..r.len() - overlap].to_vec(), alt_bases: a[..a.len() - overlap].to_vec() }
    }
    #[inline]
    pub fn end(&self) -> i64 {
        self.start + self.ref_bases.len() as i64 - 1
    }
    #[inline]
    pub fn is_snp(&self) -> bool {
        self.ref_bases.len() == 1 && self.alt_bases.len() == 1
    }
    #[inline]
    pub fn is_simple_insertion(&self) -> bool {
        self.ref_bases.len() == 1 && self.alt_bases.len() > 1
    }
    #[inline]
    pub fn is_simple_deletion(&self) -> bool {
        self.ref_bases.len() > 1 && self.alt_bases.len() == 1
    }
}

/// `EventMap.makeCompoundEvents`
fn make_compound(e1: &Event, e2: &Event) -> Event {
    if e1.is_snp() || e2.is_snp() {
        assert!(!(e1.is_snp() && e2.is_snp()), "Trying to put two overlapping SNPs in one EventMap");
        let (snp, indel) = if e1.is_snp() { (e1, e2) } else { (e2, e1) };
        if snp.ref_bases == indel.ref_bases {
            let mut alt = snp.alt_bases.clone();
            alt.extend_from_slice(&indel.alt_bases[1..]);
            Event::new(snp.start, &snp.ref_bases, &alt)
        } else {
            Event::new(snp.start, &indel.ref_bases, &snp.alt_bases)
        }
    } else {
        assert!(
            (e1.is_simple_deletion() && e2.is_simple_insertion()) || (e1.is_simple_insertion() && e2.is_simple_deletion()),
            "Can only merge single insertion with deletion (or vice versa)"
        );
        let (ins, del) = if e1.is_simple_insertion() { (e1, e2) } else { (e2, e1) };
        Event::new(e1.start, &del.ref_bases, &ins.alt_bases)
    }
}

/// `EventMap`: events of one haplotype keyed by start.
#[derive(Clone, Debug, Default)]
pub struct EventMap {
    pub events: BTreeMap<i64, Event>,
}

impl EventMap {
    fn add(&mut self, e: Event) {
        let s = e.start;
        if let Some(old) = self.events.get(&s) {
            crate::stats::inc(&crate::stats::COMPOUND_EVENTS);
            let c = make_compound(old, &e);
            self.events.insert(s, c);
        } else {
            self.events.insert(s, e);
        }
    }

    /// `EventMap.fromHaplotype(haplotype, ref, refLoc, maxMnpDistance)`
    pub fn from_haplotype(h: &Haplotype, reference: &[u8], ref_loc_start: i64, max_mnp_distance: i64) -> EventMap {
        let mut m = EventMap::default();
        let mut ref_pos = h.align_start_wrt_ref;
        if ref_pos < 0 {
            return m;
        }
        let alignment = &h.bases;
        let mut align_pos: i64 = 0;
        let n = h.cigar.len();
        let mut proposed: Vec<Event> = Vec::new();
        for (ci, &(op, len)) in h.cigar.iter().enumerate() {
            let len = len as i64;
            match op {
                b'I' => {
                    if ref_pos > 0 && ci > 0 && ci < n - 1 {
                        let start = ref_loc_start + ref_pos - 1;
                        let ref_byte = reference[(ref_pos - 1) as usize];
                        let mut ins = vec![ref_byte];
                        ins.extend_from_slice(&alignment[align_pos as usize..(align_pos + len) as usize]);
                        if all_regular(&ins) {
                            proposed.push(Event::new(start, &[ref_byte], &ins));
                        }
                    }
                    align_pos += len;
                }
                b'S' => align_pos += len,
                b'D' => {
                    if ref_pos > 0 {
                        let del = &reference[(ref_pos - 1) as usize..(ref_pos + len) as usize];
                        let start = ref_loc_start + ref_pos - 1;
                        let ref_byte = reference[(ref_pos - 1) as usize];
                        if is_regular_base(ref_byte) && all_regular(del) {
                            proposed.push(Event::new(start, del, &[ref_byte]));
                        }
                    }
                    ref_pos += len;
                }
                b'M' | b'=' | b'X' => {
                    let mut mism: std::collections::VecDeque<i64> = std::collections::VecDeque::new();
                    for off in 0..len {
                        let rb = reference[(ref_pos + off) as usize];
                        let ab = alignment[(align_pos + off) as usize];
                        if rb != ab && is_regular_base(rb) && is_regular_base(ab) {
                            mism.push_back(off);
                        }
                    }
                    while let Some(start) = mism.pop_front() {
                        let mut end = start;
                        while let Some(&p) = mism.front() {
                            if p - end <= max_mnp_distance {
                                end = p;
                                mism.pop_front();
                            } else {
                                break;
                            }
                        }
                        let r = &reference[(ref_pos + start) as usize..(ref_pos + end + 1) as usize];
                        let a = &alignment[(align_pos + start) as usize..(align_pos + end + 1) as usize];
                        proposed.push(Event::new(ref_loc_start + ref_pos + start, r, a));
                    }
                    ref_pos += len;
                    align_pos += len;
                }
                _ => panic!("Unsupported cigar operator created during SW alignment: {}", op as char),
            }
        }
        for e in proposed {
            m.add(e);
        }
        m
    }

    /// `EventMap.getOverlappingEvents(loc)`
    pub fn overlapping(&self, loc: i64) -> Vec<&Event> {
        let ov: Vec<&Event> = self.events.range(..=loc).map(|(_, e)| e).filter(|e| e.end() >= loc).collect();
        if ov.iter().any(|e| e.is_simple_insertion()) {
            ov.into_iter().filter(|e| !(e.is_simple_deletion() && e.end() == loc)).collect()
        } else {
            ov
        }
    }
}

/// `EventMap.getEventStartPositions`
pub fn event_start_positions(maps: &[EventMap]) -> BTreeSet<i64> {
    let mut s = BTreeSet::new();
    for m in maps {
        s.extend(m.events.keys().copied());
    }
    s
}

/// The merged VariantContext at one locus: ref allele first, then alts in merge order.
#[derive(Clone, Debug)]
pub struct MergedVc {
    pub start: i64,
    pub alleles: Vec<Allele>,
}

impl MergedVc {
    pub fn end(&self) -> i64 {
        self.start + self.alleles[0].len() as i64 - 1
    }
    pub fn reference(&self) -> &Allele {
        &self.alleles[0]
    }
}

/// `getVariantsFromActiveHaplotypes(loc, haplotypes, includeSpanningEvents=true)` +
/// `replaceSpanDels` + `makeMergedVariantContext`. Returns `None` when there are no events.
pub fn merged_vc_at(loc: i64, maps: &[EventMap], reference: &[u8], ref_loc_start: i64) -> Option<MergedVc> {
    merged_vc_at_with(loc, maps, reference, ref_loc_start, true)
}

/// [`merged_vc_at`] with `includeSpanningEvents` as a parameter (Mutect2 passes false: only events
/// starting at `loc`).
pub fn merged_vc_at_with(loc: i64, maps: &[EventMap], reference: &[u8], ref_loc_start: i64, include_spanning: bool) -> Option<MergedVc> {
    let mut seen: HashSet<&Event> = HashSet::new();
    // (ref allele, alt allele) of each VC in order
    let mut vcs: Vec<(Allele, Allele)> = Vec::new();
    for m in maps {
        for e in m.overlapping(loc) {
            if !include_spanning && e.start != loc {
                continue;
            }
            if seen.insert(e) {
                if e.start == loc {
                    vcs.push((Allele::new(&e.ref_bases, true), Allele::new(&e.alt_bases, false)));
                } else {
                    let rb = reference[(loc - ref_loc_start) as usize];
                    vcs.push((Allele::new(&[rb], true), Allele::span_del()));
                }
            }
        }
    }
    if vcs.is_empty() {
        return None;
    }
    // determineReferenceAllele: the longest; equal lengths must be identical
    let mut ref_allele: Option<&Allele> = None;
    for (r, _) in &vcs {
        ref_allele = match ref_allele {
            None => Some(r),
            Some(cur) if cur.len() < r.len() => Some(r),
            Some(cur) => {
                if cur.len() == r.len() {
                    assert!(cur == r, "inconsistent reference alleles at {}", loc);
                }
                Some(cur)
            }
        };
    }
    let ref_allele = ref_allele.unwrap().clone();
    // LinkedHashSet of AlleleMapper values
    let mut set: Vec<Allele> = Vec::new();
    let push = |a: Allele, set: &mut Vec<Allele>| {
        if !set.contains(&a) {
            set.push(a);
        }
    };
    for (r, a) in &vcs {
        if *r == ref_allele {
            push(r.clone(), &mut set);
            push(a.clone(), &mut set);
        } else {
            // createAlleleMapping: alts extended (spanning deletion kept), then vc.ref -> refAllele
            let extra = &ref_allele.bases[r.len()..];
            let mapped = if a.is_span_del() {
                a.clone()
            } else {
                let mut b = a.bases.clone();
                b.extend_from_slice(extra);
                Allele { bases: b, is_ref: false }
            };
            push(mapped, &mut set);
            push(ref_allele.clone(), &mut set);
        }
    }
    // VariantContext puts the reference allele first, others keep their order
    let mut alleles = Vec::with_capacity(set.len());
    alleles.push(ref_allele.clone());
    for a in set {
        if !a.is_ref {
            alleles.push(a);
        }
    }
    Some(MergedVc { start: loc, alleles })
}

/// `AssemblyBasedCallerUtils.createAlleleMapper(mergedVC, loc, haplotypes, emitSpanningDels=true)`.
/// Returns, for each allele of the result (in order), the haplotype indices supporting it. The
/// allele list is the merged VC's alleles (all are non-symbolic), plus `*` appended if needed.
pub fn create_allele_mapper(vc: &MergedVc, loc: i64, maps: &[EventMap]) -> (Vec<Allele>, Vec<Vec<usize>>) {
    let mut alleles: Vec<Allele> = vc.alleles.clone();
    let mut haps: Vec<Vec<usize>> = vec![Vec::new(); alleles.len()];
    let merged_ref = &vc.alleles[0];
    for (h, m) in maps.iter().enumerate() {
        let ov = m.overlapping(loc);
        if ov.is_empty() {
            haps[0].push(h);
            continue;
        }
        for e in ov {
            if e.start == loc {
                let alt = if e.ref_bases.len() == merged_ref.len() {
                    Allele::new(&e.alt_bases, false)
                } else if e.ref_bases.len() < merged_ref.len() {
                    let mut b = e.alt_bases.clone();
                    b.extend_from_slice(&merged_ref.bases[e.ref_bases.len()..]);
                    Allele { bases: b, is_ref: false }
                } else {
                    continue;
                };
                if let Some(i) = alleles.iter().position(|a| *a == alt) {
                    haps[i].push(h);
                }
            } else {
                let sd = Allele::span_del();
                let i = match alleles.iter().position(|a| *a == sd) {
                    Some(i) => i,
                    None => {
                        alleles.push(sd);
                        haps.push(Vec::new());
                        alleles.len() - 1
                    }
                };
                haps[i].push(h);
                break;
            }
        }
    }
    (alleles, haps)
}
