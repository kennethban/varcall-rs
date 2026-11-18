//! CIGAR utilities: port of htsjdk Cigar semantics plus GATK's CigarBuilder,
//! AlignmentUtils.{trimCigarByBases, leftAlignIndels, removeTrailingDeletions}.

/// A CIGAR element: (operator char like b'M', length).
pub type CigarElem = (u8, u32);
pub type Cigar = Vec<CigarElem>;

#[inline]
pub fn consumes_read(op: u8) -> bool {
    matches!(op, b'M' | b'I' | b'S' | b'=' | b'X')
}
#[inline]
pub fn consumes_ref(op: u8) -> bool {
    matches!(op, b'M' | b'D' | b'N' | b'=' | b'X')
}
#[inline]
pub fn is_clipping(op: u8) -> bool {
    matches!(op, b'S' | b'H')
}
#[inline]
pub fn is_indel(op: u8) -> bool {
    matches!(op, b'I' | b'D')
}
#[inline]
pub fn is_alignment(op: u8) -> bool {
    matches!(op, b'M' | b'=' | b'X')
}

pub fn reference_length(c: &[CigarElem]) -> i64 {
    c.iter().filter(|e| consumes_ref(e.0)).map(|e| e.1 as i64).sum()
}
pub fn read_length(c: &[CigarElem]) -> i64 {
    c.iter().filter(|e| consumes_read(e.0)).map(|e| e.1 as i64).sum()
}

pub fn to_string(c: &[CigarElem]) -> String {
    let mut s = String::new();
    for &(op, len) in c {
        s.push_str(&len.to_string());
        s.push(op as char);
    }
    s
}

pub fn parse(s: &str) -> Cigar {
    let mut out = Vec::new();
    let mut n: u32 = 0;
    for b in s.bytes() {
        if b.is_ascii_digit() {
            n = n * 10 + (b - b'0') as u32;
        } else {
            out.push((b, n));
            n = 0;
        }
    }
    out
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Section {
    LeftHardClip,
    LeftSoftClip,
    Middle,
    RightSoftClip,
    RightHardClip,
}

/// Port of org.broadinstitute.hellbender.utils.read.CigarBuilder
pub struct CigarBuilder {
    elems: Vec<CigarElem>,
    last_op: Option<u8>,
    section: Section,
    remove_deletions_at_ends: bool,
    leading_deletion_bases_removed: u32,
    trailing_deletion_bases_removed: u32,
    trailing_deletion_bases_removed_in_make: u32,
}

pub struct CigarBuilderResult {
    pub cigar: Cigar,
    pub leading_deletion_bases_removed: u32,
    pub trailing_deletion_bases_removed: u32,
}

impl CigarBuilder {
    pub fn new(remove_deletions_at_ends: bool) -> Self {
        CigarBuilder {
            elems: Vec::new(),
            last_op: None,
            section: Section::LeftHardClip,
            remove_deletions_at_ends,
            leading_deletion_bases_removed: 0,
            trailing_deletion_bases_removed: 0,
            trailing_deletion_bases_removed_in_make: 0,
        }
    }

    fn last_two_were_del_and_ins(&self) -> bool {
        self.last_op == Some(b'I') && self.elems.len() > 1 && self.elems[self.elems.len() - 2].0 == b'D'
    }

    pub fn add(&mut self, e: CigarElem) -> &mut Self {
        let (op, len) = e;
        if len == 0 {
            return self;
        }
        if self.remove_deletions_at_ends && op == b'D' {
            match self.last_op {
                None => {
                    self.leading_deletion_bases_removed += len;
                    return self;
                }
                Some(l) if is_clipping(l) => {
                    self.leading_deletion_bases_removed += len;
                    return self;
                }
                Some(b'I') => {
                    let n = self.elems.len();
                    if n == 1 || is_clipping(self.elems[n - 2].0) {
                        self.leading_deletion_bases_removed += len;
                        return self;
                    }
                }
                _ => {}
            }
        }
        self.advance_section(op);
        if Some(op) == self.last_op {
            let n = self.elems.len() - 1;
            self.elems[n].1 += len;
        } else if self.last_op.is_none() {
            self.elems.push(e);
            self.last_op = Some(op);
        } else if is_clipping(op) {
            let last = self.last_op.unwrap();
            if self.remove_deletions_at_ends && !consumes_read(last) && !is_clipping(last) {
                let n = self.elems.len();
                self.trailing_deletion_bases_removed += self.elems[n - 1].1;
                self.elems[n - 1] = e;
                self.last_op = Some(op);
            } else if self.remove_deletions_at_ends && self.last_two_were_del_and_ins() {
                let n = self.elems.len();
                self.trailing_deletion_bases_removed += self.elems[n - 2].1;
                self.elems[n - 2] = self.elems[n - 1];
                self.elems[n - 1] = e;
                // note: GATK does not update lastOperator here
            } else {
                self.elems.push(e);
                self.last_op = Some(op);
            }
        } else if op == b'D' && self.last_op == Some(b'I') {
            let size = self.elems.len();
            if size > 1 && self.elems[size - 2].0 == b'D' {
                self.elems[size - 2].1 += len;
            } else {
                self.elems.insert(size - 1, e);
            }
        } else {
            self.elems.push(e);
            self.last_op = Some(op);
        }
        self
    }

    fn advance_section(&mut self, op: u8) {
        if op == b'H' {
            if matches!(self.section, Section::LeftSoftClip | Section::Middle | Section::RightSoftClip) {
                self.section = Section::RightHardClip;
            }
        } else if op == b'S' {
            assert!(self.section != Section::RightHardClip, "cigar has already reached its right hard clip");
            if self.section == Section::LeftHardClip {
                self.section = Section::LeftSoftClip;
            } else if self.section == Section::Middle {
                self.section = Section::RightSoftClip;
            }
        } else {
            assert!(
                self.section != Section::RightSoftClip && self.section != Section::RightHardClip,
                "cigar has already reached right clip"
            );
            if matches!(self.section, Section::LeftHardClip | Section::LeftSoftClip) {
                self.section = Section::Middle;
            }
        }
    }

    pub fn add_all(&mut self, c: &[CigarElem]) -> &mut Self {
        for &e in c {
            self.add(e);
        }
        self
    }

    pub fn make(&mut self, _allow_empty: bool) -> Cigar {
        self.trailing_deletion_bases_removed_in_make = 0;
        if self.remove_deletions_at_ends && self.last_op == Some(b'D') {
            let n = self.elems.len();
            self.trailing_deletion_bases_removed_in_make = self.elems[n - 1].1;
            self.elems.pop();
        } else if self.remove_deletions_at_ends && self.last_two_were_del_and_ins() {
            let n = self.elems.len();
            self.trailing_deletion_bases_removed_in_make = self.elems[n - 2].1;
            self.elems.remove(n - 2);
        }
        self.elems.clone()
    }

    pub fn make_and_record(&mut self) -> CigarBuilderResult {
        let cigar = self.make(false);
        CigarBuilderResult {
            cigar,
            leading_deletion_bases_removed: self.leading_deletion_bases_removed,
            trailing_deletion_bases_removed: self.trailing_deletion_bases_removed
                + self.trailing_deletion_bases_removed_in_make,
        }
    }
}

/// Haplotype.setCigar: new CigarBuilder(false).addAll(cigar).make()
pub fn haplotype_normalize(c: &[CigarElem]) -> Cigar {
    let mut b = CigarBuilder::new(false);
    b.add_all(c);
    b.make(false)
}

/// AlignmentUtils.removeTrailingDeletions
pub fn remove_trailing_deletions(c: Cigar) -> Cigar {
    let mut c = c;
    if let Some(l) = c.last() {
        if l.0 == b'D' {
            c.pop();
        }
    }
    c
}

/// AlignmentUtils.trimCigarByBases (byReference = false)
pub fn trim_cigar_by_bases(cigar: &[CigarElem], start: i64, end: i64) -> CigarBuilderResult {
    trim_cigar(cigar, start, end, false)
}

fn trim_cigar(cigar: &[CigarElem], start: i64, end: i64, by_reference: bool) -> CigarBuilderResult {
    assert!(start >= 0 && end >= start);
    let mut nb = CigarBuilder::new(true);
    let mut element_end: i64 = 0;
    for &(op, len) in cigar {
        let element_start = element_end;
        let l = if by_reference {
            if consumes_ref(op) { len as i64 } else { 0 }
        } else if consumes_read(op) {
            len as i64
        } else {
            0
        };
        element_end = element_start + l;
        if element_end < start || (element_end == start && element_start < start) {
            continue;
        } else if element_start > end && element_end > end + 1 {
            break;
        }
        let overlap = if element_end == element_start {
            len as i64
        } else {
            (end + 1).min(element_end) - start.max(element_start)
        };
        nb.add((op, overlap as u32));
    }
    assert!(element_end > end, "cigar elements don't reach end position");
    nb.make_and_record()
}

#[derive(Clone, Copy, Debug)]
struct IndexRange {
    from: i64,
    to: i64,
}
impl IndexRange {
    fn size(&self) -> i64 {
        self.to - self.from
    }
    fn validate(&self) {
        debug_assert!(self.from <= self.to && self.from >= 0, "invalid index range");
    }
    fn shift(&mut self, s: i64) {
        self.to += s;
        self.from += s;
        self.validate();
    }
    fn shift_left(&mut self, s: i64) {
        self.shift(-s)
    }
    fn shift_start(&mut self, s: i64) {
        self.from += s;
        self.validate();
    }
    fn shift_start_left(&mut self, s: i64) {
        self.shift_start(-s)
    }
    fn shift_end(&mut self, s: i64) {
        self.to += s;
        self.validate();
    }
    fn shift_end_left(&mut self, s: i64) {
        self.shift_end(-s)
    }
}

/// AlignmentUtils.normalizeAlleles for exactly two sequences (ref, read)
fn normalize_alleles(seqs: [&[u8]; 2], bounds: &mut [IndexRange; 2], max_shift: i64, trim: bool) -> (i64, i64) {
    let mut start_shift = 0i64;
    let mut end_shift = 0i64;
    let mut min_size = bounds[0].size().min(bounds[1].size());
    let last_base_on_right_is_same = |b: &[IndexRange; 2]| -> bool {
        let x = seqs[0][(b[0].to - 1) as usize];
        seqs[1][(b[1].to - 1) as usize] == x
    };
    let first_base_on_left_is_same = |b: &[IndexRange; 2]| -> bool {
        let x = seqs[0][b[0].from as usize];
        seqs[1][b[1].from as usize] == x
    };
    let next_base_on_left_is_same = |b: &[IndexRange; 2]| -> bool {
        let x = seqs[0][(b[0].from - 1) as usize];
        seqs[1][(b[1].from - 1) as usize] == x
    };
    while trim && min_size > 0 && last_base_on_right_is_same(bounds) {
        bounds[0].shift_end_left(1);
        bounds[1].shift_end_left(1);
        min_size -= 1;
        end_shift += 1;
    }
    while trim && min_size > 0 && first_base_on_left_is_same(bounds) {
        bounds[0].shift_start(1);
        bounds[1].shift_start(1);
        min_size -= 1;
        start_shift -= 1;
    }
    while start_shift < max_shift && next_base_on_left_is_same(bounds) && last_base_on_right_is_same(bounds) {
        bounds[0].shift_left(1);
        bounds[1].shift_left(1);
        start_shift += 1;
        end_shift += 1;
    }
    (start_shift, end_shift)
}

/// AlignmentUtils.leftAlignIndels
pub fn left_align_indels(cigar: &[CigarElem], reference: &[u8], read: &[u8], read_start: i64) -> CigarBuilderResult {
    if !cigar.iter().any(|e| is_indel(e.0)) {
        return CigarBuilderResult {
            cigar: cigar.to_vec(),
            leading_deletion_bases_removed: 0,
            trailing_deletion_bases_removed: 0,
        };
    }
    let len_on_ref = |e: &CigarElem| if consumes_ref(e.0) { e.1 as i64 } else { 0 };
    let len_on_read = |e: &CigarElem| if consumes_read(e.0) { e.1 as i64 } else { 0 };
    let mut result_rtl: Vec<CigarElem> = Vec::new();
    let ref_length = reference_length(cigar);
    let mut ref_indel = IndexRange { from: read_start + ref_length, to: read_start + ref_length };
    let mut read_indel = IndexRange { from: read.len() as i64, to: read.len() as i64 };
    let ncig = cigar.len();
    for n in (0..ncig).rev() {
        let element = cigar[n];
        if is_indel(element.0) {
            ref_indel.shift_start_left(len_on_ref(&element));
            read_indel.shift_start_left(len_on_read(&element));
        } else if ref_indel.size() == 0 && read_indel.size() == 0 {
            result_rtl.push(element);
            ref_indel.shift_left(len_on_ref(&element));
            read_indel.shift_left(len_on_read(&element));
        } else {
            let max_shift = if is_alignment(element.0) { element.1 as i64 } else { 0 };
            let mut b = [ref_indel, read_indel];
            let shifts = normalize_alleles([reference, read], &mut b, max_shift, true);
            ref_indel = b[0];
            read_indel = b[1];
            result_rtl.push((b'M', shifts.1 as u32));
            let emit_indel = n == 0 || shifts.0 < max_shift || !is_alignment(element.0);
            let new_match_on_left = if shifts.0 < 0 { -shifts.0 } else { 0 };
            let remaining_on_left = if shifts.0 < 0 { element.1 as i64 } else { element.1 as i64 - shifts.0 };
            if emit_indel {
                result_rtl.push((b'D', ref_indel.size() as u32));
                result_rtl.push((b'I', read_indel.size() as u32));
                let rs = ref_indel.size();
                ref_indel.shift_end_left(rs);
                let qs = read_indel.size();
                read_indel.shift_end_left(qs);
                ref_indel.shift_left(new_match_on_left + if consumes_ref(element.0) { remaining_on_left } else { 0 });
                read_indel.shift_left(new_match_on_left + if consumes_read(element.0) { remaining_on_left } else { 0 });
            }
            result_rtl.push((b'M', new_match_on_left as u32));
            result_rtl.push((element.0, remaining_on_left as u32));
        }
    }
    result_rtl.push((b'D', ref_indel.size() as u32));
    result_rtl.push((b'I', read_indel.size() as u32));
    assert!(read_indel.from == 0, "Given cigar does not account for all bases of the read");
    let mut b = CigarBuilder::new(true);
    for e in result_rtl.iter().rev() {
        b.add(*e);
    }
    b.make_and_record()
}
