//! GATK's standard BQSR covariates (ReadGroup, QualityScore, Context, Cycle) for base substitutions,
//! and the read-group identifiers from the BAM header.

use std::collections::HashMap;

pub const MAX_QUAL: usize = 93; // QualityUtils.MAX_SAM_QUAL_SCORE
pub const QUAL_DIM: usize = MAX_QUAL + 1;
pub const MISMATCHES_CONTEXT_SIZE: usize = 2;
pub const INDELS_CONTEXT_SIZE: usize = 3;
pub const MAX_CYCLE: i32 = 500;
pub const LOW_QUAL_TAIL: u8 = 2;
pub const PRESERVE_QSCORES_LESS_THAN: u8 = 6;
const LENGTH_BITS: u32 = 4;

/// `ContextCovariate.maximumKeyValue() + 1` (context table dimension).
pub const CONTEXT_DIM: usize = {
    let len = if MISMATCHES_CONTEXT_SIZE > INDELS_CONTEXT_SIZE { MISMATCHES_CONTEXT_SIZE } else { INDELS_CONTEXT_SIZE };
    let mut key = len;
    let mut off = LENGTH_BITS;
    let mut i = 0;
    while i < len {
        key |= 3 << off;
        off += 2;
        i += 1;
    }
    key + 1
};
/// `CycleCovariate.maximumKeyValue() + 1`.
pub const CYCLE_DIM: usize = ((MAX_CYCLE as usize) << 1) + 2;

/// Read groups: `ReadGroupCovariate` keys (identifier = PU, else ID; first occurrence order).
pub struct ReadGroups {
    pub names: Vec<String>,
    pub key_by_id: HashMap<Vec<u8>, usize>,
    pub missing_pl: Vec<String>,
}

impl ReadGroups {
    pub fn from_header(text: &str) -> ReadGroups {
        let mut names: Vec<String> = Vec::new();
        let mut key_of_name: HashMap<String, usize> = HashMap::new();
        let mut key_by_id = HashMap::new();
        let mut n_header = 0;
        let mut missing_pl = Vec::new();
        for l in text.lines().filter(|l| l.starts_with("@RG\t")) {
            n_header += 1;
            let mut id = None;
            let mut pu = None;
            let mut pl = None;
            for f in l.split('\t').skip(1) {
                if let Some(v) = f.strip_prefix("ID:") {
                    id = Some(v.to_string());
                } else if let Some(v) = f.strip_prefix("PU:") {
                    pu = Some(v.to_string());
                } else if let Some(v) = f.strip_prefix("PL:") {
                    pl = Some(v.to_string());
                }
            }
            let id = id.expect("@RG without ID");
            if pl.is_none() {
                missing_pl.push(id.clone());
            }
            let name = pu.unwrap_or_else(|| id.clone());
            let k = *key_of_name.entry(name.clone()).or_insert_with(|| {
                names.push(name.clone());
                names.len() - 1
            });
            key_by_id.insert(id.into_bytes(), k);
        }
        assert!(n_header >= 1, "Number of read groups must be >= 1, but is 0");
        let _ = n_header;
        ReadGroups { names, key_by_id, missing_pl }
    }
}

#[inline]
pub fn base_index(b: u8) -> i32 {
    match b {
        b'A' | b'a' | b'*' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => -1,
    }
}

#[inline]
fn simple_complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'a' => b't',
        b'C' => b'G',
        b'c' => b'g',
        b'G' => b'C',
        b'g' => b'c',
        b'T' => b'A',
        b't' => b'a',
        x => x,
    }
}

const fn create_mask(context_size: usize) -> i32 {
    let mut mask = 0i32;
    let mut i = 0;
    while i < context_size {
        mask = (mask << 2) | 3;
        i += 1;
    }
    mask << LENGTH_BITS
}

pub fn key_from_context(dna: &[u8]) -> i32 {
    let mut key = dna.len() as i32;
    let mut off = LENGTH_BITS;
    for &b in dna {
        let bi = base_index(b);
        if bi == -1 {
            return -1;
        }
        key |= bi << off;
        off += 2;
    }
    key
}

/// `ContextCovariate.getReadContextAtEachPosition` -> `out` (cleared first).
fn read_context_at_each_position(bases: &[u8], context_size: usize, mask: i32, out: &mut Vec<i32>) {
    out.clear();
    let read_length = bases.len();
    let mut i = 1;
    while i < context_size && i <= read_length {
        out.push(-1);
        i += 1;
    }
    if read_length < context_size {
        return;
    }
    let new_base_offset = 2 * (context_size as u32 - 1) + LENGTH_BITS;
    let mut current_key = key_from_context(&bases[0..context_size]);
    out.push(current_key);
    let mut current_n_penalty: usize = 0;
    if current_key == -1 {
        current_key = 0;
        current_n_penalty = context_size - 1;
        let mut offset = new_base_offset as i32;
        loop {
            let bi = base_index(bases[current_n_penalty]);
            if bi == -1 {
                break;
            }
            current_key |= bi << offset;
            offset -= 2;
            current_n_penalty -= 1;
        }
    }
    for &b in &bases[context_size..] {
        let bi = base_index(b);
        if bi == -1 {
            current_n_penalty = context_size;
            current_key = 0;
        } else {
            current_key = (current_key >> 2) & mask;
            current_key |= bi << new_base_offset;
            current_key |= context_size as i32;
        }
        if current_n_penalty == 0 {
            out.push(current_key);
        } else {
            current_n_penalty -= 1;
            out.push(-1);
        }
    }
}

/// `ContextCovariate.contextFromKey`.
pub fn context_from_key(key: usize) -> String {
    let length = key & 15;
    let mut s = String::with_capacity(length);
    let mut mask = 48usize;
    let mut off = LENGTH_BITS;
    for _ in 0..length {
        let bi = (key & mask) >> off;
        s.push(b"ACGT"[bi] as char);
        mask <<= 2;
        off += 2;
    }
    s
}

/// `CycleCovariate.cycleFromKey`.
pub fn cycle_from_key(key: usize) -> i32 {
    let mut c = (key >> 1) as i32;
    if key & 1 != 0 {
        c = -c;
    }
    c
}

#[inline]
pub fn key_from_cycle(cycle: i32) -> usize {
    let r = cycle.abs();
    assert!(
        r <= MAX_CYCLE,
        "The maximum allowed value for the cycle is {MAX_CYCLE}, but a larger cycle ({r}) was detected"
    );
    ((r as usize) << 1) + (cycle < 0) as usize
}

/// Scratch buffers for one read.
#[derive(Default)]
pub struct Scratch {
    stranded: Vec<u8>,
    keys: Vec<i32>,
    pub context: Vec<i32>,
    pub cycle: Vec<u32>,
}

/// Context (mismatch) and cycle (substitution) covariate keys of a read with `bases`, `quals`, SAM
/// `flag`, into `s.context` / `s.cycle` (per read offset).
pub fn compute(bases: &[u8], quals: &[u8], flag: u16, s: &mut Scratch) {
    let len = bases.len();
    let neg = flag & 0x10 != 0;
    // ContextCovariate: low-quality tails written as N (ReadClipper.clipLowQualEnds, WRITE_NS),
    // stranded (reverse complemented for reverse-strand reads)
    s.stranded.clear();
    let mut right = len as i64 - 1;
    // signed comparisons, as Java's bytes: a quality byte >= 128 is negative, so low
    while right >= 0 && (quals[right as usize] as i8) <= LOW_QUAL_TAIL as i8 {
        right -= 1;
    }
    let mut left = 0i64;
    while (left as usize) < len && (quals[left as usize] as i8) <= LOW_QUAL_TAIL as i8 {
        left += 1;
    }
    s.context.clear();
    s.context.resize(len, 0);
    if left <= right {
        s.stranded.extend_from_slice(bases);
        for b in &mut s.stranded[..left as usize] {
            *b = b'N';
        }
        for b in &mut s.stranded[(right + 1) as usize..] {
            *b = b'N';
        }
        if neg {
            s.stranded.reverse();
            for b in s.stranded.iter_mut() {
                *b = simple_complement(*b);
            }
        }
        const MASK: i32 = create_mask(MISMATCHES_CONTEXT_SIZE);
        read_context_at_each_position(&s.stranded, MISMATCHES_CONTEXT_SIZE, MASK, &mut s.keys);
        for i in 0..len {
            let off = if neg { len - 1 - i } else { i };
            s.context[off] = s.keys[i];
        }
    }
    // (an all-low-quality read becomes empty: every context key is 0; all its bases are skipped)
    // CycleCovariate
    s.cycle.clear();
    let second = flag & 0x1 != 0 && flag & 0x80 != 0;
    let rof: i32 = if second { -1 } else { 1 };
    let (mut cycle, inc) = if neg { (len as i32 * rof, -rof) } else { (rof, rof) };
    for _ in 0..len {
        s.cycle.push(key_from_cycle(cycle) as u32);
        cycle += inc;
    }
}
