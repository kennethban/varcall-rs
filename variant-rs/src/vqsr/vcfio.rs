//! VCF reading for VQSR: a line reader (BGZF inflated in parallel; plain or gzip text otherwise) and the htsjdk
//! semantics VQSR depends on: `VariantContext.getType()`, filtered / variant / polymorphic records, and attribute
//! values as `CommonInfo.getAttributeAsDouble` reads them.

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::sync::mpsc::Receiver;
use std::sync::Arc;

/// Lines of a VCF (header included), without the newline.
pub enum Lines {
    Bgzf { rx: Receiver<Option<Vec<crate::bgzf::Block>>>, cur: Vec<u8>, pos: usize, done: bool },
    Text(Box<dyn BufRead + Send>),
}

fn is_bgzf(path: &str) -> bool {
    let mut h = [0u8; 18];
    let Ok(mut f) = File::open(path) else { return false };
    f.read_exact(&mut h).is_ok() && h[0] == 31 && h[1] == 139 && h[3] & 4 != 0 && h[12] == b'B' && h[13] == b'C'
}

impl Lines {
    /// A BGZF file's lines from a virtual offset (a `.tbi` seek): the block at `voff >> 16`, from byte `voff & 0xffff`.
    pub fn open_at(path: &str, voff: u64) -> Lines {
        let f = File::open(path).unwrap_or_else(|e| crate::gatk_cli::refuse(&format!("{path}: {e}")));
        let len = f.metadata().map(|m| m.len()).unwrap_or(0);
        let (coff, uoff) = (voff >> 16, (voff & 0xffff) as usize);
        let rx = crate::bgzf::spawn_block_reader(Arc::new(f), vec![(coff, len.saturating_sub(1))], true, Arc::new(Default::default()));
        let mut cur = Vec::new();
        if let Ok(Some(blocks)) = rx.recv() {
            for b in blocks {
                cur.extend_from_slice(&b.data);
            }
        }
        Lines::Bgzf { rx, pos: uoff.min(cur.len()), cur, done: false }
    }

    pub fn open(path: &str) -> Lines {
        let f = File::open(path).unwrap_or_else(|e| crate::gatk_cli::refuse(&format!("{path}: {e}")));
        if is_bgzf(path) {
            let len = f.metadata().map(|m| m.len()).unwrap_or(0);
            let rx = crate::bgzf::spawn_block_reader(Arc::new(f), vec![(0, len.saturating_sub(1))], true, Arc::new(Default::default()));
            Lines::Bgzf { rx, cur: Vec::new(), pos: 0, done: false }
        } else if path.ends_with(".gz") || path.ends_with(".bgz") {
            Lines::Text(Box::new(BufReader::with_capacity(1 << 20, flate2::read::MultiGzDecoder::new(f))))
        } else {
            Lines::Text(Box::new(BufReader::with_capacity(1 << 20, f)))
        }
    }

    /// The next line into `out` (cleared first); false at the end.
    pub fn next_line(&mut self, out: &mut Vec<u8>) -> bool {
        out.clear();
        match self {
            Lines::Text(r) => {
                if r.read_until(b'\n', out).expect("VCF read error") == 0 {
                    return false;
                }
            }
            Lines::Bgzf { rx, cur, pos, done } => loop {
                if let Some(k) = memchr::memchr(b'\n', &cur[*pos..]) {
                    out.extend_from_slice(&cur[*pos..*pos + k]);
                    *pos += k + 1;
                    return true;
                }
                out.extend_from_slice(&cur[*pos..]);
                cur.clear();
                *pos = 0;
                if *done {
                    return !out.is_empty();
                }
                match rx.recv() {
                    Ok(Some(blocks)) => {
                        for b in blocks {
                            cur.extend_from_slice(&b.data);
                        }
                    }
                    _ => *done = true,
                }
            },
        }
        if out.last() == Some(&b'\n') {
            out.pop();
        }
        if out.last() == Some(&b'\r') {
            out.pop();
        }
        true
    }
}

/// The virtual offset to read a contig from `start` (1-based) on, by a tabix index's linear index: its entry for the
/// 16 kb window of `start` is at or before every record that overlaps the window. None: no record of the contig
/// reaches `start`. Panics on a malformed index.
pub fn tbi_seek(tbi: &str, contig: &str, start: i64) -> Option<u64> {
    let mut raw = Vec::new();
    flate2::read::MultiGzDecoder::new(File::open(tbi).unwrap_or_else(|e| crate::gatk_cli::refuse(&format!("{tbi}: {e}")))).read_to_end(&mut raw).expect("tabix index");
    let mut o = 0usize;
    let i32_ = |o: &mut usize| -> i32 {
        let v = i32::from_le_bytes(raw[*o..*o + 4].try_into().unwrap());
        *o += 4;
        v
    };
    assert_eq!(&raw[0..4], b"TBI\x01", "{tbi}: not a tabix index");
    o += 4;
    let n_ref = i32_(&mut o) as usize;
    o += 24; // format, col_seq, col_beg, col_end, meta, skip
    let l_nm = i32_(&mut o) as usize;
    let names: Vec<&[u8]> = raw[o..o + l_nm].split(|&c| c == 0).filter(|n| !n.is_empty()).collect();
    let tid = names.iter().position(|n| *n == contig.as_bytes());
    o += l_nm;
    for r in 0..n_ref {
        let n_bin = i32_(&mut o);
        for _ in 0..n_bin {
            o += 4; // bin
            let n_chunk = i32_(&mut o) as usize;
            o += n_chunk * 16;
        }
        let n_intv = i32_(&mut o) as usize;
        if Some(r) == tid {
            let lin: Vec<u64> = (0..n_intv).map(|k| u64::from_le_bytes(raw[o + 8 * k..o + 8 * k + 8].try_into().unwrap())).collect();
            let w = ((start - 1).max(0) >> 14) as usize;
            if w >= lin.len() {
                return None;
            }
            // an empty window holds 0: the closest earlier window's offset is at or before it
            return Some(lin[..=w].iter().rev().copied().find(|&x| x != 0).unwrap_or(0));
        }
        o += n_intv * 8;
    }
    None
}

/// `VariantContext.Type`
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum VType {
    NoVariation,
    Snp,
    Mnp,
    Indel,
    Symbolic,
    Mixed,
}

/// htsjdk `Allele.wouldBeSymbolicAllele`
pub fn is_symbolic(a: &[u8]) -> bool {
    a.len() > 1 && (a[0] == b'<' || a[a.len() - 1] == b'>' || a.contains(&b'[') || a.contains(&b']') || a[0] == b'.' || a[a.len() - 1] == b'.')
}

/// `VariantContext.getType()` from REF and the ALT column
pub fn variant_type(refa: &[u8], alts: &[u8]) -> VType {
    if alts == b"." || alts.is_empty() {
        return VType::NoVariation;
    }
    let mut t: Option<VType> = None;
    for a in alts.split(|&c| c == b',') {
        let bt = if is_symbolic(a) {
            VType::Symbolic
        } else if a.len() == refa.len() {
            if a.len() == 1 { VType::Snp } else { VType::Mnp }
        } else {
            VType::Indel
        };
        match t {
            None => t = Some(bt),
            Some(x) if x != bt => return VType::Mixed,
            _ => {}
        }
    }
    t.unwrap_or(VType::NoVariation)
}

/// `Double.parseDouble`, then htsjdk's `parseVcfDouble` fallback (`[+-]?(INF|INFINITY|NAN)`, any case).
pub fn parse_vcf_double(s: &str) -> Option<f64> {
    let t = s.trim_matches(|c: char| c <= ' ');
    let body = t.strip_prefix(['+', '-']).unwrap_or(t);
    let lower = body.to_ascii_lowercase();
    if lower == "nan" {
        return Some(f64::NAN);
    }
    if lower == "inf" || lower == "infinity" {
        return Some(if t.starts_with('-') { f64::NEG_INFINITY } else { f64::INFINITY });
    }
    // Java accepts a trailing type suffix (1.5d, 2f); Rust's parser accepts "inf" / "nan" forms Java's does not,
    // but those were handled above
    let num = t.strip_suffix(['d', 'D', 'f', 'F']).unwrap_or(t);
    if num.is_empty() || num.bytes().any(|c| !(c.is_ascii_digit() || matches!(c, b'.' | b'e' | b'E' | b'+' | b'-'))) {
        return None;
    }
    num.parse::<f64>().ok()
}

/// One INFO value as `getAttributeAsDouble(key, NaN)` reads it: absent -> NaN; a String parsed (a number format
/// error is NaN, as `decodeAnnotation` catches it). A list or a flag is a ClassCastException in GATK: Err.
pub fn info_double(info: &[u8], key: &str) -> Result<f64, String> {
    let mut val: Option<Option<&[u8]>> = None;
    if info != b"." {
        for kv in info.split(|&c| c == b';') {
            let (k, v) = match memchr::memchr(b'=', kv) {
                Some(i) => (&kv[..i], Some(&kv[i + 1..])),
                None => (kv, None),
            };
            if k == key.as_bytes() {
                val = Some(v); // a later duplicate replaces the value, as in htsjdk's attribute map
            }
        }
    }
    match val {
        None => Ok(f64::NAN),
        Some(None) => Err(format!("INFO {key} is a flag, not a number")),
        Some(Some(v)) => {
            if v.contains(&b',') {
                return Err(format!("INFO {key} has several values"));
            }
            let s = std::str::from_utf8(v).unwrap_or("");
            Ok(parse_vcf_double(s).unwrap_or(f64::NAN))
        }
    }
}

/// The FILTER column: Some(filters) when filtered (`VariantContext.isFiltered`), None for "." / "PASS".
pub fn filters(col: &[u8]) -> Option<Vec<&[u8]>> {
    if col == b"." || col == b"PASS" {
        None
    } else {
        Some(col.split(|&c| c == b';').collect())
    }
}

/// `VariantContext.isPolymorphicInSamples()` for a variant record with sample columns (`cols[8..]`): some unfiltered
/// genotype has a called allele that is not the reference (called alleles != reference alleles).
pub fn polymorphic_in_samples(cols: &[&[u8]]) -> bool {
    if cols.len() <= 9 {
        // no samples: hasGenotypes() is false, so never monomorphic
        return true;
    }
    let format: Vec<&[u8]> = cols[8].split(|&c| c == b':').collect();
    let gt = format.iter().position(|&f| f == b"GT");
    let ft = format.iter().position(|&f| f == b"FT");
    let Some(gt) = gt else { return false };
    // called != reference alleles <=> some unfiltered genotype calls a non-reference allele
    for s in &cols[9..] {
        if let Some(ft) = ft {
            if let Some(f) = s.split(|&c| c == b':').nth(ft) {
                if f != b"." && f != b"PASS" {
                    continue; // a filtered genotype counts for nothing
                }
            }
        }
        let Some(g) = s.split(|&c| c == b':').nth(gt) else { continue };
        if g.split(|&c| c == b'/' || c == b'|').any(|a| a != b"." && a != b"0") {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn doubles() {
        assert_eq!(parse_vcf_double("1.5"), Some(1.5));
        assert_eq!(parse_vcf_double("-0.693"), Some(-0.693));
        assert_eq!(parse_vcf_double("2d"), Some(2.0));
        assert!(parse_vcf_double("NaN").unwrap().is_nan());
        assert!(parse_vcf_double("-inf").unwrap() == f64::NEG_INFINITY);
        assert_eq!(parse_vcf_double("."), None);
        assert_eq!(parse_vcf_double("1e-3"), Some(1e-3));
        assert!(info_double(b"AC=1;QD=2.5;FS=0.000", "FS").unwrap() == 0.0);
        assert!(info_double(b"AC=1", "FS").unwrap().is_nan());
        assert!(info_double(b"AC=1,2", "AC").is_err());
        assert!(info_double(b"DB;AC=1", "DB").is_err());
    }

    #[test]
    fn types() {
        assert_eq!(variant_type(b"A", b"G"), VType::Snp);
        assert_eq!(variant_type(b"A", b"G,*"), VType::Snp);
        assert_eq!(variant_type(b"AT", b"A"), VType::Indel);
        assert_eq!(variant_type(b"A", b"G,AT"), VType::Mixed);
        assert_eq!(variant_type(b"A", b"<NON_REF>"), VType::Symbolic);
        assert_eq!(variant_type(b"A", b"."), VType::NoVariation);
    }

    #[test]
    fn polymorphic() {
        let row = |s: &'static str| -> Vec<&'static [u8]> { s.split('\t').map(|x| x.as_bytes()).collect() };
        assert!(!polymorphic_in_samples(&row("1\t5\t.\tA\tG\t.\t.\t.\tGT\t0/0\t0|0\t./.")));
        assert!(polymorphic_in_samples(&row("1\t5\t.\tA\tG\t.\t.\t.\tGT\t0/0\t0|1")));
        assert!(!polymorphic_in_samples(&row("1\t5\t.\tA\tG\t.\t.\t.\tGT:FT\t0/0:PASS\t0/1:lowQ")));
        assert!(polymorphic_in_samples(&row("1\t5\t.\tA\tG\t.\t.\t.")));
    }
}
