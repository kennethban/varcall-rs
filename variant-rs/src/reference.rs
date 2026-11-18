//! Reference bases from an indexed FASTA, as GATK's `CachingIndexedFastaSequenceFile` returns them:
//! upper-cased, IUPAC ambiguity codes (R Y M K W S B D H V, any case) converted to `N`.

use std::fs::File;
use std::os::unix::fs::FileExt;

#[derive(Clone)]
pub struct FaiEntry {
    pub len: i64,
    pub offset: u64,
    pub line_bases: i64,
    pub line_width: i64,
}

pub struct Fasta {
    file: File,
    pub fai: std::collections::HashMap<String, FaiEntry>,
}

const fn iupac_table() -> [u8; 256] {
    let mut t = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        t[i] = (i as u8).to_ascii_uppercase();
        i += 1;
    }
    let codes = b"RYMKWSBDHVrymkwsbdhv";
    let mut j = 0;
    while j < codes.len() {
        t[codes[j] as usize] = b'N';
        j += 1;
    }
    t
}
static CONVERT: [u8; 256] = iupac_table();

impl Fasta {
    pub fn open(path: &str) -> Fasta {
        let file = File::open(path).unwrap_or_else(|e| panic!("{path}: {e}"));
        let text = std::fs::read_to_string(format!("{path}.fai")).expect("reference .fai index");
        let mut fai = std::collections::HashMap::new();
        for l in text.lines() {
            let f: Vec<&str> = l.split('\t').collect();
            fai.insert(
                f[0].to_string(),
                FaiEntry { len: f[1].parse().unwrap(), offset: f[2].parse().unwrap(), line_bases: f[3].parse().unwrap(), line_width: f[4].parse().unwrap() },
            );
        }
        Fasta { file, fai }
    }

    /// Bases `from..=to` (1-based, clipped to the contig) and the contig length.
    pub fn fetch(&self, contig: &str, from: i64, to: i64) -> (Vec<u8>, i64) {
        let e = self.fai.get(contig).unwrap_or_else(|| panic!("contig {contig} not in the reference .fai"));
        let from = from.max(1);
        let to = to.min(e.len);
        if to < from {
            return (Vec::new(), e.len);
        }
        let off_of = |p: i64| e.offset as i64 + (p - 1) / e.line_bases * e.line_width + (p - 1) % e.line_bases;
        let (a, b) = (off_of(from), off_of(to));
        let mut raw = vec![0u8; (b - a + 1) as usize];
        self.file.read_exact_at(&mut raw, a as u64).expect("read reference");
        let mut seq = Vec::with_capacity((to - from + 1) as usize);
        for c in raw {
            if c != b'\n' && c != b'\r' {
                seq.push(CONVERT[c as usize]);
            }
        }
        assert_eq!(seq.len() as i64, to - from + 1, "reference line layout does not match the .fai");
        (seq, e.len)
    }
}
