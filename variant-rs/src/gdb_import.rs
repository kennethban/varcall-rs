//! GenomicsDBImport's cells (`docs/specs/genomicsdb_tiledb.md` sections 5-7): each GVCF record of a sample that
//! overlaps an array's interval becomes a begin cell at (row, column) holding every attribute as GenomicsDB's loader
//! stores it (REF, ALT with `<NON_REF>` as `&`, QUAL as f32, INFO / FORMAT values as their BCF encoding, GT as
//! allele / phase pairs, the AS_ "multi-D" fields as serialized vectors), plus, for a record spanning several
//! positions, an END copy at its end column (moved before the next record of the row when they overlap). Cells of
//! all rows are merged in (column, row) order and written as one fragment per batch (`gdb_write`).
//! Refused (not verified against GATK): FILTER values other than '.', multi-sample inputs.
//!
//! The tool (`run`, `docs/specs/genomicsdb_import.md`): samples from `--sample-name-map` or `-V` (rows in sorted
//! sample order), the fields from the first sample's header (or the `-V` headers merged), `vidmap.json` /
//! `callset.json` as GenomicsDB's protobuf JSON printer writes them, `vcfheader.vcf`, one array per (merged) `-L`
//! interval, one fragment per batch (one in all with `--consolidate`), `callset.json` last.

use crate::gdb_store::{Attr, FLT_MISSING_BITS, FLT_NULL_BITS, INT_MISSING, INT_NULL};
use crate::mergevcfs::{allele, format_qual};
use crate::gdb_write::{self, FragmentWriter};
use crate::vqsr::vcfio::{self, Lines};
use std::collections::{BinaryHeap, HashMap};
use std::path::Path;

fn refuse(why: &str) -> ! {
    crate::gatk_cli::refuse(why)
}

/// GATK fails on this input too: exit with its message.
fn fail(why: &str) -> ! {
    eprintln!("variant-rs GenomicsDBImport: {why}");
    std::process::exit(1)
}

/// Element type of a vidmap field: i32, f32 or char.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Elem {
    I,
    F,
    C,
}

/// One `vidmap.json` field.
#[derive(Clone, Debug)]
pub struct Field {
    pub name: String,
    pub types: Vec<Elem>,
    /// one length descriptor per dimension ("1", "A", "R", "G", "var", "PP", ...)
    pub lengths: Vec<String>,
    pub classes: Vec<String>,
    pub delimiters: Vec<String>,
    /// a Flag: one char, 1 when present (`TileDB`'s char null 0x7F when absent)
    pub flag: bool,
}

impl Field {
    fn multid(&self) -> bool {
        self.lengths.len() > 1
    }
}

/// The vidmap: fields, the global field index (FILTER values), contigs.
pub struct Vid {
    pub fields: Vec<Field>,
    pub by_name: HashMap<String, usize>,
    /// FILTER values store these (FILTER is refused for now)
    #[allow(dead_code)]
    pub field_idx: HashMap<String, i32>,
    /// (name, length, tiledb column offset)
    pub contigs: Vec<(String, i64, i64)>,
}

fn elem(t: &str) -> Elem {
    match t {
        "Integer" | "int" => Elem::I,
        "Float" | "float" => Elem::F,
        "String" | "char" | "Character" | "Flag" | "flag" => Elem::C,
        x => refuse(&format!("vidmap field type {x} is not supported")),
    }
}

impl Vid {
    pub fn from_json(v: &serde_json::Value) -> Vid {
        let s = |v: &serde_json::Value| v.as_str().map(|x| x.to_string()).unwrap_or_else(|| v.to_string());
        let mut fields = Vec::new();
        for f in v["fields"].as_array().unwrap() {
            let list = |k: &str| f[k].as_array().map(|a| a.iter().map(s).collect::<Vec<_>>()).unwrap_or_default();
            fields.push(Field {
                name: s(&f["name"]),
                types: list("type").iter().map(|t| elem(t)).collect(),
                lengths: f["length"].as_array().map(|a| a.iter().map(|l| s(&l["variable_length_descriptor"])).collect()).unwrap_or_default(),
                classes: list("vcf_field_class"),
                delimiters: list("vcf_delimiter"),
                flag: list("type").iter().any(|t| t == "Flag" || t == "flag"),
            });
        }
        // the global field index: fields in order, `<name>_FORMAT` after an INFO+FORMAT field, the tuple elements
        // after a multi-type field
        let mut field_idx = HashMap::new();
        let mut k = 0;
        for f in &fields {
            field_idx.insert(f.name.clone(), k);
            k += 1;
            if f.classes.iter().any(|c| c == "INFO") && f.classes.iter().any(|c| c == "FORMAT") {
                field_idx.insert(format!("{}_FORMAT", f.name), k);
                k += 1;
            }
            if f.types.len() > 1 {
                for t in 0..f.types.len() {
                    field_idx.insert(format!("{}_tuple_element_{t}", f.name), k);
                    k += 1;
                }
            }
        }
        let num = |v: &serde_json::Value| -> i64 { v.as_str().map(|x| x.parse().unwrap()).or_else(|| v.as_i64()).unwrap() };
        let contigs = v["contigs"].as_array().unwrap().iter().map(|c| (s(&c["name"]), num(&c["length"]), num(&c["tiledb_column_offset"]))).collect();
        let by_name = fields.iter().enumerate().map(|(i, f)| (f.name.clone(), i)).collect();
        Vid { fields, by_name, field_idx, contigs }
    }

    /// The attribute (schema) list (`VidMapper::build_tiledb_array_schema`): END, REF, ALT, ID, QUAL, FILTER, then the
    /// INFO fields, then the FORMAT fields (`<name>_FORMAT` for one that is also INFO), multi-type fields as their
    /// tuple elements; a fixed length -> that many values, other lengths variable; multi-D fields and tuple elements
    /// as CHAR blobs.
    pub fn attrs(&self) -> Vec<Attr> {
        use gdb_write::{attr, CHAR, FLOAT32, INT32, INT64};
        let var = crate::gdb_store::TILEDB_VAR_NUM;
        let mut a = vec![attr("END", INT64, 1), attr("REF", CHAR, var), attr("ALT", CHAR, var)];
        if self.by_name.contains_key("ID") {
            a.push(attr("ID", CHAR, var));
        }
        a.push(attr("QUAL", FLOAT32, 1));
        a.push(attr("FILTER", INT32, var));
        for class in ["INFO", "FORMAT"] {
            for f in &self.fields {
                if !f.classes.iter().any(|c| c == class) || (class == "INFO" && f.name == "END") {
                    continue;
                }
                let both = f.classes.iter().any(|c| c == "INFO") && f.classes.iter().any(|c| c == "FORMAT");
                let name = if class == "FORMAT" && both { format!("{}_FORMAT", f.name) } else { f.name.clone() };
                if f.types.len() > 1 {
                    for t in 0..f.types.len() {
                        a.push(attr(&format!("{name}_tuple_element_{t}"), CHAR, var));
                    }
                    continue;
                }
                if f.multid() {
                    a.push(attr(&name, CHAR, var));
                    continue;
                }
                let typ = match f.types[0] {
                    Elem::I => INT32,
                    Elem::F => FLOAT32,
                    Elem::C => CHAR,
                };
                let n = f.lengths.first().and_then(|l| l.parse::<u32>().ok()).unwrap_or(var);
                a.push(attr(&name, typ, n));
            }
        }
        a
    }

    /// (field, tuple element, from FORMAT) of an attribute
    fn field_for_attr(&self, a: &str) -> (&Field, Option<usize>, bool) {
        if let Some(base) = a.strip_suffix("_FORMAT") {
            if let Some(&i) = self.by_name.get(base) {
                return (&self.fields[i], None, true);
            }
        }
        if let Some((base, k)) = a.rsplit_once("_tuple_element_") {
            if let (Some(&i), Ok(k)) = (self.by_name.get(base), k.parse()) {
                let f = &self.fields[i];
                let fmt = f.classes.iter().any(|c| c == "FORMAT") && !f.classes.iter().any(|c| c == "INFO");
                return (f, Some(k), fmt);
            }
        }
        let f = &self.fields[*self.by_name.get(a).unwrap_or_else(|| refuse(&format!("attribute {a} has no vidmap field")))];
        let fmt = f.classes.iter().any(|c| c == "FORMAT") && !f.classes.iter().any(|c| c == "INFO");
        (f, None, fmt)
    }
}

// ---------------------------------------------------------------- values

fn ints(v: &[i32]) -> Vec<u8> {
    v.iter().flat_map(|x| x.to_le_bytes()).collect()
}

fn floats(v: &[Option<f32>]) -> Vec<u8> {
    v.iter().flat_map(|x| x.map_or(FLT_MISSING_BITS, |f| f.to_bits()).to_le_bytes()).collect()
}

/// A VCF value as htslib's VCF parser stores it (`strtod`, then to float).
fn f32_of(tok: &str, at: &str) -> f32 {
    vcfio::parse_vcf_double(tok).unwrap_or_else(|| refuse(&format!("{at}: '{tok}' is not a number"))) as f32
}

fn i32_of(tok: &str, at: &str) -> i32 {
    tok.parse::<i32>().unwrap_or_else(|_| refuse(&format!("{at}: '{tok}' is not an integer")))
}

/// C `strtoll(tok, 0)` truncated to i32 (the multi-D parser): hex, octal, decimal; trailing junk ignored.
fn strtol0(tok: &str) -> i32 {
    let t = tok.trim_start();
    let (neg, s) = match t.as_bytes().first() {
        Some(b'-') => (true, &t[1..]),
        Some(b'+') => (false, &t[1..]),
        _ => (false, t),
    };
    let (radix, digits) = if s.len() > 1 && (s.starts_with("0x") || s.starts_with("0X")) {
        (16, &s[2..])
    } else if s.len() > 1 && s.starts_with('0') {
        (8, &s[1..])
    } else {
        (10, s)
    };
    let d: String = digits.chars().take_while(|c| c.is_digit(radix)).collect();
    let v = i64::from_str_radix(&d, radix).unwrap_or(0);
    (if neg { -v } else { v }) as i32
}

/// A multi-D field's VCF string as `GenomicsDBMultiDVectorField::parse_and_store_numeric` stores it: per tuple
/// element, `u64 S | S bytes of inner vectors | u64 N | u64 offsets[N + 1]`.
fn multid(s: &str, f: &Field) -> Vec<Vec<u8>> {
    let nt = f.types.len();
    let (outer, inner) = (f.delimiters.first().map_or("|", |d| d.as_str()), f.delimiters.get(1).map_or(",", |d| d.as_str()));
    let mut per: Vec<Vec<Vec<u8>>> = vec![Vec::new(); nt];
    for g in s.split(outer) {
        let toks: Vec<&str> = g.split(inner).collect();
        let mut vecs: Vec<Vec<u8>> = vec![Vec::new(); nt];
        let mut lens = vec![0usize; nt];
        for (j, tok) in toks.iter().enumerate() {
            let k = j % nt;
            let missing = tok.is_empty() || tok.eq_ignore_ascii_case("nan");
            let b = match f.types[k] {
                Elem::I => (if missing { INT_MISSING } else { strtol0(tok) }).to_le_bytes(),
                _ => (if missing { FLT_MISSING_BITS } else { c_strtof(tok).to_bits() }).to_le_bytes(),
            };
            vecs[k].extend_from_slice(&b);
            lens[k] += 1;
        }
        let m = *lens.iter().max().unwrap();
        for k in 0..nt {
            for _ in lens[k]..m {
                let b = if f.types[k] == Elem::I { INT_MISSING.to_le_bytes() } else { FLT_MISSING_BITS.to_le_bytes() };
                vecs[k].extend_from_slice(&b);
            }
            per[k].push(std::mem::take(&mut vecs[k]));
        }
    }
    per.into_iter()
        .map(|vs| {
            let s: usize = vs.iter().map(|v| v.len()).sum();
            let mut out = Vec::with_capacity(16 + s + 8 * (vs.len() + 1));
            out.extend_from_slice(&(s as u64).to_le_bytes());
            for v in &vs {
                out.extend_from_slice(v);
            }
            out.extend_from_slice(&(vs.len() as u64).to_le_bytes());
            let mut o = 0u64;
            out.extend_from_slice(&o.to_le_bytes());
            for v in &vs {
                o += v.len() as u64;
                out.extend_from_slice(&o.to_le_bytes());
            }
            out
        })
        .collect()
}

/// C `strtof` on a whole token (decimal; a malformed token is 0).
fn c_strtof(tok: &str) -> f32 {
    tok.trim().parse::<f32>().unwrap_or(0.0)
}

/// htslib's variant type of one ALT against REF (the subset GenomicsDB's deletion-as-interval rule needs).
fn is_indel_or_overlap(r: &str, a: &str) -> bool {
    if a == "*" {
        return true;
    }
    if r.len() == 1 && a.len() == 1 {
        return false;
    }
    if a.starts_with('<') || a.contains('[') || a.contains(']') {
        return false;
    }
    let (r, a) = (r.to_ascii_uppercase().into_bytes(), a.to_ascii_uppercase().into_bytes());
    let mut i = 0;
    while i < r.len() && i < a.len() && r[i] == a[i] {
        i += 1;
    }
    if (i < a.len() && i == r.len()) || (i < r.len() && i == a.len()) {
        return true;
    }
    if i == r.len() && i == a.len() {
        return false;
    }
    let (mut re, mut ae) = (r.len() - 1, a.len() - 1);
    while re > i && ae > i && r[re] == a[ae] {
        re -= 1;
        ae -= 1;
    }
    if ae == i {
        return re != i && r[re] == a[ae];
    }
    if re == i {
        return r[re] == a[ae];
    }
    false
}

/// One parsed GVCF record.
struct Rec<'a> {
    pos: i64,
    id: &'a str,
    refa: &'a str,
    alts: Vec<&'a str>,
    qual: &'a str,
    filter: &'a str,
    info: HashMap<&'a str, Option<&'a str>>,
    fmt: HashMap<&'a str, &'a str>,
}

impl<'a> Rec<'a> {
    fn parse(line: &'a str, path: &str) -> Rec<'a> {
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 8 {
            refuse(&format!("{path}: a VCF line with fewer than 8 columns"));
        }
        if c.len() > 10 {
            refuse(&format!("{path}: GenomicsDBImport of multi-sample VCFs is not supported"));
        }
        let mut info = HashMap::new();
        if c[7] != "." {
            for kv in c[7].split(';') {
                match kv.split_once('=') {
                    Some((k, v)) => info.insert(k, Some(v)),
                    None => info.insert(kv, None),
                };
            }
        }
        let fmt = if c.len() > 9 { c[8].split(':').zip(c[9].split(':')).collect() } else { HashMap::new() };
        Rec {
            pos: c[1].parse().unwrap_or_else(|_| refuse(&format!("{path}: bad POS {}", c[1]))),
            id: c[2],
            refa: c[3],
            alts: if c[4] == "." { vec![] } else { c[4].split(',').collect() },
            qual: c[5],
            filter: c[6],
            info,
            fmt,
        }
    }

    /// `GATKVariantContextUtils.isUnmixedMnpIgnoringNonRef`: REF of 2+ bases and every ALT `<NON_REF>` or as long.
    fn unmixed_mnp(&self) -> bool {
        self.refa.len() >= 2 && self.alts.iter().all(|a| *a == "<NON_REF>" || a.len() == self.refa.len())
    }

    /// The 1-based end GenomicsDB stores: INFO END, else a deletion's end (`treat_deletions_as_intervals`), else POS.
    fn end(&self, at: &str) -> i64 {
        if let Some(Some(e)) = self.info.get("END") {
            return e.parse().unwrap_or_else(|_| refuse(&format!("{at}: bad END {e}")));
        }
        if self.alts.iter().any(|a| is_indel_or_overlap(self.refa, a) && self.refa.len() > a.len()) {
            return self.pos + self.refa.len() as i64 - 1;
        }
        self.pos
    }

    /// Every attribute's bytes but END (schema order; END is set by the caller).
    /// `htsjdk` is GATK's default path, where the record reaches GenomicsDB re-written by htsjdk's VCF writer: QUAL
    /// rounded to 2 decimals, REF / ALT upper-cased (`--bypass-feature-reader` reads the text as it is).
    fn values(&self, vid: &Vid, attrs: &[Attr], at: &str, htsjdk: bool) -> Vec<Vec<u8>> {
        if self.filter != "." {
            refuse(&format!("{at}: FILTER {} in a GVCF: GenomicsDBImport of filtered records is not supported", self.filter));
        }
        let mut out = Vec::with_capacity(attrs.len());
        let mut cache: HashMap<&str, Vec<Vec<u8>>> = HashMap::new();
        for a in attrs {
            let v = match a.name.as_str() {
                "END" => vec![],
                "REF" => if htsjdk { allele(self.refa).into_bytes() } else { self.refa.as_bytes().to_vec() },
                "ALT" => self
                    .alts
                    .iter()
                    .map(|x| if *x == "<NON_REF>" { "&".to_string() } else if htsjdk { allele(x) } else { x.to_string() })
                    .collect::<Vec<_>>()
                    .join("|")
                    .into_bytes(),
                "ID" => if self.id == "." { vec![] } else { self.id.as_bytes().to_vec() },
                "QUAL" => {
                    let q = if htsjdk { format_qual(self.qual) } else { self.qual.to_string() };
                    if q == "." { FLT_NULL_BITS.to_le_bytes().to_vec() } else { f32_of(&q, at).to_bits().to_le_bytes().to_vec() }
                }
                "FILTER" => vec![],
                name => {
                    let (f, k, is_fmt) = vid.field_for_attr(name);
                    let val: Option<&str> = if is_fmt {
                        self.fmt.get(f.name.as_str()).copied()
                    } else {
                        if f.flag {
                            // a flag: present (bare, or with a value other than 0) -> 1, absent -> char null
                            let present = matches!(self.info.get(f.name.as_str()), Some(None)) || matches!(self.info.get(f.name.as_str()), Some(Some(v)) if *v != "0");
                            out.push(vec![if present { 1 } else { 0x7f }]);
                            continue;
                        }
                        match self.info.get(f.name.as_str()) {
                            Some(Some(v)) => Some(v),
                            Some(None) => refuse(&format!("{at}: INFO {} without a value is not supported", f.name)),
                            None => None,
                        }
                    };
                    let missing = val.is_none() || val == Some(".");
                    if f.multid() || k.is_some() {
                        if missing {
                            vec![]
                        } else {
                            let e = cache.entry(f.name.as_str()).or_insert_with(|| multid(val.unwrap(), f));
                            e[k.unwrap_or(0)].clone()
                        }
                    } else if f.name == "GT" {
                        match val {
                            None => vec![],
                            Some(".") => ints(&[-1]),
                            Some(g) => {
                                let mut v = Vec::new();
                                let mut phase = 0;
                                let mut tok = String::new();
                                let flush = |tok: &str, v: &mut Vec<i32>| v.push(if tok == "." { -1 } else { i32_of(tok, at) });
                                for ch in g.chars() {
                                    if ch == '/' || ch == '|' {
                                        flush(&tok, &mut v);
                                        tok.clear();
                                        phase = (ch == '|') as i32;
                                        v.push(phase);
                                    } else {
                                        tok.push(ch);
                                    }
                                }
                                let _ = phase;
                                flush(&tok, &mut v);
                                ints(&v)
                            }
                        }
                    } else if a.typ == gdb_write::CHAR {
                        if missing { vec![] } else { val.unwrap().split('\0').next().unwrap().as_bytes().to_vec() }
                    } else {
                        let el = f.types[0];
                        let fixed = !a.is_var();
                        if missing {
                            if fixed {
                                let one = if el == Elem::I { INT_NULL.to_le_bytes() } else { FLT_NULL_BITS.to_le_bytes() };
                                one.repeat(a.cell_val_num as usize)
                            } else {
                                vec![]
                            }
                        } else {
                            let toks: Vec<&str> = val.unwrap().split(',').collect();
                            if fixed && toks.len() != a.cell_val_num as usize {
                                refuse(&format!("{at}: {} has {} values, the header says {}", f.name, toks.len(), a.cell_val_num));
                            }
                            if el == Elem::I {
                                ints(&toks.iter().map(|t| if *t == "." || t.is_empty() { INT_MISSING } else { i32_of(t, at) }).collect::<Vec<_>>())
                            } else {
                                floats(&toks.iter().map(|t| if *t == "." || t.is_empty() { None } else { Some(f32_of(t, at)) }).collect::<Vec<_>>())
                            }
                        }
                    }
                }
            };
            out.push(v);
        }
        out
    }
}

// ---------------------------------------------------------------- cells of one row

/// The cells of one sample in an array, in column order: begin cells and END copies (section 5).
struct RowCells {
    row: i64,
    path: String,
    lines: Lines,
    line: Vec<u8>,
    contig: String,
    offset: i64,
    start: i64,
    stop: i64,
    /// an END copy waiting for the next record: (column, values with END = the begin column)
    pending: Option<(i64, Vec<Vec<u8>>)>,
    /// cells ready, in column order
    ready: std::collections::VecDeque<(i64, Vec<Vec<u8>>)>,
    /// records that begin before the interval and overlap it: only the last is kept
    before: Option<(i64, i64, Vec<Vec<u8>>)>,
    done: bool,
    last_begin: i64,
    htsjdk: bool,
}

impl RowCells {
    fn fill(&mut self, vid: &Vid, attrs: &[Attr]) {
        while self.ready.is_empty() && !self.done {
            if !self.lines.next_line(&mut self.line) {
                self.finish();
                break;
            }
            if self.line.is_empty() || self.line[0] == b'#' {
                continue;
            }
            let text = std::str::from_utf8(&self.line).unwrap_or_else(|_| refuse(&format!("{}: not UTF-8", self.path))).to_string();
            let chrom = text.split('\t').next().unwrap();
            if chrom != self.contig {
                if self.last_begin != i64::MIN {
                    // past the contig (inputs are sorted)
                    self.finish();
                }
                continue;
            }
            let rec = Rec::parse(&text, &self.path);
            let at = format!("{}:{}:{}", self.path, chrom, rec.pos);
            if rec.pos > self.stop {
                self.finish();
                break;
            }
            if self.htsjdk && rec.unmixed_mnp() {
                fail(&format!(
                    "GenomicsDBImport does not support GVCFs with MNPs. MNP found at {chrom}:{} in VCF {}",
                    rec.pos,
                    std::fs::canonicalize(&self.path).map(|p| p.display().to_string()).unwrap_or(self.path.clone())
                ));
            }
            let end = rec.end(&at);
            if end < self.start {
                continue;
            }
            let b = self.offset + rec.pos - 1;
            let e = self.offset + end - 1;
            let mut vals = rec.values(vid, attrs, &at, self.htsjdk);
            vals[0] = e.to_le_bytes().to_vec();
            if rec.pos < self.start {
                // begins before the interval: kept until the first record inside it (the latest one wins)
                self.before = Some((b, e, vals));
                self.last_begin = b;
                continue;
            }
            if let Some((bb, be, bv)) = self.before.take() {
                self.emit(bb, be, bv);
            }
            if b <= self.last_begin && self.last_begin != i64::MIN {
                refuse(&format!("{at}: records out of order or two records at one position"));
            }
            self.last_begin = b;
            self.emit(b, e, vals);
        }
    }

    fn emit(&mut self, b: i64, e: i64, vals: Vec<Vec<u8>>) {
        if let Some((pc, pv)) = self.pending.take() {
            let begin = i64::from_le_bytes(pv[0][..8].try_into().unwrap());
            if pc >= b {
                if b - 1 != begin {
                    self.ready.push_back((b - 1, pv));
                }
            } else {
                self.ready.push_back((pc, pv));
            }
        }
        if e > b {
            let mut copy = vals.clone();
            copy[0] = b.to_le_bytes().to_vec();
            self.pending = Some((e, copy));
        }
        self.ready.push_back((b, vals));
    }

    fn finish(&mut self) {
        if let Some((bb, be, bv)) = self.before.take() {
            self.emit(bb, be, bv);
        }
        if let Some((pc, pv)) = self.pending.take() {
            self.ready.push_back((pc, pv));
        }
        self.done = true;
    }
}

/// One array: the cells of `rows` (sample row, GVCF path) for contig:start-stop, merged in (column, row) order and
/// written as one fragment. Returns the number of cells.
pub fn write_fragment(array_dir: &Path, vid: &Vid, attrs: &[Attr], rows: &[(i64, String)], contig: &str, start: i64, stop: i64, htsjdk: bool) -> u64 {
    let offset = vid.contigs.iter().find(|c| c.0 == contig).unwrap_or_else(|| refuse(&format!("contig {contig} is not in the vidmap"))).2;
    let mut rcs: Vec<RowCells> = rows
        .iter()
        .map(|(row, path)| RowCells {
            row: *row,
            path: path.clone(),
            // from the index's offset for `start` (GATK queries the index too); a contig the file does not reach
            // there: nothing to read
            lines: match vcfio::tbi_seek(&format!("{path}.tbi"), contig, start) {
                Some(v) => Lines::open_at(path, v),
                None => Lines::Text(Box::new(std::io::empty())),
            },
            line: Vec::new(),
            contig: contig.to_string(),
            offset,
            start,
            stop,
            pending: None,
            ready: Default::default(),
            before: None,
            done: false,
            last_begin: i64::MIN,
            htsjdk,
        })
        .collect();
    let mut w = FragmentWriter::create(array_dir, attrs).unwrap_or_else(|e| refuse(&format!("{}: {e}", array_dir.display())));
    let mut heap = BinaryHeap::new();
    for (i, rc) in rcs.iter_mut().enumerate() {
        rc.fill(vid, attrs);
        if let Some((c, _)) = rc.ready.front() {
            heap.push(std::cmp::Reverse((*c, rc.row, i)));
        }
    }
    while let Some(std::cmp::Reverse((col, row, i))) = heap.pop() {
        let (_, vals) = rcs[i].ready.pop_front().unwrap();
        w.push(row, col, vals).unwrap_or_else(|e| refuse(&format!("{}: {e}", array_dir.display())));
        rcs[i].fill(vid, attrs);
        if let Some((c, _)) = rcs[i].ready.front() {
            heap.push(std::cmp::Reverse((*c, rcs[i].row, i)));
        }
    }
    w.finish().unwrap_or_else(|e| refuse(&format!("{}: {e}", array_dir.display())))
}

// ---------------------------------------------------------------- the tool

/// Gson's HTML-safe string escaping, quoted.
fn jstr(s: &str) -> String {
    let mut o = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            '\t' => o.push_str("\\t"),
            '\u{8}' => o.push_str("\\b"),
            '\n' => o.push_str("\\n"),
            '\r' => o.push_str("\\r"),
            '\u{c}' => o.push_str("\\f"),
            '<' | '>' | '&' | '=' | '\'' | '\u{2028}' | '\u{2029}' => o.push_str(&format!("\\u{:04x}", c as u32)),
            c if (c as u32) < 0x20 => o.push_str(&format!("\\u{:04x}", c as u32)),
            c => o.push(c),
        }
    }
    o.push('"');
    o
}

/// A JSON value as protobuf-java's `JsonFormat` printer lays it out.
enum J {
    S(String),
    Raw(String),
    L(Vec<J>),
    M(Vec<(&'static str, J)>),
}

fn jprint(v: &J, ind: usize) -> String {
    match v {
        J::S(s) => jstr(s),
        J::Raw(r) => r.clone(),
        J::L(l) => format!("[{}]", l.iter().map(|x| jprint(x, ind)).collect::<Vec<_>>().join(", ")),
        J::M(m) => {
            if m.is_empty() {
                return format!("{{\n{}}}", " ".repeat(ind));
            }
            let pad = " ".repeat(ind + 2);
            let items: Vec<String> = m.iter().map(|(k, v)| format!("{pad}{}: {}", jstr(k), jprint(v, ind + 2))).collect();
            format!("{{\n{}\n{}}}", items.join(",\n"), " ".repeat(ind))
        }
    }
}

/// One vidmap field as the importer builds it (`VidMapExtensions`, GATK's overrides).
#[derive(Clone)]
struct VidField {
    name: String,
    types: Vec<String>,
    classes: Vec<String>,
    lengths: Vec<String>,
    delimiters: Vec<String>,
    combine: Option<String>,
    vcf_name: Option<String>,
    disable_remap: bool,
}

impl VidField {
    fn json(&self) -> J {
        let list = |v: &[String]| J::L(v.iter().map(|x| J::S(x.clone())).collect());
        let mut m: Vec<(&'static str, J)> = vec![("name", J::S(self.name.clone())), ("type", list(&self.types))];
        if !self.classes.is_empty() {
            m.push(("vcf_field_class", list(&self.classes)));
        }
        m.push(("length", J::L(self.lengths.iter().map(|l| J::M(vec![("variable_length_descriptor", J::S(l.clone()))])).collect())));
        if !self.delimiters.is_empty() {
            m.push(("vcf_delimiter", list(&self.delimiters)));
        }
        if let Some(c) = &self.combine {
            m.push(("VCF_field_combine_operation", J::S(c.clone())));
        }
        if let Some(n) = &self.vcf_name {
            m.push(("vcf_name", J::S(n.clone())));
        }
        if self.disable_remap {
            m.push(("disable_remap_missing_with_non_ref", J::Raw("true".into())));
        }
        J::M(m)
    }

    /// An allele-specific field: R x var, delimiters | and , (`updateAlleleSpecificINFOFieldCombineOperation`).
    fn make_as(&mut self, types: &[&str], combine: &str) {
        self.types = types.iter().map(|t| t.to_string()).collect();
        self.lengths = vec!["R".into(), "var".into()];
        self.delimiters = vec!["|".into(), ",".into()];
        self.combine = Some(combine.into());
    }
}

/// `vidmap.json`'s fields from the header lines (`VidMapExtensions.generateVidMapFromMergedHeader` + GATK's
/// `GATKGenomicsDBUtils` overrides) and contigs (offsets: the running sum of lengths).
fn build_vidmap(lines: &[String]) -> (Vec<VidField>, Vec<(String, i64, i64)>) {
    let mut fields = vec![VidField {
        name: "ID".into(),
        types: vec!["char".into()],
        classes: vec![],
        lengths: vec!["var".into()],
        delimiters: vec![],
        combine: None,
        vcf_name: None,
        disable_remap: false,
    }];
    let mut seen: HashMap<String, Vec<(String, usize)>> = HashMap::new();
    let mut contigs = Vec::new();
    let mut offset = 0i64;
    for l in lines {
        let Some((key, id, number, typ, items)) = crate::vcfheader::compound(l) else { continue };
        let (class, types, length): (&str, Vec<String>, String) = match key.as_str() {
            "contig" => {
                let len = items.iter().find_map(|f| f.strip_prefix("length=")).map_or(0, |v| v.parse::<i64>().unwrap_or(0));
                contigs.push((id, len, offset));
                offset += len;
                continue;
            }
            "FILTER" => ("FILTER", vec!["Integer".into()], "1".into()),
            "FORMAT" | "INFO" => {
                let typ = typ.unwrap_or_else(|| refuse(&format!("header line without Type: {l}")));
                if typ == "Character" {
                    fail(&format!("Unhandled field type Character for field {id} (GenomicsDB aborts the import)"));
                }
                let number = number.unwrap_or_default();
                if key == "FORMAT" && id == "GT" {
                    ("FORMAT", vec!["int".into()], "PP".into())
                } else if typ == "String" {
                    (if key == "FORMAT" { "FORMAT" } else { "INFO" }, vec![typ], "var".into())
                } else {
                    let len = match number.as_str() {
                        "A" | "R" | "G" => number.clone(),
                        "." => "VAR".into(),
                        n if n == "0" && typ == "Flag" => "1".into(),
                        n => n.to_string(),
                    };
                    (if key == "FORMAT" { "FORMAT" } else { "INFO" }, vec![typ], len)
                }
            }
            _ => continue,
        };
        let mut f = VidField {
            name: id.clone(),
            types,
            classes: vec![class.into()],
            lengths: vec![length],
            delimiters: vec![],
            combine: None,
            vcf_name: None,
            disable_remap: false,
        };
        if class == "INFO" {
            match id.as_str() {
                "AS_RAW_BaseQRankSum" | "AS_RAW_MQRankSum" | "AS_RAW_ReadPosRankSum" => f.make_as(&["float", "int"], "histogram_sum"),
                "AS_RAW_MQ" => f.make_as(&["float"], "element_wise_sum"),
                "AS_SB_TABLE" | "AS_QUALapprox" | "AS_VarDP" => f.make_as(&["int"], "element_wise_sum"),
                _ => {}
            }
        }
        // a name seen before (`checkForDuplicateVcfFieldNames`): the same class again is skipped; else the earlier
        // entry of the first class in FILTER, INFO, FORMAT order decides: same type and length -> one entry with
        // both classes, otherwise a new entry `<ID>_<CLASS>`
        if let Some(prev) = seen.get(&id) {
            if prev.iter().any(|(c, _)| c == class) {
                continue;
            }
            let (_, k) = ["FILTER", "INFO", "FORMAT"].iter().find_map(|c| prev.iter().find(|(pc, _)| pc == c)).unwrap().clone();
            if fields[k].types == f.types && fields[k].lengths == f.lengths {
                fields[k].classes.push(class.into());
                seen.get_mut(&id).unwrap().push((class.into(), k));
                continue;
            }
            f.name = format!("{id}_{class}");
            f.vcf_name = Some(id.clone());
        }
        seen.entry(id).or_default().push((class.into(), fields.len()));
        fields.push(f);
    }
    let by = |fields: &mut Vec<VidField>, n: &str, op: &dyn Fn(&mut VidField)| {
        if let Some(f) = fields.iter_mut().find(|f| f.name == n) {
            op(f);
        }
    };
    by(&mut fields, "RAW_MQandDP", &|f| f.combine = Some("element_wise_sum".into()));
    by(&mut fields, "AS_RAW_MQ", &|f| f.make_as(&["float"], "element_wise_sum"));
    for n in ["MQ_DP", "QUALapprox", "VarDP", "TREE_SCORE", "ASSEMBLED_HAPS", "FILTERED_HAPS"] {
        by(&mut fields, n, &|f| f.combine = Some("sum".into()));
    }
    for n in ["RAW_GT_COUNT", "HAPCOMP", "HAPDOM"] {
        by(&mut fields, n, &|f| f.combine = Some("element_wise_sum".into()));
    }
    by(&mut fields, "AS_QUALapprox", &|f| f.make_as(&["int"], "element_wise_sum"));
    for n in ["AS_RAW_MQ", "AS_SB_TABLE"] {
        by(&mut fields, n, &|f| f.disable_remap = true);
    }
    (fields, contigs)
}

fn vidmap_json(fields: &[VidField], contigs: &[(String, i64, i64)]) -> String {
    let cs = contigs
        .iter()
        .map(|(n, l, o)| J::M(vec![("name", J::S(n.clone())), ("length", J::S(l.to_string())), ("tiledb_column_offset", J::S(o.to_string()))]))
        .collect();
    jprint(&J::M(vec![("fields", J::L(fields.iter().map(|f| f.json()).collect())), ("contigs", J::L(cs))]), 0)
}

fn callset_json(samples: &[(String, String)], bypass: bool) -> String {
    let cs = samples
        .iter()
        .enumerate()
        .map(|(r, (name, path))| {
            J::M(vec![
                ("sample_name", J::S(name.clone())),
                ("row_idx", J::S(r.to_string())),
                ("idx_in_file", J::S("0".into())),
                ("stream_name", J::S(if bypass { path.clone() } else { format!("{name}_stream") })),
            ])
        })
        .collect();
    jprint(&J::M(vec![("callsets", J::L(cs))]), 0)
}

/// A VCF's `##` lines (htsjdk-repaired, `##fileformat` dropped, exact duplicates once, file order) and samples.
fn read_vcf_header(path: &str) -> (Vec<String>, Vec<String>) {
    let mut lines = Lines::open(path);
    let mut line = Vec::new();
    let mut meta: Vec<String> = Vec::new();
    let mut seen = std::collections::HashSet::new();
    loop {
        if !lines.next_line(&mut line) {
            fail(&format!("{path}: no #CHROM line"));
        }
        let l = String::from_utf8_lossy(&line).into_owned();
        if !l.starts_with("##") {
            return (meta, l.split('\t').skip(9).map(|s| s.to_string()).collect());
        }
        if l.starts_with("##fileformat=") {
            continue;
        }
        let l = crate::vcfheader::repair_standard_line(l);
        if seen.insert(l.clone()) {
            meta.push(l);
        }
    }
}

/// `VCFUtils.smartMergeHeaders` over the `-V` headers (each in htsjdk's sorted order; the first line per key or
/// key-ID wins). A Number / Type conflict GATK would resolve by promoting the line is refused.
fn smart_merge(headers: &[Vec<String>]) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    let mut by_key: HashMap<String, usize> = HashMap::new();
    for h in headers {
        let mut sorted = h.clone();
        sorted.sort_by_key(|l| if l.starts_with("##contig=<") { "##contig=<".to_string() } else { l.clone() });
        for l in sorted {
            let k = crate::mergevcfs::merge_key(&l);
            match by_key.get(&k) {
                None => {
                    by_key.insert(k, out.len());
                    out.push(l);
                }
                Some(&i) => {
                    if let (Some(a), Some(b)) = (crate::vcfheader::compound(&out[i]), crate::vcfheader::compound(&l)) {
                        if (a.0 == "INFO" || a.0 == "FORMAT") && (a.2 != b.2 || a.3 != b.3) {
                            refuse(&format!("the -V headers disagree on {}: {} vs {} (GATK promotes the line); not supported", a.1, out[i], l));
                        }
                    }
                }
            }
        }
    }
    out
}

struct Args {
    workspace: String,
    sample_map: Option<String>,
    variants: Vec<String>,
    intervals: Vec<String>,
    merge_intervals: bool,
    batch_size: i64,
    consolidate: bool,
    bypass: bool,
    overwrite: bool,
    command_line: bool,
}

const USAGE: &str = "usage: variant-rs genomicsdbimport --genomicsdb-workspace-path ws (--sample-name-map map | -V x.g.vcf.gz ...) \
-L interval... [--merge-input-intervals] [--batch-size N] [--consolidate] [--bypass-feature-reader] (GATK's names)";

fn parse_args(argv: &[String]) -> Args {
    let mut a = Args {
        workspace: String::new(),
        sample_map: None,
        variants: vec![],
        intervals: vec![],
        merge_intervals: false,
        batch_size: 0,
        consolidate: false,
        bypass: false,
        overwrite: false,
        command_line: true,
    };
    let val = |i: usize| -> String { argv.get(i + 1).cloned().unwrap_or_else(|| refuse(&format!("{} needs a value\n{USAGE}", argv[i]))) };
    let flag = |i: &mut usize| -> bool {
        match argv.get(*i + 1).map(|s| s.to_ascii_lowercase()) {
            Some(v) if ["true", "false"].contains(&v.as_str()) => {
                *i += 1;
                v == "true"
            }
            _ => true,
        }
    };
    let mut i = 0;
    while i < argv.len() {
        match argv[i].as_str() {
            "--genomicsdb-workspace-path" => {
                a.workspace = val(i);
                i += 1;
            }
            "--sample-name-map" => {
                a.sample_map = Some(val(i));
                i += 1;
            }
            "-V" | "--variant" => {
                a.variants.push(val(i));
                i += 1;
            }
            "-L" | "--intervals" => {
                a.intervals.push(val(i));
                i += 1;
            }
            "--batch-size" => {
                a.batch_size = val(i).parse().unwrap_or_else(|_| refuse("--batch-size: not an integer"));
                i += 1;
            }
            "--merge-input-intervals" => a.merge_intervals = flag(&mut i),
            "--consolidate" => a.consolidate = flag(&mut i),
            "--bypass-feature-reader" => a.bypass = flag(&mut i),
            "--overwrite-existing-genomicsdb-workspace" => a.overwrite = flag(&mut i),
            "--add-output-vcf-command-line" => a.command_line = flag(&mut i),
            // no effect on the workspace's content
            "--genomicsdb-shared-posixfs-optimizations" | "--validate-sample-name-map" => {
                flag(&mut i);
            }
            "--reader-threads" | "--genomicsdb-segment-size" | "--genomicsdb-vcf-buffer-size" | "--max-num-intervals-to-import-in-parallel" | "--tmp-dir" | "-R"
            | "--reference" => i += 1,
            x => refuse(&format!("GenomicsDBImport: {x} is not supported\n{USAGE}")),
        }
        i += 1;
    }
    if a.workspace.is_empty() {
        refuse(USAGE);
    }
    if a.sample_map.is_some() == !a.variants.is_empty() {
        fail("One of --variant or --sample-name-map must be specified");
    }
    if a.intervals.is_empty() {
        fail("No intervals specified");
    }
    a
}

/// The `--sample-name-map` entries (`SampleNameMap`): sample -> file, by sample.
fn read_sample_map(path: &str, bypass: bool) -> Vec<(String, String)> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| fail(&format!("{path}: {e}")));
    let mut out: Vec<(String, String)> = Vec::new();
    for l in text.lines() {
        if l.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() != 2 && f.len() != 3 {
            fail(&format!(
                "Sample name map file must have 2 or 3 fields per line in the format:\nSample\tFile\nor:\nSample\tFile\tIndex\nbut found line: \"{l}\" with {} fields",
                f.len()
            ));
        }
        if f[0].is_empty() || f[0].trim() != f[0] {
            fail(&format!("Sample name map file must have lines in the format: Sample\tFile, but found: {l}"));
        }
        if f.len() == 3 {
            refuse("a sample name map with index files is not supported");
        }
        let file = f[1].trim().to_string();
        if let Some((_, prev)) = out.iter().find(|(s, _)| s == f[0]) {
            fail(&format!("Found two mappings for the same sample: {}\n{file}\n{prev}", f[0]));
        }
        if bypass && !file.to_ascii_lowercase().ends_with(".vcf.gz") {
            fail(&format!("Input variant files must be block compressed vcfs when using bypass-feature-reader, but {file} does not end with the standard file extension .vcf.gz"));
        }
        out.push((f[0].to_string(), file));
    }
    if out.is_empty() {
        fail("At least 1 sample is required but none were found in the sample mapping file");
    }
    out
}

fn write_text(path: &Path, text: &str) {
    std::fs::write(path, text).unwrap_or_else(|e| fail(&format!("{}: {e}", path.display())));
}

/// rapidjson's pretty printer for a flat object of integers.
fn meta_json(items: &[(&str, i64)]) -> String {
    let body: Vec<String> = items.iter().map(|(k, v)| format!("    \"{k}\": {v}")).collect();
    format!("{{\n{}\n}}", body.join(",\n"))
}

pub fn run(argv: &[String]) {
    let a = parse_args(argv);
    let t0 = std::time::Instant::now();
    // samples, rows in sorted name order (TreeMap: UTF-16 order)
    let mut samples: Vec<(String, String)> = Vec::new();
    let header_lines: Vec<String>;
    if let Some(m) = &a.sample_map {
        samples = read_sample_map(m, a.bypass);
        samples.sort_by(|x, y| x.0.encode_utf16().cmp(y.0.encode_utf16()));
        header_lines = read_vcf_header(&samples[0].1).0;
    } else {
        let mut headers = Vec::new();
        for v in &a.variants {
            let (meta, s) = read_vcf_header(v);
            if s.len() != 1 {
                fail(&format!("Input GVCF: {v} was expected to contain a single sample but actually contained {} samples.", s.len()));
            }
            if let Some((_, p)) = samples.iter().find(|(x, _)| *x == s[0]) {
                fail(&format!("Duplicate sample: {}. Sample was found in both {v} and {p}.", s[0]));
            }
            samples.push((s[0].clone(), v.clone()));
            headers.push(meta);
        }
        samples.sort_by(|x, y| x.0.encode_utf16().cmp(y.0.encode_utf16()));
        header_lines = smart_merge(&headers);
    }
    for (_, p) in &samples {
        if !(p.ends_with(".gz") || p.ends_with(".bgz")) || !Path::new(&format!("{p}.tbi")).exists() {
            refuse(&format!("{p}: GenomicsDBImport needs bgzipped, tabix-indexed GVCFs (as GATK does)"));
        }
        let (_, s) = read_vcf_header(p);
        if s.len() != 1 {
            refuse(&format!("{p}: GenomicsDBImport of a VCF with {} samples is not supported", s.len()));
        }
    }
    let (fields, contigs) = build_vidmap(&header_lines);
    if contigs.is_empty() {
        fail("The merged vcf header has no sequence dictionary. Please provide a header that contains a sequence dictionary.");
    }
    let vid_text = vidmap_json(&fields, &contigs);
    let vid = Vid::from_json(&serde_json::from_str(&vid_text).unwrap());
    let attrs = vid.attrs();

    // intervals -> arrays
    let dict = crate::intervals::Dict::new(contigs.iter().map(|c| c.0.clone()).collect(), contigs.iter().map(|c| c.1).collect());
    let mut ivs = Vec::new();
    for l in &a.intervals {
        ivs.extend(crate::intervals::parse_l_arg(l, &dict).unwrap_or_else(|e| fail(&e)));
    }
    let mut ivs = crate::intervals::sort_merge(ivs);
    if a.merge_intervals {
        let mut span: Vec<crate::intervals::Interval> = Vec::new();
        for iv in ivs {
            match span.last_mut() {
                Some(s) if s.tid == iv.tid => s.end = s.end.max(iv.end),
                _ => span.push(iv),
            }
        }
        ivs = span;
    }

    // the workspace
    let ws = Path::new(&a.workspace);
    if ws.exists() {
        if !ws.is_dir() {
            fail(&format!("Error creating GenomicsDB workspace: {} already exists and is not a directory", a.workspace));
        }
        if ws.join("__tiledb_workspace.tdb").exists() && a.overwrite {
            std::fs::remove_dir_all(ws).unwrap_or_else(|e| fail(&format!("{}: {e}", a.workspace)));
        } else if ws.join("__tiledb_workspace.tdb").exists() {
            fail(&format!("Error creating GenomicsDB workspace: {} already exists", a.workspace));
        } else {
            fail(&format!("Error creating GenomicsDB workspace: {}", a.workspace));
        }
    }
    std::fs::create_dir_all(ws).unwrap_or_else(|e| fail(&format!("Error creating GenomicsDB workspace: {}: {e}", a.workspace)));
    write_text(&ws.join("__tiledb_workspace.tdb"), "");
    write_text(&ws.join("vidmap.json"), &vid_text);
    let extra = if a.command_line { vec!["##source=GenomicsDBImport".to_string(), crate::vcfheader::command_line("GenomicsDBImport")] } else { vec![] };
    write_text(&ws.join("vcfheader.vcf"), &crate::vcfheader::write_header(&header_lines, extra, &[]));

    // batches of rows; with --consolidate (or one batch) one fragment and one meta file per array
    let n = samples.len() as i64;
    let b = if a.batch_size <= 0 { n } else { a.batch_size.min(n) };
    let batches: Vec<(i64, i64)> = if a.consolidate { vec![(0, n)] } else { (0..n).step_by(b as usize).map(|lb| (lb, (lb + b).min(n))).collect() };
    let mut total = 0u64;
    for (bi, &(lb, ub)) in batches.iter().enumerate() {
        let rows: Vec<(i64, String)> = (lb..ub).map(|r| (r, samples[r as usize].1.clone())).collect();
        for iv in &ivs {
            let contig = &dict.names[iv.tid];
            let off = contigs[iv.tid].2;
            let d = ws.join(format!("{contig}${}${}", iv.start, iv.end));
            if bi == 0 {
                std::fs::create_dir_all(d.join("genomicsdb_meta_dir")).unwrap_or_else(|e| fail(&format!("{}: {e}", d.display())));
                gdb_write::write_schema(&d.join("__array_schema.tdb"), &format!("{}/{}", a.workspace, d.file_name().unwrap().to_string_lossy()), &attrs)
                    .unwrap_or_else(|e| fail(&format!("{}: {e}", d.display())));
                write_text(&d.join(".__consolidation_lock"), "");
                write_text(&d.join("genomicsdb_meta_dir/genomicsdb_column_bounds.json"), &meta_json(&[("min_column", off + iv.start - 1), ("max_column", off + iv.end - 1)]));
            }
            // the running row bounds grow with every batch: one meta file each
            write_text(
                &d.join(format!("genomicsdb_meta_dir/genomicsdb_meta_{}.json", gdb_write::uuid())),
                &meta_json(&[("lb_row_idx", 0), ("max_valid_row_idx_in_array", ub - 1)]),
            );
            total += write_fragment(&d, &vid, &attrs, &rows, contig, iv.start, iv.end, !a.bypass);
            // fragments are ordered by their millisecond timestamp
            std::thread::sleep(std::time::Duration::from_millis(2));
        }
    }
    write_text(&ws.join("callset.json"), &callset_json(&samples, a.bypass));
    eprintln!("[genomicsdbimport] {} samples, {} arrays, {} batches, {total} cells, in {:.1}s", samples.len(), ivs.len(), batches.len(), t0.elapsed().as_secs_f64());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn htslib_types() {
        assert!(is_indel_or_overlap("GAAAA", "G"));
        assert!(is_indel_or_overlap("A", "*"));
        assert!(!is_indel_or_overlap("A", "<NON_REF>"));
        assert!(!is_indel_or_overlap("A", "G"));
        assert!(!is_indel_or_overlap("ATG", "CT"));
        assert_eq!(strtol0("0x10"), 16);
        assert_eq!(strtol0("010"), 8);
        assert_eq!(strtol0("-12"), -12);
    }
}
