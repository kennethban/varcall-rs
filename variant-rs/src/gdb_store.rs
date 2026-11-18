//! GenomicsDB 1.5.5 workspaces (as GATK's GenomicsDBImport writes them): the storage layer. Reads the TileDB fork's
//! sparse arrays cell by cell (`docs/specs/genomicsdb_tiledb.md`): `__array_schema.tdb`, live fragments (marker file,
//! ordered by timestamp), `__book_keeping.tdb.gz`, zlib tiles (+ delta-encoded coordinates and var offsets), cells in
//! (column, row) order merged over fragments (the newest fragment wins on equal coordinates).

use flate2::read::{GzDecoder, ZlibDecoder};
use std::fs::File;
use std::io::Read;
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};

pub const TILEDB_VAR_NUM: u32 = 0x7FFF_FFFF;

/// One attribute of the array schema.
#[derive(Clone, Debug)]
pub struct Attr {
    pub name: String,
    /// 0 INT32, 1 INT64, 2 FLOAT32, 4 CHAR (others not produced by GenomicsDB)
    pub typ: u8,
    pub cell_val_num: u32,
    pub compression: u8,
    pub offsets_compression: u8,
}

impl Attr {
    pub fn is_var(&self) -> bool {
        self.cell_val_num == TILEDB_VAR_NUM
    }
    pub fn type_size(&self) -> usize {
        match self.typ {
            0 | 2 | 9 => 4,
            1 | 3 | 10 => 8,
            4 | 5 | 6 => 1,
            7 | 8 => 2,
            t => panic!("TileDB type {t}"),
        }
    }
}

#[derive(Debug)]
pub struct Schema {
    pub attrs: Vec<Attr>,
    pub capacity: i64,
    pub coords_compression: u8,
}

struct Rd<'a> {
    b: &'a [u8],
    o: usize,
}

impl<'a> Rd<'a> {
    fn take(&mut self, n: usize) -> &'a [u8] {
        let s = &self.b[self.o..self.o + n];
        self.o += n;
        s
    }
    fn u8(&mut self) -> u8 {
        self.take(1)[0]
    }
    fn i32(&mut self) -> i32 {
        i32::from_le_bytes(self.take(4).try_into().unwrap())
    }
    fn i64(&mut self) -> i64 {
        i64::from_le_bytes(self.take(8).try_into().unwrap())
    }
    fn str(&mut self) -> String {
        let n = self.i32() as usize;
        String::from_utf8_lossy(self.take(n)).into_owned()
    }
}

/// `ArraySchema::deserialize` (schema version 2)
pub fn read_schema(path: &Path) -> Schema {
    let b = std::fs::read(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut r = Rd { b: &b, o: 0 };
    let tag = r.i32() as u32;
    assert_eq!(tag, 0xFFFF_FFFD, "{}: unsupported TileDB schema version", path.display());
    let _workspace = r.str();
    let _name = r.str();
    let dense = r.u8();
    assert_eq!(dense, 0, "dense TileDB arrays are not GenomicsDB arrays");
    let _tile_order = r.u8();
    let _cell_order = r.u8();
    let capacity = r.i64();
    let n = r.i32() as usize;
    let names: Vec<String> = (0..n).map(|_| r.str()).collect();
    let d = r.i32() as usize;
    assert_eq!(d, 2, "GenomicsDB arrays have 2 dimensions");
    for _ in 0..d {
        let _ = r.str();
    }
    let dom = r.i32() as usize;
    r.take(dom);
    let te = r.i32() as usize;
    r.take(te);
    let types: Vec<u8> = (0..=n).map(|_| r.u8()).collect();
    let cvn: Vec<u32> = (0..n).map(|_| r.i32() as u32).collect();
    let comp: Vec<u8> = (0..=n).map(|_| r.u8()).collect();
    let _lvl: Vec<u8> = (0..=n).map(|_| r.u8()).collect();
    let ocomp: Vec<u8> = (0..n).map(|_| r.u8()).collect();
    let _olvl: Vec<u8> = (0..n).map(|_| r.u8()).collect();
    assert_eq!(r.o, b.len(), "{}: trailing bytes", path.display());
    assert_eq!(types[n], 1, "coordinates must be INT64");
    let attrs = (0..n)
        .map(|i| Attr { name: names[i].clone(), typ: types[i], cell_val_num: cvn[i], compression: comp[i], offsets_compression: ocomp[i] })
        .collect();
    Schema { attrs, capacity, coords_compression: comp[n] }
}

/// `__book_keeping.tdb.gz`
struct BookKeeping {
    /// per attribute (coordinates last): byte offset of each tile in `<attr>.tdb`
    tile_offsets: Vec<Vec<i64>>,
    tile_var_offsets: Vec<Vec<i64>>,
    tile_var_sizes: Vec<Vec<u64>>,
    /// per tile: (col_lo, col_hi)
    mbr_cols: Vec<(i64, i64)>,
    last_tile_cell_num: i64,
}

fn read_book_keeping(path: &Path, n_attrs: usize) -> BookKeeping {
    let mut raw = Vec::new();
    GzDecoder::new(File::open(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))).read_to_end(&mut raw).unwrap();
    let mut r = Rd { b: &raw, o: 0 };
    let ned = r.i64() as usize;
    r.take(ned);
    let mbr_num = r.i64() as usize;
    let mut mbr_cols = Vec::with_capacity(mbr_num);
    for _ in 0..mbr_num {
        let (_rl, _rh, cl, ch) = (r.i64(), r.i64(), r.i64(), r.i64());
        mbr_cols.push((cl, ch));
    }
    let bc = r.i64() as usize;
    r.take(bc * 32);
    let tile_offsets = (0..=n_attrs).map(|_| {
        let k = r.i64() as usize;
        (0..k).map(|_| r.i64()).collect()
    }).collect();
    let tile_var_offsets = (0..n_attrs).map(|_| {
        let k = r.i64() as usize;
        (0..k).map(|_| r.i64()).collect()
    }).collect();
    let tile_var_sizes = (0..n_attrs).map(|_| {
        let k = r.i64() as usize;
        (0..k).map(|_| r.i64() as u64).collect()
    }).collect();
    let last_tile_cell_num = r.i64();
    BookKeeping { tile_offsets, tile_var_offsets, tile_var_sizes, mbr_cols, last_tile_cell_num }
}

fn inflate(comp: &[u8], expected: usize, what: &str) -> Vec<u8> {
    let mut out = Vec::with_capacity(expected);
    ZlibDecoder::new(comp).read_to_end(&mut out).unwrap_or_else(|e| panic!("{what}: zlib: {e}"));
    assert_eq!(out.len(), expected, "{what}: inflated {} bytes, expected {expected}", out.len());
    out
}

/// A file read slice by slice (tiles), with its size.
struct TFile {
    f: File,
    size: u64,
}

impl TFile {
    fn open(p: &Path) -> TFile {
        let f = File::open(p).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
        let size = f.metadata().unwrap().len();
        TFile { f, size }
    }
    fn slice(&self, from: u64, to: u64) -> Vec<u8> {
        let mut b = vec![0u8; (to - from) as usize];
        self.f.read_exact_at(&mut b, from).unwrap();
        b
    }
}

/// One live fragment.
pub struct Fragment {
    dir: PathBuf,
    bk: BookKeeping,
    files: Vec<TFile>,
    var_files: Vec<Option<TFile>>,
    coords: TFile,
    capacity: i64,
}

/// The decoded cells of one tile: coordinates and, per attribute, the raw value bytes of every cell.
pub struct Tile {
    pub coords: Vec<(i64, i64)>,
    /// per attribute: per cell (start, end) into `data[attr]`
    pub ranges: Vec<Vec<(usize, usize)>>,
    pub data: Vec<Vec<u8>>,
}

impl Fragment {
    fn n_tiles(&self) -> usize {
        self.bk.mbr_cols.len()
    }
    fn tile_cells(&self, t: usize) -> usize {
        if t + 1 == self.n_tiles() { self.bk.last_tile_cell_num as usize } else { self.capacity as usize }
    }
    fn tile_bytes(&self, offs: &[i64], t: usize, file: &TFile) -> Vec<u8> {
        let from = offs[t] as u64;
        let to = if t + 1 < offs.len() { offs[t + 1] as u64 } else { file.size };
        file.slice(from, to)
    }

    /// Decodes tile `t` (all attributes).
    pub fn read_tile(&self, schema: &Schema, t: usize) -> Tile {
        let n = self.tile_cells(t);
        let na = schema.attrs.len();
        // coordinates: i64 (row, column) pairs, delta-encoded with stride 2
        let raw = self.tile_bytes(&self.bk.tile_offsets[na], t, &self.coords);
        let mut v: Vec<i64> = inflate(&raw, n * 16, "coords").chunks_exact(8).map(|c| i64::from_le_bytes(c.try_into().unwrap())).collect();
        if schema.coords_compression & 0x10 != 0 {
            for i in 2..v.len() {
                v[i] = v[i].wrapping_add(v[i - 2]);
            }
        }
        let coords: Vec<(i64, i64)> = v.chunks_exact(2).map(|c| (c[0], c[1])).collect();
        let mut ranges = Vec::with_capacity(na);
        let mut data = Vec::with_capacity(na);
        for (a, at) in schema.attrs.iter().enumerate() {
            assert_eq!(at.compression & 0x0f, 1, "{}: TileDB compression {:#x} not supported", at.name, at.compression);
            if !at.is_var() {
                assert_eq!(at.compression & 0x10, 0, "{}: delta-encoded fixed attribute", at.name);
                let sz = at.cell_val_num as usize * at.type_size();
                let raw = self.tile_bytes(&self.bk.tile_offsets[a], t, &self.files[a]);
                let d = inflate(&raw, n * sz, &at.name);
                ranges.push((0..n).map(|i| (i * sz, (i + 1) * sz)).collect());
                data.push(d);
            } else {
                assert_eq!(at.offsets_compression & 0x0f, 1, "{}: offsets compression", at.name);
                let raw = self.tile_bytes(&self.bk.tile_offsets[a], t, &self.files[a]);
                let mut o: Vec<u64> = inflate(&raw, n * 8, &at.name).chunks_exact(8).map(|c| u64::from_le_bytes(c.try_into().unwrap())).collect();
                if at.offsets_compression & 0x10 != 0 {
                    for i in 1..o.len() {
                        o[i] = o[i].wrapping_add(o[i - 1]);
                    }
                }
                let vsize = self.bk.tile_var_sizes[a][t] as usize;
                let d = if vsize == 0 {
                    Vec::new()
                } else {
                    let vf = self.var_files[a].as_ref().unwrap();
                    let raw = self.tile_bytes(&self.bk.tile_var_offsets[a], t, vf);
                    inflate(&raw, vsize, &format!("{}_var", at.name))
                };
                let base = o.first().copied().unwrap_or(0);
                let r: Vec<(usize, usize)> = (0..n)
                    .map(|i| {
                        let s = (o[i] - base) as usize;
                        let e = if i + 1 < n { (o[i + 1] - base) as usize } else { vsize };
                        (s, e)
                    })
                    .collect();
                ranges.push(r);
                data.push(d);
            }
        }
        Tile { coords, ranges, data }
    }
}

/// An array (one imported interval) with its live fragments, oldest first.
pub struct Array {
    #[allow(dead_code)]
    pub name: String,
    pub schema: Schema,
    pub fragments: Vec<Fragment>,
}

pub fn open_array(dir: &Path) -> Array {
    let schema = read_schema(&dir.join("__array_schema.tdb"));
    let mut frags: Vec<(i64, PathBuf)> = Vec::new();
    for e in std::fs::read_dir(dir).unwrap() {
        let p = e.unwrap().path();
        let name = p.file_name().unwrap().to_string_lossy().into_owned();
        if !name.starts_with("__") || !p.is_dir() || !p.join("__tiledb_fragment.tdb").exists() {
            continue;
        }
        // __<uuid><tid>_<ms>: the timestamp after the first '_' following the leading "__"
        let ms: i64 = name[2..].split_once('_').map(|(_, t)| t.parse().unwrap()).unwrap_or_else(|| panic!("fragment name {name}"));
        frags.push((ms, p));
    }
    frags.sort();
    let na = schema.attrs.len();
    let fragments = frags
        .into_iter()
        .map(|(_, d)| {
            let bk = read_book_keeping(&d.join("__book_keeping.tdb.gz"), na);
            let files = schema.attrs.iter().map(|a| TFile::open(&d.join(format!("{}.tdb", a.name)))).collect();
            let var_files = schema.attrs.iter().map(|a| a.is_var().then(|| TFile::open(&d.join(format!("{}_var.tdb", a.name))))).collect();
            let coords = TFile::open(&d.join("__coords.tdb"));
            Fragment { dir: d, bk, files, var_files, coords, capacity: schema.capacity }
        })
        .collect();
    Array { name: dir.file_name().unwrap().to_string_lossy().into_owned(), schema, fragments }
}

/// A cell: coordinates and the value bytes of every attribute (schema order).
#[derive(Clone, Debug)]
pub struct Cell {
    pub row: i64,
    pub col: i64,
    pub vals: Vec<Vec<u8>>,
}

/// Streams the cells of one fragment in (column, row) order, tiles overlapping [c0, c1] only.
struct FragCursor<'a> {
    frag: &'a Fragment,
    schema: &'a Schema,
    t: usize,
    tile: Option<Tile>,
    i: usize,
    c0: i64,
    c1: i64,
}

impl<'a> FragCursor<'a> {
    fn peek_key(&mut self) -> Option<(i64, i64)> {
        loop {
            if let Some(tile) = &self.tile {
                while self.i < tile.coords.len() && tile.coords[self.i].1 < self.c0 {
                    self.i += 1;
                }
                if self.i < tile.coords.len() {
                    let (r, c) = tile.coords[self.i];
                    return if c > self.c1 { None } else { Some((c, r)) };
                }
                self.tile = None;
                self.t += 1;
            }
            while self.t < self.frag.n_tiles() && self.frag.bk.mbr_cols[self.t].1 < self.c0 {
                self.t += 1;
            }
            if self.t >= self.frag.n_tiles() || self.frag.bk.mbr_cols[self.t].0 > self.c1 {
                return None;
            }
            self.tile = Some(self.frag.read_tile(self.schema, self.t));
            self.i = 0;
        }
    }
    fn take(&mut self) -> Cell {
        let tile = self.tile.as_ref().unwrap();
        let (row, col) = tile.coords[self.i];
        let vals = (0..tile.data.len()).map(|a| {
            let (s, e) = tile.ranges[a][self.i];
            tile.data[a][s..e].to_vec()
        }).collect();
        self.i += 1;
        Cell { row, col, vals }
    }
}

/// The cells of an array with column in [c0, c1], in (column, row) order, fragments merged (the newest fragment
/// wins on equal coordinates).
pub struct CellStream<'a> {
    cur: Vec<FragCursor<'a>>,
}

impl<'a> CellStream<'a> {
    pub fn new(arr: &'a Array, c0: i64, c1: i64) -> CellStream<'a> {
        CellStream { cur: arr.fragments.iter().map(|frag| FragCursor { frag, schema: &arr.schema, t: 0, tile: None, i: 0, c0, c1 }).collect() }
    }
}

impl Iterator for CellStream<'_> {
    type Item = Cell;
    fn next(&mut self) -> Option<Cell> {
        let mut best: Option<((i64, i64), usize)> = None;
        for k in 0..self.cur.len() {
            if let Some(key) = self.cur[k].peek_key() {
                // equal keys: the later (newer) fragment wins
                if best.is_none_or(|(b, _)| key <= b) {
                    best = Some((key, k));
                }
            }
        }
        let (key, k) = best?;
        // drop the same coordinates in older fragments
        for j in 0..k {
            if self.cur[j].peek_key() == Some(key) {
                self.cur[j].take();
            }
        }
        Some(self.cur[k].take())
    }
}

/// All cells with column in [c0, c1] of the array (see `CellStream`).
pub fn for_each_cell(arr: &Array, c0: i64, c1: i64, f: impl FnMut(Cell)) {
    CellStream::new(arr, c0, c1).for_each(f);
}

#[allow(dead_code)]
impl Fragment {
    pub fn dir(&self) -> &Path {
        &self.dir
    }
}

// ------------------------------------------------------------------------------------------------
// workspace metadata and value decoding (docs/specs/genomicsdb_tiledb.md sections 1 and 7)

/// `callset.json`, `vidmap.json` and `vcfheader.vcf` of a workspace.
pub struct Workspace {
    pub dir: PathBuf,
    /// sample names by row
    pub samples: Vec<String>,
    /// (contig, length, tiledb column offset), in vidmap order
    pub contigs: Vec<(String, i64, i64)>,
    pub header: String,
    /// array directories (one per imported interval)
    pub arrays: Vec<PathBuf>,
}

fn json_i64(v: &serde_json::Value) -> i64 {
    match v {
        serde_json::Value::String(s) => s.parse().unwrap(),
        serde_json::Value::Number(n) => n.as_i64().unwrap(),
        x => panic!("not an integer: {x}"),
    }
}

impl Workspace {
    pub fn open(dir: &Path) -> Workspace {
        let rd = |f: &str| std::fs::read_to_string(dir.join(f)).unwrap_or_else(|e| panic!("{}/{f}: {e}", dir.display()));
        let cs: serde_json::Value = serde_json::from_str(&rd("callset.json")).expect("callset.json");
        let mut samples: Vec<(i64, String)> =
            cs["callsets"].as_array().unwrap().iter().map(|c| (json_i64(&c["row_idx"]), c["sample_name"].as_str().unwrap().to_string())).collect();
        samples.sort();
        for (k, (r, _)) in samples.iter().enumerate() {
            assert_eq!(*r, k as i64, "callset rows are not 0..n");
        }
        let vid: serde_json::Value = serde_json::from_str(&rd("vidmap.json")).expect("vidmap.json");
        let contigs = vid["contigs"]
            .as_array()
            .unwrap()
            .iter()
            .map(|c| (c["name"].as_str().unwrap().to_string(), json_i64(&c["length"]), json_i64(&c["tiledb_column_offset"])))
            .collect();
        let mut arrays: Vec<PathBuf> = std::fs::read_dir(dir).unwrap().map(|e| e.unwrap().path()).filter(|p| p.join("__array_schema.tdb").exists()).collect();
        arrays.sort();
        Workspace { dir: dir.to_path_buf(), samples: samples.into_iter().map(|(_, s)| s).collect(), contigs, header: rd("vcfheader.vcf"), arrays }
    }
    #[allow(dead_code)]
    /// (contig index, 1-based position) of a global column
    pub fn locate(&self, col: i64) -> (usize, i64) {
        let k = self.contigs.iter().rposition(|c| c.2 <= col).expect("column before the first contig");
        (k, col - self.contigs[k].2 + 1)
    }
    #[allow(dead_code)]
    pub fn column(&self, contig: &str, pos: i64) -> Option<i64> {
        self.contigs.iter().find(|c| c.0 == contig).map(|c| c.2 + pos - 1)
    }
}

pub const INT_NULL: i32 = 0x7FFF_FFFF;
pub const INT_MISSING: i32 = 0x8000_0000u32 as i32;
pub const FLT_NULL_BITS: u32 = 0x7F7F_FFFF;
pub const FLT_MISSING_BITS: u32 = 0x7F80_0001;

/// A decoded numeric element: a value or BCF missing ('.').
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Num {
    I(i32),
    F(f32),
    Missing,
}

/// The attribute values of a cell, by attribute name.
pub struct CellView<'a> {
    pub schema: &'a Schema,
    pub cell: &'a Cell,
}

impl<'a> CellView<'a> {
    fn raw(&self, name: &str) -> Option<&'a [u8]> {
        let i = self.schema.attrs.iter().position(|a| a.name == name)?;
        Some(&self.cell.vals[i])
    }
    pub fn end(&self) -> i64 {
        i64::from_le_bytes(self.raw("END").unwrap()[..8].try_into().unwrap())
    }
    pub fn text(&self, name: &str) -> Option<String> {
        self.raw(name).filter(|b| !b.is_empty()).map(|b| String::from_utf8_lossy(b).into_owned())
    }
    pub fn ref_allele(&self) -> String {
        self.text("REF").unwrap_or_default()
    }
    /// ALT alleles, `&` -> `<NON_REF>`
    pub fn alts(&self) -> Vec<String> {
        match self.text("ALT") {
            None => Vec::new(),
            Some(s) => s.split('|').map(|a| if a == "&" { "<NON_REF>".to_string() } else { a.to_string() }).collect(),
        }
    }
    #[allow(dead_code)]
    /// QUAL (None: missing)
    pub fn qual(&self) -> Option<f32> {
        let b = self.raw("QUAL")?;
        let bits = u32::from_le_bytes(b[..4].try_into().unwrap());
        (bits != FLT_NULL_BITS).then(|| f32::from_bits(bits))
    }
    /// A 1-D numeric field: None if the whole field is missing; elements BCF-missing as `Num::Missing`.
    pub fn nums(&self, name: &str) -> Option<Vec<Num>> {
        let i = self.schema.attrs.iter().position(|a| a.name == name)?;
        let at = &self.schema.attrs[i];
        let b = &self.cell.vals[i];
        if b.is_empty() {
            return None;
        }
        let v: Vec<Num> = match at.typ {
            0 => b.chunks_exact(4).map(|c| i32::from_le_bytes(c.try_into().unwrap())).map(|x| if x == INT_MISSING { Num::Missing } else { Num::I(x) }).collect(),
            2 => b.chunks_exact(4).map(|c| u32::from_le_bytes(c.try_into().unwrap())).map(|x| if x == FLT_MISSING_BITS { Num::Missing } else { Num::F(f32::from_bits(x)) }).collect(),
            t => panic!("{name}: TileDB type {t} is not numeric"),
        };
        // fixed-length fields: all-TileDB-null = missing field
        if !at.is_var() {
            let null = |n: &Num| match n {
                Num::I(x) => *x == INT_NULL,
                Num::F(x) => x.to_bits() == FLT_NULL_BITS,
                Num::Missing => false,
            };
            if v.iter().all(null) {
                return None;
            }
        }
        Some(v)
    }
    #[allow(dead_code)]
    /// GT: allele indices (-1 = '.') and phasing (separator before allele k is '|')
    pub fn gt(&self) -> Option<(Vec<i32>, bool)> {
        let v = self.nums("GT")?;
        let x: Vec<i32> = v.iter().map(|n| match n { Num::I(x) => *x, _ => -1 }).collect();
        let alleles: Vec<i32> = x.iter().step_by(2).copied().collect();
        let phased = x.iter().skip(1).step_by(2).any(|&p| p == 1);
        Some((alleles, phased))
    }
    /// A multi-dimensional (allele-specific) field blob: per outer entry (allele) its inner vector.
    pub fn multid(&self, name: &str, float: bool) -> Option<Vec<Vec<Num>>> {
        let b = self.raw(name)?;
        if b.is_empty() {
            return None;
        }
        let rd = |o: usize| u64::from_le_bytes(b[o..o + 8].try_into().unwrap()) as usize;
        let s = rd(0);
        let data = &b[8..8 + s];
        let n = rd(8 + s);
        let offs: Vec<usize> = (0..=n).map(|k| rd(16 + s + 8 * k)).collect();
        let el = |c: &[u8]| {
            let x = u32::from_le_bytes(c.try_into().unwrap());
            if float {
                if x == FLT_MISSING_BITS { Num::Missing } else { Num::F(f32::from_bits(x)) }
            } else if x as i32 == INT_MISSING {
                Num::Missing
            } else {
                Num::I(x as i32)
            }
        };
        Some((0..n).map(|k| data[offs[k]..offs[k + 1]].chunks_exact(4).map(el).collect()).collect())
    }
}
