//! GenomicsDB 1.5.5 workspaces, the storage layer for writing (the inverse of `gdb_store.rs`; format in
//! `docs/specs/genomicsdb_tiledb.md` sections 2-4): the array schema and fragments as GenomicsDBImport's TileDB fork
//! writes them. A fragment takes its cells in (column, row) order and stores them in tiles of `CAPACITY` cells:
//! per attribute a zlib-compressed tile of values (`<attr>.tdb`), or for a variable-length attribute a tile of
//! absolute u64 offsets, delta-encoded (`<attr>.tdb`), and a tile of the values (`<attr>_var.tdb`); the
//! coordinates (row, column) delta-encoded with stride 2 (`__coords.tdb`); then the book-keeping (gzip) and the
//! fragment marker. The compressed bytes differ from TileDB's (another deflate implementation); the cells do not.

use crate::gdb_store::{Attr, TILEDB_VAR_NUM};
use flate2::write::{GzEncoder, ZlibEncoder};
use flate2::Compression;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

/// Cells per tile (GenomicsDBImport's default; GATK does not change it).
pub const CAPACITY: usize = 1000;
const GZIP: u8 = 0x01;
const GZIP_DELTA: u8 = 0x11;

/// The attribute list's types: TileDB type codes as in `gdb_store::Attr`.
pub const INT32: u8 = 0;
pub const INT64: u8 = 1;
pub const FLOAT32: u8 = 2;
pub const CHAR: u8 = 4;

pub fn attr(name: &str, typ: u8, cell_val_num: u32) -> Attr {
    let var = cell_val_num == TILEDB_VAR_NUM;
    Attr { name: name.to_string(), typ, cell_val_num, compression: GZIP, offsets_compression: if var { GZIP_DELTA } else { 0 } }
}

fn put_str(b: &mut Vec<u8>, s: &str) {
    b.extend_from_slice(&(s.len() as i32).to_le_bytes());
    b.extend_from_slice(s.as_bytes());
}

/// `__array_schema.tdb` (version 2): sparse, column-major cells, dimensions "samples" and "position" over
/// [0, 2^63 - 2], no tile extents; `name` is the informational array name (`<workspace>/<array>`).
pub fn write_schema(path: &Path, name: &str, attrs: &[Attr]) -> std::io::Result<()> {
    let mut b = Vec::new();
    b.extend_from_slice(&0xFFFF_FFFDu32.to_le_bytes());
    put_str(&mut b, "");
    put_str(&mut b, name);
    b.extend_from_slice(&[0, 0, 1]); // sparse, tile order row-major, cell order column-major
    b.extend_from_slice(&(CAPACITY as i64).to_le_bytes());
    b.extend_from_slice(&(attrs.len() as i32).to_le_bytes());
    for a in attrs {
        put_str(&mut b, &a.name);
    }
    b.extend_from_slice(&2i32.to_le_bytes());
    put_str(&mut b, "samples");
    put_str(&mut b, "position");
    b.extend_from_slice(&32i32.to_le_bytes());
    for v in [0, i64::MAX - 1, 0, i64::MAX - 1] {
        b.extend_from_slice(&v.to_le_bytes());
    }
    b.extend_from_slice(&0i32.to_le_bytes()); // no tile extents
    b.extend(attrs.iter().map(|a| a.typ));
    b.push(INT64); // coordinates
    for a in attrs {
        b.extend_from_slice(&(a.cell_val_num as i32).to_le_bytes());
    }
    b.extend(attrs.iter().map(|a| a.compression));
    b.push(GZIP_DELTA);
    b.extend(std::iter::repeat(0xFFu8).take(attrs.len() + 1)); // levels: Z_DEFAULT_COMPRESSION (-1)
    b.extend(attrs.iter().map(|a| a.offsets_compression));
    b.extend(attrs.iter().map(|a| if a.is_var() { 0xFFu8 } else { 0 }));
    std::fs::write(path, b)
}

fn zlib(data: &[u8]) -> Vec<u8> {
    let mut e = ZlibEncoder::new(Vec::with_capacity(data.len() / 4 + 64), Compression::default());
    e.write_all(data).unwrap();
    e.finish().unwrap()
}

/// A file being written tile by tile, with the offset of each tile.
struct TileFile {
    w: BufWriter<File>,
    pos: i64,
}

impl TileFile {
    fn create(p: &Path) -> std::io::Result<TileFile> {
        Ok(TileFile { w: BufWriter::with_capacity(1 << 20, File::create(p)?), pos: 0 })
    }
    /// Writes a compressed tile; returns its offset.
    fn put(&mut self, raw: &[u8]) -> std::io::Result<i64> {
        let at = self.pos;
        let z = zlib(raw);
        self.w.write_all(&z)?;
        self.pos += z.len() as i64;
        Ok(at)
    }
}

/// One fragment being written. Cells must come in (column, row) order, at most one per coordinate.
pub struct FragmentWriter {
    tmp: PathBuf,
    done: PathBuf,
    attrs: Vec<Attr>,
    files: Vec<TileFile>,
    var_files: Vec<Option<TileFile>>,
    coords: TileFile,
    /// the tile being filled: coordinates and per attribute the cells' bytes
    tile_coords: Vec<(i64, i64)>,
    tile_vals: Vec<Vec<Vec<u8>>>,
    /// per variable-length attribute: bytes written so far to its value stream (offsets are absolute)
    var_total: Vec<u64>,
    // book-keeping
    mbrs: Vec<[i64; 4]>,
    bounds: Vec<[i64; 4]>,
    tile_offsets: Vec<Vec<i64>>,
    tile_var_offsets: Vec<Vec<i64>>,
    tile_var_sizes: Vec<Vec<u64>>,
    last: Option<(i64, i64)>,
    cells: u64,
}

/// A random (version 4) UUID, lower case.
pub fn uuid() -> String {
    let mut r = [0u8; 16];
    File::open("/dev/urandom").and_then(|mut f| std::io::Read::read_exact(&mut f, &mut r)).expect("/dev/urandom");
    r[6] = (r[6] & 0x0f) | 0x40;
    r[8] = (r[8] & 0x3f) | 0x80;
    let h: Vec<String> = r.iter().map(|b| format!("{b:02x}")).collect();
    format!("{}-{}-{}-{}-{}", h[0..4].concat(), h[4..6].concat(), h[6..8].concat(), h[8..10].concat(), h[10..16].concat())
}

/// A fragment directory name: `__<uuid><thread id>_<ms since the epoch>`.
pub fn fragment_name() -> String {
    let uuid = uuid();
    let ms = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_millis();
    format!("__{uuid}{}_{ms}", std::process::id() as u64 * 1_000_003)
}

impl FragmentWriter {
    /// Starts a fragment in `array_dir` (written as `.__name`, renamed by `finish`).
    pub fn create(array_dir: &Path, attrs: &[Attr]) -> std::io::Result<FragmentWriter> {
        let name = fragment_name();
        let tmp = array_dir.join(format!(".{name}"));
        let done = array_dir.join(&name);
        std::fs::create_dir_all(&tmp)?;
        let mut files = Vec::new();
        let mut var_files = Vec::new();
        for a in attrs {
            files.push(TileFile::create(&tmp.join(format!("{}.tdb", a.name)))?);
            var_files.push(if a.is_var() { Some(TileFile::create(&tmp.join(format!("{}_var.tdb", a.name)))?) } else { None });
        }
        let n = attrs.len();
        Ok(FragmentWriter {
            coords: TileFile::create(&tmp.join("__coords.tdb"))?,
            tmp,
            done,
            attrs: attrs.to_vec(),
            files,
            var_files,
            tile_coords: Vec::with_capacity(CAPACITY),
            tile_vals: vec![Vec::with_capacity(CAPACITY); n],
            var_total: vec![0; n],
            mbrs: vec![],
            bounds: vec![],
            tile_offsets: vec![vec![]; n + 1],
            tile_var_offsets: vec![vec![]; n],
            tile_var_sizes: vec![vec![]; n],
            last: None,
            cells: 0,
        })
    }

    /// Adds a cell: its coordinates and every attribute's bytes (schema order; a fixed attribute exactly
    /// cell_val_num values).
    pub fn push(&mut self, row: i64, col: i64, vals: Vec<Vec<u8>>) -> std::io::Result<()> {
        if let Some((r0, c0)) = self.last {
            assert!((col, row) > (c0, r0), "cells out of (column, row) order: ({row},{col}) after ({r0},{c0})");
        }
        self.last = Some((row, col));
        debug_assert_eq!(vals.len(), self.attrs.len());
        for (a, v) in vals.into_iter().enumerate() {
            debug_assert!(self.attrs[a].is_var() || v.len() == self.attrs[a].cell_val_num as usize * self.attrs[a].type_size());
            self.tile_vals[a].push(v);
        }
        self.tile_coords.push((row, col));
        self.cells += 1;
        if self.tile_coords.len() == CAPACITY {
            self.flush_tile()?;
        }
        Ok(())
    }

    fn flush_tile(&mut self) -> std::io::Result<()> {
        let n = self.tile_coords.len();
        if n == 0 {
            return Ok(());
        }
        let c = &self.tile_coords;
        let (rl, rh) = (c.iter().map(|x| x.0).min().unwrap(), c.iter().map(|x| x.0).max().unwrap());
        self.mbrs.push([rl, rh, c[0].1, c[n - 1].1]);
        self.bounds.push([c[0].0, c[0].1, c[n - 1].0, c[n - 1].1]);
        for (a, at) in self.attrs.iter().enumerate() {
            let vals = std::mem::take(&mut self.tile_vals[a]);
            if !at.is_var() {
                let raw = vals.concat();
                self.tile_offsets[a].push(self.files[a].put(&raw)?);
            } else {
                // absolute offsets into the fragment's value stream, delta-encoded within the tile
                let mut offs = Vec::with_capacity(n * 8);
                let mut prev = 0u64;
                for (i, v) in vals.iter().enumerate() {
                    let o = self.var_total[a];
                    offs.extend_from_slice(&(if i == 0 { o } else { o - prev }).to_le_bytes());
                    prev = o;
                    self.var_total[a] += v.len() as u64;
                }
                self.tile_offsets[a].push(self.files[a].put(&offs)?);
                let data = vals.concat();
                let vf = self.var_files[a].as_mut().unwrap();
                let at = if data.is_empty() { vf.pos } else { vf.put(&data)? };
                self.tile_var_offsets[a].push(at);
                self.tile_var_sizes[a].push(data.len() as u64);
            }
            self.tile_vals[a] = Vec::with_capacity(CAPACITY);
        }
        let mut raw = Vec::with_capacity(n * 16);
        for (i, &(r, col)) in c.iter().enumerate() {
            let (dr, dc) = if i == 0 { (r, col) } else { (r.wrapping_sub(c[i - 1].0), col.wrapping_sub(c[i - 1].1)) };
            raw.extend_from_slice(&dr.to_le_bytes());
            raw.extend_from_slice(&dc.to_le_bytes());
        }
        let na = self.attrs.len();
        self.tile_offsets[na].push(self.coords.put(&raw)?);
        self.tile_coords.clear();
        Ok(())
    }

    /// Writes the last tile and the book-keeping, and makes the fragment live. Returns the number of cells.
    pub fn finish(mut self) -> std::io::Result<u64> {
        let last_n = match self.tile_coords.len() {
            0 => CAPACITY as i64,
            n => n as i64,
        };
        self.flush_tile()?;
        for f in self.files.iter_mut().chain(self.var_files.iter_mut().flatten()).chain(std::iter::once(&mut self.coords)) {
            f.w.flush()?;
        }
        let mut b = Vec::new();
        let i64s = |b: &mut Vec<u8>, v: &[i64]| v.iter().for_each(|x| b.extend_from_slice(&x.to_le_bytes()));
        b.extend_from_slice(&32u64.to_le_bytes());
        i64s(&mut b, &[0, i64::MAX - 1, 0, i64::MAX - 1]);
        b.extend_from_slice(&(self.mbrs.len() as i64).to_le_bytes());
        for m in &self.mbrs {
            i64s(&mut b, m);
        }
        b.extend_from_slice(&(self.bounds.len() as i64).to_le_bytes());
        for m in &self.bounds {
            i64s(&mut b, m);
        }
        for o in &self.tile_offsets {
            b.extend_from_slice(&(o.len() as i64).to_le_bytes());
            i64s(&mut b, o);
        }
        for o in &self.tile_var_offsets {
            b.extend_from_slice(&(o.len() as i64).to_le_bytes());
            i64s(&mut b, o);
        }
        for s in &self.tile_var_sizes {
            b.extend_from_slice(&(s.len() as i64).to_le_bytes());
            s.iter().for_each(|x| b.extend_from_slice(&x.to_le_bytes()));
        }
        b.extend_from_slice(&(if self.cells == 0 { 0 } else { last_n }).to_le_bytes());
        let mut gz = GzEncoder::new(File::create(self.tmp.join("__book_keeping.tdb.gz"))?, Compression::default());
        gz.write_all(&b)?;
        gz.finish()?;
        drop(self.files);
        drop(self.var_files);
        std::fs::rename(&self.tmp, &self.done)?;
        File::create(self.done.join("__tiledb_fragment.tdb"))?;
        Ok(self.cells)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gdb_store;

    #[test]
    fn round_trip() {
        let dir = std::env::temp_dir().join(format!("gdbw_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let attrs = vec![attr("END", INT64, 1), attr("REF", CHAR, TILEDB_VAR_NUM), attr("QUAL", FLOAT32, 1), attr("GT", INT32, TILEDB_VAR_NUM)];
        write_schema(&dir.join("__array_schema.tdb"), "ws/x", &attrs).unwrap();
        let mut w = FragmentWriter::create(&dir, &attrs).unwrap();
        let mut want = Vec::new();
        for i in 0..2500i64 {
            let (row, col) = (i % 3, 100 + i / 3);
            let refa: Vec<u8> = if i % 7 == 0 { vec![] } else { b"ACGT"[..(i % 4 + 1) as usize].to_vec() };
            let gt: Vec<u8> = if i % 5 == 0 { vec![] } else { [0i32, 0, 1].iter().flat_map(|x| x.to_le_bytes()).collect() };
            let vals = vec![(col + 3).to_le_bytes().to_vec(), refa, (i as f32).to_le_bytes().to_vec(), gt];
            want.push((row, col, vals.clone()));
            w.push(row, col, vals).unwrap();
        }
        assert_eq!(w.finish().unwrap(), 2500);
        let arr = gdb_store::open_array(&dir);
        let mut got = Vec::new();
        gdb_store::for_each_cell(&arr, i64::MIN, i64::MAX, |c| got.push((c.row, c.col, c.vals)));
        assert_eq!(got, want);
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
