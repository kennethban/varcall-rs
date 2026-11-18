//! GatherBamFiles (`variant-rs gatherbamfiles`): Picard `GatherBamFiles` on BAM inputs, i.e. htsjdk
//! `BamFileIoUtils.gatherWithBlockCopying`: the first file's blocks as they are (header included, EOF block
//! dropped), then for every further file the rest of the block holding its first record, deflated again into one
//! new block, and its later blocks copied as they are; one EOF block. The output is therefore the inputs' bytes
//! except one re-deflated block per further input (Picard deflates it with its own deflater). `CREATE_INDEX`:
//! the `.bai` of the output as written (htsjdk inflates the output stream into a `BAMIndexer`); `CREATE_MD5_FILE`:
//! MD5 of the output.

use crate::bgzf::{self, EOF_BLOCK};
use crate::tabix::{self, ContigIndex};
use std::fs::File;
use std::io::{Read, Seek, SeekFrom, Write};
use std::os::unix::fs::FileExt;
use std::sync::atomic::AtomicU64;
use std::sync::Arc;

const USAGE: &str = "usage: variant-rs gatherbamfiles -I a.bam [-I b.bam ...] -O out.bam [--CREATE_INDEX true] [--CREATE_MD5_FILE true] \
[--compression-level 2]";

fn bool_arg(v: &str) -> bool {
    match v.to_ascii_lowercase().as_str() {
        "true" | "t" | "yes" | "y" => true,
        "false" | "f" | "no" | "n" => false,
        _ => panic!("not a boolean: {v}"),
    }
}

/// `BlockCompressedInputStream.checkTermination`: whether the file ends with the EOF block.
fn has_terminator(f: &File, path: &str) -> bool {
    let len = f.metadata().unwrap().len();
    if len >= EOF_BLOCK.len() as u64 {
        let mut tail = [0u8; 28];
        f.read_exact_at(&mut tail, len - 28).unwrap();
        if tail == EOF_BLOCK {
            return true;
        }
    }
    // HEALTHY (a valid last block, no terminator): nothing to skip
    assert!(len >= 28, "{path} does not have a valid GZIP block at the end of the file.");
    false
}

/// htsjdk IOUtil.unrollPaths, as Picard expands its INPUTs: a path ending in one of `exts` is an input; any other
/// is a list file, one path per line (trimmed, blank lines skipped), expanded in place, recursively, in order.
pub fn unroll_inputs(inputs: &[String], exts: &[&str]) -> Vec<String> {
    let mut out = Vec::new();
    for p in inputs {
        let name = std::path::Path::new(p).file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
        if exts.iter().any(|e| name.ends_with(e)) {
            out.push(p.clone());
        } else {
            let text = std::fs::read_to_string(p).unwrap_or_else(|e| panic!("{p}: {e} (an input without a {} extension is read as a list of paths)", exts.join("/")));
            let lines: Vec<String> = text.lines().map(|l| l.trim().to_string()).filter(|l| !l.is_empty()).collect();
            out.extend(unroll_inputs(&lines, exts));
        }
    }
    out
}

/// The size of the BGZF block at file offset `coff`.
fn block_size(f: &File, coff: u64) -> usize {
    let mut h = [0u8; 18];
    f.read_exact_at(&mut h, coff).unwrap();
    assert!(h[0] == 31 && h[1] == 139 && h[3] & 4 != 0, "not a BGZF block at {coff}");
    u16::from_le_bytes([h[16], h[17]]) as usize + 1
}

/// Copies bytes [from, to) of `f` to `out`.
fn copy_range(f: &mut File, from: u64, to: u64, out: &mut impl Write) {
    f.seek(SeekFrom::Start(from)).unwrap();
    let mut left = to - from;
    let mut buf = vec![0u8; 8 << 20];
    while left > 0 {
        let n = (left as usize).min(buf.len());
        f.read_exact(&mut buf[..n]).unwrap();
        out.write_all(&buf[..n]).unwrap();
        left -= n as u64;
    }
}

/// A writer that also hashes what it writes.
struct Out {
    w: std::io::BufWriter<File>,
    md5: Option<crate::md5::Md5>,
}

impl Write for Out {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.w.write_all(b)?;
        if let Some(m) = &mut self.md5 {
            m.update(b);
        }
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.w.flush()
    }
}

/// The `.bai` of a finished BAM: inflate it and index every record at its virtual offsets (an end at a block's end
/// is the next block's start).
fn index_bam(path: &str, bai: &str) {
    std::fs::write(bai, bai_of(path)).unwrap_or_else(|e| panic!("{bai}: {e}"));
}

/// A coordinate-sorted BAM's `.bai`, built by reading the file once. Per batch of inflated blocks, one thread finds
/// the record boundaries (two loads per record); decoding each record's position, end and virtual offsets runs in
/// parallel; the index is then fed the records in file order.
pub(crate) fn bai_of(path: &str) -> Vec<u8> {
    use rayon::prelude::*;
    let inp = crate::markdup::qname::open_input(path);
    let names = inp.header.ref_names.clone();
    let len = inp.len;
    // blocks up to the last one (the EOF block)
    let rx = bgzf::spawn_block_reader(inp.file.clone(), vec![(inp.first >> 16, len - EOF_BLOCK.len() as u64)], true, Arc::new(AtomicU64::new(0)));
    // inflated data not yet consumed, with the blocks it came from: (position in buf, file offset)
    let mut buf: Vec<u8> = Vec::new();
    let mut starts: Vec<(usize, u64)> = Vec::new();
    let mut pos = (inp.first & 0xffff) as usize; // first record in the first block
    // per contig in order of first appearance: (tid, index, mapped, unmapped); `slot[tid]` = its place
    let mut idx: Vec<(usize, ContigIndex, u64, u64)> = Vec::new();
    let mut slot = vec![usize::MAX; names.len()];
    let mut cur: Option<(usize, ContigIndex)> = None;
    let mut n_no_coor = 0u64;
    let mut recs: Vec<(usize, usize)> = Vec::new();
    loop {
        let more = match rx.recv() {
            Ok(Some(blocks)) => {
                for b in blocks {
                    if b.data.is_empty() {
                        continue;
                    }
                    starts.push((buf.len(), b.coff));
                    buf.extend_from_slice(&b.data);
                }
                true
            }
            _ => false,
        };
        // the complete records in the buffer: (position of block_size, block_size)
        recs.clear();
        while buf.len() >= pos + 4 {
            let bs = u32::from_le_bytes(buf[pos..pos + 4].try_into().unwrap()) as usize;
            if buf.len() < pos + 4 + bs {
                break;
            }
            recs.push((pos, bs));
            pos += 4 + bs;
        }
        let (buf_r, starts_r, file) = (&buf, &starts, &inp.file);
        // (tid, beg, end, mapped, start voffset, end voffset), None for records without a coordinate
        let ents: Vec<Option<(usize, i64, i64, bool, u64, u64)>> = recs
            .par_iter()
            .with_min_len(2048)
            .map(|&(p, bs)| {
                let rec = &buf_r[p + 4..p + 4 + bs];
                let (tid, pos1, flag) = crate::bam::peek_core(rec);
                if tid < 0 || pos1 < 1 {
                    return None;
                }
                // the end: in the block of the record's last byte, or, when the record ends that block, the address
                // right after it (BlockCompressedInputStream.getFilePointer at a block's end)
                let end_p = p + 4 + bs;
                let k = starts_r.partition_point(|&(s, _)| s <= end_p - 1) - 1;
                let block_end = starts_r.get(k + 1).map_or(buf_r.len(), |x| x.0);
                let end_v = if end_p < block_end {
                    (starts_r[k].1 << 16) | (end_p - starts_r[k].0) as u64
                } else {
                    (starts_r[k].1 + block_size(file, starts_r[k].1) as u64) << 16
                };
                let j = starts_r.partition_point(|&(s, _)| s <= p) - 1;
                let start_v = (starts_r[j].1 << 16) | (p - starts_r[j].0) as u64;
                let beg = pos1 - 1;
                let mapped = flag & 0x4 == 0;
                let e = if mapped { crate::bam::query_end(rec).max(beg + 1) } else { beg + 1 };
                Some((tid as usize, beg, e, mapped, start_v, end_v))
            })
            .collect();
        for x in ents {
            let Some((t, beg, e, mapped, start_v, end_v)) = x else {
                n_no_coor += 1;
                continue;
            };
            if cur.as_ref().map_or(true, |(ct, _)| *ct != t) {
                if let Some((ct, mut ci)) = cur.take() {
                    ci.finish();
                    idx[slot[ct]].1.merge(ci);
                }
                if slot[t] == usize::MAX {
                    slot[t] = idx.len();
                    idx.push((t, ContigIndex::new(&names[t]), 0, 0));
                }
                cur = Some((t, ContigIndex::new(&names[t])));
            }
            cur.as_mut().unwrap().1.push(beg, e, start_v, end_v);
            let entry = &mut idx[slot[t]];
            if mapped {
                entry.2 += 1;
            } else {
                entry.3 += 1;
            }
        }
        // drop the consumed blocks (what is left is part of one record)
        let k = starts.partition_point(|&(s, _)| s <= pos);
        if k > 1 {
            let cut = starts[k - 1].0;
            buf.drain(..cut);
            starts.drain(..k - 1);
            for s in starts.iter_mut() {
                s.0 -= cut;
            }
            pos -= cut;
        }
        if !more {
            break;
        }
    }
    if let Some((ct, mut ci)) = cur.take() {
        ci.finish();
        idx[slot[ct]].1.merge(ci);
    }
    tabix::bai_bytes(names.len(), &idx, n_no_coor)
}

pub fn run(argv: &[String]) {
    let t0 = std::time::Instant::now();
    let (mut inputs, mut out_path, mut index, mut md5, mut level) = (Vec::<String>::new(), String::new(), false, false, 2i32);
    let mut i = 0;
    while i < argv.len() {
        let v = || argv.get(i + 1).cloned().unwrap_or_else(|| panic!("{} needs a value\n{USAGE}", argv[i]));
        match argv[i].as_str() {
            "-I" | "--INPUT" => inputs.push(v()),
            "-O" | "--OUTPUT" => out_path = v(),
            "--CREATE_INDEX" => index = bool_arg(&v()),
            "--CREATE_MD5_FILE" => md5 = bool_arg(&v()),
            "--compression-level" => level = v().parse().expect("--compression-level"),
            x => panic!("unknown argument {x}\n{USAGE}"),
        }
        i += 2;
    }
    assert!(!inputs.is_empty() && !out_path.is_empty(), "{USAGE}");
    let inputs = unroll_inputs(&inputs, &[".bam", ".sam", ".cram"]);
    let mut out = Out {
        w: std::io::BufWriter::with_capacity(8 << 20, File::create(&out_path).unwrap_or_else(|e| panic!("{out_path}: {e}"))),
        md5: md5.then(crate::md5::Md5::default),
    };
    for (k, p) in inputs.iter().enumerate() {
        let inp = crate::markdup::qname::open_input(p);
        let mut f = File::open(p).unwrap_or_else(|e| panic!("{p}: {e}"));
        let term = has_terminator(&f, p);
        let len = inp.len - if term { EOF_BLOCK.len() as u64 } else { 0 };
        if k == 0 {
            copy_range(&mut f, 0, len, &mut out);
        } else {
            // the rest of the block holding the first record, as one new block; then the later blocks
            let (coff, uoff) = (inp.first >> 16, (inp.first & 0xffff) as usize);
            if coff >= len {
                continue; // no records (a header-only BAM: its first record offset is the EOF block)
            }
            let bsize = block_size(&f, coff);
            let mut raw = vec![0u8; bsize];
            f.read_exact_at(&mut raw, coff).unwrap();
            let data = bgzf::inflate_raw_block(&raw);
            if uoff < data.len() {
                out.write_all(&bgzf::deflate_block(&data[uoff..], level)).unwrap();
            }
            copy_range(&mut f, coff + bsize as u64, len, &mut out);
        }
    }
    out.write_all(&EOF_BLOCK).unwrap();
    out.flush().unwrap();
    if let Some(m) = out.md5.take() {
        std::fs::write(format!("{out_path}.md5"), m.hex()).unwrap_or_else(|e| panic!("{out_path}.md5: {e}"));
    }
    drop(out);
    if index {
        // IOUtil.basename(output) + ".bai" beside it
        let bai = format!("{}.bai", out_path.strip_suffix(".bam").unwrap_or(&out_path));
        index_bam(&out_path, &bai);
    }
    eprintln!("gatherbamfiles: {} inputs -> {out_path} ({:.1}s)", inputs.len(), t0.elapsed().as_secs_f64());
}
