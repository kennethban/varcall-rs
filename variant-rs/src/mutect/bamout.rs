//! `--bam-output`: GATK's `HaplotypeBAMWriter` (writer type CALLED_HAPLOTYPES) for Mutect2.
//!
//! Per region with calls (`writeReadsAlignedToHaplotypes`): the called haplotypes as artificial reads, in the
//! iteration order of the `HashSet` they are copied into, then every read of the likelihoods, sample by sample,
//! realigned to its best haplotype, with the tags GATK sets on the way (`HC` best-haplotype hash, `AR` alignment
//! region, `CR` callable region, `XA` supported alleles appended after any existing value, `BI` / `BD` from the
//! overlapping-pair cleaning). htsjdk then re-encodes every record's tags (sorted by binary tag, integers
//! retyped) and, as the writer is not presorted, sorts the records with `SAMRecordCoordinateComparator` (stable)
//! and writes `<out minus .bam>.bai`.

use super::call::BamoutInfo;
use super::region::SampleLikelihoods;
use hc_assembly::Haplotype;

/// A record before its final name is known (haplotype records are numbered in region order).
pub struct BRec {
    ref_id: i32,
    /// 1-based
    pos: i64,
    flag: u16,
    mapq: u8,
    /// None: haplotype record `HC_<counter>_<index>` with this index
    name: Option<Vec<u8>>,
    hap_index: usize,
    mate_ref: i32,
    /// raw BAM next_pos (0-based, -1 = none)
    mate_pos0: i32,
    tlen: i32,
    cigar: Vec<(u8, u32)>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    tags: Vec<u8>,
}

/// `Haplotype.hashCode()` = `Arrays.hashCode(bases)`.
fn hap_hash(bases: &[u8]) -> i32 {
    bases.iter().fold(1i32, |h, &b| h.wrapping_mul(31).wrapping_add(b as i8 as i32))
}

/// Bucket order of a `java.util.HashMap` table of `cap` entries filled in `order`: stable by bucket (resizes keep
/// the relative order inside a bucket).
fn by_bucket(order: &[usize], hash: &dyn Fn(usize) -> i32, cap: u32) -> Vec<usize> {
    let mut v: Vec<(u32, usize)> = order
        .iter()
        .map(|&i| {
            let h = hash(i) as u32;
            ((h ^ (h >> 16)) & (cap - 1), i)
        })
        .collect();
    v.sort_by_key(|x| x.0);
    v.into_iter().map(|x| x.1).collect()
}

/// `new HashSet<>(calledHaplotypes)` iteration order, `calledHaplotypes` being a default `HashSet` filled in
/// `insertion` order.
fn called_set_order(insertion: &[usize], haps: &[Haplotype]) -> Vec<usize> {
    let mut first: Vec<usize> = Vec::new();
    for &h in insertion {
        if !first.contains(&h) {
            first.push(h);
        }
    }
    let n = first.len();
    let hash = |i: usize| hap_hash(&haps[i].bases);
    // per bucket at most 7 entries: no treeified bins
    let mut cap1: u32 = 16;
    while n as u32 > cap1 * 3 / 4 {
        cap1 *= 2;
    }
    let order1 = by_bucket(&first, &hash, cap1);
    // HashSet(Collection c): new HashMap<>(Math.max((int) (c.size()/.75f) + 1, 16)), table = tableSizeFor
    let init = (((n as f32) / 0.75f32) as i32 + 1).max(16) as u32;
    let cap2 = init.next_power_of_two();
    let order2 = by_bucket(&order1, &hash, cap2);
    for c in [cap1, cap2] {
        let mut counts = std::collections::HashMap::new();
        for &i in &first {
            let h = hash(i) as u32;
            *counts.entry((h ^ (h >> 16)) & (c - 1)).or_insert(0) += 1;
        }
        assert!(counts.values().all(|&k| k < 8), "--bam-output: a HashSet bin with 8+ haplotypes (treeified bins are not modelled)");
    }
    order2
}

/// One decoded tag: two-character name, BAM type, value bytes (as in the record), integer value if integral.
struct Tag {
    name: [u8; 2],
    ty: u8,
    val: Vec<u8>,
    int: Option<i64>,
}

fn tag_size(t: u8) -> usize {
    match t {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        _ => 0,
    }
}

fn decode_tags(aux: &[u8]) -> Vec<Tag> {
    let mut v: Vec<Tag> = Vec::new();
    let mut o = 0;
    while o + 3 <= aux.len() {
        let name = [aux[o], aux[o + 1]];
        let ty = aux[o + 2];
        o += 3;
        let s = o;
        let mut int = None;
        match ty {
            b'Z' | b'H' => {
                while aux[o] != 0 {
                    o += 1;
                }
                o += 1;
            }
            b'B' => {
                let n = u32::from_le_bytes(aux[o + 1..o + 5].try_into().unwrap()) as usize;
                o += 5 + n * tag_size(aux[o]);
            }
            _ => {
                int = match ty {
                    b'c' => Some(aux[o] as i8 as i64),
                    b'C' => Some(aux[o] as i64),
                    b's' => Some(i16::from_le_bytes([aux[o], aux[o + 1]]) as i64),
                    b'S' => Some(u16::from_le_bytes([aux[o], aux[o + 1]]) as i64),
                    b'i' => Some(i32::from_le_bytes(aux[o..o + 4].try_into().unwrap()) as i64),
                    b'I' => Some(u32::from_le_bytes(aux[o..o + 4].try_into().unwrap()) as i64),
                    _ => None,
                };
                o += tag_size(ty);
            }
        }
        // a repeated tag replaces the earlier value
        v.retain(|t| t.name != name);
        v.push(Tag { name, ty, val: aux[s..o].to_vec(), int });
    }
    v
}

fn set_z(tags: &mut Vec<Tag>, name: &[u8; 2], s: &[u8]) {
    tags.retain(|t| &t.name != name);
    let mut val = s.to_vec();
    val.push(0);
    tags.push(Tag { name: *name, ty: b'Z', val, int: None });
}

fn set_int(tags: &mut Vec<Tag>, name: &[u8; 2], v: i64) {
    tags.retain(|t| &t.name != name);
    tags.push(Tag { name: *name, ty: b'i', val: Vec::new(), int: Some(v) });
}

fn get_z<'a>(tags: &'a [Tag], name: &[u8; 2]) -> Option<&'a [u8]> {
    tags.iter().find(|t| &t.name == name && t.ty == b'Z').map(|t| &t.val[..t.val.len() - 1])
}

/// htsjdk's encoding: sorted by binary tag (`(c1 << 8) | c0`), integers with `BinaryTagCodec.getIntegerType`.
fn encode_tags(mut tags: Vec<Tag>) -> Vec<u8> {
    tags.sort_by_key(|t| ((t.name[1] as u16) << 8) | t.name[0] as u16);
    let mut out = Vec::new();
    for t in tags {
        out.extend_from_slice(&t.name);
        match t.int {
            Some(v) => {
                if v > i32::MAX as i64 {
                    out.push(b'I');
                    out.extend_from_slice(&(v as u32).to_le_bytes());
                } else if v > 65535 {
                    out.push(b'i');
                    out.extend_from_slice(&(v as i32).to_le_bytes());
                } else if v > 32767 {
                    out.push(b'S');
                    out.extend_from_slice(&(v as u16).to_le_bytes());
                } else if v > 255 {
                    out.push(b's');
                    out.extend_from_slice(&(v as i16).to_le_bytes());
                } else if v > 127 {
                    out.push(b'C');
                    out.push(v as u8);
                } else if v >= -128 {
                    out.push(b'c');
                    out.push(v as i8 as u8);
                } else if v >= -32768 {
                    out.push(b's');
                    out.extend_from_slice(&(v as i16).to_le_bytes());
                } else {
                    out.push(b'i');
                    out.extend_from_slice(&(v as i32).to_le_bytes());
                }
            }
            None => {
                out.push(t.ty);
                out.extend_from_slice(&t.val);
            }
        }
    }
    out
}

fn interval(contig: &str, s: i64, e: i64) -> String {
    format!("{contig}:{s}-{e}")
}

/// The bamout records of one region (empty when nothing was called).
#[allow(clippy::too_many_arguments)]
pub fn region_records(
    contig: &str,
    tid: i32,
    haps: &[Haplotype],
    ref_padded_start: i64,
    span: (i64, i64),
    samples: &[SampleLikelihoods],
    info: &BamoutInfo,
) -> Vec<BRec> {
    if info.called_insertion.is_empty() {
        return Vec::new();
    }
    let cr = interval(contig, span.0, span.1);
    let mut out = Vec::new();
    for (index, h) in called_set_order(&info.called_insertion, haps).into_iter().enumerate() {
        let hap = &haps[h];
        let mut tags = Vec::new();
        set_int(&mut tags, b"HC", hap_hash(&hap.bases) as i64);
        set_z(&mut tags, b"RG", b"ArtificialHaplotypeRG");
        set_z(&mut tags, b"CR", cr.as_bytes());
        out.push(BRec {
            ref_id: tid,
            pos: ref_padded_start + hap.align_start_wrt_ref,
            flag: 0x10,
            mapq: 60,
            name: None,
            hap_index: index,
            mate_ref: -1,
            mate_pos0: -1,
            tlen: 0,
            cigar: hap.cigar.clone(),
            seq: hap.bases.clone(),
            qual: vec![b'!'; hap.bases.len()],
            tags: encode_tags(tags),
        });
    }
    for (s, sl) in samples.iter().enumerate() {
        for (r, read) in sl.reads.iter().enumerate() {
            let raw = read.bam_rec.as_ref().expect("--bam-output: read without its record");
            let aux = &raw[crate::bamwrite::aux_start(raw)..];
            let mut tags = decode_tags(aux);
            assert!(
                !tags.iter().any(|t| &t.name == b"BI" || &t.name == b"BD") || !read.bi.is_empty(),
                "--bam-output: input reads with BI/BD tags are not supported"
            );
            if !read.bi.is_empty() {
                let q33 = |q: &[u8]| -> Vec<u8> { q.iter().map(|&x| x + 33).collect() };
                set_z(&mut tags, b"BD", &q33(&read.bd));
                set_z(&mut tags, b"BI", &q33(&read.bi));
            }
            let ra = sl.realign[r];
            if ra.informative && ra.realigned {
                set_int(&mut tags, b"HC", hap_hash(&haps[ra.best].bases) as i64);
            }
            let bh = &haps[ra.best];
            set_z(&mut tags, b"AR", interval(contig, bh.loc_start, bh.loc_end).as_bytes());
            set_z(&mut tags, b"CR", cr.as_bytes());
            for &(pos, e) in &info.xa[s][r] {
                let prev = get_z(&tags, b"XA").map(|v| format!("{}, ", String::from_utf8_lossy(v))).unwrap_or_default();
                set_z(&mut tags, b"XA", format!("{prev}{contig}:{pos}={e}").as_bytes());
            }
            let ri = |o: usize| i32::from_le_bytes(raw[o..o + 4].try_into().unwrap());
            out.push(BRec {
                ref_id: ri(0),
                pos: read.start,
                flag: read.flag,
                mapq: read.mapq,
                name: Some(read.name.as_bytes().to_vec()),
                hap_index: 0,
                mate_ref: ri(20),
                mate_pos0: ri(24),
                tlen: ri(28),
                cigar: read.cigar.clone(),
                seq: read.seq.clone(),
                qual: read.qual.clone(),
                tags: encode_tags(tags),
            });
        }
    }
    out
}

/// `GenomicIndexUtil.regionToBin(beg, end)` (0-based, end exclusive)
fn reg2bin(beg: i64, end: i64) -> u16 {
    let end = end - 1;
    for (sh, off) in [(14, 4681), (17, 585), (20, 73), (23, 9), (26, 1)] {
        if beg >> sh == end >> sh {
            return (off + (beg >> sh)) as u16;
        }
    }
    0
}

fn encode(r: &BRec, name: &[u8]) -> Vec<u8> {
    let ref_len: i64 = r.cigar.iter().filter(|c| matches!(c.0, b'M' | b'D' | b'N' | b'=' | b'X')).map(|c| c.1 as i64).sum();
    let beg = r.pos - 1;
    let end = if ref_len > 0 { beg + ref_len } else { beg + 1 };
    let mut b = Vec::with_capacity(64 + name.len() + r.seq.len() * 2 + r.tags.len());
    b.extend_from_slice(&r.ref_id.to_le_bytes());
    b.extend_from_slice(&(beg as i32).to_le_bytes());
    b.push((name.len() + 1) as u8);
    b.push(r.mapq);
    b.extend_from_slice(&reg2bin(beg, end).to_le_bytes());
    b.extend_from_slice(&(r.cigar.len() as u16).to_le_bytes());
    b.extend_from_slice(&r.flag.to_le_bytes());
    b.extend_from_slice(&(r.seq.len() as i32).to_le_bytes());
    b.extend_from_slice(&r.mate_ref.to_le_bytes());
    b.extend_from_slice(&r.mate_pos0.to_le_bytes());
    b.extend_from_slice(&r.tlen.to_le_bytes());
    b.extend_from_slice(name);
    b.push(0);
    for &(op, len) in &r.cigar {
        let code = b"MIDNSHP=X".iter().position(|&c| c == op).expect("cigar op") as u32;
        b.extend_from_slice(&((len << 4) | code).to_le_bytes());
    }
    let nt = |c: u8| -> u8 { b"=ACMGRSVTWYHKDBN".iter().position(|&x| x == c.to_ascii_uppercase()).unwrap_or(15) as u8 };
    for p in r.seq.chunks(2) {
        b.push((nt(p[0]) << 4) | if p.len() > 1 { nt(p[1]) } else { 0 });
    }
    b.extend_from_slice(&r.qual);
    b.extend_from_slice(&r.tags);
    b
}

/// `HaplotypeBAMDestination`'s header: the reads' dictionary, sort order coordinate, their read groups plus
/// `ArtificialHaplotypeRG`, their programs plus `HaplotypeBAMWriter` (no comments).
pub fn header_text(source: &str) -> String {
    let mut t = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for kind in ["@SQ\t", "@RG\t"] {
        for l in source.lines().filter(|l| l.starts_with(kind)) {
            t.push_str(&crate::bamwrite::htsjdk_rg_line(l));
            t.push('\n');
        }
    }
    t.push_str("@RG\tID:ArtificialHaplotypeRG\tSM:HC\tCN:BI\n");
    for l in source.lines().filter(|l| l.starts_with("@PG\t")) {
        t.push_str(l);
        t.push('\n');
    }
    t.push_str("@PG\tID:HaplotypeBAMWriter\n");
    t
}

/// Names the haplotype records (`HC_<counter>_<index>`, counter over the run), sorts as htsjdk's writer and writes
/// the BAM and its index. `regions`: each region's records, in GATK's region order.
pub fn write(path: &str, header: &str, ref_names: &[String], ref_lens: &[i64], regions: Vec<Vec<BRec>>) {
    let mut recs: Vec<(Vec<u8>, BRec)> = Vec::new();
    let mut counter = 1u64;
    for reg in regions {
        for r in reg {
            let name = match &r.name {
                Some(n) => n.clone(),
                None => {
                    let n = format!("HC_{counter}_{}", r.hap_index).into_bytes();
                    counter += 1;
                    n
                }
            };
            recs.push((name, r));
        }
    }
    // SAMRecordCoordinateComparator; Vec::sort_by is stable, like the SortingCollection
    recs.sort_by(|(na, a), (nb, b)| {
        a.ref_id
            .cmp(&b.ref_id)
            .then(a.pos.cmp(&b.pos))
            .then_with(|| {
                let (ra, rb) = (a.flag & 0x10 != 0, b.flag & 0x10 != 0);
                if ra != rb {
                    return if ra { std::cmp::Ordering::Greater } else { std::cmp::Ordering::Less };
                }
                na.cmp(nb)
                    .then(a.flag.cmp(&b.flag))
                    .then(a.mapq.cmp(&b.mapq))
                    .then(a.mate_ref.cmp(&b.mate_ref))
                    .then(a.mate_pos0.cmp(&b.mate_pos0))
                    .then(a.tlen.cmp(&b.tlen))
            })
    });
    let mut bam = crate::bamwrite::BamFile::create(path, header, ref_names, ref_lens, 5, false);
    let mut i = 0;
    while i < recs.len() {
        let tid = recs[i].1.ref_id;
        let mut w = crate::bamwrite::RecWriter::new(5, &bam.too_large);
        while i < recs.len() && recs[i].1.ref_id == tid {
            w.push(&encode(&recs[i].1, &recs[i].0));
            i += 1;
        }
        bam.append(tid, &w.finish());
    }
    eprintln!("bam output: {path}: {} records", bam.n_records);
    bam.finish(path);
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn haplotype_hash() {
        // the truth bamout's HC_1_0 / HC_2_1 haplotypes and their HC tags
        assert_eq!(hap_hash(b"TGTCCCTGATGTGGCTGGGGTTGGTAGTGTCATTGTGGTCC"), 433870875);
        assert_eq!(hap_hash(b"TGTCCCTGATGTGGCTGGGGGTGGTAGTGTCATTGTGGTCC"), 1981276302);
    }
}
