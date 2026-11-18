//! `GVCFWriter` / `GVCFBlockCombiner` / `HomRefBlock` (GATK 4.7, default GQ bands, no floor blocks),
//! plus the writer of the unbanded records (`GATK_RAW_CALLS_VCF`), including the htsjdk aliasing of
//! the block-starting genotype's PL array (see README, "PL aliasing").

use std::io::{self, Write};

use crate::model::{gq_from_pls, SiteRecord};

/// htsjdk `VCFConstants.MAX_GENOTYPE_QUAL`
pub const MAX_GENOTYPE_QUAL: i32 = 99;

/// Default `--gvcf-gq-bands`: 1..=60, 70, 80, 90, 99.
pub fn default_gq_bands() -> Vec<i32> {
    let mut v: Vec<i32> = (1..=60).collect();
    v.extend_from_slice(&[70, 80, 90, 99]);
    v
}

/// A variant call (`HaplotypeCallerGenotypingEngine` output) as the GVCF writer sees it.
#[derive(Clone, Debug)]
pub struct CallRecord {
    pub pos: i64,
    /// VC end (start + reference allele length - 1)
    pub end: i64,
    pub rec: hc_genotype::VcfRecord,
}

impl CallRecord {
    /// `GVCFBlockCombiner.submit`: the combiner treats the record as a hom-ref site when the genotype is
    /// hom-ref (or a no-call with PL[0] == 0) and the only alt allele is <NON_REF>.
    pub fn is_hom_ref_site(&self) -> bool {
        let r = &self.rec;
        let biallelic_non_ref = r.alleles.len() == 2 && r.alleles[1].as_slice() == b"<NON_REF>";
        let hom_ref = match r.gt {
            Some((0, 0)) => true,
            None => r.pl.first() == Some(&0),
            _ => false,
        };
        biallelic_non_ref && hom_ref
    }
}

#[derive(Clone, Debug)]
pub enum Item {
    Site(SiteRecord),
    Call(CallRecord),
}

impl Item {
    pub fn pos(&self) -> i64 {
        match self {
            Item::Site(s) => s.pos,
            Item::Call(c) => c.pos,
        }
    }
}

struct Block {
    start: i64,
    end: i64,
    ref_allele: Vec<u8>,
    /// QUAL column of the starting record (a block started by a call keeps the call's QUAL)
    qual: Option<String>,
    lower: i32,
    upper: i32,
    min_pls: Vec<i32>,
    min_gq: i32,
    dps: Vec<i32>,
    /// index in the current region buffer of the record whose PL array the block aliases
    alias: Option<usize>,
}

/// GVCF writer: banded output (`out`) and, optionally, the unbanded per-position records (`raw`).
pub struct GvcfWriter<W: Write, R: Write> {
    contig: String,
    out: W,
    raw: Option<R>,
    bands: Vec<(i32, i32)>,
    block: Option<Block>,
    next_available_start: Option<i64>,
    region: Vec<Item>,
    line: Vec<u8>,
    pub n_out: u64,
    pub n_raw: u64,
}

fn push_int(b: &mut Vec<u8>, v: i64) {
    let mut buf = itoa_buf(v);
    b.append(&mut buf);
}

fn itoa_buf(v: i64) -> Vec<u8> {
    if v == 0 {
        return vec![b'0'];
    }
    let neg = v < 0;
    let mut x = v.unsigned_abs();
    let mut tmp = [0u8; 24];
    let mut i = tmp.len();
    while x > 0 {
        i -= 1;
        tmp[i] = b'0' + (x % 10) as u8;
        x /= 10;
    }
    let mut out = Vec::with_capacity(tmp.len() - i + 1);
    if neg {
        out.push(b'-');
    }
    out.extend_from_slice(&tmp[i..]);
    out
}

fn push_ints(b: &mut Vec<u8>, v: &[i32]) {
    for (i, &x) in v.iter().enumerate() {
        if i > 0 {
            b.push(b',');
        }
        push_int(b, x as i64);
    }
}

/// The data line htsjdk writes for a reference-confidence record (`GT:AD:DP:GQ:PL`).
pub fn site_line(b: &mut Vec<u8>, contig: &str, s: &SiteRecord) {
    b.extend_from_slice(contig.as_bytes());
    b.push(b'\t');
    push_int(b, s.pos);
    b.extend_from_slice(b"\t.\t");
    b.push(s.ref_base);
    b.extend_from_slice(b"\t<NON_REF>\t.\t.\t.\tGT:AD:DP:GQ:PL\t0/0:");
    push_int(b, s.ad[0] as i64);
    b.push(b',');
    push_int(b, s.ad[1] as i64);
    b.push(b':');
    push_int(b, s.dp() as i64);
    b.push(b':');
    push_int(b, s.gq.min(MAX_GENOTYPE_QUAL) as i64);
    b.push(b':');
    push_ints(b, &s.pl);
}

/// commons-math `Median` (LEGACY percentile) of ints, as `MathUtils.median` computes it.
pub fn median(v: &[i32]) -> f64 {
    let mut s: Vec<f64> = v.iter().map(|&x| x as f64).collect();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = s.len();
    if n == 1 {
        return s[0];
    }
    let pos = 0.5 * (n as f64 + 1.0);
    let fpos = pos.floor();
    let int_pos = fpos as usize;
    let dif = pos - fpos;
    if pos < 1.0 {
        return s[0];
    }
    if pos >= n as f64 {
        return s[n - 1];
    }
    let lower = s[int_pos - 1];
    let upper = s[int_pos];
    lower + dif * (upper - lower)
}

impl<W: Write, R: Write> GvcfWriter<W, R> {
    pub fn new(contig: &str, out: W, raw: Option<R>, gq_bands: &[i32]) -> Self {
        // GVCFBlockCombiner.parsePartitions
        let mut bands = Vec::new();
        let mut last = 0;
        for &v in gq_bands {
            bands.push((last, v));
            last = v;
        }
        if last <= MAX_GENOTYPE_QUAL {
            bands.push((last, MAX_GENOTYPE_QUAL + 1));
        }
        GvcfWriter {
            contig: contig.to_string(),
            out,
            raw,
            bands,
            block: None,
            next_available_start: None,
            region: Vec::new(),
            line: Vec::with_capacity(256),
            n_out: 0,
            n_raw: 0,
        }
    }

    fn band_of(&self, gq: i32) -> (i32, i32) {
        *self.bands.iter().find(|&&(lo, hi)| gq >= lo && gq < hi).expect("GQ outside every band")
    }

    fn write_block(&mut self, b: Block) -> io::Result<()> {
        let l = &mut self.line;
        l.clear();
        l.extend_from_slice(self.contig.as_bytes());
        l.push(b'\t');
        push_int(l, b.start);
        l.extend_from_slice(b"\t.\t");
        l.extend_from_slice(&b.ref_allele);
        l.extend_from_slice(b"\t<NON_REF>\t");
        l.extend_from_slice(b.qual.as_deref().unwrap_or(".").as_bytes());
        l.extend_from_slice(b"\t.\tEND=");
        push_int(l, b.end);
        // HomRefBlock.createHomRefGenotype: GT, then DP / GQ / MIN_DP / PL (DP and MIN_DP only with DPs)
        if b.dps.is_empty() {
            l.extend_from_slice(b"\tGT:GQ:PL\t0/0:");
        } else {
            l.extend_from_slice(b"\tGT:DP:GQ:MIN_DP:PL\t0/0:");
            push_int(l, hc_genotype::javamath::java_round(median(&b.dps)));
            l.push(b':');
        }
        push_int(l, b.min_gq.min(MAX_GENOTYPE_QUAL) as i64);
        l.push(b':');
        if !b.dps.is_empty() {
            push_int(l, *b.dps.iter().min().unwrap() as i64);
            l.push(b':');
        }
        push_ints(l, &b.min_pls);
        l.push(b'\n');
        self.out.write_all(l)?;
        self.n_out += 1;
        crate::stats::inc(&crate::stats::GVCF_BLOCKS);
        Ok(())
    }

    fn emit_current_block(&mut self) -> io::Result<()> {
        if let Some(b) = self.block.take() {
            self.write_block(b)?;
        }
        Ok(())
    }

    /// `GVCFBlockCombiner.submit`
    fn submit(&mut self, idx: usize) -> io::Result<()> {
        // (reference allele, GQ (-1 = missing), DP, PL, QUAL) of a hom-ref site
        type HomRef = (Vec<u8>, i32, Option<i32>, Vec<i32>, Option<String>);
        let (pos, end, hom_ref): (i64, i64, Option<HomRef>) = match &self.region[idx] {
            Item::Site(s) => (s.pos, s.pos, Some((vec![s.ref_base], s.gq, Some(s.dp()), s.pl.to_vec(), None))),
            Item::Call(c) => {
                let h = if c.is_hom_ref_site() {
                    crate::stats::inc(&crate::stats::GVCF_CALL_AS_HOM_REF_SITE);
                    let r = &c.rec;
                    Some((r.alleles[0].clone(), r.gq.unwrap_or(-1), r.dp, r.pl.clone(), Some(hc_genotype::javamath::format_qual(r.qual))))
                } else {
                    None
                };
                (c.pos, c.end, h)
            }
        };
        if let Some(b) = &self.block {
            // isContiguous: vc overlaps [block.start - 1, block.end + 1]
            if !(pos <= b.end + 1 && end >= b.start - 1) {
                self.emit_current_block()?;
            }
        }
        match hom_ref {
            Some((ref_allele, gq, dp, pl, qual)) => {
                // addHomRefSite
                if let Some(nas) = self.next_available_start {
                    if pos <= nas && end <= nas {
                        crate::stats::inc(&crate::stats::GVCF_SITES_DROPPED_UNDER_VARIANT);
                        return Ok(());
                    }
                    self.next_available_start = None;
                }
                let gqc = gq.min(MAX_GENOTYPE_QUAL);
                let mergeable = match &self.block {
                    Some(b) => gqc >= b.lower && gqc < b.upper && b.min_pls.len() == pl.len(),
                    None => false,
                };
                if mergeable {
                    let b = self.block.as_mut().unwrap();
                    if pos != b.end + 1 {
                        panic!("adding genotype at pos {} isn't contiguous with previous end {}", pos, b.end);
                    }
                    for (m, &p) in b.min_pls.iter_mut().zip(pl.iter()) {
                        *m = (*m).min(p);
                    }
                    b.min_gq = gq_from_pls(&b.min_pls);
                    b.end = pos; // getAttributeAsInt(END, start) of a single-position record
                    if let Some(d) = dp {
                        b.dps.push(d.max(0));
                    }
                    // the block's minPLs *is* the starting genotype's PL array (FastGenotype.getPL() is not a copy)
                    if let Some(i) = b.alias {
                        crate::stats::inc(&crate::stats::RAW_PL_ALIASED);
                        let mp = b.min_pls.clone();
                        match &mut self.region[i] {
                            Item::Site(s) => s.pl.copy_from_slice(&mp),
                            Item::Call(c) => c.rec.pl.copy_from_slice(&mp),
                        }
                    }
                } else {
                    self.emit_current_block()?;
                    // createNewBlock: a missing GQ counts as 0
                    let (lower, upper) = self.band_of(if gq < 0 { 0 } else { gqc });
                    let min_gq = gq_from_pls(&pl);
                    self.block = Some(Block {
                        start: pos,
                        end: pos,
                        ref_allele,
                        qual,
                        lower,
                        upper,
                        min_pls: pl,
                        min_gq,
                        dps: dp.map(|d| vec![d.max(0)]).unwrap_or_default(),
                        alias: Some(idx),
                    });
                }
            }
            None => {
                self.emit_current_block()?;
                self.next_available_start = Some(end);
                let Item::Call(c) = &self.region[idx] else { unreachable!() };
                crate::stats::inc(&crate::stats::GVCF_VARIANT_RECORDS);
                self.out.write_all(c.rec.to_line().as_bytes())?;
                self.out.write_all(b"\n")?;
                self.n_out += 1;
            }
        }
        Ok(())
    }

    /// Adds the records `callRegion` returned for one region, in order.
    pub fn add_region(&mut self, items: Vec<Item>) -> io::Result<()> {
        self.region = items;
        for i in 0..self.region.len() {
            self.submit(i)?;
        }
        // the raw writer receives the region's records after the GVCF writer has seen all of them
        if let Some(raw) = self.raw.as_mut() {
            let l = &mut self.line;
            for it in &self.region {
                l.clear();
                match it {
                    Item::Site(s) => site_line(l, &self.contig, s),
                    Item::Call(c) => l.extend_from_slice(c.rec.to_line().as_bytes()),
                }
                l.push(b'\n');
                raw.write_all(l)?;
                self.n_raw += 1;
            }
        }
        if let Some(b) = self.block.as_mut() {
            b.alias = None;
        }
        self.region.clear();
        Ok(())
    }

    /// The banded output (for in-memory writers).
    pub fn out_buffer(&self) -> &W {
        &self.out
    }
    /// The unbanded output, if any (for in-memory writers).
    pub fn raw_buffer(&self) -> &R {
        self.raw.as_ref().expect("no raw writer")
    }

    /// Returns the writers (call `finish` first).
    pub fn into_inner(self) -> (W, Option<R>) {
        (self.out, self.raw)
    }

    /// `GVCFWriter.close`: flushes the pending block.
    pub fn finish(&mut self) -> io::Result<()> {
        self.emit_current_block()?;
        self.out.flush()?;
        if let Some(r) = self.raw.as_mut() {
            r.flush()?;
        }
        Ok(())
    }
}
