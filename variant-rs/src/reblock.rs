//! ReblockGVCF (`variant-rs reblockgvcf`): GATK 4.7 `ReblockGVCF` on a single-sample, diploid HaplotypeCaller GVCF,
//! for its default arguments and WARP's (`-do-qual-approx --floor-blocks -GQB 20 -GQB 30 -GQB 40`); the options that
//! need the genotyping engine (`--drop-low-quals`) or other inputs are not supported. A literal port of
//! `ReblockGVCF.regenotypeVC` (reference blocks passed on; variants whose most likely genotype has no concrete alt
//! turned into GQ0 / subset hom-ref blocks; called variants cleaned up: uncalled alts dropped with
//! `AlleleSubsettingUtils.subsetAlleles`, trimmed, a reference block for the trimmed bases, NON_REF AD zeroed,
//! annotations copied / subset / added), of `ReblockingGVCFBlockCombiner` (the buffer of reference blocks trimmed
//! around variants) and of `GVCFBlockCombiner` / `HomRefBlock` (GQ bands, floor blocks). Records are written as
//! htsjdk's `VCFEncoder` writes decoded records.

use crate::output::OutFile;
use std::collections::BTreeMap;

const USAGE: &str = "usage: variant-rs reblockgvcf -R ref.fa -V in.g.vcf.gz -O out.g.vcf.gz [--do-qual-score-approximation] [--floor-blocks] \
[-GQB 20 -GQB 30 ...]";
const HEADER_LINES: &str = include_str!("../tables/reblock_header_lines.txt");
pub(crate) const NON_REF: &str = "<NON_REF>";
pub(crate) const SPAN_DEL: &str = "*";
pub(crate) const NO_CALL: usize = usize::MAX;
pub(crate) const MAX_GENOTYPE_QUAL: i32 = 99;
/// `GATKVariantContextUtils.SUM_GL_THRESH_NOCALL`
pub(crate) const SUM_GL_THRESH_NOCALL: f64 = -0.1;
/// `infoFieldAnnotationKeyNamesToRemove`
const INFO_TO_REMOVE: [&str; 8] = ["GVCFBlock", "HaplotypeScore", "InbreedingCoeff", "MLEAC", "MLEAF", "ExcessHet", "AS_InbreedingCoeff", "DS"];

/// code-path counters (printed with VARIANT_RS_REBLOCK_STATS=1): which parts of ReblockGVCF a run exercised
const STAT_NAMES: [&str; 14] = [
    "reference blocks in",
    "variants reblocked (no concrete alt in the best genotype)",
    "  GQ0 hom-ref (PL[0] != 0)",
    "  subset to ref + best alt (PL[0] == 0)",
    "  best alt chosen among several",
    "variants kept, nothing dropped",
    "variants kept, alleles dropped",
    "  dropped to hom-ref after subsetting",
    "  trimmed, reference block added",
    "  '*' dropped (overlap rule)",
    "genotype re-called from PLs",
    "buffer: blocks split around a variant",
    "buffer: variant dropped (buffer end before it)",
    "blocks dropped (inside the last variant)",
];
static STATS: [std::sync::atomic::AtomicU64; 14] = [const { std::sync::atomic::AtomicU64::new(0) }; 14];
fn stat(i: usize) {
    STATS[i].fetch_add(1, std::sync::atomic::Ordering::Relaxed);
}

fn refuse(why: String) -> ! {
    crate::gatk_cli::refuse(&why)
}

// ------------------------------------------------------------------------------------------------
// records

/// One sample's genotype (allele indices into the record's alleles, `NO_CALL` for '.').
#[derive(Clone, Debug)]
pub(crate) struct Gt {
    pub(crate) al: Vec<usize>,
    pub(crate) phased: bool,
    pub(crate) ad: Option<Vec<i32>>,
    pub(crate) dp: Option<i32>,
    pub(crate) gq: Option<i32>,
    pub(crate) pl: Option<Vec<i32>>,
    /// extended attributes as text
    pub(crate) ext: Vec<(String, String)>,
}

impl Gt {
    pub(crate) fn is_no_call(&self) -> bool {
        self.al.iter().all(|&a| a == NO_CALL)
    }
    pub(crate) fn is_called(&self) -> bool {
        !self.al.is_empty() && self.al.iter().all(|&a| a != NO_CALL)
    }
    /// htsjdk `Genotype.isHomRef`
    pub(crate) fn is_hom_ref(&self) -> bool {
        self.is_called() && self.al.iter().all(|&a| a == 0)
    }
    /// htsjdk `isHet` (two different called alleles) / `isHomVar` (all the same non-reference allele)
    pub(crate) fn is_het(&self) -> bool {
        self.is_called() && self.al.iter().any(|&a| a != self.al[0])
    }
    pub(crate) fn is_hom_var(&self) -> bool {
        self.is_called() && self.al[0] != 0 && self.al.iter().all(|&a| a == self.al[0])
    }
    pub(crate) fn ploidy(&self) -> usize {
        self.al.len()
    }
    pub(crate) fn remove_ext(&mut self, k: &str) {
        self.ext.retain(|(x, _)| x != k);
    }
    pub(crate) fn set_ext(&mut self, k: &str, v: String) {
        match self.ext.iter_mut().find(|(x, _)| x == k) {
            Some(e) => e.1 = v,
            None => self.ext.push((k.to_string(), v)),
        }
    }
}

#[derive(Clone, Debug)]
struct Vc {
    contig: String,
    start: i64,
    stop: i64,
    id: String,
    alleles: Vec<String>,
    /// QUAL as read (None: missing)
    qual: Option<String>,
    filter: String,
    info: Vec<(String, String)>,
    g: Gt,
}

impl Vc {
    fn attr(&self, k: &str) -> Option<&str> {
        self.info.iter().find(|(x, _)| x == k).map(|(_, v)| v.as_str())
    }
    fn attr_int(&self, k: &str, default: i64) -> i64 {
        self.attr(k).and_then(|v| v.trim().parse::<i64>().ok().or_else(|| v.trim().parse::<f64>().ok().map(|f| f as i64))).unwrap_or(default)
    }
    fn attr_double(&self, k: &str, default: f64) -> f64 {
        self.attr(k).and_then(|v| v.trim().parse::<f64>().ok()).unwrap_or(default)
    }
    fn set_attr(&mut self, k: &str, v: String) {
        match self.info.iter_mut().find(|(x, _)| x == k) {
            Some(e) => e.1 = v,
            None => self.info.push((k.to_string(), v)),
        }
    }
    fn n_alleles(&self) -> usize {
        self.alleles.len()
    }
    fn index_of(&self, a: &str) -> Option<usize> {
        self.alleles.iter().position(|x| x == a)
    }
    /// `ReblockGVCF.isHomRefBlock`: the only alt is NON_REF
    fn is_hom_ref_block(&self) -> bool {
        self.alleles.len() == 2 && self.alleles[1] == NON_REF
    }
}

pub(crate) fn is_symbolic(a: &str) -> bool {
    let b = a.as_bytes();
    a.len() > 1 && (b[0] == b'<' || b[b.len() - 1] == b'>' || a.contains('[') || a.contains(']') || b[0] == b'.' || b[b.len() - 1] == b'.')
}

/// `GATKVariantContextUtils.isConcreteAlt` for allele index `i` of `alleles`
pub(crate) fn is_concrete_alt(alleles: &[String], i: usize) -> bool {
    i != 0 && i != NO_CALL && !is_symbolic(&alleles[i]) && alleles[i] != SPAN_DEL
}

pub(crate) fn parse_ints(s: &str) -> Vec<i32> {
    s.split(',').map(|x| x.parse().unwrap_or_else(|_| panic!("not an integer list: {s}"))).collect()
}

fn parse_record(line: &str) -> Vc {
    let f: Vec<&str> = line.split('\t').collect();
    assert!(f.len() >= 10, "ReblockGVCF needs a single-sample GVCF: {line}");
    if f.len() > 10 {
        refuse("ReblockGVCF takes a single-sample GVCF".into());
    }
    let start: i64 = f[1].parse().unwrap();
    let alleles: Vec<String> = std::iter::once(f[3].to_ascii_uppercase()).chain(f[4].split(',').map(|a| if is_symbolic(a) || a == SPAN_DEL { a.to_string() } else { a.to_ascii_uppercase() })).collect();
    let info = parse_info(f[7]);
    let end = info.iter().find(|(k, _)| k == "END").map(|(_, v)| v.parse::<i64>().unwrap()).unwrap_or(start + alleles[0].len() as i64 - 1);
    let g = parse_gt(f[8], f[9]);
    if g.ploidy() != 2 {
        refuse(format!("{}:{}: ploidy {}: only diploid GVCFs are supported", f[0], start, g.ploidy()));
    }
    Vc { contig: f[0].to_string(), start, stop: end, id: f[2].to_string(), alleles, qual: if f[5] == "." { None } else { Some(f[5].to_string()) }, filter: f[6].to_string(), info, g }
}

/// One sample column under FORMAT `keys` (FT refused).
pub(crate) fn parse_gt(keys: &str, sample: &str) -> Gt {
    let vals: Vec<&str> = sample.split(':').collect();
    let mut g = Gt { al: Vec::new(), phased: false, ad: None, dp: None, gq: None, pl: None, ext: Vec::new() };
    for (k, key) in keys.split(':').enumerate() {
        let v = vals.get(k).copied().unwrap_or(".");
        match key {
            "GT" => {
                g.phased = v.contains('|');
                g.al = v.split(['/', '|']).map(|a| if a == "." { NO_CALL } else { a.parse().unwrap() }).collect();
            }
            _ if v == "." => {}
            "AD" => g.ad = Some(parse_ints(v)),
            "DP" => g.dp = Some(v.parse().unwrap()),
            "GQ" => g.gq = Some(v.parse().unwrap()),
            "PL" => g.pl = Some(parse_ints(v)),
            "FT" => refuse("genotype filters (FT): not supported".into()),
            _ => g.ext.push((key.to_string(), v.to_string())),
        }
    }
    g
}

/// INFO column -> (key, value) pairs in order ("" for a flag)
pub(crate) fn parse_info(s: &str) -> Vec<(String, String)> {
    if s == "." {
        return Vec::new();
    }
    s.split(';').filter(|s| !s.is_empty()).map(|kv| match kv.split_once('=') {
        Some((k, v)) => (k.to_string(), v.to_string()),
        None => (kv.to_string(), String::new()),
    }).collect()
}

/// htsjdk `VCFEncoder.appendGenotypeData` for decoded genotypes: FORMAT = GT then the sorted union of the keys,
/// missing values '.', trailing missing values of a sample dropped; GQ written capped at 99.
pub(crate) fn encode_genotypes(gs: &[Gt], out: &mut String) {
    let mut keys: Vec<&str> = Vec::new();
    for g in gs {
        keys.extend(g.ext.iter().map(|(k, _)| k.as_str()));
        if g.ad.is_some() {
            keys.push("AD");
        }
        if g.dp.is_some() {
            keys.push("DP");
        }
        if g.gq.is_some() {
            keys.push("GQ");
        }
        if g.pl.is_some() {
            keys.push("PL");
        }
    }
    keys.sort_unstable();
    keys.dedup();
    let ints = |x: &[i32]| x.iter().map(|i| i.to_string()).collect::<Vec<_>>().join(",");
    out.push_str("\tGT");
    for k in &keys {
        out.push(':');
        out.push_str(k);
    }
    for g in gs {
        out.push('\t');
        let gt = g.al.iter().map(|&a| if a == NO_CALL { ".".to_string() } else { a.to_string() }).collect::<Vec<_>>().join(if g.phased { "|" } else { "/" });
        out.push_str(&gt);
        let vals: Vec<Option<String>> = keys
            .iter()
            .map(|k| match *k {
                "AD" => g.ad.as_ref().map(|x| ints(x)),
                "DP" => g.dp.map(|x| x.to_string()),
                // IntGenotypeFieldAccessors: GQ written capped at 99
                "GQ" => g.gq.map(|x| x.min(MAX_GENOTYPE_QUAL).to_string()),
                "PL" => g.pl.as_ref().map(|x| ints(x)),
                _ => g.ext.iter().find(|(x, _)| x == k).map(|(_, v)| v.clone()),
            })
            .collect();
        let last = vals.iter().rposition(|v| v.is_some()).map_or(0, |i| i + 1);
        for v in &vals[..last] {
            out.push(':');
            out.push_str(v.as_deref().unwrap_or("."));
        }
    }
    out.push('\n');
}

/// htsjdk `VCFEncoder` for a decoded record: INFO sorted, FORMAT keys GT then sorted.
fn encode(v: &Vc, out: &mut String) {
    use std::fmt::Write;
    let _ = write!(out, "{}\t{}\t{}\t{}\t", v.contig, v.start, v.id, v.alleles[0]);
    if v.alleles.len() > 1 {
        out.push_str(&v.alleles[1..].join(","));
    } else {
        out.push('.');
    }
    out.push('\t');
    out.push_str(&v.qual.as_deref().map(crate::mergevcfs::format_qual).unwrap_or_else(|| ".".into()));
    out.push('\t');
    out.push_str(&v.filter);
    out.push('\t');
    if v.info.is_empty() {
        out.push('.');
    } else {
        let m: BTreeMap<&str, &str> = v.info.iter().map(|(k, x)| (k.as_str(), x.as_str())).collect();
        let s: Vec<String> = m.iter().map(|(k, x)| if x.is_empty() { k.to_string() } else { format!("{k}={x}") }).collect();
        out.push_str(&s.join(";"));
    }
    encode_genotypes(std::slice::from_ref(&v.g), out);
}

// ------------------------------------------------------------------------------------------------
// genotype likelihood arithmetic (diploid)

/// PL index of the diploid genotype (j <= k)
pub(crate) fn pl_index(j: usize, k: usize) -> usize {
    let (a, b) = if j <= k { (j, k) } else { (k, j) };
    b * (b + 1) / 2 + a
}

/// `GenotypesCache.get(2, i).asAlleleList`: allele indices of PL index i, ascending
pub(crate) fn genotype_alleles(i: usize) -> [usize; 2] {
    let mut k = 0;
    while (k + 1) * (k + 2) / 2 <= i {
        k += 1;
    }
    [i - k * (k + 1) / 2, k]
}

pub(crate) fn n_likelihoods(n_alleles: usize) -> usize {
    n_alleles * (n_alleles + 1) / 2
}

/// `MathUtils.minElementIndex` / `maxElementIndex` (first extreme)
pub(crate) fn min_index(v: &[i32]) -> usize {
    let mut b = 0;
    for i in 1..v.len() {
        if v[i] < v[b] {
            b = i;
        }
    }
    b
}

pub(crate) fn max_index_f(v: &[f64], from: usize) -> usize {
    let mut b = from;
    for i in from + 1..v.len() {
        if v[i] > v[b] {
            b = i;
        }
    }
    b
}

/// `AlleleSubsettingUtils.subsettedPLIndices(2, original, kept)`: `keep[i]` = index in the original alleles
pub(crate) fn subsetted_pl_indices(keep: &[usize]) -> Vec<usize> {
    let n = keep.len();
    let mut out = vec![0usize; n_likelihoods(n)];
    for k in 0..n {
        for j in 0..=k {
            out[pl_index(j, k)] = pl_index(keep[j], keep[k]);
        }
    }
    out
}

/// `GenotypeLikelihoods.getGQLog10FromLikelihoods`
pub(crate) fn gq_log10(chosen: usize, l: &[f64]) -> f64 {
    let mut qual = f64::NEG_INFINITY;
    for (i, &x) in l.iter().enumerate() {
        if i != chosen && x >= qual {
            qual = x;
        }
    }
    let qual = l[chosen] - qual;
    if qual < 0.0 {
        // MathUtils.normalizeFromLog10, then log10(1 - p)
        let m = l.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let s: f64 = l.iter().map(|&x| crate::jmath::pow(10.0, x - m)).sum();
        let p = crate::jmath::pow(10.0, l[chosen] - m) / s;
        crate::jmath::log10(1.0 - p)
    } else {
        -qual
    }
}

/// `GenotypeBuilder.log10PError` -> GQ
pub(crate) fn gq_of(log10: f64) -> i32 {
    hc_genotype::javamath::java_round(log10 * -10.0) as i32
}

/// `GenotypeLikelihoods.fromLog10Likelihoods(gls).getAsPLs()`
pub(crate) fn pls_of(gls: &[f64]) -> Vec<i32> {
    let m = gls.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    gls.iter().map(|&x| hc_genotype::javamath::java_round((-10.0 * (x - m)).min(i32::MAX as f64)) as i32).collect()
}

pub(crate) fn gls_of(pls: &[i32]) -> Vec<f64> {
    pls.iter().map(|&p| p as f64 / -10.0).collect()
}

#[derive(Clone, Copy, PartialEq)]
pub(crate) enum Assign {
    UsePls,
    BestMatch,
    /// GenotypeGVCFs' genotyping (PREFER_PLS) and merger (SET_TO_NO_CALL)
    PreferPls,
    SetToNoCall,
}

/// `GATKVariantContextUtils.makeGenotypeCall` (ploidy 2) on `gb`; `alleles` = the alleles to use (for NON_REF checks)
pub(crate) fn make_genotype_call(gb: &mut Gt, method: Assign, gls: Option<&[f64]>, alleles: &[String], orig: &Gt, orig_alleles: &[String]) {
    if method != Assign::SetToNoCall && (orig.is_hom_ref() || orig.is_no_call()) && orig.gq == Some(0) {
        gb.al = vec![NO_CALL; orig.ploidy().max(1)];
        if orig.dp == Some(0) {
            gb.pl = None;
            gb.dp = None;
            gb.ad = None;
            gb.gq = None;
            gb.ext.clear();
            return;
        }
    }
    // GATKVariantContextUtils.isInformative: sum of GLs < -0.1
    let informative = gls.is_some_and(|l| l.iter().sum::<f64>() < SUM_GL_THRESH_NOCALL);
    match method {
        Assign::SetToNoCall => gb.al = vec![NO_CALL; orig.ploidy().max(1)],
        Assign::PreferPls if !informative => gb.al = best_match_to_original(alleles, orig, orig_alleles),
        Assign::UsePls | Assign::PreferPls => match gls {
            None => gb.al = best_match_to_original(alleles, orig, orig_alleles),
            Some(l) => {
                let mx = max_index_f(l, 0);
                let fa = genotype_alleles(mx);
                let gq = gq_log10(mx, l);
                if fa.iter().any(|&a| alleles[a] == NON_REF) {
                    gb.al = vec![0, 0];
                    gb.pl = Some(vec![0; l.len()]);
                    gb.gq = Some(gq_of(0.0));
                } else if mx == 0 && gq > SUM_GL_THRESH_NOCALL {
                    gb.al = vec![NO_CALL; orig.ploidy().max(1)];
                } else {
                    gb.al = fa.to_vec();
                }
                if alleles.len() > 1 {
                    gb.gq = Some(gq_of(gq));
                }
            }
        },
        Assign::BestMatch => {
            if orig.gq == Some(0) && orig.pl.as_ref().is_none_or(|p| p[0] == 0) {
                gb.al = vec![NO_CALL; orig.ploidy().max(1)];
            } else {
                gb.al = best_match_to_original(alleles, orig, orig_alleles);
            }
        }
    }
}

/// `GATKVariantContextUtils.bestMatchToOriginalGT`: original alleles kept if present, else the reference
pub(crate) fn best_match_to_original(alleles: &[String], orig: &Gt, orig_alleles: &[String]) -> Vec<usize> {
    orig.al
        .iter()
        .map(|&a| if a == NO_CALL { NO_CALL } else { alleles.iter().position(|x| *x == orig_alleles[a]).unwrap_or(0) })
        .collect()
}

/// `AlleleSubsettingUtils.subsetAlleles` for one diploid genotype: `orig_alleles` -> `keep` (indices into them).
/// Allele-length FORMAT attributes (`allele_based`) other than AD / PL are dropped.
pub(crate) fn subset_alleles(g: &Gt, orig_alleles: &[String], keep: &[usize], method: Assign, allele_based: &[String]) -> Gt {
    let new_alleles: Vec<String> = keep.iter().map(|&i| orig_alleles[i].clone()).collect();
    let idx = subsetted_pl_indices(keep);
    let expected = n_likelihoods(orig_alleles.len());
    let mut new_gls: Option<Vec<f64>> = None;
    let mut new_log10_gq = f64::NEG_INFINITY;
    if let Some(pl) = &g.pl {
        if pl.len() == expected {
            let gls = gls_of(pl);
            let sub: Vec<f64> = idx.iter().map(|&i| gls[i]).collect();
            // MathUtils.scaleLogSpaceArrayForNumericalStability
            let m = sub.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let sub: Vec<f64> = sub.iter().map(|&x| x - m).collect();
            if sub.len() > 1 {
                let k = max_index_f(&sub, 0);
                new_log10_gq = gq_log10(k, &sub);
            } else {
                new_log10_gq = g.gq.unwrap_or(-1) as f64 / -10.0;
            }
            new_gls = Some(sub);
        }
    } else if let Some(gq) = g.gq {
        new_log10_gq = -0.1 * gq as f64;
    }
    let mut gb = g.clone();
    for k in ["PP", "GP", "PG"] {
        gb.remove_ext(k);
    }
    gb.ext.retain(|(k, _)| !allele_based.contains(k) || k == "AD" || k == "PL");
    gb.pl = None;
    gb.gq = None;
    if new_log10_gq != f64::NEG_INFINITY && g.gq.is_some() {
        gb.gq = Some(gq_of(new_log10_gq));
    }
    gb.pl = new_gls.as_ref().map(|l| pls_of(l));
    make_genotype_call(&mut gb, method, new_gls.as_deref(), &new_alleles, g, orig_alleles);
    if g.ext.iter().any(|(k, _)| k == "SAC") {
        refuse("SAC FORMAT attribute: not supported".into());
    }
    if let (Some(ad), true) = (&g.ad, gb.ad.is_some()) {
        gb.ad = Some(keep.iter().map(|&i| ad.get(i).copied().unwrap_or(0)).collect());
    }
    gb
}

// ------------------------------------------------------------------------------------------------
// blocks

/// `HomRefBlock`
struct Block {
    start_vc: Vc,
    lo: i32,
    hi: i32,
    min_pls: Option<Vec<i32>>,
    min_gq: i32,
    dps: Vec<i32>,
    end: i64,
}

/// `GATKVariantContextUtils.calculateGQFromPLs`
pub(crate) fn gq_from_pls(p: &[i32]) -> i32 {
    let (mut first, mut second) = (p[0], p[1]);
    if first > second {
        second = first;
        first = p[1];
    }
    for &c in &p[2..] {
        if c >= second {
            continue;
        }
        if c <= first {
            second = first;
            first = c;
        } else {
            second = c;
        }
    }
    second - first
}

impl Block {
    fn new(vc: &Vc, lo: i32, hi: i32) -> Block {
        let g = &vc.g;
        let min_pls = g.pl.clone();
        let min_gq = match &min_pls {
            Some(p) => gq_from_pls(p),
            None => g.gq.unwrap_or(-1),
        };
        let end = vc.attr_int("END", vc.start);
        Block { start_vc: vc.clone(), lo, hi, min_pls, min_gq, dps: g.dp.map(|d| vec![d.max(0)]).unwrap_or_default(), end }
    }
    fn within(&self, gq: i32) -> bool {
        gq >= self.lo && gq < self.hi
    }
    fn add(&mut self, pos: i64, new_end: i64, g: &Gt) {
        assert!(pos == self.end + 1, "adding genotype at pos {pos} isn't contiguous with previous end {}", self.end);
        let gq = g.gq.unwrap_or(-1);
        assert!(self.within(gq.min(MAX_GENOTYPE_QUAL)), "cannot add a genotype with GQ={gq}");
        match (&mut self.min_pls, &g.pl) {
            (None, p) => self.min_pls = p.clone(),
            (Some(m), Some(p)) => {
                assert_eq!(m.len(), p.len(), "trying to merge different PL array sizes");
                for (a, b) in m.iter_mut().zip(p) {
                    *a = (*a).min(*b);
                }
            }
            _ => {}
        }
        self.min_gq = match &self.min_pls {
            Some(p) => gq_from_pls(p),
            None => {
                if self.min_gq == -1 {
                    gq
                } else {
                    self.min_gq.min(gq)
                }
            }
        };
        self.end = new_end;
        if let Some(d) = g.dp {
            self.dps.push(d.max(0));
        }
    }
    /// `GVCFBlock.toVariantContext`
    fn to_vc(&self, floor: bool) -> Vc {
        let mut v = self.start_vc.clone();
        v.info = vec![("END".into(), self.end.to_string())];
        v.stop = self.end;
        let mut g = Gt { al: vec![0, 0], phased: false, ad: None, dp: None, gq: None, pl: None, ext: Vec::new() };
        if !floor {
            g.pl = self.min_pls.clone();
            g.gq = Some(self.min_gq);
            if !self.dps.is_empty() {
                g.ext.push(("MIN_DP".into(), self.dps.iter().min().unwrap().to_string()));
            }
        } else {
            g.gq = Some(self.lo);
        }
        if !self.dps.is_empty() {
            g.dp = Some(hc_genotype::javamath::java_round(hc_refconf::gvcf::median(&self.dps)) as i32);
        }
        v.g = g;
        v
    }
}

/// `GVCFBlockCombiner` + `ReblockingGVCFBlockCombiner`, writing finished records to `out`.
struct Combiner<'a> {
    bands: Vec<(i32, i32)>,
    floor: bool,
    fasta: &'a crate::reference::Fasta,
    // GVCFBlockCombiner
    block: Option<Block>,
    next_available: Option<(String, i64)>,
    // ReblockingGVCFBlockCombiner
    buffer: Vec<Vc>,
    vcf_output_end: i64,
    buffer_end: i64,
    contig: Option<String>,
    out: String,
    n_out: u64,
}

impl<'a> Combiner<'a> {
    fn ref_base(&self, contig: &str, pos: i64) -> String {
        let (b, _) = self.fasta.fetch(contig, pos, pos);
        assert_eq!(b.len(), 1, "reference base {contig}:{pos}");
        (b[0] as char).to_string()
    }
    /// `moveBuilderStart`
    fn move_start(&self, v: &mut Vc, new_start: i64) {
        v.alleles[0] = self.ref_base(&v.contig, new_start);
        v.start = new_start;
    }
    fn emit(&mut self, v: &Vc) {
        encode(v, &mut self.out);
        self.n_out += 1;
    }
    fn emit_block(&mut self) {
        if let Some(b) = self.block.take() {
            let v = b.to_vc(self.floor);
            self.emit(&v);
        }
    }
    /// `ReblockingGVCFBlockCombiner.isHomRef`
    fn is_hom_ref_r(g: &Gt) -> bool {
        match &g.pl {
            Some(p) => p[0] == 0,
            None => g.is_hom_ref() || g.is_no_call(),
        }
    }
    /// `GVCFBlockCombiner.submit`
    fn base_submit(&mut self, vc: Vc) {
        if let Some(b) = &self.block {
            // isContiguous: withinDistanceOf(block, 1)
            let contiguous = vc.contig == b.start_vc.contig && vc.start <= b.end + 1 && vc.stop >= b.start_vc.start - 1;
            if !contiguous {
                self.emit_block();
            }
        }
        let g = &vc.g;
        if (g.is_hom_ref() || (g.is_no_call() && g.pl.as_ref().is_some_and(|p| p[0] == 0))) && vc.alleles.iter().any(|a| a == NON_REF) && vc.alleles.len() == 2 {
            if let Some(done) = self.add_hom_ref_site(vc) {
                self.emit(&done);
            }
        } else {
            self.emit_block();
            self.next_available = Some((vc.contig.clone(), vc.stop));
            self.emit(&vc);
        }
    }
    /// `ReblockingGVCFBlockCombiner.addHomRefSite` (no --drop-low-quals) -> `GVCFBlockCombiner.addHomRefSite`
    fn add_hom_ref_site(&mut self, vc: Vc) -> Option<Vc> {
        let g = &vc.g;
        if !Self::is_hom_ref_r(g) {
            return None;
        }
        if g.pl.is_none() && g.gq.is_none() {
            refuse(format!("Homozygous reference genotypes must contain GQ or PL. Both are missing for hom ref genotype at {}:{}", vc.contig, vc.start));
        }
        if let Some((c, nas)) = &self.next_available {
            if vc.start <= *nas && vc.contig == *c && vc.stop <= *nas {
                stat(13);
                return None;
            }
            self.next_available = None;
        }
        let gq = g.gq.unwrap_or(-1);
        let can_merge = self.block.as_ref().is_some_and(|b| {
            b.within(gq.min(MAX_GENOTYPE_QUAL)) && b.min_pls.as_ref().is_none_or(|m| g.pl.as_ref().is_none_or(|p| p.len() == m.len()))
        });
        if can_merge {
            let end = vc.attr_int("END", vc.start);
            self.block.as_mut().unwrap().add(vc.start, end, g);
            None
        } else {
            let done = self.block.take().map(|b| b.to_vc(self.floor));
            let q = if g.gq.is_some() { gq.min(MAX_GENOTYPE_QUAL) } else { 0 };
            let &(lo, hi) = self.bands.iter().find(|(lo, hi)| q >= *lo && q < *hi).unwrap_or_else(|| panic!("GQ {q} didn't fit into any partition"));
            self.block = Some(Block::new(&vc, lo, hi));
            done
        }
    }
    /// `ReblockingGVCFBlockCombiner.submit`
    fn submit(&mut self, vc: Vc) {
        assert!(vc.start <= vc.stop, "Input variant context at position {}:{} has negative length", vc.contig, vc.start);
        match &self.contig {
            None => self.contig = Some(vc.contig.clone()),
            Some(c) if *c != vc.contig => {
                self.flush_buffer();
                self.contig = Some(vc.contig.clone());
                self.vcf_output_end = 0;
            }
            _ => {}
        }
        let hom_ref = Self::is_hom_ref_r(&vc.g);
        let mut new_block = vc.clone();
        if hom_ref && vc.start <= self.vcf_output_end {
            if vc.stop <= self.vcf_output_end {
                return;
            }
            let s = self.vcf_output_end + 1;
            self.move_start(&mut new_block, s);
        }
        let mut completed = vec![false; self.buffer.len()];
        let mut tail: Vec<Vc> = Vec::new();
        for k in 0..self.buffer.len() {
            let block_start = self.buffer[k].start;
            let variant_end = vc.stop;
            if block_start > variant_end {
                if !hom_ref {
                    self.base_submit(vc.clone());
                    self.vcf_output_end = self.vcf_output_end.max(variant_end);
                    self.finish_round(completed, tail);
                } else {
                    // add(newHomRefBlock), then removeAll(completed), addAll(tail), sort
                    self.buffer.push(new_block);
                    completed.push(false);
                    self.finish_round(completed, tail);
                }
                return;
            }
            let mut block_end = self.buffer[k].stop;
            let variant_start = vc.start;
            if block_end >= variant_start {
                block_end = self.trim_block_to_variant(&vc, k, &mut completed, &mut tail);
            }
            if block_start < variant_start && !hom_ref {
                let b = self.buffer[k].clone();
                self.base_submit(b);
                self.vcf_output_end = block_end;
                completed[k] = true;
            }
        }
        let mut k = 0;
        self.buffer.retain(|_| {
            k += 1;
            !completed[k - 1]
        });
        if hom_ref {
            if vc.is_hom_ref_block() && new_block.start < self.vcf_output_end {
                panic!("Reference positions added to buffer should not overlap positions already output to VCF. {} overlaps position {}:{} already emitted.", vc.start, vc.contig, self.vcf_output_end);
            }
            self.buffer_end = self.buffer_end.max(new_block.stop);
            self.buffer.push(new_block);
        } else if self.buffer_end >= vc.stop || self.buffer.is_empty() {
            let e = vc.stop;
            self.base_submit(vc);
            self.vcf_output_end = self.vcf_output_end.max(e);
        } else {
            stat(12);
        }
        self.buffer.extend(tail);
        self.sort_buffer();
    }
    /// the early-return path: `removeAll(completedBlocks)`, `addAll(tailBuffer)` (the sort follows)
    fn finish_round(&mut self, completed: Vec<bool>, tail: Vec<Vc>) {
        let mut k = 0;
        self.buffer.retain(|_| {
            k += 1;
            !completed[k - 1]
        });
        self.buffer.extend(tail);
        self.sort_buffer();
    }
    fn sort_buffer(&mut self) {
        // Collections.sort: stable
        self.buffer.sort_by_key(|b| b.start);
    }
    /// `trimBlockToVariant` on buffer[k]
    fn trim_block_to_variant(&mut self, vc: &Vc, k: usize, completed: &mut [bool], tail: &mut Vec<Vc>) -> i64 {
        let block_start = self.buffer[k].start;
        let variant_end = vc.stop;
        let mut block_end = self.buffer[k].stop;
        let variant_start = vc.start;
        if block_end > variant_end && block_start < variant_start {
            stat(11);
            let mut t = self.buffer[k].clone();
            self.move_start(&mut t, variant_end + 1);
            tail.push(t);
            let b = &mut self.buffer[k];
            b.stop = variant_start - 1;
            b.set_attr("END", (variant_start - 1).to_string());
            block_end = variant_start - 1;
        }
        if block_start < variant_start {
            let b = &mut self.buffer[k];
            b.set_attr("END", (variant_start - 1).to_string());
            b.stop = variant_start - 1;
        } else if vc.contig == self.buffer[k].contig && vc.start <= block_start && vc.stop >= block_end {
            completed[k] = true;
        } else {
            assert!(block_end >= variant_end + 1, "ref block end overlaps variant end");
            let mut b = self.buffer[k].clone();
            self.move_start(&mut b, variant_end + 1);
            self.buffer[k] = b;
        }
        assert!(self.buffer[k].start <= self.buffer[k].stop, "builder start follows stop");
        block_end
    }
    fn flush_buffer(&mut self) {
        let buf = std::mem::take(&mut self.buffer);
        for b in buf {
            let stop = b.stop;
            self.base_submit(b);
            self.vcf_output_end = stop;
        }
        self.buffer_end = 0;
    }
    /// `getVcfOutputEnd`: (contig, position), None before any output on the contig
    fn output_end(&self) -> Option<(&str, i64)> {
        let c = self.contig.as_deref()?;
        if self.vcf_output_end == 0 {
            None
        } else {
            Some((c, self.vcf_output_end))
        }
    }
    fn site_overlaps_buffer(&self, vc: &Vc) -> bool {
        if self.buffer.is_empty() {
            return false;
        }
        Some(vc.contig.as_str()) == self.contig.as_deref() && vc.start <= self.buffer_end && vc.start >= self.buffer[0].start
    }
    fn finish(&mut self) {
        self.flush_buffer();
        self.emit_block();
    }
}

// ------------------------------------------------------------------------------------------------
// ReblockGVCF

struct Tool<'a> {
    w: Combiner<'a>,
    qual_approx: bool,
    allele_based: Vec<String>,
}

/// (key, AS raw key, reducible, empty raw value) of the INFO annotations of StandardAnnotation + AS_StandardAnnotation
/// whose values are copied (`copyInfoAnnotations`): keys not allele specific are copied as they are.
const COPIED_KEYS: [&str; 19] = [
    "DP", "QD", "MQ", "RAW_MQandDP", "SOR", "ReadPosRankSum", "FS", "AC", "AF", "AN", "MQRankSum", "BaseQRankSum", "AS_MQ", "AS_SOR", "AS_FS",
    "AS_BaseQRankSum", "AS_ReadPosRankSum", "AS_MQRankSum", "AS_QD",
];
/// raw keys of the allele-specific annotations (subset with the alleles; the NON_REF entry set to the empty value)
const AS_KEYS: [(&str, bool, &str); 5] = [
    ("AS_RAW_MQ", true, "0.00"),
    ("AS_SB_TABLE", true, "0,0"),
    ("AS_RAW_BaseQRankSum", true, ""),
    ("AS_RAW_ReadPosRankSum", true, ""),
    ("AS_RAW_MQRankSum", true, ""),
];

impl<'a> Tool<'a> {
    /// `shouldBeReblocked`
    fn should_be_reblocked(&self, vc: &Vc) -> bool {
        let g = &vc.g;
        let Some(pls) = &g.pl else { return true };
        let fa = genotype_alleles(min_index(pls));
        !fa.iter().any(|&a| is_concrete_alt(&vc.alleles, a)) || fa.iter().any(|&a| vc.alleles[a] == NON_REF) || (g.pl.is_none() && g.gq.is_none())
    }

    fn regenotype(&mut self, vc: Vc) {
        if !vc.alleles.iter().any(|a| a == NON_REF) {
            refuse(format!("Variant Context at {}:{} does not contain a <NON-REF> allele. This tool is only intended for use with GVCFs.", vc.contig, vc.start));
        }
        if vc.is_hom_ref_block() {
            stat(0);
            if let Some((c, e)) = self.w.output_end() {
                if c == vc.contig && vc.stop <= e {
                    return;
                }
            }
            if vc.g.pl.is_none() {
                if vc.g.gq.is_some() {
                    // GATK adds such a block twice
                    self.w.submit(vc.clone());
                } else {
                    refuse(format!("Homozygous reference genotypes must contain GQ or PL. Both are missing for hom ref genotype at {}:{}", vc.contig, vc.start));
                }
            }
            self.w.submit(vc);
            return;
        }
        if self.should_be_reblocked(&vc) {
            stat(1);
            if let Some((c, e)) = self.w.output_end() {
                if c == vc.contig && vc.stop <= e {
                    return;
                }
            }
            if let Some(b) = self.low_qual_to_hom_ref(&vc) {
                self.w.submit(b);
            }
        } else if let Some(v) = self.clean_up(&vc) {
            self.w.submit(v);
        }
    }

    /// `lowQualVariantToGQ0HomRef`
    fn low_qual_to_hom_ref(&self, vc: &Vc) -> Option<Vc> {
        let (g, out_ref) = self.change_call_to_hom_ref(vc);
        let mut b = vc.clone();
        b.alleles = vec![out_ref, NON_REF.to_string()];
        b.g = g;
        if let Some((c, e)) = self.w.output_end() {
            if c == vc.contig && vc.start <= e {
                let new_start = e + 1;
                if new_start > vc.stop {
                    return None;
                }
                self.w.move_start(&mut b, new_start);
            }
        }
        b.filter = ".".into();
        b.qual = None;
        b.info = vec![("END".into(), vc.stop.to_string())];
        Some(b)
    }

    /// `changeCallToHomRefVersusNonRef`: the genotype and the output reference allele
    fn change_call_to_hom_ref(&self, vc: &Vc) -> (Gt, String) {
        let g = &vc.g;
        let mut gb;
        if g.pl.as_ref().is_none_or(|p| p[0] != 0) {
            stat(2);
            gb = g.clone();
            gb.pl = Some(vec![0; 3]);
            gb.gq = Some(0);
            gb.ad = None;
            gb.al = vec![0, 0];
            gb.ext.clear();
        } else {
            // calculateMostLikelyAlleles(vc, 2, 1, true): with one proper alt, all alleles; first non-reference is the best
            stat(3);
            let non_ref = vc.index_of(NON_REF);
            let n_proper = vc.n_alleles() - if non_ref.is_some() { 2 } else { 1 };
            let best_alt = if n_proper <= 1 {
                1
            } else {
                stat(4);
                // calculateLikelihoodSums(countAllelesWithoutHomRef = true), then the best proper alt (stable sort)
                let gls = gls_of(g.pl.as_ref().unwrap());
                let mut sums = vec![0f64; vc.n_alleles()];
                let ix = max_index_f(&gls, 1);
                let diff = (gls[ix] - gls[0]).abs();
                let ga = genotype_alleles(ix);
                for a in 1..vc.n_alleles() {
                    if ga.contains(&a) {
                        sums[a] += diff;
                    }
                }
                let mut cand: Vec<usize> = (1..vc.n_alleles()).filter(|&a| Some(a) != non_ref).collect();
                cand.sort_by(|&a, &b| sums[b].partial_cmp(&sums[a]).unwrap());
                cand[0]
            };
            let sub = subset_alleles(g, &vc.alleles, &[0, best_alt], Assign::BestMatch, &self.allele_based);
            gb = sub;
            gb.ext.clear();
            if gb.gq.is_none() {
                gb.gq = Some(0);
            }
            if gb.pl.is_none() {
                gb.pl = Some(vec![0; 3]);
            }
        }
        if vc.attr("DP").is_some() {
            let d = vc.attr_int("DP", 0) as i32;
            gb.dp = Some(d);
            gb.set_ext("MIN_DP", d.to_string());
        } else if let Some(ad) = &g.ad {
            let d: i32 = ad.iter().sum();
            gb.dp = Some(d);
            gb.set_ext("MIN_DP", d.to_string());
        }
        let in_ref = &vc.alleles[0];
        let out_ref = if in_ref.len() > 1 || g.al.iter().any(|&a| a == NO_CALL || vc.alleles[a] == SPAN_DEL) { in_ref[..1].to_string() } else { in_ref.clone() };
        gb.al = vec![0, 0];
        (gb, out_ref)
    }

    /// `getCalledGenotype`
    fn called_genotype(&self, vc: &Vc) -> Gt {
        let g = &vc.g;
        let pls = g.pl.as_ref().unwrap_or_else(|| panic!("Cannot verify called genotype without likelihoods or posteriors. Error at {}:{}", vc.contig, vc.start));
        let fa = genotype_alleles(min_index(pls));
        let mismatch = !g.al.iter().all(|a| fa.contains(a));
        if g.is_no_call() || mismatch {
            stat(10);
            let mut gb = g.clone();
            let gls = gls_of(pls);
            make_genotype_call(&mut gb, Assign::UsePls, Some(&gls), &vc.alleles, g, &vc.alleles);
            gb
        } else {
            g.clone()
        }
    }

    /// `getAllelesToDrop` (indices into vc.alleles)
    fn alleles_to_drop(&self, vc: &Vc, called: &Gt) -> Vec<usize> {
        let mut d: Vec<usize> = (1..vc.n_alleles()).filter(|&a| !is_symbolic(&vc.alleles[a]) && !called.al.contains(&a)).collect();
        if let Some(sd) = vc.index_of(SPAN_DEL) {
            if called.al.contains(&sd) && (self.w.site_overlaps_buffer(vc) || self.w.output_end().is_none_or(|(_, e)| e < vc.start)) {
                stat(9);
                d.push(sd);
            }
        }
        d
    }

    /// `cleanUpHighQualityVariant`
    fn clean_up(&mut self, vc: &Vc) -> Option<Vc> {
        let genotype = self.called_genotype(vc);
        let mut b = vc.clone();
        b.info.clear();
        b.g = genotype.clone();
        let drop = self.alleles_to_drop(vc, &genotype);
        let mut relevant: Vec<usize> = vec![0; vc.n_alleles()];
        stat(if drop.is_empty() { 5 } else { 6 });
        if !drop.is_empty() {
            let keep: Vec<usize> = (0..vc.n_alleles()).filter(|a| !drop.contains(a)).collect();
            let sub = subset_alleles(&vc.g, &vc.alleles, &keep, Assign::UsePls, &self.allele_based);
            if sub.is_hom_ref() || sub.gq.is_none() || sub.al.contains(&NO_CALL) {
                stat(7);
                if let Some(h) = self.low_qual_to_hom_ref(vc) {
                    self.w.submit(h);
                }
                return None;
            }
            b.g = sub;
            b.alleles = keep.iter().map(|&i| vc.alleles[i].clone()).collect();
            // trimAlleles(vc, false, true)
            trim_reverse(&mut b);
            relevant = keep.clone();
            let depth = if vc.attr("DP").is_some() { vc.attr_int("DP", 0) as i32 } else { genotype.dp.unwrap_or(0) };
            self.add_ref_block_if_necessary(vc, &drop, &b, depth);
        }
        // removeNonRefADs (in place: the annotations below see the zeroed AD too)
        let nr = b.index_of(NON_REF);
        if let (Some(ad), Some(nr)) = (b.g.ad.as_mut(), nr) {
            if ad.len() >= nr && ad[nr] > 0 {
                let c = ad[nr];
                ad[nr] = 0;
                b.g.dp = Some(match b.g.dp {
                    Some(d) => d - c,
                    None => ad.iter().sum(),
                });
            }
        }
        let attrs = self.compose_annotations(vc, &relevant, &b);
        b.info = attrs;
        b.filter = ".".into();
        Some(b)
    }

    /// `addRefBlockIfNecessary`
    fn add_ref_block_if_necessary(&mut self, orig: &Vc, drop: &[usize], trimmed: &Vc, depth: i32) {
        let old_len = orig.alleles[0].len() as i64;
        let new_len = trimmed.alleles[0].len() as i64;
        let out_end = self.w.output_end().map_or(-1, |(_, e)| e);
        if new_len >= old_len {
            return;
        }
        let Some(pls) = &orig.g.pl else { return };
        // shortest dropped alt (not *), first of equal length
        let shortest = drop.iter().copied().filter(|&a| orig.alleles[a] != SPAN_DEL).min_by_key(|&a| orig.alleles[a].len())
            .unwrap_or_else(|| panic!("No shortest ALT at {} across alleles", orig.start));
        let idx = subsetted_pl_indices(&[0, shortest]);
        let sub: Vec<i32> = idx.iter().map(|&i| pls[i]).collect();
        let m = *sub.iter().min().unwrap();
        let mut lk: Vec<i32> = sub.iter().map(|&x| x - m).collect();
        if lk[0] != 0 {
            // as GATK's loop: in place, so lk[0] is already 0 for the later entries (only the first changes)
            for i in 0..lk.len() {
                lk[i] = (lk[i] - lk[0]).max(0);
            }
        }
        let ref_start = (orig.stop - (old_len - new_len)).max(out_end) + 1;
        let new_ref = self.w.ref_base(&orig.contig, ref_start);
        // MathUtils.secondSmallestMinusSmallest(lk, 0)
        let gq = {
            let mut s = lk.clone();
            s.sort_unstable();
            if s.len() < 2 { 0 } else { s[1] - s[0] }
        };
        if ref_start > out_end && orig.stop > out_end {
            stat(8);
            let start = ref_start.max(out_end + 1);
            let v = Vc {
                contig: orig.contig.clone(),
                start,
                stop: orig.stop,
                id: ".".into(),
                alleles: vec![new_ref, NON_REF.into()],
                qual: None,
                filter: ".".into(),
                info: vec![("END".into(), orig.stop.to_string())],
                g: Gt { al: vec![0, 0], phased: false, ad: None, dp: Some(depth), gq: Some(gq), pl: Some(lk), ext: Vec::new() },
            };
            self.w.submit(v);
        }
    }

    /// `composeUpdatedAnnotations`
    fn compose_annotations(&self, vc: &Vc, relevant: &[usize], updated: &Vc) -> Vec<(String, String)> {
        let mut d: Vec<(String, String)> = Vec::new();
        let put = |d: &mut Vec<(String, String)>, k: &str, v: String| match d.iter_mut().find(|(x, _)| x == k) {
            Some(e) => e.1 = v,
            None => d.push((k.to_string(), v)),
        };
        // updateMQAnnotations
        if vc.attr("RAW_MQandDP").is_none() {
            if vc.attr("RAW_MQ").is_some() {
                refuse("the deprecated RAW_MQ annotation: not supported".into());
            }
            let mq = vc.attr_double("MQ", 60.0);
            let dp = vc.attr_int("DP", 0);
            let raw = hc_genotype::javamath::java_round(mq * mq * dp as f64) as i32;
            put(&mut d, "RAW_MQandDP", format!("{raw},{dp}"));
        }
        // copyInfoAnnotations
        let subsetting = relevant.len() < vc.n_alleles();
        if subsetting && vc.attr("AS_QD").is_some() {
            refuse(format!("{}:{}: AS_QD with alleles to drop: not supported", vc.contig, vc.start));
        }
        for k in COPIED_KEYS {
            if let Some(v) = vc.attr(k) {
                put(&mut d, k, v.to_string());
            }
        }
        for (k, reducible, empty) in AS_KEYS {
            let Some(v) = vc.attr(k) else { continue };
            if !subsetting {
                put(&mut d, k, v.to_string());
                continue;
            }
            let vals: Vec<&str> = if reducible { v.split('|').collect() } else { v.split(',').collect() };
            let offset = if reducible { 0 } else { 1 };
            let mut sub: Vec<String> = (offset..relevant.len())
                .map(|i| {
                    let old = relevant[i];
                    if old >= vals.len() + offset { String::new() } else { vals[old - offset].to_string() }
                })
                .collect();
            if vc.alleles[*relevant.last().unwrap()] == NON_REF {
                let l = sub.len() - 1;
                sub[l] = empty.to_string();
            }
            put(&mut d, k, sub.join(if reducible { "|" } else { "," }));
        }
        let ug = &updated.g;
        if self.qual_approx {
            if let Some(pl) = &ug.pl {
                put(&mut d, "QUALapprox", pl[0].to_string());
                let mut var_dp = qd_depth(ug);
                if var_dp == 0 {
                    var_dp = d.iter().find(|(k, _)| k == "DP").map(|(_, v)| v.parse().unwrap()).unwrap_or(1);
                }
                put(&mut d, "VarDP", var_dp.to_string());
                // AS_QUALapprox / AS_VarDP
                let mut quals = Vec::new();
                for a in 1..updated.n_alleles() {
                    if updated.alleles[a] == NON_REF || updated.alleles[a] == SPAN_DEL {
                        quals.push("0".to_string());
                        continue;
                    }
                    let s = subset_alleles(ug, &updated.alleles, &[0, a], Assign::BestMatch, &self.allele_based);
                    quals.push(s.pl.map_or("0".to_string(), |p| p[0].to_string()));
                }
                put(&mut d, "AS_QUALapprox", format!("|{}", quals.join("|")));
                if let Some(ad) = &ug.ad {
                    let mut depths = vec![0i32; ad.len()];
                    if ug.is_het() || ug.is_hom_var() {
                        let tot: i32 = ad.iter().sum();
                        if tot - ad[0] > 1 {
                            for (x, y) in depths.iter_mut().zip(ad) {
                                *x += y;
                            }
                        }
                    }
                    put(&mut d, "AS_VarDP", depths.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("|"));
                }
            }
        } else {
            for k in ["AS_VarDP", "QUALapprox"] {
                if let Some(v) = vc.attr(k) {
                    put(&mut d, k, v.to_string());
                }
            }
        }
        let has_ref = ug.al.contains(&0);
        put(&mut d, "RAW_GT_COUNT", if has_ref { "0,1,0" } else { "0,0,1" }.into());
        d
    }
}

/// `QualByDepth.getDepth` for one genotype (no likelihoods)
fn qd_depth(g: &Gt) -> i32 {
    if !g.is_het() && !g.is_hom_var() {
        return 0;
    }
    if let Some(ad) = &g.ad {
        let tot: i32 = ad.iter().sum();
        if tot != 0 {
            return if tot - ad[0] > 1 { tot } else { tot };
        }
    }
    g.dp.unwrap_or(0)
}

/// `GATKVariantContextUtils.trimAlleles(vc, false, true)`: common trailing bases of the non-symbolic alleles
/// removed (one kept when an allele would become empty); nothing when an allele is one base long (but '*').
fn trim_reverse(v: &mut Vc) {
    if reverse_trim(&mut v.alleles) > 0 {
        v.stop = v.start + v.alleles[0].len() as i64 - 1;
    }
}

/// `trimAlleles(vc, false, true)` on the alleles; returns the number of bases clipped from the end.
pub(crate) fn reverse_trim(alleles: &mut [String]) -> usize {
    if alleles.len() <= 1 || alleles.iter().any(|a| a.len() == 1 && a != SPAN_DEL) {
        return 0;
    }
    let seqs: Vec<&[u8]> = alleles.iter().filter(|a| !is_symbolic(a) && *a != SPAN_DEL).map(|a| a.as_bytes()).collect();
    // AlignmentUtils.normalizeAlleles(maxShift 0, trim): end trim first, then start trim (counted only for the
    // empty-allele rule)
    let mut ranges: Vec<(usize, usize)> = seqs.iter().map(|s| (0, s.len())).collect();
    let mut min_size = ranges.iter().map(|r| r.1 - r.0).min().unwrap();
    let mut end_trim = 0usize;
    while min_size > 0 && seqs.iter().zip(&ranges).all(|(s, r)| s[r.1 - 1] == seqs[0][ranges[0].1 - 1]) {
        ranges.iter_mut().for_each(|r| r.1 -= 1);
        min_size -= 1;
        end_trim += 1;
    }
    let mut start_trim = 0usize;
    while min_size > 0 && seqs.iter().zip(&ranges).all(|(s, r)| s[r.0] == seqs[0][ranges[0].0]) {
        ranges.iter_mut().for_each(|r| r.0 += 1);
        min_size -= 1;
        start_trim += 1;
    }
    let empty = ranges.iter().any(|r| r.1 == r.0);
    let end_clip = if empty && start_trim == 0 { end_trim - 1 } else { end_trim };
    if end_clip == 0 {
        return 0;
    }
    let seqs_owned: Vec<bool> = alleles.iter().map(|a| !is_symbolic(a) && a != SPAN_DEL).collect();
    for (a, t) in alleles.iter_mut().zip(seqs_owned) {
        if t {
            a.truncate(a.len() - end_clip);
        }
    }
    end_clip
}

// ------------------------------------------------------------------------------------------------
// driver

fn bands_of(gqb: &[i32]) -> Vec<(i32, i32)> {
    let mut bands = Vec::new();
    let mut last = 0;
    for &v in gqb {
        // GVCFBlockCombiner.parsePartitions: increasing, distinct, in (0, 100]
        if v <= last || v > MAX_GENOTYPE_QUAL + 1 {
            refuse(format!("GQBands are malformed: {gqb:?}"));
        }
        bands.push((last, v));
        last = v;
    }
    if last <= MAX_GENOTYPE_QUAL {
        bands.push((last, MAX_GENOTYPE_QUAL + 1));
    }
    bands
}

pub fn run(argv: &[String]) {
    let t0 = std::time::Instant::now();
    let (mut input, mut output, mut reference) = (None::<String>, None::<String>, None::<String>);
    let (mut qual_approx, mut floor, mut gqb) = (false, false, Vec::<i32>::new());
    let mut i = 0;
    while i < argv.len() {
        let v = || argv.get(i + 1).cloned().unwrap_or_else(|| panic!("{} needs a value\n{USAGE}", argv[i]));
        match argv[i].as_str() {
            "--do-qual-score-approximation" | "-do-qual-approx" => {
                qual_approx = true;
                i += 1;
                continue;
            }
            "--floor-blocks" => {
                floor = true;
                i += 1;
                continue;
            }
            "-V" | "--variant" => input = Some(v()),
            "-O" | "--output" => output = Some(v()),
            "-R" | "--reference" => reference = Some(v()),
            "-GQB" | "--gvcf-gq-bands" => gqb.push(v().parse().expect("-GQB")),
            x => panic!("unknown argument {x}\n{USAGE}"),
        }
        i += 2;
    }
    let (input, output, reference) = (input.expect(USAGE), output.expect(USAGE), reference.expect(USAGE));
    if gqb.is_empty() {
        gqb = vec![20, 100];
    }
    let fasta = crate::reference::Fasta::open(&reference);
    let mut lines = crate::mergevcfs::open_lines(&input);

    // header: the input's lines without GVCFBlock and single-sample-uninformative INFO lines, plus ReblockGVCF's
    let mut hdr: Vec<String> = Vec::new();
    let chrom = loop {
        let l = lines.next().unwrap_or_else(|| panic!("{input}: no #CHROM line"));
        if l.starts_with("#CHROM") {
            break l;
        }
        hdr.push(l);
    };
    if chrom.split('\t').count() != 10 {
        refuse(format!("ReblockGVCF can take multiple input GVCFs, but they must be non-overlapping shards from the same sample. {input} has {} samples", chrom.split('\t').count().saturating_sub(9)));
    }
    let mut allele_based: Vec<String> = Vec::new();
    let mut set: Vec<String> = Vec::new();
    for l in &hdr {
        if l.starts_with("##fileformat=") {
            continue;
        }
        if let Some(r) = l.strip_prefix("##FORMAT=<ID=") {
            let id = r.split(',').next().unwrap();
            let num = r.split(",Number=").nth(1).and_then(|x| x.split(',').next()).unwrap_or("");
            if ["A", "R", "G"].contains(&num) {
                allele_based.push(id.to_string());
            }
        }
        if l.starts_with("##GVCFBlock") {
            continue;
        }
        if let Some(r) = l.strip_prefix("##INFO=<ID=") {
            let id = r.split(',').next().unwrap();
            if INFO_TO_REMOVE.contains(&id) || id == "RAW_MQ" {
                if id == "RAW_MQ" {
                    refuse("the deprecated RAW_MQ annotation: not supported".into());
                }
                continue;
            }
        }
        set.push(l.clone());
    }
    for l in HEADER_LINES.lines() {
        set.push(l.to_string());
    }
    for &(lo, hi) in &bands_of(&gqb) {
        set.push(format!("##GVCFBlock{lo}-{hi}=minGQ={lo}(inclusive),maxGQ={hi}(exclusive)"));
    }
    let cmdline = argv.join(" ");
    set.push(format!(
        "##variant-rsCommandLine=<ID=ReblockGVCF,CommandLine=\"reblockgvcf {}\",Version=\"{}\",Date=\"{}\">",
        cmdline.replace('"', "\\\""),
        env!("CARGO_PKG_VERSION"),
        crate::now_string()
    ));
    let key = |l: &String| -> String {
        if l.starts_with("##contig=<") {
            "##contig=<".into()
        } else if let Some(rest) = l.strip_prefix("##variant-rsCommandLine") {
            format!("##GATKCommandLine{rest}")
        } else {
            l.clone()
        }
    };
    // a set: identical lines once; contig lines keep their order
    let mut seen = std::collections::HashSet::new();
    set.retain(|l| seen.insert(l.clone()));
    set.sort_by_key(|l| key(l));
    let mut header = String::from("##fileformat=VCFv4.2\n");
    for l in &set {
        header.push_str(l);
        header.push('\n');
    }
    header.push_str(&chrom);
    header.push('\n');
    let gz = crate::output::block_compressed(&output);
    let mut out = OutFile::create(&output, &header, if gz { Some(2) } else { None }).unwrap_or_else(|e| panic!("{e}"));

    let mut tool = Tool {
        w: Combiner {
            bands: bands_of(&gqb),
            floor,
            fasta: &fasta,
            block: None,
            next_available: None,
            buffer: Vec::new(),
            vcf_output_end: 0,
            buffer_end: 0,
            contig: None,
            out: String::with_capacity(8 << 20),
            n_out: 0,
        },
        qual_approx,
        allele_based,
    };
    let mut n_in = 0u64;
    for l in lines {
        if l.is_empty() {
            continue;
        }
        n_in += 1;
        let vc = parse_record(&l);
        tool.regenotype(vc);
        if tool.w.out.len() >= 4 << 20 {
            out.add_text(tool.w.out.as_bytes()).unwrap();
            tool.w.out.clear();
        }
    }
    tool.w.finish();
    out.add_text(tool.w.out.as_bytes()).unwrap();
    out.finish().unwrap();
    if std::env::var_os("VARIANT_RS_REBLOCK_STATS").is_some() {
        for (n, s) in STAT_NAMES.iter().zip(&STATS) {
            eprintln!("  {:>9}  {n}", s.load(std::sync::atomic::Ordering::Relaxed));
        }
    }
    eprintln!("reblockgvcf: {n_in} records in, {} out ({:.1}s)", tool.w.n_out, t0.elapsed().as_secs_f64());
}
