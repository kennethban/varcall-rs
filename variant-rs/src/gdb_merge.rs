//! GenomicsDB 1.5.5 query / merge as GATK GenotypeGVCFs receives it from `-V gendb://ws` (`docs/specs/genomicsdb_merge.md`,
//! a port of the validated prototype `merge_proto.py`): the walker's query windows, `scan_and_operate` segmentation
//! (calls per sample, replacement of overlapped calls, 1-bp stepping while a deletion / MNV is active),
//! `BroadCombinedGVCFOperator` (spanning deletions reduced to `*`, REF / ALT merge, AD / PL / AS remapping, INFO
//! combine operations, FORMAT, the GQ0 / PL0 rule) and the htslib text htsjdk parses. Records are produced as VCF
//! lines (what htsjdk would write for the received VariantContexts) for `genotypegvcfs::parse_mrec`.

use crate::gdb_store::{Array, Cell, CellStream, CellView, Num, Workspace};

/// `VariantWalkerBase.DEFAULT_DRIVING_VARIANTS_LOOKAHEAD_BASES`: shard size and feature-cache lookahead
pub const SHARD: i64 = 100_000;
const MAX_ALTS: usize = 50;
const MAX_GENOTYPE_COUNT: u64 = 1024;
const NON_REF: &str = "<NON_REF>";

fn refuse(why: String) -> ! {
    crate::gatk_cli::refuse(&why)
}

// ------------------------------------------------------------------------------------------------
// numbers as htslib prints them

/// C `%g` (6 significant digits)
fn fmt_g(d: f64) -> String {
    if d == 0.0 {
        return "0".into();
    }
    let e = d.abs().log10().floor() as i32;
    let s = format!("{:.5e}", d);
    let (m, ex) = s.split_once('e').unwrap();
    let ex: i32 = ex.parse().unwrap();
    let _ = e;
    if ex < -4 || ex >= 6 {
        let m = if m.contains('.') { m.trim_end_matches('0').trim_end_matches('.').to_string() } else { m.to_string() };
        format!("{m}e{}{:02}", if ex < 0 { '-' } else { '+' }, ex.abs())
    } else {
        let prec = (5 - ex).max(0) as usize;
        let t = format!("{:.*}", prec, d);
        if t.contains('.') { t.trim_end_matches('0').trim_end_matches('.').to_string() } else { t }
    }
}

/// htslib `kputd` (a float32 promoted to double)
pub fn kputd(d: f64) -> String {
    if d == 0.0 {
        return if d.is_sign_negative() { "-0".into() } else { "0".into() };
    }
    if d.is_nan() {
        return "nan".into();
    }
    let neg = d < 0.0;
    let d = d.abs();
    if !(0.0001..=999999.0).contains(&d) {
        return format!("{}{}", if neg { "-" } else { "" }, fmt_g(d));
    }
    let mut i = (d * 10000000000.0) as u64;
    i += if d < 0.001 {
        5
    } else if d < 0.01 {
        50
    } else if d < 0.1 {
        500
    } else if d < 1.0 {
        5000
    } else if d < 10.0 {
        50000
    } else if d < 100.0 {
        500000
    } else if d < 1000.0 {
        5000000
    } else if d < 10000.0 {
        50000000
    } else if d < 100000.0 {
        500000000
    } else {
        5000000000
    };
    let digits = i.to_string();
    let p = digits.len();
    let mut s = if p <= 10 {
        format!("0.{}{}", "0".repeat(10 - p), &digits[..6.min(p)])
    } else {
        let mut t = format!("{}.{}", &digits[..p - 10], &digits[p - 10..]);
        t.truncate(7);
        if t.ends_with('.') {
            t.pop();
        }
        t
    };
    if s.contains('.') {
        while s.ends_with('0') {
            s.pop();
        }
        if s.ends_with('.') {
            s.push('0');
        }
    }
    format!("{}{s}", if neg { "-" } else { "" })
}

fn fmt3(x: f64) -> String {
    format!("{:.3}", x)
}

// ------------------------------------------------------------------------------------------------
// vid map fields

#[derive(Clone, Debug)]
struct Field {
    name: String,
    op: String,
    remap_missing_with_non_ref: bool,
    two_d: bool,
    is_float: bool,
    /// first length descriptor is "A" (alt-only)
    a_length: bool,
}

fn load_fields(ws: &Workspace) -> (Vec<Field>, Vec<String>) {
    let v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(ws.dir.join("vidmap.json")).unwrap()).unwrap();
    let mut info = Vec::new();
    let mut format = Vec::new();
    for f in v["fields"].as_array().unwrap() {
        let name = f["name"].as_str().unwrap().to_string();
        let classes: Vec<&str> = f["vcf_field_class"].as_array().map(|a| a.iter().map(|x| x.as_str().unwrap()).collect()).unwrap_or_default();
        let lengths: Vec<String> = f["length"].as_array().map(|a| a.iter().map(|l| l["variable_length_descriptor"].as_str().unwrap_or("").to_string()).collect()).unwrap_or_default();
        let types: Vec<String> = f["type"].as_array().map(|a| a.iter().map(|x| x.as_str().unwrap().to_lowercase()).collect()).unwrap_or_default();
        let op = f["VCF_field_combine_operation"].as_str().map(|s| s.to_string()).unwrap_or_else(|| {
            match name.as_str() {
                "BaseQRankSum" | "ClippingRankSum" | "MQRankSum" | "ReadPosRankSum" | "MQ" | "MQ0" | "ExcessHet" => "median",
                "RAW_MQ" => "sum",
                "DP" => "DP",
                _ => "unknown",
            }
            .to_string()
        });
        let fld = Field {
            name: name.clone(),
            op,
            remap_missing_with_non_ref: !f["disable_remap_missing_with_non_ref"].as_bool().unwrap_or(false),
            two_d: lengths.len() == 2,
            is_float: types.first().is_some_and(|t| t == "float"),
            a_length: lengths.first().is_some_and(|l| l == "A"),
        };
        if classes.contains(&"INFO") && name != "END" && name != "DP" && ["sum", "median", "mean", "element_wise_sum", "histogram_sum", "concatenate"].contains(&fld.op.as_str()) {
            if ["mean", "concatenate"].contains(&fld.op.as_str()) {
                refuse(format!("GenomicsDB INFO field {name} with combine operation {}: not supported", fld.op));
            }
            info.push(fld.clone());
        }
        if classes.contains(&"FORMAT") {
            format.push(name);
        }
    }
    (info, format)
}

// ------------------------------------------------------------------------------------------------
// calls

/// An INFO value of a call, by field kind
#[derive(Clone, Debug)]
enum IV {
    None,
    Scalar(f64),
    Vec(Vec<Option<f64>>),
    TwoD(Vec<Vec<Option<f64>>>),
    Hist(Vec<Vec<(Option<f32>, Option<i64>)>>),
}

/// A stored record (begin cell) of one sample, positions 1-based within the contig.
#[derive(Clone)]
struct Rec {
    row: usize,
    begin: i64,
    end: i64,
    cell: Cell,
}

/// A record filled into a sample's slot (fresh state per query).
#[derive(Clone)]
struct Call {
    row: usize,
    begin: i64,
    end: i64,
    ref_: String,
    alts: Vec<String>,
    del_or_mnv: bool,
    info_invalid: bool,
    gt: Option<Vec<(Option<i64>, bool)>>,
    ad: Option<Vec<Option<i64>>>,
    pl: Option<Vec<Option<i64>>>,
    sb: Option<Vec<Option<i64>>>,
    gq: Option<i64>,
    min_dp: Option<i64>,
    ps: Option<i64>,
    dp_format: Option<i64>,
    dp_info: Option<i64>,
    pgt: Option<String>,
    pid: Option<String>,
    id: Option<String>,
    info: Vec<IV>,
}

fn is_symbolic(a: &str) -> bool {
    a == NON_REF || a == "*" || (!a.is_empty() && ((a.starts_with('<') && a.ends_with('>')) || a.contains('[') || a.contains(']')))
}

fn is_deletion(r: &str, a: &str) -> bool {
    r.len() > 1 && (a == "*" || (!is_symbolic(a) && a.len() < r.len()))
}

fn is_mnv(r: &str, a: &str) -> bool {
    r.len() > 1 && !is_symbolic(a) && a.len() == r.len() && r.bytes().zip(a.bytes()).filter(|(x, y)| x != y).count() > 1
}

/// VCF genotype index of an allele combination (any ploidy)
fn gt_index(alleles: &[i64]) -> usize {
    let mut a = alleles.to_vec();
    a.sort_unstable();
    let mut idx = 0u64;
    for (k, &x) in a.iter().enumerate() {
        idx += binom(x as u64 + k as u64, k as u64 + 1);
    }
    idx as usize
}

fn binom(n: u64, k: u64) -> u64 {
    if k > n {
        return 0;
    }
    let mut r: u64 = 1;
    for i in 0..k {
        r = r * (n - i) / (i + 1);
    }
    r
}

fn num_opt(n: &Num) -> Option<f64> {
    match n {
        Num::I(x) => Some(*x as f64),
        Num::F(x) => Some(*x as f64),
        Num::Missing => None,
    }
}

fn ints(v: Option<Vec<Num>>) -> Option<Vec<Option<i64>>> {
    v.map(|v| v.iter().map(|n| num_opt(n).map(|x| x as i64)).collect())
}

fn int1(v: Option<Vec<Num>>) -> Option<i64> {
    v.and_then(|v| v.first().and_then(num_opt)).map(|x| x as i64)
}

struct Ctx {
    info: Vec<Field>,
    has_fmt: Vec<(&'static str, bool)>,
    n: usize,
}

impl Ctx {
    fn has(&self, k: &str) -> bool {
        self.has_fmt.iter().any(|(x, h)| *x == k && *h)
    }
}

fn make_call(rec: &Rec, schema: &crate::gdb_store::Schema, ctx: &Ctx) -> Call {
    let v = CellView { schema, cell: &rec.cell };
    let ref_ = v.ref_allele();
    let alts = v.alts();
    let del_or_mnv = alts.iter().any(|a| is_deletion(&ref_, a) || is_mnv(&ref_, a));
    let gt = v.nums("GT").map(|x| {
        let raw: Vec<i64> = x.iter().map(|n| num_opt(n).map_or(-1, |y| y as i64)).collect();
        let mut out = vec![(if raw[0] < 0 { None } else { Some(raw[0]) }, false)];
        let mut k = 1;
        while k + 1 < raw.len() {
            out.push((if raw[k + 1] < 0 { None } else { Some(raw[k + 1]) }, raw[k] == 1));
            k += 2;
        }
        out
    });
    let info = ctx
        .info
        .iter()
        .map(|f| {
            if f.op == "histogram_sum" {
                let b = v.multid(&format!("{}_tuple_element_0", f.name), true);
                let c = v.multid(&format!("{}_tuple_element_1", f.name), false);
                match (b, c) {
                    (Some(b), Some(c)) => IV::Hist(
                        b.iter()
                            .zip(&c)
                            .map(|(bs, cs)| {
                                (0..bs.len().max(cs.len()))
                                    .map(|i| {
                                        let bin = bs.get(i).and_then(|n| match n {
                                            Num::F(x) if !x.is_nan() => Some(*x),
                                            _ => None,
                                        });
                                        let cnt = cs.get(i).and_then(num_opt).map(|x| x as i64);
                                        (bin, cnt)
                                    })
                                    .collect()
                            })
                            .collect(),
                    ),
                    _ => IV::None,
                }
            } else if f.two_d {
                match v.multid(&f.name, f.is_float) {
                    Some(m) => IV::TwoD(m.iter().map(|e| e.iter().map(num_opt).collect()).collect()),
                    None => IV::None,
                }
            } else if f.op == "element_wise_sum" {
                match v.nums(&f.name) {
                    Some(x) => IV::Vec(x.iter().map(num_opt).collect()),
                    None => IV::None,
                }
            } else {
                match v.nums(&f.name).and_then(|x| x.first().and_then(num_opt)) {
                    Some(x) => IV::Scalar(x),
                    None => IV::None,
                }
            }
        })
        .collect();
    Call {
        row: rec.row,
        begin: rec.begin,
        end: rec.end,
        del_or_mnv,
        info_invalid: false,
        gt,
        ad: if ctx.has("AD") { ints(v.nums("AD")) } else { None },
        pl: if ctx.has("PL") { ints(v.nums("PL")) } else { None },
        sb: if ctx.has("SB") { ints(v.nums("SB")) } else { None },
        gq: if ctx.has("GQ") { int1(v.nums("GQ")) } else { None },
        min_dp: if ctx.has("MIN_DP") { int1(v.nums("MIN_DP")) } else { None },
        ps: if ctx.has("PS") { int1(v.nums("PS")) } else { None },
        dp_format: int1(v.nums("DP_FORMAT")),
        dp_info: int1(v.nums("DP")),
        pgt: if ctx.has("PGT") { v.text("PGT") } else { None },
        pid: if ctx.has("PID") { v.text("PID") } else { None },
        id: v.text("ID"),
        ref_,
        alts,
        info,
    }
}

// ------------------------------------------------------------------------------------------------
// remapping (merged -> input LUT `inv`, `input_nr` = the sample's NON_REF index)

type Inv = std::collections::HashMap<usize, usize>;

fn remap_alleles(data: &[Option<i64>], inv: &Inv, nall: usize, input_nr: Option<usize>, alt_only: bool) -> Vec<Option<i64>> {
    let rng = if alt_only { 1..nall } else { 0..nall };
    rng.map(|j| {
        let ij = match inv.get(&j) {
            Some(&x) => x,
            None => match input_nr {
                Some(x) => x,
                None => return None,
            },
        };
        let ix = if alt_only { ij - 1 } else { ij };
        data.get(ix).copied().flatten()
    })
    .collect()
}

fn remap_genotypes(data: &[Option<i64>], inv: &Inv, nall: usize, input_nr: Option<usize>, ploidy: usize) -> Option<Vec<Option<i64>>> {
    let get = |j: usize| inv.get(&j).copied().or(input_nr);
    match ploidy {
        0 => None,
        1 => Some((0..nall).map(|j| get(j).and_then(|ij| data.get(ij).copied().flatten())).collect()),
        2 => {
            let mut out = vec![None; nall * (nall + 1) / 2];
            for j in 0..nall {
                for k in j..nall {
                    let gi = k * (k + 1) / 2 + j;
                    let (Some(ij), Some(ik)) = (get(j), get(k)) else { continue };
                    let (a, b) = if ij <= ik { (ij, ik) } else { (ik, ij) };
                    let igi = b * (b + 1) / 2 + a;
                    out[gi] = data.get(igi).copied().flatten();
                }
            }
            Some(out)
        }
        p => refuse(format!("GenomicsDB: PL remapping for ploidy {p}: not supported")),
    }
}

fn remap_2d<T: Clone>(data: &[Vec<T>], inv: &Inv, nall: usize, input_nr: Option<usize>, alt_only: bool) -> Vec<Vec<T>> {
    let rng = if alt_only { 1..nall } else { 0..nall };
    rng.map(|j| {
        let ij = match inv.get(&j) {
            Some(&x) => x,
            None => match input_nr {
                Some(x) => x,
                None => return Vec::new(),
            },
        };
        let ix = if alt_only { ij - 1 } else { ij };
        data.get(ix).cloned().unwrap_or_default()
    })
    .collect()
}

// ------------------------------------------------------------------------------------------------
// BroadCombinedGVCFOperator

/// `handle_deletions`: spanning deletions / MNVs reduced to REF `N`, ALT `*`(, NON_REF)
fn handle_deletions(s: i64, calls: &mut [Option<Call>], ctx: &Ctx) {
    for c in calls.iter_mut().flatten() {
        if !c.del_or_mnv || s <= c.begin {
            continue;
        }
        c.info_invalid = true;
        if c.alts[0] == "*" && (c.alts.len() == 1 || (c.alts.len() == 2 && c.alts[1] == NON_REF)) {
            continue;
        }
        let ploidy = c.gt.as_ref().map_or(0, |g| g.len());
        let mut lowest: i64 = -1;
        let mut lowest_pl = i64::MAX;
        let mut non_ref: Option<usize> = None;
        for (i, a) in c.alts.iter().enumerate() {
            let ai = i + 1;
            if is_deletion(&c.ref_, a) || is_mnv(&c.ref_, a) {
                if lowest < 0 {
                    lowest = ai as i64;
                }
                if let Some(pl) = &c.pl {
                    let gi = gt_index(&vec![ai as i64; ploidy]);
                    if let Some(Some(v)) = pl.get(gi) {
                        if *v < lowest_pl {
                            lowest_pl = *v;
                            lowest = ai as i64;
                        }
                    }
                }
            } else if a == NON_REF {
                non_ref = Some(ai);
            }
        }
        let lowest = lowest as usize;
        let mut lut: std::collections::HashMap<usize, usize> = [(0, 0), (lowest, 1)].into_iter().collect();
        let mut inv: Inv = [(0, 0), (1, lowest)].into_iter().collect();
        if let Some(nr) = non_ref {
            lut.insert(nr, 2);
            inv.insert(2, nr);
        }
        let nred = if non_ref.is_some() { 3 } else { 2 };
        if let Some(ad) = &c.ad {
            c.ad = Some(remap_alleles(ad, &inv, nred, non_ref, false));
        }
        if let Some(pl) = &c.pl {
            c.pl = remap_genotypes(pl, &inv, nred, non_ref, ploidy);
        }
        for (k, f) in ctx.info.iter().enumerate() {
            if !f.two_d {
                continue;
            }
            c.info[k] = match &c.info[k] {
                IV::TwoD(d) => IV::TwoD(remap_2d(d, &inv, nred, non_ref, f.a_length)),
                IV::Hist(d) => IV::Hist(remap_2d(d, &inv, nred, non_ref, f.a_length)),
                x => x.clone(),
            };
        }
        if let Some(gt) = &c.gt {
            c.gt = Some(
                gt.iter()
                    .map(|&(a, ph)| (a.map(|a| lut.get(&(a as usize)).map(|&x| x as i64).or(if non_ref.is_some() { Some(2) } else { None })).flatten(), ph))
                    .collect(),
            );
        }
        c.ref_ = "N".into();
        c.alts = if non_ref.is_some() { vec!["*".into(), NON_REF.into()] } else { vec!["*".into()] };
    }
}

/// One emitted record `[s, e]` as a VCF line (htsjdk normalized: INFO sorted, FORMAT GT then sorted).
fn operate(s: i64, e: i64, calls: &mut [Option<Call>], ctx: &Ctx, contig: &str, fasta: &crate::reference::Fasta) -> String {
    handle_deletions(s, calls, ctx);
    // merge_reference_allele
    let mut mref = String::from("N");
    for c in calls.iter().flatten() {
        if c.begin < s {
            continue;
        }
        let r = &c.ref_;
        if r.len() > mref.len() {
            mref = if mref.starts_with('N') { r.clone() } else { format!("{mref}{}", &r[mref.len()..]) };
        } else if mref.starts_with('N') && !r.starts_with('N') {
            mref = r.clone();
        }
    }
    // merge_alt_alleles
    let mut malts: Vec<String> = Vec::new();
    let mut luts: Vec<Option<std::collections::HashMap<usize, usize>>> = vec![None; ctx.n];
    let mut nonref_in: Vec<Option<usize>> = vec![None; ctx.n];
    let mut nr_exists = false;
    for c in calls.iter().flatten() {
        let mut lut: std::collections::HashMap<usize, usize> = [(0, 0)].into_iter().collect();
        let suffix = if c.ref_.len() < mref.len() { mref[c.ref_.len()..].to_string() } else { String::new() };
        for (i, a) in c.alts.iter().enumerate() {
            let ai = i + 1;
            if a == NON_REF {
                nonref_in[c.row] = Some(ai);
                nr_exists = true;
                continue;
            }
            let a2 = if !suffix.is_empty() && !is_symbolic(a) { format!("{a}{suffix}") } else { a.clone() };
            let idx = match malts.iter().position(|x| *x == a2) {
                Some(p) => p + 1,
                None => {
                    malts.push(a2);
                    malts.len()
                }
            };
            lut.insert(ai, idx);
        }
        luts[c.row] = Some(lut);
    }
    if nr_exists {
        malts.push(NON_REF.into());
        let nri = malts.len();
        for c in calls.iter().flatten() {
            if let Some(nr) = nonref_in[c.row] {
                luts[c.row].as_mut().unwrap().insert(nr, nri);
            }
        }
    }
    let nall = malts.len() + 1;
    let ref_block_only = mref.len() == 1 && malts.len() == 1 && malts[0] == NON_REF;
    let remap = !ref_block_only;
    let too_many = malts.len() > MAX_ALTS;
    if mref == "N" {
        let (b, _) = fasta.fetch(contig, s, s);
        mref = if !b.is_empty() && b"ACGT".contains(&b[0]) { (b[0] as char).to_string() } else { "N".into() };
    }
    // inverse LUTs (merged -> input; with duplicates the last added pair wins, LUT iteration order: the input
    // index ascending as the prototype's dict insertion order, NON_REF last)
    let mut invs: Vec<Inv> = vec![Inv::new(); ctx.n];
    for c in calls.iter().flatten() {
        let lut = luts[c.row].as_ref().unwrap();
        let mut keys: Vec<usize> = lut.keys().copied().filter(|k| Some(*k) != nonref_in[c.row]).collect();
        keys.sort_unstable();
        if let Some(nr) = nonref_in[c.row] {
            if lut.contains_key(&nr) {
                keys.push(nr);
            }
        }
        let mut d = Inv::new();
        for k in keys {
            d.insert(lut[&k], k);
        }
        invs[c.row] = d;
    }
    // remapped GT / AD / PL per call
    struct R {
        gt: Option<Vec<(Option<i64>, bool)>>,
        ad: Option<Vec<Option<i64>>>,
        pl: Option<Vec<Option<i64>>>,
    }
    let mut rs: Vec<Option<R>> = (0..ctx.n).map(|_| None).collect();
    for c in calls.iter().flatten() {
        let lut = luts[c.row].as_ref().unwrap();
        let nr_in = nonref_in[c.row];
        let gt = match (&c.gt, remap) {
            (Some(g), true) => Some(
                g.iter()
                    .map(|&(a, ph)| (a.and_then(|a| lut.get(&(a as usize)).map(|&x| x as i64).or(if nr_exists { Some(nall as i64 - 1) } else { None })), ph))
                    .collect(),
            ),
            (g, _) => g.clone(),
        };
        let ploidy = c.gt.as_ref().map_or(0, |g| g.len());
        let nr_arg = if nr_exists { nr_in } else { None };
        let ad = match &c.ad {
            Some(v) if remap => Some(remap_alleles(v, &invs[c.row], nall, nr_arg, false)),
            v => v.clone(),
        };
        let pl = match &c.pl {
            Some(v) if remap => {
                if too_many || binom((nall - 1 + ploidy) as u64, ploidy as u64) > MAX_GENOTYPE_COUNT {
                    None
                } else {
                    remap_genotypes(v, &invs[c.row], nall, nr_arg, ploidy)
                }
            }
            v => v.clone(),
        };
        rs[c.row] = Some(R { gt, ad, pl });
    }
    // ---------------- INFO
    let mut info: Vec<(String, String)> = Vec::new();
    if e > s {
        info.push(("END".into(), e.to_string()));
    }
    for (k, f) in ctx.info.iter().enumerate() {
        let vals: Vec<(&Call, &IV)> = calls.iter().flatten().filter(|c| !c.info_invalid && !matches!(c.info[k], IV::None)).map(|c| (c, &c.info[k])).collect();
        if vals.is_empty() {
            continue;
        }
        match f.op.as_str() {
            "sum" => {
                let good: Vec<f64> = vals.iter().filter_map(|(_, v)| if let IV::Scalar(x) = v { Some(*x) } else { None }).collect();
                if !good.is_empty() {
                    let t: f64 = good.iter().sum();
                    info.push((f.name.clone(), if f.is_float { kputd(t as f32 as f64) } else { (t as i64).to_string() }));
                }
            }
            "median" => {
                let mut good: Vec<f64> = vals.iter().filter_map(|(_, v)| if let IV::Scalar(x) = v { Some(*x) } else { None }).collect();
                if !good.is_empty() {
                    good.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    info.push((f.name.clone(), kputd(good[good.len() / 2])));
                }
            }
            "element_wise_sum" if !f.two_d => {
                let mut res: Vec<Option<f64>> = Vec::new();
                for (_, v) in &vals {
                    let IV::Vec(v) = v else { continue };
                    for (i, x) in v.iter().enumerate() {
                        if i >= res.len() {
                            res.push(*x);
                        } else if let Some(x) = x {
                            res[i] = Some(res[i].map_or(*x, |r| r + x));
                        }
                    }
                }
                let t: Vec<String> = res.iter().map(|x| match x {
                    None => ".".to_string(),
                    Some(x) => if f.is_float { kputd(*x as f32 as f64) } else { (*x as i64).to_string() },
                }).collect();
                info.push((f.name.clone(), t.join(",")));
            }
            "element_wise_sum" => {
                let mut res: Vec<Vec<Option<f64>>> = Vec::new();
                let mut any_valid = false;
                for (c, v) in &vals {
                    let IV::TwoD(v) = v else { continue };
                    let v = if remap {
                        let nr = if nr_exists && f.remap_missing_with_non_ref { nonref_in[c.row] } else { None };
                        remap_2d(v, &invs[c.row], nall, nr, f.a_length)
                    } else {
                        v.clone()
                    };
                    while res.len() < v.len() {
                        res.push(Vec::new());
                    }
                    for (i, ent) in v.iter().enumerate() {
                        let tgt = &mut res[i];
                        while tgt.len() < ent.len() {
                            tgt.push(None);
                        }
                        for (j, x) in ent.iter().enumerate() {
                            let Some(x) = x else { continue };
                            any_valid = true;
                            tgt[j] = Some(match tgt[j] {
                                None => *x,
                                Some(t) => if f.is_float { (t as f32 + *x as f32) as f64 } else { t + x },
                            });
                        }
                    }
                }
                if any_valid {
                    let t: Vec<String> = res
                        .iter()
                        .map(|ent| ent.iter().map(|x| match x {
                            None => String::new(),
                            Some(x) => if f.is_float { fmt3(*x) } else { (*x as i64).to_string() },
                        }).collect::<Vec<_>>().join(","))
                        .collect();
                    info.push((f.name.clone(), t.join("|")));
                }
            }
            "histogram_sum" => {
                let mut res: Vec<Vec<(f32, i64)>> = Vec::new();
                for (c, v) in &vals {
                    let IV::Hist(v) = v else { continue };
                    let v = if remap {
                        let nr = if nr_exists && f.remap_missing_with_non_ref { nonref_in[c.row] } else { None };
                        remap_2d(v, &invs[c.row], nall, nr, false)
                    } else {
                        v.clone()
                    };
                    while res.len() < v.len() {
                        res.push(Vec::new());
                    }
                    for (i, ent) in v.iter().enumerate() {
                        for &(b, cnt) in ent {
                            let (Some(b), Some(cnt)) = (b, cnt) else { continue };
                            match res[i].iter_mut().find(|(x, _)| *x == b) {
                                Some(e) => e.1 += cnt,
                                None => res[i].push((b, cnt)),
                            }
                        }
                    }
                }
                let t: Vec<String> = res
                    .iter_mut()
                    .map(|h| {
                        h.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                        h.iter().map(|(b, c)| format!("{},{c}", fmt3(*b as f64))).collect::<Vec<_>>().join(",")
                    })
                    .collect();
                info.push((f.name.clone(), t.join("|")));
            }
            op => refuse(format!("GenomicsDB combine operation {op}: not supported")),
        }
    }
    // ---------------- FORMAT
    let n = ctx.n;
    let by: Vec<Option<&Call>> = {
        let mut v: Vec<Option<&Call>> = vec![None; n];
        for c in calls.iter().flatten() {
            v[c.row] = Some(c);
        }
        v
    };
    // name -> per sample text (None = missing)
    let mut fmt: Vec<(&str, Vec<Option<String>>)> = Vec::new();
    let list = |v: &Option<Vec<Option<i64>>>| v.as_ref().map(|x| x.iter().map(|y| y.map_or(".".to_string(), |z| z.to_string())).collect::<Vec<_>>().join(","));
    let mut gts: Option<Vec<Option<Vec<(Option<i64>, bool)>>>> = None;
    if (0..n).any(|i| by[i].is_some() && rs[i].as_ref().unwrap().gt.is_some()) {
        gts = Some((0..n).map(|i| if by[i].is_some() { rs[i].as_ref().unwrap().gt.clone() } else { None }).collect());
    }
    let mut gq_zero: Vec<usize> = Vec::new();
    let mut pl_zero: Vec<usize> = Vec::new();
    for k in ["AD", "GQ", "MIN_DP", "PGT", "PID", "PL", "PS", "SB"] {
        if !ctx.has(k) || (k == "PL" && too_many) {
            continue;
        }
        let col: Vec<Option<String>> = (0..n)
            .map(|i| {
                let c = by[i]?;
                let r = rs[i].as_ref().unwrap();
                match k {
                    "AD" => list(&r.ad),
                    "PL" => list(&r.pl),
                    "SB" => list(&c.sb),
                    "GQ" => c.gq.map(|x| x.to_string()),
                    "MIN_DP" => c.min_dp.map(|x| x.to_string()),
                    "PS" => c.ps.map(|x| x.to_string()),
                    "PGT" => c.pgt.clone(),
                    _ => c.pid.clone(),
                }
            })
            .collect();
        if col.iter().all(|v| v.is_none()) {
            continue;
        }
        if k == "GQ" {
            gq_zero = (0..n).filter(|&i| by[i].is_some_and(|c| c.gq == Some(0))).collect();
        }
        if k == "PL" {
            pl_zero = (0..n).filter(|&i| rs[i].as_ref().is_some_and(|r| r.pl.as_ref().is_some_and(|p| !p.is_empty() && p[0] == Some(0)))).collect();
        }
        fmt.push((k, col));
    }
    if let Some(g) = gts.as_mut() {
        for &i in &gq_zero {
            if pl_zero.contains(&i) {
                if let Some(x) = g[i].as_mut() {
                    *x = vec![(None, false); x.len()];
                }
            }
        }
    }
    // DP
    let dpf: Vec<Option<i64>> = (0..n).map(|i| by[i].and_then(|c| c.dp_format)).collect();
    let dpi: Vec<Option<i64>> = (0..n).map(|i| by[i].and_then(|c| c.dp_info)).collect();
    let has_min_dp = fmt.iter().any(|(k, _)| *k == "MIN_DP");
    let mdp: Vec<Option<i64>> = (0..n).map(|i| if has_min_dp { by[i].and_then(|c| c.min_dp) } else { None }).collect();
    if dpf.iter().any(|v| v.is_some()) || dpi.iter().any(|v| v.is_some()) {
        let mut tot = 0i64;
        for i in 0..n {
            if let Some(v) = dpi[i].or(mdp[i]).or(dpf[i]) {
                tot += v;
            }
        }
        if dpf.iter().any(|v| v.is_some()) {
            fmt.push(("DP", dpf.iter().map(|v| v.map(|x| x.to_string())).collect()));
        }
        if tot > 0 && !ref_block_only {
            info.push(("DP".into(), tot.to_string()));
        }
    }
    // ID: the active calls' IDs. GenomicsDB keeps them as an unordered set, so several come in no fixed order; sorted
    // here. No reader uses them: GenotypeGVCFs drops the input ID (its records are rebuilt), and
    // CreateSomaticPanelOfNormals reads only the alleles and AD.
    let mut ids: Vec<&str> = calls.iter().flatten().filter_map(|c| c.id.as_deref()).flat_map(|i| i.split(';')).collect();
    ids.sort_unstable();
    ids.dedup();
    let id = if ids.is_empty() { ".".to_string() } else { ids.join(";") };
    // line
    info.sort();
    let infos = if info.is_empty() { ".".to_string() } else { info.iter().map(|(k, v)| format!("{k}={v}")).collect::<Vec<_>>().join(";") };
    fmt.sort_by(|a, b| a.0.cmp(b.0));
    let mut line = format!("{contig}\t{s}\t{id}\t{mref}\t{}\t.\t.\t{infos}\tGT", if malts.is_empty() { ".".to_string() } else { malts.join(",") });
    for (k, _) in &fmt {
        line.push(':');
        line.push_str(k);
    }
    for i in 0..n {
        line.push('\t');
        let gt = gts.as_ref().and_then(|g| g[i].clone());
        match gt {
            None => line.push('.'),
            Some(g) => {
                let ph = g.iter().any(|(_, p)| *p);
                let t: Vec<String> = g.iter().map(|(a, _)| a.map_or(".".to_string(), |x| x.to_string())).collect();
                line.push_str(&t.join(if ph { "|" } else { "/" }));
            }
        }
        let mut vals: Vec<String> = fmt.iter().map(|(_, col)| col[i].clone().unwrap_or_else(|| ".".into())).collect();
        while vals.last().is_some_and(|v| v.chars().all(|c| c == '.' || c == ',')) {
            vals.pop();
        }
        for v in vals {
            line.push(':');
            line.push_str(&v);
        }
    }
    line
}

// ------------------------------------------------------------------------------------------------
// the scan (VariantQueryProcessor::scan_and_operate), streamed over the array

pub struct Scanner<'a> {
    arr: &'a Array,
    stream: CellStream<'a>,
    peeked: Option<Rec>,
    /// per row: its last record beginning before the current window
    prev: Vec<Option<Rec>>,
    ctx: Ctx,
    contig: String,
    offset: i64,
}

impl<'a> Scanner<'a> {
    pub fn new(ws: &Workspace, arr: &'a Array, contig: &str) -> Scanner<'a> {
        Self::with_fields(ws, arr, contig, false)
    }

    /// `other_format`: FORMAT fields other than GVCF ones (Mutect2's AF, F1R2, ...) are left out of the records
    /// instead of refused (for callers that read only the alleles, GT and AD).
    pub fn with_fields(ws: &Workspace, arr: &'a Array, contig: &str, other_format: bool) -> Scanner<'a> {
        let (info, format) = load_fields(ws);
        for f in &format {
            if !other_format && !["GT", "AD", "DP", "GQ", "MIN_DP", "PGT", "PID", "PL", "PS", "SB"].contains(&f.as_str()) {
                refuse(format!("GenomicsDB FORMAT field {f}: not supported"));
            }
        }
        let has_fmt = ["AD", "PL", "SB", "GQ", "MIN_DP", "PS", "PGT", "PID", "GT"].iter().map(|k| (*k, format.iter().any(|f| f == k))).collect();
        let offset = ws.contigs.iter().find(|c| c.0 == contig).unwrap().2;
        Scanner {
            arr,
            stream: CellStream::new(arr, i64::MIN, i64::MAX),
            peeked: None,
            prev: vec![None; ws.samples.len()],
            ctx: Ctx { info, has_fmt, n: ws.samples.len() },
            contig: contig.to_string(),
            offset,
        }
    }

    /// The next begin cell (END copies skipped), as a record with 1-based positions.
    fn next_rec(&mut self) -> Option<Rec> {
        if let Some(r) = self.peeked.take() {
            return Some(r);
        }
        for cell in self.stream.by_ref() {
            let end = CellView { schema: &self.arr.schema, cell: &cell }.end();
            if end < cell.col {
                continue; // END copy
            }
            return Some(Rec { row: cell.row as usize, begin: cell.col - self.offset + 1, end: end - self.offset + 1, cell });
        }
        None
    }

    /// The records of one query [a, b] (1-based positions), as VCF lines in order.
    pub fn scan(&mut self, a: i64, b: i64, fasta: &crate::reference::Fasta, out: &mut Vec<String>) {
        // records beginning before a: the latest per row
        let mut win: Vec<Rec> = Vec::new();
        while let Some(r) = self.next_rec() {
            if r.begin < a {
                let row = r.row;
                self.prev[row] = Some(r);
            } else if r.begin <= b {
                win.push(r);
            } else {
                self.peeked = Some(r);
                break;
            }
        }
        let n = self.ctx.n;
        let schema = &self.arr.schema;
        let mut calls: Vec<Option<Call>> = vec![None; n];
        // gt_get_column(a): the first stored cell >= a of each row (a begin cell at a, or the END copy of the previous
        // record, truncated to the next begin - 1 and absent when that is its own begin)
        for r in 0..n {
            let first = win.iter().find(|x| x.row == r);
            if let Some(f) = first {
                if f.begin == a {
                    calls[r] = Some(make_call(f, schema, &self.ctx));
                    continue;
                }
            }
            if let Some(p) = &self.prev[r] {
                let ce = match first {
                    Some(f) => p.end.min(f.begin - 1),
                    None => p.end,
                };
                if ce > p.begin && ce >= a {
                    let mut c = make_call(p, schema, &self.ctx);
                    c.end = ce;
                    calls[r] = Some(c);
                }
            }
        }
        let mut cur: Option<i64> = if calls.iter().any(|c| c.is_some()) { Some(a) } else { None };
        let mut ndel = calls.iter().flatten().filter(|c| c.del_or_mnv).count();
        let ctx = &self.ctx;
        let contig = self.contig.clone();
        let ranges = |calls: &mut Vec<Option<Call>>, mut cur: i64, nxt: i64, ndel: &mut usize, out: &mut Vec<String>| -> i64 {
            while cur < nxt {
                let top = match calls.iter().flatten().map(|c| c.end).min() {
                    Some(t) => t,
                    None => break,
                };
                let mut m = if top < nxt - 1 { top } else { nxt - 1 };
                if *ndel > 0 {
                    m = cur;
                }
                out.push(operate(cur, m, calls, ctx, &contig, fasta));
                for c in calls.iter_mut() {
                    if c.as_ref().is_some_and(|x| x.end == m) {
                        if c.as_ref().unwrap().del_or_mnv {
                            *ndel -= 1;
                        }
                        *c = None;
                    }
                }
                cur = m + 1;
            }
            cur
        };
        for r in win.iter().filter(|r| r.begin > a) {
            let cb = r.begin;
            let c0 = *cur.get_or_insert(cb);
            if cb != c0 {
                ranges(&mut calls, c0, cb, &mut ndel, out);
                cur = Some(cb);
            }
            if let Some(old) = &calls[r.row] {
                if old.end >= cb {
                    let ref_block = old.ref_.len() == 1 && old.alts.len() == 1 && old.alts[0] == NON_REF;
                    if !old.del_or_mnv && !ref_block {
                        refuse(format!("GenomicsDB: Unhandled overlapping variants at {}:{cb}", self.contig));
                    }
                    if old.del_or_mnv {
                        ndel -= 1;
                    }
                }
            }
            let c = make_call(r, schema, ctx);
            if c.del_or_mnv {
                ndel += 1;
            }
            calls[r.row] = Some(c);
        }
        if let Some(c0) = cur {
            ranges(&mut calls, c0, b + 1, &mut ndel, out);
        }
        // the window's records become the rows' previous records
        for r in win {
            let row = r.row;
            self.prev[row] = Some(r);
        }
    }
}

/// The received VariantContexts of one traversal interval [S, E] on `contig` (VCF lines), window by window.
pub fn interval_records(ws: &Workspace, arrays: &mut [(i64, i64, Scanner)], contig: &str, s: i64, e: i64, fasta: &crate::reference::Fasta, mut each: impl FnMut(&str)) {
    let _ = ws;
    let mut a = s;
    let mut out = Vec::new();
    while a <= e {
        let shard_end = (a + SHARD - 1).min(e);
        let win_end = shard_end + SHARD;
        out.clear();
        for (lo, hi, sc) in arrays.iter_mut() {
            let (qa, qb) = (a.max(*lo), win_end.min(*hi));
            if qa <= qb {
                sc.scan(qa, qb, fasta, &mut out);
            }
        }
        for line in &out {
            let pos: i64 = line.split('\t').nth(1).unwrap().parse().unwrap();
            if pos <= e {
                each(line);
            }
        }
        let _ = contig;
        a += 2 * SHARD;
    }
}

/// The records a walker receives from a whole workspace (no -L): `GenomicsDBFeatureReader.iterator()` queries every
/// array over its whole range [start, end], arrays in `ChrArrayFolderComparator` order (the contig name as
/// `String.compareTo`, then the start).
pub fn whole_workspace_records(ws: &Workspace, fasta: &crate::reference::Fasta, other_format: bool, mut each: impl FnMut(&str)) {
    let mut arrs: Vec<(String, i64, i64, &std::path::PathBuf)> = Vec::new();
    for p in &ws.arrays {
        let name = p.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
        let f: Vec<&str> = name.split('$').collect();
        match (f.len(), f.get(1).and_then(|x| x.parse::<i64>().ok()), f.get(2).and_then(|x| x.parse::<i64>().ok())) {
            (3, Some(lo), Some(hi)) => arrs.push((f[0].to_string(), lo, hi, p)),
            _ => refuse(format!("-V gendb:// without -L: array {name} is not named contig$start$end; not supported")),
        }
    }
    arrs.sort_by(|a, b| a.0.encode_utf16().cmp(b.0.encode_utf16()).then(a.1.cmp(&b.1)));
    let mut out = Vec::new();
    for (contig, lo, hi, p) in &arrs {
        let a = crate::gdb_store::open_array(p);
        let mut sc = Scanner::with_fields(ws, &a, contig, other_format);
        out.clear();
        sc.scan(*lo, *hi, fasta, &mut out);
        for l in &out {
            each(l);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn htslib_kputd() {
        for (x, s) in [(2.0, "2.0"), (0.0, "0"), (1.383, "1.383"), (-1.0, "-1.0"), (0.5, "0.5"), (60.0, "60.0"), (-0.431, "-0.431"), (36.45, "36.45")] {
            assert_eq!(kputd(x as f32 as f64), s, "{x}");
        }
    }
}
