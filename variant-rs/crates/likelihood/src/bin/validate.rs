//! Validates hc_likelihood against GATK dumps.
//!
//! usage: validate <gatk_assembly.jsonl> <gatk_genotyping.jsonl> [--batch N] [--limit N] [--detail N] [--threads N] [--mates mates.jsonl] [--cpu avx512|avx2]
//!
//! Joins the two dumps by region start, runs prepare (events, trimming, haplotype trimming,
//! genotyping reads, PairHMM inputs) per region, batches the PairHMM of `--batch` regions per GPU
//! launch, runs finish, and compares every field with GATK.

use hc_likelihood::{cigar_string, Haplotype, Params, Prepared, Read, RegionInput, RegionResult, Span, StageTimes};
use rayon::prelude::*;
use serde::de::{Deserializer, Visitor};
use serde::Deserialize;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::time::{Duration, Instant};

#[derive(Deserialize, Clone)]
struct JRead {
    name: String,
    flag: u16,
    start: i64,
    end: i64,
    cigar: String,
    seq: String,
    qual: String,
    mapq: u8,
    mate_start: i64,
    tlen: i64,
}

#[derive(Deserialize, Clone)]
struct JHap {
    bases: String,
    cigar: String,
    align_start_wrt_ref: i64,
    loc_start: i64,
    loc_end: i64,
    is_ref: bool,
    score: Option<f64>,
    kmer: u32,
}

#[derive(Deserialize)]
struct JAsm {
    start: i64,
    end: i64,
    padded_start: i64,
    padded_end: i64,
    final_reads: Vec<JRead>,
    ref_padded_start: i64,
    ref_padded: String,
    haplotypes: Vec<JHap>,
}

#[derive(Deserialize)]
struct JHmmRead {
    bases: String,
    qual: Vec<i16>,
    ins: Vec<i16>,
    del: Vec<i16>,
    gcp: Vec<i16>,
}

/// a JSON number, or "NaN"/"Infinity"/"-Infinity"
#[derive(Clone, Copy)]
struct JNum(f64);
impl<'de> Deserialize<'de> for JNum {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct V;
        impl<'de> Visitor<'de> for V {
            type Value = JNum;
            fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                f.write_str("number or NaN/Infinity string")
            }
            fn visit_f64<E>(self, v: f64) -> Result<JNum, E> {
                Ok(JNum(v))
            }
            fn visit_i64<E>(self, v: i64) -> Result<JNum, E> {
                Ok(JNum(v as f64))
            }
            fn visit_u64<E>(self, v: u64) -> Result<JNum, E> {
                Ok(JNum(v as f64))
            }
            fn visit_str<E: serde::de::Error>(self, v: &str) -> Result<JNum, E> {
                match v {
                    "NaN" => Ok(JNum(f64::NAN)),
                    "Infinity" => Ok(JNum(f64::INFINITY)),
                    "-Infinity" => Ok(JNum(f64::NEG_INFINITY)),
                    _ => Err(E::custom("bad number string")),
                }
            }
        }
        d.deserialize_any(V)
    }
}

#[derive(Deserialize)]
struct JGeno {
    start: i64,
    exit: Option<String>,
    trim_variation_present: Option<bool>,
    variant_start: Option<i64>,
    variant_end: Option<i64>,
    variant_padded_start: Option<i64>,
    variant_padded_end: Option<i64>,
    geno_start: Option<i64>,
    geno_end: Option<i64>,
    geno_padded_start: Option<i64>,
    geno_padded_end: Option<i64>,
    assembly_variation_present: Option<bool>,
    genotyping_reads: Option<Vec<JRead>>,
    filtered_reads: Option<Vec<JRead>>,
    trimmed_haplotypes: Option<Vec<JHap>>,
    lk_haplotypes: Option<Vec<JHap>>,
    lk_reads: Option<Vec<JRead>>,
    lk_filtered_reads: Option<Vec<JRead>>,
    lk_matrix: Option<Vec<Vec<JNum>>>,
    hmm_haps: Option<Vec<String>>,
    hmm_reads: Option<Vec<JHmmRead>>,
}

#[derive(Deserialize)]
struct JMates {
    start: i64,
    mate_other_contig: Vec<bool>,
}

struct Region {
    span: Span,
    padded: Span,
    reads: Vec<Read>,
    haps: Vec<Haplotype>,
    ref_padded: Vec<u8>,
    ref_padded_start: i64,
}

fn to_read(j: &JRead) -> Read {
    Read {
        name: j.name.clone(),
        flag: j.flag,
        start: j.start,
        end: j.end,
        cigar: hc_assembly::cigar::parse(&j.cigar),
        seq: j.seq.as_bytes().to_vec(),
        qual: j.qual.bytes().map(|b| b - 33).collect(),
        mapq: j.mapq,
        mate_start: j.mate_start,
        tlen: j.tlen,
        ref_id: 0,
        mate_ref_id: 0,
        bi: Vec::new(),
        bd: Vec::new(),
        sample: 0,
        oa_other_contig: false,
        bam_rec: None,
    }
}

fn to_hap(j: &JHap) -> Haplotype {
    Haplotype {
        bases: j.bases.as_bytes().to_vec(),
        cigar: hc_assembly::cigar::parse(&j.cigar),
        align_start_wrt_ref: j.align_start_wrt_ref,
        loc_start: j.loc_start,
        loc_end: j.loc_end,
        is_ref: j.is_ref,
        score: j.score.unwrap_or(f64::NAN),
        kmer: j.kmer,
    }
}

fn read_eq(r: &Read, j: &JRead) -> bool {
    r.name == j.name
        && r.flag == j.flag
        && r.start == j.start
        && r.end == j.end
        && cigar_string(&r.cigar) == j.cigar
        && r.seq == j.seq.as_bytes()
        && r.qual.len() == j.qual.len()
        && r.qual.iter().zip(j.qual.bytes()).all(|(&a, b)| a + 33 == b)
        && r.mapq == j.mapq
        && r.mate_start == j.mate_start
        && r.tlen == j.tlen
}

fn reads_eq(r: &[Read], j: &[JRead]) -> bool {
    r.len() == j.len() && r.iter().zip(j).all(|(a, b)| read_eq(a, b))
}

fn hap_eq(h: &Haplotype, j: &JHap) -> bool {
    let score_ok = match j.score {
        None => !h.score.is_finite(),
        Some(x) => h.score == x,
    };
    h.bases == j.bases.as_bytes()
        && cigar_string(&h.cigar) == j.cigar
        && h.align_start_wrt_ref == j.align_start_wrt_ref
        && h.loc_start == j.loc_start
        && h.loc_end == j.loc_end
        && h.is_ref == j.is_ref
        && score_ok
        && h.kmer == j.kmer
}

fn haps_eq(h: &[Haplotype], j: &[JHap]) -> bool {
    h.len() == j.len() && h.iter().zip(j).all(|(a, b)| hap_eq(a, b))
}

fn read_lines(path: &str, keep: impl Fn(&str) -> bool) -> Vec<String> {
    let f = std::fs::File::open(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let rd = BufReader::with_capacity(1 << 24, f);
    let mut v = Vec::new();
    for l in rd.lines() {
        let l = l.unwrap();
        if keep(&l) {
            v.push(l);
        }
    }
    v
}

#[derive(Default)]
struct Counts {
    n: usize,
    c: Vec<(&'static str, usize, usize)>, // name, ok, total
}
impl Counts {
    fn add(&mut self, name: &'static str, ok: bool) -> bool {
        if let Some(e) = self.c.iter_mut().find(|e| e.0 == name) {
            e.1 += ok as usize;
            e.2 += 1;
        } else {
            self.c.push((name, ok as usize, 1));
        }
        ok
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: validate <gatk_assembly.jsonl> <gatk_genotyping.jsonl> [--batch N] [--limit N] [--detail N] [--threads N] [--mates mates.jsonl] [--cpu avx512|avx2]");
        std::process::exit(2);
    }
    let mut batch = 256usize;
    let mut limit = usize::MAX;
    let mut detail = 5usize;
    let mut threads = 8usize;
    let mut mates_path: Option<String> = None;
    let mut cpu: Option<hc_likelihood::cpu::Flavor> = None;
    let mut i = 3;
    while i < args.len() {
        let v = || args[i + 1].parse::<usize>().unwrap();
        match args[i].as_str() {
            "--batch" => batch = v(),
            "--limit" => limit = v(),
            "--detail" => detail = v(),
            "--threads" => threads = v(),
            "--mates" => mates_path = Some(args[i + 1].clone()),
            "--cpu" => {
                cpu = Some(match args[i + 1].as_str() {
                    "avx512" => hc_likelihood::cpu::Flavor::Avx512,
                    "avx2" => hc_likelihood::cpu::Flavor::Avx2,
                    x => panic!("--cpu avx512|avx2, got {x}"),
                })
            }
            a => panic!("unknown argument {a}"),
        }
        i += 2;
    }
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();

    // ---------------- load dumps ----------------
    let t0 = Instant::now();
    let glines = read_lines(&args[2], |l| l.contains("\"exit\":"));
    let glines: Vec<String> = glines.into_iter().take(limit).collect();
    let geno: Vec<JGeno> = glines.par_iter().map(|l| serde_json::from_str(l).expect("bad genotyping json")).collect();
    drop(glines);
    let starts: std::collections::HashSet<i64> = geno.iter().map(|g| g.start).collect();
    let alines = read_lines(&args[1], |l| {
        let s = &l[l.find("\"start\":").unwrap() + 8..];
        let e = s.find(',').unwrap();
        starts.contains(&s[..e].parse::<i64>().unwrap())
    });
    let asm: HashMap<i64, Region> = alines
        .par_iter()
        .map(|l| {
            let a: JAsm = serde_json::from_str(l).expect("bad assembly json");
            (
                a.start,
                Region {
                    span: Span { start: a.start, end: a.end },
                    padded: Span { start: a.padded_start, end: a.padded_end },
                    reads: a.final_reads.iter().map(to_read).collect(),
                    haps: a.haplotypes.iter().map(to_hap).collect(),
                    ref_padded: a.ref_padded.into_bytes(),
                    ref_padded_start: a.ref_padded_start,
                },
            )
        })
        .collect();
    drop(alines);
    let load = t0.elapsed();
    println!("loaded {} regions with an exit ({} assembly records) in {:.2} s", geno.len(), asm.len(), load.as_secs_f64());
    // optional per-read "mate on another contig" flags (synthetic fuzz data only)
    let mates: HashMap<i64, Vec<bool>> = match &mates_path {
        None => HashMap::new(),
        Some(path) => read_lines(path, |_| true)
            .iter()
            .map(|l| {
                let m: JMates = serde_json::from_str(l).unwrap();
                (m.start, m.mate_other_contig)
            })
            .collect(),
    };
    let regions: Vec<&Region> = geno.iter().map(|g| asm.get(&g.start).expect("region missing from the assembly dump")).collect();
    let p = Params::default();
    let inputs: Vec<RegionInput> = regions
        .iter()
        .map(|r| RegionInput {
            span: r.span,
            padded_span: r.padded,
            reads: &r.reads,
            mate_on_other_contig: mates.get(&r.span.start).map(|v| &v[..]).unwrap_or(&[]),
            haplotypes: &r.haps,
            ref_padded: &r.ref_padded,
            ref_padded_start: r.ref_padded_start,
        })
        .collect();

    // ---------------- phase 1 ----------------
    let mut st = StageTimes::default();
    let t1 = Instant::now();
    let preps: Vec<Prepared> = inputs.iter().map(|inp| hc_likelihood::prepare_with(inp, &p, Some(&mut st))).collect();
    let prep_wall_1t = t1.elapsed();
    let t1p = Instant::now();
    let preps_par: Vec<Prepared> = inputs.par_iter().map(|inp| hc_likelihood::prepare_with(inp, &p, None)).collect();
    let prep_wall_par = t1p.elapsed();
    drop(preps_par);

    // ---------------- PairHMM (GPU, or the CPU fallback with --cpu) ----------------
    let mut gpu = Duration::ZERO;
    let mut raws: Vec<Vec<f64>> = Vec::with_capacity(preps.len());
    let mut launches = 0;
    let mut pairs = 0usize;
    let mut cells = 0f64;
    for p in preps.iter().filter(|p| p.needs_hmm()) {
        pairs += p.reads.len() * p.haplotypes.len();
        let rl: usize = p.reads.iter().map(|r| r.seq.len()).sum();
        let hl: usize = p.haplotypes.iter().map(|h| h.bases.len()).sum();
        cells += rl as f64 * hl as f64;
    }
    let hmm_label;
    if let Some(flavor) = cpu {
        let t = Instant::now();
        let chunks: Vec<Vec<Vec<f64>>> = preps
            .par_chunks(64)
            .map(|c| {
                let refs: Vec<&Prepared> = c.iter().collect();
                hc_likelihood::compute_likelihoods_cpu(&refs, flavor)
            })
            .collect();
        gpu = t.elapsed();
        raws = chunks.into_iter().flatten().collect();
        hmm_label = format!("PairHMM CPU {:?} ({} threads; {} pairs, {:.3e} cells)", flavor, threads, pairs, cells);
    } else {
        assert!(hc_cudahmm::init(false), "no usable CUDA device");
        // warm-up launch (context creation, table upload)
        {
            let first: Vec<&Prepared> = preps.iter().filter(|p| p.needs_hmm()).take(1).collect();
            let _ = hc_likelihood::compute_likelihoods_gpu(&first);
        }
        for chunk in preps.chunks(batch) {
            let refs: Vec<&Prepared> = chunk.iter().collect();
            let t = Instant::now();
            let r = hc_likelihood::compute_likelihoods_gpu(&refs);
            gpu += t.elapsed();
            launches += 1;
            raws.extend(r);
        }
        hmm_label = format!("PairHMM GPU ({} launches of <= {} regions; {} pairs, {:.3e} cells)", launches, batch, pairs, cells);
    }

    // PairHMM inputs vs GATK's batch (before finish consumes the prepared regions)
    let mut hmm_ok = 0;
    let mut hmm_total = 0;
    let mut hmm_shown = 0;
    for (g, pr) in geno.iter().zip(&preps) {
        if let (Some(hh), Some(hr)) = (&g.hmm_haps, &g.hmm_reads) {
            hmm_total += 1;
            let grp = pr.hmm_group();
            let ok = hh.len() == grp.haps.len()
                && hh.iter().zip(&grp.haps).all(|(a, b)| a.as_bytes() == *b)
                && hr.len() == grp.reads.len()
                && hr.iter().zip(&grp.reads).all(|(a, b)| {
                    let eq = |x: &[i16], y: &[u8]| x.len() == y.len() && x.iter().zip(y).all(|(&u, &v)| u == v as i16);
                    a.bases.as_bytes() == b.bases && eq(&a.qual, b.quals) && eq(&a.ins, b.ins) && eq(&a.del, b.del) && eq(&a.gcp, b.gcp)
                });
            hmm_ok += ok as usize;
            if !ok && hmm_shown < detail {
                hmm_shown += 1;
                println!("HMM INPUT MISMATCH region {}", g.start);
            }
        }
    }

    // ---------------- phase 2 ----------------
    let t2 = Instant::now();
    let results: Vec<RegionResult> = preps.into_iter().zip(raws).map(|(pr, raw)| hc_likelihood::finish_with(pr, &raw, &p, Some(&mut st))).collect();
    let finish_wall = t2.elapsed();

    // ---------------- code-path coverage of the real data ----------------
    {
        use hc_likelihood::{events, trim};
        let mut k: Vec<(&str, usize)> = Vec::new();
        let bump = |k: &mut Vec<(&'static str, usize)>, name: &'static str, n: usize| {
            if let Some(e) = k.iter_mut().find(|e| e.0 == name) { e.1 += n } else { k.push((name, n)) }
        };
        for (inp, g) in inputs.iter().zip(&geno) {
            let ev = events::variation_events(inp.haplotypes, inp.ref_padded, inp.ref_padded_start, 0);
            for e in &ev {
                bump(&mut k, if e.is_snp() { "events: SNP" } else if e.is_simple_insertion() { "events: insertion" } else if e.is_simple_deletion() { "events: deletion" } else { "events: other (MNP/complex)" }, 1);
                if e.is_indel() && e.ref_bases[0] != e.alt_bases[0] { bump(&mut k, "events: compound SNP+indel", 1); }
                if e.ref_bases.len() > 1 && e.alt_bases.len() > 1 && e.is_indel() { bump(&mut k, "events: compound deletion+insertion", 1); }
                if e.is_indel() && trim::overlaps(inp.span.start, inp.span.end, e.start, e.end()) {
                    let off = (inp.padded_span.start - inp.ref_padded_start) as usize;
                    let ctx = &inp.ref_padded[off..off + (inp.padded_span.end - inp.padded_span.start + 1) as usize];
                    match trim::event_tandem_repeat(ctx, inp.padded_span.start, e) {
                        Some((n, u)) => { bump(&mut k, "trim: indel with STR padding", 1); if n * u as i64 + 75 > 75 + 50 { bump(&mut k, "trim: STR padding > 125", 1); } }
                        None => bump(&mut k, "trim: indel with plain padding", 1),
                    }
                }
            }
            if let Some(vp) = g.variant_padded_start {
                let span = Span { start: vp, end: g.variant_padded_end.unwrap() };
                let mut seen = std::collections::HashSet::new();
                for h in inp.haplotypes {
                    match trim::trim_haplotype(h, span) {
                        None => bump(&mut k, "trimTo: haplotype dropped (ends in deletion / all insertion)", 1),
                        Some(t) => { if !seen.insert(t.bases.clone()) { bump(&mut k, if h.is_ref { "trimTo: reference merged with an earlier haplotype" } else { "trimTo: haplotype merged into an earlier one" }, 1); } }
                    }
                }
                let mut scratch = Vec::new();
                for r in inp.reads {
                    let mut c = r.clone();
                    let ok = hc_likelihood::clip::hard_clip_to_region(&mut c, span.start, g.variant_padded_end.unwrap(), &mut scratch);
                    if !ok { bump(&mut k, "reads: dropped by hardClipToRegion", 1); }
                    else if c.seq.len() != r.seq.len() { bump(&mut k, "reads: clipped", 1); if c.seq.len() < 10 { bump(&mut k, "reads: stub (< 10 bases) removed", 1); } }
                }
            }
        }
        println!("code paths exercised by this input:");
        for (n, c) in &k { println!("  {:<60} {}", n, c); }
        println!();
    }

    // ---------------- compare ----------------
    let mut c = Counts::default();
    c.c.push(("hmm_reads (bases, qual, ins, del, gcp; batch order)", hmm_ok, hmm_total));
    let mut max_diff = 0f64;
    let mut n_vals = 0usize;
    let mut n_gt_1e5 = 0usize;
    let mut n_gt_1e4 = 0usize;
    let mut n_exact = 0usize;
    let mut shown = 0;
    let mut near_threshold = 0usize;
    let (mut n_filtered, mut n_poor, mut n_evidence) = (0usize, 0usize, 0usize);
    for ((g, res), reg) in geno.iter().zip(&results).zip(&regions) {
        c.n += 1;
        let mut bad: Vec<&'static str> = Vec::new();
        let mut chk = |c: &mut Counts, name: &'static str, ok: bool| {
            if !c.add(name, ok) {
                bad.push(name);
            }
        };
        let gexit = g.exit.as_deref().unwrap();
        chk(&mut c, "exit", gexit == res.exit.as_str());
        chk(&mut c, "trim_variation_present", g.trim_variation_present == Some(res.variant_span.is_some()));
        if let Some(vs) = g.variant_start {
            let gv = (vs, g.variant_end.unwrap(), g.variant_padded_start.unwrap(), g.variant_padded_end.unwrap());
            let rv = res.variant_span.zip(res.variant_padded_span).map(|(a, b)| (a.start, a.end, b.start, b.end));
            chk(&mut c, "variant spans", Some(gv) == rv);
        }
        if let Some(gs) = g.geno_start {
            // regionForGenotyping == the trimmed variant region
            let gv = (gs, g.geno_end.unwrap(), g.geno_padded_start.unwrap(), g.geno_padded_end.unwrap());
            let rv = res.variant_span.zip(res.variant_padded_span).map(|(a, b)| (a.start, a.end, b.start, b.end));
            chk(&mut c, "geno spans", Some(gv) == rv);
            chk(&mut c, "assembly_variation_present", g.assembly_variation_present == Some(res.assembly_variation_present));
        }
        // genotyping reads = evidence + poorly modelled (in evidence order) when likelihoods were computed
        let (haps, all_reads): (&[Haplotype], Vec<Read>) = match &res.likelihoods {
            Some(lk) => {
                // reconstruct the pre-filter order: merge by original order is not stored, so compare via GATK's lists
                (&lk.haplotypes, Vec::new())
            }
            None => (&res.haplotypes, res.reads.clone()),
        };
        if let Some(th) = &g.trimmed_haplotypes {
            chk(&mut c, "trimmed haplotypes (all fields, order)", haps_eq(haps, th));
        }
        n_filtered += res.filtered_reads.len();
        if let Some(lk) = &res.likelihoods {
            n_poor += lk.filtered_reads.len();
            n_evidence += lk.reads.len();
        }
        if let Some(fr) = &g.filtered_reads {
            chk(&mut c, "filtered_reads (filterNonPassingReads)", reads_eq(&res.filtered_reads, fr));
        }
        if let (Some(gr), None) = (&g.genotyping_reads, &res.likelihoods) {
            chk(&mut c, "genotyping_reads", reads_eq(&all_reads, gr));
        }
        if let Some(lk) = &res.likelihoods {
            let gr = g.genotyping_reads.as_ref().unwrap();
            // genotyping reads: evidence + poorly modelled interleaved in original order
            let mut merged: Vec<&Read> = Vec::with_capacity(gr.len());
            {
                let (mut a, mut b) = (0, 0);
                for j in gr {
                    if a < lk.reads.len() && read_eq(&lk.reads[a], j) {
                        merged.push(&lk.reads[a]);
                        a += 1;
                    } else if b < lk.filtered_reads.len() && read_eq(&lk.filtered_reads[b], j) {
                        merged.push(&lk.filtered_reads[b]);
                        b += 1;
                    }
                }
                let ok = merged.len() == gr.len() && a == lk.reads.len() && b == lk.filtered_reads.len();
                chk(&mut c, "genotyping_reads", ok);
            }
            chk(&mut c, "lk_haplotypes", haps_eq(&lk.haplotypes, g.lk_haplotypes.as_ref().unwrap()));
            chk(&mut c, "lk_reads (evidence order)", reads_eq(&lk.reads, g.lk_reads.as_ref().unwrap()));
            chk(&mut c, "lk_filtered_reads", reads_eq(&lk.filtered_reads, g.lk_filtered_reads.as_ref().unwrap()));
            let m = g.lk_matrix.as_ref().unwrap();
            let shape_ok = m.len() == lk.haplotypes.len() && m.iter().all(|row| row.len() == lk.reads.len());
            chk(&mut c, "lk_matrix shape", shape_ok);
            if shape_ok {
                let mut region_max = 0f64;
                for (h, row) in m.iter().enumerate() {
                    for (r, x) in row.iter().enumerate() {
                        let d = (x.0 - lk.get(h, r)).abs();
                        let d = if x.0 == lk.get(h, r) { 0.0 } else { d };
                        region_max = region_max.max(d);
                        n_vals += 1;
                        n_gt_1e5 += (d > 1e-5) as usize;
                        n_gt_1e4 += (d > 1e-4) as usize;
                        n_exact += (x.0.to_bits() == lk.get(h, r).to_bits()) as usize;
                        if std::env::var("SHOWDIFF").is_ok() && x.0 != lk.get(h, r) { eprintln!("diff region {} h{} r{} gatk {:e} rust {:e} best {:e}", g.start, h, r, x.0, lk.get(h, r), (0..lk.haplotypes.len()).map(|hh| lk.get(hh, r)).fold(f64::NEG_INFINITY, f64::max)); }
                    }
                }
                max_diff = max_diff.max(region_max);
                chk(&mut c, "lk_matrix within 1e-4", region_max <= 1e-4);
            }
            // discrete outcomes near a threshold: reads whose best likelihood is within 1e-4 of the
            // filterPoorlyModeledEvidence threshold
            for (r, read) in lk.reads.iter().enumerate() {
                let best = (0..lk.haplotypes.len()).map(|h| lk.get(h, r)).fold(f64::NEG_INFINITY, f64::max);
                let thr = hc_likelihood::min_true_likelihood(&read.seq, &p); // static threshold: length only
                if (best - thr).abs() < 1e-4 {
                    near_threshold += 1;
                }
            }
        }
        if !bad.is_empty() && shown < detail {
            shown += 1;
            println!("MISMATCH region {} (exit {}): {:?}", g.start, gexit, bad);
            let _ = reg;
        }
    }


    println!("\nregions compared: {}", c.n);
    let mut by_exit: HashMap<&str, usize> = HashMap::new();
    for g in &geno {
        *by_exit.entry(g.exit.as_deref().unwrap()).or_default() += 1;
    }
    println!("GATK exits: {:?}", by_exit);
    for (name, ok, tot) in &c.c {
        println!("  {:<55} {:>6} / {:<6}{}", name, ok, tot, if ok == tot { "" } else { "   <-- MISMATCH" });
    }
    println!("\nreads: {} kept as evidence, {} poorly modelled (filterPoorlyModeledEvidence), {} removed by filterNonPassingReads", n_evidence, n_poor, n_filtered);
    println!("lk_matrix: {} values, max |diff| = {:.3e}, > 1e-5: {}, > 1e-4: {}, bit-identical: {}", n_vals, max_diff, n_gt_1e5, n_gt_1e4, n_exact);
    println!("reads kept whose best likelihood is within 1e-4 of the poorly-modelled threshold: {}", near_threshold);
    println!("\ntimings (single thread unless noted; JSON parsing excluded):");
    println!("  load+parse dumps ({} threads)           {:>8.3} s", threads, load.as_secs_f64());
    println!("  events (getVariationEvents)             {:>8.3} s", st.events.as_secs_f64());
    println!("  trim (AssemblyRegionTrimmer)            {:>8.3} s", st.trim.as_secs_f64());
    println!("  trimTo (haplotype trimming)             {:>8.3} s", st.trim_haplotypes.as_secs_f64());
    println!("  genotyping reads (clip/sort/filter)     {:>8.3} s", st.genotyping_reads.as_secs_f64());
    println!("  PairHMM inputs (PCR model, caps)        {:>8.3} s", st.hmm_inputs.as_secs_f64());
    println!("  prepare total, 1 thread (wall)          {:>8.3} s", prep_wall_1t.as_secs_f64());
    println!("  prepare total, {} threads (wall)         {:>8.3} s", threads, prep_wall_par.as_secs_f64());
    println!("  {} {:>8.3} s", hmm_label, gpu.as_secs_f64());
    println!("  finish (normalize + filter)             {:>8.3} s", st.finish.as_secs_f64());
    println!("  finish total (wall)                     {:>8.3} s", finish_wall.as_secs_f64());
    if cpu.is_none() {
        hc_cudahmm::report();
    }
}
