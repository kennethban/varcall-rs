//! Validate the Rust assembler against GATK's dumped assembly results.
//!
//! usage: validate <gatk_assembly.jsonl> [--limit N] [--detail N] [--only IDX] [--verbose]

use hc_assembly::{cigar, Assembler, Haplotype, Read, RegionInfo};
use serde::Deserialize;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::time::{Duration, Instant};

#[derive(Deserialize)]
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
struct JRegion {
    contig: String,
    start: i64,
    end: i64,
    padded_start: i64,
    padded_end: i64,
    final_reads: Vec<JRead>,
    ref_padded_start: i64,
    ref_padded_end: i64,
    ref_padded: String,
    ref_haplotype: JHap,
    haplotypes: Vec<JHap>,
    #[serde(default)]
    error: Option<String>,
}

type Key = (String, String, i64, i64, i64, bool);

fn jkey(h: &JHap) -> Key {
    (h.bases.clone(), h.cigar.clone(), h.align_start_wrt_ref, h.loc_start, h.loc_end, h.is_ref)
}
fn rkey(h: &Haplotype) -> Key {
    (
        String::from_utf8_lossy(&h.bases).into_owned(),
        cigar::to_string(&h.cigar),
        h.align_start_wrt_ref,
        h.loc_start,
        h.loc_end,
        h.is_ref,
    )
}

fn score_matches(g: Option<f64>, r: f64) -> bool {
    match g {
        None => !r.is_finite(),
        Some(x) => {
            if !r.is_finite() {
                return false;
            }
            if x == r {
                return true;
            }
            ((x - r).abs() / x.abs().max(r.abs())) <= 1e-9
        }
    }
}

fn fmt_score(s: f64) -> String {
    if s.is_finite() {
        format!("{}", s)
    } else {
        "null".to_string()
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args.get(1).expect("usage: validate <jsonl> [--limit N] [--detail N] [--only IDX] [--verbose]");
    let mut limit = usize::MAX;
    let mut detail = 5usize;
    let mut only: Option<usize> = None;
    let mut verbose = false;
    let mut i = 2;
    while i < args.len() {
        match args[i].as_str() {
            "--limit" => {
                limit = args[i + 1].parse().unwrap();
                i += 1;
            }
            "--detail" => {
                detail = args[i + 1].parse().unwrap();
                i += 1;
            }
            "--only" => {
                only = Some(args[i + 1].parse().unwrap());
                i += 1;
            }
            "--verbose" => verbose = true,
            _ => panic!("unknown arg {}", args[i]),
        }
        i += 1;
    }

    let f = std::fs::File::open(path).expect("open");
    let reader = BufReader::with_capacity(1 << 22, f);
    let mut assembler = Assembler::default();
    std::panic::set_hook(Box::new(|_| {}));
    let mut gatk_errors = 0usize;
    let mut both_errors = 0usize;
    let mut error_mismatch = 0usize;

    let mut n_regions = 0usize;
    let mut n_set_match = 0usize;
    let mut n_order_match = 0usize;
    let mut n_full_match = 0usize;
    let mut tot_gatk = 0usize;
    let mut tot_rust = 0usize;
    let mut missing = 0usize;
    let mut extra = 0usize;
    let mut kmer_mismatch = 0usize;
    let mut score_mismatch = 0usize;
    let mut ref_mismatch = 0usize;
    let mut n_detail = 0usize;
    let mut asm_time = Duration::ZERO;
    let mut reads_with_softclip = 0usize;
    let mut mismatching_regions: Vec<usize> = Vec::new();

    for (idx, line) in reader.lines().enumerate() {
        if idx >= limit {
            break;
        }
        if let Some(o) = only {
            if idx != o {
                continue;
            }
        }
        let line = line.expect("read line");
        let reg: JRegion = serde_json::from_str(&line).expect("parse json");
        let reads: Vec<Read> = reg
            .final_reads
            .iter()
            .map(|r| {
                if r.cigar.contains('S') {
                    reads_with_softclip += 1;
                }
                Read {
                    name: r.name.clone(),
                    flag: r.flag,
                    start: r.start,
                    end: r.end,
                    cigar: cigar::parse(&r.cigar),
                    seq: r.seq.as_bytes().to_vec(),
                    qual: r.qual.bytes().map(|b| b - 33).collect(),
                    mapq: r.mapq,
                    mate_start: r.mate_start,
                    tlen: r.tlen,
                    ref_id: 0,
                    mate_ref_id: 0,
                    bi: Vec::new(),
                    bd: Vec::new(),
                    sample: 0,
                    oa_other_contig: false,
                    bam_rec: None,
                }
            })
            .collect();
        let region = RegionInfo { start: reg.start, end: reg.end, padded_start: reg.padded_start, padded_end: reg.padded_end };
        let t0 = Instant::now();
        let res = assembler.try_assemble_region(&reads, region, reg.ref_padded.as_bytes(), reg.ref_padded_start, reg.ref_padded_end);
        asm_time += t0.elapsed();
        n_regions += 1;
        if reg.error.is_some() {
            gatk_errors += 1;
        }
        let res = match (res, &reg.error) {
            (Err(e), Some(ge)) => {
                both_errors += 1;
                if verbose {
                    println!("--- region #{}: both failed. GATK: {} / Rust: {}", idx, ge, e);
                }
                continue;
            }
            (Err(e), None) => {
                error_mismatch += 1;
                mismatching_regions.push(idx);
                println!("--- region #{}: Rust failed ({}) but GATK did not", idx, e);
                continue;
            }
            (Ok(_), Some(ge)) => {
                error_mismatch += 1;
                mismatching_regions.push(idx);
                println!("--- region #{}: GATK failed ({}) but Rust did not", idx, ge);
                continue;
            }
            (Ok(r), None) => r,
        };

        // reference haplotype
        let rref = &res.ref_haplotype;
        if jkey(&reg.ref_haplotype) != rkey(rref) || reg.ref_haplotype.kmer != rref.kmer {
            ref_mismatch += 1;
        }

        let g: Vec<Key> = reg.haplotypes.iter().map(jkey).collect();
        let r: Vec<Key> = res.haplotypes.iter().map(rkey).collect();
        tot_gatk += g.len();
        tot_rust += r.len();
        let mut gcount: HashMap<&Key, i64> = HashMap::new();
        for k in &g {
            *gcount.entry(k).or_insert(0) += 1;
        }
        let mut rcount: HashMap<&Key, i64> = HashMap::new();
        for k in &r {
            *rcount.entry(k).or_insert(0) += 1;
        }
        let mut miss = 0usize;
        let mut ext = 0usize;
        for (k, &c) in &gcount {
            let rc = *rcount.get(k).unwrap_or(&0);
            if c > rc {
                miss += (c - rc) as usize;
            }
        }
        for (k, &c) in &rcount {
            let gc = *gcount.get(k).unwrap_or(&0);
            if c > gc {
                ext += (c - gc) as usize;
            }
        }
        missing += miss;
        extra += ext;
        let set_match = miss == 0 && ext == 0;
        let order_match = g == r;
        // kmer/score on matched haplotypes
        let mut rmap: HashMap<&Key, &Haplotype> = HashMap::new();
        for (k, h) in r.iter().zip(res.haplotypes.iter()) {
            rmap.entry(k).or_insert(h);
        }
        let mut km = 0usize;
        let mut sm = 0usize;
        for (k, jh) in g.iter().zip(reg.haplotypes.iter()) {
            if let Some(rh) = rmap.get(k) {
                if rh.kmer != jh.kmer {
                    km += 1;
                }
                if !score_matches(jh.score, rh.score) {
                    sm += 1;
                }
            }
        }
        kmer_mismatch += km;
        score_mismatch += sm;
        if set_match {
            n_set_match += 1;
            if order_match {
                n_order_match += 1;
                if km == 0 && sm == 0 {
                    n_full_match += 1;
                }
            }
        }
        if !set_match || !order_match || km > 0 || sm > 0 {
            mismatching_regions.push(idx);
            if n_detail < detail || verbose {
                n_detail += 1;
                println!(
                    "--- region #{} {}:{}-{} (padded {}-{}) reads={} gatk={} rust={} missing={} extra={} order_match={} kmer_mm={} score_mm={}",
                    idx,
                    reg.contig,
                    reg.start,
                    reg.end,
                    reg.padded_start,
                    reg.padded_end,
                    reads.len(),
                    g.len(),
                    r.len(),
                    miss,
                    ext,
                    order_match,
                    km,
                    sm
                );
                let gk: HashMap<&Key, usize> = g.iter().enumerate().map(|(i, k)| (k, i)).collect();
                let rk: HashMap<&Key, usize> = r.iter().enumerate().map(|(i, k)| (k, i)).collect();
                for (i, (k, jh)) in g.iter().zip(reg.haplotypes.iter()).enumerate() {
                    let tag = match rk.get(k) {
                        Some(&ri) => format!("rust#{}", ri),
                        None => "MISSING".to_string(),
                    };
                    let rinfo = rmap
                        .get(k)
                        .map(|h| format!(" rust_kmer={} rust_score={}", h.kmer, fmt_score(h.score)))
                        .unwrap_or_default();
                    println!(
                        "  G{:<3} {:<8} len={} cigar={} ref={} kmer={} score={}{}",
                        i,
                        tag,
                        k.0.len(),
                        k.1,
                        k.5,
                        jh.kmer,
                        jh.score.map(|s| s.to_string()).unwrap_or("null".into()),
                        rinfo
                    );
                }
                for (i, (k, rh)) in r.iter().zip(res.haplotypes.iter()).enumerate() {
                    if !gk.contains_key(k) {
                        println!(
                            "  R{:<3} EXTRA    len={} cigar={} ref={} kmer={} score={}",
                            i,
                            k.0.len(),
                            k.1,
                            k.5,
                            rh.kmer,
                            fmt_score(rh.score)
                        );
                    }
                }
            }
        }
    }

    println!("==================== summary ====================");
    println!("regions:                         {}", n_regions);
    println!("haplotype set matches exactly:   {}", n_set_match);
    println!("  ...and order matches:          {}", n_order_match);
    println!("  ...and kmer+score match too:   {}", n_full_match);
    println!("reference haplotype mismatches:  {}", ref_mismatch);
    println!("total haplotypes GATK / Rust:    {} / {}", tot_gatk, tot_rust);
    println!("missing (in GATK, not Rust):     {}", missing);
    println!("extra (in Rust, not GATK):       {}", extra);
    println!("kmer mismatches (matched haps):  {}", kmer_mismatch);
    println!("score mismatches (matched haps): {}", score_mismatch);
    println!("reads with soft clips in input:  {}", reads_with_softclip);
    if gatk_errors > 0 || error_mismatch > 0 {
        println!("GATK exceptions / also in Rust:  {} / {}", gatk_errors, both_errors);
        println!("exception mismatches:            {}", error_mismatch);
    }
    println!("assembly time (excl. JSON):      {:.3} s", asm_time.as_secs_f64());
    let pr = hc_assembly::prof::report() + &hc_assembly::prof::hits_report();
    if !pr.is_empty() {
        print!("{}", pr);
    }
    if !mismatching_regions.is_empty() {
        let show: Vec<String> = mismatching_regions.iter().take(50).map(|x| x.to_string()).collect();
        println!("mismatching region indices (first 50): {}", show.join(","));
    }
}
