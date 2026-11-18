//! Validator: runs the genotyping stage on GATK's dumped inputs and compares against GATK.
use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::time::Instant;

use hc_assembly::cigar;
use hc_genotype::{Haplotype, Read};
use serde::Deserialize;

#[derive(Deserialize)]
struct JRead {
    name: String,
    flag: u16,
    start: i64,
    end: i64,
    cigar: String,
    seq: String,
    qual: String,
    mapq: i64,
    mate_start: i64,
    tlen: i64,
}

#[derive(Deserialize)]
struct JHap {
    bases: String,
    cigar: String,
    align_start_wrt_ref: i64,
    loc_start: i64,
    loc_end: i64,
    is_ref: bool,
    score: Option<serde_json::Value>,
    kmer: u32,
}

#[derive(Deserialize)]
struct JRegion {
    start: i64,
    end: i64,
    exit: Option<String>,
    geno_start: Option<i64>,
    geno_end: Option<i64>,
    filtered_reads: Option<Vec<JRead>>,
    lk_haplotypes: Option<Vec<JHap>>,
    lk_reads: Option<Vec<JRead>>,
    lk_filtered_reads: Option<Vec<JRead>>,
    lk_matrix: Option<Vec<Vec<serde_json::Value>>>,
    realigned_reads: Option<Vec<JRead>>,
    n_calls: i64,
    ref_padded: Option<String>,
    ref_padded_start: Option<i64>,
    ref_padded_end: Option<i64>,
}

#[derive(Deserialize)]
struct JAsm {
    start: i64,
    ref_padded_start: i64,
    ref_padded_end: i64,
    ref_padded: String,
}

fn num(v: &serde_json::Value) -> f64 {
    match v {
        serde_json::Value::Number(n) => n.as_f64().unwrap(),
        serde_json::Value::String(s) => match s.as_str() {
            "NaN" => f64::NAN,
            "Infinity" => f64::INFINITY,
            "-Infinity" => f64::NEG_INFINITY,
            _ => panic!("bad number {}", s),
        },
        serde_json::Value::Null => f64::NAN,
        _ => panic!("bad number"),
    }
}

fn conv_read(r: &JRead) -> Read {
    Read {
        name: r.name.clone(),
        flag: r.flag,
        start: r.start,
        end: r.end,
        cigar: cigar::parse(&r.cigar),
        seq: r.seq.as_bytes().to_vec(),
        qual: r.qual.bytes().map(|q| q - 33).collect(),
        mapq: r.mapq as u8,
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
}

fn conv_hap(h: &JHap) -> Haplotype {
    Haplotype {
        bases: h.bases.as_bytes().to_vec(),
        cigar: cigar::parse(&h.cigar),
        align_start_wrt_ref: h.align_start_wrt_ref,
        loc_start: h.loc_start,
        loc_end: h.loc_end,
        is_ref: h.is_ref,
        score: h.score.as_ref().map(num).unwrap_or(f64::NAN),
        kmer: h.kmer,
    }
}


fn load_vcf(path: &str) -> Vec<String> {
    let out = std::process::Command::new("bash")
        .arg("-c")
        .arg(format!("{} '{}' | grep -v '^#'", if path.ends_with(".gz") { "zcat" } else { "cat" }, path))
        .output()
        .expect("zcat");
    String::from_utf8(out.stdout).unwrap().lines().map(|l| l.to_string()).collect()
}

const COLS: [&str; 10] = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"];

fn kv_map(s: &str, sep: char) -> Vec<(String, String)> {
    s.split(sep)
        .map(|f| match f.split_once('=') {
            Some((k, v)) => (k.to_string(), v.to_string()),
            None => (f.to_string(), String::new()),
        })
        .collect()
}

#[derive(Default)]
struct Stats {
    lines_cmp: u64,
    identical: u64,
    col_mismatch: [u64; 10],
    info_key: std::collections::BTreeMap<String, u64>,
    fmt_key: std::collections::BTreeMap<String, u64>,
    shown: usize,
}

fn compare_line(st: &mut Stats, got: &str, exp: &str, show: usize) {
    st.lines_cmp += 1;
    if got == exp {
        st.identical += 1;
        return;
    }
    let g: Vec<&str> = got.split('\t').collect();
    let e: Vec<&str> = exp.split('\t').collect();
    for c in 0..10 {
        if g.get(c) != e.get(c) {
            st.col_mismatch[c] += 1;
        }
    }
    if g.get(7) != e.get(7) {
        let gm: HashMap<String, String> = kv_map(g[7], ';').into_iter().collect();
        let em: HashMap<String, String> = kv_map(e[7], ';').into_iter().collect();
        let mut keys: Vec<&String> = gm.keys().chain(em.keys()).collect();
        keys.sort();
        keys.dedup();
        for k in keys {
            if gm.get(k) != em.get(k) {
                *st.info_key.entry(k.clone()).or_default() += 1;
            }
        }
    }
    if g.get(8) != e.get(8) || g.get(9) != e.get(9) {
        let gk: Vec<&str> = g[8].split(':').collect();
        let gv: Vec<&str> = g[9].split(':').collect();
        let ek: Vec<&str> = e[8].split(':').collect();
        let ev: Vec<&str> = e[9].split(':').collect();
        let gm: HashMap<&str, &str> = gk.iter().copied().zip(gv.iter().copied()).collect();
        let em: HashMap<&str, &str> = ek.iter().copied().zip(ev.iter().copied()).collect();
        let mut keys: Vec<&str> = gk.iter().chain(ek.iter()).copied().collect();
        keys.sort();
        keys.dedup();
        for k in keys {
            if gm.get(k) != em.get(k) {
                *st.fmt_key.entry(k.to_string()).or_default() += 1;
            }
        }
    }
    if st.shown < show {
        st.shown += 1;
        println!("MISMATCH\n  rust: {}\n  gatk: {}", got, exp);
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    // defaults: the dumps under $VARCALL_DATA (tests/harness/README.md)
    let data = std::env::var("VARCALL_DATA").unwrap_or_else(|_| ".".to_string());
    let mut geno_path = format!("{data}/geno/gatk_genotyping.jsonl");
    let mut asm_path = format!("{data}/asm/gatk_assembly.jsonl");
    let mut vcf_path = format!("{data}/geno/gatk.vcf.gz");
    let mut show = 10usize;
    let mut out_vcf: Option<String> = None;
    let mut erc = false;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--geno" => { geno_path = args[i + 1].clone(); i += 1; }
            "--asm" => { asm_path = args[i + 1].clone(); i += 1; }
            "--vcf" => { vcf_path = args[i + 1].clone(); i += 1; }
            "--show" => { show = args[i + 1].parse().unwrap(); i += 1; }
            "--out" => { out_vcf = Some(args[i + 1].clone()); i += 1; }
            "--erc" => { erc = true; }
            _ => panic!("unknown argument {}", args[i]),
        }
        i += 1;
    }
    let t0 = Instant::now();
    let mut asm: HashMap<i64, JAsm> = HashMap::new();
    if asm_path != "-" {
        for line in BufReader::with_capacity(1 << 22, std::fs::File::open(&asm_path).unwrap()).lines() {
            let a: JAsm = serde_json::from_str(&line.unwrap()).unwrap();
            asm.insert(a.start, a);
        }
    }
    let gatk_lines = load_vcf(&vcf_path);
    eprintln!("loaded {} assembly regions, {} GATK VCF records in {:.1}s", asm.len(), gatk_lines.len(), t0.elapsed().as_secs_f64());

    let mut engine = hc_genotype::GenotypingEngine::new(if erc {
        hc_genotype::EngineConfig::new_erc("chr20", 64444167)
    } else {
        hc_genotype::EngineConfig::new("chr20", 64444167)
    });
    let mut sw = hc_assembly::sw::SwAligner::new();
    let (mut n_reads, mut n_start_ok, mut n_cigar_ok, mut n_both_ok, mut n_changed) = (0u64, 0u64, 0u64, 0u64, 0u64);
    let mut shown_realign = 0;
    let (mut realign_time, mut geno_time) = (0f64, 0f64);
    let mut n_regions = 0;
    let (mut regions_count_ok, mut regions_count_bad) = (0u64, 0u64);
    let mut st = Stats::default();
    let mut gi = 0usize; // cursor into GATK lines (sorted by position)
    let mut n_rust_lines = 0u64;
    let mut all_rust_lines: Vec<String> = Vec::new();
    let mut count_bad_shown = 0;
    for line in BufReader::with_capacity(1 << 22, std::fs::File::open(&geno_path).unwrap()).lines() {
        let line = line.unwrap();
        let reg: JRegion = serde_json::from_str(&line).unwrap();
        // GATK records of this region: POS within [start, end]
        let gpos = |l: &str| -> i64 { l.split('\t').nth(1).unwrap().parse().unwrap() };
        while gi < gatk_lines.len() && gpos(&gatk_lines[gi]) < reg.start {
            gi += 1;
        }
        let mut gend = gi;
        while gend < gatk_lines.len() && gpos(&gatk_lines[gend]) <= reg.end {
            gend += 1;
        }
        let gatk_region: Vec<&String> = gatk_lines[gi..gend].iter().collect();
        gi = gend;
        if reg.exit.as_deref() != Some("genotyped") {
            // (in -ERC GVCF dumps n_calls counts the reference-confidence records too)
            if !gatk_region.is_empty() || (reg.n_calls != 0 && !erc) {
                println!("region {} not genotyped but has {} GATK records / n_calls {}", reg.start, gatk_region.len(), reg.n_calls);
            }
            continue;
        }
        n_regions += 1;
        let inline;
        let a = match &reg.ref_padded {
            Some(rp) => {
                inline = JAsm { start: reg.start, ref_padded_start: reg.ref_padded_start.unwrap(), ref_padded_end: reg.ref_padded_end.unwrap(), ref_padded: rp.clone() };
                &inline
            }
            None => &asm[&reg.start],
        };
        let haps: Vec<Haplotype> = reg.lk_haplotypes.as_ref().unwrap().iter().map(conv_hap).collect();
        let orig: Vec<Read> = reg.lk_reads.as_ref().unwrap().iter().map(conv_read).collect();
        let filtered: Vec<Read> = reg.filtered_reads.as_ref().unwrap().iter().map(conv_read).collect();
        let lk_filtered: Vec<Read> = reg.lk_filtered_reads.as_ref().map(|v| v.iter().map(conv_read).collect()).unwrap_or_default();
        let matrix: Vec<Vec<f64>> =
            reg.lk_matrix.as_ref().unwrap().iter().map(|row| row.iter().map(num).collect()).collect();
        let mut reads = orig.clone();
        let t = Instant::now();
        hc_genotype::realign::realign_reads_to_best_haplotype(&mut sw, &mut reads, &haps, &matrix, a.ref_padded_start);
        realign_time += t.elapsed().as_secs_f64();
        let exp = reg.realigned_reads.as_ref().unwrap();
        assert_eq!(exp.len(), reads.len());
        for ((r, e), o) in reads.iter().zip(exp.iter()).zip(orig.iter()) {
            n_reads += 1;
            let s_ok = r.start == e.start;
            let c_ok = cigar::to_string(&r.cigar) == e.cigar;
            n_start_ok += s_ok as u64;
            n_cigar_ok += c_ok as u64;
            n_both_ok += (s_ok && c_ok) as u64;
            n_changed += (o.start != e.start || cigar::to_string(&o.cigar) != e.cigar) as u64;
            if !(s_ok && c_ok) && shown_realign < 10 {
                shown_realign += 1;
                println!("realign mismatch region {} read {}: got {} {} expected {} {}", reg.start, r.name, r.start, cigar::to_string(&r.cigar), e.start, e.cigar);
            }
        }
        let input = hc_genotype::RegionInput {
            haplotypes: &haps,
            reads: &orig,
            matrix: &matrix,
            filtered_reads: &filtered,
            filtered_evidence: &lk_filtered,
            ref_padded: a.ref_padded.as_bytes(),
            ref_padded_start: a.ref_padded_start,
            span_start: reg.geno_start.unwrap(),
            span_end: reg.geno_end.unwrap(),
        };
        let t = Instant::now();
        let recs = if erc {
            engine.genotype_region_erc(&input, &reads).into_iter().map(|c| c.rec).collect()
        } else {
            engine.genotype_region(&input, &reads)
        };
        geno_time += t.elapsed().as_secs_f64();
        let lines: Vec<String> = recs.iter().map(|r| r.to_line()).collect();
        n_rust_lines += lines.len() as u64;
        if (erc || lines.len() as i64 == reg.n_calls) && gatk_region.len() == lines.len() {
            regions_count_ok += 1;
        } else {
            regions_count_bad += 1;
            if count_bad_shown < 10 {
                count_bad_shown += 1;
                println!("COUNT region {}-{}: rust {} n_calls {} gatk vcf {}", reg.start, reg.end, lines.len(), reg.n_calls, gatk_region.len());
                for l in &lines { println!("  rust: {}", l); }
                for l in &gatk_region { println!("  gatk: {}", l); }
            }
        }
        // pair up by position
        let mut used = vec![false; gatk_region.len()];
        for l in &lines {
            let p: i64 = l.split('\t').nth(1).unwrap().parse().unwrap();
            if let Some(k) = (0..gatk_region.len()).find(|&k| !used[k] && gpos(gatk_region[k]) == p) {
                used[k] = true;
                compare_line(&mut st, l, gatk_region[k], show);
            }
        }
        let _ = (reg.end, a.ref_padded_end);
        if out_vcf.is_some() {
            all_rust_lines.extend(lines);
        }
    }
    if let Some(p) = out_vcf {
        std::fs::write(&p, all_rust_lines.join("\n") + "\n").unwrap();
    }
    println!("== realignment ({} genotyped regions)", n_regions);
    println!("reads: {}  (changed by realignment: {})  start match: {}  cigar match: {}  both: {}", n_reads, n_changed, n_start_ok, n_cigar_ok, n_both_ok);
    println!("== calls per region");
    println!("regions with rust count == {}GATK VCF count: {}   differing: {}", if erc { "" } else { "n_calls == " }, regions_count_ok, regions_count_bad);
    println!("records: rust {}  GATK VCF {}", n_rust_lines, gatk_lines.len());
    println!("== VCF lines (paired by region and POS)");
    println!("compared: {}  identical: {}  unpaired GATK: {}  unpaired rust: {}", st.lines_cmp, st.identical, gatk_lines.len() as u64 - st.lines_cmp, n_rust_lines - st.lines_cmp);
    for c in 0..10 {
        if st.col_mismatch[c] > 0 {
            println!("  column {:7} mismatches: {}", COLS[c], st.col_mismatch[c]);
        }
    }
    for (k, v) in &st.info_key {
        println!("  INFO {:15} mismatches: {}", k, v);
    }
    for (k, v) in &st.fmt_key {
        println!("  FORMAT {:13} mismatches: {}", k, v);
    }
    println!("QD jitter draws: {}", engine.n_qd_jitter);
    println!("code paths: {}", hc_genotype::stats::report());
    println!("== time: realignment {:.2}s, genotyping+annotation+formatting {:.2}s, total wall {:.1}s", realign_time, geno_time, t0.elapsed().as_secs_f64());
}
