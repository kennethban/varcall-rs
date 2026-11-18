//! Validate hc_preproc::finalize_region against a GATK assembly dump (JSONL).
//!
//! usage: validate <dump.jsonl> [max_detail=10]
//!
//! Only the fields needed here are deserialized (serde ignores the rest, e.g. haplotypes).

use hc_preproc::{cigar_to_string, finalize_region_owned, parse_cigar, Read};
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

#[derive(Deserialize)]
struct JRegion {
    contig: String,
    start: i64,
    end: i64,
    padded_start: i64,
    padded_end: i64,
    raw_reads: Vec<JRead>,
    final_reads: Vec<JRead>,
}

fn to_read(j: &JRead) -> Read {
    Read {
        name: j.name.clone(),
        flag: j.flag,
        start: j.start,
        end: j.end,
        cigar: parse_cigar(&j.cigar),
        seq: j.seq.as_bytes().to_vec(),
        qual: j.qual.bytes().map(|b| b - 33).collect(),
        mapq: j.mapq,
        mate_start: j.mate_start,
        tlen: j.tlen,
        ..Default::default()
    }
}

fn qstr(q: &[u8]) -> String {
    q.iter().map(|&b| (b + 33) as char).collect()
}

fn show(r: &Read) -> String {
    format!(
        "{} flag={} {}-{} {} mapq={} mate={} tlen={}\n      seq={}\n      qual={}",
        r.name,
        r.flag,
        r.start,
        r.end,
        cigar_to_string(&r.cigar),
        r.mapq,
        r.mate_start,
        r.tlen,
        String::from_utf8_lossy(&r.seq),
        qstr(&r.qual)
    )
}

#[derive(Default)]
struct Stats {
    regions: u64,
    regions_exact: u64,
    raw_reads: u64,
    expected_reads: u64,
    produced_reads: u64,
    exact: u64,
    missing: u64,
    extra: u64,
    order: u64,
    start: u64,
    end: u64,
    cigar: u64,
    seq: u64,
    qual: u64,
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: validate <dump.jsonl> [max_detail]");
        std::process::exit(2);
    }
    let max_detail: usize = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(10);
    let f = std::fs::File::open(&args[1]).expect("open dump");
    let rdr = BufReader::with_capacity(1 << 22, f);

    let mut st = Stats::default();
    let mut details = 0usize;
    let mut t_parse = Duration::ZERO;
    let mut t_conv = Duration::ZERO;
    let mut t_fin = Duration::ZERO;
    let t_all = Instant::now();

    for line in rdr.lines() {
        let line = line.expect("read line");
        if line.is_empty() {
            continue;
        }
        let t0 = Instant::now();
        let reg: JRegion = serde_json::from_str(&line).expect("parse json");
        let t1 = Instant::now();
        let raw: Vec<Read> = reg.raw_reads.iter().map(to_read).collect();
        let expected: Vec<Read> = reg.final_reads.iter().map(to_read).collect();
        let t2 = Instant::now();
        let got = finalize_region_owned(raw, reg.padded_start, reg.padded_end);
        let t3 = Instant::now();
        t_parse += t1 - t0;
        t_conv += t2 - t1;
        t_fin += t3 - t2;

        st.regions += 1;
        st.raw_reads += reg.raw_reads.len() as u64;
        st.expected_reads += expected.len() as u64;
        st.produced_reads += got.len() as u64;

        // match reads by (name, flag, occurrence)
        let key = |v: &[Read]| -> Vec<(String, u16, u32)> {
            let mut cnt: HashMap<(&str, u16), u32> = HashMap::new();
            v.iter()
                .map(|r| {
                    let c = cnt.entry((r.name.as_str(), r.flag)).or_insert(0);
                    *c += 1;
                    (r.name.clone(), r.flag, *c)
                })
                .collect()
        };
        let ke = key(&expected);
        let kg = key(&got);
        let gi: HashMap<&(String, u16, u32), usize> = kg.iter().enumerate().map(|(i, k)| (k, i)).collect();
        let ei: HashMap<&(String, u16, u32), usize> = ke.iter().enumerate().map(|(i, k)| (k, i)).collect();

        let mut region_ok = true;
        let mut msgs: Vec<String> = Vec::new();
        // order check over the common reads
        let common_e: Vec<&(String, u16, u32)> = ke.iter().filter(|k| gi.contains_key(k)).collect();
        let common_g: Vec<&(String, u16, u32)> = kg.iter().filter(|k| ei.contains_key(k)).collect();
        if common_e != common_g {
            st.order += 1;
            region_ok = false;
            msgs.push("  ORDER differs".to_string());
        }
        for (i, k) in ke.iter().enumerate() {
            let e = &expected[i];
            match gi.get(k) {
                None => {
                    st.missing += 1;
                    region_ok = false;
                    msgs.push(format!("  MISSING read: {}", show(e)));
                }
                Some(&j) => {
                    let g = &got[j];
                    let mut bad = Vec::new();
                    if g.start != e.start {
                        st.start += 1;
                        bad.push("start");
                    }
                    if g.end != e.end {
                        st.end += 1;
                        bad.push("end");
                    }
                    if g.cigar != e.cigar {
                        st.cigar += 1;
                        bad.push("cigar");
                    }
                    if g.seq != e.seq {
                        st.seq += 1;
                        bad.push("seq");
                    }
                    if g.qual != e.qual {
                        st.qual += 1;
                        bad.push("qual");
                    }
                    if bad.is_empty() {
                        st.exact += 1;
                    } else {
                        region_ok = false;
                        msgs.push(format!(
                            "  MISMATCH [{}]\n    expected: {}\n    got:      {}",
                            bad.join(","),
                            show(e),
                            show(g)
                        ));
                    }
                }
            }
        }
        for (j, k) in kg.iter().enumerate() {
            if !ei.contains_key(k) {
                st.extra += 1;
                region_ok = false;
                msgs.push(format!("  EXTRA read: {}", show(&got[j])));
            }
        }
        if region_ok {
            st.regions_exact += 1;
        } else if details < max_detail {
            details += 1;
            println!(
                "Region {}:{}-{} padded {}-{} (raw {}, expected {}, got {})",
                reg.contig,
                reg.start,
                reg.end,
                reg.padded_start,
                reg.padded_end,
                reg.raw_reads.len(),
                expected.len(),
                got.len()
            );
            for m in msgs.iter().take(6) {
                println!("{}", m);
            }
            if msgs.len() > 6 {
                println!("  ... {} more", msgs.len() - 6);
            }
            // print raw read(s) for the first problem read
            if let Some(m) = msgs.iter().find(|m| !m.contains("ORDER")) {
                for rr in reg.raw_reads.iter().filter(|rr| m.contains(&rr.name)).take(2) {
                    println!("    raw:      {}", show(&to_read(rr)));
                }
            }
        }
    }

    let total = t_all.elapsed();
    println!("==== summary ====");
    println!("regions:              {} ({} fully exact)", st.regions, st.regions_exact);
    println!("raw reads:            {}", st.raw_reads);
    println!("expected final reads: {}", st.expected_reads);
    println!("produced final reads: {}", st.produced_reads);
    println!("exact matches:        {}", st.exact);
    println!("missing reads:        {}", st.missing);
    println!("extra reads:          {}", st.extra);
    println!("regions w/ order diff:{}", st.order);
    println!("start mismatches:     {}", st.start);
    println!("end mismatches:       {}", st.end);
    println!("cigar mismatches:     {}", st.cigar);
    println!("seq mismatches:       {}", st.seq);
    println!("qual mismatches:      {}", st.qual);
    let pct = if st.expected_reads > 0 { 100.0 * st.exact as f64 / st.expected_reads as f64 } else { 100.0 };
    println!("read match rate:      {:.6}%", pct);
    if st.regions_exact == st.regions && st.exact == st.expected_reads && st.extra == 0 {
        println!("RESULT: 100% MATCH");
    } else {
        println!("RESULT: MISMATCHES PRESENT");
    }
    println!(
        "time: total {:.3}s | json parse {:.3}s | struct conversion {:.3}s | finalize_region {:.3}s",
        total.as_secs_f64(),
        t_parse.as_secs_f64(),
        t_conv.as_secs_f64(),
        t_fin.as_secs_f64()
    );
}
