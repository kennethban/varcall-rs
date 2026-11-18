//! Checks `calculate_ref_confidence` against GATK's `ReferenceConfidenceModel.calculateRefConfidence` on the
//! synthetic regions written by `scratch/fuzz/RefconfFuzz.java` (a development probe, not in the repository).
//!
//! usage: refconf_fuzzcheck <fuzz.txt> [--show N]

use hc_assembly::Read;
use hc_refconf::{calculate_ref_confidence, CallSpan, Entry};

fn parse_cigar(s: &str) -> Vec<(u8, u32)> {
    let mut v = Vec::new();
    let mut n = 0u32;
    for b in s.bytes() {
        if b.is_ascii_digit() {
            n = n * 10 + (b - b'0') as u32;
        } else {
            v.push((b, n));
            n = 0;
        }
    }
    v
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let show: usize = a.iter().position(|x| x == "--show").map(|i| a[i + 1].parse().unwrap()).unwrap_or(5);
    let text = std::fs::read_to_string(&a[1]).unwrap();
    let (mut regions, mut ok_regions, mut positions, mut mism) = (0, 0, 0u64, 0u64);
    let mut shown = 0;
    let mut lines = text.lines();
    let c0 = hc_refconf::stats::all().iter().map(|c| c.1.load(std::sync::atomic::Ordering::Relaxed)).collect::<Vec<_>>();
    let _ = c0;
    while let Some(l) = lines.next() {
        if !l.starts_with("REGION") {
            continue;
        }
        let f: Vec<i64> = l.split(' ').skip(1).map(|x| x.parse().unwrap()).collect();
        let (s0, e0, ps, _pe) = (f[0], f[1], f[2], f[3]);
        let refb = lines.next().unwrap().strip_prefix("REF ").unwrap().as_bytes().to_vec();
        let mut reads: Vec<Read> = Vec::new();
        let mut calls: Vec<CallSpan> = Vec::new();
        let mut expected: Vec<String> = Vec::new();
        for l in lines.by_ref() {
            if l == "END" {
                break;
            }
            let t: Vec<&str> = l.split(' ').collect();
            match t[0] {
                "READ" => {
                    let cigar = parse_cigar(t[4]);
                    let start: i64 = t[3].parse().unwrap();
                    let rl: i64 = cigar.iter().filter(|e| matches!(e.0, b'M' | b'D' | b'N' | b'=' | b'X')).map(|e| e.1 as i64).sum();
                    reads.push(Read {
                        name: t[1].to_string(),
                        flag: t[2].parse().unwrap(),
                        start,
                        end: start + rl - 1,
                        cigar,
                        seq: t[5].as_bytes().to_vec(),
                        qual: t[6].bytes().map(|q| q - 33).collect(),
                        mapq: t[7].parse().unwrap(),
                        mate_start: t[8].parse().unwrap(),
                        tlen: t[9].parse().unwrap(),
                        ref_id: 0,
                        mate_ref_id: 0,
                        bi: Vec::new(),
                        bd: Vec::new(),
                        sample: 0,
                        oa_other_contig: false,
                        bam_rec: None,
                    });
                }
                "CALL" => calls.push(CallSpan { start: t[1].parse().unwrap(), end: t[2].parse().unwrap() }),
                "OUT" => expected.push(l.to_string()),
                _ => panic!("bad line {l}"),
            }
        }
        regions += 1;
        let rr: Vec<&Read> = reads.iter().collect();
        let got: Vec<String> = calculate_ref_confidence(&refb, ps, s0, e0, &rr, &calls)
            .into_iter()
            .enumerate()
            .map(|(i, e)| {
                let pos = s0 + i as i64;
                match e {
                    Entry::Call(k) => format!("OUT {} C {}", pos, k),
                    Entry::Site(s) => format!(
                        "OUT {} {} {} {} {} {} {},{},{}",
                        s.pos, s.ref_base as char, s.ad[0], s.ad[1], s.dp(), s.gq, s.pl[0], s.pl[1], s.pl[2]
                    ),
                }
            })
            .collect();
        let mut bad = got.len() != expected.len();
        for (g, x) in got.iter().zip(expected.iter()) {
            positions += 1;
            if g != x {
                mism += 1;
                bad = true;
                if shown < show {
                    println!("region {s0}-{e0}: GATK `{x}` rust `{g}`");
                    shown += 1;
                }
            }
        }
        if !bad {
            ok_regions += 1;
        }
    }
    println!("regions: {regions}  identical: {ok_regions}  positions: {positions}  mismatching positions: {mism}");
    println!("paths: {}", hc_refconf::stats::report());
    if ok_regions != regions {
        std::process::exit(1);
    }
}
