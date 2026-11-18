//! Compares the data lines of two (g)VCF files (plain or gzip), joined on POS (one record per POS on one
//! contig, as in HaplotypeCaller GVCFs and the unbanded per-position records).
//!
//! usage: gvcfcmp <gatk.vcf[.gz]> <rust.vcf[.gz]> [--show N] [--only-variants]
//!
//! Prints: lines compared / identical / missing on either side, mismatches by column and by INFO / FORMAT
//! key, and the first N differing pairs.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

fn open(p: &str) -> Box<dyn BufRead> {
    let f = File::open(p).unwrap_or_else(|e| panic!("{p}: {e}"));
    let r: Box<dyn Read> = if p.ends_with(".gz") { Box::new(flate2::read::MultiGzDecoder::new(f)) } else { Box::new(f) };
    Box::new(BufReader::with_capacity(1 << 20, r))
}

struct Lines {
    it: Box<dyn BufRead>,
    buf: String,
}

impl Lines {
    fn next(&mut self) -> Option<(i64, String)> {
        loop {
            self.buf.clear();
            if self.it.read_line(&mut self.buf).unwrap() == 0 {
                return None;
            }
            if self.buf.starts_with('#') {
                continue;
            }
            let l = self.buf.trim_end_matches(['\n', '\r']).to_string();
            let pos: i64 = l.split('\t').nth(1).unwrap().parse().unwrap();
            return Some((pos, l));
        }
    }
}

fn kv(s: &str) -> BTreeMap<String, String> {
    if s == "." {
        return BTreeMap::new();
    }
    s.split(';')
        .map(|x| match x.split_once('=') {
            Some((k, v)) => (k.to_string(), v.to_string()),
            None => (x.to_string(), String::new()),
        })
        .collect()
}

fn fmt_map(keys: &str, vals: &str) -> BTreeMap<String, String> {
    let v: Vec<&str> = vals.split(':').collect();
    keys.split(':').enumerate().map(|(i, k)| (k.to_string(), v.get(i).unwrap_or(&"").to_string())).collect()
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 3 {
        eprintln!("usage: gvcfcmp <gatk.vcf[.gz]> <rust.vcf[.gz]> [--show N] [--only-variants]");
        std::process::exit(2);
    }
    let mut show = 10usize;
    let mut only_var = false;
    let mut i = 3;
    while i < a.len() {
        match a[i].as_str() {
            "--show" => {
                show = a[i + 1].parse().unwrap();
                i += 2;
            }
            "--only-variants" => {
                only_var = true;
                i += 1;
            }
            x => panic!("unknown arg {x}"),
        }
    }
    let is_var = |l: &str| l.split('\t').nth(4).map_or(false, |alt| alt != "<NON_REF>");
    let mut g = Lines { it: open(&a[1]), buf: String::new() };
    let mut r = Lines { it: open(&a[2]), buf: String::new() };
    let (mut ng, mut nr) = (g.next(), r.next());
    let (mut compared, mut identical, mut only_g, mut only_r) = (0u64, 0u64, 0u64, 0u64);
    let (mut var_compared, mut var_identical) = (0u64, 0u64);
    let cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"];
    let mut by_col = [0u64; 10];
    let mut by_info: BTreeMap<String, u64> = BTreeMap::new();
    let mut by_fmt: BTreeMap<String, u64> = BTreeMap::new();
    let mut shown = 0usize;
    let mut examples_missing: Vec<String> = Vec::new();
    loop {
        match (&ng, &nr) {
            (None, None) => break,
            (Some((pg, lg)), Some((pr, lr))) if pg == pr => {
                if !only_var || is_var(lg) || is_var(lr) {
                    compared += 1;
                    let v = is_var(lg) || is_var(lr);
                    if v {
                        var_compared += 1;
                    }
                    if lg == lr {
                        identical += 1;
                        if v {
                            var_identical += 1;
                        }
                    } else {
                        let fg: Vec<&str> = lg.split('\t').collect();
                        let fr: Vec<&str> = lr.split('\t').collect();
                        for c in 0..10 {
                            if fg.get(c) != fr.get(c) {
                                by_col[c] += 1;
                            }
                        }
                        let (ig, ir) = (kv(fg.get(7).unwrap_or(&".")), kv(fr.get(7).unwrap_or(&".")));
                        for k in ig.keys().chain(ir.keys()) {
                            if ig.get(k) != ir.get(k) {
                                *by_info.entry(k.clone()).or_default() += 1;
                            }
                        }
                        let (mg, mr) = (
                            fmt_map(fg.get(8).unwrap_or(&""), fg.get(9).unwrap_or(&"")),
                            fmt_map(fr.get(8).unwrap_or(&""), fr.get(9).unwrap_or(&"")),
                        );
                        for k in mg.keys().chain(mr.keys()) {
                            if mg.get(k) != mr.get(k) {
                                *by_fmt.entry(k.clone()).or_default() += 1;
                            }
                        }
                        if shown < show {
                            println!("- GATK: {lg}\n+ rust: {lr}");
                            shown += 1;
                        }
                    }
                }
                ng = g.next();
                nr = r.next();
            }
            (Some((pg, lg)), Some((pr, _))) if pg < pr => {
                if !only_var || is_var(lg) {
                    only_g += 1;
                    if examples_missing.len() < show {
                        examples_missing.push(format!("only GATK: {lg}"));
                    }
                }
                ng = g.next();
            }
            (Some(_), Some((_, lr))) => {
                if !only_var || is_var(lr) {
                    only_r += 1;
                    if examples_missing.len() < show {
                        examples_missing.push(format!("only rust: {lr}"));
                    }
                }
                nr = r.next();
            }
            (Some((_, lg)), None) => {
                if !only_var || is_var(lg) {
                    only_g += 1;
                    if examples_missing.len() < show {
                        examples_missing.push(format!("only GATK: {lg}"));
                    }
                }
                ng = g.next();
            }
            (None, Some((_, lr))) => {
                if !only_var || is_var(lr) {
                    only_r += 1;
                    if examples_missing.len() < show {
                        examples_missing.push(format!("only rust: {lr}"));
                    }
                }
                nr = r.next();
            }
        }
    }
    for e in &examples_missing {
        println!("{e}");
    }
    println!("lines compared (same POS): {compared}   identical: {identical}   differing: {}", compared - identical);
    println!("  of which with a real ALT on either side: {var_compared}   identical: {var_identical}");
    println!("only in GATK: {only_g}   only in rust: {only_r}");
    let mut s = String::new();
    for c in 0..10 {
        if by_col[c] > 0 {
            s.push_str(&format!(" {}={}", cols[c], by_col[c]));
        }
    }
    println!("differing lines by column:{}", if s.is_empty() { " none".into() } else { s });
    println!("differing INFO keys: {:?}", by_info);
    println!("differing FORMAT keys: {:?}", by_fmt);
    if compared == identical && only_g == 0 && only_r == 0 {
        println!("RESULT: IDENTICAL");
    } else {
        println!("RESULT: DIFFERENT");
        std::process::exit(1);
    }
}
