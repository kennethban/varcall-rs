//! Compares the Rust numerics against Java test vectors (scratch/probes/NumProbe.java, a development probe not in the repository).
use std::io::{BufRead, BufReader};

use hc_genotype::annotate;
use hc_genotype::javamath::{fast_log, format_fixed};

fn bits(s: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(s, 16).unwrap())
}

fn main() {
    let path = std::env::args().nth(1).expect("tsv");
    let mut fmt = (0u64, 0u64);
    let mut fns = [(0u64, 0u64); 4];
    let mut sb = (0u64, 0u64, 0u64);
    let mut mwu = (0u64, 0u64, 0u64);
    let mut shown = 0;
    for line in BufReader::new(std::fs::File::open(&path).unwrap()).lines() {
        let line = line.unwrap();
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "fmt" => {
                let d = bits(f[1]);
                fmt.0 += 1;
                let ok = format_fixed(d, 2) == f[2] && format_fixed(d, 3) == f[3] && format_fixed(d, 4) == f[4];
                if !ok {
                    fmt.1 += 1;
                    if shown < 10 {
                        shown += 1;
                        println!("fmt {:e}: java {} {} {} rust {} {} {}", d, f[2], f[3], f[4], format_fixed(d, 2), format_fixed(d, 3), format_fixed(d, 4));
                    }
                }
            }
            "fn" => {
                let x = bits(f[1]);
                let y = bits(f[5]);
                let r = [fast_log(x), x.ln(), x.log10(), 10f64.powf(y)];
                let j = [bits(f[2]), bits(f[3]), bits(f[4]), bits(f[6])];
                for k in 0..4 {
                    fns[k].0 += 1;
                    if r[k].to_bits() != j[k].to_bits() {
                        fns[k].1 += 1;
                    }
                }
            }
            "sb" => {
                let t = [[f[1].parse().unwrap(), f[2].parse().unwrap()], [f[3].parse().unwrap(), f[4].parse().unwrap()]];
                sb.0 += 1;
                if annotate::fisher_strand(t) != f[5] {
                    sb.1 += 1;
                    if shown < 20 {
                        shown += 1;
                        println!("FS {:?}: java {} rust {}", t, f[5], annotate::fisher_strand(t));
                    }
                }
                if annotate::strand_odds_ratio(t) != f[6] {
                    sb.2 += 1;
                }
            }
            "mwu" => {
                let p = |s: &str| -> Vec<f64> { s.split(',').map(|x| x.parse::<f64>().unwrap()).collect() };
                let (s1, s2) = (p(f[1]), p(f[2]));
                let z = annotate::mann_whitney_z(&s1, &s2);
                mwu.0 += 1;
                let jz = bits(f[3]);
                if z.to_bits() != jz.to_bits() && !(z.is_nan() && jz.is_nan()) {
                    mwu.1 += 1;
                    if shown < 30 {
                        shown += 1;
                        println!("MWU {:?} {:?}: java {} rust {}", s1, s2, jz, z);
                    }
                }
                if format_fixed(z, 3) != f[4] {
                    mwu.2 += 1;
                }
            }
            "const" => {
                println!("LOG10_E java==rust: {}", bits(f[2]) == std::f64::consts::LOG10_E);
                println!("log10(MIN_VALUE) java==rust: {}", bits(f[4]) == f64::from_bits(1).log10());
                println!("SQRT2 java==rust: {}", bits(f[6]) == 2f64.sqrt());
            }
            _ => {}
        }
    }
    println!("format_fixed: {} values, {} mismatching (any of %.2f/%.3f/%.4f)", fmt.0, fmt.1);
    let names = ["FastMath.log vs ln", "Math.log vs ln", "Math.log10 vs log10", "Math.pow(10,y) vs powf"];
    for k in 0..4 {
        println!("{}: {} values, {} differ in the last bit", names[k], fns[k].0, fns[k].1);
    }
    println!("FisherStrand: {} tables, {} FS strings differ, {} SOR strings differ", sb.0, sb.1, sb.2);
    println!("MannWhitneyU: {} tests, {} z values differ bitwise, {} \"%.3f\" strings differ", mwu.0, mwu.1, mwu.2);
}
