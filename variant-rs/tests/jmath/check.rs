// Reads JMathGen.java records from stdin and compares them with hc-rs/assembly/src/jmath.rs, bit for bit.
// Built and run by run.sh (plain rustc; jmath.rs depends only on fdlibm.rs).
#[path = "../../crates/assembly/src/jmath.rs"]
#[allow(dead_code)]
mod jmath;
#[path = "../../crates/assembly/src/fdlibm.rs"]
mod fdlibm;

use std::io::{BufReader, Read};

const NAMES: [&str; 6] = ["exp", "log", "log10", "pow", "log1p", "expm1"];
const NF: usize = NAMES.len();

fn main() {
    let mut r = BufReader::with_capacity(1 << 20, std::io::stdin().lock());
    let mut rec = [0u8; 26];
    let mut total = [[0u64; 16]; NF];
    let mut bad = [[0u64; 16]; NF];
    let mut nan_payload = [0u64; NF];
    let mut libm_diff = [[0u64; 16]; NF];
    let mut shown = 0;
    while r.read_exact(&mut rec).is_ok() {
        let (fun, cat) = (rec[0] as usize, rec[1] as usize);
        let u = |o: usize| u64::from_le_bytes(rec[o..o + 8].try_into().unwrap());
        let (x, y, want) = (f64::from_bits(u(2)), f64::from_bits(u(10)), u(18));
        let (got, libm) = match fun {
            0 => (jmath::exp(x), x.exp()),
            1 => (jmath::log(x), x.ln()),
            2 => (jmath::log10(x), x.log10()),
            3 => (jmath::pow(x, y), x.powf(y)),
            4 => (jmath::log1p(x), x.ln_1p()),
            _ => (jmath::expm1(x), x.exp_m1()),
        };
        total[fun][cat] += 1;
        let w = f64::from_bits(want);
        if libm.to_bits() != want && !(libm.is_nan() && w.is_nan()) {
            libm_diff[fun][cat] += 1;
        }
        if got.to_bits() == want {
            continue;
        }
        if got.is_nan() && w.is_nan() {
            nan_payload[fun] += 1;
            continue;
        }
        bad[fun][cat] += 1;
        if shown < 30 {
            shown += 1;
            eprintln!(
                "MISMATCH {}({:e} [{:#018x}], {:e} [{:#018x}]) cat {}: java {:e} [{:#018x}] rust {:e} [{:#018x}]",
                NAMES[fun], x, x.to_bits(), y, y.to_bits(), cat, w, want, got, got.to_bits()
            );
        }
    }
    let mut fail = false;
    for f in 0..NF {
        let t: u64 = total[f].iter().sum();
        let b: u64 = bad[f].iter().sum();
        fail |= b > 0 || t == 0; // no input at all (generator failed) is a failure too
        let per: Vec<String> = (0..16).filter(|&c| total[f][c] > 0).map(|c| format!("{}:{}/{}", c, bad[f][c], total[f][c])).collect();
        let lper: Vec<String> = (0..16).filter(|&c| total[f][c] > 0).map(|c| format!("{}:{}", c, libm_diff[f][c])).collect();
        let l: u64 = libm_diff[f].iter().sum();
        println!(
            "{:5}  inputs {:>9}  mismatches {:>6}  (NaN payload only {})  by category {}\n       [libm would differ on {}; by category {}]",
            NAMES[f], t, b, nan_payload[f], per.join(" "), l, lper.join(" ")
        );
    }
    if fail {
        std::process::exit(1);
    }
}
