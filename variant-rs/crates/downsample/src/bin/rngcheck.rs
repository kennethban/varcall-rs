//! Compares hc_downsample::JavaRandom with java.util.Random: reads the output of
//! probes/RngProbe.java (lines "I bound value" / "G value-bits" / "D value-bits") from stdin and
//! replays the same call sequence.
use hc_downsample::JavaRandom;
use std::io::BufRead;

fn main() {
    let mut rng = JavaRandom::gatk();
    let (mut n, mut bad) = (0u64, 0u64);
    for line in std::io::stdin().lock().lines() {
        let line = line.unwrap();
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.is_empty() { continue; }
        n += 1;
        let ok = match f[0] {
            "S" => { rng = JavaRandom::new(f[1].parse().unwrap()); true }
            "I" => { let b: i32 = f[1].parse().unwrap(); let v: i32 = f[2].parse().unwrap(); rng.next_int(b) == v }
            "G" => { let v = u64::from_str_radix(f[1], 16).unwrap(); rng.next_gaussian().to_bits() == v }
            "D" => { let v = u64::from_str_radix(f[1], 16).unwrap(); rng.next_double().to_bits() == v }
            _ => panic!("bad line {line}"),
        };
        if !ok { bad += 1; if bad < 10 { eprintln!("MISMATCH at call {n}: {line}"); } }
    }
    println!("rngcheck: {n} calls, {bad} mismatches");
    std::process::exit(if bad == 0 { 0 } else { 1 });
}
