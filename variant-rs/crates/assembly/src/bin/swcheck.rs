//! Compare the Rust Smith-Waterman / calculateCigar ports against GATK results recorded by
//! scratch/SwProbe.java, a development probe not in the repository.   usage: swcheck <probe.tsv>
use hc_assembly::cigar;
use hc_assembly::cigar_utils;
use hc_assembly::sw::{self, Overhang, FixedRefAligner, SwAligner};
use std::io::{BufRead, BufReader};

fn main() {
    let path = std::env::args().nth(1).expect("usage: swcheck <probe.tsv>");
    let f = BufReader::new(std::fs::File::open(path).unwrap());
    let mut aligner = SwAligner::new();
    let mut cache = FixedRefAligner::new();
    let mut pad = (Vec::new(), Vec::new());
    let (mut n_sw, mut bad_sw, mut n_cc, mut bad_cc, mut n_cache, mut bad_cache) = (0, 0, 0, 0, 0, 0);
    for line in f.lines() {
        let line = line.unwrap();
        let t: Vec<&str> = line.split('\t').collect();
        let params = if t[1] == "0" { sw::STANDARD_NGS } else { sw::NEW_SW_PARAMETERS };
        let strat = match t[2] {
            "SOFTCLIP" => Overhang::SoftClip,
            "INDEL" => Overhang::Indel,
            "LEADING_INDEL" => Overhang::LeadingIndel,
            _ => Overhang::Ignore,
        };
        let (r, a) = (t[3].as_bytes(), t[4].as_bytes());
        if t[0] == "SW" {
            n_sw += 1;
            let res = aligner.align(r, a, &params, strat);
            let got = format!("{}\t{}", cigar::to_string(&res.cigar), res.offset);
            let exp = format!("{}\t{}", t[5], t[6]);
            if got != exp {
                bad_sw += 1;
                if bad_sw <= 5 {
                    println!("SW mismatch {:?} {}: expected {} got {}\n  ref={}\n  alt={}", strat, t[1], exp, got, t[3], t[4]);
                }
            }
            if strat == Overhang::SoftClip {
                n_cache += 1;
                cache.set_reference(r, &params);
                let res2 = cache.align(a);
                let got2 = format!("{}\t{}", cigar::to_string(&res2.cigar), res2.offset);
                if got2 != exp {
                    bad_cache += 1;
                    if bad_cache <= 5 {
                        println!("cached SW mismatch: expected {} got {}", exp, got2);
                    }
                }
            }
        } else {
            n_cc += 1;
            let c = cigar_utils::calculate_cigar(&mut aligner, r, a, &params, strat, &mut pad);
            let got = c.map(|c| cigar::to_string(&c)).unwrap_or_else(|| "null".into());
            if got != t[5] {
                bad_cc += 1;
                if bad_cc <= 5 {
                    println!("calculateCigar mismatch {:?}: expected {} got {}\n  ref={}\n  alt={}", strat, t[5], got, t[3], t[4]);
                }
            }
        }
    }
    println!("SW cases: {} mismatches: {}", n_sw, bad_sw);
    println!("FixedRefAligner SOFTCLIP SW cases: {} mismatches: {}", n_cache, bad_cache);
    println!("calculateCigar cases: {} mismatches: {}", n_cc, bad_cc);
}
