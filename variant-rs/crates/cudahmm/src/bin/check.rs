//! Validates the Rust->CUDA path against GATK's AVX PairHMM results recorded in a
//! GATK_PAIRHMM_DUMP file (see variant-rs/cuda-pairhmm/native/cuda_pairhmm.cu). usage: check dump.bin [regions_per_launch]
use hc_cudahmm::*;
use std::time::Instant;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let data = std::fs::read(&a[1]).unwrap();
    let merge: usize = a.get(2).map(|s| s.parse().unwrap()).unwrap_or(64);
    let mut o = 0usize;
    let ri = |o: &mut usize| { let v = i32::from_be_bytes(data[*o..*o + 4].try_into().unwrap()); *o += 4; v as usize };
    struct B { haps: Vec<Vec<u8>>, reads: Vec<[Vec<u8>; 5]>, exp: Vec<f64> }
    let mut batches = vec![];
    while o < data.len() {
        let (nr, nh) = (ri(&mut o), ri(&mut o));
        let rbytes = |o: &mut usize| { let l = ri(o); let v = data[*o..*o + l].to_vec(); *o += l; v };
        let haps: Vec<Vec<u8>> = (0..nh).map(|_| rbytes(&mut o)).collect();
        let reads: Vec<[Vec<u8>; 5]> = (0..nr).map(|_| [rbytes(&mut o), rbytes(&mut o), rbytes(&mut o), rbytes(&mut o), rbytes(&mut o)]).collect();
        let exp = (0..nr * nh).map(|_| { let v = f64::from_be_bytes(data[o..o + 8].try_into().unwrap()); o += 8; v }).collect();
        batches.push(B { haps, reads, exp });
    }
    assert!(init(false), "no CUDA device");
    let t = Instant::now();
    let (mut maxd, mut n, mut big) = (0f64, 0usize, 0usize); let (mut same, mut lowdiff, mut low) = (0usize, 0usize, 0usize);
    for chunk in batches.chunks(merge) {
        let groups: Vec<HmmGroup> = chunk.iter().map(|b| HmmGroup {
            reads: b.reads.iter().map(|r| HmmRead { bases: &r[0], quals: &r[1], ins: &r[2], del: &r[3], gcp: &r[4] }).collect(),
            haps: b.haps.iter().map(|h| h.as_slice()).collect(),
        }).collect();
        for (b, res) in chunk.iter().zip(compute(&groups)) {
            for (k, (x, y)) in b.exp.iter().zip(res).enumerate() { let d = (x - y).abs(); maxd = maxd.max(d); if d > 1e-3 { big += 1; } n += 1; if x.to_bits() == y.to_bits() { same += 1; } else if std::env::var("SHOW").is_ok() && lowdiff < 12 { lowdiff += 1; let r = &b.reads[k / b.haps.len()]; let h = &b.haps[k % b.haps.len()]; eprintln!("diff: gkl={:.9} cuda={:.9} d={:.2e} rlen={} hlen={} readN={} hapN={} minq={} min_gcp={}", x, y, d, r[0].len(), h.len(), r[0].contains(&b'N'), h.contains(&b'N'), r[1].iter().min().unwrap(), r[4].iter().min().unwrap()); } if *x < -28.0 + 36.0 - 36.0 - 0.0 && *x < -28.0 { low += 1; if d > 0.0 { lowdiff += 1; } } }
        }
    }
    println!("bit-identical={} of {}; values < -28 (double-fallback range): {} of which differ: {}", same, n, low, lowdiff);
    println!("batches={} pairs={} max|diff vs GKL AVX|={:.3e} n(>1e-3)={} time={:.3}s (regions/launch={})", batches.len(), n, maxd, big, t.elapsed().as_secs_f64(), merge);
    report();
}
