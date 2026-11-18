//! CPU PairHMM throughput: scalar `cpu::compute` vs lane-batched `cpu::compute_batched` (bit-identical check).
//!
//! usage: hmmbench [n_groups] [reads_per_group] [haps_per_group] [seed]
use hc_likelihood::{cpu, HmmGroup, HmmRead};
use std::time::Instant;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

struct OwnedRead {
    bases: Vec<u8>,
    quals: Vec<u8>,
    ins: Vec<u8>,
    del: Vec<u8>,
    gcp: Vec<u8>,
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let arg = |i: usize, d: u64| a.get(i).map_or(d, |s| s.parse().unwrap());
    let (n_groups, n_reads, n_haps, seed) = (arg(1, 200) as usize, arg(2, 100) as usize, arg(3, 8) as usize, arg(4, 7));
    let mut rng = Rng(seed | 1);
    const B: &[u8] = b"ACGT";
    let mut owned: Vec<(Vec<OwnedRead>, Vec<Vec<u8>>)> = Vec::new();
    for _ in 0..n_groups {
        let hl = 300 + rng.below(200) as usize;
        let base: Vec<u8> = (0..hl).map(|_| B[rng.below(4) as usize]).collect();
        let mut haps = vec![base.clone()];
        for _ in 1..n_haps {
            let mut h = base.clone();
            for _ in 0..1 + rng.below(4) {
                let p = rng.below(h.len() as u64) as usize;
                match rng.below(3) {
                    0 => h[p] = B[rng.below(4) as usize],
                    1 => h.insert(p, B[rng.below(4) as usize]),
                    _ => {
                        h.remove(p);
                    }
                }
            }
            if std::env::var("NOISE").is_ok() && rng.below(5) == 0 {
                let p = rng.below(h.len() as u64) as usize;
                h[p] = b'N';
            }
            haps.push(h);
        }
        let mut reads = Vec::new();
        for _ in 0..n_reads {
            let src = &haps[rng.below(haps.len() as u64) as usize];
            let rl = if rng.below(5) == 0 { 60 + rng.below(90) as usize } else { 150 }.min(src.len());
            let st = rng.below((src.len() - rl + 1) as u64) as usize;
            let noise = std::env::var("NOISE").is_ok();
            let random = noise && rng.below(25) == 0;
            let mut bases: Vec<u8> = if random { (0..rl).map(|_| B[rng.below(4) as usize]).collect() } else { src[st..st + rl].to_vec() };
            for b in bases.iter_mut() {
                if rng.below(100) == 0 {
                    *b = B[rng.below(4) as usize];
                }
                if noise && rng.below(500) == 0 {
                    *b = b'N';
                }
            }
            let wild = noise && rng.below(10) == 0;
            let q = |rng: &mut Rng, d: u8| if wild { rng.below(256) as u8 } else { d };
            let quals: Vec<u8> = (0..rl).map(|_| { let d = if rng.below(10) == 0 { 2 + rng.below(20) as u8 } else { 20 + rng.below(22) as u8 }; q(&mut rng, d) }).collect();
            let ins: Vec<u8> = (0..rl).map(|_| { let d = if rng.below(20) == 0 { 10 + rng.below(35) as u8 } else { 45 }; q(&mut rng, d) }).collect();
            let del: Vec<u8> = (0..rl).map(|k| if wild { rng.below(256) as u8 } else { ins[k] }).collect();
            let gcp: Vec<u8> = (0..rl).map(|_| q(&mut rng, 10)).collect();
            reads.push(OwnedRead { bases, quals, ins, del, gcp });
        }
        owned.push((reads, haps));
    }
    let groups: Vec<HmmGroup> = owned
        .iter()
        .map(|(rs, hs)| HmmGroup {
            reads: rs.iter().map(|r| HmmRead { bases: &r.bases, quals: &r.quals, ins: &r.ins, del: &r.del, gcp: &r.gcp }).collect(),
            haps: hs.iter().map(|h| h.as_slice()).collect(),
        })
        .collect();
    let cells: f64 = groups.iter().map(|g| g.reads.iter().map(|r| r.bases.len()).sum::<usize>() as f64 * g.haps.iter().map(|h| h.len()).sum::<usize>() as f64).sum();
    let t = Instant::now();
    let s = cpu::compute(&groups, cpu::Flavor::Avx512);
    let ts = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let b = cpu::compute_batched(&groups, cpu::Flavor::Avx512);
    let tb = t.elapsed().as_secs_f64();
    let same = s.iter().zip(&b).all(|(x, y)| x.iter().zip(y).all(|(p, q)| p.to_bits() == q.to_bits()));
    println!("{:.2} Gcells; scalar {:.2}s ({:.3} Gcells/s); batched {:.2}s ({:.3} Gcells/s); speedup {:.1}x; bit-identical {}", cells / 1e9, ts, cells / 1e9 / ts, tb, cells / 1e9 / tb, ts / tb, same);
}
