//! Bitwise comparison of the CUDA PairHMM against the exact CPU port of GKL AVX-512
//! (`cpu::compute`, which is bit-identical to GATK's default PairHMM).
//!
//! usage: hmmcmp [n_groups] [seed]
//! Synthetic regions: haplotypes 100–600 bp, reads 30–300 bp sampled from a haplotype with
//! substitutions/indels and realistic qualities; some reads are random (to force the
//! double-precision fallback).
use hc_likelihood::{cpu, HmmGroup, HmmRead};

extern "C" {
    fn cuda_pairhmm_debug_cells(rb: *const u8, rq: *const u8, ri: *const u8, rd: *const u8, rc: *const u8, r: i32, hb: *const u8, h: i32, cells: *mut f32, raw: *mut f32) -> i32;
}

fn cell_diff(r: &HmmRead, h: &[u8]) {
    let (cc, craw) = cpu::debug_cells_f32(r, h);
    let mut gc = vec![0f32; cc.len()];
    let mut graw = 0f32;
    assert_eq!(unsafe { cuda_pairhmm_debug_cells(r.bases.as_ptr(), r.quals.as_ptr(), r.ins.as_ptr(), r.del.as_ptr(), r.gcp.as_ptr(), r.bases.len() as i32, h.as_ptr(), h.len() as i32, gc.as_mut_ptr(), &mut graw) }, 1);
    let w = h.len() + 1;
    let mut nd = 0;
    for i in 1..=r.bases.len() {
        for j in 1..=h.len() {
            let o = 3 * (i * w + j);
            for (k, nm) in ["M", "X", "Y"].iter().enumerate() {
                let (a, b) = (gc[o + k], cc[o + k]);
                if a.to_bits() != b.to_bits() {
                    if nd < 6 {
                        let ow = |ii: usize, jj: usize| 3 * (ii * w + jj);
                        let (d, u, l) = (ow(i - 1, j - 1), ow(i - 1, j), ow(i, j - 1));
                        println!("  first cells: {nm}[{i}][{j}] gpu {a:e} ({:#010x}) cpu {b:e} ({:#010x}); cpu diag M,X,Y {:e} {:e} {:e}; up M,X {:e} {:e}; left M,Y {:e} {:e}; q {} ins {} del {} gcp {} rb {} hb {}",
                            a.to_bits(), b.to_bits(), cc[d], cc[d+1], cc[d+2], cc[u], cc[u+1], cc[l], cc[l+2],
                            r.quals[i-1], r.ins[i-1], r.del[i-1], r.gcp[i-1], r.bases[i-1] as char, h[j-1] as char);
                    }
                    nd += 1;
                }
            }
        }
    }
    println!("  raw gpu {graw:e} ({:#010x}) cpu {craw:e} ({:#010x}); differing cell values {nd}", graw.to_bits(), craw.to_bits());
}

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
    let args: Vec<String> = std::env::args().collect();
    let n_groups: usize = args.get(1).map(|s| s.parse().unwrap()).unwrap_or(200);
    let seed: u64 = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(12345);
    let mut rng = Rng(seed | 1);
    const B: &[u8] = b"ACGT";
    let mut groups_owned: Vec<(Vec<OwnedRead>, Vec<Vec<u8>>)> = Vec::new();
    for _ in 0..n_groups {
        // HAPMIN/HAPSPAN: haplotype length range (default 100..600)
        let hmin: usize = std::env::var("HAPMIN").ok().and_then(|v| v.parse().ok()).unwrap_or(100);
        let hspan: u64 = std::env::var("HAPSPAN").ok().and_then(|v| v.parse().ok()).unwrap_or(500);
        let hl = hmin + rng.below(hspan) as usize;
        let base: Vec<u8> = (0..hl).map(|_| B[rng.below(4) as usize]).collect();
        let nh = 1 + rng.below(6) as usize;
        let mut haps = vec![base.clone()];
        for _ in 1..nh {
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
            if std::env::var("NO_N").is_err() && rng.below(10) == 0 {
                let p = rng.below(h.len() as u64) as usize;
                h[p] = b'N';
            }
            haps.push(h);
        }
        let nr = 1 + rng.below(40) as usize;
        let mut reads = Vec::new();
        for _ in 0..nr {
            let src = &haps[rng.below(nh as u64) as usize];
            let rl = (30 + rng.below(271) as usize).min(src.len());
            let st = rng.below((src.len() - rl + 1) as u64) as usize;
            let random = rng.below(25) == 0;
            let mut bases: Vec<u8> = if random { (0..rl).map(|_| B[rng.below(4) as usize]).collect() } else { src[st..st + rl].to_vec() };
            let err = rng.below(4);
            for b in bases.iter_mut() {
                if rng.below(1000) < err * 10 {
                    *b = B[rng.below(4) as usize];
                }
                if std::env::var("NO_N").is_err() && rng.below(2000) == 0 {
                    *b = b'N';
                }
            }
            let quals: Vec<u8> = (0..rl).map(|_| if rng.below(10) == 0 { 2 + rng.below(20) as u8 } else { 20 + rng.below(22) as u8 }).collect();
            let ins: Vec<u8> = (0..rl).map(|_| if rng.below(20) == 0 { 10 + rng.below(35) as u8 } else { 45 }).collect();
            let del = ins.clone();
            let wild = std::env::var("WILD").is_ok();
            let gcp: Vec<u8> = (0..rl).map(|_| if wild { rng.below(256) as u8 } else { 10 }).collect();
            let (quals, ins, del) = if wild {
                ((0..rl).map(|_| rng.below(256) as u8).collect(), (0..rl).map(|_| rng.below(256) as u8).collect(), (0..rl).map(|_| rng.below(256) as u8).collect())
            } else {
                (quals, ins.clone(), del)
            };
            reads.push(OwnedRead { bases, quals, ins, del, gcp });
        }
        groups_owned.push((reads, haps));
    }
    let groups: Vec<HmmGroup> = groups_owned
        .iter()
        .map(|(rs, hs)| HmmGroup {
            reads: rs.iter().map(|r| HmmRead { bases: &r.bases, quals: &r.quals, ins: &r.ins, del: &r.del, gcp: &r.gcp }).collect(),
            haps: hs.iter().map(|h| h.as_slice()).collect(),
        })
        .collect();
    assert!(hc_cudahmm::init(false));
    let gpu = hc_cudahmm::compute(&groups);
    let cpu = cpu::compute(&groups, cpu::Flavor::Avx512);
    let lim = (1e-28f64).log10() - (2f64.powi(120)).log10();
    let (mut n, mut same, mut nd, mut same_d, mut maxd) = (0usize, 0usize, 0usize, 0usize, 0f64);
    let mut shown = 0;
    for (gi, (g, c)) in gpu.iter().zip(&cpu).enumerate() {
        let nh = groups[gi].haps.len();
        for (k, (a, b)) in g.iter().zip(c).enumerate() {
            n += 1;
            let dbl = *b < lim + 1.0; // roughly: double-fallback pairs
            nd += dbl as usize;
            if a.to_bits() == b.to_bits() {
                same += 1;
                same_d += dbl as usize;
            } else {
                maxd = maxd.max((a - b).abs());
                if std::env::var("BISECT").is_ok() && shown < 30 {
                    let r = &groups[gi].reads[k / nh];
                    let h = groups[gi].haps[k % nh];
                    let mut first = None;
                    for len in 1..=r.bases.len() {
                        let rr = HmmRead { bases: &r.bases[..len], quals: &r.quals[..len], ins: &r.ins[..len], del: &r.del[..len], gcp: &r.gcp[..len] };
                        let g1 = [HmmGroup { reads: vec![rr], haps: vec![h] }];
                        let a1 = hc_cudahmm::compute(&g1)[0][0];
                        let b1 = cpu::compute(&g1, cpu::Flavor::Avx512)[0][0];
                        if a1.to_bits() != b1.to_bits() { first = Some(len); break; }
                    }
                    println!("pair g{gi} r{} h{}: R {} H {} first differing prefix row {:?}", k / nh, k % nh, r.bases.len(), h.len(), first);
                }
                if shown < 8 {
                    if std::env::var("CELLS").is_ok() { cell_diff(&groups[gi].reads[k / nh], groups[gi].haps[k % nh]); }
                    shown += 1;
                    println!("diff g{gi} r{} h{}: gpu {a:.17e} cpu {b:.17e} (readlen {} haplen {})", k / nh, k % nh, groups[gi].reads[k / nh].bases.len(), groups[gi].haps[k % nh].len());
                }
            }
        }
    }
    println!("pairs {n}: bit-identical {same} ({:.3}%), max |diff| {maxd:.3e}; low-likelihood (double path) {nd}, identical {same_d}", 100.0 * same as f64 / n as f64);
}
