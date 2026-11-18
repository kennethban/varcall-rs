//! Replays a Smith-Waterman call trace recorded with `VARIANT_RS_SW_DUMP` through the reference scalar aligners in
//! `old_sw.rs` (the validated port of the Java aligner) and the current ones, checks every result is identical, and times both.
//!
//! usage: swbench <trace> [threads]   (threads > 1: the trace is split round-robin, each thread replays its share)
#[allow(dead_code)]
mod old_sw;

use hc_assembly::sw::{self, Overhang, FixedRefAligner, SwAligner, SwParameters};
use std::time::Instant;

struct SCall {
    p: SwParameters,
    strat: Overhang,
    reference: Vec<u8>,
    alt: Vec<u8>,
}
struct Block {
    p: SwParameters,
    reference: Vec<u8>,
    batches: Vec<Vec<Vec<u8>>>,
}

fn params(s: &str) -> SwParameters {
    let v: Vec<i32> = s.split(',').map(|x| x.parse().unwrap()).collect();
    SwParameters { match_value: v[0], mismatch_penalty: v[1], gap_open_penalty: v[2], gap_extend_penalty: v[3] }
}

fn load(path: &str) -> (Vec<SCall>, Vec<Block>) {
    let data = std::fs::read(path).unwrap();
    let (mut s, mut b) = (Vec::new(), Vec::<Block>::new());
    for line in data.split(|&c| c == b'\n').filter(|l| !l.is_empty()) {
        let t: Vec<&[u8]> = line.split(|&c| c == b'\t').collect();
        let st = |x: &[u8]| std::str::from_utf8(x).unwrap().to_string();
        match t[0] {
            b"S" => s.push(SCall {
                p: params(&st(t[1])),
                strat: match t[2] {
                    b"SoftClip" => Overhang::SoftClip,
                    b"Indel" => Overhang::Indel,
                    b"LeadingIndel" => Overhang::LeadingIndel,
                    _ => Overhang::Ignore,
                },
                reference: t[3].to_vec(),
                alt: t[4].to_vec(),
            }),
            b"R" => b.push(Block { p: params(&st(t[1])), reference: t[2].to_vec(), batches: Vec::new() }),
            b"B" => b.last_mut().unwrap().batches.push(Vec::new()),
            b"A" => b.last_mut().unwrap().batches.last_mut().unwrap().push(t[1].to_vec()),
            _ => panic!("bad line"),
        }
    }
    (s, b)
}

type Res = Vec<(Vec<(u8, u32)>, i32)>;

fn to_old(p: &SwParameters) -> old_sw::SwParameters {
    old_sw::SwParameters {
        match_value: p.match_value,
        mismatch_penalty: p.mismatch_penalty,
        gap_open_penalty: p.gap_open_penalty,
        gap_extend_penalty: p.gap_extend_penalty,
    }
}
fn old_strat(s: Overhang) -> old_sw::Overhang {
    match s {
        Overhang::SoftClip => old_sw::Overhang::SoftClip,
        Overhang::Indel => old_sw::Overhang::Indel,
        Overhang::LeadingIndel => old_sw::Overhang::LeadingIndel,
        Overhang::Ignore => old_sw::Overhang::Ignore,
    }
}

fn s_old(calls: &[&SCall]) -> Res {
    let mut a = old_sw::SwAligner::new();
    calls.iter().map(|c| {
        let r = a.align(&c.reference, &c.alt, &to_old(&c.p), old_strat(c.strat));
        (r.cigar, r.offset)
    }).collect()
}
fn s_new(calls: &[&SCall]) -> Res {
    let mut a = SwAligner::new();
    calls.iter().map(|c| {
        let r = a.align(&c.reference, &c.alt, &c.p, c.strat);
        (r.cigar, r.offset)
    }).collect()
}
fn h_old(blocks: &[&Block]) -> Res {
    let mut a = old_sw::PrefixCachedSw::new();
    let mut out = Vec::new();
    for b in blocks {
        a.set_reference(&b.reference, &to_old(&b.p));
        for batch in &b.batches {
            for alt in batch {
                let r = a.align(alt);
                out.push((r.cigar, r.offset));
            }
        }
    }
    out
}
fn h_new(blocks: &[&Block]) -> Res {
    let mut a = FixedRefAligner::new();
    let mut out = Vec::new();
    for b in blocks {
        a.set_reference(&b.reference, &b.p);
        for batch in &b.batches {
            let alts: Vec<&[u8]> = batch.iter().map(|x| &x[..]).collect();
            for r in a.align_batch(&alts) {
                out.push((r.cigar, r.offset));
            }
        }
    }
    out
}

/// Runs `f` on `threads` round-robin shares of `items`; returns (seconds, results in share order).
fn timed<T: Sync>(items: &[T], threads: usize, f: impl Fn(&[&T]) -> Res + Sync) -> (f64, Vec<Res>) {
    let shares: Vec<Vec<&T>> = (0..threads).map(|t| items.iter().skip(t).step_by(threads).collect()).collect();
    let t0 = Instant::now();
    let out = std::thread::scope(|s| {
        let hs: Vec<_> = shares.iter().map(|sh| s.spawn(|| f(sh))).collect();
        hs.into_iter().map(|h| h.join().unwrap()).collect::<Vec<_>>()
    });
    (t0.elapsed().as_secs_f64(), out)
}

fn compare(name: &str, a: &[Res], b: &[Res]) {
    let (mut n, mut bad) = (0usize, 0usize);
    for (x, y) in a.iter().zip(b) {
        assert_eq!(x.len(), y.len());
        for (p, q) in x.iter().zip(y) {
            n += 1;
            if p != q {
                bad += 1;
                if bad <= 5 {
                    eprintln!("{name} MISMATCH: old {:?} new {:?}", p, q);
                }
            }
        }
    }
    println!("{name}: {n} alignments, {}", if bad == 0 { "all identical".to_string() } else { format!("{bad} DIFFERENT") });
}

/// `swbench --random N`: N random pairs per strategy and parameter set (lengths 1-60 over a small alphabet, with
/// related pairs carrying indels), old vs new.
fn random_test(n: usize) {
    let mut seed = 0x2545F4914F6CDD1Du64;
    let mut rnd = move |k: u64| {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed % k
    };
    let ps = [sw::STANDARD_NGS, sw::NEW_SW_PARAMETERS, SwParameters { match_value: 10, mismatch_penalty: -15, gap_open_penalty: -30, gap_extend_penalty: -5 }];
    let strats = [Overhang::SoftClip, Overhang::Indel, Overhang::LeadingIndel, Overhang::Ignore];
    let (mut a_old, mut a_new) = (old_sw::SwAligner::new(), SwAligner::new());
    let (mut tot, mut bad) = (0, 0);
    for it in 0..n {
        let alpha: &[u8] = if it % 3 == 0 { b"AC" } else { b"ACGTN" };
        let la = 1 + rnd(60) as usize;
        let r: Vec<u8> = (0..la).map(|_| alpha[rnd(alpha.len() as u64) as usize]).collect();
        let a: Vec<u8> = if it % 2 == 0 {
            (0..1 + rnd(60) as usize).map(|_| alpha[rnd(alpha.len() as u64) as usize]).collect()
        } else {
            // a mutated copy of the reference: substitutions, insertions, deletions
            let mut a = Vec::new();
            let mut i = 0;
            while i < r.len() {
                match rnd(12) {
                    0 => a.push(alpha[rnd(alpha.len() as u64) as usize]),
                    1 => i += 1 + rnd(8) as usize,
                    2 => {
                        for _ in 0..1 + rnd(8) {
                            a.push(alpha[rnd(alpha.len() as u64) as usize]);
                        }
                    }
                    _ => {
                        a.push(r[i]);
                        i += 1;
                    }
                }
            }
            if a.is_empty() {
                a.push(b'A');
            }
            a
        };
        for p in &ps {
            for &st in &strats {
                for (x, y) in [(&r, &a), (&a, &r)] {
                    tot += 1;
                    let o = a_old.align(x, y, &to_old(p), old_strat(st));
                    let w = a_new.align(x, y, p, st);
                    if (o.cigar.clone(), o.offset) != (w.cigar.clone(), w.offset) {
                        bad += 1;
                        if bad <= 5 {
                            eprintln!("MISMATCH {:?} {:?} ref={} alt={}: old {:?} {} new {:?} {}", p, st, String::from_utf8_lossy(x), String::from_utf8_lossy(y), o.cigar, o.offset, w.cigar, w.offset);
                        }
                    }
                }
            }
        }
    }
    println!("random: {tot} alignments, {}", if bad == 0 { "all identical".to_string() } else { format!("{bad} DIFFERENT") });
    // FixedRefAligner batches: alternates derived from a few stems (shared prefixes, some prefixes of others),
    // sorted or not, several batches per reference; each result vs the old aligner from scratch
    let (mut tot, mut bad) = (0, 0);
    let mut fr = FixedRefAligner::new();
    for it in 0..n / 10 {
        let alpha: &[u8] = if it % 3 == 0 { b"AC" } else { b"ACGT" };
        let r: Vec<u8> = (0..1 + rnd(80) as usize).map(|_| alpha[rnd(alpha.len() as u64) as usize]).collect();
        let p = ps[it % ps.len()];
        fr.set_reference(&r, &p);
        for _ in 0..1 + rnd(3) {
            let stems: Vec<Vec<u8>> = (0..1 + rnd(3)).map(|_| (0..rnd(40) as usize).map(|_| alpha[rnd(alpha.len() as u64) as usize]).collect()).collect();
            let mut batch: Vec<Vec<u8>> = (0..1 + rnd(12))
                .map(|_| {
                    let st = &stems[rnd(stems.len() as u64) as usize];
                    let mut a = st[..rnd(st.len() as u64 + 1) as usize].to_vec();
                    for _ in 0..rnd(30) {
                        a.push(alpha[rnd(alpha.len() as u64) as usize]);
                    }
                    if a.is_empty() {
                        a.push(b'G');
                    }
                    a
                })
                .collect();
            if rnd(4) != 0 {
                batch.sort();
            }
            let alts: Vec<&[u8]> = batch.iter().map(|x| &x[..]).collect();
            for (a, w) in alts.iter().zip(fr.align_batch(&alts)) {
                tot += 1;
                let o = a_old.align(&r, a, &to_old(&p), old_sw::Overhang::SoftClip);
                if (o.cigar.clone(), o.offset) != (w.cigar.clone(), w.offset) {
                    bad += 1;
                    if bad <= 5 {
                        eprintln!("BATCH MISMATCH ref={} alt={}: old {:?} {} new {:?} {}", String::from_utf8_lossy(&r), String::from_utf8_lossy(a), o.cigar, o.offset, w.cigar, w.offset);
                    }
                }
            }
        }
    }
    println!("random FixedRefAligner batches: {tot} alignments, {}", if bad == 0 { "all identical".to_string() } else { format!("{bad} DIFFERENT") });
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args[1] == "--random" {
        return random_test(args[2].parse().unwrap());
    }
    let threads: usize = args.get(2).map(|x| x.parse().unwrap()).unwrap_or(1);
    let t0 = Instant::now();
    let (mut s, mut b) = load(&args[1]);
    // SWB_EVERY=k: replay every k-th reference block and call only (quicker experiments)
    if let Some(k) = std::env::var("SWB_EVERY").ok().map(|x| x.parse::<usize>().unwrap()) {
        s = s.into_iter().step_by(k).collect();
        b = b.into_iter().step_by(k).collect();
    }
    println!("loaded {} SwAligner calls, {} references ({:.1} s)", s.len(), b.len(), t0.elapsed().as_secs_f64());
    let _ = sw::STANDARD_NGS;
    let (to, ro) = timed(&b, threads, h_old);
    let (tn, rn) = timed(&b, threads, h_new);
    compare("haplotype-to-reference (FixedRefAligner)", &ro, &rn);
    println!("  {threads} thread(s): old {to:.2} s, new {tn:.2} s ({:.2}x)", to / tn);
    let (to, ro) = timed(&s, threads, s_old);
    let (tn, rn) = timed(&s, threads, s_new);
    compare("other calls (SwAligner)", &ro, &rn);
    println!("  {threads} thread(s): old {to:.2} s, new {tn:.2} s ({:.2}x)", to / tn);
    if std::env::var_os("VARIANT_RS_ASM_PROFILE").is_some() {
        println!("{}", hc_assembly::prof::hits_report());
    }
}
