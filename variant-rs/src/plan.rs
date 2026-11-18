//! Shard planner (`variant-rs plan`): cuts the genome into load-balanced shard files for `--shard`.
//!
//! Cost model: the BAM index. A BAI's linear index holds, for every 16 kb window, the smallest file
//! offset of the reads overlapping it, so the compressed bytes between consecutive windows estimate how
//! much read data a window holds, without reading the BAM. A small per-base term keeps read-free sequence
//! from being free. HaplotypeCaller time also depends on local complexity, which reads alone do not
//! capture; the scheduler compensates by running many more shards than workers, largest first (LPT).
//!
//! Partition: contigs in dictionary order. A contig costing well above the target per-shard cost is
//! split into near-equal-cost pieces; each cut is moved to the middle of an N run of >= 1 kb when one
//! lies near the ideal cut (no read or active region crosses an N gap, so the cut changes nothing but
//! the RNG domain), otherwise to the window with the fewest reads nearby. Smaller contigs are grouped,
//! in order, into shards of about the target cost. Shards are contiguous in dictionary order, so their
//! outputs concatenate into a sorted VCF/GVCF.
//!
//! Each shard is one GATK JVM's worth of work (`-L shardNNN.intervals`), so `variant-rs --shard` output
//! equals GATK run on the same shard files.

use crate::bam::Bai;
use crate::intervals::Dict;
use crate::reference::Fasta;
use std::io::Write;

pub const WIN: i64 = 16384;

/// Estimated cost of every 16 kb window of BAM reference `tid` (index units: compressed bytes).
pub fn window_costs(bai: &Bai, tid: usize, len: i64) -> Vec<f64> {
    let n = ((len + WIN - 1) / WIN) as usize;
    let Some(ri) = bai.refs.get(tid) else { return vec![0.0; n] };
    let end_off = ri.bins.values().flatten().map(|c| c.1 >> 16).max().unwrap_or(0);
    // linear[w] (compressed part); empty windows inherit the next non-empty offset
    let mut off = vec![end_off; n + 1];
    let mut next = end_off;
    for w in (0..n).rev() {
        let v = ri.linear.get(w).copied().unwrap_or(0) >> 16;
        if v != 0 {
            next = v;
        }
        off[w] = next;
    }
    (0..n).map(|w| off[w + 1].saturating_sub(off[w]) as f64).collect()
}

/// Per-window cost with a per-base floor (`floor` = fraction of the mean cost per window).
fn with_floor(c: &[f64], mean_per_win: f64, floor: f64) -> Vec<f64> {
    c.iter().map(|&x| x + floor * mean_per_win).collect()
}

/// Cost of 1-based inclusive `[s, e]` from window costs (pro rata at the edges).
pub fn interval_cost(costs: &[f64], s: i64, e: i64) -> f64 {
    if costs.is_empty() || e < s {
        return 0.0;
    }
    let (a, b) = ((s - 1) / WIN, (e - 1) / WIN);
    let mut t = 0.0;
    for w in a..=b {
        let ws = w * WIN + 1;
        let we = ws + WIN - 1;
        let ov = (e.min(we) - s.max(ws) + 1).max(0) as f64;
        t += costs.get(w as usize).copied().unwrap_or(0.0) * ov / WIN as f64;
    }
    t
}

/// N runs of at least `min_len` bp, 1-based inclusive.
fn n_runs(seq: &[u8], min_len: i64) -> Vec<(i64, i64)> {
    let mut v = Vec::new();
    let mut i = 0usize;
    while i < seq.len() {
        if seq[i] == b'N' {
            let s = i;
            while i < seq.len() && seq[i] == b'N' {
                i += 1;
            }
            if (i - s) as i64 >= min_len {
                v.push((s as i64 + 1, i as i64));
            }
        } else {
            i += 1;
        }
    }
    v
}

pub struct PlanArgs {
    pub bam: String,
    pub fasta: String,
    pub out_dir: String,
    pub pieces: usize,
    pub floor: f64,
    pub min_gap: i64,
}

fn usage() -> ! {
    eprintln!(
        "usage: variant-rs plan -I in.bam -R ref.fa -o OUT_DIR [--pieces N (96)] [--min-gap BP (1000)] [--base-cost F (0.05)]\n\
         writes OUT_DIR/shardNNN.intervals (contiguous, in dictionary order) and OUT_DIR/plan.tsv; run them with\n\
         variant-rs ... $(for f in OUT_DIR/shard*.intervals; do echo --shard $f; done) -O out.vcf.gz"
    );
    std::process::exit(2)
}

pub fn parse(args: &[String]) -> PlanArgs {
    let mut p = PlanArgs { bam: String::new(), fasta: String::new(), out_dir: String::new(), pieces: 96, floor: 0.05, min_gap: 1000 };
    let mut i = 0;
    while i < args.len() {
        let v = || args.get(i + 1).cloned().unwrap_or_else(|| usage());
        match args[i].as_str() {
            "-I" => p.bam = v(),
            "-R" => p.fasta = v(),
            "-o" | "-O" => p.out_dir = v(),
            "--pieces" => p.pieces = v().parse().unwrap(),
            "--base-cost" => p.floor = v().parse().unwrap(),
            "--min-gap" => p.min_gap = v().parse().unwrap(),
            _ => usage(),
        }
        i += 2;
    }
    if p.bam.is_empty() || p.fasta.is_empty() || p.out_dir.is_empty() || p.pieces == 0 {
        usage();
    }
    p
}

/// The BAM's `.bai` (`x.bam.bai` or `x.bai`), if there is one.
pub fn find_bai(bam: &str) -> Option<String> {
    [format!("{bam}.bai"), bam.strip_suffix(".bam").map(|s| format!("{s}.bai")).unwrap_or_default()]
        .into_iter()
        .find(|p| std::path::Path::new(p).is_file())
}

/// Cost-model input shared with the scheduler: per BAM tid, window costs including the floor.
pub struct CostModel {
    pub costs: Vec<Vec<f64>>,
}

impl CostModel {
    /// `raw`: per reference, read data per 16 kb window (`Input::window_costs`).
    pub fn new(raw: Vec<Vec<f64>>, floor: f64) -> CostModel {
        let (tot, nw): (f64, usize) = raw.iter().fold((0.0, 0), |(a, b), v| (a + v.iter().sum::<f64>(), b + v.len()));
        let mean = if nw > 0 { tot / nw as f64 } else { 0.0 };
        CostModel { costs: raw.iter().map(|v| with_floor(v, mean, floor)).collect() }
    }
    pub fn cost(&self, tid: usize, intervals: &[(i64, i64)]) -> f64 {
        intervals.iter().map(|&(s, e)| interval_cost(&self.costs[tid], s, e)).sum()
    }
}

pub fn run(args: &[String]) {
    let a = parse(args);
    let (input, header) = crate::input::Input::open(&a.bam, &a.fasta);
    let ref_dict = Dict::from_reference(&a.fasta);
    let fasta = Fasta::open(&a.fasta);
    let model = CostModel::new(input.window_costs(&header.ref_lens), a.floor);
    // contigs: the BAM's, in reference-dictionary order (as variant-rs without -L)
    let mut contigs: Vec<(usize, usize)> = header
        .ref_names
        .iter()
        .enumerate()
        .map(|(t, n)| (*ref_dict.index.get(n).unwrap_or_else(|| panic!("BAM contig {n} not in the reference")), t))
        .collect();
    contigs.sort();
    let total: f64 = contigs.iter().map(|&(_, t)| model.costs[t].iter().sum::<f64>()).sum();
    let target = total / a.pieces as f64;

    // shards: list of (contig name, start, end)
    let mut shards: Vec<Vec<(String, i64, i64)>> = Vec::new();
    let mut costs: Vec<f64> = Vec::new();
    let mut group: Vec<(String, i64, i64)> = Vec::new();
    let mut gcost = 0.0;
    let mut n_split = 0;
    let mut n_gap_cuts = 0;
    for &(rt, t) in &contigs {
        let name = &ref_dict.names[rt];
        let len = ref_dict.lens[rt];
        let wc = &model.costs[t];
        let c: f64 = wc.iter().sum();
        // every piece at least one 16 kb window long (a short, very deep contig such as chrM gets fewer pieces)
        let m = ((c / target).round() as usize).min((len / WIN) as usize);
        if m >= 2 {
            if !group.is_empty() {
                shards.push(std::mem::take(&mut group));
                costs.push(gcost);
                gcost = 0.0;
            }
            n_split += 1;
            let (seq, _) = fasta.fetch(name, 1, len);
            let gaps = n_runs(&seq, a.min_gap);
            drop(seq);
            // cumulative cost at window boundaries
            let mut cum = vec![0.0; wc.len() + 1];
            for (w, &x) in wc.iter().enumerate() {
                cum[w + 1] = cum[w] + x;
            }
            let pos_at = |q: f64| -> i64 {
                let w = cum.partition_point(|&x| x < q).max(1) - 1;
                let frac = if wc[w.min(wc.len() - 1)] > 0.0 { (q - cum[w]) / wc[w.min(wc.len() - 1)] } else { 0.0 };
                ((w as f64 + frac.clamp(0.0, 1.0)) * WIN as f64) as i64
            };
            let mut cuts: Vec<i64> = Vec::new(); // last base of each piece except the final one
            for i in 1..m {
                let ideal = pos_at(c * i as f64 / m as f64).clamp(1, len - 1);
                let radius = (len / m as i64) / 8; // allow +/- 1/8 of a piece length
                // this cut in [lo, hi]: a window after the previous cut, a window for each piece still to come
                let lo = cuts.last().copied().unwrap_or(0) + WIN;
                let hi = len - (m - i) as i64 * WIN;
                let gap = gaps
                    .iter()
                    .map(|&(s, e)| (s + e) / 2)
                    .filter(|&mid| (mid - ideal).abs() <= radius && mid > lo && mid < len - WIN && mid <= hi)
                    .min_by_key(|&mid| (mid - ideal).abs());
                let cut = match gap {
                    Some(mid) => {
                        n_gap_cuts += 1;
                        mid
                    }
                    None => {
                        // least-read 16 kb window within +/- radius/4; cut in its middle
                        let r = (radius / 4).max(WIN);
                        let (w0, w1) = (((ideal - r).max(lo) / WIN) as usize, (((ideal + r) / WIN) as usize).min(wc.len() - 1));
                        let best = (w0..=w1.max(w0)).min_by(|&x, &y| wc[x].partial_cmp(&wc[y]).unwrap().then((x as i64 * WIN - ideal).abs().cmp(&(y as i64 * WIN - ideal).abs()))).unwrap();
                        (best as i64 * WIN + WIN / 2).clamp(lo, hi)
                    }
                };
                cuts.push(cut);
            }
            let mut s = 1;
            for &cut in cuts.iter().chain(std::iter::once(&len)) {
                shards.push(vec![(name.clone(), s, cut)]);
                costs.push(interval_cost(wc, s, cut));
                s = cut + 1;
            }
        } else {
            group.push((name.clone(), 1, len));
            gcost += c;
            if gcost >= target * 0.9 {
                shards.push(std::mem::take(&mut group));
                costs.push(gcost);
                gcost = 0.0;
            }
        }
    }
    if !group.is_empty() {
        shards.push(group);
        costs.push(gcost);
    }

    std::fs::create_dir_all(&a.out_dir).unwrap();
    let mut tsv = std::fs::File::create(format!("{}/plan.tsv", a.out_dir)).unwrap();
    writeln!(tsv, "shard\tcontigs\tfirst\tlast\tbases\test_cost_fraction").unwrap();
    for (i, (sh, &c)) in shards.iter().zip(&costs).enumerate() {
        let name = format!("shard{i:03}");
        let mut f = std::io::BufWriter::new(std::fs::File::create(format!("{}/{name}.intervals", a.out_dir)).unwrap());
        for (n, s, e) in sh {
            writeln!(f, "{n}:{s}-{e}").unwrap();
        }
        let bases: i64 = sh.iter().map(|(_, s, e)| e - s + 1).sum();
        let (fa, la) = (&sh[0], sh.last().unwrap());
        writeln!(tsv, "{name}\t{}\t{}:{}-{}\t{}:{}-{}\t{bases}\t{:.5}", sh.len(), fa.0, fa.1, fa.2, la.0, la.1, la.2, c / total).unwrap();
    }
    let maxc = costs.iter().cloned().fold(0.0, f64::max);
    eprintln!(
        "{} shards from {} contigs ({} split, {} cuts in N gaps); target cost {:.4} of total, largest shard {:.4}",
        shards.len(),
        contigs.len(),
        n_split,
        n_gap_cuts,
        1.0 / a.pieces as f64,
        maxc / total
    );
}
