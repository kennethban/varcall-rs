//! Stage S3: `SomaticGenotypingEngine.callMutations` (tumor-normal, default arguments) with the
//! `StandardMutectAnnotation` group, `phaseCalls`, the event-count annotations (ECNT / ECNTH) and the
//! VCF data line as htsjdk writes it.
//!
//! Input is what stage S2 produces for a region: the trimmed haplotypes and, per sample (sorted sample
//! order), the realigned reads with their natural-log likelihoods.

use super::germline::Germline;
use super::region::SampleLikelihoods;
use super::somatic;
use hc_assembly::{Haplotype, Read};
use hc_genotype::annotate::{base_quality_at, read_position};
use hc_genotype::events::{create_allele_mapper, event_start_positions, merged_vc_at_with, Event, EventMap};
use hc_genotype::javamath::format_vcf_double;
use hc_genotype::likelihoods::{best_allele_t, BestAllele};
use std::collections::{HashMap, HashSet};
use std::f64::consts::{LN_10, LOG10_E};

/// `AlleleLikelihoods.NATURAL_LOG_INFORMATIVE_THRESHOLD`
const LN_INFORMATIVE: f64 = 0.2 * LN_10;
/// `M2ArgumentCollection`: normal LOD 2.2, informative read overlap margin 2, global read mismapping rate Q45
/// (the emission LOD is in `M2Opts`).
const NORMAL_LOD: f64 = 2.2;
const OVERLAP_MARGIN: i64 = 2;
const MISMAPPING_QUAL: f64 = 45.0;
/// `OrientationBiasReadCounts.MINIMUM_BASE_QUALITY`
const F1R2_MIN_BQ: u8 = 20;
/// `StrandBiasUtils.MIN_COUNT`
const AS_SB_MIN_COUNT: i32 = 2;

pub struct CallCtx<'a> {
    pub contig: &'a str,
    pub contig_len: i64,
    /// per sorted sample: true for normal samples
    pub is_normal: &'a [bool],
    pub germline: Option<&'a Germline>,
    pub pon: Option<&'a super::germline::Pon>,
    /// `af-of-alleles-not-in-resource` default (1e-6 with a normal, 5e-8 tumor-only)
    pub default_af: f64,
    pub opts: super::M2Opts,
    /// `--bam-output`: the contig's reference index in the reads' header
    pub bamout_tid: Option<i32>,
}

/// What `callMutations` needs from a region.
pub struct RegionCall<'a> {
    pub haplotypes: &'a [Haplotype],
    pub samples: &'a [SampleLikelihoods],
    pub ref_padded: &'a [u8],
    pub ref_padded_start: i64,
    /// `regionForGenotyping.getSpan()`
    pub span: (i64, i64),
    /// reference bases of the original region's padded span (the walker's `ReferenceContext`)
    pub ref_window: &'a [u8],
    pub ref_window_start: i64,
}

#[derive(Clone, Debug, Default)]
pub struct M2Genotype {
    /// allele indices (into the call's alleles)
    pub gt: Vec<usize>,
    pub ad: Vec<i32>,
    pub af: Vec<f64>,
    pub dp: i32,
    pub f1r2: Vec<i32>,
    pub f2r1: Vec<i32>,
    pub fad: Vec<i32>,
    pub sb: [i32; 4],
    /// (PGT, PID, PS)
    pub phase: Option<(&'static str, String, i64)>,
}

#[derive(Clone, Debug)]
pub struct M2Call {
    pub pos: i64,
    /// trimmed alleles, reference first
    pub alleles: Vec<Vec<u8>>,
    /// INFO fields except ECNT / ECNTH, already formatted
    pub info: Vec<(&'static str, String)>,
    /// ECNTH values (None: no assembled event found -> ".")
    pub ecnth: Option<Vec<i32>>,
    pub ecnt: usize,
    pub genotypes: Vec<M2Genotype>,
}

#[derive(Clone)]
struct Fragment {
    /// all reads of the name group (evidence order): their likelihoods are summed
    group: Vec<usize>,
    /// `Fragment.getReads()` (at most two)
    reads: Vec<usize>,
    start: i64,
    end: i64,
}

/// `AlleleLikelihoods.groupEvidence(GATKRead::getName, Fragment::createAndAvoidFailure)` for one sample.
fn group_fragments(reads: &[Read]) -> Vec<Fragment> {
    let names: Vec<&str> = reads.iter().map(|r| r.name.as_str()).collect();
    let mut by_name: HashMap<&str, Vec<usize>> = HashMap::new();
    for (i, n) in names.iter().enumerate() {
        by_name.entry(n).or_default().push(i);
    }
    somatic::java_hashmap_order(&names)
        .into_iter()
        .map(|first| {
            let group = by_name[names[first]].clone();
            let reads_used: Vec<usize> = if group.len() <= 2 {
                group.clone()
            } else {
                let non_supp: Vec<usize> = group.iter().copied().filter(|&i| reads[i].flag & (0x400 | 0x100 | 0x800) == 0).collect();
                if non_supp.len() > 2 {
                    non_supp[..2].to_vec()
                } else if non_supp.is_empty() {
                    group[..1].to_vec()
                } else {
                    non_supp
                }
            };
            let (start, end) = if reads_used.len() == 1 {
                let r = &reads[reads_used[0]];
                (r.start.min(r.end), r.start.max(r.end))
            } else {
                let (a, b) = (&reads[reads_used[0]], &reads[reads_used[1]]);
                let s = a.start.min(b.start);
                let e = a.end.max(b.end);
                (s.min(e), s.max(e))
            };
            Fragment { group, reads: reads_used, start, end }
        })
        .collect()
}

#[inline]
fn best(col: &[f64], pr: &[f64]) -> BestAllele {
    best_allele_t(col, Some(pr), LN_INFORMATIVE)
}

/// `MathUtils.median(int[])`: commons-math `Median` (LEGACY estimation) then `FastMath.round`.
fn median(v: &mut [i32]) -> i32 {
    v.sort_unstable();
    let n = v.len();
    let pos = 0.5 * (n as f64 + 1.0);
    let x = if pos < 1.0 {
        v[0] as f64
    } else if pos >= n as f64 {
        v[n - 1] as f64
    } else {
        let fpos = pos.floor();
        let ip = fpos as usize;
        let d = pos - fpos;
        let lo = v[ip - 1] as f64;
        let hi = v[ip] as f64;
        lo + d * (hi - lo)
    };
    hc_assembly::m2math::java_round(x) as i32
}

/// `GATKVariantContextUtils.trimAlleles(vc, true, true)` on allele bases (no symbolic alleles here).
fn trim_alleles(start: i64, alleles: &[Vec<u8>]) -> (i64, Vec<Vec<u8>>) {
    if alleles.len() <= 1 || alleles.iter().any(|a| a.len() == 1) {
        return (start, alleles.to_vec());
    }
    let n = alleles.len();
    let mut from = vec![0i64; n];
    let mut to: Vec<i64> = alleles.iter().map(|a| a.len() as i64).collect();
    let (mut start_shift, mut end_shift) = (0i64, 0i64);
    let mut min_size = (0..n).map(|i| to[i] - from[i]).min().unwrap();
    let last_same = |to: &[i64]| (0..n).all(|i| alleles[i][(to[i] - 1) as usize] == alleles[0][(to[0] - 1) as usize]);
    let first_same = |from: &[i64]| (0..n).all(|i| alleles[i][from[i] as usize] == alleles[0][from[0] as usize]);
    while min_size > 0 && last_same(&to) {
        to.iter_mut().for_each(|t| *t -= 1);
        min_size -= 1;
        end_shift += 1;
    }
    while min_size > 0 && first_same(&from) {
        from.iter_mut().for_each(|f| *f += 1);
        min_size -= 1;
        start_shift -= 1;
    }
    // maxShift = 0: the left-shift loop never runs
    let end_trim = end_shift;
    let start_trim = -start_shift;
    let empty = (0..n).any(|i| to[i] - from[i] == 0);
    let end_clip = if empty && start_trim == 0 { end_trim - 1 } else { end_trim };
    let start_clip = if empty && start_trim > 0 { start_trim - 1 } else { start_trim };
    if start_clip == 0 && end_clip == 0 {
        return (start, alleles.to_vec());
    }
    let out = alleles.iter().map(|a| a[start_clip as usize..a.len() - end_clip as usize].to_vec()).collect();
    (start + start_clip, out)
}

/// `GATKVariantContextUtils.findRepeatedSubstring`
fn find_repeated_substring(b: &[u8]) -> usize {
    let mut rep = 1;
    while rep <= b.len() {
        let unit = &b[..rep];
        let mut all = true;
        let mut s = rep;
        while s < b.len() {
            if s + rep > b.len() || &b[s..s + rep] != unit {
                all = false;
                break;
            }
            s += rep;
        }
        if all {
            return rep;
        }
        rep += 1;
    }
    rep
}

/// `findNumberOfRepetitions(repeatUnit, testString, leadingRepeats = true)`
fn leading_repeats(unit: &[u8], test: &[u8]) -> i32 {
    if test.is_empty() {
        return 0;
    }
    let ul = unit.len();
    let mut n = 0;
    let mut s = 0usize;
    while s + ul <= test.len() {
        if &test[s..s + ul] == unit {
            n += 1;
        } else {
            return n;
        }
        s += ul;
    }
    n
}

/// `GATKVariantContextUtils.getNumTandemRepeatUnits(refAllele, altAlleles, refBasesStartingAtVCWithoutPad)`
fn tandem_repeat_units(ref_allele: &[u8], alts: &[Vec<u8>], context: &[u8]) -> Option<(Vec<i32>, Vec<u8>)> {
    if alts.iter().all(|a| a.len() == ref_allele.len()) {
        return None;
    }
    let rb = &ref_allele[1..];
    let mut unit_out = Vec::new();
    let mut lengths = Vec::new();
    for a in alts {
        let ab = &a[1..];
        let long_b = if ab.len() > rb.len() { ab } else { rb };
        let ul = find_repeated_substring(long_b);
        let mut unit = long_b.to_vec();
        unit.resize(ul, 0); // Arrays.copyOf pads with zeros
        let in_ref = leading_repeats(&unit, rb);
        let mut r_ctx = rb.to_vec();
        r_ctx.extend_from_slice(context);
        let mut a_ctx = ab.to_vec();
        a_ctx.extend_from_slice(context);
        let c0 = leading_repeats(&unit, &r_ctx) - in_ref;
        let c1 = leading_repeats(&unit, &a_ctx) - in_ref;
        if c0 == 0 || c1 == 0 {
            return None;
        }
        if lengths.is_empty() {
            lengths.push(c0);
        }
        lengths.push(c1);
        unit_out = unit;
    }
    Some((lengths, unit_out))
}

fn join<T: ToString>(v: &[T]) -> String {
    v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
}

fn usable(r: &Read) -> bool {
    r.mapq != 0 && r.mapq != 255
}

/// `SomaticGenotypingEngine.callMutations` for one region.
/// What `--bam-output` needs from `callMutations`.
#[derive(Default)]
pub struct BamoutInfo {
    /// haplotypes in the order they are added to `calledHaplotypes` (repeats included)
    pub called_insertion: Vec<usize>,
    /// per sample, per read: `annotateReadLikelihoodsWithSupportedAlleles` entries (call position, allele index)
    pub xa: Vec<Vec<Vec<(i64, usize)>>>,
}

/// `SomaticGenotypingEngine.callMutations` for one region; also fills `bamout` when given.
pub fn call_mutations_ext(ctx: &CallCtx, rc: &RegionCall, mut bamout: Option<&mut BamoutInfo>) -> Vec<M2Call> {
    let haps = rc.haplotypes;
    let nh = haps.len();
    let ns = rc.samples.len();
    let has_normal = ctx.is_normal.iter().any(|&x| x);
    if let Some(b) = bamout.as_deref_mut() {
        b.xa = rc.samples.iter().map(|s| vec![Vec::new(); s.reads.len()]).collect();
    }
    let maps: Vec<EventMap> = haps.iter().map(|h| EventMap::from_haplotype(h, rc.ref_padded, rc.ref_padded_start, ctx.opts.max_mnp_distance as _)).collect();
    let starts: Vec<i64> = event_start_positions(&maps).into_iter().filter(|&l| rc.span.0 <= l && l <= rc.span.1).collect();

    // logReadLikelihoods.normalizeLikelihoods(qualToLogErrorProb(45), true)
    let cap = MISMAPPING_QUAL * (-LN_10 / 10.0);
    let mut lk: Vec<Vec<Vec<f64>>> = rc.samples.iter().map(|s| s.matrix.clone()).collect();
    if nh > 1 {
        for (s, m) in lk.iter_mut().enumerate() {
            for r in 0..rc.samples[s].reads.len() {
                let mut b = m[0][r];
                for h in 1..nh {
                    if m[h][r] > b {
                        b = m[h][r];
                    }
                }
                let worst = b + cap;
                for h in 0..nh {
                    if m[h][r] < worst {
                        m[h][r] = worst;
                    }
                }
            }
        }
    }
    // fragments
    let frags: Vec<Vec<Fragment>> = rc.samples.iter().map(|s| group_fragments(&s.reads)).collect();
    let flk: Vec<Vec<Vec<f64>>> = (0..ns)
        .map(|s| {
            (0..nh)
                .map(|h| {
                    frags[s]
                        .iter()
                        .map(|f| {
                            let mut v = 0.0;
                            for &r in &f.group {
                                v += lk[s][h][r];
                            }
                            v
                        })
                        .collect()
                })
                .collect()
        })
        .collect();

    let emission = ctx.opts.emission_log_odds;
    let normal_thr = NORMAL_LOD * LN_10;
    let mut potential: HashSet<Event> = HashSet::new();
    let mut calls: Vec<M2Call> = Vec::new();
    let mut called_haps = vec![false; nh];

    for loc in starts {
        let merged = match merged_vc_at_with(loc, &maps, rc.ref_padded, rc.ref_padded_start, false) {
            Some(m) => m,
            None => continue,
        };
        let (alleles, haps_of) = create_allele_mapper(&merged, loc, &maps);
        let n_merged = merged.alleles.len();
        let na = alleles.len();
        let lo = (loc - OVERLAP_MARGIN).max(1);
        let hi = (merged.end() + OVERLAP_MARGIN).min(ctx.contig_len);
        let marg = |vals: &[Vec<f64>], col: usize, a: usize| -> f64 {
            let mut m = f64::NEG_INFINITY;
            for &h in &haps_of[a] {
                let v = vals[h][col];
                if v > m || v.is_nan() {
                    m = v;
                }
            }
            m
        };
        // fragment allele likelihoods, retained evidence only: fm[s][a][k]
        let mut kept_f: Vec<Vec<usize>> = Vec::with_capacity(ns);
        let mut fm: Vec<Vec<Vec<f64>>> = Vec::with_capacity(ns);
        for s in 0..ns {
            let kept: Vec<usize> = (0..frags[s].len()).filter(|&f| frags[s][f].start <= hi && lo <= frags[s][f].end).collect();
            fm.push((0..na).map(|a| kept.iter().map(|&f| marg(&flk[s], f, a)).collect()).collect());
            kept_f.push(kept);
        }
        let combine = |want_normal: bool| -> (Vec<Vec<f64>>, usize) {
            let mut m = vec![Vec::new(); na];
            for s in 0..ns {
                if ctx.is_normal[s] == want_normal {
                    for a in 0..na {
                        m[a].extend_from_slice(&fm[s][a]);
                    }
                }
            }
            let n = m[0].len();
            (m, n)
        };
        let (tm, tn) = combine(false);
        let (nm, nn) = combine(true);
        let tumor_lods = somatic::somatic_log_odds(&tm, tn);
        let normal_lods = somatic::diploid_alt_log_odds(&nm, nn);
        let normal_artifact_lods = somatic::somatic_log_odds(&nm, nn);

        let tumor_alts: Vec<usize> = (1..n_merged).filter(|&a| tumor_lods[a] > emission).collect();
        // considerGermlineActive = genotypeGermlineSites && rng.nextDouble() < fraction (1.0: always)
        let germline_active = ctx.opts.genotype_germline;
        let to_genotype = tumor_alts.iter().filter(|&&a| !has_normal || germline_active || normal_lods[a] > normal_thr).count();
        for a in 1..n_merged {
            if tumor_lods[a] > emission && (!has_normal || normal_lods[a] > normal_thr) {
                potential.insert(Event::new(loc, &merged.alleles[0].bases, &alleles[a].bases));
            }
        }
        if to_genotype == 0 {
            continue;
        }
        let emit: Vec<usize> = std::iter::once(0).chain(tumor_alts.iter().copied()).collect();
        let emit_bases: Vec<Vec<u8>> = emit.iter().map(|&a| alleles[a].bases.clone()).collect();

        let mut info: Vec<(&'static str, String)> = Vec::new();
        // POPAF
        let afs = match ctx.germline {
            Some(g) => g.alt_allele_frequencies(ctx.contig, loc, &emit_bases, ctx.default_af),
            None => vec![ctx.default_af; emit.len() - 1],
        };
        info.push(("POPAF", afs.iter().map(|&x| format_vcf_double(-crate::jmath::log10(x))).collect::<Vec<_>>().join(",")));
        info.push(("TLOD", tumor_alts.iter().map(|&a| format_vcf_double(tumor_lods[a] * LOG10_E)).collect::<Vec<_>>().join(",")));
        if has_normal {
            info.push(("NALOD", tumor_alts.iter().map(|&a| format_vcf_double(normal_artifact_lods[a] * LOG10_E)).collect::<Vec<_>>().join(",")));
            info.push(("NLOD", tumor_alts.iter().map(|&a| format_vcf_double(normal_lods[a] * LOG10_E)).collect::<Vec<_>>().join(",")));
        }
        if ctx.pon.is_some_and(|p| p.starts_at(ctx.contig, loc)) {
            info.push(("PON", String::new()));
        }

        // genotypes (addGenotypes): GT and AF from the fragment likelihoods subset to the emitted alleles
        let ne = emit.len();
        let mut gts: Vec<M2Genotype> = (0..ns)
            .map(|s| {
                let sub: Vec<Vec<f64>> = emit.iter().map(|&a| fm[s][a].clone()).collect();
                let afm = somatic::allele_fractions_mean(&sub, kept_f[s].len());
                M2Genotype {
                    gt: if ctx.is_normal[s] { vec![0, 0] } else { (0..ne).collect() },
                    af: afm[1..].to_vec(),
                    ..Default::default()
                }
            })
            .collect();

        let (pos, trimmed) = trim_alleles(loc, &emit_bases);
        let vc_end = pos + trimmed[0].len() as i64 - 1;
        let pr: Vec<f64> = (0..ne).map(|e| if e == 0 { 1.0 } else { 0.0 }).collect();

        // read-level likelihoods for annotation: retained reads, emitted alleles
        let mut all_best: Vec<(usize, usize, Option<usize>)> = Vec::new(); // (sample, read, informative best allele)
        let mut col = vec![0f64; ne];
        let mut n_reads_total = 0usize;
        for s in 0..ns {
            let reads = &rc.samples[s].reads;
            let (mut ad, mut dp, mut sb) = (vec![0i32; ne], 0i32, [0i32; 4]);
            let mut n_kept = 0;
            for (r, rd) in reads.iter().enumerate() {
                if !(rd.start <= hi && lo <= rd.end) {
                    continue;
                }
                n_kept += 1;
                for (e, &a) in emit.iter().enumerate() {
                    col[e] = marg(&lk[s], r, a);
                }
                let b = best(&col, &pr);
                let inf = if b.is_informative() { b.allele } else { None };
                if let Some(e) = inf {
                    ad[e] += 1;
                    dp += 1;
                    let rev = rd.flag & 0x10 != 0;
                    sb[if e == 0 { 0 } else { 2 } + rev as usize] += 1;
                }
                all_best.push((s, r, inf));
            }
            n_reads_total += n_kept;
            let g = &mut gts[s];
            g.ad = ad;
            g.dp = if n_kept == 0 { 0 } else { dp };
            g.sb = if sb.iter().sum::<i32>() > 0 { sb } else { [0; 4] };
            // fragment-level: FAD, F1R2 / F2R1
            let (mut fad, mut f1r2, mut f2r1) = (vec![0i32; ne], vec![0i32; ne], vec![0i32; ne]);
            for (k, &f) in kept_f[s].iter().enumerate() {
                for (e, &a) in emit.iter().enumerate() {
                    col[e] = fm[s][a][k];
                }
                let b = best(&col, &pr);
                if !b.is_informative() {
                    continue;
                }
                let e = b.allele.unwrap();
                fad[e] += 1;
                if let Some(bo) = bamout.as_deref_mut() {
                    for &r in &frags[s][f].reads {
                        bo.xa[s][r].push((pos, e));
                    }
                }
                let first = &reads[frags[s][f].reads[0]];
                if usable(first) && base_quality_at(first, pos).map_or(0, |q| q as u8) >= F1R2_MIN_BQ {
                    let f2r1_read = (first.flag & 0x10 != 0) == (first.flag & 0x40 != 0);
                    if f2r1_read { f2r1[e] += 1 } else { f1r2[e] += 1 }
                }
            }
            g.fad = fad;
            g.f1r2 = f1r2;
            g.f2r1 = f2r1;
        }

        // INFO annotations over all samples' reads
        // AS_SB_TABLE: per sample, only samples with more than MIN_COUNT informative reads
        let mut as_sb = vec![[0i32; 2]; ne];
        for s in 0..ns {
            let mut t = vec![[0i32; 2]; ne];
            for &(ss, r, b) in &all_best {
                if ss == s {
                    if let Some(e) = b {
                        t[e][(rc.samples[s].reads[r].flag & 0x10 != 0) as usize] += 1;
                    }
                }
            }
            if t.iter().map(|x| x[0] + x[1]).sum::<i32>() > AS_SB_MIN_COUNT {
                for e in 0..ne {
                    as_sb[e][0] += t[e][0];
                    as_sb[e][1] += t[e][1];
                }
            }
        }
        info.push(("AS_SB_TABLE", as_sb.iter().map(|x| format!("{},{}", x[0], x[1])).collect::<Vec<_>>().join("|")));
        if n_reads_total > 0 {
            info.push(("DP", n_reads_total.to_string()));
        }
        // PerAlleleAnnotation: MBQ, MFRL, MMQ, MPOS
        let (mut mbq, mut mfrl, mut mmq, mut mpos) = (vec![Vec::new(); ne], vec![Vec::new(); ne], vec![Vec::new(); ne], vec![Vec::new(); ne]);
        for &(s, r, b) in &all_best {
            let e = match b {
                Some(e) => e,
                None => continue,
            };
            let rd = &rc.samples[s].reads[r];
            if !usable(rd) {
                continue;
            }
            if !(pos < rd.start || rd.end < pos) {
                if let Some(q) = base_quality_at(rd, pos) {
                    mbq[e].push(hc_assembly::m2math::java_round(q) as i32);
                }
                if let Some(p) = read_position(rd, pos, vc_end) {
                    mpos[e].push(p as i32);
                }
            }
            mfrl[e].push(rd.tlen.abs() as i32);
            mmq[e].push(rd.mapq as i32);
        }
        let agg = |v: &mut Vec<Vec<i32>>, empty: i32, include_ref: bool| -> String {
            let vals: Vec<i32> =
                (if include_ref { 0 } else { 1 }..ne).map(|e| if v[e].is_empty() { empty } else { median(&mut v[e]) }).collect();
            join(&vals)
        };
        info.push(("MBQ", agg(&mut mbq, 0, true)));
        info.push(("MFRL", agg(&mut mfrl, 0, true)));
        info.push(("MMQ", agg(&mut mmq, 60, true)));
        info.push(("MPOS", agg(&mut mpos, 50, false)));
        if ctx.opts.original_alignment {
            // OriginalAlignment: informative reads whose best allele is the alt with the largest TLOD (first on ties)
            // and whose OA tag names another contig
            let mut top = 0;
            for k in 1..tumor_alts.len() {
                if tumor_lods[tumor_alts[k]] > tumor_lods[tumor_alts[top]] {
                    top = k;
                }
            }
            let n = all_best.iter().filter(|&&(s, r, b)| b == Some(top + 1) && rc.samples[s].reads[r].oa_other_contig).count();
            info.push(("OCM", n.to_string()));
        }
        // TandemRepeat (indels only)
        if trimmed[1..].iter().all(|a| a.len() != trimmed[0].len()) {
            let idx = (pos + 1 - rc.ref_window_start) as usize;
            let context = if idx <= rc.ref_window.len() { &rc.ref_window[idx..] } else { &[][..] };
            if let Some((lens, unit)) = tandem_repeat_units(&trimmed[0], &trimmed[1..], context) {
                info.push(("RPA", join(&lens)));
                info.push(("RU", String::from_utf8_lossy(&unit).into_owned()));
                info.push(("STR", String::new()));
            }
        }
        for &a in &emit {
            for &h in &haps_of[a] {
                called_haps[h] = true;
                if let Some(bo) = bamout.as_deref_mut() {
                    bo.called_insertion.push(h);
                }
            }
        }
        calls.push(M2Call { pos, alleles: trimmed, info, ecnth: None, ecnt: 0, genotypes: gts });
    }

    // phaseCalls
    let phases = phase_calls(&calls, &called_haps, &maps);
    for (c, p) in calls.iter_mut().zip(phases) {
        if let Some((p01, pid, ps)) = p {
            for g in c.genotypes.iter_mut() {
                // phaseVC: a het genotype is reversed when the allele at the phased index is not site-specific
                let het = g.gt.iter().any(|&x| x != g.gt[0]);
                let idx = if p01 { 1 } else { 0 };
                if het && g.gt[idx] == 0 {
                    g.gt.reverse();
                }
                g.phase = Some((if p01 { "0|1" } else { "1|0" }, pid.clone(), ps));
            }
        }
    }

    // event counts: support of each haplotype = number of reads whose best haplotype it is
    let mut support = vec![0i32; nh];
    let hpr: Vec<f64> = haps.iter().map(|h| if h.is_ref { 1.0 } else { 0.0 }).collect();
    let mut hcol = vec![0f64; nh];
    for s in 0..ns {
        for r in 0..rc.samples[s].reads.len() {
            for h in 0..nh {
                hcol[h] = lk[s][h][r];
            }
            if let Some(h) = best(&hcol, &hpr).allele {
                support[h] += 1;
            }
        }
    }
    let mut haps_by_event: HashMap<&Event, Vec<usize>> = HashMap::new();
    for (h, m) in maps.iter().enumerate() {
        for e in m.events.values() {
            haps_by_event.entry(e).or_default().push(h);
        }
    }
    let ecnt = potential.len();
    for c in calls.iter_mut() {
        let mut v: Vec<i32> = Vec::new();
        for alt in &c.alleles[1..] {
            let ev = Event::new(c.pos, &c.alleles[0], alt);
            if let Some(hs) = haps_by_event.get(&ev) {
                // stable sort by decreasing support: the first haplotype with the maximum
                let mut bh = hs[0];
                for &h in hs {
                    if support[h] > support[bh] {
                        bh = h;
                    }
                }
                v.push(maps[bh].events.values().filter(|e| potential.contains(*e)).count() as i32);
            }
        }
        c.ecnth = if v.is_empty() { None } else { Some(v) };
        c.ecnt = ecnt;
    }
    calls
}

fn is_site_specific_alt(a: &[u8]) -> bool {
    a != b"*" && a.first() != Some(&b'<')
}

/// `AssemblyBasedCallerUtils.phaseCalls`: per call Some((PHASE_01, PID, PS)).
fn phase_calls(calls: &[M2Call], called_haps: &[bool], maps: &[EventMap]) -> Vec<Option<(bool, String, i64)>> {
    let nc = calls.len();
    let nh = called_haps.len();
    let hap_sets: Vec<Vec<bool>> = calls
        .iter()
        .map(|c| {
            let ss: Vec<&Vec<u8>> = c.alleles[1..].iter().filter(|a| is_site_specific_alt(a)).collect();
            if ss.len() != 1 {
                return vec![false; nh];
            }
            let alt = ss[0];
            (0..nh).map(|h| called_haps[h] && maps[h].events.values().any(|e| e.start == c.pos && &e.alt_bases == alt)).collect()
        })
        .collect();
    let count = |s: &[bool]| s.iter().filter(|&&x| x).count();
    let contains_all = |a: &[bool], b: &[bool]| a.iter().zip(b).all(|(&x, &y)| x || !y);
    let mut union = vec![false; nh];
    for s in &hap_sets {
        for h in 0..nh {
            union[h] |= s[h];
        }
    }
    let total = count(&union);
    let mut mapping: Vec<Option<(usize, bool)>> = vec![None; nc];
    let mut unique = 0usize;
    let mut out: Vec<Option<(bool, String, i64)>> = vec![None; nc];
    for i in 0..nc.saturating_sub(1) {
        let hw = &hap_sets[i];
        let nw = count(hw);
        if nw == 0 {
            continue;
        }
        let call_on_all = nw == total;
        let mut available = hw.clone();
        for j in i + 1..nc {
            let hc = &hap_sets[j];
            let ncomp = count(hc);
            if ncomp == 0 {
                continue;
            }
            let comp_on_all = ncomp == total;
            if (nw == ncomp && contains_all(hw, hc)) || (call_on_all && contains_all(&available, hc)) || comp_on_all {
                if mapping[i].is_none() {
                    if mapping[j].is_some() {
                        return out; // un-phasable: GATK clears the mapping
                    }
                    mapping[i] = Some((unique, true));
                    mapping[j] = Some((unique, true));
                    for h in 0..nh {
                        available[h] &= hc[h];
                    }
                    unique += 1;
                } else if mapping[j].is_none() {
                    mapping[j] = mapping[i];
                }
            } else if nw + ncomp == total && (0..nh).all(|h| !(hw[h] && hc[h])) {
                if mapping[i].is_none() {
                    if mapping[j].is_some() {
                        return out;
                    }
                    mapping[i] = Some((unique, true));
                    mapping[j] = Some((unique, false));
                    unique += 1;
                } else if mapping[j].is_none() {
                    let (g, p01) = mapping[i].unwrap();
                    mapping[j] = Some((g, !p01));
                }
            }
        }
    }
    for group in 0..unique {
        let idx: Vec<usize> = (0..nc).filter(|&k| mapping[k].map_or(false, |(g, _)| g == group)).collect();
        assert!(idx.len() >= 2, "phase group with fewer than 2 members");
        let first = &calls[idx[0]];
        let pid = format!("{}_{}_{}", first.pos, String::from_utf8_lossy(&first.alleles[0]), String::from_utf8_lossy(&first.alleles[1]));
        for &k in &idx {
            out[k] = Some((mapping[k].unwrap().1, pid.clone(), first.pos));
        }
    }
    out
}

/// The VCF data line (htsjdk `VCFEncoder`): INFO keys sorted, FORMAT keys GT first then sorted.
pub fn to_line(contig: &str, c: &M2Call) -> String {
    let mut info = c.info.clone();
    info.push(("ECNT", c.ecnt.to_string()));
    info.push(("ECNTH", c.ecnth.as_ref().map_or(".".to_string(), |v| join(v))));
    info.sort_by(|a, b| a.0.cmp(b.0));
    let info_s: Vec<String> = info.iter().map(|(k, v)| if *k == "STR" || *k == "PON" { k.to_string() } else { format!("{k}={v}") }).collect();
    let alts: Vec<String> = c.alleles[1..].iter().map(|a| String::from_utf8_lossy(a).into_owned()).collect();
    let phased = c.genotypes.iter().any(|g| g.phase.is_some());
    let mut keys = vec!["GT", "AD", "AF", "DP", "F1R2", "F2R1", "FAD"];
    if phased {
        keys.extend(["PGT", "PID", "PS"]);
    }
    keys.push("SB");
    let mut s = format!(
        "{}\t{}\t.\t{}\t{}\t.\t.\t{}\t{}",
        contig,
        c.pos,
        String::from_utf8_lossy(&c.alleles[0]),
        alts.join(","),
        info_s.join(";"),
        keys.join(":")
    );
    for g in &c.genotypes {
        let sep = if g.phase.is_some() { "|" } else { "/" };
        let mut v = vec![
            g.gt.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(sep),
            join(&g.ad),
            g.af.iter().map(|&x| format_vcf_double(x)).collect::<Vec<_>>().join(","),
            g.dp.to_string(),
            join(&g.f1r2),
            join(&g.f2r1),
            join(&g.fad),
        ];
        if let Some((pgt, pid, ps)) = &g.phase {
            v.push(pgt.to_string());
            v.push(pid.clone());
            v.push(ps.to_string());
        } else if phased {
            v.extend([".".to_string(), ".".to_string(), ".".to_string()]);
        }
        v.push(join(&g.sb));
        s.push('\t');
        s.push_str(&v.join(":"));
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn trims() {
        assert_eq!(trim_alleles(10, &[b"ACG".to_vec(), b"TCG".to_vec()]), (10, vec![b"A".to_vec(), b"T".to_vec()]));
        assert_eq!(trim_alleles(10, &[b"AT".to_vec(), b"A".to_vec()]), (10, vec![b"AT".to_vec(), b"A".to_vec()]));
        assert_eq!(trim_alleles(10, &[b"ATT".to_vec(), b"AT".to_vec(), b"GTT".to_vec()]), (10, vec![b"AT".to_vec(), b"A".to_vec(), b"GT".to_vec()]));
        assert_eq!(trim_alleles(10, &[b"ATTC".to_vec(), b"ATC".to_vec()]), (10, vec![b"AT".to_vec(), b"A".to_vec()]));
        let mut v = vec![30, 33];
        assert_eq!(median(&mut v), 32);
        let mut v = vec![5, 1, 9];
        assert_eq!(median(&mut v), 5);
    }
}
