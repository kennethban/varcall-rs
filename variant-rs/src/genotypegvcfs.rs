//! GenotypeGVCFs (`variant-rs genotypegvcfs`): GATK 4.7 `GenotypeGVCFs` on a GVCF file (single sample, or several
//! samples combined by CombineGVCFs), diploid. By-variant traversal: every input record on its own goes through
//! `ReferenceConfidenceVariantContextMerger.merge` (NON_REF dropped, PL / AD remapped, genotypes set to no-call or
//! the best match of the original, raw annotations combined), `GenotypeGVCFsEngine.regenotypeVC`
//! (`MinimalGenotypingEngine.calculateGenotypes`: the allele-frequency EM over all samples, QUAL, alleles to output,
//! PREFER_PLS genotypes, MLEAC / MLEAF / AS_QUAL), the annotations finalized and recomputed, and the genotypes
//! cleaned up. A GenomicsDB workspace (`gendb://`): `gdb_store` reads it, `gdb_merge` produces the records GATK's
//! GenomicsDB reader gives the walker.

use crate::reblock::{
    genotype_alleles, gls_of, make_genotype_call, max_index_f, n_likelihoods, parse_gt, parse_info, pl_index, subset_alleles, Assign, Gt,
    NON_REF, NO_CALL, SPAN_DEL,
};
use hc_genotype::mathutils::{log10_sum_log10, log10_sum_log10_2, normalize_log10_in_place};

pub(crate) fn refuse(why: String) -> ! {
    crate::gatk_cli::refuse(&why)
}

/// A (merged) multi-sample record.
#[derive(Clone, Debug)]
pub(crate) struct MRec {
    pub contig: String,
    pub start: i64,
    pub stop: i64,
    pub id: String,
    pub alleles: Vec<String>,
    /// QUAL: log10 P(error) when set
    pub log10_p_error: Option<f64>,
    pub filter: String,
    pub info: Vec<(String, String)>,
    pub gs: Vec<Gt>,
}

impl MRec {
    pub fn attr(&self, k: &str) -> Option<&str> {
        self.info.iter().find(|(x, _)| x == k).map(|(_, v)| v.as_str())
    }
    pub fn attr_int(&self, k: &str, default: i64) -> i64 {
        self.attr(k).and_then(|v| v.trim().parse::<i64>().ok()).unwrap_or(default)
    }
}

pub(crate) fn parse_mrec(line: &str, n_samples: usize) -> MRec {
    let f: Vec<&str> = line.split('\t').collect();
    assert_eq!(f.len(), 9 + n_samples, "record with {} columns, header has {n_samples} samples: {}", f.len(), &line[..line.len().min(200)]);
    let start: i64 = f[1].parse().unwrap();
    let up = |a: &str| if crate::reblock::is_symbolic(a) || a == SPAN_DEL { a.to_string() } else { a.to_ascii_uppercase() };
    let alleles: Vec<String> = std::iter::once(f[3].to_ascii_uppercase()).chain(f[4].split(',').filter(|a| *a != ".").map(up)).collect();
    let info = parse_info(f[7]);
    let end = info.iter().find(|(k, _)| k == "END").map(|(_, v)| v.parse::<i64>().unwrap()).unwrap_or(start + alleles[0].len() as i64 - 1);
    let gs: Vec<Gt> = f[9..].iter().map(|s| parse_gt(f[8], s)).collect();
    for g in &gs {
        // a haploid no-call: GenomicsDB's "." for a sample without a call
        if g.ploidy() == 1 && g.is_no_call() && g.pl.is_none() {
            continue;
        }
        if g.ploidy() != 2 {
            refuse(format!("{}:{}: ploidy {}: only diploid GVCFs are supported", f[0], start, g.ploidy()));
        }
    }
    let qual = if f[5] == "." { None } else { Some(f[5].to_string()) };
    MRec {
        contig: f[0].to_string(),
        start,
        stop: end,
        id: f[2].to_string(),
        alleles,
        log10_p_error: qual.as_ref().map(|q| q.parse::<f64>().unwrap() / -10.0),
        filter: f[6].to_string(),
        info,
        gs,
    }
}

// ------------------------------------------------------------------------------------------------
// ReferenceConfidenceVariantContextMerger.merge, one record at its own position

/// The merged genotypes and alleles; the attributes are combined by the annotation layer (`annotations`).
pub(crate) fn merge_genotypes(vc: &MRec) -> (Vec<String>, Vec<Gt>) {
    // collectTargetAlleles: the reference, then every called alt but NON_REF (the '*' / mixed rules only re-add
    // alleles that are already there for a record at its own position)
    let mut target: Vec<String> = vec![vc.alleles[0].clone()];
    for a in &vc.alleles[1..] {
        if a != NON_REF && !target.contains(a) {
            target.push(a.clone());
        }
    }
    let non_ref = vc.alleles.iter().position(|a| a == NON_REF).unwrap_or_else(|| {
        refuse(format!(
            "The list of input alleles must contain <NON_REF> as an allele but that is not the case at position {}; please use the Haplotype Caller with gVCF output to generate appropriate records",
            vc.start
        ))
    });
    // getIndexesOfRelevantAllelesForGVCF: target allele -> remapped index, NON_REF when absent
    let relevant: Vec<usize> = target.iter().enumerate().map(|(i, a)| if i == 0 { 0 } else { vc.alleles.iter().position(|x| x == a).unwrap_or(non_ref) }).collect();
    let new_to_old: Vec<usize> = (0..n_likelihoods(target.len())).map(|i| {
        let [j, k] = genotype_alleles(i);
        pl_index(relevant[j], relevant[k])
    }).collect();
    let mut gs = Vec::with_capacity(vc.gs.len());
    for g in &vc.gs {
        let mut gb = g.clone();
        if let Some(pl) = &g.pl {
            gb.pl = Some(new_to_old.iter().map(|&i| pl[i]).collect());
        }
        if let Some(ad) = &g.ad {
            // remapRLengthList, filler 0
            gb.ad = Some(relevant.iter().map(|&i| ad.get(i).copied().unwrap_or(0)).collect());
        }
        if exclude_from_annotations(g) {
            gb.al = vec![NO_CALL; g.ploidy()];
        }
        // callGTAlleles: GenotypeUtils.shouldBeCalled
        let method = if should_be_called(g) { Assign::BestMatch } else { Assign::SetToNoCall };
        let gls = g.pl.as_ref().map(|p| gls_of(p));
        make_genotype_call(&mut gb, method, gls.as_deref(), &target, g, &vc.alleles);
        gs.push(gb);
    }
    (target, gs)
}

/// `GenotypeGVCFsEngine.excludeFromAnnotations`
pub(crate) fn exclude_from_annotations(g: &Gt) -> bool {
    (g.is_hom_ref() || g.is_no_call()) && g.dp.is_none_or(|d| d == 0) && g.gq == Some(0)
}

/// `GenotypeUtils.shouldBeCalled`: informative PLs (not all zero) or GQ > 0
fn should_be_called(g: &Gt) -> bool {
    let non_informative = g.pl.as_ref().is_none_or(|p| p.iter().all(|&x| x == 0));
    !non_informative || g.gq.is_some_and(|q| q > 0)
}

/// `calculateVCDepth`: INFO DP, else the samples' MIN_DP / DP
pub(crate) fn vc_depth(vc: &MRec) -> i64 {
    if vc.attr("DP").is_some() {
        return vc.attr_int("DP", 0);
    }
    vc.gs
        .iter()
        .map(|g| match g.ext.iter().find(|(k, _)| k == "MIN_DP") {
            Some((_, v)) => v.parse::<i64>().unwrap(),
            None => g.dp.unwrap_or(0) as i64,
        })
        .sum()
}

// ------------------------------------------------------------------------------------------------
// GenotypingEngine.calculateGenotypes (MinimalGenotypingEngine, EMIT_VARIANTS_ONLY, PREFER_PLS)

pub(crate) struct Call {
    /// output alleles (indices into the merged alleles), the reference first
    pub alleles: Vec<usize>,
    pub log10_conf: f64,
    pub gs: Vec<Gt>,
    pub mleac: Vec<i32>,
    pub mleaf: Vec<f64>,
    pub as_qual: Vec<i32>,
}

pub(crate) struct Genotyper {
    ref_pseudo: f64,
    snp_pseudo: f64,
    indel_pseudo: f64,
    standard_conf: f64,
    max_alt: usize,
    as_calcs: bool,
    /// recordDeletions: (contig, start, end) of emitted deletions
    upstream: Vec<(String, i64, i64)>,
}

/// `GenotypeUtils.makeApproximateDiploidLog10LikelihoodsFromGQ`: PLs 0 / GQ / 10*GQ
fn approx_gls_from_gq(gq: i32, n_alleles: usize) -> Vec<f64> {
    // PLOIDY_2_HOM_VAR_SCALE_FACTOR = round(30 / -10.0 / log10(0.5))
    let scale = hc_genotype::javamath::java_round(30.0 / -10.0 / 0.5f64.log10()) as i32;
    let pls: Vec<i32> = (0..n_likelihoods(n_alleles))
        .map(|i| {
            let [a, b] = genotype_alleles(i);
            if i == 0 {
                0
            } else if a == 0 || b == 0 {
                gq
            } else {
                scale * gq
            }
        })
        .collect();
    gls_of(&pls)
}

/// The genotype's log10 likelihoods for the AF calculation (`genotypeIsUsableForAFCalculation` /
/// `log10NormalizedGenotypePosteriors`), None if it is skipped.
fn af_gls(g: &Gt, n_alleles: usize) -> Option<Vec<f64>> {
    if let Some(pl) = &g.pl {
        return Some(gls_of(pl));
    }
    if g.is_hom_ref() && g.gq.is_some() {
        return Some(approx_gls_from_gq(g.gq.unwrap(), n_alleles));
    }
    None
}

fn log10_posteriors(gls: &[f64], log10_af: &[f64]) -> Vec<f64> {
    let mut post: Vec<f64> = gls
        .iter()
        .enumerate()
        .map(|(g, &l)| {
            let [i, j] = genotype_alleles(g);
            // log10CombinationCount + likelihood + sum count * log10 AF
            let comb = if i == j { 0.0 } else { std::f64::consts::LOG10_2 };
            let s = if i == j { 0.0 + 2.0 * log10_af[i] } else { 0.0 + 1.0 * log10_af[i] + 1.0 * log10_af[j] };
            comb + l + s
        })
        .collect();
    normalize_log10_in_place(&mut post);
    post
}

pub(crate) struct AfResult {
    pub mle: Vec<i32>,
    pub log10_p_no_variant: f64,
    /// per allele (index 0 unused)
    pub log10_p_absent: Vec<f64>,
}

impl Genotyper {
    pub fn new(standard_conf: f64, max_alt: usize, as_calcs: bool) -> Genotyper {
        // heterozygosity 0.001, indel heterozygosity 1.25E-4, heterozygosity stdev 0.01
        let snp_het = 0.001f64;
        let ref_pseudo = snp_het / 0.01f64.powf(2.0);
        Genotyper { ref_pseudo, snp_pseudo: snp_het * ref_pseudo, indel_pseudo: 1.25E-4 * ref_pseudo, standard_conf, max_alt, as_calcs, upstream: Vec::new() }
    }

    /// `AlleleFrequencyCalculator.calculate` over the samples
    pub fn af(&self, alleles: &[String], gs: &[Gt]) -> AfResult {
        let n = alleles.len();
        let ref_len = alleles[0].len();
        let prior: Vec<f64> = (0..n)
            .map(|a| if a == 0 { self.ref_pseudo } else if alleles[a].len() == ref_len { self.snp_pseudo } else { self.indel_pseudo })
            .collect();
        let samples: Vec<Vec<f64>> = gs.iter().filter_map(|g| af_gls(g, n)).collect();
        let mut counts = vec![0.0f64; n];
        let flat = -crate::jmath::log10(n as f64);
        let mut log10_af = vec![flat; n];
        loop {
            let mut log10_res = vec![f64::NEG_INFINITY; n];
            for gls in &samples {
                let post = log10_posteriors(gls, &log10_af);
                for (g, &p) in post.iter().enumerate() {
                    let [i, j] = genotype_alleles(g);
                    if i == j {
                        log10_res[i] = log10_sum_log10_2(log10_res[i], p + crate::jmath::log10(2.0));
                    } else {
                        log10_res[i] = log10_sum_log10_2(log10_res[i], p + crate::jmath::log10(1.0));
                        log10_res[j] = log10_sum_log10_2(log10_res[j], p + crate::jmath::log10(1.0));
                    }
                }
            }
            let new_counts: Vec<f64> = log10_res.iter().map(|&x| crate::jmath::pow(10.0, x)).collect();
            let max_diff = counts.iter().zip(&new_counts).map(|(a, b)| (a - b).abs()).fold(f64::NEG_INFINITY, f64::max);
            counts = new_counts;
            let pp: Vec<f64> = prior.iter().zip(&counts).map(|(a, b)| a + b).collect();
            let s: f64 = pp.iter().sum();
            log10_af = pp.iter().map(|&x| crate::jmath::log10(x / s)).collect();
            if !(max_diff > 0.1) {
                break;
            }
        }
        let sd = alleles.iter().position(|a| a == SPAN_DEL);
        let mut p_absent = vec![0.0f64; n];
        let mut p_no_var = 0.0f64;
        for gls in &samples {
            let post = log10_posteriors(gls, &log10_af);
            match sd {
                None => p_no_var += post[0],
                Some(sd) => {
                    let v = [post[pl_index(0, 0)], post[pl_index(0, sd)], post[pl_index(sd, sd)]];
                    p_no_var += 0f64.min(log10_sum_log10(&v));
                }
            }
            if n == 2 && sd.is_none() {
                continue;
            }
            let mut bufs: Vec<Vec<f64>> = vec![Vec::new(); n];
            for (g, &p) in post.iter().enumerate() {
                let [i, j] = genotype_alleles(g);
                for (a, b) in bufs.iter_mut().enumerate() {
                    if a != i && a != j {
                        b.push(p);
                    }
                }
            }
            for a in 0..n {
                p_absent[a] += log10_sum_log10(&bufs[a]).min(0.0);
            }
        }
        if n == 2 && sd.is_none() {
            p_absent[1] = p_no_var;
        }
        AfResult { mle: counts[1..].iter().map(|&x| hc_genotype::javamath::java_round(x) as i32).collect(), log10_p_no_variant: p_no_var, log10_p_absent: p_absent }
    }

    fn covered_by_deletion(&self, contig: &str, start: i64) -> bool {
        self.upstream.iter().any(|(c, s, e)| c == contig && *s < start && start <= *e)
    }

    /// `calculateGenotypes(vc)`; None when the site is not emitted
    pub fn calculate(&mut self, vc: &MRec) -> Option<Call> {
        let n = vc.alleles.len();
        // cannotBeGenotyped
        let ok = n <= 50 && vc.gs.iter().filter(|g| !(g.is_no_call() || g.is_hom_ref())).all(|g| g.pl.is_some()) && vc.gs.iter().any(|g| g.pl.is_some());
        if !ok || vc.gs.is_empty() {
            return None;
        }
        // more alts than --max-alternate-alleles: the most likely ones (calculateMostLikelyAlleles, BEST_MATCH subset)
        let (reduced, reduced_gs): (Vec<usize>, Vec<Gt>) = if self.max_alt < n - 1 {
            let keep = most_likely_alleles(vc, self.max_alt);
            let g2 = vc.gs.iter().map(|g| subset_alleles(g, &vc.alleles, &keep, Assign::BestMatch, &[])).collect();
            (keep, g2)
        } else {
            ((0..n).collect(), vc.gs.clone())
        };
        let red_alleles: Vec<String> = reduced.iter().map(|&i| vc.alleles[i].clone()).collect();
        let af = self.af(&red_alleles, &reduced_gs);
        // calculateOutputAlleleSubset
        let threshold = self.standard_conf * -0.1;
        let mut outs: Vec<usize> = Vec::new(); // positions in `reduced`
        let mut mono = true;
        for k in 1..reduced.len() {
            let plausible = af.log10_p_absent[k] + 1.0e-10 < threshold;
            let spurious = red_alleles[k] == SPAN_DEL && !self.covered_by_deletion(&vc.contig, vc.start);
            if plausible && !spurious {
                outs.push(k);
            }
            mono &= !(plausible && !spurious);
        }
        let log10_conf = af.log10_p_no_variant + 0.0;
        let phred = (-10.0 * log10_conf) + 0.0;
        // passesEmitThreshold (EMIT_VARIANTS_ONLY): not monomorphic and QUAL >= stand-call-conf
        if !(!mono && phred >= self.standard_conf) {
            return None;
        }
        if outs.len() == 1 && red_alleles[outs[0]] == SPAN_DEL {
            return None;
        }
        let mut out_alleles = vec![0usize];
        out_alleles.extend(outs.iter().map(|&k| reduced[k]));
        // recordDeletions
        self.upstream.retain(|(c, _, e)| *c == vc.contig && *e >= vc.start);
        let ref_len = vc.alleles[0].len() as i64;
        for &a in &out_alleles {
            let del = ref_len - vc.alleles[a].len() as i64;
            if del > 0 {
                self.upstream.push((vc.contig.clone(), vc.start, vc.start + del));
            }
        }
        let gs: Vec<Gt> = vc.gs.iter().map(|g| subset_alleles(g, &vc.alleles, &out_alleles, Assign::PreferPls, &[])).collect();
        let mleac: Vec<i32> = outs.iter().map(|&k| af.mle[k - 1]).collect();
        let an = gs.iter().flat_map(|g| g.al.iter()).filter(|&&a| a != NO_CALL).count() as f64;
        let mleaf = mleac.iter().map(|&c| 1.0f64.min(c as f64 / an)).collect();
        let mut as_qual = Vec::new();
        if self.as_calcs {
            if reduced.len() > 2 {
                for &k in &outs {
                    as_qual.push(hc_genotype::javamath::java_round(af.log10_p_absent[k] * -10.0) as i32);
                }
            } else {
                as_qual.push(hc_genotype::javamath::java_round(af.log10_p_no_variant * -10.0) as i32);
            }
        }
        Some(Call { alleles: out_alleles, log10_conf, gs, mleac, mleaf, as_qual })
    }
}

/// `AlleleSubsettingUtils.calculateMostLikelyAlleles(vc, 2, numAltAllelesToKeep, ensureReturnContainsAlt = false)`
fn most_likely_alleles(vc: &MRec, keep: usize) -> Vec<usize> {
    let n = vc.alleles.len();
    let non_ref = vc.alleles.iter().position(|a| a == NON_REF);
    let n_proper = n - if non_ref.is_some() { 2 } else { 1 };
    if keep >= n_proper {
        return (0..n).collect();
    }
    let mut sums = vec![0f64; n];
    for g in &vc.gs {
        let Some(pl) = &g.pl else { continue };
        let gls = gls_of(pl);
        let ix = max_index_f(&gls, 0);
        let diff = (gls[ix] - gls[0]).abs();
        let ga = genotype_alleles(ix);
        for a in 1..n {
            if ga.contains(&a) {
                sums[a] += diff;
            }
        }
    }
    let mut cand: Vec<usize> = (1..n).filter(|&a| Some(a) != non_ref).collect();
    cand.sort_by(|&a, &b| sums[b].partial_cmp(&sums[a]).unwrap());
    cand.truncate(keep);
    (0..n).filter(|&i| i == 0 || Some(i) == non_ref || cand.contains(&i)).collect()
}

// ------------------------------------------------------------------------------------------------
// driver

const USAGE: &str = "usage: variant-rs genotypegvcfs -R ref.fa -V in.g.vcf.gz -O out.vcf.gz [-D dbsnp.vcf.gz] [-G StandardAnnotation \
-G AS_StandardAnnotation] [-L interval ...] [--only-output-calls-starting-in-intervals] [--merge-input-intervals]";

pub fn run(argv: &[String]) {
    let t0 = std::time::Instant::now();
    let (mut input, mut output, mut reference, mut dbsnp) = (None::<String>, None::<String>, None::<String>, None::<String>);
    let (mut groups, mut ivs, mut starts_in, mut merge_ivs) = (Vec::<String>::new(), Vec::<String>::new(), false, false);
    let mut i = 0;
    while i < argv.len() {
        let v = || argv.get(i + 1).cloned().unwrap_or_else(|| panic!("{} needs a value\n{USAGE}", argv[i]));
        match argv[i].as_str() {
            "--only-output-calls-starting-in-intervals" => {
                starts_in = true;
                i += 1;
                continue;
            }
            "--merge-input-intervals" => {
                merge_ivs = true;
                i += 1;
                continue;
            }
            "-V" | "--variant" => input = Some(v()),
            "-O" | "--output" => output = Some(v()),
            "-R" | "--reference" => reference = Some(v()),
            "-D" | "--dbsnp" => dbsnp = Some(v()),
            "-G" | "--annotation-group" => groups.push(v()),
            "-L" | "--intervals" => ivs.push(v()),
            x => panic!("unknown argument {x}\n{USAGE}"),
        }
        i += 2;
    }
    let (input, output, reference) = (input.expect(USAGE), output.expect(USAGE), reference.expect(USAGE));
    let gendb: Option<crate::gdb_store::Workspace> = input.strip_prefix("gendb://").map(|p| crate::gdb_store::Workspace::open(std::path::Path::new(p)));
    let as_annotations = groups.iter().any(|g| g == "AS_StandardAnnotation");
    let dict = crate::intervals::Dict::from_reference(&reference);
    // traversal intervals (-L, merged; --merge-input-intervals: one spanning interval per contig) and the output
    // filter (STARTS_IN: the user's intervals)
    let mut user: Vec<crate::intervals::Interval> = Vec::new();
    for l in &ivs {
        user.extend(crate::intervals::parse_l_arg(l, &dict).unwrap_or_else(|e| panic!("{e}")));
    }
    let user = crate::intervals::sort_merge(user);
    let traversal: Vec<crate::intervals::Interval> = if merge_ivs {
        let mut t: Vec<crate::intervals::Interval> = Vec::new();
        for iv in &user {
            match t.last_mut() {
                Some(l) if l.tid == iv.tid => l.end = l.end.max(iv.end),
                _ => t.push(crate::intervals::Interval { tid: iv.tid, start: iv.start, end: iv.end }),
            }
        }
        t
    } else {
        user.iter().map(|x| crate::intervals::Interval { tid: x.tid, start: x.start, end: x.end }).collect()
    };
    let overlaps = |set: &[crate::intervals::Interval], c: &str, s: i64, e: i64| -> bool {
        let Some(&t) = dict.index.get(c) else { return false };
        set.iter().any(|iv| iv.tid == t && iv.start <= e && s <= iv.end)
    };
    let mut lines: Box<dyn Iterator<Item = String>> = Box::new(std::iter::empty());
    let mut hdr: Vec<String> = Vec::new();
    let chrom = match &gendb {
        // the workspace's merged header plus the PASS filter line GenomicsDB's reader adds; samples in row order
        Some(ws) => {
            hdr.extend(ws.header.lines().filter(|l| l.starts_with("##")).map(|l| l.to_string()));
            hdr.push("##FILTER=<ID=PASS,Description=\"All filters passed\">".into());
            format!("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}", ws.samples.join("\t"))
        }
        None => {
            lines = crate::mergevcfs::open_lines(&input);
            loop {
                let l = lines.next().unwrap_or_else(|| panic!("{input}: no #CHROM line"));
                if l.starts_with("#CHROM") {
                    break l;
                }
                hdr.push(l);
            }
        }
    };
    let n_samples = chrom.split('\t').count() - 9;
    // output samples sorted (TreeSet of names)
    let names: Vec<&str> = chrom.split('\t').skip(9).collect();
    let mut order: Vec<usize> = (0..n_samples).collect();
    order.sort_by(|&x, &y| names[x].cmp(names[y]));
    let types: std::collections::HashMap<String, String> = hdr
        .iter()
        .filter_map(|l| l.strip_prefix("##INFO=<ID="))
        .filter_map(|r| {
            let id = r.split(',').next()?;
            let t = r.split(",Type=").nth(1)?.split([',', '>']).next()?;
            Some((id.to_string(), t.to_string()))
        })
        .collect();
    // header (setupVCFWriter): the input's lines without GVCFBlock, plus the annotation / genotyping lines
    let mut set: Vec<String> = hdr.iter().filter(|l| !l.starts_with("##fileformat=") && !l.starts_with("##GVCFBlock")).cloned().collect();
    for l in include_str!("../tables/genotypegvcfs_header_lines.txt").lines() {
        let (cond, line) = l.split_once('\t').unwrap();
        if cond == "all" || (cond == "as" && as_annotations) || (cond == "dbsnp" && dbsnp.is_some()) {
            set.push(line.to_string());
        }
    }
    set.push(format!(
        "##variant-rsCommandLine=<ID=GenotypeGVCFs,CommandLine=\"genotypegvcfs {}\",Version=\"{}\",Date=\"{}\">",
        argv.join(" ").replace('"', "\\\""),
        env!("CARGO_PKG_VERSION"),
        crate::now_string()
    ));
    let mut seen = std::collections::HashSet::new();
    set.retain(|l| seen.insert(l.clone()));
    let key = |l: &String| -> String {
        if l.starts_with("##contig=<") {
            "##contig=<".into()
        } else if let Some(rest) = l.strip_prefix("##variant-rsCommandLine") {
            format!("##GATKCommandLine{rest}")
        } else {
            l.clone()
        }
    };
    set.sort_by_key(|l| key(l));
    let mut header = String::from("##fileformat=VCFv4.2\n");
    for l in &set {
        header.push_str(l);
        header.push('\n');
    }
    let mut sorted_names: Vec<&str> = names.clone();
    sorted_names.sort();
    header.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for n in &sorted_names {
        header.push('\t');
        header.push_str(n);
    }
    header.push('\n');
    let contigs: std::collections::HashSet<String> = dict.names.iter().cloned().collect();
    let mut ann = Annotator {
        as_ann: as_annotations,
        types,
        rng: hc_genotype::javamath::JavaRandom::new(47382911),
        dbsnp: dbsnp.as_ref().map(|p| crate::dbsnp::Dbsnp::load(p, &contigs)),
    };
    let gz = crate::output::block_compressed(&output);
    let mut out = crate::output::OutFile::create(&output, &header, if gz { Some(2) } else { None }).unwrap_or_else(|e| panic!("{e}"));
    let mut gt = Genotyper::new(30.0, 6, as_annotations);
    let mut buf = String::with_capacity(8 << 20);
    let (mut n_in, mut n_out) = (0u64, 0u64);
    let is_gendb = gendb.is_some();
    let mut dump = std::env::var_os("VARIANT_RS_GDB_DUMP").filter(|p| !p.is_empty()).map(|p| std::io::BufWriter::new(std::fs::File::create(p).unwrap()));
    let mut handle = |l: &str| {
        if l.is_empty() {
            return;
        }
        if let Some(d) = dump.as_mut() {
            use std::io::Write;
            writeln!(d, "{l}").unwrap();
        }
        n_in += 1;
        let rec = parse_mrec(l, n_samples);
        if !is_gendb && !user.is_empty() && !overlaps(&traversal, &rec.contig, rec.start, rec.stop) {
            return;
        }
        if !rec.alleles[1..].iter().any(|a| a != NON_REF) {
            return;
        }
        let Some((o, attrs)) = genotype_record(&rec, &mut gt, &mut ann) else { return };
        if starts_in && !user.is_empty() && !overlaps(&user, &o.contig, o.start, o.start) {
            return;
        }
        encode(&o, &attrs, &order, &mut buf);
        n_out += 1;
        if buf.len() >= 4 << 20 {
            out.add_text(buf.as_bytes()).unwrap();
            buf.clear();
        }
    };
    match &gendb {
        None => {
            for l in lines {
                handle(&l);
            }
        }
        Some(ws) if ivs.is_empty() => {
            // no -L (as nf-core/sarek runs it on each workspace): GenomicsDBFeatureReader.iterator() queries every
            // array over its whole range [start, end], arrays in ChrArrayFolderComparator order (contig name as
            // String.compareTo, then start), and the walker takes the records in that order
            let fasta = crate::reference::Fasta::open(&reference);
            let mut arrs: Vec<(String, i64, i64, &std::path::PathBuf)> = Vec::new();
            for p in &ws.arrays {
                let name = p.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
                let f: Vec<&str> = name.split('$').collect();
                match (f.len(), f.get(1).and_then(|x| x.parse::<i64>().ok()), f.get(2).and_then(|x| x.parse::<i64>().ok())) {
                    (3, Some(lo), Some(hi)) => arrs.push((f[0].to_string(), lo, hi, p)),
                    _ => refuse(format!("-V gendb:// without -L: array {name} is not named contig$start$end; not supported")),
                }
            }
            arrs.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
            let mut out_lines = Vec::new();
            for (contig, lo, hi, p) in &arrs {
                let a = crate::gdb_store::open_array(p);
                let mut sc = crate::gdb_merge::Scanner::new(ws, &a, contig);
                out_lines.clear();
                sc.scan(*lo, *hi, &fasta, &mut out_lines);
                for l in &out_lines {
                    handle(l);
                }
            }
        }
        Some(ws) => {
            let fasta = crate::reference::Fasta::open(&reference);
            for iv in &traversal {
                let contig = &dict.names[iv.tid];
                // the arrays (imported intervals contig$start$end) overlapping the traversal interval, in order
                let mut arrs: Vec<(i64, i64, crate::gdb_store::Array)> = ws
                    .arrays
                    .iter()
                    .filter_map(|p| {
                        let name = p.file_name()?.to_string_lossy().into_owned();
                        let mut f = name.rsplitn(3, '$');
                        let (hi, lo, c) = (f.next()?.parse::<i64>().ok()?, f.next()?.parse::<i64>().ok()?, f.next()?);
                        (c == contig && lo <= iv.end && iv.start <= hi).then(|| (lo, hi, crate::gdb_store::open_array(p)))
                    })
                    .collect();
                arrs.sort_by_key(|x| x.0);
                let mut scanners: Vec<(i64, i64, crate::gdb_merge::Scanner)> =
                    arrs.iter().map(|(lo, hi, a)| (*lo, *hi, crate::gdb_merge::Scanner::new(ws, a, contig))).collect();
                crate::gdb_merge::interval_records(ws, &mut scanners, contig, iv.start, iv.end, &fasta, |l| handle(l));
            }
        }
    }
    drop(handle);
    if let Some(mut d) = dump {
        use std::io::Write;
        d.flush().unwrap();
    }
    out.add_text(buf.as_bytes()).unwrap();
    out.finish().unwrap();
    eprintln!("genotypegvcfs: {n_in} records in, {n_out} out ({:.1}s)", t0.elapsed().as_secs_f64());
}

// ------------------------------------------------------------------------------------------------
// annotations (docs/specs/genotypegvcfs_annotations.md: VariantAnnotatorEngine with StandardAnnotation [+ AS_StandardAnnotation])

/// An attribute value as htsjdk holds it (the writer formats Double with formatVCFDouble)
#[derive(Clone, Debug)]
enum Val {
    S(String),
    I(i64),
    D(f64),
    IL(Vec<i64>),
    DL(Vec<f64>),
    Flag,
}

impl Val {
    fn text(&self) -> String {
        use hc_genotype::javamath::format_vcf_double as vd;
        match self {
            Val::S(s) => s.clone(),
            Val::I(x) => x.to_string(),
            Val::D(x) => vd(*x),
            Val::IL(v) => v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","),
            Val::DL(v) => v.iter().map(|&x| vd(x)).collect::<Vec<_>>().join(","),
            Val::Flag => String::new(),
        }
    }
}

type Attrs = Vec<(String, Val)>;

fn put(a: &mut Attrs, k: &str, v: Val) {
    match a.iter_mut().find(|(x, _)| x == k) {
        Some(e) => e.1 = v,
        None => a.push((k.to_string(), v)),
    }
}

fn get<'a>(a: &'a Attrs, k: &str) -> Option<&'a Val> {
    a.iter().find(|(x, _)| x == k).map(|(_, v)| v)
}

fn remove(a: &mut Attrs, k: &str) {
    a.retain(|(x, _)| x != k);
}

fn f(x: f64, n: usize) -> String {
    hc_genotype::javamath::format_fixed(x, n)
}

const RAW_AS: [&str; 8] = ["AS_RAW_BaseQRankSum", "AS_SB_TABLE", "AS_RAW_MQRankSum", "AS_QUALapprox", "AS_QUAL", "AS_VarDP", "AS_RAW_MQ", "AS_RAW_ReadPosRankSum"];

pub(crate) struct Annotator {
    pub as_ann: bool,
    /// INFO header types of the input ("Integer", "Float", ...)
    pub types: std::collections::HashMap<String, String>,
    pub rng: hc_genotype::javamath::JavaRandom,
    pub dbsnp: Option<crate::dbsnp::Dbsnp>,
}

/// Histogram (bin 0.1, BIN_EPSILON 0.01): bin -> count
type Hist = std::collections::BTreeMap<i64, i64>;

fn bin_of(d: f64) -> i64 {
    hc_genotype::javamath::java_round(((d + 0.01 * 0.1) / 0.1).floor())
}

/// `AS_RankSumTest.parseRawDataString`: per allele of `alleles` a histogram
fn parse_rank_raw(raw: &str, alleles: &[String]) -> Vec<Hist> {
    let mut out = vec![Hist::new(); alleles.len()];
    let raw = if raw.starts_with('[') { &raw[1..raw.len() - 1] } else { raw };
    let mut toks: Vec<&str> = raw.split('|').collect();
    while toks.last().is_some_and(|t| t.is_empty()) {
        toks.pop();
    }
    for (i, t) in toks.iter().enumerate() {
        if i >= alleles.len() {
            break;
        }
        let p: Vec<&str> = t.split(',').collect();
        let mut j = 0;
        while j < p.len() {
            if !p[j].is_empty() {
                let v: f64 = p[j].trim().parse().unwrap();
                if !v.is_nan() && j + 1 < p.len() && !p[j + 1].is_empty() {
                    let c: i64 = p[j + 1].trim().parse().unwrap();
                    *out[i].entry(bin_of(v)).or_insert(0) += c;
                }
            }
            j += 2;
        }
    }
    out
}

fn hist_text(h: &Hist) -> String {
    h.iter().map(|(&k, &c)| format!("{},{}", f(k as f64 * 0.1, 1), c)).collect::<Vec<_>>().join(",")
}

/// `Histogram.median`
fn hist_median(h: &Hist) -> Option<f64> {
    let n: i64 = h.values().sum();
    if n == 0 {
        return None;
    }
    let mid = (n + 1) / 2;
    let mut counter = 0;
    let mut first: Option<f64> = None;
    for (&k, &c) in h {
        counter += c;
        if counter > mid {
            return Some(match first {
                None => k as f64 * 0.1,
                Some(fk) => (fk + k as f64) / 2.0 * 0.1,
            });
        }
        if counter == mid {
            if n % 2 == 1 {
                return Some(k as f64 * 0.1);
            }
            first = Some(k as f64);
        }
    }
    None
}

/// `AS_StrandBiasTest.parseRawDataString`: per allele its [fwd, rev] list (None when the field count is wrong)
fn parse_sb_raw(raw: &str, n: usize) -> Vec<Vec<i64>> {
    let raw: String = if raw.starts_with('[') { raw[1..raw.len() - 1].chars().filter(|c| !c.is_whitespace()).collect() } else { raw.to_string() };
    let fields: Vec<&str> = raw.split('|').collect();
    assert_eq!(fields.len(), n, "AS_SB_TABLE {raw}: {} fields for {n} alleles", fields.len());
    fields.iter().map(|fl| fl.split(',').filter(|t| !t.is_empty()).map(|t| t.trim().parse().unwrap()).collect()).collect()
}

/// `split("\\|")` (trailing empty strings removed), per allele: None for an empty token
fn parse_pipe_doubles(raw: &str, n: usize) -> Vec<Option<f64>> {
    let mut toks: Vec<&str> = raw.split('|').collect();
    while toks.last().is_some_and(|t| t.is_empty()) {
        toks.pop();
    }
    (0..n).map(|i| toks.get(i).filter(|t| !t.is_empty()).map(|t| t.trim().parse().unwrap())).collect()
}

fn sb_of(g: &Gt) -> Option<[i64; 4]> {
    let (_, v) = g.ext.iter().find(|(k, _)| k == "SB")?;
    let x: Vec<i64> = v.split(',').map(|t| t.trim().parse().unwrap()).collect();
    Some([x[0], x[1], x[2], x[3]])
}

/// `StrandBiasTest.getTableFromSamples(genotypes, minCount)`
fn table_from_samples(gs: &[Gt], min_count: i64) -> Option<[[i32; 2]; 2]> {
    let mut sb = [0i64; 4];
    let mut found = false;
    for g in gs {
        if let Some(d) = sb_of(g) {
            found = true;
            if d[0] + d[1] + d[2] + d[3] > min_count {
                for k in 0..4 {
                    sb[k] += d[k];
                }
            }
        }
    }
    found.then(|| [[sb[0] as i32, sb[1] as i32], [sb[2] as i32, sb[3] as i32]])
}

/// `AS_QualByDepth.getAlleleDepths`
fn allele_depths(gs: &[Gt]) -> Option<Vec<i64>> {
    let n = gs.iter().find_map(|g| g.ad.as_ref().map(|a| a.len()))?;
    let mut d = vec![0i64; n];
    for g in gs {
        if !(g.is_het() || g.is_hom_var()) {
            continue;
        }
        if let Some(ad) = &g.ad {
            let tot: i64 = ad.iter().map(|&x| x as i64).sum();
            if tot - ad[0] as i64 > 1 {
                for (x, &y) in d.iter_mut().zip(ad) {
                    *x += y as i64;
                }
            }
        }
    }
    Some(d)
}

/// `QualByDepth.getDepth` (likelihoods null)
fn qd_depth(gs: &[Gt]) -> i64 {
    let (mut depth, mut restricted) = (0i64, 0i64);
    for g in gs {
        if !(g.is_het() || g.is_hom_var()) {
            continue;
        }
        if let Some(ad) = &g.ad {
            let tot: i64 = ad.iter().map(|&x| x as i64).sum();
            if tot != 0 {
                if tot - ad[0] as i64 > 1 {
                    restricted += tot;
                }
                depth += tot;
                continue;
            }
        }
        if let Some(d) = g.dp {
            depth += d as i64;
        }
    }
    if restricted > 0 { restricted } else { depth }
}

impl Annotator {
    /// `QualByDepth.fixTooHighQD`
    fn fix_qd(&mut self, qd: f64) -> f64 {
        if qd < 35.0 { qd } else { 30.0 + self.rng.next_gaussian() * 3.0 }
    }

    /// merger attributes: raw keys (combined per the spec) and typed non-raw values
    fn merge_attrs(&self, rec: &MRec, target: &[String], depth: i64) -> Attrs {
        let reducible = |k: &str| k == "RAW_MQandDP" || (self.as_ann && RAW_AS.contains(&k));
        let mut raw: std::collections::HashMap<&str, &str> = std::collections::HashMap::new();
        let mut nonraw: Vec<(String, Vec<Val>)> = Vec::new();
        for (k, v) in &rec.info {
            if ["TLOD", "SQ", "POPAF", "MMQ", "MBQ", "MPOS", "MFRL"].contains(&k.as_str()) {
                continue;
            }
            if reducible(k) {
                raw.insert(k, v);
                continue;
            }
            let mut vals = Vec::new();
            // a comma-separated value is a List whose toString is "[a, b]": tokens after the first keep a space
            let toks: Vec<String> = if !v.contains(',') {
                vec![v.clone()]
            } else {
                let mut t: Vec<String> = Vec::new();
                for (i, p) in v.split(',').enumerate() {
                    t.push(if i == 0 { p.to_string() } else { format!(" {p}") });
                }
                t
            };
            let typ = self.types.get(k).map(|s| s.as_str());
            let flag = v.is_empty();
            for t in toks {
                let parsed = match typ {
                    _ if flag => None,
                    Some("Integer") => t.parse::<i64>().ok().map(Val::I),
                    Some("Float") => t.trim().parse::<f64>().ok().or_else(|| if t.trim() == "NaN" { Some(f64::NAN) } else { None }).map(Val::D),
                    None => if t.contains('.') { t.trim().parse::<f64>().ok().map(Val::D) } else { t.parse::<i64>().ok().map(Val::I) },
                    _ => None,
                };
                match parsed {
                    Some(x) => vals.push(x),
                    None => break,
                }
            }
            match nonraw.iter_mut().find(|(x, _)| x == k) {
                Some(e) => e.1 = vals,
                None => nonraw.push((k.clone(), vals)),
            }
        }
        let mut a: Attrs = Vec::new();
        let remapped = &rec.alleles;
        let idx_in_remapped = |al: &str| remapped.iter().position(|x| x == al);
        if self.as_ann {
            for key in ["AS_RAW_BaseQRankSum", "AS_RAW_MQRankSum", "AS_RAW_ReadPosRankSum"] {
                if let Some(r) = raw.get(key) {
                    let hs = parse_rank_raw(r, remapped);
                    let mut s = String::new();
                    for (i, al) in target.iter().enumerate() {
                        if i != 0 {
                            s.push('|');
                        }
                        if let Some(h) = idx_in_remapped(al).map(|j| &hs[j]) {
                            if !h.is_empty() {
                                s.push_str(&hist_text(h));
                            }
                        }
                    }
                    put(&mut a, key, Val::S(s));
                }
            }
            if let Some(r) = raw.get("AS_SB_TABLE") {
                let sb = parse_sb_raw(r, remapped.len());
                let s = target
                    .iter()
                    .map(|al| match idx_in_remapped(al) {
                        Some(j) => sb[j].iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","),
                        None => "0,0".to_string(),
                    })
                    .collect::<Vec<_>>()
                    .join("|");
                put(&mut a, "AS_SB_TABLE", Val::S(s));
            }
            if let Some(r) = raw.get("AS_RAW_MQ") {
                let v = parse_pipe_doubles(r, remapped.len());
                let s = target.iter().map(|al| f(idx_in_remapped(al).and_then(|j| v[j]).unwrap_or(0.0), 2)).collect::<Vec<_>>().join("|");
                put(&mut a, "AS_RAW_MQ", Val::S(s));
            }
        }
        if let Some(r) = raw.get("RAW_MQandDP") {
            let t: Vec<i64> = r.trim().trim_matches(|c| c == '[' || c == ']').split(',').map(|x| x.trim().parse().unwrap()).collect();
            assert_eq!(t.len(), 2, "RAW_MQandDP {r}");
            put(&mut a, "RAW_MQandDP", Val::S(format!("{},{}", t[0], t[1])));
        }
        for (k, vals) in nonraw {
            if vals.is_empty() {
                continue;
            }
            let v = if vals.len() == 1 {
                vals[0].clone()
            } else {
                // Utils.getMedianValue: upper median
                let mut v = vals;
                v.sort_by(|x, y| match (x, y) {
                    (Val::I(a), Val::I(b)) => a.cmp(b),
                    (Val::D(a), Val::D(b)) => crate::jmath::java_cmp(*a, *b),
                    _ => std::cmp::Ordering::Equal,
                });
                v[v.len() / 2].clone()
            };
            put(&mut a, &k, v);
        }
        if depth > 0 {
            put(&mut a, "DP", Val::S(depth.to_string()));
        }
        for k in ["AC", "AF", "AN", "MLEAC", "MLEAF", "END", "ECNT", "ECNTH"] {
            remove(&mut a, k);
        }
        a
    }

    /// `finalizeAnnotations`: `vc_alleles` / `gs` = the regenotyped VC, `target` = the merged alleles
    fn finalize(&mut self, a: &mut Attrs, target: &[String], vc_alleles: &[String], gs: &[Gt], as_qual: &[i32]) {
        let alt_of = |al: &str| target.iter().position(|x| x == al).unwrap();
        let text = |a: &Attrs, k: &str| get(a, k).map(|v| v.text());
        if self.as_ann {
            // AS_BaseQRankSum, AS_FS, AS_MQRankSum, AS_QD, AS_MQ, AS_ReadPosRankSum, AS_SOR (class-name order)
            let rank = |a: &mut Attrs, raw: &str, key: &str| {
                if let Some(r) = text(a, raw) {
                    let hs = parse_rank_raw(&r, target);
                    let s = vc_alleles[1..].iter().map(|al| hist_median(&hs[alt_of(al)]).map_or(".".to_string(), |m| f(m, 3))).collect::<Vec<_>>().join(",");
                    put(a, key, Val::S(s));
                }
            };
            rank(a, "AS_RAW_BaseQRankSum", "AS_BaseQRankSum");
            let sb = text(a, "AS_SB_TABLE").map(|r| parse_sb_raw(&r, target.len()));
            if let Some(sb) = &sb {
                let s = vc_alleles[1..]
                    .iter()
                    .map(|al| {
                        let (r, x) = (&sb[0], &sb[alt_of(al)]);
                        if x.is_empty() {
                            ".".to_string()
                        } else {
                            hc_genotype::annotate::fisher_strand([[r[0] as i32, r[1] as i32], [x[0] as i32, x[1] as i32]])
                        }
                    })
                    .collect::<Vec<_>>()
                    .join(",");
                put(a, "AS_FS", Val::S(s));
            }
            rank(a, "AS_RAW_MQRankSum", "AS_MQRankSum");
            if !as_qual.is_empty() {
                if let Some(depth) = allele_depths(gs) {
                    assert_eq!(as_qual.len(), vc_alleles.len() - 1, "AS_QUAL entries");
                    let mut v = Vec::new();
                    for (i, &q) in as_qual.iter().enumerate() {
                        let qd = q as f64 / (depth[i + 1] as f64 + depth[0] as f64);
                        v.push(f(self.fix_qd(qd), 2));
                    }
                    put(a, "AS_QD", Val::S(v.join(",")));
                }
            }
            if let Some(r) = text(a, "AS_RAW_MQ") {
                let per = parse_pipe_doubles(&r, target.len());
                let s = vc_alleles[1..]
                    .iter()
                    .enumerate()
                    .map(|(k, al)| {
                        let ad_sum: i64 = gs.iter().filter_map(|g| g.ad.as_ref().map(|ad| ad[k + 1] as i64)).sum();
                        match per[alt_of(al)] {
                            None => ".".to_string(),
                            Some(v) => f((v / ad_sum as f64).sqrt(), 2),
                        }
                    })
                    .collect::<Vec<_>>()
                    .join(",");
                put(a, "AS_MQ", Val::S(s));
            }
            rank(a, "AS_RAW_ReadPosRankSum", "AS_ReadPosRankSum");
            if let Some(sb) = &sb {
                let s = vc_alleles[1..]
                    .iter()
                    .map(|al| {
                        let (r, x) = (&sb[0], &sb[alt_of(al)]);
                        if x.is_empty() {
                            ".".to_string()
                        } else {
                            hc_genotype::annotate::strand_odds_ratio([[r[0] as i32, r[1] as i32], [x[0] as i32, x[1] as i32]])
                        }
                    })
                    .collect::<Vec<_>>()
                    .join(",");
                put(a, "AS_SOR", Val::S(s));
            }
            for k in RAW_AS {
                remove(a, k);
            }
        }
        if let Some(r) = text(a, "RAW_MQandDP") {
            let t: Vec<i64> = r.split(',').map(|x| x.trim().parse().unwrap()).collect();
            put(a, "MQ", Val::S(f((t[0] as f64 / t[1] as f64).sqrt(), 2)));
        }
        remove(a, "RAW_MQandDP");
        for k in ["AS_QUAL", "QUALapprox", "VarDP", "RAW_GT_COUNT"] {
            remove(a, k);
        }
    }

    /// `annotateContext` (likelihoods null): AN / AC / AF, ExcessHet, FS, QD, SOR (InbreedingCoeff and
    /// AS_InbreedingCoeff need 10 samples), then dbSNP
    fn annotate(&mut self, a: &mut Attrs, o: &mut MRec, log10_conf: f64) {
        let gs = &o.gs;
        let n_alt = o.alleles.len() - 1;
        // AS_InbreedingCoeff (HeterozygosityCalculator; 10 samples or more)
        if self.as_ann {
            if let Some(v) = as_inbreeding_coeff(gs, o.alleles.len()) {
                put(a, "AS_InbreedingCoeff", Val::S(v));
            }
        }
        // ChromosomeCounts
        let an: i64 = gs.iter().flat_map(|g| g.al.iter()).filter(|&&x| x != NO_CALL).count() as i64;
        if an > 0 {
            put(a, "AN", Val::I(an));
            let ac: Vec<i64> = (1..=n_alt).map(|k| gs.iter().flat_map(|g| g.al.iter()).filter(|&&x| x == k).count() as i64).collect();
            let af: Vec<f64> = ac.iter().map(|&c| c as f64 / an as f64).collect();
            if n_alt == 1 {
                put(a, "AC", Val::I(ac[0]));
                put(a, "AF", Val::D(af[0]));
            } else {
                put(a, "AC", Val::IL(ac));
                put(a, "AF", Val::DL(af));
            }
        }
        // ExcessHet
        let counts = genotype_counts(gs, o.alleles.len(), true);
        if let Some(v) = excess_het(gs, counts) {
            put(a, "ExcessHet", Val::S(v));
        }
        // FisherStrand
        if let Some(t) = table_from_samples(gs, 2) {
            put(a, "FS", Val::S(hc_genotype::annotate::fisher_strand(t)));
        }
        // InbreedingCoeff (10 samples or more)
        if gs.len() >= 10 {
            let n = gs.iter().filter(|g| (g.is_called() && (g.pl.is_some() || g.gq.is_some())) || g.pl.is_some()).count();
            if n >= 10 {
                let [refs, hets, homs] = genotype_counts(gs, o.alleles.len(), false);
                let p = (2.0 * refs + hets) / (2.0 * (refs + hets + homs));
                let q = 1.0 - p;
                let ic = 1.0 - (hets / (2.0 * p * q * n as f64));
                put(a, "InbreedingCoeff", Val::S(f(ic, 4)));
            }
        }
        // QualByDepth
        let depth = qd_depth(gs);
        if depth != 0 {
            let qual = -10.0 * log10_conf;
            let qd = self.fix_qd(qual / depth as f64);
            put(a, "QD", Val::S(f(qd, 2)));
        }
        // StrandOddsRatio
        if let Some(t) = table_from_samples(gs, 0) {
            put(a, "SOR", Val::S(hc_genotype::annotate::strand_odds_ratio(t)));
        }
        // dbSNP: ID and DB
        if let Some(db) = &self.dbsnp {
            let ids = db.ids(&o.contig, o.start, &o.alleles);
            if !ids.is_empty() {
                let rs = ids.join(";");
                if o.id == "." {
                    o.id = rs;
                } else if !o.id.contains(&rs) {
                    o.id = format!("{};{rs}", o.id);
                }
                put(a, "DB", Val::Flag);
            }
        }
    }
}

/// `GenotypeUtils.computeDiploidGenotypeCounts`: [refs, hets, homs] (ExcessHet rounds each genotype's
/// contribution, InbreedingCoeff does not)
fn genotype_counts(gs: &[Gt], n_alleles: usize, round: bool) -> [f64; 3] {
    let mut counts = [0f64; 3]; // ref, het, hom
    for g in gs {
        let called = g.is_called();
        if !(g.pl.is_some() || (called && g.gq.is_some())) {
            continue;
        }
        let Some(pl) = &g.pl else {
            if g.is_hom_ref() {
                let gq = g.gq.unwrap_or(0);
                if round {
                    counts[0] += 1.0;
                } else if gq == 0 {
                    counts[0] += 1.0 / 3.0;
                    counts[1] += 1.0 / 3.0;
                    counts[2] += 1.0 / 3.0;
                } else {
                    // QualityUtils.qualToProb
                    let p = 1.0 - crate::jmath::pow(10.0, gq as f64 / -10.0);
                    counts[0] += p;
                    counts[1] += 1.0 - p;
                }
                continue;
            }
            refuse("ExcessHet for a called genotype without PL that is not hom-ref: not supported".into());
        };
        let norm = hc_genotype::mathutils::normalize_from_log10_to_linear(&gls_of(pl));
        let (r, h, v);
        if n_alleles > 2 {
            let mut m = 0;
            for i in 1..norm.len() {
                if norm[i] > norm[m] {
                    m = i;
                }
            }
            let [a1, a2] = genotype_alleles(m);
            if a1 != 0 && a2 != 0 {
                counts[2] += 1.0;
                continue;
            }
            let (mut max_l, mut het_i, mut var_i) = (norm[1], 1, 2);
            for a in 1..n_alleles {
                let hi = pl_index(0, a);
                if norm[hi] > max_l {
                    max_l = norm[hi];
                    het_i = hi;
                    var_i = pl_index(a, a);
                }
            }
            let b = [norm[0], norm[het_i], norm[var_i]];
            let s: f64 = hc_genotype::mathutils::sum(&b);
            r = b[0] / s;
            h = b[1] / s;
            v = b[2] / s;
        } else {
            r = norm[0];
            h = norm[1];
            v = norm[2];
        }
        if round {
            let fr = |d: f64| if d > 0.0 { (d + 0.5) as i32 } else { (d - 0.5) as i32 };
            counts[0] += fr(r) as f64;
            if r != h {
                counts[1] += fr(h) as f64;
            }
            if v != h {
                counts[2] += fr(v) as f64;
            }
        } else {
            counts[0] += r;
            counts[1] += h;
            counts[2] += v;
        }
    }
    counts
}

/// `ExcessHet.annotate`: the exact test on the genotype counts
fn excess_het(gs: &[Gt], counts: [f64; 3]) -> Option<String> {
    let samples = gs.iter().filter(|g| g.is_called() && (g.pl.is_some() || g.gq.is_some())).count();
    if gs.is_empty() || samples < 1 {
        return None;
    }
    let p = hc_genotype::annotate::exact_test(counts[1] as i32, counts[0] as i32, counts[2] as i32);
    let eh = if p < 10e-60 { 160.0 } else { -10.0 * crate::jmath::log10(p) + 0.0 };
    Some(f(eh, 4))
}

/// `AS_InbreedingCoeff.annotate` with `HeterozygosityCalculator` (unrounded likelihood counts), per alt "%.4f"
fn as_inbreeding_coeff(gs: &[Gt], n_alleles: usize) -> Option<String> {
    let mut het = vec![0f64; n_alleles];
    let mut ac = vec![0f64; n_alleles];
    let mut samples = 0usize;
    for g in gs {
        if !(g.is_called() && (g.pl.is_some() || g.gq.is_some())) {
            continue;
        }
        samples += 1;
        let gls = match &g.pl {
            Some(pl) => gls_of(pl),
            None => approx_gls_from_gq(g.gq.unwrap(), n_alleles),
        };
        let l = hc_genotype::mathutils::normalize_from_log10_to_linear(&gls);
        for k in 1..n_alleles {
            for i in 0..n_alleles {
                if i == k {
                    ac[k] += 2.0 * l[pl_index(k, k)];
                    continue;
                }
                let ab = pl_index(i.min(k), i.max(k));
                het[k] += l[ab];
                ac[k] += l[ab];
                ac[0] += l[ab];
            }
            ac[0] += 2.0 * l[0];
        }
    }
    if samples < 10 {
        return None;
    }
    let an = gs.iter().flat_map(|g| g.al.iter()).filter(|&&x| x != NO_CALL).count() as f64;
    let v: Vec<String> = (1..n_alleles)
        .map(|k| {
            let ic = if n_alleles == 2 {
                let ref_af = ac[0] / (ac[k] + ac[0]);
                let alt_af = 1.0 - ref_af;
                1.0 - (het[k] / (2.0 * ref_af * alt_af * samples as f64))
            } else {
                let alt_af = ac[k] / an;
                1.0 - (het[k] / (2.0 * (1.0 - alt_af) * alt_af * samples as f64))
            };
            f(ic, 4)
        })
        .collect();
    Some(v.join(","))
}

/// merge + regenotypeVC + annotations + cleanup for one input record
fn genotype_record(rec: &MRec, gt: &mut Genotyper, ann: &mut Annotator) -> Option<(MRec, Attrs)> {
    let (alleles, gs) = merge_genotypes(rec);
    let depth = vc_depth(rec);
    let merged = MRec {
        contig: rec.contig.clone(),
        start: rec.start,
        stop: rec.start + alleles[0].len() as i64 - 1,
        id: {
            let ids: Vec<&str> = rec.id.split(',').filter(|x| *x != ".").collect();
            if ids.is_empty() { ".".into() } else { ids.join(",") }
        },
        alleles: alleles.clone(),
        log10_p_error: None,
        filter: ".".into(),
        info: Vec::new(),
        gs,
    };
    let mut attrs = ann.merge_attrs(rec, &alleles, depth);
    if !(merged.alleles.len() > 1 && depth > 0) {
        return None;
    }
    let call = gt.calculate(&merged)?;
    let mut o = merged.clone();
    o.alleles = call.alleles.iter().map(|&a| merged.alleles[a].clone()).collect();
    // isProperlyPolymorphic
    let alts = &o.alleles[1..];
    let poly = if alts.is_empty() { false } else if alts.len() == 1 { alts[0] != SPAN_DEL } else { !(alts[0] == SPAN_DEL && alts[1] == NON_REF) };
    if !poly {
        return None;
    }
    o.gs = call.gs;
    o.log10_p_error = Some(call.log10_conf);
    put(&mut attrs, "MLEAC", Val::IL(call.mleac.iter().map(|&x| x as i64).collect()));
    put(&mut attrs, "MLEAF", Val::DL(call.mleaf.clone()));
    ann.finalize(&mut attrs, &alleles, &o.alleles, &o.gs, &call.as_qual);
    if crate::reblock::reverse_trim(&mut o.alleles) > 0 {
        o.stop = o.start + o.alleles[0].len() as i64 - 1;
    }
    // isPolymorphicInSamples
    if !o.gs.iter().any(|g| g.al.iter().any(|&x| x != NO_CALL && x != 0)) {
        return None;
    }
    for g in o.gs.iter_mut() {
        if exclude_from_annotations(g) {
            g.al = vec![NO_CALL; g.ploidy()];
        }
    }
    ann.annotate(&mut attrs, &mut o, call.log10_conf);
    // cleanupGenotypeAnnotations(createRefGTs = false, keepSB = false)
    let n_alleles = o.alleles.len();
    for g in o.gs.iter_mut() {
        let mut depth = g.dp.unwrap_or(0) as i64;
        if let Some(pos) = g.ext.iter().position(|(k, _)| k == "MIN_DP") {
            depth = g.ext[pos].1.parse().unwrap();
            g.dp = Some(depth as i32);
            g.ext.remove(pos);
        }
        g.remove_ext("SB");
        if g.is_hom_var() && g.ext.iter().any(|(k, _)| k == "PGT") {
            g.set_ext("PGT", "1|1".into());
        }
        if g.ad.is_none() && depth > 0 {
            let mut ad = vec![0i32; n_alleles];
            ad[0] = depth as i32;
            g.ad = Some(ad);
        }
    }
    Some((o, attrs))
}

fn encode(o: &MRec, a: &Attrs, order: &[usize], out: &mut String) {
    use std::fmt::Write;
    let qual = o.log10_p_error.map(|l| hc_genotype::javamath::format_qual(-10.0 * l)).unwrap_or_else(|| ".".into());
    let _ = write!(out, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t", o.contig, o.start, o.id, o.alleles[0], o.alleles[1..].join(","), qual, o.filter);
    let m: std::collections::BTreeMap<&str, String> = a.iter().map(|(k, v)| (k.as_str(), v.text())).collect();
    if m.is_empty() {
        out.push('.');
    } else {
        out.push_str(&m.iter().map(|(k, v)| if v.is_empty() { k.to_string() } else { format!("{k}={v}") }).collect::<Vec<_>>().join(";"));
    }
    let gs: Vec<Gt> = order.iter().map(|&i| o.gs[i].clone()).collect();
    crate::reblock::encode_genotypes(&gs, out);
}
