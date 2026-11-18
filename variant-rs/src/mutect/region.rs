//! Mutect2 region processing (stage S2): `Mutect2Engine.callRegion` from `removeUnmarkedDuplicates` to the
//! read likelihoods after realignment (natural log), i.e. everything `SomaticGenotypingEngine.callMutations` takes.
//!
//! * reads: the region's reads with their qualities at `callRegion` time (after the in-place
//!   `cleanOverlappingReadPairs` of this and earlier regions; BI/BD then exist, capped at 20 in mate overlaps)
//! * `removeUnmarkedDuplicates`, `finalizeRegion` without overlap correction, local assembly with
//!   `AdaptiveChainPruner` (prune factor 0) and per-sample read threading, events with `maxMnpDistance = 1`,
//!   trimming, `removeReadStubs` (< 10 bases; no `filterNonPassingReads`), PairHMM per sample with the reads'
//!   BI/BD as gap-open penalties, `normalizeLikelihoods` + `filterPoorlyModeledEvidence` per sample,
//!   `switchToNaturalLog`, `realignReadsToTheirBestHaplotype`.

use crate::frontend::{FilledRegion, M2Quals, RefSeq};
use hc_assembly::{Haplotype, Read};
use std::collections::HashMap;
use std::sync::Arc;

pub const REF_PAD_FOR_ASSEMBLY: i64 = 500;
pub const HUGE_FRAGMENT_LENGTH: i64 = 1_000_000;

/// One sample's read likelihoods (natural log, `[haplotype][read]`) after realignment.
pub struct SampleLikelihoods {
    pub reads: Vec<Read>,
    pub filtered: Vec<Read>,
    pub matrix: Vec<Vec<f64>>,
    /// per read: best haplotype, informative, realigned (`realignReadsToTheirBestHaplotype`)
    pub realign: Vec<hc_genotype::realign::Realigned>,
}

pub struct RegionOut {
    pub exit: &'static str,
    pub raw_reads: Vec<Read>,
    pub dedup_reads: Vec<Read>,
    pub final_reads: Vec<Read>,
    pub ref_padded: (i64, i64),
    pub haplotypes: Vec<Haplotype>,
    pub events: Vec<hc_likelihood::events::Event>,
    pub variant_region: Option<(i64, i64, i64, i64)>,
    pub trimmed: Vec<Haplotype>,
    pub genotyping_reads: Vec<Read>,
    /// in `samples` (sorted) order
    pub likelihoods: Option<(Vec<Haplotype>, Vec<SampleLikelihoods>)>,
}

pub struct RegionCtx<'a> {
    pub refseq: &'a RefSeq,
    /// sample index -> rank in the sorted sample list (GATK `samplesList` order)
    pub sample_rank: &'a [usize],
    pub n_samples: usize,
    /// CUDA PairHMM (bit-identical to the CPU port of GKL's AVX-512 kernel)
    pub use_gpu: bool,
}

fn to_asm(r: hc_preproc::Read) -> Read {
    Read {
        name: r.name,
        flag: r.flag,
        start: r.start,
        end: r.end,
        cigar: r.cigar,
        seq: r.seq,
        qual: r.qual,
        mapq: r.mapq,
        mate_start: r.mate_start,
        tlen: r.tlen,
        ref_id: r.ref_id,
        mate_ref_id: r.mate_ref_id,
        bi: r.bi,
        bd: r.bd,
        sample: r.sample,
        oa_other_contig: r.oa_other_contig,
        bam_rec: r.bam_rec,
    }
}

/// `Mutect2Engine.removeUnmarkedDuplicates` on the region's reads (in region order).
fn remove_unmarked_duplicates(reads: Vec<hc_preproc::Read>) -> Vec<hc_preproc::Read> {
    let unclipped_start = |r: &hc_preproc::Read| -> i64 {
        let mut s = r.start;
        for &(op, len) in &r.cigar {
            if op == b'S' || op == b'H' {
                s -= len as i64;
            } else {
                break;
            }
        }
        s
    };
    // groups by (sample, +/- unclipped start), in encounter order
    let mut groups: HashMap<(u16, i64), Vec<usize>> = HashMap::new();
    for (i, r) in reads.iter().enumerate() {
        let paired = r.flag & 0x1 != 0;
        let mate_unmapped = r.flag & 0x8 != 0 || r.mate_start == 0;
        if !(paired && !mate_unmapped && (r.mate_ref_id != r.ref_id || r.tlen.abs() > HUGE_FRAGMENT_LENGTH)) {
            continue;
        }
        let sign = if r.flag & 0x40 != 0 { 1 } else { -1 };
        groups.entry((r.sample, sign * unclipped_start(r))).or_default().push(i);
    }
    let mut drop = vec![false; reads.len()];
    for list in groups.values() {
        let mut by_contig: HashMap<i32, Vec<usize>> = HashMap::new();
        for &i in list {
            by_contig.entry(reads[i].mate_ref_id).or_default().push(i);
        }
        for cr in by_contig.values() {
            let skip = if cr.len() > list.len() / 2 { 1 } else { 0 };
            for &i in &cr[skip..] {
                drop[i] = true;
            }
        }
    }
    reads.into_iter().zip(drop).filter(|x| !x.1).map(|x| x.0).collect()
}

/// Accumulated nanoseconds per stage of `process_region` (+ genotyping), for `VARIANT_RS_TRACE_STAGES`.
pub static STAGE_NS: [std::sync::atomic::AtomicU64; 7] = [const { std::sync::atomic::AtomicU64::new(0) }; 7];
pub const STAGE_NAMES: [&str; 7] = ["preprocess", "assembly", "prepare", "pairhmm", "finish", "realign", "genotype"];
pub fn stage_add(i: usize, t: std::time::Instant) {
    STAGE_NS[i].fetch_add(t.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
}

pub fn process_region(
    fr: &FilledRegion,
    quals: &[Option<Arc<M2Quals>>],
    ctx: &RegionCtx,
    asm: &mut hc_assembly::Assembler,
    sw: &mut hc_assembly::sw::SwAligner,
    lp: &hc_likelihood::Params,
) -> RegionOut {
    let contig_len = ctx.refseq.contig_len;
    let t = std::time::Instant::now();
    // reads with their qualities at callRegion time
    let raw: Vec<hc_preproc::Read> = fr
        .reads
        .iter()
        .zip(quals)
        .map(|(r, q)| {
            let mut p = r.to_preproc();
            if let Some(q) = q {
                p.qual = q.qual.clone();
                p.bi = q.bi.clone();
                p.bd = q.bd.clone();
            }
            p
        })
        .collect();
    let raw_reads: Vec<Read> = raw.iter().cloned().map(to_asm).collect();
    let dedup = remove_unmarked_duplicates(raw);
    let dedup_reads: Vec<Read> = dedup.iter().cloned().map(to_asm).collect();
    let fp = hc_preproc::FinalizeParams { correct_overlapping_base_qualities: false, ..Default::default() };
    let finals: Vec<Read> = hc_preproc::finalize_region_owned_with(dedup, fr.padded_start, fr.padded_end, &fp).into_iter().map(to_asm).collect();
    let (rps, rpe) = ((fr.padded_start - REF_PAD_FOR_ASSEMBLY).max(1), (fr.padded_end + REF_PAD_FOR_ASSEMBLY).min(contig_len));
    let ref_padded = ctx.refseq.window(rps, rpe);
    stage_add(0, t);
    let t = std::time::Instant::now();
    let info = hc_assembly::RegionInfo { start: fr.start, end: fr.end, padded_start: fr.padded_start, padded_end: fr.padded_end };
    let a = asm.try_assemble_region(&finals, info, ref_padded, rps, rpe).unwrap_or_else(|e| panic!("assembly failed at {}-{}: {e}", fr.start, fr.end));
    stage_add(1, t);
    let t = std::time::Instant::now();
    let events = hc_likelihood::events::variation_events(&a.haplotypes, ref_padded, rps, lp.max_mnp_distance);
    let mut out = RegionOut {
        exit: "",
        raw_reads,
        dedup_reads,
        final_reads: finals.clone(),
        ref_padded: (rps, rpe),
        haplotypes: a.haplotypes.clone(),
        events,
        variant_region: None,
        trimmed: Vec::new(),
        genotyping_reads: Vec::new(),
        likelihoods: None,
    };
    let prep = hc_likelihood::prepare_with(
        &hc_likelihood::RegionInput {
            span: hc_likelihood::Span { start: fr.start, end: fr.end },
            padded_span: hc_likelihood::Span { start: fr.padded_start, end: fr.padded_end },
            reads: &finals,
            mate_on_other_contig: &[],
            haplotypes: &a.haplotypes,
            ref_padded,
            ref_padded_start: rps,
        },
        lp,
        None,
    );
    if let (Some(v), Some(vp)) = (prep.variant_span, prep.variant_padded_span) {
        out.variant_region = Some((v.start, v.end, vp.start, vp.end));
    }
    stage_add(2, t);
    out.trimmed = prep.haplotypes.clone();
    out.genotyping_reads = prep.reads.clone();
    match prep.exit {
        hc_likelihood::Exit::NoVariationAfterTrim => {
            out.exit = "no_variation_after_trim";
            return out;
        }
        hc_likelihood::Exit::NoAssemblyVariation => {
            out.exit = "no_assembly_variation";
            return out;
        }
        _ => {}
    }
    // per sample (sorted order): PairHMM, normalization, poorly-modelled filter, natural log, realignment
    let mut samples = Vec::with_capacity(ctx.n_samples);
    for rank in 0..ctx.n_samples {
        let keep: Vec<bool> = prep.reads.iter().map(|r| ctx.sample_rank[r.sample as usize] == rank).collect();
        let sub = prep.subset(&keep);
        let t = std::time::Instant::now();
        let raw = if ctx.use_gpu {
            hc_likelihood::compute_likelihoods_gpu(&[&sub]).pop().unwrap()
        } else {
            hc_likelihood::compute_likelihoods_cpu(&[&sub], hc_likelihood::cpu::Flavor::Avx512).pop().unwrap()
        };
        stage_add(3, t);
        let t = std::time::Instant::now();
        let res = hc_likelihood::finish_with(sub, &raw, lp, None);
        let lk = res.likelihoods.expect("likelihoods");
        let n = lk.reads.len();
        let mut mat: Vec<Vec<f64>> = if n == 0 { vec![Vec::new(); lk.haplotypes.len()] } else { lk.values.chunks(n).map(|c| c.to_vec()).collect() };
        for row in mat.iter_mut() {
            for x in row.iter_mut() {
                *x *= std::f64::consts::LN_10; // MathUtils.log10ToLog
            }
        }
        stage_add(4, t);
        let t = std::time::Instant::now();
        let mut reads = lk.reads;
        // AlleleLikelihoods.NATURAL_LOG_INFORMATIVE_THRESHOLD = log10ToLog(0.2)
        let realign = hc_genotype::realign::realign_reads_to_best_haplotype_t(sw, &mut reads, &lk.haplotypes, &mat, rps, 0.2 * std::f64::consts::LN_10);
        stage_add(5, t);
        samples.push(SampleLikelihoods { reads, filtered: lk.filtered_reads, matrix: mat, realign });
    }
    out.likelihoods = Some((prep.haplotypes.clone(), samples));
    out.exit = "genotyped";
    out
}

// ------------------------------------------------------------------------------------------------
// JSON dump in the schema of GATK's M2Dumper (validation)
// ------------------------------------------------------------------------------------------------

fn json_str(s: &[u8]) -> String {
    let mut o = String::with_capacity(s.len() + 2);
    for &c in s {
        match c {
            b'"' => o.push_str("\\\""),
            b'\\' => o.push_str("\\\\"),
            c => o.push(c as char),
        }
    }
    o
}

fn quals33(q: &[u8], n: usize, default: u8) -> String {
    if q.is_empty() {
        return std::iter::repeat((default + 33) as char).take(n).collect();
    }
    json_str(&q.iter().map(|&b| b + 33).collect::<Vec<u8>>())
}

fn reads_json(reads: &[Read], samples: &[String]) -> String {
    let mut s = String::from("[");
    for (i, r) in reads.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        s.push_str(&format!(
            "{{\"name\":\"{}\",\"sample\":\"{}\",\"flag\":{},\"start\":{},\"end\":{},\"cigar\":\"{}\",\"seq\":\"{}\",\"qual\":\"{}\",\"bi\":\"{}\",\"bd\":\"{}\",\"mapq\":{}}}",
            r.name,
            samples[r.sample as usize],
            r.flag,
            r.start,
            r.end,
            hc_likelihood::cigar_string(&r.cigar),
            String::from_utf8_lossy(&r.seq),
            quals33(&r.qual, r.seq.len(), 0),
            quals33(&r.bi, r.seq.len(), 45),
            quals33(&r.bd, r.seq.len(), 45),
            r.mapq
        ));
    }
    s.push(']');
    s
}

fn hap_json(h: &Haplotype) -> String {
    let score = if h.score.is_finite() { java_double(h.score) } else { "null".into() };
    format!(
        "{{\"bases\":\"{}\",\"cigar\":\"{}\",\"align_start_wrt_ref\":{},\"loc_start\":{},\"loc_end\":{},\"is_ref\":{},\"score\":{},\"kmer\":{}}}",
        String::from_utf8_lossy(&h.bases),
        hc_likelihood::cigar_string(&h.cigar),
        h.align_start_wrt_ref,
        h.loc_start,
        h.loc_end,
        h.is_ref,
        score,
        h.kmer
    )
}

fn haps_json(hs: &[Haplotype]) -> String {
    format!("[{}]", hs.iter().map(hap_json).collect::<Vec<_>>().join(","))
}

/// Java `Double.toString` (shortest repr; the comparer parses numbers, so only the value matters).
fn java_double(x: f64) -> String {
    if x.is_nan() {
        "\"NaN\"".into()
    } else if x.is_infinite() {
        if x > 0.0 { "\"Infinity\"".into() } else { "\"-Infinity\"".into() }
    } else {
        format!("{x:?}")
    }
}

pub fn dump_json(fr: &FilledRegion, o: &RegionOut, samples: &[String]) -> String {
    let mut s = format!(
        "{{\"start\":{},\"end\":{},\"padded_start\":{},\"padded_end\":{},\"raw_reads\":{},\"dedup_reads\":{},\"final_reads\":{},\"assembly\":{{\"ref_padded_start\":{},\"ref_padded_end\":{},\"haplotypes\":{}}},\"events\":[{}]",
        fr.start,
        fr.end,
        fr.padded_start,
        fr.padded_end,
        reads_json(&o.raw_reads, samples),
        reads_json(&o.dedup_reads, samples),
        reads_json(&o.final_reads, samples),
        o.ref_padded.0,
        o.ref_padded.1,
        haps_json(&o.haplotypes),
        o.events.iter().map(|e| format!("[{},\"{}\",\"{}\"]", e.start, String::from_utf8_lossy(&e.ref_bases), String::from_utf8_lossy(&e.alt_bases))).collect::<Vec<_>>().join(",")
    );
    if let Some((a, b, c, d)) = o.variant_region {
        s.push_str(&format!(",\"variant_region\":[{a},{b},{c},{d}]"));
    }
    if o.exit != "no_variation_after_trim" {
        s.push_str(&format!(",\"trimmed\":{{\"ref_padded_start\":{},\"ref_padded_end\":{},\"haplotypes\":{}}}", o.ref_padded.0, o.ref_padded.1, haps_json(&o.trimmed)));
    }
    if o.exit == "genotyped" {
        s.push_str(&format!(",\"genotyping_reads\":{}", reads_json(&o.genotyping_reads, samples)));
        let (haps, per) = o.likelihoods.as_ref().unwrap();
        let mut sorted: Vec<(usize, &String)> = samples.iter().enumerate().collect();
        sorted.sort_by(|a, b| a.1.cmp(b.1));
        s.push_str(&format!(",\"likelihoods\":{{\"haplotypes\":{},\"samples\":[", haps_json(haps)));
        for (k, sl) in per.iter().enumerate() {
            if k > 0 {
                s.push(',');
            }
            let m: Vec<String> = sl.matrix.iter().map(|row| format!("[{}]", row.iter().map(|&x| java_double(x)).collect::<Vec<_>>().join(","))).collect();
            s.push_str(&format!(
                "{{\"sample\":\"{}\",\"reads\":{},\"filtered\":{},\"matrix\":[{}]}}",
                sorted[k].1,
                reads_json(&sl.reads, samples),
                reads_json(&sl.filtered, samples),
                m.join(",")
            ));
        }
        s.push_str("]}");
    }
    s.push_str(&format!(",\"exit\":\"{}\"}}", o.exit));
    s
}
