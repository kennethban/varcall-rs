//! Per-contig back end: batches of filled regions through finalizeRegion, assembly, likelihoods
//! (PairHMM on the GPU, shared by all contigs through `hc_cudahmm`'s lock), realignment and genotyping
//! (sequential, in genomic order: the engine carries cross-region state), then VCF text or, in GVCF
//! mode, reference confidence and GVCF banding on a writer thread.

use crate::bgzf::BgzfWriter;
use crate::frontend::{FilledRegion, RefSeq, REGION_PADDING};
use crate::tabix::{ContigIndex, VcfIndexer};
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use std::sync::mpsc::Receiver;
use std::sync::Arc;
use std::time::Instant;
use hc_downsample::JavaRandom;

pub const ASM_BATCH: usize = 256; // active regions per batch
const REF_PAD_FOR_ASSEMBLY: i64 = 500;

#[derive(Default, Clone, Debug)]
pub struct Stats {
    pub convert: f64,
    pub preproc: f64,
    pub asm: f64,
    pub prep: f64,
    pub hmm: f64,
    pub finish: f64,
    pub realign: f64,
    pub geno: f64,
    pub refconf: f64,
    pub write: f64,
    pub regions: u64,
    pub active: u64,
    pub genotyped: u64,
    pub calls: u64,
    pub records: u64,
    pub haplotypes: u64,
    pub final_reads: u64,
    pub qd_jitter: u64,
    pub qd_drawn: u64,
    pub ds_draws: u64,
    pub ds_groups: u64,
    pub ds_kept: u64,
    pub ds_in: u64,
    pub out_lines: u64,
    pub writer_busy: f64,
    pub deflate: f64,
}

impl Stats {
    pub fn add(&mut self, o: &Stats) {
        self.convert += o.convert;
        self.preproc += o.preproc;
        self.asm += o.asm;
        self.prep += o.prep;
        self.hmm += o.hmm;
        self.finish += o.finish;
        self.realign += o.realign;
        self.geno += o.geno;
        self.refconf += o.refconf;
        self.write += o.write;
        self.regions += o.regions;
        self.active += o.active;
        self.genotyped += o.genotyped;
        self.calls += o.calls;
        self.records += o.records;
        self.haplotypes += o.haplotypes;
        self.final_reads += o.final_reads;
        self.qd_jitter += o.qd_jitter;
        self.qd_drawn += o.qd_drawn;
        self.ds_draws += o.ds_draws;
        self.ds_groups += o.ds_groups;
        self.ds_kept += o.ds_kept;
        self.ds_in += o.ds_in;
        self.out_lines += o.out_lines;
        self.writer_busy += o.writer_busy;
        self.deflate += o.deflate;
    }
}

/// A finished per-contig output file: BGZF (with block sizes and a tabix index in pseudo offsets) or plain.
pub struct Part {
    pub path: PathBuf,
    pub block_sizes: Option<Vec<u32>>,
    pub index: Vec<ContigIndex>,
}

/// Output writer of a part: bgzip + tabix bookkeeping, or plain text.
pub enum PartWriter {
    Gz(BgzfWriter<BufWriter<File>, VcfIndexer>, PathBuf),
    Plain(BufWriter<File>, PathBuf),
}

impl PartWriter {
    pub fn create(path: PathBuf, level: Option<i32>) -> io::Result<PartWriter> {
        let f = BufWriter::with_capacity(1 << 20, File::create(&path)?);
        Ok(match level {
            Some(l) => PartWriter::Gz(BgzfWriter::new(f, l, VcfIndexer::default(), true), path),
            None => PartWriter::Plain(f, path),
        })
    }
    pub fn finish(self) -> io::Result<(Part, f64)> {
        match self {
            PartWriter::Gz(w, path) => {
                let (_, sizes, mut idx, secs) = w.finish()?;
                idx.finish();
                Ok((Part { path, block_sizes: Some(sizes), index: idx.contigs }, secs))
            }
            PartWriter::Plain(mut w, path) => {
                w.flush()?;
                Ok((Part { path, block_sizes: None, index: Vec::new() }, 0.0))
            }
        }
    }
}

impl Write for PartWriter {
    fn write(&mut self, b: &[u8]) -> io::Result<usize> {
        match self {
            PartWriter::Gz(w, _) => w.write(b),
            PartWriter::Plain(w, _) => w.write(b),
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// HaplotypeCaller options that change the records (shared by every contig).
pub struct HcOpts {
    /// `--gvcf-gq-bands`
    pub gq_bands: Vec<i32>,
    /// `--dbsnp`
    pub dbsnp: Option<crate::dbsnp::Dbsnp>,
    /// `-contamination` (0: off)
    pub contamination: f64,
    /// `-G AS_StandardAnnotation`
    pub as_annotations: bool,
    /// number of samples (several: VCF mode)
    pub n_samples: usize,
}

/// Contamination downsampling draws from GATK's `Well19937c` in genotyping order: a contig takes it before its first
/// genotyping (after every earlier contig of the domain is done) and hands it on at its end.
fn hold_contam_rng(cfg: &BackendCfg, engine: &mut hc_genotype::GenotypingEngine) {
    if cfg.opts.contamination > 0.0 && engine.contam_rng.is_none() {
        engine.contam_rng = Some(cfg.contam.take(cfg.domain_idx));
    }
}

fn release_contam_rng(cfg: &BackendCfg, engine: &mut hc_genotype::GenotypingEngine) {
    match engine.contam_rng.take() {
        Some(r) => cfg.contam.finish(cfg.domain_idx, r),
        None => cfg.contam.pass(cfg.domain_idx),
    }
}

pub struct BackendCfg {
    pub contig: String,
    pub tid: i32,
    pub contig_len: i64,
    pub refseq: Arc<RefSeq>,
    pub use_gpu: bool,
    pub erc: bool,
    /// GVCF mode: banded output part (and optional raw records part)
    pub gvcf_part: Option<(PathBuf, Option<i32>)>,
    pub raw_part: Option<PathBuf>,
    pub domain: Arc<crate::rng::DomainRng>,
    pub domain_idx: usize,
    pub opts: Arc<HcOpts>,
    /// the domain's contamination-downsampling generator
    pub contam: Arc<crate::rng::DomainRng<hc_genotype::well::Well19937c>>,
}

pub struct BackendOut {
    /// VCF mode: the data lines
    pub vcf_text: Vec<u8>,
    /// false: QD jitter placeholders remain (the contig never got the generator; drawn at output)
    pub vcf_resolved: bool,
    pub gvcf: Option<Part>,
    pub raw: Option<PathBuf>,
    pub stats: Stats,
}

/// Converts and finalizes (`finalizeRegion`) the reads of each region, in parallel.
fn finalize_regions(regions: &[FilledRegion], keep: &[bool]) -> Vec<Vec<hc_assembly::Read>> {
    regions
        .par_iter()
        .zip(keep.par_iter())
        .map(|(fr, &k)| {
            if !k {
                return Vec::new();
            }
            crate::cpuprof::measure(0, || {
                let raw: Vec<hc_preproc::Read> = fr.reads.iter().map(|r| r.to_preproc()).collect();
                hc_preproc::finalize_region_owned(raw, fr.padded_start, fr.padded_end).into_iter().map(to_asm_read).collect()
            })
        })
        .collect()
}

fn to_asm_read(r: hc_preproc::Read) -> hc_assembly::Read {
    hc_assembly::Read {
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

fn pair_hmm(preps: &[&hc_likelihood::Prepared], use_gpu: bool) -> Vec<Vec<f64>> {
    if use_gpu {
        crate::cpuprof::measure(3, || hc_likelihood::compute_likelihoods_gpu(preps))
    } else {
        let chunk = (preps.len() / (rayon::current_num_threads() * 4)).max(1);
        preps
            .par_chunks(chunk)
            .flat_map_iter(|c| crate::cpuprof::measure(3, || hc_likelihood::compute_likelihoods_cpu(c, hc_likelihood::cpu::Flavor::Avx512)))
            .collect()
    }
}

/// Assembly through realignment for the regions with `run[i]`; returns per region the likelihood result and
/// the realigned evidence (only for genotyped regions).
#[allow(clippy::type_complexity)]
fn assemble_and_score(
    regions: &[FilledRegion],
    finals: &[Vec<hc_assembly::Read>],
    run: &[bool],
    refseq: &RefSeq,
    contig_len: i64,
    use_gpu: bool,
    st: &mut Stats,
) -> (Vec<Option<hc_likelihood::RegionResult>>, Vec<Option<(Vec<Vec<f64>>, Vec<hc_assembly::Read>)>>) {
    let rp = |fr: &FilledRegion| ((fr.padded_start - REF_PAD_FOR_ASSEMBLY).max(1), (fr.padded_end + REF_PAD_FOR_ASSEMBLY).min(contig_len));
    // 1. local assembly
    let t = Instant::now();
    let asm: Vec<Option<Result<hc_assembly::AssemblyResult, String>>> = regions
        .par_iter()
        .zip(finals.par_iter())
        .zip(run.par_iter())
        .map_init(
            || hc_assembly::Assembler::new(hc_assembly::AssemblerConfig::default()),
            |a, ((fr, f), &go)| {
                if !go {
                    return None;
                }
                let (rps, rpe) = rp(fr);
                let info = hc_assembly::RegionInfo { start: fr.start, end: fr.end, padded_start: fr.padded_start, padded_end: fr.padded_end };
                Some(crate::cpuprof::measure(1, || a.try_assemble_region(f, info, refseq.window(rps, rpe), rps, rpe)))
            },
        )
        .collect();
    st.asm += t.elapsed().as_secs_f64();
    for (fr, r) in regions.iter().zip(asm.iter()) {
        match r {
            Some(Ok(a)) => st.haplotypes += a.haplotypes.len() as u64,
            // GATK rethrows an assembly exception (AssemblyBasedCallerUtils.assembleReads) and the run fails; so
            // does this one, before any output is finalized (the panic hook exits): skipping the region gave a
            // complete-looking VCF missing its calls
            Some(Err(e)) => panic!("assembly failed at {}-{}: {e}", fr.start, fr.end),
            None => {}
        }
    }
    // 2. trimming + genotyping reads + PairHMM inputs
    let t = Instant::now();
    let preps: Vec<Option<hc_likelihood::Prepared>> = regions
        .par_iter()
        .zip(finals.par_iter())
        .zip(asm.par_iter())
        .map(|((fr, f), a)| {
            let a = a.as_ref()?.as_ref().ok()?;
            let (rps, rpe) = rp(fr);
            Some(crate::cpuprof::measure(2, || hc_likelihood::prepare(&hc_likelihood::RegionInput {
                span: hc_likelihood::Span { start: fr.start, end: fr.end },
                padded_span: hc_likelihood::Span { start: fr.padded_start, end: fr.padded_end },
                reads: f,
                mate_on_other_contig: &[], // taken from the reads' ref_id / mate_ref_id
                haplotypes: &a.haplotypes,
                ref_padded: refseq.window(rps, rpe),
                ref_padded_start: rps,
            })))
        })
        .collect();
    drop(asm);
    st.prep += t.elapsed().as_secs_f64();
    // 3. PairHMM for the whole batch in one call (one GPU launch)
    let t = Instant::now();
    let refs: Vec<&hc_likelihood::Prepared> = preps.iter().flatten().collect();
    let raws = pair_hmm(&refs, use_gpu);
    st.hmm += t.elapsed().as_secs_f64();
    // 4. likelihood post-processing + realignment of reads to their best haplotype
    let t = Instant::now();
    let mut raw_it = raws.into_iter();
    let paired: Vec<Option<(hc_likelihood::Prepared, Vec<f64>)>> = preps.into_iter().map(|p| p.map(|p| (p, raw_it.next().unwrap()))).collect();
    let finished: Vec<Option<hc_likelihood::RegionResult>> = paired.into_par_iter().map(|o| o.map(|(p, r)| crate::cpuprof::measure(4, || hc_likelihood::finish(p, &r)))).collect();
    st.finish += t.elapsed().as_secs_f64();
    let t = Instant::now();
    let realigned: Vec<Option<(Vec<Vec<f64>>, Vec<hc_assembly::Read>)>> = finished
        .par_iter()
        .zip(regions.par_iter())
        .map_init(hc_assembly::sw::SwAligner::new, |sw, (res, fr)| crate::cpuprof::measure(5, || {
            let res = res.as_ref()?;
            if res.exit != hc_likelihood::Exit::Genotyped {
                return None;
            }
            let lk = res.likelihoods.as_ref()?;
            let n = lk.reads.len();
            let mat: Vec<Vec<f64>> = if n == 0 { vec![Vec::new(); lk.haplotypes.len()] } else { lk.values.chunks(n).map(|c| c.to_vec()).collect() };
            let mut reads = lk.reads.clone();
            hc_genotype::realign::realign_reads_to_best_haplotype(sw, &mut reads, &lk.haplotypes, &mat, rp(fr).0);
            Some((mat, reads))
        }))
        .collect();
    st.realign += t.elapsed().as_secs_f64();
    (finished, realigned)
}

/// VCF mode: one batch of active regions.
/// A batch after its per-region work (finalizeRegion, assembly, likelihoods, realignment): independent of every
/// other batch, so several batches can be in this stage at once.
pub struct Computed {
    regions: Vec<FilledRegion>,
    finals: Vec<Vec<hc_assembly::Read>>,
    finished: Vec<Option<hc_likelihood::RegionResult>>,
    realigned: Vec<Option<(Vec<Vec<f64>>, Vec<hc_assembly::Read>)>>,
    st: Stats,
}

/// The per-region stages of a batch (VCF: every region; GVCF: assembly only for active regions with reads).
fn compute_batch(regions: Vec<FilledRegion>, cfg: &BackendCfg) -> Computed {
    let mut st = Stats::default();
    if regions.is_empty() {
        return Computed { regions, finals: Vec::new(), finished: Vec::new(), realigned: Vec::new(), st };
    }
    let t = Instant::now();
    let keep = vec![true; regions.len()];
    let finals = finalize_regions(&regions, &keep);
    st.preproc += t.elapsed().as_secs_f64();
    let run: Vec<bool> = if cfg.erc {
        // regions that go through assembly: active with at least one read (region.size() > 0)
        regions.iter().map(|fr| fr.active && !fr.reads.is_empty()).collect()
    } else {
        st.final_reads += finals.iter().map(|f| f.len() as u64).sum::<u64>();
        finals.iter().map(|_| true).collect()
    };
    let (finished, realigned) = assemble_and_score(&regions, &finals, &run, &cfg.refseq, cfg.contig_len, cfg.use_gpu, &mut st);
    Computed { regions, finals, finished, realigned, st }
}

/// VCF mode: genotyping, annotation and formatting of a computed batch, in genomic order.
fn genotype_vcf_batch(c: Computed, cfg: &BackendCfg, engine: &mut hc_genotype::GenotypingEngine, text: &mut QdText, st: &mut Stats) {
    st.add(&c.st);
    let Computed { regions, finished, realigned, .. } = c;
    if regions.is_empty() {
        return;
    }
    // genotyping + annotation + VCF formatting, in genomic order
    let t = Instant::now();
    let rp0 = |fr: &FilledRegion| (fr.padded_start - REF_PAD_FOR_ASSEMBLY).max(1);
    let rp1 = |fr: &FilledRegion| (fr.padded_end + REF_PAD_FOR_ASSEMBLY).min(cfg.contig_len);
    for ((fr, res), ra) in regions.iter().zip(finished.iter()).zip(realigned.iter()) {
        let (Some(res), Some((mat, reads))) = (res, ra) else { continue };
        let lk = res.likelihoods.as_ref().unwrap();
        let vs = res.variant_span.unwrap();
        let input = hc_genotype::RegionInput {
            haplotypes: &lk.haplotypes,
            reads: &lk.reads,
            matrix: mat,
            filtered_reads: &res.filtered_reads,
            filtered_evidence: &lk.filtered_reads,
            ref_padded: cfg.refseq.window(rp0(fr), rp1(fr)),
            ref_padded_start: rp0(fr),
            span_start: vs.start,
            span_end: vs.end,
        };
        let mut recs = engine.genotype_region(&input, reads);
        if let Some(db) = &cfg.opts.dbsnp {
            recs.iter_mut().for_each(|r| db.annotate(&cfg.contig, r));
        }
        st.genotyped += 1;
        st.records += recs.len() as u64;
        let mut t = Vec::with_capacity(recs.len() * 256);
        for r in &recs {
            t.extend_from_slice(r.to_line().as_bytes());
            t.push(b'\n');
        }
        text.add(fr.id, t);
    }
    st.geno += t.elapsed().as_secs_f64();
    st.regions += regions.len() as u64;
    st.active += regions.len() as u64;
}

/// What `callRegion` does with a region in reference-confidence mode.
enum RefPlan {
    /// `referenceModelForNoVariation(region)`
    Whole,
    /// left flank + `calculateRefConfidence` over the genotyping span with the calls + right flank
    Called { calls: Vec<hc_genotype::erc::ErcCall>, vspan: hc_likelihood::Span, vpadded: hc_likelihood::Span, reads: Vec<hc_assembly::Read> },
}

/// `AssemblyRegion.trim(span, padding = 100)` of a flank: reads hard-clipped to the new padded span, sorted.
fn flank_reads(reads: &[hc_assembly::Read], ps: i64, pe: i64) -> Vec<hc_assembly::Read> {
    let mut scratch = Vec::new();
    let mut out: Vec<hc_assembly::Read> = Vec::new();
    for r in reads {
        let mut c = r.clone();
        if hc_likelihood::clip::hard_clip_to_region(&mut c, ps, pe, &mut scratch) && !c.seq.is_empty() && hc_likelihood::clip::coord_overlaps(ps, pe, c.start, c.end) {
            out.push(c);
        }
    }
    out.sort_by(hc_likelihood::clip::read_coordinate_cmp);
    out
}

/// GVCF mode: one batch of regions (active or not), in order. Returns the records of every region.
/// GVCF mode: genotyping (in genomic order) and reference confidence of a computed batch.
fn genotype_erc_batch(c: Computed, cfg: &BackendCfg, engine: &mut hc_genotype::GenotypingEngine, st: &mut Stats) -> Vec<Vec<hc_refconf::Item>> {
    st.add(&c.st);
    let Computed { regions, finals, finished, realigned, .. } = c;
    if regions.is_empty() {
        return Vec::new();
    }
    let contig_len = cfg.contig_len;
    let refseq = &*cfg.refseq;
    let rp = |fr: &FilledRegion| ((fr.padded_start - REF_PAD_FOR_ASSEMBLY).max(1), (fr.padded_end + REF_PAD_FOR_ASSEMBLY).min(contig_len));
    // genotyping in genomic order (cross-region state), deciding what each region emits
    let t = Instant::now();
    let mut plans: Vec<RefPlan> = Vec::with_capacity(regions.len());
    for ((fr, res), ra) in regions.iter().zip(finished.into_iter()).zip(realigned.into_iter()) {
        let (Some(res), Some((mat, reads))) = (res, ra) else {
            plans.push(RefPlan::Whole);
            continue;
        };
        let lk = res.likelihoods.as_ref().unwrap();
        let (vs, vp) = (res.variant_span.unwrap(), res.variant_padded_span.unwrap());
        let (rps, rpe) = rp(fr);
        let input = hc_genotype::RegionInput {
            haplotypes: &lk.haplotypes,
            reads: &lk.reads,
            matrix: &mat,
            filtered_reads: &res.filtered_reads,
            filtered_evidence: &lk.filtered_reads,
            ref_padded: refseq.window(rps, rpe),
            ref_padded_start: rps,
            span_start: vs.start,
            span_end: vs.end,
        };
        let mut calls = engine.genotype_region_erc(&input, &reads);
        if let Some(db) = &cfg.opts.dbsnp {
            calls.iter_mut().for_each(|c| db.annotate(&cfg.contig, &mut c.rec));
        }
        st.genotyped += 1;
        st.calls += calls.len() as u64;
        // containsCalls: some genotype is called
        if calls.iter().any(|c| c.called) {
            plans.push(RefPlan::Called { calls, vspan: vs, vpadded: vp, reads });
        } else {
            plans.push(RefPlan::Whole);
        }
    }
    st.geno += t.elapsed().as_secs_f64();
    // reference confidence (parallel over regions)
    let t = Instant::now();
    let items: Vec<Vec<hc_refconf::Item>> = regions
        .par_iter()
        .zip(finals.par_iter())
        .zip(plans.par_iter())
        .map(|((j, f), plan)| {
            let sites = |reads: &[hc_assembly::Read], s: i64, e: i64, ps: i64, pe: i64| -> Vec<hc_refconf::Item> {
                hc_refconf::ref_model_for_no_variation(reads, s, e, ps, refseq.window(ps, pe)).into_iter().map(hc_refconf::Item::Site).collect()
            };
            match plan {
                RefPlan::Whole => sites(f, j.start, j.end, j.padded_start, j.padded_end),
                RefPlan::Called { calls, vspan, vpadded, reads } => {
                    let mut v: Vec<hc_refconf::Item> = Vec::with_capacity((j.end - j.start + 1) as usize);
                    // trim(flank, assemblyRegionPadding = 100): padded span = flank +/- 100 within the contig, intersected
                    let flank = |s: i64, e: i64| -> (i64, i64) { ((s - REGION_PADDING).max(1).max(j.padded_start), (e + REGION_PADDING).min(contig_len).min(j.padded_end)) };
                    if j.start < vspan.start {
                        let (s, e) = (j.start, vspan.start - 1);
                        let (ps, pe) = flank(s, e);
                        v.extend(sites(&flank_reads(f, ps, pe), s, e, ps, pe));
                    }
                    let spans: Vec<hc_refconf::CallSpan> = calls.iter().map(|c| hc_refconf::CallSpan { start: c.rec.pos, end: c.end }).collect();
                    let rr: Vec<&hc_assembly::Read> = reads.iter().collect();
                    for e in hc_refconf::calculate_ref_confidence(refseq.window(vpadded.start, vpadded.end), vpadded.start, vspan.start, vspan.end, &rr, &spans) {
                        v.push(match e {
                            hc_refconf::Entry::Site(s) => hc_refconf::Item::Site(s),
                            hc_refconf::Entry::Call(k) => hc_refconf::Item::Call(hc_refconf::CallRecord { pos: calls[k].rec.pos, end: calls[k].end, rec: calls[k].rec.clone() }),
                        });
                    }
                    if vspan.end < j.end {
                        let (s, e) = (vspan.end + 1, j.end);
                        let (ps, pe) = flank(s, e);
                        v.extend(sites(&flank_reads(f, ps, pe), s, e, ps, pe));
                    }
                    v
                }
            }
        })
        .collect();
    st.refconf += t.elapsed().as_secs_f64();
    st.regions += regions.len() as u64;
    st.active += regions.iter().filter(|j| j.active).count() as u64;
    items
}

/// What the compute stage passes to the in-order genotyping stage.
enum Staged {
    Done(Computed),
    Fwd(BackMsg),
}

/// Runs `consume` on the front end's messages in order, with every batch's per-region stages already done; up to
/// `VARIANT_RS_BACK_STAGES` (default 3) batches are in those stages at once (round-robin over compute threads, results read
/// back in the same order).
fn staged(cfg: &BackendCfg, rx: Receiver<BackMsg>, mut consume: impl FnMut(Staged)) {
    let k: usize = std::env::var("VARIANT_RS_BACK_STAGES").ok().and_then(|v| v.parse().ok()).unwrap_or(3).max(1);
    std::thread::scope(|scope| {
        let mut work = Vec::new();
        let mut results = Vec::new();
        for _ in 0..k {
            let (wtx, wrx) = std::sync::mpsc::sync_channel::<BackMsg>(1);
            let (rtx, rrx) = std::sync::mpsc::sync_channel::<Staged>(1);
            work.push(wtx);
            results.push(rrx);
            std::thread::Builder::new().name("back-compute".into()).spawn_scoped(scope, move || {
                for m in wrx {
                    let out = match m {
                        BackMsg::Batch(b) => Staged::Done(crate::cpuprof::measure(9, || compute_batch(b, cfg))),
                        other => Staged::Fwd(other),
                    };
                    if rtx.send(out).is_err() {
                        break;
                    }
                }
            }).unwrap();
        }
        std::thread::Builder::new().name("back-dispatch".into()).spawn_scoped(scope, move || {
            for (seq, m) in rx.into_iter().enumerate() {
                if work[seq % k].send(m).is_err() {
                    break;
                }
            }
        }).unwrap();
        let mut seq = 0usize;
        while let Ok(m) = results[seq % k].recv() {
            consume(m);
            seq += 1;
        }
    });
}

/// VCF text of a contig whose QD jitter draws are not all made yet: `done` is final text, `pending` the
/// regions (in order) whose placeholders still wait for the generator.
#[derive(Default)]
pub struct QdText {
    pub done: Vec<u8>,
    pending: std::collections::VecDeque<(usize, Vec<u8>)>,
    pub draws: u64,
}

impl QdText {
    fn add(&mut self, id: usize, t: Vec<u8>) {
        if self.pending.is_empty() && memchr::memchr(1, &t).is_none() {
            self.done.extend_from_slice(&t);
        } else {
            self.pending.push_back((id, t));
        }
    }
    /// Makes the QD draws of the regions with id < `before`.
    fn resolve(&mut self, before: usize, rng: &mut JavaRandom) {
        while let Some((id, _)) = self.pending.front() {
            if *id >= before {
                break;
            }
            let (_, t) = self.pending.pop_front().unwrap();
            self.draws += crate::rng::resolve_qd(&t, rng, &mut self.done);
        }
    }
    pub fn has_pending(&self) -> bool {
        !self.pending.is_empty()
    }
    /// Remaining text, placeholders unresolved.
    pub fn into_raw(mut self) -> Vec<u8> {
        for (_, t) in self.pending.drain(..) {
            self.done.extend_from_slice(&t);
        }
        self.done
    }
}

pub enum BackMsg {
    Batch(Vec<FilledRegion>),
    /// make the QD draws of the regions with id < `before`, then hand the generator to the front end
    Resolve { before: usize, rng: Option<JavaRandom>, reply: std::sync::mpsc::SyncSender<JavaRandom> },
    /// end of the contig (with the generator if the front end holds it)
    End { rng: Option<JavaRandom> },
}

/// Runs the back end of one contig over the messages the front end sends.
pub fn run(cfg: BackendCfg, rx: Receiver<BackMsg>) -> BackendOut {
    let mut st = Stats::default();
    if !cfg.erc {
        let mut ecfg = hc_genotype::EngineConfig::new(&cfg.contig, cfg.contig_len);
        ecfg.n_samples = cfg.opts.n_samples;
        ecfg.defer_qd_jitter = true;
        ecfg.contamination = cfg.opts.contamination;
        ecfg.as_annotations = cfg.opts.as_annotations;
        let mut engine = hc_genotype::GenotypingEngine::new(ecfg);
        let mut text = QdText::default();
        let mut held: Option<JavaRandom> = None;
        let (domain, idx) = (cfg.domain.clone(), cfg.domain_idx);
        let mut end_rng = None;
        staged(&cfg, rx, |m| match m {
            Staged::Done(c) => {
                hold_contam_rng(&cfg, &mut engine);
                crate::cpuprof::measure(6, || genotype_vcf_batch(c, &cfg, &mut engine, &mut text, &mut st))
            }
            Staged::Fwd(BackMsg::Resolve { before, rng, reply }) => {
                let mut r = rng.or(held.take()).unwrap_or_else(|| domain.take(idx));
                text.resolve(before, &mut r);
                reply.send(r).unwrap();
            }
            Staged::Fwd(BackMsg::End { rng }) => end_rng = rng.or(held.take()),
            Staged::Fwd(BackMsg::Batch(_)) => unreachable!(),
        });
        release_contam_rng(&cfg, &mut engine);
        st.qd_jitter = engine.n_qd_jitter;
        // the generator is done with this contig once every placeholder is drawn
        let resolved = match end_rng {
            Some(mut r) => {
                text.resolve(usize::MAX, &mut r);
                domain.finish(idx, r);
                true
            }
            None if !text.has_pending() => {
                domain.pass(idx);
                true
            }
            None => false,
        };
        st.qd_drawn = text.draws;
        return BackendOut { vcf_text: if resolved { text.done } else { text.into_raw() }, vcf_resolved: resolved, gvcf: None, raw: None, stats: st };
    }
    // GVCF: banding + output on a writer thread
    let (part_path, level) = cfg.gvcf_part.clone().expect("GVCF part");
    let out = PartWriter::create(part_path, level).expect("create GVCF part");
    let raw = cfg.raw_part.as_ref().map(|p| BufWriter::with_capacity(1 << 20, File::create(p).expect("create raw part")));
    let (tx, wrx) = std::sync::mpsc::sync_channel::<Vec<Vec<hc_refconf::Item>>>(4);
    let contig = cfg.contig.clone();
    let bands = cfg.opts.gq_bands.clone();
    let writer = std::thread::Builder::new()
        .name("gvcf-writer".into())
        .spawn(move || crate::cpuprof::measure(8, || {
            let mut gw = hc_refconf::GvcfWriter::new(&contig, out, raw, &bands);
            let mut busy = 0f64;
            for batch in wrx {
                let t = Instant::now();
                for it in batch {
                    gw.add_region(it).unwrap();
                }
                busy += t.elapsed().as_secs_f64();
            }
            let t = Instant::now();
            gw.finish().unwrap();
            let n_out = gw.n_out;
            let (out, raw) = gw.into_inner();
            if let Some(mut r) = raw {
                r.flush().unwrap();
            }
            let (part, deflate) = out.finish().unwrap();
            busy += t.elapsed().as_secs_f64();
            (part, n_out, busy, deflate)
        }))
        .unwrap();
    let mut ecfg = hc_genotype::EngineConfig::new_erc(&cfg.contig, cfg.contig_len);
    ecfg.contamination = cfg.opts.contamination;
    ecfg.as_annotations = cfg.opts.as_annotations;
    let mut engine = hc_genotype::GenotypingEngine::new(ecfg);
    staged(&cfg, rx, |m| {
        let Staged::Done(c) = m else { return };
        hold_contam_rng(&cfg, &mut engine);
        let items = crate::cpuprof::measure(7, || genotype_erc_batch(c, &cfg, &mut engine, &mut st));
        let t = Instant::now();
        tx.send(items).unwrap();
        st.write += t.elapsed().as_secs_f64();
    });
    release_contam_rng(&cfg, &mut engine);
    drop(tx);
    let (part, n_out, busy, deflate) = writer.join().unwrap();
    st.out_lines += n_out;
    st.writer_busy += busy;
    st.deflate += deflate;
    let _ = cfg.tid;
    BackendOut { vcf_text: Vec::new(), vcf_resolved: true, gvcf: Some(part), raw: cfg.raw_part.clone(), stats: st }
}
