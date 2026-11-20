//! `variant-rs mutect2`: port of GATK 4.7 Mutect2 (tumor-normal, default arguments), in stages:
//! S1 front end (reads, filters, downsampling, isActive, assembly regions, callable sites) — this file and
//! `front.rs`; S2 region processing (`region.rs`); S3 somatic genotyping, annotation and VCF output
//! (`call.rs`, `somatic.rs`, `germline.rs`).

pub mod bamout;
pub mod call;
pub mod cluster;
pub mod cm_linalg;
pub mod contam;
pub mod orientation;
pub mod pon;
pub mod fmc;
pub mod fmc_math;
pub mod front;
pub mod gathers;
pub mod germline;
pub mod math;
pub mod region;
pub mod somatic;

use crate::active::{Region, Tables};
use crate::frontend::{FilledRegion, M2Quals, RefSeq};
use std::sync::Arc;
use crate::input::Input;
use crate::intervals::{self, Dict, Interval};
use front::{M2FilterCounts, M2FrontResult, M2Input, M2Sink};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

const USAGE: &str = "usage: variant-rs mutect2 -R ref.fa -I tumor.bam -I normal.bam -normal SAMPLE [-L interval|file]... \
[--germline-resource gnomad.vcf] [-pon pon.vcf] [--f1r2-tar-gz f1r2.tar.gz] -O out.vcf [--pairhmm auto|cuda|cpu] \
[--assembly-region-out regions.igv] [--dump-activity file] [--max-reads-per-alignment-start N (50)]  (env VARIANT_RS_THREADS: region workers)";

struct Args {
    reference: String,
    inputs: Vec<String>,
    normals: Vec<String>,
    intervals: Vec<String>,
    germline: Option<String>,
    f1r2: Option<String>,
    pon: Option<String>,
    pairhmm: String,
    out: Option<String>,
    region_out: Option<String>,
    dump: Option<String>,
    dump_s2: Option<String>,
    max_reads_per_start: usize,
    genotype_germline: bool,
    genotype_pon: bool,
    af_not_in_resource: Option<f64>,
    mitochondria: bool,
    annotations: Vec<String>,
    read_filters: Vec<String>,
    max_mnp_distance: i32,
    bam_output: Option<String>,
}

/// Mutect2 options that change the calls (M2ArgumentCollection).
#[derive(Clone, Copy, Debug)]
pub struct M2Opts {
    /// `--genotype-germline-sites`
    pub genotype_germline: bool,
    /// `--genotype-pon-sites`
    pub genotype_pon: bool,
    /// `getInitialLogOdds()` (natural log)
    pub initial_log_odds: f64,
    /// `getEmissionLogOdds()` (natural log)
    pub emission_log_odds: f64,
    /// `--max-population-af`
    pub max_population_af: f64,
    /// the `OriginalAlignment` annotation (OCM; added by `--mitochondria-mode` unless `--annotation` is given)
    pub original_alignment: bool,
    /// `--max-mnp-distance`
    pub max_mnp_distance: i32,
}

fn parse_args(argv: &[String]) -> Args {
    let mut a = Args { reference: String::new(), inputs: vec![], normals: vec![], intervals: vec![], germline: None, f1r2: None, pon: None, pairhmm: "auto".into(), out: None, region_out: None, dump: None, dump_s2: None, max_reads_per_start: 50, genotype_germline: false, genotype_pon: false, af_not_in_resource: None, mitochondria: false, annotations: vec![], read_filters: vec![], max_mnp_distance: 1, bam_output: None };
    let mut i = 0;
    while i < argv.len() {
        let val = |i: usize| -> String { argv.get(i + 1).cloned().unwrap_or_else(|| panic!("{} needs a value\n{USAGE}", argv[i])) };
        match argv[i].as_str() {
            "-R" | "--reference" => a.reference = val(i),
            "-I" | "--input" => a.inputs.push(val(i)),
            "-normal" | "--normal-sample" => a.normals.push(val(i)),
            "-L" | "--intervals" => a.intervals.push(val(i)),
            "--germline-resource" => a.germline = Some(val(i)),
            "--f1r2-tar-gz" => a.f1r2 = Some(val(i)),
            "-pon" | "--panel-of-normals" => a.pon = Some(val(i)),
            "--pairhmm" => a.pairhmm = val(i),
            "-O" | "--output" => a.out = Some(val(i)),
            "--assembly-region-out" => a.region_out = Some(val(i)),
            "--dump-activity" => a.dump = Some(val(i)),
            "--dump-s2" => a.dump_s2 = Some(val(i)),
            "--max-reads-per-alignment-start" => a.max_reads_per_start = val(i).parse().expect("--max-reads-per-alignment-start"),
            "--genotype-germline-sites" => {
                a.genotype_germline = true;
                i += 1;
                continue;
            }
            "--genotype-pon-sites" => {
                a.genotype_pon = true;
                i += 1;
                continue;
            }
            "--mitochondria-mode" => {
                a.mitochondria = true;
                i += 1;
                continue;
            }
            "-A" | "--annotation" => a.annotations.push(val(i)),
            "-bamout" | "--bam-output" => a.bam_output = Some(val(i)),
            "-RF" | "--read-filter" => a.read_filters.push(val(i)),
            "-mnp-dist" | "--max-mnp-distance" => a.max_mnp_distance = val(i).parse().expect("--max-mnp-distance"),
            "-default-af" | "--af-of-alleles-not-in-resource" => a.af_not_in_resource = Some(val(i).parse().expect("--af-of-alleles-not-in-resource")),
            x => panic!("unknown argument {x}\n{USAGE}"),
        }
        i += 2;
    }
    assert!(!a.reference.is_empty() && !a.inputs.is_empty(), "{USAGE}");
    a
}

/// `@RG` ID -> SM of a SAM header text.
fn read_groups(text: &str) -> Vec<(String, String)> {
    let mut v = Vec::new();
    for l in text.lines().filter(|l| l.starts_with("@RG")) {
        let (mut id, mut sm) = (None, None);
        for f in l.split('\t').skip(1) {
            if let Some(x) = f.strip_prefix("ID:") {
                id = Some(x.to_string());
            } else if let Some(x) = f.strip_prefix("SM:") {
                sm = Some(x.to_string());
            }
        }
        if let (Some(id), Some(sm)) = (id, sm) {
            v.push((id, sm));
        }
    }
    v
}

/// Java `Double.toString` for a non-negative integer-valued double.
fn java_double_int(v: u64) -> String {
    if v < 10_000_000 {
        return format!("{v}.0");
    }
    let s = format!("{:e}", v as f64); // shortest digits, e.g. 1.2345e7
    let (m, e) = s.split_once('e').unwrap();
    let m = if m.contains('.') { m.to_string() } else { format!("{m}.0") };
    format!("{m}E{e}")
}

/// One contig being called: its index in output order, name and reference.
struct ContigData {
    idx: usize,
    name: String,
    refseq: Arc<RefSeq>,
}

/// One active region handed to the workers: (contig, sequence number, region, base qualities as of `callRegion`).
type Job = (Arc<ContigData>, usize, FilledRegion, Vec<Option<Arc<M2Quals>>>);
/// A worker's result for one region: (contig index, sequence number, exit, --dump-s2 line, VCF lines).
type Done = (usize, usize, &'static str, Option<String>, Vec<String>, Vec<bamout::BRec>);

/// Sink of the dry pass (no regions are made).
struct NoSink;
impl M2Sink for NoSink {
    fn region(&mut self, _: FilledRegion, _: Vec<Option<Arc<M2Quals>>>) {}
}

/// Hands each active region to the stage S2 / S3 worker pool as GATK's walker reaches it (regions are
/// independent once their read qualities are snapshotted; results are put back in region order).
struct RegionSink {
    tx: std::sync::mpsc::SyncSender<Job>,
    next: usize,
    contig: Arc<ContigData>,
}

impl M2Sink for RegionSink {
    fn region(&mut self, fr: FilledRegion, quals: Vec<Option<Arc<M2Quals>>>) {
        if !fr.active || fr.reads.is_empty() {
            return;
        }
        // VARIANT_RS_TRACE_REGIONS=1: per-region hand-off times (profiling)
        if std::env::var_os("VARIANT_RS_TRACE_REGIONS").is_some() {
            eprintln!("region {}-{} reads={} t={:.1}", fr.padded_start, fr.padded_end, fr.reads.len(), crate::mutect::trace_t());
        }
        self.tx.send((self.contig.clone(), self.next, fr, quals)).expect("region workers stopped");
        self.next += 1;
    }
}

/// Stage S2 (assembly, likelihoods) and S3 (somatic genotyping) of one region.
fn process_job(
    job: (usize, FilledRegion, Vec<Option<Arc<M2Quals>>>),
    ctx: &region::RegionCtx,
    call_ctx: &call::CallCtx,
    samples: &[String],
    asm: &mut hc_assembly::Assembler,
    sw: &mut hc_assembly::sw::SwAligner,
    lp: &hc_likelihood::Params,
    dump: bool,
) -> Done {
    let (id, fr, quals) = job;
    let t0 = std::time::Instant::now();
    let o = region::process_region(&fr, &quals, ctx, asm, sw, lp);
    if std::env::var_os("VARIANT_RS_TRACE_REGIONS").is_some() {
        eprintln!("work {}-{} reads={} s={:.2}", fr.padded_start, fr.padded_end, fr.reads.len(), t0.elapsed().as_secs_f64());
    }
    let d = if dump { Some(region::dump_json(&fr, &o, samples)) } else { None };
    let mut lines = Vec::new();
    let mut brecs = Vec::new();
    if let (Some((haps, per)), Some((vs, ve, _, _))) = (&o.likelihoods, o.variant_region) {
        let refseq = ctx.refseq;
        let rc = call::RegionCall {
            haplotypes: haps,
            samples: per,
            ref_padded: refseq.window(o.ref_padded.0, o.ref_padded.1),
            ref_padded_start: o.ref_padded.0,
            span: (vs, ve),
            ref_window: refseq.window(fr.padded_start, fr.padded_end),
            ref_window_start: fr.padded_start,
        };
        let t = std::time::Instant::now();
        let mut info = call::BamoutInfo::default();
        let want_bamout = call_ctx.bamout_tid.is_some();
        for c in call::call_mutations_ext(call_ctx, &rc, want_bamout.then_some(&mut info)) {
            lines.push(call::to_line(call_ctx.contig, &c));
        }
        if let Some(tid) = call_ctx.bamout_tid {
            brecs = bamout::region_records(call_ctx.contig, tid, haps, o.ref_padded.0, (vs, ve), per, &info);
        }
        region::stage_add(6, t);
    }
    (0, id, o.exit, d, lines, brecs)
}

/// Mutect2's assembler (adaptive pruning); `--mitochondria-mode` recovers all dangling branches and prunes with a
/// log-odds threshold of ln(1e-4) (`DEFAULT_MITO_PRUNING_LOG_ODDS_THRESHOLD`).
fn new_assembler(mito: bool) -> hc_assembly::Assembler {
    let mut ap = hc_assembly::rtgraph::AdaptivePruning::default();
    if mito {
        ap.log_odds_threshold = -4.0 * std::f64::consts::LN_10;
    }
    hc_assembly::Assembler::new(hc_assembly::AssemblerConfig {
        min_prune_factor: 0,
        adaptive_pruning: Some(ap),
        recover_all_dangling_branches: mito,
        ..Default::default()
    })
}

pub fn run(argv: &[String]) {
    let t0 = std::time::Instant::now();
    let args = parse_args(argv);
    // inputs and samples (sample order: first appearance over the inputs)
    let mut samples: Vec<String> = Vec::new();
    let mut inputs: Vec<M2Input> = Vec::new();
    let mut names_per_input: Vec<Vec<String>> = Vec::new();
    let mut contig_lines: Vec<(String, i64, Option<String>)> = Vec::new();
    let mut headers: Vec<crate::bam::Header> = Vec::new();
    for p in &args.inputs {
        let (input, header) = Input::open(p, &args.reference);
        if args.bam_output.is_some() {
            headers.push(crate::bam::Header { text: header.text.clone(), ref_names: header.ref_names.clone(), ref_lens: header.ref_lens.clone() });
        }
        if contig_lines.is_empty() {
            let asm = header.sq_assembly();
            contig_lines = header.ref_names.iter().zip(&header.ref_lens).map(|(n, &l)| (n.clone(), l, asm.get(n).cloned())).collect();
        }
        let mut rg_sample = HashMap::new();
        for (id, sm) in read_groups(&header.text) {
            let k = samples.iter().position(|s| *s == sm).unwrap_or_else(|| {
                samples.push(sm.clone());
                samples.len() - 1
            });
            rg_sample.insert(id.into_bytes(), k as u16);
        }
        names_per_input.push(header.ref_names.clone());
        inputs.push(M2Input {
            input,
            rg_sample,
            ref_lens: header.ref_lens.clone(),
            ref_names: header.ref_names.clone(),
            mate_same_contig: args.read_filters.iter().any(|f| f == "MateOnSameContigOrNoMappedMateReadFilter"),
            mate_mapped: args.read_filters.iter().any(|f| f == "MateUnmappedAndUnmappedReadFilter"),
            keep_records: args.bam_output.is_some(),
            merged: None,
        });
    }
    // --bam-output with several inputs: GATK reads them through htsjdk's header merger
    let merged = (args.bam_output.is_some() && headers.len() > 1).then(|| {
        let hs: Vec<&crate::bam::Header> = headers.iter().collect();
        Arc::new(crate::markdup::qname::merge_header_list(&hs, "coordinate"))
    });
    if let Some(m) = &merged {
        for (i, inp) in inputs.iter_mut().enumerate() {
            inp.merged = Some((m.clone(), i));
        }
    }
    for n in &args.normals {
        assert!(samples.contains(n), "normal sample {n} is not in the BAM headers: {samples:?}");
    }
    let is_tumor: Vec<bool> = samples.iter().map(|s| !args.normals.contains(s)).collect();
    let n_tumor = is_tumor.iter().filter(|&&t| t).count();
    assert_eq!(n_tumor, 1, "exactly one tumor sample is supported for now (samples {samples:?})");
    eprintln!("samples: {}", samples.iter().zip(&is_tumor).map(|(s, t)| format!("{s} ({})", if *t { "tumor" } else { "normal" })).collect::<Vec<_>>().join(", "));

    let ref_dict = Dict::from_reference(&args.reference);
    let ivs: Vec<Interval> = if args.intervals.is_empty() {
        names_per_input[0]
            .iter()
            .map(|n| {
                let t = *ref_dict.index.get(n).unwrap_or_else(|| panic!("BAM contig {n} not in the reference"));
                Interval { tid: t, start: 1, end: ref_dict.lens[t] }
            })
            .collect()
    } else {
        let mut v = Vec::new();
        for l in &args.intervals {
            v.extend(intervals::parse_l_arg(l, &ref_dict).unwrap_or_else(|e| panic!("-L {l}: {e}")));
        }
        v
    };
    let ivs = intervals::sort_merge(ivs);
    let fasta = crate::reference::Fasta::open(&args.reference);
    let tables = Tables::load();
    let math = math::M2Math::new();
    let rng = hc_downsample::JavaRandom::gatk();
    let mut dump = args.dump.as_ref().map(|p| BufWriter::with_capacity(1 << 20, File::create(p).unwrap()));
    let mut dump_s2 = args.dump_s2.as_ref().map(|p| BufWriter::with_capacity(1 << 22, File::create(p).unwrap()));
    // GATK samplesList: sorted sample names
    let mut sorted: Vec<usize> = (0..samples.len()).collect();
    sorted.sort_by(|&a, &b| samples[a].cmp(&samples[b]));
    let mut sample_rank = vec![0usize; samples.len()];
    for (rank, &s) in sorted.iter().enumerate() {
        sample_rank[s] = rank;
    }
    let mut exits: std::collections::BTreeMap<&'static str, u64> = Default::default();
    let use_gpu = args.pairhmm != "cpu" && hc_cudahmm::init(false);
    assert!(args.pairhmm != "cuda" || use_gpu, "--pairhmm cuda requested but no CUDA device (use --pairhmm cpu or auto)");
    // VARIANT_RS_THREADS: region workers (default: all cores)
    let n_workers: usize = std::env::var("VARIANT_RS_THREADS").ok().and_then(|v| v.parse().ok()).unwrap_or_else(|| std::thread::available_parallelism().map_or(4, |n| n.get()));
    let mut f1r2 = args.f1r2.as_ref().map(|_| orientation::F1R2Counts::new(&samples[is_tumor.iter().position(|&t| t).unwrap()]));
    let is_normal_sorted: Vec<bool> = sorted.iter().map(|&s| !is_tumor[s]).collect();
    let germline = args.germline.as_ref().map(|p| {
        let contigs: std::collections::HashSet<String> = ivs.iter().map(|iv| ref_dict.names[iv.tid].clone()).collect();
        germline::Germline::load(p, &contigs)
    });
    let pon = args.pon.as_ref().map(|p| {
        let contigs: std::collections::HashSet<String> = ivs.iter().map(|iv| ref_dict.names[iv.tid].clone()).collect();
        germline::Pon::load(p, &contigs)
    });
    // M2ArgumentCollection.getDefaultAlleleFrequency (a negative --af-of-alleles-not-in-resource means unset)
    // --mitochondria-mode (ModeArgumentUtils: only arguments not given on the command line): default AF 4e-3,
    // emission and initial LOD 0, recover all dangling branches, pruning LOD ln(1e-4), annotation OriginalAlignment
    let default_af = match args.af_not_in_resource {
        Some(x) => if x >= 0.0 { x } else if args.normals.is_empty() { 5e-8 } else { 1e-6 },
        None if args.mitochondria => 4e-3,
        None => if args.normals.is_empty() { 5e-8 } else { 1e-6 },
    };
    for a in &args.annotations {
        assert!(a == "StrandBiasBySample" || a == "OriginalAlignment", "--annotation {a} is not supported");
    }
    for f in &args.read_filters {
        assert!(f == "MateOnSameContigOrNoMappedMateReadFilter" || f == "MateUnmappedAndUnmappedReadFilter", "--read-filter {f} is not supported");
    }
    let opts = M2Opts {
        genotype_germline: args.genotype_germline,
        genotype_pon: args.genotype_pon,
        initial_log_odds: if args.mitochondria { 0.0 } else { 2.0 * std::f64::consts::LN_10 },
        emission_log_odds: if args.mitochondria { 0.0 } else { 3.0 * std::f64::consts::LN_10 },
        max_population_af: 0.01,
        original_alignment: args.annotations.iter().any(|a| a == "OriginalAlignment") || (args.mitochondria && args.annotations.is_empty()),
        max_mnp_distance: args.max_mnp_distance,
    };
    let mut vcf_lines: Vec<String> = Vec::new();

    let mut all_regions: Vec<(String, Vec<Region>)> = Vec::new();
    let mut fc = M2FilterCounts::default();
    let (mut callable, mut loci, mut active_loci, mut reads_in, mut reads_kept, mut groups, mut assigned) = (0u64, 0u64, 0u64, 0u64, 0u64, 0u64, 0u64);
    // Work units ("contigs" below): each contig's intervals, split into chunks at gaps of at least VARIANT_RS_SPLIT_GAP
    // (default 5,000 bp; 0 = no splitting), except where a read crosses the gap: no read belongs to two chunks (a
    // chunk reads the records overlapping its intervals padded by REGION_PADDING). The activity profile is flushed at every interval gap anyway, so a
    // chunk behaves exactly as that stretch of the contig did. Chunks run concurrently (VARIANT_RS_CONTIG_THREADS front ends;
    // 1 with --dump-activity). The MutectDownsampler RNG runs through the chunks in order and its draws only depend on
    // the drawing groups' sizes: a dry pass counts them per chunk, a replay gives each chunk's starting RNG state.
    let n_front_wanted: usize = if dump.is_some() {
        1
    } else {
        std::env::var("VARIANT_RS_CONTIG_THREADS").ok().and_then(|v| v.parse().ok()).unwrap_or((n_workers * 3 / 8).max(1))
    };
    let split_gap: i64 = if n_front_wanted == 1 { 0 } else { std::env::var("VARIANT_RS_SPLIT_GAP").ok().and_then(|v| v.parse().ok()).unwrap_or(5000) };
    // whether a record of any input overlaps both a padded end at `a` and a padded start at `b` (> a + 1): such a
    // record covers a + 1, so the index query of that one locus finds it
    let crosses = |tids: &[usize], a: i64, b: i64| -> bool {
        inputs.iter().zip(tids).any(|(inp, &tid)| {
            let mut s = inp.input.stream(tid, &[(a + 1, a + 1)], false, Arc::new(std::sync::atomic::AtomicU64::new(0)));
            let mut body = Vec::new();
            while s.next_record(&mut body) {
                let (rid, rstart, _) = crate::bam::peek_core(&body);
                if rid == tid as i32 && rstart <= a && crate::bam::query_end(&body) >= b {
                    return true;
                }
            }
            false
        })
    };
    let pad = crate::frontend::REGION_PADDING;
    let by_contig = intervals::group_by_contig(&ivs);
    let total_loci: i64 = by_contig.iter().flat_map(|(_, l)| l.iter()).map(|&(a, b)| b - a + 1).sum();
    let target = (total_loci / (n_front_wanted as i64 * 4)).max(1);
    let mut contigs: Vec<(usize, Vec<(i64, i64)>, Vec<usize>)> = Vec::new();
    for (rt, list) in by_contig {
        let cname = &ref_dict.names[rt];
        let tids: Vec<usize> = names_per_input
            .iter()
            .map(|names| names.iter().position(|n| n == cname).unwrap_or_else(|| panic!("contig {cname} is not in every BAM header")))
            .collect();
        let mut cur: Vec<(i64, i64)> = Vec::new();
        let mut loci = 0i64;
        for iv in list {
            if let Some(last) = cur.last() {
                if split_gap > 0 && loci >= target && iv.0 - last.1 - 1 >= split_gap && !crosses(&tids, last.1 + pad, iv.0 - pad) {
                    contigs.push((rt, std::mem::take(&mut cur), tids.clone()));
                    loci = 0;
                }
            }
            loci += iv.1 - iv.0 + 1;
            cur.push(iv);
        }
        contigs.push((rt, cur, tids));
    }
    // the whole contig's reference (palindrome clipping looks at the adaptor boundary, which may be far from the read),
    // shared by the chunks of a contig that are in flight
    let ref_cache: std::sync::Mutex<HashMap<usize, std::sync::Weak<RefSeq>>> = std::sync::Mutex::new(HashMap::new());
    let fetch = |rt: usize, _list: &[(i64, i64)]| -> Arc<RefSeq> {
        let mut c = ref_cache.lock().unwrap();
        if let Some(r) = c.get(&rt).and_then(|w| w.upgrade()) {
            return r;
        }
        let (seq, contig_len) = fasta.fetch(&ref_dict.names[rt], 1, ref_dict.lens[rt]);
        let r = Arc::new(RefSeq { from: 1, seq, contig_len });
        c.insert(rt, Arc::downgrade(&r));
        r
    };
    let n_front = n_front_wanted.clamp(1, contigs.len().max(1));
    // two-stage front ends (read stage on its own thread) while the machine has room for them (VARIANT_RS_TWO_STAGE=0/1)
    let two_stage = match std::env::var("VARIANT_RS_TWO_STAGE").ok().as_deref() {
        Some("0") => false,
        Some(_) => true,
        None => n_front * 4 <= n_workers,
    };
    // BAM decoding threads per input and front end (VARIANT_RS_DECODE_THREADS)
    let decode_threads: usize = std::env::var("VARIANT_RS_DECODE_THREADS").ok().and_then(|v| v.parse().ok()).unwrap_or((n_workers / (n_front * inputs.len()).max(1)).clamp(1, 8));
    let tumor_name = samples[is_tumor.iter().position(|&t| t).unwrap()].clone();
    // dry-pass results: per contig the drawing groups' sizes (None = not counted yet)
    let dry_sizes: std::sync::Mutex<Vec<Option<Vec<usize>>>> = std::sync::Mutex::new(vec![None; contigs.len()]);
    let dry_cv = std::sync::Condvar::new();
    let need_dry = n_front > 1 && args.max_reads_per_start > 0;
    if n_front > 1 {
        eprintln!("front ends: {n_front} over {} chunks", contigs.len());
    }
    // RNG state at the start of contig k: wait for the dry counts of contigs 0..k, replay their draws
    let start_state = |k: usize| -> hc_downsample::JavaRandomState {
        let mut r = rng.clone();
        if !need_dry {
            return r.state();
        }
        let mut g = dry_sizes.lock().unwrap();
        while g[..k].iter().any(|x| x.is_none()) {
            g = dry_cv.wait(g).unwrap();
        }
        for sz in g[..k].iter() {
            for &n in sz.as_ref().unwrap() {
                for t in args.max_reads_per_start + 1..=n {
                    r.next_int(t as i32);
                }
            }
        }
        r.state()
    };
    let lp = hc_likelihood::Params { max_mnp_distance: args.max_mnp_distance as _, mutect: true, ..Default::default() };
    let mito = args.mitochondria;
    let want_bam = args.bam_output.is_some();
    let bam_names: Vec<String> = headers.first().map_or(Vec::new(), |h| h.ref_names.clone());
    let want_dump = dump_s2.is_some();
    type FrontOut = (M2FrontResult, Option<orientation::F1R2Counts>);
    let outs: std::sync::Mutex<Vec<Option<FrontOut>>> = std::sync::Mutex::new((0..contigs.len()).map(|_| None).collect());
    let dump_slot = std::sync::Mutex::new(dump.take());
    let chained = std::sync::Mutex::new(rng.clone());
    let mut done: Vec<Done> = std::thread::scope(|scope| {
        let (tx, rx) = std::sync::mpsc::sync_channel::<Job>(2 * n_workers);
        let rx = Arc::new(std::sync::Mutex::new(rx));
        let (dtx, drx) = std::sync::mpsc::channel::<Done>();
        for _ in 0..n_workers {
            let (rx, dtx, lp, samples) = (rx.clone(), dtx.clone(), &lp, &samples);
            let (sample_rank, is_normal_sorted, germline, pon) = (&sample_rank, &is_normal_sorted, germline.as_ref(), pon.as_ref());
            let bam_names = &bam_names;
            scope.spawn(move || {
                let (mut asm, mut sw) = (new_assembler(mito), hc_assembly::sw::SwAligner::new());
                loop {
                    let job = rx.lock().unwrap().recv();
                    let Ok((cd, id, fr, q)) = job else { break };
                    let ctx = region::RegionCtx { refseq: &cd.refseq, sample_rank, n_samples: samples.len(), use_gpu };
                    let call_ctx = call::CallCtx { contig: &cd.name, contig_len: cd.refseq.contig_len, is_normal: is_normal_sorted, germline, pon, default_af, opts, bamout_tid: want_bam.then(|| bam_names.iter().position(|n| *n == cd.name).unwrap() as i32) };
                    let mut d = process_job((id, fr, q), &ctx, &call_ctx, samples, &mut asm, &mut sw, lp, want_dump);
                    d.0 = cd.idx;
                    dtx.send(d).unwrap();
                }
            });
        }
        drop(dtx);
        // dry pass (contigs in order; the last contig's counts are never needed)
        if need_dry {
            let next_dry = Arc::new(std::sync::atomic::AtomicUsize::new(0));
            for _ in 0..(n_workers / 4).clamp(1, contigs.len()) {
                let next_dry = next_dry.clone();
                let (contigs, dry_sizes, dry_cv, fetch) = (&contigs, &dry_sizes, &dry_cv, &fetch);
                let (inputs, is_tumor, tables, math, pon, ref_dict, max_reads) = (&inputs, &is_tumor, &tables, &math, pon.as_ref(), &ref_dict, args.max_reads_per_start);
                let germline = germline.as_ref();
                scope.spawn(move || loop {
                    let k = next_dry.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if k + 1 >= contigs.len() {
                        break;
                    }
                    let (rt, list, tids) = &contigs[k];
                    let refseq = fetch(*rt, list);
                    let mut dummy = hc_downsample::JavaRandom::gatk();
                    let (res, _) = front::run(inputs, tids, is_tumor, &ref_dict.names[*rt], list, &refseq, tables, math, max_reads, &mut dummy, None, NoSink, None, pon, germline, opts, true, decode_threads, false);
                    dry_sizes.lock().unwrap()[k] = Some(res.group_sizes);
                    dry_cv.notify_all();
                });
            }
        }
        let next = std::sync::atomic::AtomicUsize::new(0);
        std::thread::scope(|fs| {
            for _ in 0..n_front {
                let tx = tx.clone();
                let (next, outs, contigs, start_state, dump_slot, chained, fetch, tumor_name) = (&next, &outs, &contigs, &start_state, &dump_slot, &chained, &fetch, &tumor_name);
                let (inputs, is_tumor, tables, math, pon, ref_dict) = (&inputs, &is_tumor, &tables, &math, pon.as_ref(), &ref_dict);
                let germline = germline.as_ref();
                let (want_f1r2, max_reads) = (args.f1r2.is_some(), args.max_reads_per_start);
                fs.spawn(move || loop {
                    let k = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    let Some((rt, list, tids)) = contigs.get(k) else { break };
                    let t_contig = std::time::Instant::now();
                    let cd = Arc::new(ContigData { idx: k, name: ref_dict.names[*rt].clone(), refseq: fetch(*rt, list) });
                    let mut f1r2 = want_f1r2.then(|| orientation::F1R2Counts::new(tumor_name));
                    // one front end: the RNG continues from the previous contig; several: replayed start states
                    let mut chained_guard = if n_front == 1 { Some(chained.lock().unwrap()) } else { None };
                    let mut own = (n_front > 1).then(|| hc_downsample::JavaRandom::from_state(start_state(k)));
                    let rng: &mut hc_downsample::JavaRandom = match (own.as_mut(), chained_guard.as_mut()) {
                        (Some(r), _) => r,
                        (None, Some(g)) => g,
                        _ => unreachable!(),
                    };
                    // the activity dump only exists with a single front end
                    let mut dump_guard = (n_front == 1).then(|| dump_slot.lock().unwrap());
                    let (res, sink) = front::run(
                        inputs, tids, is_tumor, &cd.name, list, &cd.refseq, tables, math, max_reads, rng, dump_guard.as_mut().and_then(|g| g.as_mut()),
                        RegionSink { tx: tx.clone(), next: 0, contig: cd.clone() }, f1r2.as_mut(), pon, germline, opts, false, decode_threads, two_stage,
                    );
                    drop(dump_guard);
                    drop(sink);
                    if std::env::var_os("VARIANT_RS_TRACE_STAGES").is_some() {
                        eprintln!("contig {} front {:.1}s (ended at {:.1}s)", cd.name, t_contig.elapsed().as_secs_f64(), crate::mutect::trace_t());
                    }
                    outs.lock().unwrap()[k] = Some((res, f1r2));
                });
            }
        });
        drop(tx);
        drx.iter().collect::<Vec<Done>>()
    });
    if let Some(d) = dump_slot.into_inner().unwrap() {
        dump = Some(d);
    }
    done.sort_by_key(|d| (d.0, d.1));
    let mut bam_regions: Vec<Vec<bamout::BRec>> = Vec::new();
    for (_, _, exit, d, lines, brecs) in done {
        bam_regions.push(brecs);
        *exits.entry(exit).or_default() += 1;
        if let (Some(w), Some(d)) = (dump_s2.as_mut(), d) {
            writeln!(w, "{d}").unwrap();
        }
        vcf_lines.extend(lines);
    }
    if let Some(path) = &args.bam_output {
        // the reads' header: the one input's, or htsjdk's merge of them (ReadsPathDataSource)
        let src = match &merged {
            Some(m) => m.text.clone(),
            None => headers[0].text.clone(),
        };
        bamout::write(path, &bamout::header_text(&src), &headers[0].ref_names, &headers[0].ref_lens, bam_regions);
    }
    for (k, o) in outs.into_inner().unwrap().into_iter().enumerate() {
        let (res, c) = o.expect("contig not processed");
        if let (Some(all), Some(c)) = (f1r2.as_mut(), c) {
            all.merge(c);
        }
        fc.add(&res.filters);
        callable += res.callable;
        loci += res.n_loci;
        active_loci += res.n_active_loci;
        reads_in += res.reads_in;
        reads_kept += res.reads_kept;
        groups += res.drawing_groups;
        assigned += res.assigned;
        all_regions.push((ref_dict.names[contigs[k].0].clone(), res.regions));
    }
    if let Some(mut d) = dump {
        d.flush().unwrap();
    }
    if let (Some(c), Some(p)) = (&f1r2, &args.f1r2) {
        c.write_tar_gz(p);
    }
    if let Some(mut d) = dump_s2 {
        d.flush().unwrap();
    }
    eprintln!("active regions: {exits:?}");
    if front::ft_on() {
        for (n, c) in front::FT_NAMES.iter().zip(&front::FT) {
            eprintln!("front time {n}: {:.1} s", c.load(std::sync::atomic::Ordering::Relaxed) as f64 * 1e-9);
        }
        eprintln!("front pairs skipped {} adjusted {}", front::FTC[0].load(std::sync::atomic::Ordering::Relaxed), front::FTC[1].load(std::sync::atomic::Ordering::Relaxed));
    }
    if std::env::var_os("VARIANT_RS_TRACE_STAGES").is_some() {
        for (n, c) in region::STAGE_NAMES.iter().zip(&region::STAGE_NS) {
            eprintln!("stage {n}: {:.1} s", c.load(std::sync::atomic::Ordering::Relaxed) as f64 * 1e-9);
        }
    }
    if let Some(p) = &args.region_out {
        let mut w = BufWriter::new(File::create(p).unwrap());
        writeln!(w, "#track graphType=line").unwrap();
        writeln!(w, "Chromosome\tStart\tEnd\tFeature\tAssemblyRegions").unwrap();
        for (c, regs) in &all_regions {
            for r in regs {
                writeln!(w, "{}\t{}\t{}\tend-marker\t0.00000", c, r.start - 1, r.start).unwrap();
                writeln!(w, "{}\t{}\t{}\tsize={}\t{}", c, r.start - 1, r.end, r.end - r.start + 1, if r.active { "1.00000" } else { "-1.00000" }).unwrap();
            }
        }
        w.flush().unwrap();
    }
    if let Some(o) = &args.out {
        let sorted_names: Vec<String> = sorted.iter().map(|&k| samples[k].clone()).collect();
        let tumors: Vec<String> = sorted_names.iter().zip(&is_normal_sorted).filter(|x| !*x.1).map(|x| x.0.clone()).collect();
        let cmdline: String = std::env::args().collect::<Vec<_>>().join(" ");
        // `.vcf.gz`: bgzip + .tbi, as GATK (--create-output-variant-index true)
        let level = crate::output::block_compressed(o).then_some(5);
        let header = m2_header(&contig_lines, &sorted_names, &args.normals, &tumors, &cmdline, &crate::now_string());
        let mut w = crate::output::OutFile::create(o, &header, level).unwrap_or_else(|e| panic!("{e}"));
        let mut body = String::with_capacity(vcf_lines.iter().map(|l| l.len() + 1).sum());
        for l in &vcf_lines {
            body.push_str(l);
            body.push('\n');
        }
        w.add_text(body.as_bytes()).unwrap();
        w.finish().unwrap();
        eprintln!("{}: {} records", o, vcf_lines.len());
        // MutectStats table
        std::fs::write(format!("{o}.stats"), format!("statistic\tvalue\ncallable\t{}\n", java_double_int(callable))).unwrap();
    }
    let n_regions: usize = all_regions.iter().map(|(_, r)| r.len()).sum();
    let n_act: usize = all_regions.iter().map(|(_, r)| r.iter().filter(|x| x.active).count()).sum();
    eprintln!(
        "reads: total={} passed={} filtered: mapq={} mq255={} unmapped={} secondary={} duplicate={} qcfail={} chimericOA={} zeroRefLen={} short={} badCigar={} malformed={}; palindrome-clipped={}",
        fc.total, fc.passed, fc.mapq, fc.mapq_unavailable, fc.unmapped, fc.secondary, fc.duplicate, fc.qcfail, fc.chimeric_oa, fc.zero_ref_len, fc.short, fc.bad_cigar, fc.malformed, fc.palindrome_clipped
    );
    eprintln!("downsampling (max {} reads per start): {reads_in} in, {reads_kept} kept, {groups} drawing groups", args.max_reads_per_start);
    eprintln!("loci={loci} active={active_loci} callable={callable} regions={n_regions} (active {n_act}) read->region assignments={assigned} in {:.1}s", t0.elapsed().as_secs_f64());
    if std::env::var_os("VARIANT_RS_ASM_PROFILE").is_some() {
        eprintln!("assembly events:\n{}", hc_assembly::prof::hits_report());
    }
}

const M2_HEADER_LINES: &str = include_str!("../../tables/m2_header_lines.txt");

/// Mutect2's VCF header as htsjdk writes it: `##fileformat` first, then the other lines sorted by their
/// text (contig lines in dictionary order), then `#CHROM` with the samples in sorted order.
/// `##GATKCommandLine` is replaced at its position by `##variant-rsCommandLine`.
fn m2_header(contigs: &[(String, i64, Option<String>)], samples: &[String], normals: &[String], tumors: &[String], cmdline: &str, date: &str) -> String {
    let quote = |v: &str| if v.contains(',') || v.contains(' ') { format!("\"{}\"", v.replace('"', "\\\"")) } else { v.to_string() };
    // (sort key, contig index, line)
    let mut lines: Vec<(String, usize, String)> = M2_HEADER_LINES.lines().map(|l| (l[2..].to_string(), 0, l.to_string())).collect();
    lines.push((
        "GATKCommandLine".into(),
        0,
        format!(
            "##variant-rsCommandLine=<ID=variant-rs,CommandLine=\"{}\",Version=\"{}\",Date=\"{}\">",
            cmdline.replace('"', "\\\""),
            env!("CARGO_PKG_VERSION"),
            date
        ),
    ));
    for (i, (n, l, a)) in contigs.iter().enumerate() {
        let mut t = format!("##contig=<ID={},length={}", quote(n), l);
        if let Some(a) = a {
            t.push_str(&format!(",assembly={}", quote(a)));
        }
        t.push('>');
        lines.push(("contig=<ID=".into(), i, t));
    }
    for n in normals {
        lines.push((format!("normal_sample={n}"), 0, format!("##normal_sample={n}")));
    }
    for t in tumors {
        lines.push((format!("tumor_sample={t}"), 0, format!("##tumor_sample={t}")));
    }
    lines.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut out = String::from("##fileformat=VCFv4.2\n");
    for (_, _, l) in lines {
        out.push_str(&l);
        out.push('\n');
    }
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for s in samples {
        out.push('\t');
        out.push_str(s);
    }
    out.push('\n');
    out
}

pub(crate) fn trace_t() -> f64 {
    static T0: std::sync::OnceLock<std::time::Instant> = std::sync::OnceLock::new();
    T0.get_or_init(std::time::Instant::now).elapsed().as_secs_f64()
}
