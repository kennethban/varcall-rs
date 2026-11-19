//! Rust port of GATK 4.7 HaplotypeCaller (default arguments; one sample in VCF or `-ERC GVCF` mode, several samples
//! called jointly in VCF mode).
//!
//! Whole-genome driver: BAM + BAI, any number of intervals (or all contigs), one GATK read shard per
//! contig (`AssemblyRegionWalker.makeReadShards`), several contigs processed concurrently, all sharing the
//! GPU PairHMM and the rayon pool, output assembled in dictionary order as bgzipped VCF/GVCF + tabix.
//!
//! Run-wide state: GATK's only state that crosses contigs is the global `java.util.Random`
//! (`Utils.getRandomGenerator()`, seed 47382911), drawn by the QD jitter (and by positional downsampling).
//! It is reset per *domain*: the whole run (single-run mode, = one GATK JVM over all intervals) or each
//! shard file (`--shard`, = one GATK JVM per shard file). QD jitter draws are deferred to output time, where
//! they are drawn in output order from the domain's generator, which is GATK's draw order.

mod active;
mod backend;
mod bam;
mod bamwrite;
mod bgzf;
mod bqsr;
mod cpuprof;
mod cram;
mod dbsnp;
mod filter;
mod frontend;
mod gatherbam;
mod gdb_merge;
mod gdb_store;
mod genotypegvcfs;
mod gatk_cli;
mod input;
mod intervals;
pub use hc_assembly::jmath; // bit-exact java.lang.Math (HotSpot x86_64)
mod markdup;
mod md5;
mod mergevcfs;
mod mutect;
mod output;
mod plan;
mod reblock;
mod reference;
mod rng;
mod sortsam;
mod tabix;

use active::{FilterCounts, Region, Tables};
use frontend::{FilledRegion, RefSeq};
use intervals::{Dict, Interval};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

const USAGE: &str = "usage: variant-rs <GATK tool> [GATK arguments]  (HaplotypeCaller, BaseRecalibrator, ApplyBQSR, \
MarkDuplicates; `variant-rs <tool> --help` lists what is supported; GATK_FALLBACK=/path/to/gatk runs GATK for the rest)\n\
   or: variant-rs -I in.bam -R ref.fa [-L interval|file]... [--shard file.intervals]... \
-O out.vcf[.gz] [-ERC GVCF] [--threads N] [--contig-workers N] [--pairhmm auto|cuda|cpu] \
[--shard-output-dir DIR] [--assembly-region-out regions.igv] [--gvcf-raw raw.vcf] [--compression-level N] \
[--scatter-intervals (each -L interval is its own GATK run)] \
[--tmp-dir DIR] [--max-reads-per-alignment-start N (50; 0 = off)] [--max-workers-per-shard N (2)] \
[--dump-activity file] [--vcf-header file] [--bqsr-recal-file recal.table (ApplyBQSR on the fly)]";

struct Args {
    /// `-I` (several: merged as GATK does; their samples are called jointly)
    bams: Vec<String>,
    fasta: String,
    intervals: Vec<String>,
    shards: Vec<String>,
    out: Option<String>,
    erc: bool,
    region_out: Option<String>,
    gvcf_raw: Option<String>,
    shard_dir: Option<String>,
    vcf_header: Option<String>,
    dump: Option<String>,
    pairhmm: String,
    threads: usize,
    workers: usize,
    level: i32,
    tmp: Option<String>,
    max_reads_per_start: usize,
    max_workers_per_shard: usize,
    /// `--scatter-intervals`: every `-L` interval is its own GATK run (as nf-core/sarek scatters HaplotypeCaller)
    scatter: bool,
    bqsr_table: Option<String>,
    padding: i64,
    excludes: Vec<String>,
    excl_padding: i64,
    dbsnp: Option<String>,
    gq_bands: Vec<i32>,
    contamination: f64,
    as_annotations: bool,
}

fn parse_args(a: &[String]) -> Args {
    let mut r = Args {
        bams: Vec::new(),
        fasta: String::new(),
        intervals: vec![],
        shards: vec![],
        out: None,
        erc: false,
        region_out: None,
        gvcf_raw: None,
        shard_dir: None,
        vcf_header: None,
        dump: None,
        pairhmm: "auto".into(),
        threads: 1,
        workers: 0,
        level: 5,
        tmp: None,
        max_reads_per_start: 50,
        max_workers_per_shard: 2,
        scatter: false,
        bqsr_table: None,
        padding: 0,
        excludes: vec![],
        excl_padding: 0,
        dbsnp: None,
        gq_bands: vec![],
        contamination: 0.0,
        as_annotations: false,
    };
    let mut i = 1;
    let val = |i: usize| -> String { a.get(i + 1).cloned().unwrap_or_else(|| panic!("{} needs a value\n{USAGE}", a[i])) };
    while i < a.len() {
        match a[i].as_str() {
            "-I" | "--input" => r.bams.push(val(i)),
            "-R" | "--reference" => r.fasta = val(i),
            "-L" | "--intervals" => r.intervals.push(val(i)),
            "--shard" => r.shards.push(val(i)),
            "-ip" | "--interval-padding" => r.padding = val(i).parse().unwrap(),
            "-XL" | "--exclude-intervals" => r.excludes.push(val(i)),
            "-ixp" | "--interval-exclusion-padding" => r.excl_padding = val(i).parse().unwrap(),
            "-D" | "--dbsnp" => r.dbsnp = Some(val(i)),
            "-contamination" | "--contamination-fraction-to-filter" => r.contamination = val(i).parse().unwrap(),
            // HaplotypeCaller's default groups are StandardAnnotation and StandardHCAnnotation; naming them again
            // changes nothing (GATKAnnotationPluginDescriptor merges by class name)
            "-G" | "--annotation-group" => match val(i).as_str() {
                "StandardAnnotation" | "StandardHCAnnotation" => {}
                "AS_StandardAnnotation" => r.as_annotations = true,
                g => panic!("--annotation-group {g} is not supported"),
            },
            "-GQB" | "--gvcf-gq-bands" => r.gq_bands.push(val(i).parse().unwrap_or_else(|_| panic!("--gvcf-gq-bands: not an integer"))),
            "-O" | "--output" => {
                let v = val(i);
                // legacy: `-O regions.igv` was the assembly-region output
                if v.ends_with(".igv") { r.region_out = Some(v) } else { r.out = Some(v) }
            }
            "--vcf" => r.out = Some(val(i)),
            "--gvcf" => {
                r.out = Some(val(i));
                r.erc = true;
            }
            "-ERC" | "--emit-ref-confidence" => {
                let v = val(i);
                r.erc = match v.as_str() {
                    "GVCF" => true,
                    "NONE" => false,
                    _ => panic!("only -ERC GVCF or NONE is supported"),
                };
            }
            "--assembly-region-out" => r.region_out = Some(val(i)),
            "--gvcf-raw" => r.gvcf_raw = Some(val(i)),
            "--shard-output-dir" => r.shard_dir = Some(val(i)),
            "--vcf-header" => r.vcf_header = Some(val(i)),
            "--dump-activity" => r.dump = Some(val(i)),
            "--pairhmm" => r.pairhmm = val(i),
            "--bqsr-recal-file" | "-bqsr" => r.bqsr_table = Some(val(i)),
            "--threads" => r.threads = val(i).parse().unwrap(),
            "--contig-workers" => r.workers = val(i).parse().unwrap(),
            "--compression-level" => r.level = val(i).parse().unwrap(),
            "--tmp-dir" => r.tmp = Some(val(i)),
            "--max-reads-per-alignment-start" => r.max_reads_per_start = val(i).parse().unwrap(),
            "--max-workers-per-shard" => r.max_workers_per_shard = val(i).parse().unwrap(),
            "--scatter-intervals" => {
                r.scatter = true;
                i += 1;
                continue;
            }
            "-h" | "--help" => {
                println!("{USAGE}");
                std::process::exit(0);
            }
            x => panic!("unknown argument {x}\n{USAGE}"),
        }
        i += 2;
    }
    if r.bams.is_empty() || r.fasta.is_empty() {
        eprintln!("{USAGE}");
        std::process::exit(2);
    }
    if let Err(e) = output::check_gq_bands(&r.gq_bands) {
        eprintln!("variant-rs: --gvcf-gq-bands: {e}");
        std::process::exit(2);
    }
    if r.as_annotations && !r.erc {
        // HaplotypeCallerEngine.initialize (UserException, exit code 2)
        eprintln!("A USER ERROR has occurred: Allele-specific annotations are not yet supported in the VCF mode");
        std::process::exit(2);
    }
    if !r.shards.is_empty() && !r.intervals.is_empty() {
        panic!("use either -L (single run) or --shard (one GATK run per shard file), not both");
    }
    if r.scatter && r.intervals.is_empty() {
        panic!("--scatter-intervals needs -L");
    }
    r
}

/// One contig of one domain: a GATK read shard.
struct Job {
    domain: usize,
    didx: usize, // index within the domain (RNG order)
    tid: usize, // BAM reference index
    name: String,
    len: i64,
    intervals: Vec<(i64, i64)>,
}

struct Domain {
    name: String,
    first: usize,
    last: usize, // inclusive job range; first > last for an empty domain
    rng: Arc<rng::DomainRng>,
    contam: Arc<rng::DomainRng<hc_genotype::well::Well19937c>>,
}

struct ContigResult {
    order: usize,
    vcf_text: Vec<u8>,
    vcf_resolved: bool,
    ds: hc_downsample::DownsampleStats,
    gvcf: Option<backend::Part>,
    raw: Option<PathBuf>,
    regions: Vec<Region>,
    filters: FilterCounts,
    loci: u64,
    active_loci: u64,
    assigned: u64,
    stats: backend::Stats,
    front_secs: f64,
    wall: f64,
    finished_at: f64,
}

fn domains_name(d: &[(String, Vec<Interval>)], i: usize) -> &str {
    &d[i].0
}

pub(crate) fn now_string() -> String {
    // UTC date without extra dependencies
    let secs = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_secs() as i64;
    let (days, rem) = (secs.div_euclid(86400), secs.rem_euclid(86400));
    let (y, m, d) = civil_from_days(days);
    format!("{y:04}-{m:02}-{d:02}T{:02}:{:02}:{:02}Z", rem / 3600, rem % 3600 / 60, rem % 60)
}

fn civil_from_days(z: i64) -> (i64, u32, u32) {
    let z = z + 719468;
    let era = z.div_euclid(146097);
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    (if m <= 2 { y + 1 } else { y }, m, d)
}

fn main() {
    // A panic on any worker thread (GPU, decoder, back end) ends the run with an error, before any output is
    // finalized; otherwise the threads waiting on the dead one would wait forever. Assembly's own panics are
    // caught where they happen (GATK assembly exceptions) and turned into a run failure by the caller, as GATK's
    // rethrow does.
    let default_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        default_hook(info);
        if !hc_assembly::panic_is_caught() {
            eprintln!("variant-rs: fatal error on thread {:?}, exiting", std::thread::current().name().unwrap_or("?"));
            std::process::exit(101);
        }
    }));
    // the `gatk` launcher's own `--java-options "..."` before the tool name (`gatk --java-options "-Xmx4g" Tool ...`,
    // as nf-core modules write it): no effect here (a GATK fallback gets the original arguments)
    let mut argv: Vec<String> = gatk_cli::strip_launcher_options(std::env::args().collect());
    // the `gatk` launcher's own requests: `--version` (nf-core records it with `sed -n '/GATK.*v/s/.*v//p'`), and
    // `--list` / `--help` / no arguments, which GATK answers (GATK_FALLBACK) when it is there
    match argv.get(1).map(|s| s.as_str()) {
        Some("--version") => {
            println!("The Genome Analysis Toolkit (GATK) v{}", gatk_cli::GATK_VERSION);
            println!("variant-rs {} (Rust/CUDA port of the {} tools)", env!("CARGO_PKG_VERSION"), gatk_cli::GATK_VERSION);
            return;
        }
        Some("--list") | Some("--help") | Some("-h") | None if std::env::var_os("GATK_FALLBACK").is_some() => gatk_cli::run_gatk_instead(),
        _ => {}
    }
    // `variant-rs <GatkTool> ...`: GATK's command line, translated to the tool's own arguments (or run by GATK)
    if argv.len() > 1 && !gatk_cli::is_tool_name(&argv[1]) && argv[1].starts_with(|c: char| c.is_ascii_uppercase()) && std::env::var_os("VARIANT_RS_DRY_RUN").is_some_and(|v| !v.is_empty() && v != "0") {
        println!("REFUSE {} is not supported by variant-rs (supported: {})", argv[1], gatk_cli::tool_names().join(", "));
        return;
    }
    if argv.len() > 1 && gatk_cli::is_tool_name(&argv[1]) {
        let translated = gatk_cli::dispatch(&argv);
        argv = std::iter::once(argv[0].clone()).chain(translated).collect();
    }
    if argv.get(1).map(|s| s.as_str()) == Some("mutect2") {
        mutect::run(&argv[2..]);
        return;
    }
    if argv.get(1).map(|s| s.as_str()) == Some("learnreadorientationmodel") {
        mutect::orientation::run_lrom(&argv[2..]);
        return;
    }
    if argv.get(1).map(|s| s.as_str()) == Some("getpileupsummaries") {
        mutect::contam::run_gps(&argv[2..]);
        return;
    }
    if argv.get(1).map(|s| s.as_str()) == Some("calculatecontamination") {
        mutect::contam::run_cc(&argv[2..]);
        return;
    }
    if argv.get(1).map(|s| s.as_str()) == Some("filtermutectcalls") {
        mutect::fmc::run(&argv[2..]);
        return;
    }
    if argv.get(1).map(|s| s.as_str()) == Some("filter") {
        filter::run(&argv[2..]);
        return;
    }
    if argv.get(1).map(|s| s.as_str()) == Some("mergemutectstats") {
        mutect::gathers::run_merge_stats(&argv[2..]);
        return;
    }
    if argv.get(1).map(|s| s.as_str()) == Some("gatherpileupsummaries") {
        mutect::gathers::run_gather_pileups(&argv[2..]);
        return;
    }
    if argv.get(1).map(|s| s.as_str()) == Some("gatherbqsrreports") {
        bqsr::gather::run(&argv[2..]);
        return;
    }
    if argv.get(1).map(|s| s.as_str()) == Some("gatherbamfiles") {
        gatherbam::run(&argv[2..]);
        return;
    }
    if argv.get(1).map(|s| s.as_str()) == Some("mergevcfs") {
        mergevcfs::run(&argv[2..]);
        return;
    }
    if argv.get(1).map(|s| s.as_str()) == Some("reblockgvcf") {
        reblock::run(&argv[2..]);
        return;
    }
    if argv.get(1).map(|s| s.as_str()) == Some("gdbdump") {
        // debugging: every cell of every array of a GenomicsDB workspace, attribute bytes in hex
        let ws = std::path::Path::new(&argv[2]);
        let mut arrays: Vec<_> = std::fs::read_dir(ws).unwrap().map(|e| e.unwrap().path()).filter(|p| p.join("__array_schema.tdb").exists()).collect();
        arrays.sort();
        use std::io::Write;
        let mut out = std::io::BufWriter::new(std::io::stdout().lock());
        for a in arrays {
            let arr = gdb_store::open_array(&a);
            gdb_store::for_each_cell(&arr, i64::MIN, i64::MAX, |c| {
                let h: Vec<String> = c.vals.iter().map(|v| v.iter().map(|b| format!("{b:02x}")).collect()).collect();
                writeln!(out, "{} {} {}", c.row, c.col, h.join(" ")).unwrap();
            });
        }
        return;
    }
    if argv.get(1).map(|s| s.as_str()) == Some("genotypegvcfs") {
        genotypegvcfs::run(&argv[2..]);
        return;
    }
    if argv.get(1).map(|s| s.as_str()) == Some("sortsam") {
        sortsam::run(&argv[2..]);
        return;
    }
    if argv.get(1).map(|s| s.as_str()) == Some("markdup") {
        markdup::run(&argv[2..]);
        return;
    }
    if argv.get(1).map(|s| s.as_str()) == Some("applybqsr") {
        bqsr::apply::run(&argv[2..]);
        return;
    }
    if argv.get(1).map(|s| s.as_str()) == Some("bqsr") {
        bqsr::run(&argv[2..]);
        return;
    }
    if argv.get(1).map(|s| s.as_str()) == Some("plan") {
        plan::run(&argv[2..]);
        return;
    }
    let args = parse_args(&argv);
    let t_all = Instant::now();
    rayon::ThreadPoolBuilder::new().num_threads(args.threads).thread_name(|_| "rayon".into()).build_global().unwrap();
    let tables = Tables::load();

    // inputs
    let (inputs, headers): (Vec<input::Input>, Vec<bam::Header>) = args.bams.iter().map(|b| input::Input::open(b, &args.fasta)).unzip();
    let header = &headers[0];
    // ReadUtils.getSamplesFromHeader over the merged header (a TreeSet: sorted names); several samples are called
    // jointly (VCF mode only, as in GATK)
    let samples: Vec<String> = headers.iter().flat_map(|h| h.sample_names()).collect::<std::collections::BTreeSet<String>>().into_iter().collect();
    let rg_sample: Option<Vec<std::collections::HashMap<Vec<u8>, u16>>> = (samples.len() > 1).then(|| {
        headers
            .iter()
            .map(|h| h.read_group_samples().into_iter().map(|(id, sm)| (id.into_bytes(), samples.iter().position(|s| *s == sm).unwrap() as u16)).collect())
            .collect()
    });
    if samples.len() > 1 {
        if args.erc {
            gatk_cli::refuse("-ERC GVCF with several samples (GATK: \"Can only be used in single sample mode currently\")");
        }
        if samples.len() > u16::MAX as usize {
            gatk_cli::refuse("too many samples");
        }
        eprintln!("{} samples: {}", samples.len(), samples.join(", "));
    }
    for (b, h) in args.bams.iter().zip(&headers).skip(1) {
        if h.ref_names != header.ref_names || h.ref_lens != header.ref_lens {
            gatk_cli::refuse(&format!("{b}: its sequence dictionary differs from {}'s; several inputs need the same one", args.bams[0]));
        }
    }
    if inputs.len() > 1 && args.bqsr_table.is_some() {
        gatk_cli::refuse("--bqsr-recal-file takes one input");
    }
    let bqsr: Option<Arc<bqsr::apply::Recalibrator>> = args.bqsr_table.as_ref().map(|t| Arc::new(bqsr::apply::Recalibrator::new(t, &header.text)));
    let fasta = reference::Fasta::open(&args.fasta);
    let ref_dict = Dict::from_reference(&args.fasta);
    let bam_dict = Dict::new(header.ref_names.clone(), header.ref_lens.clone());
    for (n, &l) in bam_dict.names.iter().zip(&bam_dict.lens) {
        if let Some(&t) = ref_dict.index.get(n) {
            assert_eq!(ref_dict.lens[t], l, "contig {n}: length differs between BAM and reference");
        }
    }

    // domains and jobs: one job per (domain, contig), in output order
    // IntervalArgumentCollection.parseIntervals: -L padded by --interval-padding, minus -XL padded by
    // --interval-exclusion-padding (-XL alone: the whole dictionary minus -XL)
    let exclude: Vec<Interval> = {
        let mut x = Vec::new();
        for l in &args.excludes {
            x.extend(intervals::parse_l_arg(l, &ref_dict).unwrap_or_else(|e| panic!("{e}")));
        }
        intervals::pad(&mut x, args.excl_padding, &ref_dict);
        intervals::sort_merge(x)
    };
    let finish_ivs = |mut iv: Vec<Interval>| -> Vec<Interval> {
        intervals::pad(&mut iv, args.padding, &ref_dict);
        let iv = intervals::sort_merge(iv);
        if exclude.is_empty() {
            return iv;
        }
        let out = intervals::subtract(&iv, &exclude);
        if out.is_empty() {
            eprintln!("A USER ERROR has occurred: Argument -L,-XL has a bad value: {:?}, {:?}. The intervals specified for exclusion with -XL removed all territory specified by -L.", args.intervals, args.excludes);
            std::process::exit(2);
        }
        out
    };
    let mut domain_ivs: Vec<(String, Vec<Interval>)> = Vec::new();
    if !args.shards.is_empty() {
        for s in &args.shards {
            let iv = intervals::read_interval_file(s, &ref_dict).unwrap_or_else(|e| panic!("{e}"));
            let name = std::path::Path::new(s).file_stem().unwrap().to_string_lossy().to_string();
            domain_ivs.push((name, finish_ivs(iv)));
        }
    } else if !args.intervals.is_empty() {
        let mut all = Vec::new();
        for l in &args.intervals {
            all.extend(intervals::parse_l_arg(l, &ref_dict).unwrap_or_else(|e| panic!("{e}")));
        }
        if args.scatter {
            // one GATK run per interval (sarek: one HaplotypeCaller task per interval BED), in dictionary order.
            // Overlapping intervals are merged (separate runs would call the overlap twice); abutting ones stay
            // separate runs, as sarek's tasks are.
            all.sort();
            let mut merged: Vec<Interval> = Vec::with_capacity(all.len());
            for i in all {
                match merged.last_mut() {
                    Some(l) if l.tid == i.tid && i.start <= l.end => l.end = l.end.max(i.end),
                    _ => merged.push(i),
                }
            }
            for iv in merged {
                domain_ivs.push((format!("{}_{}-{}", ref_dict.names[iv.tid], iv.start, iv.end), finish_ivs(vec![iv])));
            }
        } else {
            domain_ivs.push(("all".into(), finish_ivs(all)));
        }
    } else if !exclude.is_empty() {
        let all: Vec<Interval> = bam_dict.names.iter().filter_map(|n| ref_dict.index.get(n)).map(|&t| Interval { tid: t, start: 1, end: ref_dict.lens[t] }).collect();
        domain_ivs.push(("all".into(), finish_ivs(all)));
    } else {
        // no -L: every contig of the reads' sequence dictionary, in its order (AssemblyRegionWalker.onStartup uses
        // getAllIntervalsForReference(reads dictionary); GATK accepts a BAM whose @SQ order differs from the .dict)
        let all: Vec<Interval> = bam_dict
            .names
            .iter()
            .map(|n| {
                let t = *ref_dict.index.get(n).unwrap_or_else(|| panic!("BAM contig {n} not in the reference"));
                Interval { tid: t, start: 1, end: ref_dict.lens[t] }
            })
            .collect();
        domain_ivs.push(("all".into(), all));
    }
    // contig rank in traversal order: the reads' order without -L / -XL / --shard, else the reference's
    let reads_order = args.shards.is_empty() && args.intervals.is_empty() && args.excludes.is_empty();
    let rank = |name: &str| if reads_order { bam_dict.index[name] } else { ref_dict.index[name] };
    let mut jobs: Vec<Job> = Vec::new();
    let mut domains: Vec<Domain> = Vec::new();
    for (di, (name, ivs)) in domain_ivs.iter().enumerate() {
        let first = jobs.len();
        for (rt, list) in intervals::group_by_contig(ivs) {
            let cname = &ref_dict.names[rt];
            let tid = *bam_dict.index.get(cname).unwrap_or_else(|| panic!("contig {cname} is not in the BAM header"));
            jobs.push(Job { domain: di, didx: jobs.len() - first, tid, name: cname.clone(), len: ref_dict.lens[rt], intervals: list });
        }
        let n = jobs.len() - first;
        domains.push(Domain { name: name.clone(), first, last: jobs.len().wrapping_sub(1), rng: Arc::new(rng::DomainRng::new(n)), contam: Arc::new(rng::DomainRng::with(n, hc_genotype::well::Well19937c::gatk())) });
    }
    if args.out.is_some() {
        // the merged output is the jobs' outputs in job order: they must be sorted and disjoint
        // (a contig may be split across shards, e.g. by `variant-rs plan`)
        for w in jobs.windows(2) {
            let (a, b) = (&w[0], &w[1]);
            let (ra, rb) = (rank(&a.name), rank(&b.name));
            let ok = ra < rb || (ra == rb && a.intervals.last().unwrap().1 < b.intervals[0].0);
            if !ok {
                panic!("shards {} and {} are out of order or overlap on {}: the merged output would not be sorted",
                       domains_name(&domain_ivs, a.domain), domains_name(&domain_ivs, b.domain), b.name);
            }
        }
    }
    if args.dump.is_some() && jobs.len() != 1 {
        panic!("--dump-activity needs a single contig");
    }
    eprintln!("{} contig job(s) in {} domain(s), {} interval bases", jobs.len(), domains.len(), jobs.iter().flat_map(|j| j.intervals.iter()).map(|i| i.1 - i.0 + 1).sum::<i64>());

    let hc_opts = Arc::new(backend::HcOpts {
        gq_bands: if args.gq_bands.is_empty() { hc_refconf::default_gq_bands() } else { args.gq_bands.clone() },
        dbsnp: args.dbsnp.as_ref().map(|p| dbsnp::Dbsnp::load(p, &jobs.iter().map(|j| j.name.clone()).collect())),
        contamination: args.contamination,
        as_annotations: args.as_annotations,
        n_samples: samples.len().max(1),
    });

    // outputs
    let gz = |p: &str| output::block_compressed(p);
    let level_for = |p: &str| if gz(p) { Some(args.level) } else { None };
    let contig_lines: Vec<output::ContigLine> = {
        let asm = header.sq_assembly();
        bam_dict.names.iter().zip(&bam_dict.lens).map(|(n, &l)| output::ContigLine { name: n.clone(), len: l, assembly: asm.get(n).cloned() }).collect()
    };
    let header_samples: Vec<String> = if samples.is_empty() { vec!["SAMPLE".into()] } else { samples.clone() };
    let cmdline: String = std::env::args().collect::<Vec<_>>().join(" ");
    let date = now_string();
    let header_text = |template: &str| -> String {
        match &args.vcf_header {
            Some(h) => std::fs::read_to_string(h).unwrap(),
            None => output::make_header(template, &contig_lines, &header_samples, &cmdline, &date),
        }
    };
    // header lines that depend on the options
    let mut extra: Vec<String> = Vec::new();
    let mut drop_lines: Vec<&str> = Vec::new();
    if args.erc && !args.gq_bands.is_empty() {
        drop_lines.push("##GVCFBlock");
        extra.extend(output::gvcf_block_lines(&args.gq_bands));
    }
    if args.dbsnp.is_some() {
        extra.push(dbsnp::DB_HEADER.to_string());
    }
    if args.as_annotations {
        extra.extend(output::AS_GVCF_HEADER.iter().map(|l| l.to_string()));
    }
    let template_owned = output::template_with(if args.erc { output::GVCF_TEMPLATE } else { output::VCF_TEMPLATE }, &drop_lines, &extra);
    let template = template_owned.as_str();
    let main_header = header_text(template);
    // all parts share one compression setting: the main output's (or the shard outputs' if no main output)
    let part_level: Option<i32> = match (&args.out, &args.shard_dir) {
        (Some(o), _) => level_for(o),
        (None, Some(_)) => Some(args.level),
        (None, None) => None,
    };
    if let (Some(o), Some(_)) = (&args.out, &args.shard_dir) {
        if !gz(o) {
            panic!("with --shard-output-dir the main output must be .gz too");
        }
    }
    let tmp_dir: PathBuf = match &args.tmp {
        Some(t) => PathBuf::from(t).join(format!("variant-rs.{}", std::process::id())),
        None => {
            let base = args.out.clone().or_else(|| args.shard_dir.as_ref().map(|d| format!("{d}/x"))).unwrap_or_else(|| "variant-rs".into());
            PathBuf::from(format!("{base}.parts.{}", std::process::id()))
        }
    };
    if args.erc {
        std::fs::create_dir_all(&tmp_dir).unwrap();
    }

    let use_gpu = args.pairhmm != "cpu" && hc_cudahmm::init(false);
    if args.pairhmm == "cuda" && !use_gpu {
        panic!("--pairhmm cuda requested but no CUDA device (use --pairhmm cpu or auto)");
    }
    // contig workers: half the threads up to 8, and a quarter of them beyond (up to 32). With the GPU PairHMM's engine
    // pool more workers keep the device busier: 30x whole genome, L40S, 64 threads: 8 workers 160 s, 16 workers 121 s
    // (128 threads, 32 workers: 114 s); on 32 threads (RTX 5070 Ti) the CPU stages set the pace and 16 workers gain
    // nothing (403 s against 401 s) for 6.6 GB more memory, so up to 32 threads the count is as before.
    let n_workers = if args.workers > 0 { args.workers } else { (args.threads / 2).clamp(1, 8).max((args.threads / 4).min(32)) }.min(jobs.len().max(1));
    // two-stage front end (read stage on its own thread) when the threads allow it; VARIANT_RS_TWO_STAGE=0/1 overrides
    let two_stage = match std::env::var("VARIANT_RS_TWO_STAGE").ok().as_deref() {
        Some("0") => false,
        Some(_) => true,
        None => n_workers * 4 <= args.threads,
    };
    eprintln!(
        "threads {} | contig workers {} | front end {} | PairHMM {}",
        args.threads,
        n_workers,
        if two_stage { "two-stage" } else { "single-thread" },
        if use_gpu { "CUDA" } else { "CPU (GKL AVX-512 port)" }
    );

    // ---------------------------------------------------------------------------------------------
    // workers. A worker claims the next contig of a domain (contigs of one domain in order, at most
    // --max-workers-per-shard workers per domain unless there is a single domain); longest estimated first.
    // ---------------------------------------------------------------------------------------------
    // Longest-processing-time-first: domains by estimated duration (BAM-index cost of their intervals,
    // divided by the number of workers the domain can use at once), largest first. (Starting domains with more
    // contigs first, as their duration is not predicted by the cost, was tried and did not help: 2026-11-14
    // NOTES.)
    let per_dom_cap = if domains.len() == 1 { usize::MAX } else { args.max_workers_per_shard.max(1) };
    let model = {
        // read data per window, summed over the inputs
        let mut c = inputs[0].window_costs(&header.ref_lens);
        for inp in &inputs[1..] {
            for (a, b) in c.iter_mut().zip(inp.window_costs(&header.ref_lens)) {
                a.iter_mut().zip(b).for_each(|(x, y)| *x += y);
            }
        }
        plan::CostModel::new(c, 0.05)
    };
    let mut dom_order: Vec<usize> = (0..domains.len()).filter(|&d| domains[d].first <= domains[d].last && domains[d].last != usize::MAX).collect();
    let dom_est: Vec<f64> = (0..domains.len())
        .map(|d| {
            if !dom_order.contains(&d) {
                return 0.0;
            }
            let js = &jobs[domains[d].first..=domains[d].last];
            js.iter().map(|j| model.cost(j.tid, &j.intervals)).sum::<f64>() / (js.len().min(per_dom_cap) as f64)
        })
        .collect();
    dom_order.sort_by(|&a, &b| dom_est[b].partial_cmp(&dom_est[a]).unwrap());
    // (next unclaimed job per domain, active workers per domain)
    let sched: Mutex<(Vec<usize>, Vec<usize>)> = Mutex::new((domains.iter().map(|d| d.first).collect(), vec![0; domains.len()]));
    // time each domain's first contig was claimed (for the "domains" report)
    let dom_start: Mutex<Vec<f64>> = Mutex::new(vec![f64::NAN; domains.len()]);
    let claim = |prev_dom: Option<usize>| -> Option<usize> {
        let mut s = sched.lock().unwrap();
        if let Some(d) = prev_dom {
            s.1[d] -= 1;
        }
        for &d in &dom_order {
            if s.0[d] <= domains[d].last && s.1[d] < per_dom_cap {
                let k = s.0[d];
                s.0[d] += 1;
                s.1[d] += 1;
                let mut st = dom_start.lock().unwrap();
                if st[d].is_nan() {
                    st[d] = t_all.elapsed().as_secs_f64();
                }
                return Some(k);
            }
        }
        None
    };
    let inflate_ns = Arc::new(AtomicU64::new(0));
    let (res_tx, res_rx) = std::sync::mpsc::channel::<ContigResult>();
    let dump_file: Mutex<Option<BufWriter<File>>> = Mutex::new(args.dump.as_ref().map(|p| BufWriter::with_capacity(1 << 20, File::create(p).unwrap())));
    let failed = std::sync::atomic::AtomicBool::new(false);
    let _ = AtomicUsize::new(0);

    std::thread::scope(|scope| {
        for _ in 0..n_workers {
            let res_tx = res_tx.clone();
            let (jobs, domains, inputs, fasta, tables, header, claim) = (&jobs, &domains, &inputs, &fasta, &tables, header, &claim);
            let hc_opts = &hc_opts;
            let (inflate_ns, tmp_dir, dump_file, failed) = (&inflate_ns, &tmp_dir, &dump_file, &failed);
            let raw_wanted = args.gvcf_raw.is_some();
            let dump_on = args.dump.is_some();
            let erc = args.erc;
            let cap = args.max_reads_per_start;
            let bqsr = &bqsr;
            let (rg_sample, samples) = (&rg_sample, &samples);
            std::thread::Builder::new().name("contig-front".into()).spawn_scoped(scope, move || {
                let mut prev: Option<usize> = None;
                loop {
                    let Some(k) = claim(prev) else { return };
                    if failed.load(Ordering::Relaxed) {
                        return;
                    }
                    let job = &jobs[k];
                    prev = Some(job.domain);
                    let domain = domains[job.domain].rng.clone();
                    let t0 = Instant::now();
                    // reference: the intervals +/- 700 (regions are within intervals; assembly needs the padded
                    // region +/- 500, region padding is 100)
                    let (lo, hi) = (job.intervals[0].0 - 700, job.intervals.last().unwrap().1 + 700);
                    let (seq, clen) = fasta.fetch(&job.name, lo, hi);
                    assert_eq!(clen, job.len);
                    let refseq = Arc::new(RefSeq { from: lo.max(1), seq, contig_len: clen });
                    let safe = job.name.replace(|c: char| !c.is_ascii_alphanumeric() && c != '_' && c != '-', "_");
                    let cfg = backend::BackendCfg {
                        contig: job.name.clone(),
                        tid: job.tid as i32,
                        contig_len: clen,
                        refseq: refseq.clone(),
                        use_gpu,
                        erc,
                        gvcf_part: if erc { Some((tmp_dir.join(format!("{k:05}.{safe}.part")), part_level)) } else { None },
                        raw_part: if erc && raw_wanted { Some(tmp_dir.join(format!("{k:05}.{safe}.raw"))) } else { None },
                        domain: domain.clone(),
                        domain_idx: job.didx,
                        opts: hc_opts.clone(),
                        contam: domains[job.domain].contam.clone(),
                    };
                    // front -> back end queue (batches of FilledRegion); VARIANT_RS_BACK_QUEUE overrides the depth
                    let depth: usize = std::env::var("VARIANT_RS_BACK_QUEUE").ok().and_then(|v| v.parse().ok()).unwrap_or(2);
                    let (btx, brx) = std::sync::mpsc::sync_channel::<backend::BackMsg>(depth);
                    let back = std::thread::Builder::new().name(format!("back-{k}")).spawn(move || backend::run(cfg, brx)).unwrap();
                    let ctx = frontend::FrontCtx {
                        inputs,
                        ref_lens: &header.ref_lens,
                        tables,
                        parallel_inflate: true,
                        inflate_ns: inflate_ns.clone(),
                        bqsr: bqsr.clone(),
                        rg_sample: rg_sample.as_deref(),
                        samples,
                    };
                    let mut dump_guard = if dump_on { Some(dump_file.lock().unwrap()) } else { None };
                    let dump = dump_guard.as_mut().and_then(|g| g.as_mut());
                    let tf = Instant::now();
                    let sink = Sink { erc, tx: btx, batch: Vec::new(), n_active: 0, rng: None, resolved: 0, domain: domain.clone(), didx: job.didx };
                    let (fr, mut sink) = frontend::run(&ctx, job.tid, &job.name, &job.intervals, &refseq, cap, dump, sink, two_stage);
                    sink.send_batch();
                    if erc {
                        match sink.rng.take() {
                            Some(r) => domain.finish(job.didx, r),
                            None => domain.pass(job.didx),
                        }
                    } else {
                        sink.tx.send(backend::BackMsg::End { rng: sink.rng.take() }).unwrap();
                    }
                    drop(sink);
                    drop(dump_guard);
                    let front_secs = tf.elapsed().as_secs_f64();
                    let out = match back.join() {
                        Ok(o) => o,
                        Err(_) => {
                            failed.store(true, Ordering::Relaxed);
                            panic!("back end of {} failed", job.name);
                        }
                    };
                    res_tx
                        .send(ContigResult {
                            order: k,
                            vcf_text: out.vcf_text,
                            vcf_resolved: out.vcf_resolved,
                            ds: fr.ds,
                            gvcf: out.gvcf,
                            raw: out.raw,
                            regions: fr.regions,
                            filters: fr.filters,
                            loci: fr.n_loci,
                            active_loci: fr.n_active_loci,
                            assigned: fr.assigned,
                            stats: out.stats,
                            front_secs,
                            wall: t0.elapsed().as_secs_f64(),
                            finished_at: t_all.elapsed().as_secs_f64(),
                        })
                        .unwrap();
                }
            }).unwrap();
        }
        drop(res_tx);

        // -----------------------------------------------------------------------------------------
        // results: QD draws of contigs that never held the generator (in domain order), then output in
        // job order
        // -----------------------------------------------------------------------------------------
        let mut main_out = args.out.as_ref().map(|p| output::OutFile::create(p, &main_header, level_for(p)).unwrap());
        let mut raw_out = args.gvcf_raw.as_ref().map(|p| BufWriter::with_capacity(4 << 20, File::create(p).unwrap()));
        if let Some(r) = raw_out.as_mut() {
            r.write_all(header_text(output::RAW_TEMPLATE).as_bytes()).unwrap();
        }
        let mut shard_out: Option<output::OutFile> = None;
        let mut unresolved: BTreeMap<usize, ContigResult> = BTreeMap::new();
        let mut ready: BTreeMap<usize, ContigResult> = BTreeMap::new();
        let mut next = 0usize;
        let mut regions_out: Vec<(usize, Vec<Region>)> = Vec::new();
        let mut tot = backend::Stats::default();
        let mut fc = FilterCounts::default();
        let mut ds_tot = hc_downsample::DownsampleStats::default();
        let (mut loci, mut active_loci, mut assigned, mut front_secs) = (0u64, 0u64, 0u64, 0f64);
        let mut draws = 0u64;
        let mut slowest: Vec<(f64, String)> = Vec::new();
        let mut dom_done: Vec<f64> = vec![0.0; domains.len()];
        let open_domain = |d: usize, shard_out: &mut Option<output::OutFile>| {
            if let Some(dir) = &args.shard_dir {
                std::fs::create_dir_all(dir).unwrap();
                let ext = if args.erc { "g.vcf.gz" } else { "vcf.gz" };
                let p = format!("{dir}/{}.{ext}", domains[d].name);
                *shard_out = Some(output::OutFile::create(&p, &main_header, Some(args.level)).unwrap());
            }
        };
        let write = |r: ContigResult, main_out: &mut Option<output::OutFile>, shard_out: &mut Option<output::OutFile>, raw_out: &mut Option<BufWriter<File>>| {
            let job = &jobs[r.order];
            let dom = &domains[job.domain];
            if r.order == dom.first {
                open_domain(job.domain, shard_out);
            }
            if !args.erc {
                if let Some(o) = main_out.as_mut() {
                    o.add_text(&r.vcf_text).unwrap();
                }
                if let Some(o) = shard_out.as_mut() {
                    o.add_text(&r.vcf_text).unwrap();
                }
            } else if let Some(p) = &r.gvcf {
                if let Some(o) = main_out.as_mut() {
                    o.add_part(p).unwrap();
                }
                if let Some(o) = shard_out.as_mut() {
                    o.add_part(p).unwrap();
                }
                std::fs::remove_file(&p.path).ok();
            }
            if let (Some(p), Some(w)) = (&r.raw, raw_out.as_mut()) {
                std::io::copy(&mut File::open(p).unwrap(), w).unwrap();
                std::fs::remove_file(p).ok();
            }
            if r.order == dom.last {
                if let Some(o) = shard_out.take() {
                    o.finish().unwrap();
                }
            }
        };
        for r in res_rx {
            if r.vcf_resolved {
                ready.insert(r.order, r);
            } else {
                unresolved.insert(r.order, r);
            }
            // contigs whose QD draws wait for the generator: draw them when their turn comes
            loop {
                let mut progressed = false;
                let keys: Vec<usize> = unresolved.keys().copied().collect();
                for k in keys {
                    let job = &jobs[k];
                    if let Some(mut g) = domains[job.domain].rng.try_take(job.didx) {
                        let mut r = unresolved.remove(&k).unwrap();
                        let mut text = Vec::with_capacity(r.vcf_text.len() + 4096);
                        let n = rng::resolve_qd(&r.vcf_text, &mut g, &mut text);
                        r.stats.qd_drawn += n;
                        r.vcf_text = text;
                        domains[job.domain].rng.finish(job.didx, g);
                        ready.insert(k, r);
                        progressed = true;
                    }
                }
                if !progressed {
                    break;
                }
            }
            while let Some(mut r) = ready.remove(&next) {
                tot.add(&r.stats);
                fc.add(&r.filters);
                loci += r.loci;
                active_loci += r.active_loci;
                assigned += r.assigned;
                front_secs += r.front_secs;
                draws += r.stats.qd_drawn;
                ds_tot.submitted += r.ds.submitted;
                ds_tot.kept += r.ds.kept;
                ds_tot.drawing_groups += r.ds.drawing_groups;
                ds_tot.draws += r.ds.draws;
                slowest.push((r.wall, jobs[r.order].name.clone()));
                let dd = &mut dom_done[jobs[r.order].domain];
                *dd = dd.max(r.finished_at);
                regions_out.push((r.order, std::mem::take(&mut r.regions)));
                write(r, &mut main_out, &mut shard_out, &mut raw_out);
                next += 1;
            }
        }
        if next < jobs.len() {
            panic!("only {next} of {} contig jobs finished", jobs.len());
        }
        // shard files of domains without any job
        for (d, dom) in domains.iter().enumerate() {
            if dom.first > dom.last || dom.last == usize::MAX {
                open_domain(d, &mut shard_out);
                if let Some(o) = shard_out.take() {
                    o.finish().unwrap();
                }
            }
        }
        let main_lines = main_out.as_ref().map_or(0, |o| o.lines);
        let main_deflate = main_out.as_ref().map_or(0.0, |o| o.deflate_secs);
        if let Some(o) = main_out.take() {
            o.finish().unwrap();
        }
        if let Some(mut r) = raw_out.take() {
            r.flush().unwrap();
        }
        if args.erc {
            std::fs::remove_dir_all(&tmp_dir).ok();
        }
        if let Some(mut d) = dump_file.lock().unwrap().take() {
            d.flush().unwrap();
        }

        // assembly regions (IGV format of GATK --assembly-region-out)
        if let Some(p) = &args.region_out {
            let mut w = BufWriter::new(File::create(p).unwrap());
            writeln!(w, "#track graphType=line").unwrap();
            writeln!(w, "Chromosome\tStart\tEnd\tFeature\tAssemblyRegions").unwrap();
            for (k, regs) in &regions_out {
                let c = &jobs[*k].name;
                for r in regs {
                    writeln!(w, "{}\t{}\t{}\tend-marker\t0.00000", c, r.start - 1, r.start).unwrap();
                    writeln!(w, "{}\t{}\t{}\tsize={}\t{}", c, r.start - 1, r.end, r.end - r.start + 1, if r.active { "1.00000" } else { "-1.00000" }).unwrap();
                }
            }
            w.flush().unwrap();
        }

        // report
        let n_regions: usize = regions_out.iter().map(|(_, r)| r.len()).sum();
        let n_act: usize = regions_out.iter().map(|(_, r)| r.iter().filter(|x| x.active).count()).sum();
        let s = &tot;
        eprintln!("reads: total={} passed={} filtered: mapq={} unmapped={} secondary={} duplicate={} qcfail={} zeroRefLen={} badCigar={} malformed={}",
                  fc.total, fc.passed, fc.mapq, fc.unmapped, fc.secondary, fc.duplicate, fc.qcfail, fc.zero_ref_len, fc.bad_cigar, fc.malformed);
        eprintln!("downsampling (max {} reads per start): {} reads in, {} kept, {} drawing groups, {} nextInt draws",
                  args.max_reads_per_start, ds_tot.submitted, ds_tot.kept, ds_tot.drawing_groups, ds_tot.draws);
        eprintln!("loci={} with isActive>0: {}  regions={} (active {})  read->region assignments={}", loci, active_loci, n_regions, n_act, assigned);
        crate::cpuprof::report();
        hc_assembly::sw::dump::flush();
        if args.erc {
            eprintln!("GVCF: regions {} (active {}, genotyped {}), calls {}, GVCF lines {}", s.regions, s.active, s.genotyped, s.calls, s.out_lines);
        } else {
            eprintln!("VCF: regions genotyped {}, records {}, haplotypes {}, QD jitter draws {}", s.genotyped, s.records, s.haplotypes, draws);
        }
        if crate::mutect::front::ft_on() {
            for (n, c) in crate::frontend::HFT_NAMES.iter().zip(&crate::frontend::HFT) {
                eprintln!("front time {n}: {:.1} s", c.load(std::sync::atomic::Ordering::Relaxed) as f64 * 1e-9);
            }
        }
        eprintln!("summed stage wall times (contigs overlap): front end {:.1}s (BGZF inflate {:.1}s) | finalize {:.1}s | assembly {:.1}s | prep {:.1}s | PairHMM {:.1}s | finish {:.1}s | realign {:.1}s | genotyping {:.1}s | ref-confidence {:.1}s | GVCF writer busy {:.1}s (deflate {:.1}s) | output deflate {:.1}s",
                  front_secs, inflate_ns.load(Ordering::Relaxed) as f64 * 1e-9, s.preproc, s.asm, s.prep, s.hmm, s.finish, s.realign, s.geno, s.refconf, s.writer_busy, s.deflate, main_deflate);
        slowest.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        let top: Vec<String> = slowest.iter().take(5).map(|(t, n)| format!("{n} {t:.1}s")).collect();
        eprintln!("slowest contigs: {}", top.join(", "));
        if domains.len() > 1 {
            let st = dom_start.lock().unwrap();
            let v: Vec<String> = domains.iter().zip(&dom_done).zip(st.iter()).map(|((d, t), s)| format!("{} {:.0}-{:.0}s", d.name, s, t)).collect();
            eprintln!("domains started-finished at: {}", v.join(", "));
        }
        if use_gpu {
            hc_cudahmm::report();
        }
        if main_lines > 0 {
            eprintln!("output data lines: {main_lines}");
        }
        eprintln!("time: total={:.3}s (threads={}, contig workers={})", t_all.elapsed().as_secs_f64(), args.threads, n_workers);
    });
}

/// The worker side of a contig: batches filled regions for the back end, and fetches the domain's random
/// generator for the downsampler's draws (in VCF mode through the back end, which must first make the QD
/// draws GATK makes before them).
struct Sink {
    erc: bool,
    tx: std::sync::mpsc::SyncSender<backend::BackMsg>,
    batch: Vec<FilledRegion>,
    n_active: usize,
    rng: Option<hc_downsample::JavaRandom>,
    resolved: usize,
    domain: Arc<rng::DomainRng>,
    didx: usize,
}

impl Sink {
    fn send_batch(&mut self) {
        if !self.batch.is_empty() {
            let t = crate::frontend::hft_now();
            self.tx.send(backend::BackMsg::Batch(std::mem::take(&mut self.batch))).unwrap();
            crate::frontend::hft_add(10, t);
            self.n_active = 0;
        }
    }
}

impl frontend::FrontSink for Sink {
    fn region(&mut self, r: FilledRegion) {
        if !self.erc && !r.active {
            return;
        }
        self.n_active += r.active as usize;
        self.batch.push(r);
        if self.n_active >= backend::ASM_BATCH || self.batch.len() >= 8 * backend::ASM_BATCH {
            self.send_batch();
        }
    }

    fn rng_for_draw(&mut self, resolved_before: usize) -> &mut hc_downsample::JavaRandom {
        if self.erc {
            // GVCF mode makes no QD draws: only the order of contigs matters
            if self.rng.is_none() {
                self.rng = Some(self.domain.take(self.didx));
            }
        } else if self.rng.is_none() || resolved_before > self.resolved {
            self.send_batch();
            let (reply, rx) = std::sync::mpsc::sync_channel(1);
            let t = crate::frontend::hft_now();
            self.tx.send(backend::BackMsg::Resolve { before: resolved_before, rng: self.rng.take(), reply }).unwrap();
            self.rng = Some(rx.recv().expect("back end failed"));
            crate::frontend::hft_add(9, t);
            self.resolved = resolved_before;
        }
        self.rng.as_mut().unwrap()
    }
}
