//! align-rs: bwa-mem2 2.2.1 `mem`, output identical (see docs/notes/align-rs.md and docs/architecture.md).
//!
//! Per batch (`-K` bases, as fastmap.cpp reads it): seeds, chains and extended regions per read (512-read chunks in
//! parallel, as bwa-mem2's worker_bwt / worker_aln), then the batch's insert-size distribution, then pairing and
//! SAM per pair (worker_sam), written in input order.

mod bam;
mod chain;
mod ext;
mod fastq;
#[cfg(gpu)]
mod gpu;
mod index;
mod ksort;
mod ksw;
mod prof;
mod sam;
mod smem;
#[cfg(gpu)]
mod tune;

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use rayon::prelude::*;
use std::io::Write;

/// ALIGN_RS_PROFILE=1: summed thread time per stage (ns) and wall time of the sequential parts
static PROF: [std::sync::atomic::AtomicU64; 9] = [const { std::sync::atomic::AtomicU64::new(0) }; 9];
static DROP_PENDING: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
static DROP_MAX: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
const PROF_NAMES: [&str; 9] = ["wall: seeding (SMEM)", "wall: SA + chaining", "wall: extension (rounds + replay + dedup)", "pairing + SAM", "wall: FASTQ read + gunzip", "wall: seed/chain/extend phase", "wall: pestat + SAM phase", "wall: output write", "dropper busy"];
fn vb<T>(v: &Vec<T>) -> usize {
    v.capacity() * std::mem::size_of::<T>()
}
fn mbsz(what: &str, n_processed: i64, b: usize) {
    if n_processed == 2 * 666668 && std::env::var("ALIGN_RS_PROFILE").is_ok() {
        eprintln!("[align-rs] size {what}: {} MB", b >> 20);
    }
}
fn padd(i: usize, t: std::time::Instant) {
    PROF[i].fetch_add(t.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
}

/// ALIGN_RS_PROFILE: (stage, batch, start s, end s) of every stage run, printed as a timeline at the end
static TIMELINE: std::sync::Mutex<Vec<(&'static str, usize, f64, f64)>> = std::sync::Mutex::new(Vec::new());
static T_START: std::sync::OnceLock<std::time::Instant> = std::sync::OnceLock::new();
fn tl(stage: &'static str, batch: usize, start: std::time::Instant) {
    let t0 = *T_START.get().unwrap();
    TIMELINE.lock().unwrap().push((stage, batch, (start - t0).as_secs_f64(), t0.elapsed().as_secs_f64()));
}

const BATCH_SIZE: usize = 512;
/// GPU seeding: room for seeds per read (more: the CPU redoes the read)
#[cfg_attr(not(gpu), allow(dead_code))] // GPU build only
const SMEM_CAP: usize = 256;

/// batches a stage may finish ahead of the next one (ALIGN_RS_SLACK_<STAGE>=n, ALIGN_RS_SLACK=n for all; default 0)
fn slack(stage: &str) -> usize {
    let v = |k: &str| std::env::var(k).ok().and_then(|v| v.parse().ok());
    v(&format!("ALIGN_RS_SLACK_{}", stage.to_uppercase())).or_else(|| v("ALIGN_RS_SLACK")).unwrap_or(0)
}

/// a batch's identity (lane, first read), as gpu::batch_key
fn batch_key(lane: usize, n_processed: i64) -> u64 {
    ((lane as u64) << 40) | n_processed as u64
}

/// A raw pointer for disjoint parallel writes
#[derive(Clone, Copy)]
struct SendPtr<T>(*mut T);
unsafe impl<T> Send for SendPtr<T> {}
unsafe impl<T> Sync for SendPtr<T> {}
impl<T> SendPtr<T> {
    fn get(self) -> *mut T {
        self.0
    }
} // reads per kernel call (kthread.cpp chunks: fixed boundaries)

struct Args {
    threads: usize,
    chunk: i64,
    reference: String,
    r1: String,
    r2: Option<String>,
    rg: Option<String>,
    softclip: bool,
    dump: Option<String>,
    max_occ: i64,
    out: Option<String>,
    gpu: bool,
    level: u32,
    /// --lanes FILE: one line per lane, `out.bam<TAB>@RG line (escaped as for -R)<TAB>R1[<TAB>R2]`
    lanes: Option<String>,
}

fn parse(a: &[String]) -> Args {
    let mut r = Args { threads: 1, chunk: 10_000_000, reference: String::new(), r1: String::new(), r2: None, rg: None, softclip: false, dump: None, max_occ: 500, out: None, gpu: false, level: 1, lanes: None };
    let mut pos = Vec::new();
    let mut i = 0;
    let mut chunk_fixed = None;
    while i < a.len() {
        let v = |i: usize| a.get(i + 1).cloned().unwrap_or_else(|| panic!("{} needs a value", a[i]));
        match a[i].as_str() {
            "-t" => {
                r.threads = v(i).parse().unwrap();
                i += 1;
            }
            "-K" => {
                chunk_fixed = Some(v(i).parse::<i64>().unwrap());
                i += 1;
            }
            "-R" => {
                r.rg = Some(v(i));
                i += 1;
            }
            "-Y" => r.softclip = true,
            "--gpu" => r.gpu = true,
            "-o" => {
                r.out = Some(v(i));
                i += 1;
            }
            "-l" => {
                r.level = v(i).parse().unwrap();
                i += 1;
            }
            "--dump" => {
                r.dump = Some(v(i));
                i += 1;
            }
            "--lanes" => {
                r.lanes = Some(v(i));
                i += 1;
            }
            x if x.starts_with('-') => panic!("unsupported option {x}"),
            x => pos.push(x.to_string()),
        }
        i += 1;
    }
    r.chunk = chunk_fixed.unwrap_or(r.chunk * r.threads as i64);
    r.reference = pos[0].clone();
    if r.lanes.is_none() {
        r.r1 = pos.get(1).cloned().expect("no FASTQ given");
        r.r2 = pos.get(2).cloned();
    }
    r
}

/// One input: a FASTQ (pair), its read group and output
struct Lane {
    r1: String,
    r2: Option<String>,
    rg_line: Option<String>,
    sopt: sam::SamOpt,
    /// -o (SAM text when not .bam; stdout when None)
    out: Option<String>,
    /// the @PG command line: the run's own, or with --lanes the single-lane command giving the same output
    cl: String,
}

impl Lane {
    fn new(args: &Args, r1: String, r2: Option<String>, rg: Option<&str>, out: Option<String>, cl: String) -> Lane {
        let (rg_line, rg_id) = match rg {
            Some(s) => {
                let (l, id) = set_rg(s);
                (Some(l), Some(id))
            }
            None => (None, None),
        };
        Lane { r1, r2, rg_line, sopt: sam::SamOpt::bwa_default(args.softclip, rg_id), out, cl }
    }

    fn bam(&self) -> bool {
        self.out.as_ref().is_some_and(|p| p.ends_with(".bam"))
    }

    /// bwa_print_sam_hdr
    fn header(&self, bns: &index::Bns) -> String {
        let mut h = String::new();
        for a in &bns.anns {
            h.push_str(&format!("@SQ\tSN:{}\tLN:{}", a.name, a.len));
            h.push_str(if a.is_alt { "\tAH:*\n" } else { "\n" });
        }
        if let Some(l) = &self.rg_line {
            h.push_str(l);
            h.push('\n');
        }
        h.push_str(&format!("@PG\tID:align-rs\tPN:align-rs\tVN:0.1.0 (bwa-mem2 2.2.1 output)\tCL:{}\n", self.cl));
        h
    }
}

/// The lanes of this run: the command line's, or --lanes FILE's (BAM outputs only)
fn lanes(args: &Args, argv: &[String]) -> Vec<Lane> {
    let Some(f) = &args.lanes else {
        return vec![Lane::new(args, args.r1.clone(), args.r2.clone(), args.rg.as_deref(), args.out.clone(), argv.join(" "))];
    };
    assert!(args.out.is_none() && args.rg.is_none() && args.dump.is_none(), "--lanes: -o, -R and --dump are per lane (in the file)");
    let text = std::fs::read_to_string(f).unwrap_or_else(|e| panic!("{f}: {e}"));
    let mut v = Vec::new();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let c: Vec<&str> = line.split('\t').collect();
        assert!(c.len() == 3 || c.len() == 4, "{f}: expected out.bam, @RG line, R1[, R2]: {line}");
        assert!(c[0].ends_with(".bam"), "{f}: lane outputs are .bam: {}", c[0]);
        let r2 = c.get(3).filter(|x| !x.is_empty()).map(|x| x.to_string());
        // the command that writes this lane alone, byte for byte the same file
        let cl = format!(
            "{} mem -t {} -K {}{}{} -R {} -l {} {} {}{} -o {}",
            argv[0],
            args.threads,
            args.chunk,
            if args.softclip { " -Y" } else { "" },
            if args.gpu { " --gpu" } else { "" },
            c[1],
            args.level,
            args.reference,
            c[2],
            r2.as_ref().map_or(String::new(), |x| format!(" {x}")),
            c[0]
        );
        v.push(Lane::new(args, c[2].to_string(), r2, Some(c[1]), Some(c[0].to_string()), cl));
    }
    assert!(!v.is_empty(), "{f}: no lanes");
    // two lanes writing one BAM would truncate the first one's spill file under its final pass
    let canon = |p: &str| {
        let p = std::path::absolute(p).unwrap();
        let dir = p.parent().and_then(|d| std::fs::canonicalize(d).ok());
        dir.zip(p.file_name()).map_or(p.clone(), |(d, n)| d.join(n))
    };
    let mut outs: Vec<std::path::PathBuf> = v.iter().map(|l| canon(l.out.as_deref().unwrap())).collect();
    outs.sort();
    if let Some(w) = outs.windows(2).find(|w| w[0] == w[1]) {
        eprintln!("[align-rs] {f}: {} is the output of more than one lane", w[0].display());
        std::process::exit(1);
    }
    v
}

/// bwa_escape + bwa_set_rg: the header line and its ID
fn set_rg(s: &str) -> (String, String) {
    assert!(s.starts_with("@RG"), "the read group line is not started with @RG");
    let mut line = String::new();
    let mut it = s.chars();
    while let Some(c) = it.next() {
        if c == '\\' {
            match it.next() {
                Some('t') => line.push('\t'),
                Some('n') => line.push('\n'),
                Some('r') => line.push('\r'),
                Some('\\') => line.push('\\'),
                _ => {}
            }
        } else {
            line.push(c);
        }
    }
    let p = line.find("\tID:").expect("no ID at the read group line") + 4;
    let id: String = line[p..].chars().take_while(|&c| c != '\t' && c != '\n').collect();
    (line, id)
}

#[cfg(gpu)]
fn gpu_setup(gpu: bool, o: &ext::ExtOpt, bns: &index::Bns, fmi: &index::Fmi, ref_string: &[u8]) -> Option<gpu::Gpu> {
    if !gpu {
        return None;
    }
    // the device's 2-bit reference must be what .0123 holds (bwa-mem2 builds .0123 from .pac)
    // positions inside the reference only (a fixed probe past the end panicked on any reference under ~500 Mb)
    for &c in [0i64, 12345, bns.l_pac - 1, bns.l_pac, bns.l_pac * 2 - 1, 1_000_000_007].iter().filter(|&&c| c < bns.l_pac * 2) {
        let b = if c < bns.l_pac { bns.pac_base(c) } else { 3 - bns.pac_base(bns.l_pac * 2 - 1 - c) };
        assert_eq!(b, ref_string[c as usize], ".pac and .0123 disagree at {c}");
    }
    let g = gpu::Gpu::new(o, bns, fmi);
    if g.is_none() {
        eprintln!("[align-rs] no usable CUDA device: everything on the CPU");
    }
    g
}

fn main() {
    // a panic on any thread ends the process at once (abort, exit 134): unwinding would let the other stages run on
    // and finalize BAMs that look complete (a lane cut short, the lanes not reached written empty)
    let report = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        report(info);
        std::process::abort();
    }));
    let argv: Vec<String> = std::env::args().collect();
    if argv.len() < 2 || (argv[1] != "mem" && argv[1] != "tune") {
        eprintln!("usage: align-rs mem [-t N] [-K N] [-Y] [-R RG] [-o out.sam | -o out.bam [-l level]] [--dump FILE] ref.fa R1.fq[.gz] [R2.fq[.gz]]");
        eprintln!("       align-rs mem [-t N] [-K N] [-Y] [-l level] --lanes LANES.tsv ref.fa   (LANES.tsv: out.bam<TAB>@RG...<TAB>R1[<TAB>R2] per lane)");
        eprintln!("       align-rs tune [-t N] [-K N] ref.fa R1.fq[.gz] [R2.fq[.gz]]   (GPU settings measured on the first batch, saved for `mem`)");
        std::process::exit(1);
    }
    // tune: the first batch only, its output discarded, the GPU calls' inputs captured (tune.rs)
    let tuning = argv[1] == "tune";
    let mut args = parse(&argv[2..]);
    if tuning {
        assert!(args.lanes.is_none() && args.out.is_none() && args.dump.is_none(), "tune: no -o, --lanes or --dump");
        args.gpu = true;
        args.out = Some("/dev/null".to_string());
    }
    let max_batches = if tuning { 1 } else { usize::MAX };
    {
        // mimalloc returns freed memory to the OS after 10 ms (MADV_DONTNEED); each batch frees and re-allocates
        // gigabytes (the extension state alone ~2 GB), and every page faulted in again: ~10M faults on 1/6 of a
        // lane, 12% of the CPU time in the kernel. Purging with MADV_FREE instead keeps the pages mapped for the
        // next batch unless the kernel needs the memory (faults 10.6M -> 2.0M, 5% faster).
        // ALIGN_RS_MI_DECOMMIT=1 restores decommitting; ALIGN_RS_MI_PURGE=ms sets the purge delay (-1: never).
        extern "C" {
            fn mi_option_set(option: i32, value: std::ffi::c_long);
        }
        const MI_OPTION_PURGE_DECOMMITS: i32 = 5;
        const MI_OPTION_PURGE_DELAY: i32 = 15;
        if let Some(v) = std::env::var("ALIGN_RS_MI_PURGE").ok().and_then(|v| v.parse::<i64>().ok()) {
            unsafe { mi_option_set(MI_OPTION_PURGE_DELAY, v as std::ffi::c_long) };
        }
        let decommit = std::env::var("ALIGN_RS_MI_DECOMMIT").ok().and_then(|v| v.parse::<i64>().ok()).unwrap_or(0);
        unsafe { mi_option_set(MI_OPTION_PURGE_DECOMMITS, decommit as std::ffi::c_long) };
    }
    if std::env::var("ALIGN_RS_PTRACER").is_ok() {
        // profiling: let a sampler (gdb) attach although yama's ptrace_scope is 1 (PR_SET_PTRACER, PR_SET_PTRACER_ANY)
        extern "C" {
            fn prctl(option: i32, arg2: u64, arg3: u64, arg4: u64, arg5: u64) -> i32;
        }
        unsafe { prctl(0x59616d61, u64::MAX, 0, 0, 0) };
    }
    rayon::ThreadPoolBuilder::new().num_threads(args.threads).thread_name(|i| format!("aln-{i}")).build_global().unwrap();
    let t0 = std::time::Instant::now();
    T_START.set(t0).unwrap();
    let fmi = index::Fmi::load(&args.reference);
    let bns = index::Bns::load(&args.reference);
    let ref_map = index::Map::open(&format!("{}.0123", args.reference));
    assert_eq!(ref_map.bytes().len() as i64, bns.l_pac * 2, ".0123 length");
    let ref_string = ref_map.bytes();
    let opt = smem::Opt { min_seed_len: 19, split_factor: 1.5, split_width: 10, max_mem_intv: 20 };
    let copt = chain::ChainOpt { w: 100, max_chain_gap: 10000, max_occ: args.max_occ, mask_level: 0.5, drop_ratio: 0.5, min_seed_len: 19, min_chain_weight: 0, max_chain_extend: 1 << 30, a: 1 };
    let lanes = lanes(&args, &argv);
    let sopt = &lanes[0].sopt; // what does not depend on the read group
    if let Ok(p) = std::env::var("ALIGN_RS_DUMP_PAIR") {
        let _ = sam::PDUMP.set(std::sync::Mutex::new(std::fs::File::create(p).unwrap()));
    }

    // -o x.bam: coordinate-sorted BAM + x.bam.bai (bam.rs), per lane; otherwise SAM text (one lane)
    let mut out: Box<dyn Write + Send> = match &lanes[0].out {
        _ if lanes[0].bam() => Box::new(std::io::sink()),
        Some(p) => Box::new(std::io::BufWriter::with_capacity(16 << 20, std::fs::File::create(p).unwrap())),
        None => Box::new(std::io::BufWriter::with_capacity(16 << 20, std::io::stdout())),
    };
    if !lanes[0].bam() {
        out.write_all(lanes[0].header(&bns).as_bytes()).unwrap();
    }
    let ref_lens: Vec<i64> = bns.anns.iter().map(|a| a.len as i64).collect();
    let level = args.level;
    // before the GPU is set up: a tune starts from the device's defaults, not from a previous tune's cache
    #[cfg(gpu)]
    tune::ON.store(tuning, std::sync::atomic::Ordering::Relaxed);
    #[cfg(gpu)]
    let gpu_dev = gpu_setup(args.gpu, &sopt.e, &bns, &fmi, ref_string);
    #[cfg(not(gpu))]
    if args.gpu {
        eprintln!("[align-rs] built without CUDA (make -C cuda): everything on the CPU");
    }
    let cpu_ext = ext::CpuExtender { ref_string };
    #[cfg(gpu)]
    let gpu_ext = gpu_dev.as_ref().map(|g| gpu::GpuExtender { gpu: g, cpu: ext::CpuExtender { ref_string } });
    #[cfg(gpu)]
    let extender: &dyn ext::Extender = match &gpu_ext {
        Some(g) => g,
        None => &cpu_ext,
    };
    #[cfg(not(gpu))]
    let extender: &dyn ext::Extender = &cpu_ext;
    let mut dump = args.dump.as_ref().map(|p| std::io::BufWriter::with_capacity(8 << 20, std::fs::File::create(p).unwrap()));

    // Five stages, one batch each at a time (bounded channels): reading + decompressing; seeds, SA, chains;
    // extension; insert sizes + mate rescue; pairing + SAM + writing. Batch boundaries and read ids are fixed by the
    // reader, so the overlap changes nothing in the output.
    let chunk_size = args.chunk;
    let threads = args.threads;
    let (fmi, bns, opt, copt, lanes, ref_lens) = (&fmi, &bns, &opt, &copt, &lanes, &ref_lens);
    // pairing + SAM and mate rescue each have their own pool: in a shared rayon pool, workers busy with one stage's
    // parallel loops pick up another stage's injected jobs only when they run out of local work, so the stages
    // alternate instead of overlapping (rescue's job building waited for the previous batch's SAM text)
    let sam_pool_ = rayon::ThreadPoolBuilder::new().num_threads(threads).thread_name(|i| format!("sam-{i}")).build().unwrap();
    let rescue_pool_ = rayon::ThreadPoolBuilder::new().num_threads(threads).thread_name(|i| format!("rsc-{i}")).build().unwrap();
    let (sam_pool, sam_pool2) = (&rescue_pool_, &sam_pool_);
    std::thread::scope(|sc| {
        // every message carries its lane: batches and read ids restart per lane, as separate runs would number them
        let (rtx, rrx) = std::sync::mpsc::sync_channel::<(usize, Vec<fastq::Seq>, Vec<Vec<u8>>, i64)>(slack("read"));
        let (atx, arx) = std::sync::mpsc::sync_channel::<(usize, Vec<fastq::Seq>, Vec<Vec<u8>>, Vec<Vec<ext::Reg>>, i64)>(slack("extend"));
        sc.spawn(move || {
            let mut sent = 0usize;
            'lanes: for (li, lane) in lanes.iter().enumerate() {
                let mut f1 = fastq::Fastq::open(&lane.r1);
                let mut f2 = lane.r2.as_ref().map(|p| fastq::Fastq::open(p));
                let mut n_processed: i64 = 0;
                let mut empty_batches = 0;
                loop {
                    let tr = std::time::Instant::now();
                    let batch = fastq::read_batch(&mut f1, f2.as_mut(), chunk_size);
                    if batch.is_empty() {
                        // bwa-mem2 reads through kt_pipeline's 2 workers (process(.., no_mt_io ? 1 : 2)), each of
                        // which stops at its first empty batch: the input ends at the second one (a batch can be
                        // empty before the end, cut short by a malformed first record)
                        empty_batches += 1;
                        if empty_batches == 2 {
                            break;
                        }
                        continue;
                    }
                    // extension jobs keep query positions in 16 bits; bwa-mem2 itself overflows int16 positions
                    // from 32,768 bp on: refuse such reads plainly instead of failing deep in a stage
                    if let Some(r) = batch.iter().find(|r| r.seq.len() > u16::MAX as usize) {
                        eprintln!("[align-rs] read {} is {} bp: reads over {} bp are not supported", String::from_utf8_lossy(&r.name), r.seq.len(), u16::MAX);
                        std::process::exit(1);
                    }
                    let seqs: Vec<Vec<u8>> = batch.iter().map(|s| s.seq.iter().map(|&c| fastq::nt4(c)).collect()).collect();
                    padd(4, tr);
                    let n = seqs.len() as i64;
                    tl("read", (n_processed / 666668) as usize, tr);
                    if rtx.send((li, batch, seqs, n_processed)).is_err() {
                        break 'lanes;
                    }
                    n_processed += n;
                    sent += 1;
                    if sent >= max_batches {
                        break 'lanes;
                    }
                }
            }
            prof::stage_done("read + gunzip");
        });
        // final pass of each lane's BAM (sort + write + index), one lane at a time, while the next lane aligns
        let (ftx, frx) = std::sync::mpsc::channel::<bam::BamSink>();
        sc.spawn(move || {
            for b in frx {
                let tf = std::time::Instant::now();
                let (n, path) = (b.n_records, b.path().to_string());
                b.finish(threads);
                eprintln!("[align-rs] sorted BAM + index: {path}: {n} records, final pass {:.1} s", tf.elapsed().as_secs_f64());
            }
            prof::stage_done("BAM final pass");
        });
        let (dropper, drop_rx) = std::sync::mpsc::channel::<Box<dyn Send>>();
        sc.spawn(move || {
            for x in drop_rx {
                let t = std::time::Instant::now();
                drop(x);
                padd(8, t);
                DROP_PENDING.fetch_sub(1, std::sync::atomic::Ordering::Relaxed);
            }
            prof::stage_done("dropper");
        });
        #[cfg(gpu)]
        let gpu_sam: Option<&gpu::Gpu> = gpu_dev.as_ref().filter(|_| std::env::var("ALIGN_RS_GPU_MATESW").map_or(true, |v| v != "0"));
        // mate-rescue stage: insert-size statistics and the batch's rescue alignments (GPU), so that they overlap
        // the previous batch's pairing + SAM text
        let dropper_r = dropper.clone();
        let (mtx, mrx) = std::sync::mpsc::sync_channel::<(usize, Vec<fastq::Seq>, Vec<Vec<u8>>, Vec<Vec<ext::Reg>>, [sam::Pestat; 4], Vec<Vec<sam::PreSw>>, i64)>(slack("rescue"));
        sc.spawn(move || sam_pool.install(|| {
            #[cfg_attr(not(gpu), allow(unused_variables))] // used by the GPU mate rescue
            let dropper = dropper_r;
            for (li, batch, seqs, regs, n_processed) in arx {
                let paired = lanes[li].r2.is_some();
                let tq = std::time::Instant::now();
                let bi = (n_processed / 666668) as usize;
                if let Some(d) = dump.as_mut() {
                    for (i, rv) in regs.iter().enumerate() {
                        for (ri, p) in rv.iter().enumerate() {
                            writeln!(d, "REG\t{}\t{ri}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:08x}", n_processed + i as i64, p.rid, p.rb, p.re, p.qb, p.qe, p.score, p.truesc, p.sub, p.csub, p.w, p.seedcov, p.seedlen0, p.n_comp, p.is_alt, p.frac_rep.to_bits()).unwrap();
                        }
                    }
                }
                let mut pes = [sam::Pestat::default(); 4];
                #[allow(unused_mut)]
                let mut pre: Vec<Vec<sam::PreSw>> = Vec::new();
                if paired {
                    pes = prof::cpu(prof::C_PESTAT, || sam::pestat(&sopt, bns.l_pac, &regs));
                    if std::env::var("ALIGN_RS_DEBUG").is_ok() {
                        for (d, p) in pes.iter().enumerate() {
                            eprintln!("[align-rs][PE] dir {d}: failed={} low={} high={} avg={:.2} std={:.2}", p.failed, p.low, p.high, p.avg, p.std);
                        }
                    }
                    // mate-rescue alignments precomputed on the GPU (mem_matesw_batch_pre), looked up by sam_pe
                    #[cfg(gpu)]
                    if let Some(g) = gpu_sam {
                        let jobs: Vec<Vec<(u8, u16, u8, sam::SwSpec)>> = regs
                            .par_chunks(2)
                            .enumerate()
                            .map(|(p, a)| prof::cpu(prof::C_RJOB, || {
                                let mut v = Vec::new();
                                sam::matesw_jobs(&sopt, &bns, &pes, a, [seqs[2 * p].len(), seqs[2 * p + 1].len()], [(2 * p) as u32, (2 * p + 1) as u32], &mut v);
                                v.retain(|x| x.3.xtra & ksw::KSW_XBYTE != 0);
                                v
                            }))
                            .collect();
                        let mut first = Vec::with_capacity(jobs.len() + 1);
                        let mut k = 0usize;
                        for v in &jobs {
                            first.push(k);
                            k += v.len();
                        }
                        let mut flat: Vec<gpu::MJob> = Vec::with_capacity(k);
                        let fp = SendPtr(flat.as_mut_ptr());
                        jobs.par_iter().zip(first.par_iter()).for_each(|(v, &f)| {
                            let fp = fp;
                            for (m, x) in v.iter().enumerate() {
                                // SAFETY: each pair writes its own range first..first+len of the k slots
                                unsafe { *fp.0.add(f + m) = gpu::MJob { read: x.3.read, qlen: seqs[x.3.read as usize].len() as i32, rev: x.3.rev as i32, tlen: x.3.tlen, rb: x.3.rb, xtra: x.3.xtra, pad: 0, boff: 0 } };
                            }
                        });
                        unsafe { flat.set_len(k) };
                        tl(" rescue.jobs", bi, tq);
                        let tg = std::time::Instant::now();
                        let res = g.matesw(&seqs, gpu::batch_key(li, n_processed), &flat);
                        tl(" rescue.gpu", bi, tg);
                        if prof::on() {
                            prof::count(prof::N_MSW, flat.len() as u64);
                            prof::count(prof::N_MCELLS, flat.iter().map(|j| j.qlen as u64 * j.tlen as u64).sum());
                            prof::count(prof::N_MTLEN, flat.iter().map(|j| j.tlen as u64).sum());
                            prof::count(prof::N_MPASS, res.iter().filter(|r| r.ok != 0 && r.score >= sopt.min_seed_len * sopt.e.a).count() as u64);
                            let nd = res.iter().filter(|r| r.ok == 0).count();
                            if nd > 0 {
                                eprintln!("[align-rs] GPU mate SW: {} jobs, {} not done (ok = 0)", res.len(), nd);
                            }
                        }
                        pre = jobs
                            .par_iter()
                            .zip(first.par_iter())
                            .map(|(v, &f)| prof::cpu(prof::C_RUNPACK, || {
                                v.iter()
                                    .zip(&res[f..f + v.len()])
                                    .filter(|(_, r)| r.ok != 0)
                                    .map(|(x, r)| (x.0, x.1, x.2, ksw::Kswr { score: r.score, te: r.te, qe: r.qe, score2: r.score2, te2: r.te2, tb: r.tb, qb: r.qb }))
                                    .collect()
                            }))
                            .collect();
{ let k = DROP_PENDING.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1; DROP_MAX.fetch_max(k, std::sync::atomic::Ordering::Relaxed); }
                        dropper.send(Box::new((jobs, flat, res))).unwrap();
                    }
                }
                tl("rescue", bi, tq);
                let tw = std::time::Instant::now();
                mtx.send((li, batch, seqs, regs, pes, pre, n_processed)).unwrap();
                tl("rescue->sam wait", bi, tw);
            }
            if let Some(mut d) = dump {
                d.flush().unwrap();
            }
            prof::stage_done("rescue");
        }));
        let dropper3 = dropper.clone();
        let writer = sc.spawn(move || sam_pool2.install(|| {
            let dropper = dropper3;
            // the lane being written, and the next lane not yet started (lanes without reads get empty BAMs)
            let mut sink: Option<bam::BamSink> = None;
            let mut next_lane = 0usize;
            let new_sink = |l: &Lane| bam::BamSink::new(l.out.as_ref().unwrap(), &l.header(bns), ref_lens, level);
            for (li, batch, seqs, mut regs, pes, pre, n_processed) in mrx {
                let ts = std::time::Instant::now();
                let (lane, sopt, paired) = (&lanes[li], &lanes[li].sopt, lanes[li].r2.is_some());
                if lane.bam() && li >= next_lane {
                    if let Some(b) = sink.take() {
                        ftx.send(b).unwrap();
                    }
                    for l in &lanes[next_lane..li] {
                        ftx.send(new_sink(l)).unwrap();
                    }
                    sink = Some(new_sink(lane));
                    next_lane = li + 1;
                }
                let reads: Vec<sam::Read> = batch.iter().zip(&seqs).map(|(b, s)| sam::Read { name: &b.name, seq: s, qual: &b.qual }).collect();
                let text: Vec<Vec<u8>> = if paired {
                    let pairs: Vec<[Vec<ext::Reg>; 2]> = regs.chunks_mut(2).map(|c| [std::mem::take(&mut c[0]), std::mem::take(&mut c[1])]).collect();
                    let pre = &pre;
                    pairs
                        .into_par_iter()
                        .enumerate()
                        .map(|(p, a)| {
                            let id = (n_processed >> 1) + p as i64;
                            let t4 = std::time::Instant::now();
                            let (s0, s1) = prof::cpu(prof::C_SAMPE, || sam::sam_pe(&sopt, &bns, &pes, id, [&reads[2 * p], &reads[2 * p + 1]], a, pre.get(p).map(|v| v.as_slice())));
                            padd(3, t4);
                            let mut s0 = s0;
                            s0.extend_from_slice(&s1);
                            s0
                        })
                        .collect()
                } else {
                    regs.into_par_iter().enumerate().map(|(i, a)| sam::sam_se(&sopt, &bns, n_processed + i as i64, &reads[i], a)).collect()
                };
                padd(6, ts);
                tl("sam", (n_processed / 666668) as usize, ts);
                let tw = std::time::Instant::now();
                match sink.as_mut() {
                    Some(b) => b.add_batch(&text),
                    None => {
                        for t in &text {
                            out.write_all(t).unwrap();
                        }
                    }
                }
                padd(7, tw);
                tl("write", (n_processed / 666668) as usize, tw);
                drop(reads);
                let n = seqs.len();
{ let k = DROP_PENDING.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1; DROP_MAX.fetch_max(k, std::sync::atomic::Ordering::Relaxed); }
                mbsz("SAM text", n_processed, vb(&text) + text.iter().map(|t| t.capacity()).sum::<usize>());
                mbsz("rescue pre", n_processed, vb(&pre) + pre.iter().map(vb).sum::<usize>());
                dropper.send(Box::new((text, batch, seqs, pre))).unwrap();
                eprintln!("[align-rs] batch: {} reads, {:.2}s", n, t0.elapsed().as_secs_f64());
            }
            out.flush().unwrap();
            if lanes[0].bam() {
                if let Some(b) = sink.take() {
                    ftx.send(b).unwrap();
                }
                for l in &lanes[next_lane..] {
                    ftx.send(new_sink(l)).unwrap();
                }
            }
            drop(ftx);
            prof::stage_done("sam + write");
        }));
        // seeds, SA positions, chains: its own stage, so that batch N's extension overlaps batch N+1's seeding
        let (ctx, crx) = std::sync::mpsc::sync_channel::<(usize, Vec<fastq::Seq>, Vec<Vec<u8>>, Vec<Vec<chain::Chain>>, i64, std::time::Instant)>(slack("chain"));
        // seeding and SA + chaining are stages of their own (each half of what used to be one stage)
        let (stx, srx) = std::sync::mpsc::sync_channel::<(usize, Vec<fastq::Seq>, Vec<Vec<u8>>, Vec<Vec<index::Smem>>, i64, std::time::Instant)>(slack("seed"));
        let dropper2b = dropper.clone();
        #[cfg(gpu)]
        let gpu_dev = &gpu_dev;
        sc.spawn(move || {
            for (li, batch, seqs, n_processed) in rrx {
        // seeds, chains, regions: 512-read chunks
        let chunks: Vec<(usize, usize)> = (0..seqs.len()).step_by(BATCH_SIZE).map(|st| (st, (st + BATCH_SIZE).min(seqs.len()))).collect();
        let tp = std::time::Instant::now();
        // seeds per read
        let t1 = std::time::Instant::now();
        #[allow(unused_mut)]
        let mut smems: Option<Vec<Vec<index::Smem>>> = None;
        #[cfg(gpu)]
        if let Some(g) = gpu_dev.as_ref().filter(|g| g.seeding && g.set_reads(0, &seqs, gpu::batch_key(li, n_processed))) {
            const CAP: usize = SMEM_CAP;
            let tg = std::time::Instant::now();
            let split_len = (opt.min_seed_len as f64 * opt.split_factor + 0.499) as i32;
            let (out, cnt) = g.smem(seqs.len(), CAP, opt.min_seed_len, split_len, opt.split_width, opt.max_mem_intv);
            tl(" seed.gpu", (n_processed / 666668) as usize, tg);
            prof::count(prof::N_SEEDREDO, cnt.iter().filter(|&&c| c < 0).count() as u64);
            // where each read's seeds start in the packed output
            let mut soff = Vec::with_capacity(cnt.len());
            let mut acc = 0usize;
            for &c in &cnt {
                soff.push(acc);
                acc += c.max(0) as usize;
            }
            smems = Some(
                (0..seqs.len())
                    .into_par_iter()
                    .map(|i| prof::cpu(prof::C_SEED, || {
                        if cnt[i] < 0 {
                            return smem::collect(&fmi, &opt, &seqs[i], (i % BATCH_SIZE) as u32);
                        }
                        let mut v: Vec<index::Smem> = out[soff[i]..soff[i] + cnt[i] as usize]
                            .iter()
                            .map(|x| index::Smem { rid: (i % BATCH_SIZE) as u32, m: x.m as u32, n: x.n as u32, k: x.k, l: 0, s: x.s })
                            .collect();
                        v.sort_by_key(|s| ((s.m as u64) << 32) | s.n as u64);
                        v
                    }))
                    .collect(),
            );
        }
        let smems: Vec<Vec<index::Smem>> = smems.unwrap_or_else(|| {
            chunks.par_iter().flat_map_iter(|&(st, ed)| prof::cpu(prof::C_SEED, || (st..ed).map(|i| smem::collect(&fmi, &opt, &seqs[i], (i - st) as u32)).collect::<Vec<_>>())).collect()
        });
        padd(0, t1);
        tl(" seed (all)", (n_processed / 666668) as usize, t1);
        tl("seed", (n_processed / 666668) as usize, tp);
        let tw = std::time::Instant::now();
        if stx.send((li, batch, seqs, smems, n_processed, tp)).is_err() {
            break;
        }
        tl("seed->chain wait", (n_processed / 666668) as usize, tw);
            }
            prof::stage_done("seed");
        });
        sc.spawn(move || {
            let dropper = dropper2b;
            for (li, batch, seqs, smems, n_processed, _tp) in srx {
        let chunks: Vec<(usize, usize)> = (0..seqs.len()).step_by(BATCH_SIZE).map(|st| (st, (st + BATCH_SIZE).min(seqs.len()))).collect();
        let tp = std::time::Instant::now();
        // SA positions, per read (mem_chain_seeds: no read of a 512-read chunk with a single seed gets chains)
        let t2 = std::time::Instant::now();
        let eligible: Vec<bool> = (0..seqs.len())
            .map(|i| {
                let (st, ed) = ((i / BATCH_SIZE) * BATCH_SIZE, ((i / BATCH_SIZE + 1) * BATCH_SIZE).min(seqs.len()));
                let _ = (st, ed);
                !smems[i].is_empty() && (seqs[i].len() as i32) >= opt.min_seed_len
            })
            .collect();
        let chunk_ok: Vec<bool> = chunks.iter().map(|&(st, ed)| (st..ed).map(|i| smems[i].len()).sum::<usize>() > 1).collect();
        // SA positions of all reads in one flat buffer, filled in parallel (starts[i]..starts[i+1]: read i's)
        let counts: Vec<usize> = (0..seqs.len())
            .into_par_iter()
            .map(|i| if eligible[i] && chunk_ok[i / BATCH_SIZE] { prof::cpu(prof::C_SAPOS, || chain::sa_count(&copt, &smems[i])) } else { 0 })
            .collect();
        let mut starts = Vec::with_capacity(seqs.len() + 1);
        let mut acc = 0usize;
        for &c in &counts {
            starts.push(acc);
            acc += c;
        }
        starts.push(acc);
        let mut flat: Vec<i64> = Vec::with_capacity(acc);
        {
            let base = SendPtr(flat.as_mut_ptr());
            (0..seqs.len()).into_par_iter().filter(|&i| counts[i] > 0).for_each(|i| prof::cpu(prof::C_SAPOS, || {
                let mut v = Vec::with_capacity(counts[i]);
                chain::sa_positions(&copt, &smems[i], &mut v);
                debug_assert_eq!(v.len(), counts[i]);
                // SAFETY: disjoint ranges of a buffer with capacity acc
                unsafe { std::ptr::copy_nonoverlapping(v.as_ptr(), base.get().add(starts[i]), v.len()) };
            }));
            unsafe { flat.set_len(acc) };
        }
        prof::count(prof::N_SAPOS, flat.len() as u64);
        if prof::fine() {
            let mut u = flat.clone();
            u.par_sort_unstable();
            u.dedup();
            prof::count(prof::N_SAUNIQ, u.len() as u64);
        }
        tl(" sa.positions", (n_processed / 666668) as usize, t2);
        let t2g = std::time::Instant::now();
        #[allow(unused_mut)]
        let mut coords: Option<Vec<i64>> = None;
        #[cfg(gpu)]
        if let Some(g) = gpu_dev.as_ref().filter(|g| g.sa) {
            coords = g.sa(&flat);
        }
        let coords: Vec<i64> = coords.unwrap_or_else(|| flat.par_chunks(4096).flat_map_iter(|c| prof::cpu(prof::C_SACPU, || c.iter().map(|&p| fmi.sa_entry(p)).collect::<Vec<_>>())).collect());
        tl(" sa.lookup", (n_processed / 666668) as usize, t2g);
        let t2c = std::time::Instant::now();
        let chains: Vec<Vec<chain::Chain>> = (0..seqs.len())
            .into_par_iter()
            .map(|i| prof::cpu(prof::C_CHAIN, || {
                if !(eligible[i] && chunk_ok[i / BATCH_SIZE]) {
                    return Vec::new();
                }
                let l_seq = seqs[i].len() as i32;
                let c = prof::cpu(prof::C_CSEEDS, || chain::chain_seeds(&bns, &copt, &smems[i], &coords[starts[i]..starts[i + 1]], l_seq, (i % BATCH_SIZE) as i32));
                let mut c = prof::cpu(prof::C_CFLT, || chain::chain_flt(&copt, c));
                chain::flt_chained_seeds(&copt, &sopt.e, &bns, &seqs[i], &mut c);
                c
            }))
            .collect();
        padd(1, t2);
        tl(" chain", (n_processed / 666668) as usize, t2c);
{ let k = DROP_PENDING.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1; DROP_MAX.fetch_max(k, std::sync::atomic::Ordering::Relaxed); }
        mbsz("smems", n_processed, vb(&smems) + smems.iter().map(vb).sum::<usize>());
        mbsz("SA positions + coords", n_processed, vb(&flat) + vb(&coords));
        dropper.send(Box::new((smems, flat, coords))).unwrap();
                tl("sa+chain", (n_processed / 666668) as usize, tp);
                if ctx.send((li, batch, seqs, chains, n_processed, tp)).is_err() {
                    break;
                }
            }
            prof::stage_done("SA + chain");
        });
        for (li, batch, seqs, chains, n_processed, tp) in crx {
            let t3 = std::time::Instant::now();
            let regs = ext::extend_batch(&sopt.e, &bns, extender, &chains, &seqs, batch_key(li, n_processed));
            mbsz("fastq records", n_processed, vb(&batch) + batch.iter().map(|x| vb(&x.name) + vb(&x.comment) + vb(&x.seq) + vb(&x.qual)).sum::<usize>());
            mbsz("seqs", n_processed, vb(&seqs) + seqs.iter().map(vb).sum::<usize>());
            mbsz("chains", n_processed, vb(&chains) + chains.iter().map(|c| vb(c) + c.iter().map(|x| vb(&x.seeds)).sum::<usize>()).sum::<usize>());
            mbsz("regs", n_processed, vb(&regs) + regs.iter().map(vb).sum::<usize>());
            padd(2, t3);
            tl("extend", (n_processed / 666668) as usize, t3);
{ let k = DROP_PENDING.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1; DROP_MAX.fetch_max(k, std::sync::atomic::Ordering::Relaxed); }
            dropper.send(Box::new(chains)).unwrap();
            padd(5, tp);
            let tq = std::time::Instant::now();
            atx.send((li, batch, seqs, regs, n_processed)).unwrap();
            tl("extend->sam wait", (n_processed / 666668) as usize, tq);
        }
        drop(atx);
        prof::stage_done("extend (main thread)");
        drop(dropper);
        writer.join().unwrap();
    });
    if tuning {
        #[cfg(gpu)]
        {
            let Some(g) = gpu_dev.as_ref() else {
                eprintln!("[align-rs tune] no usable GPU: nothing to tune");
                std::process::exit(1);
            };
            let sp = tune::SeedParams { cap: SMEM_CAP, min_seed_len: opt.min_seed_len, split_len: (opt.min_seed_len as f64 * opt.split_factor + 0.499) as i32, split_width: opt.split_width, max_mem_intv: opt.max_mem_intv };
            if let Err(e) = tune::run(g, &sp) {
                eprintln!("[align-rs tune] FAILED: {e}");
                std::process::exit(1);
            }
        }
        #[cfg(not(gpu))]
        {
            eprintln!("[align-rs tune] built without CUDA: nothing to tune");
            std::process::exit(1);
        }
    }
    if std::env::var("ALIGN_RS_PROFILE").is_ok() {
        for (n, c) in PROF_NAMES.iter().zip(&PROF) {
            eprintln!("[align-rs] {n}: {:.1} s", c.load(std::sync::atomic::Ordering::Relaxed) as f64 * 1e-9);
        }
        let e: Vec<f64> = ext::EPROF.iter().map(|c| c.load(std::sync::atomic::Ordering::Relaxed) as f64).collect();
        eprintln!("[align-rs] extension: setup {:.1} s, job build {:.1} s, run {:.1} s ({} jobs), apply {:.1} s, replay+dedup {:.1} s", e[0] * 1e-9, e[1] * 1e-9, e[2] * 1e-9, e[5] as u64, e[3] * 1e-9, e[4] * 1e-9);
        let sp: Vec<f64> = sam::SPROF.iter().map(|c| c.load(std::sync::atomic::Ordering::Relaxed) as f64 * 1e-9).collect();
        eprintln!("[align-rs] mate SW calls {} (u8 {}), cells {:.3e}", sam::KSW_CALLS[0].load(std::sync::atomic::Ordering::Relaxed), sam::KSW_CALLS[1].load(std::sync::atomic::Ordering::Relaxed), sam::KSW_CALLS[2].load(std::sync::atomic::Ordering::Relaxed) as f64);
        eprintln!("[align-rs] pairing + SAM split (thread s): mate SW {:.1}, CIGAR/NM/MD {:.1}, SAM text {:.1}", sp[0], sp[1], sp[3]);
        #[cfg(gpu)]
        if let Some(g) = &gpu_dev {
            g.report();
        }
        let mut tlv = TIMELINE.lock().unwrap().clone();
        tlv.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap());
        for (st, b, s0, s1) in &tlv {
            eprintln!("[align-rs] timeline {:>16} batch {:3} {:8.2} -> {:8.2} ({:6.2} s)", st, b, s0, s1, s1 - s0);
        }
        eprintln!("[align-rs] dropper backlog max {}", DROP_MAX.load(std::sync::atomic::Ordering::Relaxed));
        prof::report();
        eprintln!("[align-rs] total wall {:.1} s", t0.elapsed().as_secs_f64());
    }
}
