//! ALIGN_RS_PROFILE accounting (off otherwise): thread CPU time per piece of work and event counters, in per-thread
//! slots (no shared cache line between threads), summed at the end. Thread CPU time, unlike wall time inside a
//! parallel loop, does not grow when other stages' threads take the cores.

// Some slots are filled only by the GPU code.
#![cfg_attr(not(gpu), allow(dead_code))]

use std::sync::atomic::{AtomicU64, Ordering::Relaxed};
use std::sync::{Arc, Mutex, OnceLock};

#[repr(C)]
struct Ts {
    s: i64,
    ns: i64,
}
extern "C" {
    fn clock_gettime(clk: i32, ts: *mut Ts) -> i32;
}
const CLOCK_THREAD_CPUTIME_ID: i32 = 3;

/// this thread's CPU time (ns)
pub fn tcpu() -> u64 {
    let mut t = Ts { s: 0, ns: 0 };
    unsafe { clock_gettime(CLOCK_THREAD_CPUTIME_ID, &mut t) };
    t.s as u64 * 1_000_000_000 + t.ns as u64
}

pub fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var("ALIGN_RS_PROFILE").is_ok())
}

/// ALIGN_RS_PROFILE=2: also the thread CPU time per piece of work (a system call per timed item under WSL, which
/// slows the timed stages) and the SA position repeat count
pub fn fine() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var("ALIGN_RS_PROFILE").is_ok_and(|v| v == "2"))
}

/// CPU slots
pub const CPU_NAMES: [&str; 26] = [
    "seed: GPU output unpack / CPU seeding",
    "SA positions",
    "SA lookups (CPU)",
    "chaining",
    "ext: setup",
    "ext: job build",
    "ext: apply",
    "ext: replay + dedup",
    "rescue: pestat",
    "rescue: job build",
    "rescue: unpack",
    "sam_pe (all of it)",
    "  sam_pe: CIGAR/NM/MD (gen_cigar2)",
    "  sam_pe: SAM text",
    "BAM: parse (sam_parse1)",
    "BAM: deflate",
    "  sam_pe: mate rescue (CPU part)",
    "  sam_pe: mark_primary_se",
    "  sam_pe: pair",
    "  sam_pe: XA (gen_alt, incl. its CIGARs)",
    "  sam_pe: primary reg2aln (incl. CIGAR)",
    "  chaining: chain_seeds",
    "  chaining: chain_flt",
    "  replay: skip test",
    "  replay: sort_dedup_patch",
    "  setup: fetch_window",
];
pub const C_SEED: usize = 0;
pub const C_SAPOS: usize = 1;
pub const C_SACPU: usize = 2;
pub const C_CHAIN: usize = 3;
pub const C_ESETUP: usize = 4;
#[allow(dead_code)] // the slot is no longer timed
pub const C_EJOB: usize = 5;
pub const C_EAPPLY: usize = 6;
pub const C_EREPLAY: usize = 7;
pub const C_PESTAT: usize = 8;
pub const C_RJOB: usize = 9;
pub const C_RUNPACK: usize = 10;
pub const C_SAMPE: usize = 11;
pub const C_CIGAR: usize = 12;
pub const C_SAMTEXT: usize = 13;
pub const C_BPARSE: usize = 14;
pub const C_BDEFL: usize = 15;
pub const C_PRESCUE: usize = 16;
pub const C_PMARK: usize = 17;
pub const C_PPAIR: usize = 18;
pub const C_PXA: usize = 19;
pub const C_PALN: usize = 20;
pub const C_CSEEDS: usize = 21;
pub const C_CFLT: usize = 22;
pub const C_RSKIP: usize = 23;
pub const C_RDEDUP: usize = 24;
pub const C_FETCHW: usize = 25;

/// counters
pub const N_NAMES: [&str; 31] = [
    "ext: regions (seeds) extended",
    "ext: regions kept by the replay",
    "ext: jobs run",
    "ext: jobs of kept regions",
    "ext: GPU jobs redone on the CPU",
    "ext: reads with 0 kept regions",
    "ext: reads with 1",
    "ext: reads with 2",
    "ext: reads with 3",
    "ext: reads with 4",
    "ext: reads with 5-8",
    "ext: reads with 9-16",
    "ext: reads with 17+",
    "SA: positions looked up",
    "SA: distinct positions (per batch)",
    "mate SW: jobs",
    "mate SW: cells (qlen x tlen)",
    "mate SW: jobs with score >= min_seed_len",
    "mate SW: sum of tlen",
    "CIGAR: gen_cigar2 calls",
    "CIGAR: of them banded DP (not the ungapped shortcut)",
    "CIGAR: DP cells (qlen x band)",
    "XA: alternative hits written",
    "seed: reads redone on the CPU",
    "region sorts",
    "  of them with unequal ties (ks_introsort)",
    "GPU seeding: aln_gpu_smem calls (us)",
    "GPU seeding: fetch (us)",
    "GPU seeding: seeds",
    "GPU set_reads, all slots (us)",
    "GPU seeding: output Vec growth (us)",
];
pub const N_REG: usize = 0;
pub const N_KEPT: usize = 1;
pub const N_JOBS: usize = 2;
pub const N_KJOBS: usize = 3;
pub const N_REDO: usize = 4;
pub const N_HIST: usize = 5;
pub const N_SAPOS: usize = 13;
pub const N_SAUNIQ: usize = 14;
pub const N_MSW: usize = 15;
pub const N_MCELLS: usize = 16;
pub const N_MPASS: usize = 17;
pub const N_MTLEN: usize = 18;
pub const N_CIGAR: usize = 19;
pub const N_CIGDP: usize = 20;
pub const N_CIGCELLS: usize = 21;
pub const N_XA: usize = 22;
pub const N_SEEDREDO: usize = 23;
pub const N_SORTS: usize = 24;
pub const N_SORTSLOW: usize = 25;
pub const N_GSMEM: usize = 26;
pub const N_GFETCH: usize = 27;
pub const N_GSEEDS: usize = 28;
pub const N_GSETREADS: usize = 29;
pub const N_GRESIZE: usize = 30;

const NS: usize = CPU_NAMES.len() + N_NAMES.len();
type Slots = [AtomicU64; NS];

static ALL: Mutex<Vec<Arc<Slots>>> = Mutex::new(Vec::new());
thread_local! {
    static MINE: Arc<Slots> = {
        let s: Arc<Slots> = Arc::new([const { AtomicU64::new(0) }; NS]);
        ALL.lock().unwrap().push(s.clone());
        s
    };
}

fn add(i: usize, v: u64) {
    MINE.with(|s| s[i].fetch_add(v, Relaxed));
}

/// runs f, adding this thread's CPU time to slot i (profiling on)
#[inline]
pub fn cpu<T>(i: usize, f: impl FnOnce() -> T) -> T {
    if !fine() {
        return f();
    }
    let t = tcpu();
    let r = f();
    add(i, tcpu().saturating_sub(t));
    r
}

pub fn add_cpu(i: usize, ns: u64) {
    add(i, ns);
}

#[inline]
pub fn count(i: usize, v: u64) {
    if on() {
        add(CPU_NAMES.len() + i, v);
    }
}

/// thread CPU time of the named stage threads (recorded by each before it ends)
static STAGES: Mutex<Vec<(String, f64)>> = Mutex::new(Vec::new());
pub fn stage_done(name: &str) {
    if on() {
        STAGES.lock().unwrap().push((name.to_string(), tcpu() as f64 * 1e-9));
    }
}

/// CPU time (s) of this process's live threads, by name with trailing digits removed (pool threads)
fn task_cpu() -> Vec<(String, f64)> {
    let mut m: std::collections::BTreeMap<String, f64> = Default::default();
    if let Ok(rd) = std::fs::read_dir("/proc/self/task") {
        for e in rd.flatten() {
            let p = e.path();
            let comm = std::fs::read_to_string(p.join("comm")).unwrap_or_default().trim().to_string();
            let stat = std::fs::read_to_string(p.join("stat")).unwrap_or_default();
            // fields after the ")": state is field 3; utime / stime are fields 14 / 15
            if let Some(rest) = stat.rsplit_once(')').map(|x| x.1) {
                let f: Vec<&str> = rest.split_whitespace().collect();
                if f.len() > 13 {
                    let t = (f[11].parse::<f64>().unwrap_or(0.0) + f[12].parse::<f64>().unwrap_or(0.0)) / 100.0;
                    let key = comm.trim_end_matches(|c: char| c.is_ascii_digit()).to_string();
                    *m.entry(key).or_default() += t;
                }
            }
        }
    }
    m.into_iter().collect()
}

pub fn report() {
    if !on() {
        return;
    }
    let mut tot = [0u64; NS];
    for s in ALL.lock().unwrap().iter() {
        for i in 0..NS {
            tot[i] += s[i].load(Relaxed);
        }
    }
    eprintln!("[align-rs] thread CPU time (s), by piece of work:");
    for (i, n) in CPU_NAMES.iter().enumerate() {
        eprintln!("[align-rs]   {:<40} {:9.1}", n, tot[i] as f64 * 1e-9);
    }
    eprintln!("[align-rs] thread CPU time (s), by thread:");
    for (n, t) in STAGES.lock().unwrap().iter() {
        eprintln!("[align-rs]   stage thread {:<27} {:9.1}", n, t);
    }
    for (n, t) in task_cpu() {
        eprintln!("[align-rs]   live threads {:<27} {:9.1}", n, t);
    }
    eprintln!("[align-rs] counters:");
    for (i, n) in N_NAMES.iter().enumerate() {
        eprintln!("[align-rs]   {:<52} {:>15}", n, tot[CPU_NAMES.len() + i]);
    }
}

/// ALIGN_RS_TRACE=FILE: one line per GPU call, per extension round and CPU step, and per stage run, for working
/// out where the device waits. Tab-separated: event, thread, start µs, end µs (since the run started), details.
static TRACE: OnceLock<Option<Mutex<std::io::BufWriter<std::fs::File>>>> = OnceLock::new();
static TRACE_T0: OnceLock<std::time::Instant> = OnceLock::new();
fn trace_file() -> Option<&'static Mutex<std::io::BufWriter<std::fs::File>>> {
    TRACE.get_or_init(|| std::env::var_os("ALIGN_RS_TRACE").map(|p| Mutex::new(std::io::BufWriter::new(std::fs::File::create(p).expect("ALIGN_RS_TRACE: cannot create the file"))))).as_ref()
}
pub fn tr_on() -> bool {
    trace_file().is_some()
}
/// the trace's time origin (the run's start)
pub fn tr_start(t0: std::time::Instant) {
    let _ = TRACE_T0.set(t0);
}
pub fn tr(ev: &str, t0: std::time::Instant, t1: std::time::Instant, detail: std::fmt::Arguments) {
    if let Some(f) = trace_file() {
        use std::io::Write;
        let z = *TRACE_T0.get_or_init(std::time::Instant::now);
        let us = |t: std::time::Instant| t.saturating_duration_since(z).as_micros();
        let th = std::thread::current();
        let line = format!("{ev}\t{}\t{}\t{}\t{detail}\n", th.name().map_or_else(|| format!("{:?}", th.id()), str::to_string), us(t0), us(t1));
        f.lock().unwrap().write_all(line.as_bytes()).unwrap();
    }
}
pub fn tr_flush() {
    if let Some(f) = trace_file() {
        use std::io::Write;
        f.lock().unwrap().flush().unwrap();
    }
}
