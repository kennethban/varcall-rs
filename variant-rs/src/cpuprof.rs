//! `VARIANT_RS_TRACE_CPU=1`: thread CPU time (CLOCK_THREAD_CPUTIME_ID) summed per HaplotypeCaller back-end stage. Each
//! sample covers one task on one thread (a region on a rayon worker, a batch on the genotyping thread), so the sums
//! are CPU seconds, unlike the stage wall times.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

#[repr(C)]
struct Timespec {
    tv_sec: i64,
    tv_nsec: i64,
}
extern "C" {
    fn clock_gettime(clk: i32, ts: *mut Timespec) -> i32;
}
const CLOCK_THREAD_CPUTIME_ID: i32 = 3;

/// CPU time of the calling thread, in nanoseconds.
pub fn thread_cpu_ns() -> u64 {
    let mut ts = Timespec { tv_sec: 0, tv_nsec: 0 };
    // SAFETY: valid pointer to a timespec
    unsafe { clock_gettime(CLOCK_THREAD_CPUTIME_ID, &mut ts) };
    ts.tv_sec as u64 * 1_000_000_000 + ts.tv_nsec as u64
}

pub fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("VARIANT_RS_TRACE_CPU").is_some())
}

pub const N: usize = 12;
pub const NAMES: [&str; N] = [
    "finalizeRegion (read conversion + preprocessing)",
    "assembly",
    "prepare (trimming, genotyping reads, PairHMM inputs)",
    "PairHMM call (GPU: pack + wait on the calling thread; CPU: compute)",
    "finish (likelihood post-processing)",
    "realignment to best haplotype",
    "VCF genotyping + annotation + formatting (genotyping thread)",
    "GVCF genotyping + reference confidence (genotyping thread)",
    "GVCF writer thread (banding, formatting, deflate)",
    "batch drop / hand-off on compute threads (compute_batch outside the stages)",
    "-",
    "-",
];
pub static CPU: [AtomicU64; N] = [const { AtomicU64::new(0) }; N];
pub static COUNT: [AtomicU64; N] = [const { AtomicU64::new(0) }; N];

/// Runs `f`, adding its thread CPU time to stage `i` when profiling is on.
#[inline]
pub fn measure<T>(i: usize, f: impl FnOnce() -> T) -> T {
    if !on() {
        return f();
    }
    let t = thread_cpu_ns();
    let r = f();
    CPU[i].fetch_add(thread_cpu_ns() - t, Ordering::Relaxed);
    COUNT[i].fetch_add(1, Ordering::Relaxed);
    r
}

pub fn report() {
    if !on() {
        return;
    }
    let tot: u64 = CPU.iter().map(|c| c.load(Ordering::Relaxed)).sum();
    eprintln!("back-end CPU by stage (thread CPU seconds, summed over tasks):");
    for i in 0..N {
        let c = CPU[i].load(Ordering::Relaxed);
        if c == 0 {
            continue;
        }
        eprintln!("  cpu {:>9.1} s {:>5.1}%  tasks {:>9}  {}", c as f64 / 1e9, 100.0 * c as f64 / tot as f64, COUNT[i].load(Ordering::Relaxed), NAMES[i]);
    }
    eprintln!("  cpu {:>9.1} s total", tot as f64 / 1e9);
    if std::env::var_os("VARIANT_RS_ASM_PROFILE").is_some() {
        eprintln!("assembly phases (VARIANT_RS_ASM_PROFILE, summed over threads; nested scopes overlap):\n{}assembly events:\n{}", hc_assembly::prof::report(), hc_assembly::prof::hits_report());
    }
}

/// `VARIANT_RS_TRACE_CPU=1` for MarkDuplicates and BaseRecalibrator: thread CPU time by step, summed over threads.
/// A [`Laps`] charges the time since its last lap to a step; its sums go into [`LAP`] when it is dropped.
pub const NL: usize = 12;
pub static LAP: [AtomicU64; NL] = [const { AtomicU64::new(0) }; NL];

pub struct Laps {
    on: bool,
    last: u64,
    ns: [u64; NL],
}

impl Laps {
    pub fn new() -> Laps {
        let on = on();
        Laps { on, last: if on { thread_cpu_ns() } else { 0 }, ns: [0; NL] }
    }
    /// Charges the thread CPU time since the last lap to step `i`.
    #[inline]
    pub fn lap(&mut self, i: usize) {
        if self.on {
            let t = thread_cpu_ns();
            self.ns[i] += t - self.last;
            self.last = t;
        }
    }
}

impl Drop for Laps {
    fn drop(&mut self) {
        for (a, &n) in LAP.iter().zip(&self.ns) {
            if n > 0 {
                a.fetch_add(n, Ordering::Relaxed);
            }
        }
    }
}

/// The steps' CPU seconds (`names`: step index and name), when profiling is on.
pub fn report_laps(title: &str, names: &[(usize, &str)]) {
    if !on() {
        return;
    }
    let tot: u64 = names.iter().map(|&(i, _)| LAP[i].load(Ordering::Relaxed)).sum();
    eprintln!("{title} CPU by step (thread CPU seconds, summed over threads):");
    for &(i, n) in names {
        let c = LAP[i].load(Ordering::Relaxed);
        eprintln!("  cpu {:>8.1} s {:>5.1}%  {n}", c as f64 / 1e9, 100.0 * c as f64 / tot.max(1) as f64);
    }
    eprintln!("  cpu {:>8.1} s total", tot as f64 / 1e9);
}
