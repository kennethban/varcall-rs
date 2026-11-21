//! `align-rs tune`: per-GPU settings measured on real work. The aligner runs on the first batch of the given reads
//! while the GPU calls' inputs are captured (the seeding batch, every extension round's jobs, the mate-SW jobs);
//! each candidate setting then replays them, one call at a time, timing kernels only. Every candidate's results
//! are compared with the first's: settings must not change them, and the tune stops if one does. The fastest
//! settings go to a cache file that `align-rs mem` loads at start (the environment still overrides it).

use crate::ext::{Job, Round};
use crate::gpu::{Gpu, MJob};
use std::sync::atomic::{AtomicBool, Ordering::Relaxed};
use std::sync::Mutex;

/// capturing (set by `align-rs tune` before the batch runs)
pub static ON: AtomicBool = AtomicBool::new(false);

#[derive(Default)]
pub struct Capture {
    pub seed_reads: Option<Vec<Vec<u8>>>,
    pub ext_reads: Option<Vec<Vec<u8>>>,
    pub ext: Vec<(Round, Vec<Job>)>,
    pub msw_reads: Option<Vec<Vec<u8>>>,
    pub msw: Vec<Vec<MJob>>,
}

pub static CAPTURE: Mutex<Option<Capture>> = Mutex::new(None);

pub fn capture(f: impl FnOnce(&mut Capture)) {
    if ON.load(Relaxed) {
        f(CAPTURE.lock().unwrap().get_or_insert_with(Capture::default));
    }
}

/// candidate values per setting (seeding lanes per read, extension blocks per SM, mate-SW register cap)
const SMEM_G: [i32; 5] = [1, 2, 4, 8, 16];
const EXT_BPSM: [i32; 11] = [0, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16];
const MSW_MINB: [i32; 3] = [12, 18, 24];
const REPS: usize = 3;
/// a value must beat the device default by this fraction to replace it
const TIE: f64 = 0.03;

/// Seeding parameters of the run (as main passes them to Gpu::smem)
pub struct SeedParams {
    pub cap: usize,
    pub min_seed_len: i32,
    pub split_len: i32,
    pub split_width: i64,
    pub max_mem_intv: i64,
}

/// The cache file: $ALIGN_RS_TUNE_CACHE, else $XDG_CACHE_HOME/align-rs/gpu-settings.tsv, else ~/.cache/align-rs/...
pub fn cache_path() -> Option<std::path::PathBuf> {
    if let Ok(p) = std::env::var("ALIGN_RS_TUNE_CACHE") {
        return Some(p.into());
    }
    let base = std::env::var("XDG_CACHE_HOME").map(std::path::PathBuf::from).or_else(|_| std::env::var("HOME").map(|h| std::path::Path::new(&h).join(".cache"))).ok()?;
    Some(base.join("align-rs").join("gpu-settings.tsv"))
}

/// One cached line: device name, sm_XY, driver, seeding lanes, extension blocks/SM, mate-SW cap, date
pub struct Cached {
    pub settings: [i32; 3],
    pub driver: i32,
    pub date: String,
}

pub fn load(name: &str, sm: &str) -> Option<Cached> {
    let text = std::fs::read_to_string(cache_path()?).ok()?;
    let mut found = None;
    for l in text.lines() {
        let c: Vec<&str> = l.split('\t').collect();
        if c.len() >= 7 && c[0] == name && c[1] == sm {
            let p = |i: usize| c[i].parse::<i32>().ok();
            if let (Some(d), Some(a), Some(b), Some(m)) = (p(2), p(3), p(4), p(5)) {
                found = Some(Cached { settings: [a, b, m], driver: d, date: c[6].to_string() });
            }
        }
    }
    found
}

fn save(name: &str, sm: &str, driver: i32, s: [i32; 3]) -> std::io::Result<std::path::PathBuf> {
    let p = cache_path().ok_or_else(|| std::io::Error::other("no cache location (set ALIGN_RS_TUNE_CACHE or HOME)"))?;
    if let Some(d) = p.parent() {
        std::fs::create_dir_all(d)?;
    }
    let old = std::fs::read_to_string(&p).unwrap_or_default();
    let mut out: String = old.lines().filter(|l| !l.starts_with(&format!("{name}\t{sm}\t"))).map(|l| format!("{l}\n")).collect();
    let date = std::process::Command::new("date").arg("+%Y-%m-%d").output().ok().map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string()).unwrap_or_default();
    out.push_str(&format!("{name}\t{sm}\t{driver}\t{}\t{}\t{}\t{date}\n", s[0], s[1], s[2]));
    std::fs::write(&p, out)?;
    Ok(p)
}

/// min over REPS of the kernel seconds of `kernel` (index into Gpu::kernel_times) that `f` takes; f's result
fn time_kernel<T>(g: &Gpu, kernel: usize, mut f: impl FnMut() -> T) -> (f64, T) {
    let mut best = f64::MAX;
    let mut last = None;
    for _ in 0..REPS {
        let t0 = g.kernel_times()[kernel];
        let r = f();
        best = best.min(g.kernel_times()[kernel] - t0);
        last = Some(r);
    }
    (best, last.unwrap())
}

fn pick<T: PartialEq>(g: &Gpu, what: &str, unit: &str, cands: &[i32], current: i32, kernel: usize, set: impl Fn(i32), run: impl Fn() -> T) -> Result<i32, String> {
    let mut times = Vec::new();
    let mut first: Option<T> = None;
    for &c in cands {
        set(c);
        let (t, r) = time_kernel(g, kernel, &run);
        match &first {
            None => first = Some(r),
            Some(f) if *f != r => return Err(format!("{what} {c}: results differ from {what} {}", cands[0])),
            _ => {}
        }
        times.push(t);
    }
    let fastest = (0..cands.len()).min_by(|&a, &b| times[a].partial_cmp(&times[b]).unwrap()).unwrap();
    let cur = cands.iter().position(|&c| c == current);
    // the device default stays unless another value is clearly faster (run-to-run noise is 1-3% here)
    let best = match cur {
        Some(c) if times[fastest] > times[c] * (1.0 - TIE) => c,
        _ => fastest,
    };
    eprintln!("[align-rs tune] {what} ({unit}): results identical for all {} values", cands.len());
    for (i, &c) in cands.iter().enumerate() {
        eprintln!(
            "[align-rs tune]   {c:>3}: {:7.3} s{}{}",
            times[i],
            if i == best { "  <- chosen" } else if i == fastest { "  (fastest, within the tie margin)" } else { "" },
            if Some(i) == cur { "  (device default)" } else { "" }
        );
    }
    set(cands[best]);
    Ok(cands[best])
}

/// Replays the captured batch under each candidate setting, keeps the fastest, writes the cache
pub fn run(g: &Gpu, sp: &SeedParams) -> Result<(), String> {
    ON.store(false, Relaxed);
    let cap = CAPTURE.lock().unwrap().take().ok_or("nothing captured (no GPU work in the batch?)")?;
    let defaults = g.settings();
    let mut chosen = defaults;
    let (name, sm, driver) = g.device();
    eprintln!("[align-rs tune] {name} ({sm}, driver {driver}); each value timed {REPS}x, kernels only, one call at a time");
    // seeding
    if let (true, Some(reads)) = (g.seeding, cap.seed_reads.as_ref()) {
        let key = u64::MAX - 1;
        chosen[0] = pick(
            g,
            "seeding lanes per read",
            &format!("{} reads", reads.len()),
            &SMEM_G,
            defaults[0],
            1,
            |c| g.set_settings([c, -1, -1]),
            || {
                assert!(g.set_reads(0, reads, key), "GPU read upload failed");
                g.smem(reads.len(), sp.cap, sp.min_seed_len, sp.split_len, sp.split_width, sp.max_mem_intv)
            },
        )?;
    }
    // extension
    if let Some(reads) = cap.ext_reads.as_ref() {
        let key = u64::MAX - 2;
        let n: usize = cap.ext.iter().map(|x| x.1.len()).sum();
        chosen[1] = pick(
            g,
            "extension blocks per SM",
            &format!("{} rounds, {n} jobs", cap.ext.len()),
            &EXT_BPSM,
            defaults[1],
            0,
            |c| g.set_settings([-1, c, -1]),
            || {
                assert!(g.set_reads(1, reads, key), "GPU read upload failed");
                cap.ext.iter().map(|(rd, jobs)| g.extend_raw(*rd, jobs)).collect::<Vec<_>>()
            },
        )?;
    }
    // mate SW: the register cap is the striped kernel's (ALIGN_RS_MSW=striped); the wavefront kernel's is fixed
    if g.msw_wave() {
        eprintln!("[align-rs tune] mate SW: the wavefront kernel (fixed register cap), not tuned; ALIGN_RS_MSW=striped tunes the striped kernel's cap");
    } else if let Some(reads) = cap.msw_reads.as_ref() {
        let key = u64::MAX - 3;
        let n: usize = cap.msw.iter().map(|x| x.len()).sum();
        chosen[2] = pick(
            g,
            "mate-SW blocks per SM (register cap)",
            &format!("{n} jobs"),
            &MSW_MINB,
            defaults[2],
            3,
            |c| g.set_settings([-1, -1, c]),
            || cap.msw.iter().map(|jobs| g.matesw(reads, key, jobs)).collect::<Vec<_>>(),
        )?;
    }
    let p = save(&name, &sm, driver, chosen).map_err(|e| format!("writing the cache: {e}"))?;
    eprintln!(
        "[align-rs tune] saved to {}: seeding {} lanes per read, extension {} blocks/SM, mate SW {} (device defaults {} / {} / {}); `align-rs mem` uses them on this GPU (ALIGN_RS_SMEM_G, ALIGN_RS_EXT_BPSM, ALIGN_RS_MSW_MINB override)",
        p.display(),
        chosen[0],
        chosen[1],
        chosen[2],
        defaults[0],
        defaults[1],
        defaults[2]
    );
    Ok(())
}
