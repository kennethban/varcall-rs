//! Optional lightweight phase profiler (enabled with env VARIANT_RS_ASM_PROFILE=1). Totals are process-wide (all threads).
use std::sync::Mutex;
use std::time::Instant;

static ACC: Mutex<Vec<(&'static str, u64, u64)>> = Mutex::new(Vec::new());
static HITS: Mutex<Vec<(&'static str, u64)>> = Mutex::new(Vec::new());

pub fn enabled() -> bool {
    static EN: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *EN.get_or_init(|| std::env::var_os("VARIANT_RS_ASM_PROFILE").is_some())
}

pub struct Guard {
    name: &'static str,
    t0: Option<Instant>,
}

/// Times the enclosing scope (wall time of the calling thread; within one rayon task this is its CPU time
/// unless the thread is preempted).
#[inline]
pub fn scope(name: &'static str) -> Guard {
    Guard { name, t0: if enabled() { Some(Instant::now()) } else { None } }
}

impl Drop for Guard {
    fn drop(&mut self) {
        if let Some(t0) = self.t0 {
            let d = t0.elapsed().as_nanos() as u64;
            let mut a = ACC.lock().unwrap();
            if let Some(x) = a.iter_mut().find(|x| x.0 == self.name) {
                x.1 += d;
                x.2 += 1;
            } else {
                a.push((self.name, d, 1));
            }
        }
    }
}

pub fn report() -> String {
    let a = ACC.lock().unwrap();
    let mut s = String::new();
    for (n, d, c) in a.iter() {
        s.push_str(&format!("  {:<28} {:>9.3} s  ({} calls)\n", n, *d as f64 / 1e9, c));
    }
    s
}

/// Count an event (only when VARIANT_RS_ASM_PROFILE is set); used for code-path coverage reporting.
#[inline]
pub fn hit(name: &'static str) {
    add(name, 1);
}

/// Add `v` to a counter (only when VARIANT_RS_ASM_PROFILE is set).
#[inline]
pub fn add(name: &'static str, v: u64) {
    if !enabled() {
        return;
    }
    let mut h = HITS.lock().unwrap();
    if let Some(x) = h.iter_mut().find(|x| x.0 == name) {
        x.1 += v;
    } else {
        h.push((name, v));
    }
}

pub fn hits_report() -> String {
    let h = HITS.lock().unwrap();
    let mut s = String::new();
    for (n, c) in h.iter() {
        s.push_str(&format!("  {:<40} {}\n", n, c));
    }
    s
}
