//! GPU kernels (cuda/aln_ext.cu): seeding, suffix-array lookups and seed extension. What runs on the device
//! depends on its free memory, checked at start:
//!   * the 2-bit reference (.pac, 0.8 GB for GRCh38): extension;
//!   * + the FM-index checkpoints (6.4 GB): seeding;
//!   * + the sampled suffix array (4.0 GB): SA lookups.
//! Each kernel is an exact port of the CPU code (integer only); work a kernel does not take (reads longer than 256
//! bases, more seeds than the per-read room) is redone on the CPU.

use crate::ext::{CpuExtender, ExtOpt, Extender, Job, JobOut, Round};
use crate::index::{Bns, Fmi};

#[repr(C)]
#[derive(Clone, Copy, Default, PartialEq)]
pub struct OutSmem {
    pub k: i64,
    pub s: i64,
    pub m: i32,
    pub n: i32,
}

#[repr(C)]
#[derive(Clone, Copy, Default)]
pub struct MJob {
    pub read: u32,
    pub qlen: i32,
    pub rev: i32,
    pub tlen: i32,
    pub rb: i64,
    pub xtra: i32,
    pub pad: i32,
    pub boff: i64,
}

#[repr(C)]
#[derive(Clone, Copy, Default, PartialEq)]
pub struct MOut {
    pub score: i32,
    pub te: i32,
    pub qe: i32,
    pub score2: i32,
    pub te2: i32,
    pub tb: i32,
    pub qb: i32,
    pub ok: i32,
}

extern "C" {
    fn aln_gpu_matesw(jobs: *const MJob, n: i32, scratch: i64, out: *mut MOut) -> i32;
    fn aln_gpu_report3();
    fn aln_gpu_init(pac: *const u8, pac_bytes: i64, l_pac: i64, mat: *const i8, o_del: i32, e_del: i32, o_ins: i32, e_ins: i32, zdrop: i32) -> i32;
    fn aln_gpu_report();
    fn aln_gpu_report2();
    fn aln_gpu_set_reads(slot: i32, key: u64, reads: *const u8, reads_bytes: i64, off: *const i64, n_reads: i32) -> i32;
    fn aln_gpu_reads_cached(slot: i32, key: u64) -> i32;
    fn aln_gpu_extend(jobs: *const Job, n_jobs: i32, w: i32, end_bonus: i32, dir: i32, out: *mut JobOut) -> i32;
    fn aln_gpu_ext_take(out: *mut f64);
    fn aln_gpu_msw_take(out: *mut f64);
    fn aln_gpu_load_fmi(cp: *const u8, n_cp: i64, ms: *const i8, ls: *const u32, n_sa: i64, count: *const i64, sentinel: i64, with_sa: i32, reserve: i64) -> i32;
    fn aln_gpu_free_mem() -> i64;
    fn aln_gpu_smem(r0: i32, n: i32, cap: i32, min_seed_len: i32, split_len: i32, split_width: i64, max_mem_intv: i64, cnt: *mut i32) -> i64;
    fn aln_gpu_smem_fetch(out: *mut OutSmem, total: i64) -> i32;
    fn aln_gpu_sa(pos: *const i64, n: i64, out: *mut i64) -> i32;
    fn aln_gpu_set_settings(smem_g: i32, ext_bpsm: i32, msw_minb: i32) -> i32;
    fn aln_gpu_get_settings(out: *mut i32);
    fn aln_gpu_kernel_times(out: *mut f64);
    fn aln_gpu_device(name: *mut u8, cap: i32, major: *mut i32, minor: *mut i32, driver: *mut i32) -> i32;
}

/// device memory kept free for the per-batch buffers (reads, seeds, SA positions, extension and mate-SW jobs):
/// every kernel runs in chunks of the sizes below, so these stay ~1.2 GB whatever the batch
const RESERVE: i64 = 2 << 30;
const SMEM_CHUNK: usize = 50_000; // x cap 256 x 24 bytes: 307 MB of per-read rows on the device
const SA_CHUNK: usize = 8 << 20;
const EXT_CHUNK: usize = 4 << 20;
const MSW_CHUNK: usize = 64 << 10; // ~1,100 per-row slots (4 bytes) each: ~280 MB of scratch

/// A GPU call failed (device memory exhausted, typically: another process took it, or a large -K): the work it
/// had is redone on the CPU, which gives the same results. Said once per stage.
fn failed(stage: usize, what: &str) {
    static SAID: [std::sync::atomic::AtomicBool; 5] = [const { std::sync::atomic::AtomicBool::new(false) }; 5];
    if !SAID[stage].swap(true, std::sync::atomic::Ordering::Relaxed) {
        eprintln!("[align-rs] GPU {what} failed (see [align-rs-gpu] above): that work runs on the CPU (same results), here and whenever it fails again");
    }
}

/// A batch's identity for the device read buffers: its lane and first read
pub fn batch_key(lane: usize, n_processed: i64) -> u64 {
    ((lane as u64) << 40) | n_processed as u64
}

pub struct Gpu {
    pub seeding: bool,
    pub sa: bool,
    /// per stage (0 seeding, 1 extension, 2 mate SW): the batch whose reads it points at on the device
    reads_key: [std::sync::Mutex<u64>; 3],
    msw_lock: std::sync::Mutex<()>,
}

impl Gpu {
    /// None when there is no CUDA device
    pub fn new(o: &ExtOpt, bns: &Bns, fmi: &Fmi) -> Option<Gpu> {
        let pac = bns.pac_bytes();
        let rc = unsafe { aln_gpu_init(pac.as_ptr(), pac.len() as i64, bns.l_pac, o.mat.as_ptr(), o.o_del, o.e_del, o.o_ins, o.e_ins, o.zdrop) };
        if rc != 0 {
            return None;
        }
        let (cp, n_cp, ms, ls, n_sa) = fmi.raw();
        let mut g = Gpu { seeding: false, sa: false, reads_key: [std::sync::Mutex::new(u64::MAX), std::sync::Mutex::new(u64::MAX), std::sync::Mutex::new(u64::MAX)], msw_lock: std::sync::Mutex::new(()) };
        for with_sa in [1, 0] {
            let rc = unsafe { aln_gpu_load_fmi(cp as *const u8, n_cp as i64, ms, ls, n_sa as i64, fmi.count.as_ptr(), fmi.sentinel_index, with_sa, RESERVE) };
            assert!(rc >= 0, "GPU: loading the FM index failed");
            if rc == 0 {
                g.seeding = true;
                g.sa = with_sa == 1;
                break;
            }
        }
        // settings measured by `align-rs tune` on this GPU model, unless set in the environment
        let (name, sm, driver) = g.device();
        let tuning = crate::tune::ON.load(std::sync::atomic::Ordering::Relaxed);
        if let Some(c) = crate::tune::load(&name, &sm).filter(|_| !tuning) {
            let env = ["ALIGN_RS_SMEM_G", "ALIGN_RS_EXT_BPSM", "ALIGN_RS_MSW_MINB"];
            let mut s = [-1; 3];
            for i in 0..3 {
                if std::env::var(env[i]).is_err() {
                    s[i] = c.settings[i];
                }
            }
            g.set_settings(s);
            let cur = g.settings();
            eprintln!(
                "[align-rs] GPU settings from `align-rs tune` ({}, {}{}): seeding {} lanes per read, extension {} blocks/SM, mate SW cap {}",
                crate::tune::cache_path().map(|p| p.display().to_string()).unwrap_or_default(),
                c.date,
                if c.driver != driver { format!(", tuned with driver {}, now {driver}", c.driver) } else { String::new() },
                cur[0],
                cur[1],
                cur[2]
            );
        }
        let free = unsafe { aln_gpu_free_mem() };
        eprintln!(
            "[align-rs] GPU: extension{}{}; {:.1} GB device memory free",
            if g.seeding { ", seeding" } else { "" },
            if g.sa { ", SA lookups" } else { "" },
            free as f64 / (1u64 << 30) as f64
        );
        Some(g)
    }

    /// Points a stage (0 seeding, 1 extension, 2 mate SW) at batch `key`'s reads on the device, uploading them
    /// unless an earlier stage did. `key`: batch_key(lane, first read). False if the upload failed (the stage then
    /// runs the batch on the CPU).
    #[must_use]
    pub fn set_reads(&self, slot: usize, reads: &[Vec<u8>], key: u64) -> bool {
        if slot == 0 {
            crate::tune::capture(|c| {
                c.seed_reads.get_or_insert_with(|| reads.to_vec());
            });
        }
        let mut g = self.reads_key[slot].lock().unwrap();
        if *g == key || unsafe { aln_gpu_reads_cached(slot as i32, key) } == 1 {
            *g = key;
            return true;
        }
        let t0 = std::time::Instant::now();
        let total: usize = reads.iter().map(|r| r.len()).sum();
        let mut off = Vec::with_capacity(reads.len() + 1);
        let mut buf = Vec::with_capacity(total);
        for r in reads {
            off.push(buf.len() as i64);
            buf.extend_from_slice(r);
        }
        off.push(buf.len() as i64);
        let rc = unsafe { aln_gpu_set_reads(slot as i32, key, buf.as_ptr(), buf.len() as i64, off.as_ptr(), reads.len() as i32) };
        if rc != 0 {
            *g = u64::MAX;
            failed(4, "read upload");
            return false;
        }
        *g = key;
        crate::prof::count(crate::prof::N_GSETREADS, t0.elapsed().as_micros() as u64);
        true
    }

    /// Seeds of every read set by set_reads: (seeds packed in read order, count per read; -1 = not done, the CPU
    /// redoes it, and it has none in the packed seeds)
    pub fn smem(&self, n_reads: usize, cap: usize, min_seed_len: i32, split_len: i32, split_width: i64, max_mem_intv: i64) -> (Vec<OutSmem>, Vec<i32>) {
        let mut out: Vec<OutSmem> = Vec::with_capacity(n_reads * 16); // ~12 per 151-base read
        let mut cnt = vec![0i32; n_reads];
        let mut r0 = 0;
        while r0 < n_reads {
            let n = (n_reads - r0).min(SMEM_CHUNK);
            let t = std::time::Instant::now();
            let total = unsafe { aln_gpu_smem(r0 as i32, n as i32, cap as i32, min_seed_len, split_len, split_width, max_mem_intv, cnt[r0..].as_mut_ptr()) };
            if total < 0 {
                failed(0, "seeding");
                cnt[r0..r0 + n].fill(-1);
                r0 += n;
                continue;
            }
            crate::prof::count(crate::prof::N_GSMEM, t.elapsed().as_micros() as u64);
            let t = std::time::Instant::now();
            let at = out.len();
            out.resize(at + total as usize, OutSmem::default());
            crate::prof::count(crate::prof::N_GRESIZE, t.elapsed().as_micros() as u64);
            let t = std::time::Instant::now();
            if unsafe { aln_gpu_smem_fetch(out[at..].as_mut_ptr(), total) } != 0 {
                failed(0, "seeding");
                out.truncate(at);
                cnt[r0..r0 + n].fill(-1);
            }
            crate::prof::count(crate::prof::N_GFETCH, t.elapsed().as_micros() as u64);
            crate::prof::count(crate::prof::N_GSEEDS, total as u64);
            r0 += n;
        }
        (out, cnt)
    }

    /// None if the device failed (the caller looks them up on the CPU)
    pub fn sa(&self, pos: &[i64]) -> Option<Vec<i64>> {
        let mut out = vec![0i64; pos.len()];
        for (p, o) in pos.chunks(SA_CHUNK).zip(out.chunks_mut(SA_CHUNK)) {
            if unsafe { aln_gpu_sa(p.as_ptr(), p.len() as i64, o.as_mut_ptr()) } != 0 {
                failed(1, "SA lookups");
                return None;
            }
        }
        Some(out)
    }

    /// Mate-rescue alignments (u8 jobs) of `reads` (the SAM stage's batch)
    pub fn matesw(&self, reads: &[Vec<u8>], key: u64, jobs: &[MJob]) -> Vec<MOut> {
        crate::tune::capture(|c| {
            c.msw_reads.get_or_insert_with(|| reads.to_vec());
            c.msw.push(jobs.to_vec());
        });
        let t_lock = std::time::Instant::now();
        let _g = self.msw_lock.lock().unwrap();
        let t_up = std::time::Instant::now();
        // ok = 0: not done, the CPU aligns it
        let mut out = vec![MOut::default(); jobs.len()];
        if !self.set_reads(2, reads, key) {
            return out;
        }
        crate::prof::tr("gmsw0", t_lock, std::time::Instant::now(), format_args!("{key}\t{}\t{}", jobs.len(), (t_up - t_lock).as_micros()));
        unsafe {
            for (ci, (j, o)) in jobs.chunks(MSW_CHUNK).zip(out.chunks_mut(MSW_CHUNK)).enumerate() {
                let t_prep = std::time::Instant::now();
                let mut jj = j.to_vec();
                let mut off = 0i64;
                for x in jj.iter_mut() {
                    x.boff = off;
                    off += x.tlen as i64 + 1;
                }
                let t_call = std::time::Instant::now();
                let rc = aln_gpu_matesw(jj.as_ptr(), jj.len() as i32, off, o.as_mut_ptr());
                if crate::prof::tr_on() {
                    // CPU preparation, then the device call (upload / kernel / download by events)
                    let mut d = [0f64; 3];
                    aln_gpu_msw_take(d.as_mut_ptr());
                    let t_end = std::time::Instant::now();
                    let cells: i64 = jj.iter().map(|x| x.qlen as i64 * x.tlen as i64).sum();
                    crate::prof::tr("gmsw", t_prep, t_end, format_args!("{key}\t{ci}\t{}\t{cells}\t{}\t{}\t{:.0}\t{:.0}\t{:.0}", jj.len(), (t_call - t_prep).as_micros(), (t_end - t_call).as_micros(), d[0] * 1e6, d[1] * 1e6, d[2] * 1e6));
                }
                if rc != 0 {
                    failed(3, "mate SW");
                    o.fill(MOut::default());
                    continue;
                }
                // ALIGN_RS_DUMP_MSW=FILE: the first full chunk's reads, jobs and results, for cuda/bench_msw.cu
                static DUMPED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
                if jj.len() == MSW_CHUNK && !DUMPED.swap(true, std::sync::atomic::Ordering::Relaxed) {
                    if let Ok(p) = std::env::var("ALIGN_RS_DUMP_MSW") {
                        use std::io::Write;
                        let mut f = std::io::BufWriter::new(std::fs::File::create(p).unwrap());
                        let bytes = |x: &[u8]| x.to_vec();
                        let as_u8 = |p: *const u8, n: usize| std::slice::from_raw_parts(p, n);
                        let total: usize = reads.iter().map(|r| r.len()).sum();
                        f.write_all(&(reads.len() as u64).to_le_bytes()).unwrap();
                        f.write_all(&(total as u64).to_le_bytes()).unwrap();
                        let mut acc = 0i64;
                        for r in reads {
                            f.write_all(&acc.to_le_bytes()).unwrap();
                            acc += r.len() as i64;
                        }
                        f.write_all(&acc.to_le_bytes()).unwrap();
                        for r in reads {
                            f.write_all(&bytes(r)).unwrap();
                        }
                        f.write_all(&(jj.len() as u64).to_le_bytes()).unwrap();
                        f.write_all(&(off as u64).to_le_bytes()).unwrap();
                        f.write_all(as_u8(jj.as_ptr() as *const u8, jj.len() * std::mem::size_of::<MJob>())).unwrap();
                        f.write_all(as_u8(o.as_ptr() as *const u8, o.len() * std::mem::size_of::<MOut>())).unwrap();
                    }
                }
            }
        }
        out
    }

    /// (seeding lanes per read, extension blocks per SM, mate-SW register cap)
    pub fn settings(&self) -> [i32; 3] {
        let mut s = [0i32; 3];
        unsafe { aln_gpu_get_settings(s.as_mut_ptr()) };
        s
    }

    /// values < 0 are kept
    pub fn set_settings(&self, s: [i32; 3]) {
        assert_eq!(unsafe { aln_gpu_set_settings(s[0], s[1], s[2]) }, 0, "invalid GPU settings {s:?}");
    }

    /// kernel seconds so far: extension, seeding, SA, mate SW
    pub fn kernel_times(&self) -> [f64; 4] {
        let mut t = [0f64; 4];
        unsafe { aln_gpu_kernel_times(t.as_mut_ptr()) };
        t
    }

    /// (name, "sm_XY", driver version)
    pub fn device(&self) -> (String, String, i32) {
        let mut name = vec![0u8; 256];
        let (mut ma, mut mi, mut dr) = (0, 0, 0);
        unsafe { aln_gpu_device(name.as_mut_ptr(), 256, &mut ma, &mut mi, &mut dr) };
        let n = name.iter().position(|&c| c == 0).unwrap_or(0);
        (String::from_utf8_lossy(&name[..n]).to_string(), format!("sm_{ma}{mi}"), dr)
    }

    /// one round's extension jobs on the device (reads set in slot 1); jobs it leaves are INT_MIN
    pub fn extend_raw(&self, rd: Round, jobs: &[Job]) -> Vec<JobOut> {
        let mut out = vec![JobOut::default(); jobs.len()];
        for (j, o) in jobs.chunks(EXT_CHUNK).zip(out.chunks_mut(EXT_CHUNK)) {
            let rc = unsafe { aln_gpu_extend(j.as_ptr(), j.len() as i32, rd.w, rd.end_bonus, rd.dir, o.as_mut_ptr()) };
            assert!(rc >= 0, "GPU extension failed");
        }
        out
    }

    pub fn report(&self) {
        unsafe {
            aln_gpu_report();
            aln_gpu_report2();
            aln_gpu_report3();
        }
    }
}

pub struct GpuExtender<'a> {
    pub gpu: &'a Gpu,
    pub cpu: CpuExtender<'a>,
}

impl Extender for GpuExtender<'_> {
    fn run(&self, o: &ExtOpt, rd: Round, jobs: &[Job], reads: &[Vec<u8>], key: u64) -> Vec<JobOut> {
        use rayon::prelude::*;
        crate::tune::capture(|c| {
            c.ext_reads.get_or_insert_with(|| reads.to_vec());
            c.ext.push((rd, jobs.to_vec()));
        });
        // the device orders each chunk's jobs by size (cuda/aln_ext.cu) and returns results in job order
        let t_run = std::time::Instant::now();
        let mut out = vec![JobOut::default(); jobs.len()];
        let on_device = self.gpu.set_reads(1, reads, key);
        let t_dev = std::time::Instant::now();
        // jobs left for the CPU (the kernel's count, or a whole chunk when the device failed)
        let mut left = 0usize;
        for (j, o) in jobs.chunks(EXT_CHUNK).zip(out.chunks_mut(EXT_CHUNK)) {
            let rc = if on_device { unsafe { aln_gpu_extend(j.as_ptr(), j.len() as i32, rd.w, rd.end_bonus, rd.dir, o.as_mut_ptr()) } } else { -1 };
            if rc < 0 {
                if on_device {
                    failed(2, "extension");
                }
                o.iter_mut().for_each(|x| x.score = i32::MIN);
                left += o.len();
            } else {
                left += rc as usize;
            }
        }
        // jobs the kernel left (score INT_MIN): CPU. Looked for only when there are some: the scan over all the
        // round's results (a parallel loop) waited ~80 ms a batch for the thread pool at 64 threads, finding none
        let t_redo = std::time::Instant::now();
        let redo: Vec<usize> = if left == 0 { Vec::new() } else { (0..jobs.len()).into_par_iter().filter(|&i| out[i].score == i32::MIN).collect() };
        assert_eq!(redo.len(), left, "GPU extension: jobs left for the CPU");
        crate::prof::count(crate::prof::N_REDO, redo.len() as u64);
        if !redo.is_empty() {
            let sub: Vec<Job> = redo.iter().map(|&i| jobs[i]).collect();
            let r = self.cpu.run(o, rd, &sub, reads, key);
            for (k, &i) in redo.iter().enumerate() {
                out[i] = r[k];
            }
        }
        if crate::prof::tr_on() {
            // reads upload, device call (and in it, upload / sort / kernel / download by events), CPU redo
            let mut d = [0f64; 4];
            unsafe { aln_gpu_ext_take(d.as_mut_ptr()) };
            let us = |a: std::time::Instant, b: std::time::Instant| (b - a).as_micros();
            let t_end = std::time::Instant::now();
            crate::prof::tr("gext", t_run, t_end, format_args!("{key}\t{}\t{}\t{}\t{:.0}\t{:.0}\t{:.0}\t{:.0}\t{}\t{}", jobs.len(), us(t_run, t_dev), us(t_dev, t_redo), d[0] * 1e6, d[1] * 1e6, d[2] * 1e6, d[3] * 1e6, redo.len(), us(t_redo, t_end)));
        }
        out
    }
}
