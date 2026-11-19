//! Safe wrapper over the CUDA PairHMM C ABI (cuda-pairhmm/native/cuda_pairhmm.cu).
//! Numerics are bit-identical to Intel GKL's AVX-512 PairHMM (float pass + double fallback),
//! i.e. GATK's default PairHMM.
use std::sync::Mutex;

#[cfg(not(no_cuda))]
extern "C" {
    fn cuda_pairhmm_init(use_double: i32) -> i32;
    fn cuda_pairhmm_compute(n_groups: i32, group_n_reads: *const i32, group_n_haps: *const i32,
        read_bases: *const *const u8, read_quals: *const *const u8, read_ins: *const *const u8,
        read_del: *const *const u8, read_gcp: *const *const u8, read_lens: *const i32,
        hap_bases: *const *const u8, hap_lens: *const i32, out: *mut f64) -> i32;
    fn cuda_pairhmm_report();
    fn cuda_pairhmm_max_hap_len() -> i32;
}

// Built without libcudapairhmm.so: no GPU (init() is false, so nothing else is called)
#[cfg(no_cuda)]
unsafe fn cuda_pairhmm_init(_use_double: i32) -> i32 {
    0
}
#[cfg(no_cuda)]
#[allow(clippy::too_many_arguments)]
unsafe fn cuda_pairhmm_compute(_: i32, _: *const i32, _: *const i32, _: *const *const u8, _: *const *const u8, _: *const *const u8,
    _: *const *const u8, _: *const *const u8, _: *const i32, _: *const *const u8, _: *const i32, _: *mut f64) -> i32 {
    unreachable!("CUDA PairHMM not built")
}
#[cfg(no_cuda)]
unsafe fn cuda_pairhmm_report() {}
#[cfg(no_cuda)]
unsafe fn cuda_pairhmm_max_hap_len() -> i32 {
    0
}

/// PairHMM inputs of one read (all slices have the read's length).
#[derive(Clone, Copy)]
pub struct HmmRead<'a> { pub bases: &'a [u8], pub quals: &'a [u8], pub ins: &'a [u8], pub del: &'a [u8], pub gcp: &'a [u8] }

/// One region: reads x haplotypes.
pub struct HmmGroup<'a> { pub reads: Vec<HmmRead<'a>>, pub haps: Vec<&'a [u8]> }

static LOCK: Mutex<bool> = Mutex::new(false);

/// Initializes the GPU engine (idempotent). Returns false if no CUDA device is usable.
pub fn init(use_double: bool) -> bool {
    let mut g = LOCK.lock().unwrap();
    if !*g { *g = unsafe { cuda_pairhmm_init(use_double as i32) } == 1; }
    *g
}

/// Computes all groups in one GPU launch. Returns, per group, a row-major (read x haplotype)
/// matrix of log10 likelihoods (same layout GATK's VectorLoglessPairHMM uses).
/// Callers may run at once: each takes one of the engine's streams for the whole call
/// (VARIANT_RS_GPU_STREAMS, default 8), so one packs its batch while another's kernels run.
pub fn compute(groups: &[HmmGroup]) -> Vec<Vec<f64>> {
    assert!(*LOCK.lock().unwrap(), "hc_cudahmm::init() not called or failed");
    let (mut nr, mut nh) = (Vec::with_capacity(groups.len()), Vec::with_capacity(groups.len()));
    let (mut rb, mut rq, mut ri, mut rd, mut rc, mut rl) = (vec![], vec![], vec![], vec![], vec![], vec![]);
    let (mut hb, mut hl) = (vec![], vec![]);
    let mut total = 0usize;
    for g in groups {
        nr.push(g.reads.len() as i32); nh.push(g.haps.len() as i32);
        for r in &g.reads {
            debug_assert!(r.quals.len() == r.bases.len() && r.ins.len() == r.bases.len() && r.del.len() == r.bases.len() && r.gcp.len() == r.bases.len());
            rb.push(r.bases.as_ptr()); rq.push(r.quals.as_ptr()); ri.push(r.ins.as_ptr());
            rd.push(r.del.as_ptr()); rc.push(r.gcp.as_ptr()); rl.push(r.bases.len() as i32);
        }
        for h in &g.haps { hb.push(h.as_ptr()); hl.push(h.len() as i32); }
        total += g.reads.len() * g.haps.len();
    }
    let mut out = vec![0f64; total];
    if total > 0 {
        let ok = unsafe { cuda_pairhmm_compute(groups.len() as i32, nr.as_ptr(), nh.as_ptr(), rb.as_ptr(), rq.as_ptr(),
            ri.as_ptr(), rd.as_ptr(), rc.as_ptr(), rl.as_ptr(), hb.as_ptr(), hl.as_ptr(), out.as_mut_ptr()) };
        assert_eq!(ok, 1, "CUDA PairHMM failed");
    }
    let mut res = Vec::with_capacity(groups.len());
    let mut o = 0;
    for g in groups { let n = g.reads.len() * g.haps.len(); res.push(out[o..o + n].to_vec()); o += n; }
    res
}

/// Longest haplotype the GPU accepts (shared-memory limit for one warp per block: 4,223 bases with
/// 99 KB). Groups with a longer haplotype must go to the CPU PairHMM.
pub fn max_hap_len() -> usize {
    let _g = LOCK.lock().unwrap();
    unsafe { cuda_pairhmm_max_hap_len() }.max(0) as usize
}

/// Prints the engine's timing/throughput summary to stderr.
pub fn report() { let _g = LOCK.lock().unwrap(); unsafe { cuda_pairhmm_report() } }
