// align-rs GPU seed extension: bwa-mem2's scalarBandedSWA (= bwa's ksw_extend2 with end_bonus), one thread per job.
// Integer-only, the same recurrence as the CPU port (src/ext.rs ksw_extend2), so results are identical by
// construction. The reference lives on the device as bwa's 2-bit .pac (forward strand); the reverse strand is
// complemented on the fly, as bwa-mem2's .0123 holds it.
#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <stdint.h>
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <atomic>

// Run-time settings, chosen in aln_gpu_init from the device (gpu_settings) unless set in the environment. None
// changes any result: they only spread the same work differently.
static int g_smem_g = 4;     // ALIGN_RS_SMEM_G: seeding lanes per read (1: a thread per read)
static int g_ext_bpsm = 6;   // ALIGN_RS_EXT_BPSM: extension blocks per SM (0: a block per 128 jobs)
static int g_msw_minb = 18;  // ALIGN_RS_MSW_MINB: k_matesw's blocks per SM the registers are capped for (12, 18, 24)
static int g_msw_wave = 1;   // mate SW: k_matesw_wave (ALIGN_RS_MSW=striped: k_matesw)
static int g_n_sm = 0;
static int gpu_settings(int dev);
#define MAXQ 256  // longest query a job may have; longer jobs come back with score INT_MIN (CPU recomputes)

struct Job {  // = align-rs ext::Job (repr(C)); band, end bonus and direction are the round's (kernel arguments)
    int64_t rstart;
    uint32_t read, rlen;
    uint16_t qstart, qlen, h0, pad;
};
static_assert(sizeof(Job) == 24, "Job layout");
struct JobOut { int32_t score, qle, tle, gtle, gscore, max_off; };

__constant__ int8_t c_mat[25];
__constant__ int32_t c_par[5];  // o_del, e_del, o_ins, e_ins, zdrop

static uint8_t *d_pac = 0;
static int64_t g_l_pac = 0;
static uint8_t *d_reads = 0;  static size_t cap_reads = 0;
static int64_t *d_off = 0;    static size_t cap_off = 0;
static Job *d_jobs = 0;       static size_t cap_jobs = 0;
static double t_h2d = 0, t_kernel = 0, t_d2h = 0, t_smem = 0, t_sa = 0;
static int g_n_reads = 0;
static cudaEvent_t ev[5];
static JobOut *d_out = 0;     static size_t cap_out = 0;

#define CK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { fprintf(stderr, "[align-rs-gpu] %s: %s\n", #x, cudaGetErrorString(e_)); return -1; } } while (0)

__device__ __forceinline__ int ref_base(const uint8_t *pac, int64_t l_pac, int64_t c) {
    if (c < l_pac) return (pac[c >> 2] >> ((~c & 3) << 1)) & 3;
    int64_t f = (l_pac << 1) - 1 - c;
    return 3 - ((pac[f >> 2] >> ((~f & 3) << 1)) & 3);
}

// thread t runs job order[t] (jobs of similar size side by side, so that a warp's threads finish together) and
// writes its result in the job's own slot
// The jobs the ungapped shortcut took (k_extkey: *nug of them) are sorted last and skipped.
__global__ void k_extend(const Job *jobs, const uint32_t *order, int n, const int *nug, const uint8_t *reads, const int64_t *off, const uint8_t *pac, int64_t l_pac, JobOut *out, int *nleft, int jw, int jend_bonus, int jdir) {
  n -= *nug;
  // grid-stride: a fixed number of blocks per SM (aln_gpu_extend), so the threads' DP rows fit L1
  for (int t0 = blockIdx.x * blockDim.x + threadIdx.x; t0 < n; t0 += gridDim.x * blockDim.x) {
    const int t = order[t0];
    const Job jb = jobs[t];
    const int qlen = jb.qlen, tlen = (int)jb.rlen;
    JobOut r;
    if (qlen > MAXQ) { r.score = INT_MIN; r.qle = r.tle = r.gtle = r.gscore = r.max_off = 0; out[t] = r; atomicAdd(nleft, 1); continue; }
    const int o_del = c_par[0], e_del = c_par[1], o_ins = c_par[2], e_ins = c_par[3], zdrop = c_par[4];
    const int oe_del = o_del + e_del, oe_ins = o_ins + e_ins;
    // per column one word: H (bits 0-13), E (14-27) and the query base (28-30), so a cell is one local load
    // and one store (the row lives in local memory, which spills to L2: the kernel is bound by L2 traffic).
    // H and E stay below h0 + qlen * max(mat); jobs where that could reach 2^14 go to the CPU.
    uint32_t eh[MAXQ + 1];
#define EH_H(x) ((int)((x) & 0x3fff))
#define EH_E(x) ((int)(((x) >> 14) & 0x3fff))
#define EH_Q(x) ((int)((x) >> 28))
#define EH_QM 0xf0000000u
#define EH(h, e) ((uint32_t)(h) | ((uint32_t)(e) << 14))
    {
        int mx = 0;
        for (int i = 0; i < 25; ++i) mx = mx > c_mat[i] ? mx : c_mat[i];
        if (jb.h0 + qlen * mx >= (1 << 14)) { r.score = INT_MIN; r.qle = r.tle = r.gtle = r.gscore = r.max_off = 0; out[t] = r; atomicAdd(nleft, 1); continue; }
    }
    const uint8_t *rd = reads + off[jb.read];
    for (int i = 0; i < qlen; ++i) eh[i] = (uint32_t)rd[jb.qstart + (int64_t)i * jdir] << 28;
    eh[qlen] = 0;
    const int h0 = jb.h0, end_bonus = jend_bonus;
    int w = jw;
    {
        int h = h0;
        eh[0] |= EH(h0, 0);
        if (qlen >= 1) { h = h0 > oe_ins ? h0 - oe_ins : 0; eh[1] |= EH(h, 0); }
        for (int j = 2; j <= qlen && h > e_ins; ++j) { h -= e_ins; eh[j] |= EH(h, 0); }
    }
    int max = 0;
    for (int i = 0; i < 25; ++i) max = max > c_mat[i] ? max : c_mat[i];
    int max_ins = (int)((double)(qlen * max + end_bonus - o_ins) / e_ins + 1.);
    max_ins = max_ins > 1 ? max_ins : 1;
    w = w < max_ins ? w : max_ins;
    int max_del = (int)((double)(qlen * max + end_bonus - o_del) / e_del + 1.);
    max_del = max_del > 1 ? max_del : 1;
    w = w < max_del ? w : max_del;
    max = h0;
    int max_i = -1, max_j = -1, max_ie = -1, gscore = -1, max_off = 0;
    int beg = 0, end = qlen;
    for (int i = 0; i < tlen; ++i) {
        int f = 0, h1, m = 0, mj = -1;
        const int tb = ref_base(pac, l_pac, jb.rstart + (int64_t)i * jdir);
        const int8_t *mrow = c_mat + tb * 5;
        if (beg < i - w) beg = i - w;
        if (end > i + w + 1) end = i + w + 1;
        if (end > qlen) end = qlen;
        if (beg == 0) { h1 = h0 - (o_del + e_del * (i + 1)); if (h1 < 0) h1 = 0; } else h1 = 0;
        int j;
        for (j = beg; j < end; ++j) {
            uint32_t x = eh[j];
            int M = EH_H(x), e = EH_E(x), h, tt;
            M = M ? M + mrow[EH_Q(x)] : 0;
            h = M > e ? M : e;
            h = h > f ? h : f;
            mj = m > h ? mj : j;
            m = m > h ? m : h;
            tt = M - oe_del; tt = tt > 0 ? tt : 0;
            e -= e_del; e = e > tt ? e : tt;
            eh[j] = EH(h1, e) | (x & EH_QM);  // H(i,j-1) for the next row, E(i+1,j)
            h1 = h;
            tt = M - oe_ins; tt = tt > 0 ? tt : 0;
            f -= e_ins; f = f > tt ? f : tt;
        }
        eh[end] = EH(h1, 0) | (eh[end] & EH_QM);
        if (j == qlen) { max_ie = gscore > h1 ? max_ie : i; gscore = gscore > h1 ? gscore : h1; }
        if (m == 0) break;
        if (m > max) {
            max = m; max_i = i; max_j = mj;
            int d = mj - i; d = d < 0 ? -d : d;
            max_off = max_off > d ? max_off : d;
        } else if (zdrop > 0) {
            if (i - max_i > mj - max_j) { if (max - m - ((i - max_i) - (mj - max_j)) * e_del > zdrop) break; }
            else { if (max - m - ((mj - max_j) - (i - max_i)) * e_ins > zdrop) break; }
        }
        for (j = beg; j < end && (eh[j] & ~EH_QM) == 0; ++j);
        beg = j;
        for (j = end; j >= beg && (eh[j] & ~EH_QM) == 0; --j);
        end = j + 2 < qlen ? j + 2 : qlen;
    }
    r.score = max; r.qle = max_j + 1; r.tle = max_i + 1; r.gtle = max_ie + 1; r.gscore = gscore; r.max_off = max_off;
    out[t] = r;
  }
}

// ALIGN_RS_GPU_SERIAL=1 (profiling): one wrapper at a time on the device, so that each kernel's time is its own cost
#include <mutex>
static std::mutex g_serial_mu;
static bool g_serial() { static int s = -1; if (s < 0) s = getenv("ALIGN_RS_GPU_SERIAL") != nullptr; return s; }
#define SERIAL std::unique_lock<std::mutex> serial_lk(g_serial_mu, std::defer_lock); if (g_serial()) serial_lk.lock()
// kernel-only seconds per kind: 0 extension, 1 seeding, 2 SA, 3 mate SW
static double t_kern[4] = {0, 0, 0, 0};
static std::mutex g_tk_mu;
static void add_kern(int k, cudaEvent_t a, cudaEvent_t b) { float ms_; cudaEventElapsedTime(&ms_, a, b); std::lock_guard<std::mutex> g(g_tk_mu); t_kern[k] += ms_ * 1e-3; }

// timing events, created once per calling thread (creating and destroying them per call cost ~8 ms a
// destroy under WSL: 13.7 s of stage time on 1/6 of a lane)
struct ThreadEvents { cudaEvent_t e[4]; bool init = false; };
static thread_local ThreadEvents t_ev;
static cudaEvent_t *tev() {
    if (!t_ev.init) { for (int i = 0; i < 4; ++i) cudaEventCreate(&t_ev.e[i]); t_ev.init = true; }
    return t_ev.e;
}

// ALIGN_RS_GPU_FAIL_GROW=N (tests): every per-batch allocation after the first N fails, as when device memory runs out
static long long grow_fail_after() {
    static long long v = [] { const char *e = getenv("ALIGN_RS_GPU_FAIL_GROW"); return e ? atoll(e) : -1LL; }();
    return v;
}
static std::atomic<long long> g_grows{0};
template <typename T> static int grow(T **p, size_t *cap, size_t n) {
    if (n <= *cap) return 0;
    // a failed allocation leaves no buffer (not the freed one with its old capacity)
    if (*p) cudaFree(*p);
    *p = 0;
    *cap = 0;
    size_t c = n + n / 4 + 1024;
    long long lim = grow_fail_after();
    if (lim >= 0 && g_grows.fetch_add(1) >= lim) {
        fprintf(stderr, "[align-rs-gpu] cudaMalloc: out of memory (ALIGN_RS_GPU_FAIL_GROW)\n");
        return -1;
    }
    CK(cudaMalloc((void **)p, c * sizeof(T)));
    *cap = c;
    return 0;
}

extern "C" int aln_gpu_init(const uint8_t *pac, int64_t pac_bytes, int64_t l_pac, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int zdrop) {
    int n = 0;
    if (cudaGetDeviceCount(&n) != cudaSuccess || n == 0) return -1;
    int dev = 0;
    CK(cudaGetDevice(&dev));
    if (gpu_settings(dev)) return -2;  // no kernels for this device in the build
    // block instead of spinning while waiting for the device: four stage threads call in, and a spinning thread
    // holds a core the CPU stages need
    CK(cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync));
    CK(cudaMalloc((void **)&d_pac, pac_bytes));
    CK(cudaMemcpy(d_pac, pac, pac_bytes, cudaMemcpyHostToDevice));
    g_l_pac = l_pac;
    int32_t par[5] = {o_del, e_del, o_ins, e_ins, zdrop};
    // the mate-SW kernel computes scores from bwa_fill_scmat's form: match / mismatch, N (row and column 4) NS
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j) {
            int8_t want = (i > 3 || j > 3) ? mat[4] : (i == j ? mat[0] : mat[1]);
            if (mat[i * 5 + j] != want) { fprintf(stderr, "[align-rs-gpu] scoring matrix not of bwa's form\n"); return -1; }
        }
    CK(cudaMemcpyToSymbol(c_mat, mat, 25));
    CK(cudaMemcpyToSymbol(c_par, par, sizeof(par)));
    CK(cudaDeviceSetLimit(cudaLimitStackSize, 8192));
    for (int i = 0; i < 5; ++i) CK(cudaEventCreate(&ev[i]));
    return 0;
}

// A batch's reads (concatenated, 0..4) and their offsets (n_reads + 1 entries), per slot: 0 seeding,
// 1 extension, 2 mate SW. The pipeline stages run on different batches at once, so each has its own copy.
static uint8_t *d_rs[3] = {0, 0, 0};
static int64_t *d_os[3] = {0, 0, 0};
// Batches' reads on the device, keyed by batch (lane << 40 | first read), uploaded once and used by seeding,
// extension and mate SW in turn (each stage points at one: d_rs / d_os[stage]). A buffer no stage points at
// is replaced oldest first; a batch replaced too early is uploaded again.
#define RB_N 6
struct ReadBuf { uint64_t key; uint64_t age; uint8_t *rs; size_t cap_rs; int64_t *os; size_t cap_os; };
static ReadBuf g_rb[RB_N];
static int g_stage_rb[3] = {-1, -1, -1};
static uint64_t g_rb_clock = 0;
static std::mutex g_rb_mu;
// 1 when batch `key`'s reads are on the device (stage now points at them), 0 otherwise
extern "C" int aln_gpu_reads_cached(int stage, uint64_t key) {
    std::lock_guard<std::mutex> lk(g_rb_mu);
    for (int i = 0; i < RB_N; ++i)
        if (g_rb[i].rs && g_rb[i].key == key) {
            g_stage_rb[stage] = i; d_rs[stage] = g_rb[i].rs; d_os[stage] = g_rb[i].os;
            return 1;
        }
    return 0;
}
extern "C" int aln_gpu_set_reads(int slot, uint64_t key, const uint8_t *reads, int64_t reads_bytes, const int64_t *off, int n_reads) {
    SERIAL;
    std::lock_guard<std::mutex> lk(g_rb_mu);
    // the oldest buffer no other stage points at (unused ones first)
    int pick = -1;
    for (int i = 0; i < RB_N; ++i) {
        bool used = false;
        for (int s = 0; s < 3; ++s) used |= s != slot && g_stage_rb[s] == i;
        if (used) continue;
        if (pick < 0 || !g_rb[i].rs || (g_rb[pick].rs && g_rb[i].age < g_rb[pick].age)) pick = i;
    }
    if (pick < 0) return -1;
    ReadBuf &b = g_rb[pick];
    b.key = ~0ull;
    if (grow(&b.rs, &b.cap_rs, (size_t)reads_bytes + 1) || grow(&b.os, &b.cap_os, (size_t)n_reads + 1)) return -1;
    CK(cudaMemcpy(b.rs, reads, reads_bytes, cudaMemcpyHostToDevice));
    CK(cudaMemcpy(b.os, off, (n_reads + 1) * sizeof(int64_t), cudaMemcpyHostToDevice));
    b.key = key; b.age = ++g_rb_clock;
    g_stage_rb[slot] = pick; d_rs[slot] = b.rs; d_os[slot] = b.os;
    return 0;
}

static uint32_t *d_ekey = 0; static size_t cap_ekey = 0;
static uint32_t *d_eidx = 0; static size_t cap_eidx = 0;
static uint8_t *d_etmp = 0;  static size_t cap_etmp = 0;
// ksw_extend2's result read off the main diagonal when no path with a gap can change it (bwa-mem3's ungapped
// shortcut; the conditions and why they suffice: align-rs/src/ext.rs `ungapped_extend`, the same walk): the
// diagonal's total loss (a - score per cell) below min(o_del + e_del, o_ins + e_ins), tlen >= qlen, every diagonal
// value > 0, no z-drop on the diagonal. false: the job needs the DP.
__device__ bool ungapped_ext(const Job &jb, const uint8_t *rd, const uint8_t *pac, int64_t l_pac, int dir, JobOut &r) {
    const int qlen = jb.qlen;
    if (qlen < 1 || (int64_t)jb.rlen < qlen) return false;
    int a = 0;
    for (int i = 0; i < 25; ++i) a = a > c_mat[i] ? a : c_mat[i];
    const int oe = min(c_par[0] + c_par[1], c_par[2] + c_par[3]), zdrop = c_par[4];
    int d = jb.h0, mx = d, mi = -1, loss = 0;
    for (int i = 0; i < qlen; ++i) {
        const int sc = c_mat[rd[jb.qstart + (int64_t)i * dir] * 5 + ref_base(pac, l_pac, jb.rstart + (int64_t)i * dir)];
        loss += a - sc;
        if (loss >= oe) return false;
        d += sc;
        if (d <= 0) return false;
        if (d > mx) { mx = d; mi = i; }
        else if (zdrop > 0 && mx - d > zdrop) return false;
    }
    r.score = mx; r.qle = mi + 1; r.tle = mi + 1; r.gtle = qlen; r.gscore = d; r.max_off = 0;
    return true;
}

// The sort key (qlen << 16 | rlen, clamped), and the ungapped shortcut (ALIGN_RS_UNGAPPED, `mode`): 1, a job it
// takes gets its result here and the largest key (sorted last, k_extend skips it); 2 (check), its result goes to
// ug[] and the job still runs the DP (k_extcheck compares).
__global__ void k_extkey(const Job *jobs, int n, uint32_t *key, uint32_t *idx, int mode, const uint8_t *reads, const int64_t *off, const uint8_t *pac, int64_t l_pac, int dir, JobOut *out, JobOut *ug, int *nug) {
    int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n) return;
    const Job jb = jobs[t];
    uint32_t k = ((uint32_t)jb.qlen << 16) | min(jb.rlen, 65535u);
    if (mode) {
        JobOut r;
        const bool ok = ungapped_ext(jb, reads + off[jb.read], pac, l_pac, dir, r);
        if (mode == 1) {
            if (ok) { out[t] = r; k = 0xffffffffu; atomicAdd(nug, 1); }
        } else {
            if (!ok) r.score = INT_MIN;
            ug[t] = r;
        }
    }
    key[t] = k;
    idx[t] = t;
}
// ALIGN_RS_UNGAPPED=check: the shortcut's results against the DP's
__global__ void k_extcheck(const JobOut *out, const JobOut *ug, int n, int *cnt) {
    int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n || ug[t].score == INT_MIN) return;
    const JobOut a = out[t], b = ug[t];
    atomicAdd(&cnt[0], 1);
    if (a.score != b.score || a.qle != b.qle || a.tle != b.tle || a.gtle != b.gtle || a.gscore != b.gscore || a.max_off != b.max_off) atomicAdd(&cnt[1], 1);
}
static int ungapped_mode() {
    static int m = [] { const char *e = getenv("ALIGN_RS_UNGAPPED"); return e && !strcmp(e, "0") ? 0 : e && !strcmp(e, "check") ? 2 : 1; }();
    return m;
}
static JobOut *d_ug = 0; static size_t cap_ug = 0;
static int *d_nug = 0;
static std::atomic<long long> g_ug_jobs{0}, g_ug_checked{0};

static thread_local double t_ext_call[4] = {0, 0, 0, 0};
static int *d_nleft = 0;
// -1 if the device failed, else the number of jobs the kernel left (score INT_MIN) for the CPU
extern "C" int aln_gpu_extend(const Job *jobs, int n_jobs, int w, int end_bonus, int dir, JobOut *out) {
    if (n_jobs == 0) return 0;
    SERIAL;
    if (grow(&d_jobs, &cap_jobs, (size_t)n_jobs) || grow(&d_out, &cap_out, (size_t)n_jobs)) return -1;
    cudaEventRecord(ev[0]);
    CK(cudaMemcpy(d_jobs, jobs, (size_t)n_jobs * sizeof(Job), cudaMemcpyHostToDevice));
    cudaEventRecord(ev[1]);
    int bs = 128, nb = (n_jobs + bs - 1) / bs;
    // order by size: key qlen << 16 | rlen (clamped), radix-sorted with the job indexes
    if (grow(&d_ekey, &cap_ekey, (size_t)n_jobs * 2) || grow(&d_eidx, &cap_eidx, (size_t)n_jobs * 2)) return -1;
    const int ugm = ungapped_mode();
    if (!d_nug) CK(cudaMalloc((void **)&d_nug, 2 * sizeof(int)));
    CK(cudaMemsetAsync(d_nug, 0, 2 * sizeof(int), cudaStreamPerThread));
    if (ugm == 2 && grow(&d_ug, &cap_ug, (size_t)n_jobs)) return -1;
    k_extkey<<<nb, bs>>>(d_jobs, n_jobs, d_ekey, d_eidx, ugm, d_rs[1], d_os[1], d_pac, g_l_pac, dir, d_out, d_ug, d_nug);
    size_t tb = 0;
    CK(cub::DeviceRadixSort::SortPairs(nullptr, tb, d_ekey, d_ekey + n_jobs, d_eidx, d_eidx + n_jobs, n_jobs, 0, 32, cudaStreamPerThread));
    if (grow(&d_etmp, &cap_etmp, tb)) return -1;
    CK(cub::DeviceRadixSort::SortPairs(d_etmp, tb, d_ekey, d_ekey + n_jobs, d_eidx, d_eidx + n_jobs, n_jobs, 0, 32, cudaStreamPerThread));
    cudaEventRecord(ev[4]);
    // a fixed number of blocks per SM, grid-strided: each thread's DP row lives in local memory, and with more
    // threads resident the rows no longer fit L1
    if (g_ext_bpsm > 0 && nb > g_ext_bpsm * g_n_sm) nb = g_ext_bpsm * g_n_sm;
    if (!d_nleft) CK(cudaMalloc((void **)&d_nleft, sizeof(int)));
    CK(cudaMemsetAsync(d_nleft, 0, sizeof(int), cudaStreamPerThread));
    k_extend<<<nb, bs>>>(d_jobs, d_eidx + n_jobs, n_jobs, d_nug, d_rs[1], d_os[1], d_pac, g_l_pac, d_out, d_nleft, w, end_bonus, dir);
    CK(cudaGetLastError());
    if (ugm == 2) k_extcheck<<<(n_jobs + bs - 1) / bs, bs>>>(d_out, d_ug, n_jobs, d_nug);
    cudaEventRecord(ev[2]);
    CK(cudaMemcpy(out, d_out, (size_t)n_jobs * sizeof(JobOut), cudaMemcpyDeviceToHost));
    int nleft = 0, nug[2] = {0, 0};
    CK(cudaMemcpy(&nleft, d_nleft, sizeof(int), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(nug, d_nug, 2 * sizeof(int), cudaMemcpyDeviceToHost));
    if (ugm == 2) {
        // check: nug[0] jobs the shortcut took, nug[1] of them differ from the DP
        if (nug[1]) { fprintf(stderr, "[align-rs-gpu] ALIGN_RS_UNGAPPED=check: %d of %d shortcut results differ from the DP\n", nug[1], nug[0]); abort(); }
        g_ug_checked += nug[0];
    } else g_ug_jobs += nug[0];
    cudaEventRecord(ev[3]);
    cudaEventSynchronize(ev[3]);
    float a, b, c, srt;
    cudaEventElapsedTime(&a, ev[0], ev[1]); cudaEventElapsedTime(&b, ev[1], ev[2]); cudaEventElapsedTime(&c, ev[2], ev[3]);
    cudaEventElapsedTime(&srt, ev[1], ev[4]);
    t_h2d += a * 1e-3; t_kernel += b * 1e-3; t_d2h += c * 1e-3;
    t_kern[0] += b * 1e-3;
    t_ext_call[0] += a * 1e-3; t_ext_call[1] += srt * 1e-3; t_ext_call[2] += (b - srt) * 1e-3; t_ext_call[3] += c * 1e-3;
    return nleft;
}
// ALIGN_RS_TRACE: this thread's extension device time since the last call (s): upload, sort, kernel, download
extern "C" void aln_gpu_ext_take(double *out) {
    for (int i = 0; i < 4; ++i) { out[i] = t_ext_call[i]; t_ext_call[i] = 0; }
}

extern "C" void aln_gpu_report() {
    fprintf(stderr, "[align-rs-gpu] extension: H2D %.2f s, kernel %.2f s, D2H %.2f s; ungapped shortcut: %lld jobs%s\n", t_h2d, t_kernel, t_d2h,
            ungapped_mode() == 2 ? g_ug_checked.load() : g_ug_jobs.load(), ungapped_mode() == 2 ? " checked against the DP, all identical" : ungapped_mode() == 0 ? " (off)" : "");
}

// ------------------------------------------------------------------------------------------------
// FM index on the device: seeding (smem.rs `collect`) and suffix-array lookups (index.rs `sa_entry`)
// ------------------------------------------------------------------------------------------------
struct CpOcc { int64_t cnt[4]; uint64_t one_hot[4]; };
struct DSmem { int64_t k, l, s; int32_t m, n; };
struct OutSmem { int64_t k; int64_t s; int32_t m, n; };  // = align-rs gpu::OutSmem

__constant__ int64_t c_count[5];
__constant__ int64_t c_sentinel;
static CpOcc *d_cp = 0;
static int8_t *d_sa_ms = 0;
static uint32_t *d_sa_ls = 0;
static OutSmem *d_smem = 0;  static size_t cap_smem = 0;
static int32_t *d_cnt = 0;   static size_t cap_cnt = 0;
static int64_t *d_pos = 0;   static size_t cap_pos = 0;

__device__ __forceinline__ int64_t occ(const CpOcc *cp, int64_t pp, int c) {
    const CpOcc *x = cp + (pp >> 6);
    int y = (int)(pp & 63);
    uint64_t mask = y ? (~0ULL << (64 - y)) : 0ULL;
    return x->cnt[c] + __popcll(x->one_hot[c] & mask);
}

__device__ __forceinline__ DSmem backward_ext(const CpOcc *cp, DSmem sm, int a) {
    int64_t k[4], l[4], s[4];
    int64_t sp = sm.k, ep = sm.k + sm.s;
    for (int b = 0; b < 4; ++b) {
        int64_t os = occ(cp, sp, b), oe = occ(cp, ep, b);
        k[b] = c_count[b] + os;
        s[b] = oe - os;
    }
    int64_t so = (sm.k <= c_sentinel && sm.k + sm.s > c_sentinel) ? 1 : 0;
    l[3] = sm.l + so; l[2] = l[3] + s[3]; l[1] = l[2] + s[2]; l[0] = l[1] + s[1];
    sm.k = k[a]; sm.l = l[a]; sm.s = s[a];
    return sm;
}

#define SMEM_MAXL 256

struct SmemOut { OutSmem *o; int n, cap; bool over; };
__device__ __forceinline__ void push(SmemOut &o, const DSmem &x) {
    if (o.n >= o.cap) { o.over = true; return; }
    OutSmem v; v.k = x.k; v.s = x.s; v.m = x.m; v.n = x.n;
    o.o[o.n++] = v;
}

// getSMEMsOnePosOneThread for one read and start x; returns next_x
__device__ int one_pos(const CpOcc *cp, const uint8_t *q, int len, int x, int64_t min_intv, int min_seed_len, SmemOut &out, DSmem *prev) {
    int next_x = x + 1;
    int a = q[x];
    if (a >= 4) return next_x;
    DSmem sm; sm.m = x; sm.n = x; sm.k = c_count[a]; sm.l = c_count[3 - a]; sm.s = c_count[a + 1] - c_count[a];
    int np = 0;
    int j;
    for (j = x + 1; j < len; ++j) {
        a = q[j];
        next_x = j + 1;
        if (a < 4) {
            DSmem t = sm; t.k = sm.l; t.l = sm.k;
            DSmem n_ = backward_ext(cp, t, 3 - a);
            DSmem ns = n_; ns.k = n_.l; ns.l = n_.k; ns.n = j;
            if (ns.s != sm.s) prev[np++] = sm;
            if (ns.s < min_intv) { next_x = j; break; }
            sm = ns;
        } else break;
    }
    if (sm.s >= min_intv) prev[np++] = sm;
    for (int p = 0; p < np / 2; ++p) { DSmem t = prev[p]; prev[p] = prev[np - p - 1]; prev[np - p - 1] = t; }
    for (int jj = x - 1; jj >= 0; --jj) {
        a = q[jj];
        if (a > 3) break;
        int nc = 0;
        int64_t curr_s = -1;  // `int curr_s` in bwa-mem2
        int p;
        for (p = 0; p < np; ++p) {
            DSmem s0 = prev[p];
            DSmem ns = backward_ext(cp, s0, a); ns.m = jj;
            if (ns.s < min_intv && (s0.n - s0.m + 1) >= min_seed_len) { push(out, s0); break; }
            if (ns.s >= min_intv && ns.s != curr_s) { curr_s = (int64_t)(int32_t)ns.s; prev[nc++] = ns; break; }
        }
        p++;
        for (; p < np; ++p) {
            DSmem s0 = prev[p];
            DSmem ns = backward_ext(cp, s0, a); ns.m = jj;
            if (ns.s >= min_intv && ns.s != curr_s) { curr_s = (int64_t)(int32_t)ns.s; prev[nc++] = ns; }
        }
        np = nc;
        if (nc == 0) break;
    }
    if (np != 0) {
        DSmem s0 = prev[0];
        if ((s0.n - s0.m + 1) >= min_seed_len) push(out, s0);
    }
    return next_x;
}

#ifndef SMEM_BS
#define SMEM_BS 64
#endif
#ifndef SMEM_MINB
#define SMEM_MINB 1
#endif
__global__ void __launch_bounds__(SMEM_BS, SMEM_MINB) k_smem(int n_reads, const uint8_t *reads, const int64_t *off, const CpOcc *cp, int cap, OutSmem *outs, int32_t *cnt,
                       int min_seed_len, int split_len, int64_t split_width, int64_t max_mem_intv) {
    int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= n_reads) return;
    const uint8_t *q = reads + off[r];
    int len = (int)(off[r + 1] - off[r]);
    if (len > SMEM_MAXL) { cnt[r] = -1; return; }
    DSmem prev[SMEM_MAXL + 1];
    SmemOut out; out.o = outs + (int64_t)r * cap; out.n = 0; out.cap = cap; out.over = false;
    // round 1
    int x = 0;
    while (x < len && !out.over) x = one_pos(cp, q, len, x, 1, min_seed_len, out, prev);
    // round 2
    int n1 = out.n;
    for (int i = 0; i < n1 && !out.over; ++i) {
        OutSmem p = out.o[i];
        int start = p.m, end = p.n + 1;
        if (end - start < split_len || p.s > split_width) continue;
        int pos = (end + start) >> 1;
        int64_t mi = (int64_t)(int32_t)(p.s + 1);
        one_pos(cp, q, len, pos, mi, min_seed_len, out, prev);
    }
    // round 3 (bwtSeedStrategyAllPosOneThread)
    if (max_mem_intv > 0) {
        int msl = min_seed_len + 1;
        x = 0;
        while (x < len && !out.over) {
            int next_x = x + 1;
            int a = q[x];
            if (a < 4) {
                DSmem sm; sm.m = x; sm.n = x; sm.k = c_count[a]; sm.l = c_count[3 - a]; sm.s = c_count[a + 1] - c_count[a];
                for (int j = x + 1; j < len; ++j) {
                    next_x = j + 1;
                    a = q[j];
                    if (a < 4) {
                        DSmem t = sm; t.k = sm.l; t.l = sm.k;
                        DSmem n_ = backward_ext(cp, t, 3 - a);
                        DSmem ns = n_; ns.k = n_.l; ns.l = n_.k; ns.n = j;
                        sm = ns;
                        if (sm.s < max_mem_intv && (sm.n - sm.m + 1) >= msl) {
                            if (sm.s > 0) push(out, sm);
                            break;
                        }
                    } else break;
                }
            }
            x = next_x;
        }
    }
    cnt[r] = out.over ? -1 : out.n;
}

// ------------------------------------------------------------------------------------------------
// Seeding with a group of G lanes per read (k_smem: one thread per read left 3 of 32 lanes busy, a few
// repetitive reads per warp running ~10x longer than the rest). The sequential parts (forward walks, rounds 2
// and 3) run in lockstep on all lanes of the group (identical loads); in the backward phase the candidates are
// extended G at a time, one per lane, and every lane then replays bwa's selection over them in order (values
// shuffled from their lanes): the same decisions, pushes and in-place compaction as one_pos. Groups are
// persistent: each takes the next read from a counter; its `prev` array lives in its own scratch slot.
// ------------------------------------------------------------------------------------------------
__device__ __forceinline__ void gpush(SmemOut &o, const DSmem &x, bool writer) {
    if (o.n >= o.cap) { o.over = true; return; }
    if (writer) { OutSmem v; v.k = x.k; v.s = x.s; v.m = x.m; v.n = x.n; o.o[o.n] = v; }
    o.n++;
}

template <int G>
__device__ int one_pos_g(const CpOcc *cp, const uint8_t *q, int len, int x, int64_t min_intv, int min_seed_len, SmemOut &out, DSmem *prev, int lane, unsigned mask) {
    int next_x = x + 1;
    int a = q[x];
    if (a >= 4) return next_x;
    DSmem sm; sm.m = x; sm.n = x; sm.k = c_count[a]; sm.l = c_count[3 - a]; sm.s = c_count[a + 1] - c_count[a];
    int np = 0;
    for (int j = x + 1; j < len; ++j) {
        a = q[j];
        next_x = j + 1;
        if (a < 4) {
            DSmem t = sm; t.k = sm.l; t.l = sm.k;
            DSmem n_ = backward_ext(cp, t, 3 - a);
            DSmem ns = n_; ns.k = n_.l; ns.l = n_.k; ns.n = j;
            if (ns.s != sm.s) { if (lane == 0) prev[np] = sm; np++; }
            if (ns.s < min_intv) { next_x = j; break; }
            sm = ns;
        } else break;
    }
    if (sm.s >= min_intv) { if (lane == 0) prev[np] = sm; np++; }
    __syncwarp(mask);
    for (int p = lane; p < np / 2; p += G) { DSmem t = prev[p]; prev[p] = prev[np - p - 1]; prev[np - p - 1] = t; }
    __syncwarp(mask);
    for (int jj = x - 1; jj >= 0; --jj) {
        a = q[jj];
        if (a > 3) break;
        int nc = 0;
        int64_t curr_s = -1;  // `int curr_s` in bwa-mem2
        bool first_loop = true;  // one_pos's first loop (up to its break), then its second
        for (int base = 0; base < np; base += G) {
            const int p = base + lane;
            DSmem s0, ns;
            s0.k = s0.l = s0.s = 0; s0.m = s0.n = 0; ns = s0;
            if (p < np) { s0 = prev[p]; ns = backward_ext(cp, s0, a); ns.m = jj; }
            __syncwarp(mask);  // the chunk is read before any of it is overwritten
            const int lim = np - base < G ? np - base : G;
            for (int i = 0; i < lim; ++i) {
                const int64_t nss = __shfl_sync(mask, ns.s, i, G);
                if (first_loop) {
                    const int len0 = __shfl_sync(mask, s0.n - s0.m + 1, i, G);
                    if (nss < min_intv && len0 >= min_seed_len) { gpush(out, s0, lane == i); first_loop = false; continue; }
                    if (nss >= min_intv && nss != curr_s) {
                        curr_s = (int64_t)(int32_t)nss;
                        if (lane == i) prev[nc] = ns;
                        nc++; first_loop = false;
                    }
                } else if (nss >= min_intv && nss != curr_s) {
                    curr_s = (int64_t)(int32_t)nss;
                    if (lane == i) prev[nc] = ns;
                    nc++;
                }
            }
        }
        __syncwarp(mask);
        np = nc;
        if (nc == 0) break;
    }
    if (np != 0) {
        DSmem s0 = prev[0];
        if ((s0.n - s0.m + 1) >= min_seed_len) gpush(out, s0, lane == 0);
    }
    __syncwarp(mask);
    return next_x;
}

template <int G>
__global__ void __launch_bounds__(128) k_smem_g(int n_reads, const uint8_t *reads, const int64_t *off, const CpOcc *cp, int cap, OutSmem *outs, int32_t *cnt,
                       int min_seed_len, int split_len, int64_t split_width, int64_t max_mem_intv, DSmem *scratch, int *next_read) {
    const int lane = threadIdx.x & (G - 1);
    const unsigned mask = G == 32 ? 0xffffffffu : (((1u << G) - 1) << ((threadIdx.x & 31) & ~(G - 1)));
    const int64_t gid = ((int64_t)blockIdx.x * blockDim.x + threadIdx.x) / G;
    DSmem *prev = scratch + gid * (SMEM_MAXL + 1);
    for (;;) {
        int r = 0;
        if (lane == 0) r = atomicAdd(next_read, 1);
        r = __shfl_sync(mask, r, 0, G);
        if (r >= n_reads) return;
        const uint8_t *q = reads + off[r];
        int len = (int)(off[r + 1] - off[r]);
        if (len > SMEM_MAXL) { if (lane == 0) cnt[r] = -1; continue; }
        SmemOut out; out.o = outs + (int64_t)r * cap; out.n = 0; out.cap = cap; out.over = false;
        // round 1
        int x = 0;
        while (x < len && !out.over) x = one_pos_g<G>(cp, q, len, x, 1, min_seed_len, out, prev, lane, mask);
        // round 2
        int n1 = out.n;
        for (int i = 0; i < n1 && !out.over; ++i) {
            OutSmem p = out.o[i];
            int start = p.m, end = p.n + 1;
            if (end - start < split_len || p.s > split_width) continue;
            int pos = (end + start) >> 1;
            int64_t mi = (int64_t)(int32_t)(p.s + 1);
            one_pos_g<G>(cp, q, len, pos, mi, min_seed_len, out, prev, lane, mask);
        }
        // round 3 (bwtSeedStrategyAllPosOneThread)
        if (max_mem_intv > 0) {
            int msl = min_seed_len + 1;
            x = 0;
            while (x < len && !out.over) {
                int next_x = x + 1;
                int a = q[x];
                if (a < 4) {
                    DSmem sm; sm.m = x; sm.n = x; sm.k = c_count[a]; sm.l = c_count[3 - a]; sm.s = c_count[a + 1] - c_count[a];
                    for (int j = x + 1; j < len; ++j) {
                        next_x = j + 1;
                        a = q[j];
                        if (a < 4) {
                            DSmem t = sm; t.k = sm.l; t.l = sm.k;
                            DSmem n_ = backward_ext(cp, t, 3 - a);
                            DSmem ns = n_; ns.k = n_.l; ns.l = n_.k; ns.n = j;
                            sm = ns;
                            if (sm.s < max_mem_intv && (sm.n - sm.m + 1) >= msl) {
                                if (sm.s > 0) gpush(out, sm, lane == 0);
                                break;
                            }
                        } else break;
                    }
                }
                x = next_x;
            }
        }
        if (lane == 0) cnt[r] = out.over ? -1 : out.n;
        __syncwarp(mask);
    }
}

// k_smem_g's launch: as many groups as fit on the device at once, each with a scratch slot
static DSmem *d_sscr = 0; static size_t cap_sscr = 0;
static int *d_snext = 0;
template <int G>
static int launch_smem_g(int n, const uint8_t *rs, const int64_t *os, int cap, int min_seed_len, int split_len, int64_t split_width, int64_t max_mem_intv) {
    static int blocks = 0;
    const int bs = 128;
    if (!blocks) {
        int per_sm = 0, dev = 0, n_sm = 0;
        CK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&per_sm, k_smem_g<G>, bs, 0));
        CK(cudaGetDevice(&dev));
        CK(cudaDeviceGetAttribute(&n_sm, cudaDevAttrMultiProcessorCount, dev));
        blocks = per_sm * n_sm;
    }
    size_t groups = (size_t)blocks * bs / G;
    if (grow(&d_sscr, &cap_sscr, groups * (SMEM_MAXL + 1))) return -1;
    if (!d_snext) CK(cudaMalloc((void **)&d_snext, sizeof(int)));
    CK(cudaMemsetAsync(d_snext, 0, sizeof(int), cudaStreamPerThread));
    k_smem_g<G><<<blocks, bs>>>(n, rs, os, d_cp, cap, d_smem, d_cnt, min_seed_len, split_len, split_width, max_mem_intv, d_sscr, d_snext);
    CK(cudaGetLastError());
    return 0;
}
// ALIGN_RS_SMEM_G: lanes per read (4, 8, 16, 32; 1 = k_smem, a thread per read)
static int smem_g() { return g_smem_g; }

__device__ __forceinline__ int64_t sa_sample(const int8_t *ms, const uint32_t *ls, int64_t i) {
    return ((int64_t)ms[i] << 32) + (int64_t)ls[i];
}

__global__ void k_sa(const int64_t *pos, int n, const CpOcc *cp, const int8_t *ms, const uint32_t *ls, int64_t *out) {
    int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n) return;
    int64_t p = pos[t];
    if ((p & 7) == 0) { out[t] = sa_sample(ms, ls, p >> 3); return; }
    int64_t sp = p, offset = 0;
    for (;;) {
        const CpOcc *x = cp + (sp >> 6);
        int y = 63 - (int)(sp & 63);
        int b;
        if ((x->one_hot[0] >> y) & 1) b = 0;
        else if ((x->one_hot[1] >> y) & 1) b = 1;
        else if ((x->one_hot[2] >> y) & 1) b = 2;
        else if ((x->one_hot[3] >> y) & 1) b = 3;
        else { out[t] = 0; return; }  // the sentinel: call_one_step returns 0
        sp = c_count[b] + occ(cp, sp, b);
        offset++;
        if ((sp & 7) == 0) { out[t] = sa_sample(ms, ls, sp >> 3) + offset; return; }
    }
}

// Loads the FM index if it fits in free device memory (keeping `reserve` bytes free): 0 loaded,
// 1 not enough memory (nothing loaded), -1 error. `with_sa` 0: the checkpoints only.
extern "C" int aln_gpu_load_fmi(const void *cp, int64_t n_cp, const int8_t *ms, const uint32_t *ls, int64_t n_sa, const int64_t *count, int64_t sentinel, int with_sa, int64_t reserve) {
    size_t fr = 0, tot = 0;
    CK(cudaMemGetInfo(&fr, &tot));
    size_t need = (size_t)n_cp * sizeof(CpOcc) + (with_sa ? (size_t)n_sa * 5 : 0);
    if (need + (size_t)reserve > fr) return 1;
    CK(cudaMalloc((void **)&d_cp, (size_t)n_cp * sizeof(CpOcc)));
    CK(cudaMemcpy(d_cp, cp, (size_t)n_cp * sizeof(CpOcc), cudaMemcpyHostToDevice));
    if (with_sa) {
        CK(cudaMalloc((void **)&d_sa_ms, (size_t)n_sa));
        CK(cudaMemcpy(d_sa_ms, ms, (size_t)n_sa, cudaMemcpyHostToDevice));
        CK(cudaMalloc((void **)&d_sa_ls, (size_t)n_sa * 4));
        CK(cudaMemcpy(d_sa_ls, ls, (size_t)n_sa * 4, cudaMemcpyHostToDevice));
    }
    CK(cudaMemcpyToSymbol(c_count, count, 5 * sizeof(int64_t)));
    CK(cudaMemcpyToSymbol(c_sentinel, &sentinel, sizeof(int64_t)));
    return 0;
}

extern "C" int64_t aln_gpu_free_mem() {
    size_t fr = 0, tot = 0;
    if (cudaMemGetInfo(&fr, &tot) != cudaSuccess) return -1;
    return (int64_t)fr;
}

// Seeds of the reads set by aln_gpu_set_reads: up to `cap` per read; cnt[r] = -1 when the read did not fit
// (longer than SMEM_MAXL or more seeds than cap): the CPU redoes it.
// k_smem's output, per read in rows of `cap`, packed: seeds of read r at off[r] (exclusive sum of the counts,
// reads not done counting 0)
static OutSmem *d_spack = 0;  static size_t cap_spack = 0;
static int64_t *d_soff = 0;   static size_t cap_soff = 0;
static int32_t *d_sc2 = 0;    static size_t cap_sc2 = 0;
static uint8_t *d_stmp = 0;   static size_t cap_stmp = 0;
__global__ void k_scount(const int32_t *cnt, int n, int32_t *c2) {
    int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r < n) c2[r] = cnt[r] > 0 ? cnt[r] : 0;
}
__global__ void k_spack(const OutSmem *rows, int cap, const int32_t *c2, const int64_t *off, int n, OutSmem *pack) {
    int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= n) return;
    const OutSmem *src = rows + (int64_t)r * cap;
    OutSmem *dst = pack + off[r];
    for (int i = 0; i < c2[r]; ++i) dst[i] = src[i];
}

// Seeds of reads r0..r0+n (set by aln_gpu_set_reads): up to `cap` per read; cnt[r] = -1 when the read did not fit
// (longer than SMEM_MAXL or more seeds than cap): the CPU redoes it. The seeds stay on the device, packed; returns
// their number (aln_gpu_smem_fetch copies them), -1 on error.
extern "C" int64_t aln_gpu_smem(int r0, int n, int cap, int min_seed_len, int split_len, int64_t split_width, int64_t max_mem_intv, int32_t *cnt) {
    if (n == 0) return 0;
    SERIAL;
    if (grow(&d_smem, &cap_smem, (size_t)n * cap) || grow(&d_cnt, &cap_cnt, (size_t)n) || grow(&d_sc2, &cap_sc2, (size_t)n) || grow(&d_soff, &cap_soff, (size_t)n + 1)) return -1;
    cudaEvent_t e0 = tev()[0], e1 = tev()[1], ek = tev()[2]; cudaEventRecord(e0);
    int bs = SMEM_BS, nb = (n + bs - 1) / bs;
    switch (smem_g()) {
    case 2: if (launch_smem_g<2>(n, d_rs[0], d_os[0] + r0, cap, min_seed_len, split_len, split_width, max_mem_intv)) return -1; break;
    case 4: if (launch_smem_g<4>(n, d_rs[0], d_os[0] + r0, cap, min_seed_len, split_len, split_width, max_mem_intv)) return -1; break;
    case 8: if (launch_smem_g<8>(n, d_rs[0], d_os[0] + r0, cap, min_seed_len, split_len, split_width, max_mem_intv)) return -1; break;
    case 16: if (launch_smem_g<16>(n, d_rs[0], d_os[0] + r0, cap, min_seed_len, split_len, split_width, max_mem_intv)) return -1; break;
    case 32: if (launch_smem_g<32>(n, d_rs[0], d_os[0] + r0, cap, min_seed_len, split_len, split_width, max_mem_intv)) return -1; break;
    default:
        k_smem<<<nb, bs>>>(n, d_rs[0], d_os[0] + r0, d_cp, cap, d_smem, d_cnt, min_seed_len, split_len, split_width, max_mem_intv);
        CK(cudaGetLastError());
    }
    cudaEventRecord(ek);
    k_scount<<<(n + 255) / 256, 256>>>(d_cnt, n, d_sc2);
    size_t tb = 0;
    CK(cub::DeviceScan::ExclusiveSum(nullptr, tb, d_sc2, d_soff, n, cudaStreamPerThread));
    if (grow(&d_stmp, &cap_stmp, tb)) return -1;
    CK(cub::DeviceScan::ExclusiveSum(d_stmp, tb, d_sc2, d_soff, n, cudaStreamPerThread));
    int64_t off_last = 0; int32_t c_last = 0;
    CK(cudaMemcpy(cnt, d_cnt, (size_t)n * sizeof(int32_t), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(&off_last, d_soff + n - 1, sizeof(int64_t), cudaMemcpyDeviceToHost));
    c_last = cnt[n - 1] > 0 ? cnt[n - 1] : 0;
    int64_t total = off_last + c_last;
    if (grow(&d_spack, &cap_spack, (size_t)total + 1)) return -1;
    k_spack<<<(n + 127) / 128, 128>>>(d_smem, cap, d_sc2, d_soff, n, d_spack);
    CK(cudaGetLastError());
    cudaEventRecord(e1); cudaEventSynchronize(e1);
    float ms_; cudaEventElapsedTime(&ms_, e0, e1); t_smem += ms_ * 1e-3;
    add_kern(1, e0, ek);
    return total;
}

// the packed seeds of the last aln_gpu_smem call
extern "C" int aln_gpu_smem_fetch(OutSmem *out, int64_t total) {
    if (total == 0) return 0;
    SERIAL;
    CK(cudaMemcpy(out, d_spack, (size_t)total * sizeof(OutSmem), cudaMemcpyDeviceToHost));
    return 0;
}

extern "C" int aln_gpu_sa(const int64_t *pos, int64_t n, int64_t *out) {
    if (n == 0) return 0;
    int64_t *d_out64 = 0;
    SERIAL;
    if (grow(&d_pos, &cap_pos, (size_t)n * 2)) return -1;
    d_out64 = d_pos + n;
    cudaEvent_t e0 = tev()[0], e1 = tev()[1]; cudaEventRecord(e0);
    CK(cudaMemcpy(d_pos, pos, (size_t)n * sizeof(int64_t), cudaMemcpyHostToDevice));
    int bs = 128; int64_t nb = (n + bs - 1) / bs;
    cudaEvent_t ek0 = tev()[2], ek1 = tev()[3]; cudaEventRecord(ek0);
    k_sa<<<(unsigned)nb, bs>>>(d_pos, (int)n, d_cp, d_sa_ms, d_sa_ls, d_out64);
    CK(cudaGetLastError());
    cudaEventRecord(ek1);
    CK(cudaMemcpy(out, d_out64, (size_t)n * sizeof(int64_t), cudaMemcpyDeviceToHost));
    cudaEventRecord(e1); cudaEventSynchronize(e1);
    float ms_; cudaEventElapsedTime(&ms_, e0, e1); t_sa += ms_ * 1e-3;
    add_kern(2, ek0, ek1);
    return 0;
}

extern "C" void aln_gpu_report2() {
    fprintf(stderr, "[align-rs-gpu] seeding %.2f s, SA lookups %.2f s (incl. transfers)\n", t_smem, t_sa);
}

// ------------------------------------------------------------------------------------------------
// Mate rescue: ksw.cpp ksw_align2 with KSW_XBYTE (ksw_u8), 16 threads per alignment, thread l = byte lane l of
// the SSE vectors: _mm_slli_si128(x, 1) is a shuffle up by one lane, the row maximum a reduction, the lazy-F exit
// a vote. Saturating u8 arithmetic as the SSE instructions do. Forward pass, then (KSW_XSTART) the reverse pass
// that finds the start, as ksw_align2 runs them.
// ------------------------------------------------------------------------------------------------
struct MJob { uint32_t read; int32_t qlen, rev, tlen; int64_t rb; int32_t xtra, pad; int64_t boff; };
struct MOut { int32_t score, te, qe, score2, te2, tb, qb, ok; };

#define MS_SLEN 16   // queries up to 256 bases
// per-row maxima for score2: room for one per target row (b entry: imax << 20 | row)

__device__ __forceinline__ uint8_t adds8(uint8_t a, uint8_t b) { unsigned s = a + b; return s > 255 ? 255 : s; }
__device__ __forceinline__ uint8_t subs8(uint8_t a, uint8_t b) { return a > b ? a - b : 0; }
__device__ __forceinline__ uint8_t max8(uint8_t a, uint8_t b) { return a > b ? a : b; }

struct MRes { int score, te, qe, score2, te2; bool over; };

// one ksw_u8 pass; qbase(k) / tbase(i): the query / target bases; lane = 0..15, mask = the half-warp.
// SL = slen (the query's 16-byte segments) is a template parameter, so that the per-lane rows (H, E, the saved row,
// the query bases) are register arrays: with a runtime slen they lived in local memory, a load or store per cell.
// The score of query base q against target base t is bwa's matrix (bwa_fill_scmat: match A, mismatch B, N against
// anything NS; aln_gpu_init checks the matrix has that form); positions past the query end score 0, as ksw_qinit
// pads them.
// half-warp (16-lane segment) max / min: redux.sync from sm_80, shuffles before
__device__ __forceinline__ unsigned red_max16(unsigned mask, unsigned v) {
#if __CUDA_ARCH__ >= 800 && !defined(ALIGN_RS_NO_REDUX)
    return __reduce_max_sync(mask, v);
#else
    for (int o = 8; o; o >>= 1) v = max(v, (unsigned)__shfl_xor_sync(mask, v, o, 16));
    return v;
#endif
}
__device__ __forceinline__ unsigned red_min16(unsigned mask, unsigned v) {
#if __CUDA_ARCH__ >= 800 && !defined(ALIGN_RS_NO_REDUX)
    return __reduce_min_sync(mask, v);
#else
    for (int o = 8; o; o >>= 1) v = min(v, (unsigned)__shfl_xor_sync(mask, v, o, 16));
    return v;
#endif
}

__device__ __forceinline__ int adds8i(int a, int b) { return min(a + b, 255); }
__device__ __forceinline__ int subs8i(int a, int b) { return max(a - b, 0); }
__device__ __forceinline__ int max8i(int a, int b) { return max(a, b); }

template <int SL, typename QF, typename TF>
__device__ MRes ksw_u8_pass_t(int qlen, int tlen, QF qbase, TF tbase, int xtra, int lane, unsigned mask, uint32_t *b, int bcap) {
    const int o_del = c_par[0], e_del = c_par[1], o_ins = c_par[2], e_ins = c_par[3];
    // shift = 256 - min(mat) as uint8 (mat min -4 -> 4); max = max(mat)
    int mn = 127, mx = 0;
    for (int a = 0; a < 25; ++a) { mn = c_mat[a] < mn ? c_mat[a] : mn; mx = c_mat[a] > mx ? c_mat[a] : mx; }
    const uint8_t shift = (uint8_t)(256 - (int)(uint8_t)(int8_t)mn);
    const int qmax = mx;
    const int sA = (uint8_t)(int8_t)(c_mat[0] + shift), sB = (uint8_t)(int8_t)(c_mat[1] + shift), sN = (uint8_t)(int8_t)(c_mat[4] + shift), s0 = shift;
    // per query position: its base (255 past the end) and its score against a target base other than itself and
    // other than N (mismatch; N against anything; 0 past the end), all less the shift
    int qv[SL], base[SL];
#pragma unroll
    for (int j = 0; j < SL; ++j) {
        int k = j + lane * SL;
        qv[j] = k >= qlen ? 255 : qbase(k);
        base[j] = (qv[j] == 255 ? s0 : (qv[j] > 3 ? sN : sB)) - shift;
    }
    const int sAs = sA - shift, sNs = sN - shift, s0s = 0;
    // H never reaches 255 when Q * max + max + shift <= 255 (max(h + s - shift, 0) is then u8's saturating result)
    const bool sat = SL * 16 * qmax + qmax + shift > 255;
    const int minsc = (xtra & 0x40000) ? (xtra & 0xffff) : 0x10000;
    const int endsc = (xtra & 0x20000) ? (xtra & 0xffff) : 0x10000;
    const int oe_del = (uint8_t)(o_del + e_del), ed = (uint8_t)e_del, oe_ins = (uint8_t)(o_ins + e_ins), ei = (uint8_t)e_ins;
    // u8 values held in ints (no byte packing): saturating add / subtract as min(x, 255) / max(x, 0)
    int H0[SL], H1[SL], E[SL];
#pragma unroll
    for (int j = 0; j < SL; ++j) { H0[j] = 0; E[j] = 0; H1[j] = 0; }
    // qe as ksw_u8 takes it from the saved row (the smallest position of the row's largest value), chosen when
    // the row is saved instead of saving it (10 registers fewer)
    int qe_save = 0;
    int n_b = 0, te = -1, gmax = 0;
    uint32_t last_b = 0;   // b[n_b - 1]
    bool over = false;
    int tc_next = tlen > 0 ? tbase(0) : 0;
    for (int i = 0; i < tlen; ++i) {
        const int tc = tc_next;
        if (i + 1 < tlen) tc_next = tbase(i + 1);   // issued a row ahead
        int f = 0, mxv = 0;
        int h = __shfl_up_sync(mask, H0[SL - 1], 1, 16);
        if (lane == 0) h = 0;
        // the first sweep; tc != 4 && !sat is the common case
#define MS_CELL(SC, SATURATE) { \
            if (SATURATE) { h = subs8i(adds8i(h, (SC) + shift), shift); } else h = max(h + (SC), 0); \
            int e = E[j]; \
            h = max(h, e); \
            h = max(h, f); \
            mxv = max(mxv, h); \
            H1[j] = h; \
            E[j] = max(e - ed, max(h - oe_del, 0)); \
            f = max(f - ei, max(h - oe_ins, 0)); \
            h = H0[j]; }
        if (tc < 4 && !sat) {
#pragma unroll
            for (int j = 0; j < SL; ++j) MS_CELL(qv[j] == tc ? sAs : base[j], false)
        } else if (tc < 4) {
#pragma unroll
            for (int j = 0; j < SL; ++j) MS_CELL(qv[j] == tc ? sAs : base[j], true)
        } else {
#pragma unroll
            for (int j = 0; j < SL; ++j) MS_CELL(qv[j] == 255 ? s0s : sNs, sat)
        }
#undef MS_CELL
        for (int k = 0; k < 16; ++k) {
            int fs = __shfl_up_sync(mask, f, 1, 16);
            f = lane == 0 ? 0 : fs;
            bool done = false;
#pragma unroll
            for (int j = 0; j < SL; ++j) {
                int hh = max8i(H1[j], f);
                H1[j] = hh;
                hh = subs8i(hh, oe_ins);
                f = subs8i(f, ei);
                unsigned all = __ballot_sync(mask, f <= hh);
                if ((all & mask) == mask) { done = true; break; }
            }
            if (done) break;
        }
        const int imax = (int)red_max16(mask, (unsigned)mxv);
        if (imax >= minsc) {
            if (lane == 0) {
                if (n_b == 0 || (int)(last_b & 0xfffff) + 1 != i) {
                    last_b = (uint32_t)imax << 20 | (uint32_t)i;
                    if (n_b < bcap) b[n_b] = last_b;
                    else over = true;
                    n_b++;
                } else if ((int)(last_b >> 20) < imax) {
                    last_b = (uint32_t)imax << 20 | (uint32_t)i;
                    if (n_b <= bcap) b[n_b - 1] = last_b;
                }
            } else {
                // keep n_b in step on the other lanes (only lane 0 reads b)
                if (n_b == 0) n_b = 1;
            }
        }
        if (imax > gmax) {
            gmax = imax; te = i;
            int bv = -1, bp = 0x7fffffff;
#pragma unroll
            for (int j = 0; j < SL; ++j) {
                int v = H1[j], p = j + lane * SL;
                if (v > bv || (v == bv && p < bp)) { bv = v; bp = p; }
            }
            int mv = (int)red_max16(mask, (unsigned)bv);
            qe_save = (int)red_min16(mask, (unsigned)(bv == mv ? bp : 0x7fffffff));
            if (gmax + shift >= 255 || gmax >= endsc) break;
        }
#pragma unroll
        for (int j = 0; j < SL; ++j) H0[j] = H1[j];
    }
    MRes r; r.score = gmax + shift < 255 ? gmax : 255; r.te = te; r.qe = -1; r.score2 = -1; r.te2 = -1; r.over = false;
    if (r.score != 255) {
        // qe: the smallest query position among the cells of the saved row with the maximal value
        r.qe = te >= 0 ? qe_save : 0;
        if (lane == 0 && n_b > 0) {
            int w = (r.score + qmax - 1) / qmax;
            int low = te - w, high = te + w;
            int nb = n_b < bcap ? n_b : bcap;
            for (int k = 0; k < nb; ++k) {
                int e = (int)(b[k] & 0xfffff);
                if ((e < low || e > high) && (int)(b[k] >> 20) > r.score2) { r.score2 = (int)(b[k] >> 20); r.te2 = e; }
            }
            r.over = over;
        }
    }
    return r;
}

template <typename QF, typename TF>
__device__ MRes ksw_u8_pass(int qlen, int tlen, QF qbase, TF tbase, int xtra, int lane, unsigned mask, uint32_t *b, int bcap) {
    switch ((qlen + 15) / 16) {
#define MS_CASE(n) case n: return ksw_u8_pass_t<n>(qlen, tlen, qbase, tbase, xtra, lane, mask, b, bcap);
        MS_CASE(1) MS_CASE(2) MS_CASE(3) MS_CASE(4) MS_CASE(5) MS_CASE(6) MS_CASE(7) MS_CASE(8)
        MS_CASE(9) MS_CASE(10) MS_CASE(11) MS_CASE(12) MS_CASE(13) MS_CASE(14) MS_CASE(15) MS_CASE(16)
#undef MS_CASE
    }
    MRes r; r.score = 0; r.te = -1; r.qe = -1; r.score2 = -1; r.te2 = -1; r.over = true;  // longer query: CPU
    return r;
}

// 64 threads (4 jobs) per block, registers capped for MINB blocks per SM: the kernel waits on its shuffles and
// votes row after row, so more resident warps beat fewer registers (some spilling in the long-query variants).
// RTX 5070 Ti: MINB 18 (56 registers) 164 -> 306 GCUPS on a lane's chunk (cuda/bench_msw.cu), results identical.
// Variants 12 / 18 / 24, one chosen per device (gpu_settings).
#define MSW_BS 64
template <int MINB>
__global__ void __launch_bounds__(MSW_BS, MINB) k_matesw(const MJob *jobs, int n, const uint8_t *reads, const int64_t *off, const uint8_t *pac, int64_t l_pac, MOut *out, uint32_t *bscratch) {
    int gt = blockIdx.x * blockDim.x + threadIdx.x;
    int job = gt >> 4, lane = threadIdx.x & 15;
    unsigned mask = 0xffffu << (threadIdx.x & 16);
    if (job >= n) return;  // both halves of a warp: whole half-warps leave together
    const MJob jb = jobs[job];
    const uint8_t *rd = reads + off[jb.read];
    const int lq = jb.qlen;
    auto qf = [&](int k) -> int { int c = jb.rev ? rd[lq - 1 - k] : rd[k]; return jb.rev ? (c < 4 ? 3 - c : 4) : c; };
    auto tf = [&](int i) -> int { return ref_base(pac, l_pac, jb.rb + i); };
    uint32_t *b = bscratch + jb.boff;
    const int bcap = jb.tlen + 1;
    MRes r = ksw_u8_pass(lq, jb.tlen, qf, tf, jb.xtra, lane, mask, b, bcap);
    MOut o; o.score = r.score; o.te = r.te; o.qe = r.qe; o.score2 = r.score2; o.te2 = r.te2; o.tb = -1; o.qb = -1; o.ok = 1;
    if (lane == 0 && r.over) o.ok = 0;
    const int xtra = jb.xtra;
    if (!((xtra & 0x80000) == 0 || ((xtra & 0x40000) && r.score < (xtra & 0xffff)))) {
        // reverse pass: the query's first qe+1 bases reversed, the target's first te+1 reversed (the rest as is)
        const int qe = r.qe, te = r.te;
        auto qf2 = [&](int k) -> int { return qf(qe - k); };
        auto tf2 = [&](int i) -> int { return i <= te ? tf(te - i) : tf(i); };
        MRes rr = ksw_u8_pass(qe + 1, jb.tlen, qf2, tf2, 0x20000 | r.score, lane, mask, b, bcap);
        if (r.score == rr.score) { o.tb = r.te - rr.te; o.qb = r.qe - rr.qe; }
    }
    if (lane == 0) out[job] = o;
}


// ------------------------------------------------------------------------------------------------
// Mate rescue as a wavefront (the default; ALIGN_RS_MSW=striped runs k_matesw): the same ksw_u8 results, 16 threads per
// job. Thread t owns query segment t (positions t*SL .. t*SL+SL-1, as SSE lane t does) in registers and runs one
// target row behind thread t-1, which hands it the row's boundary by shuffle: the final H at its last position (the
// next row's diagonal), the F leaving it, the row maximum so far, and (when the row may beat the best so far) the
// row's largest final H and its first position. No lazy-F loop, no per-row vote or reduction (the layout of
// CUDASW++4.0 / gpuPairHMM). bwa's striped semantics are kept: the first sweep's F restarts at each segment, and the
// row maximum and E take the values before the lazy-F pass (hf); H takes the full F (hF). So two F chains: fs within
// the segment, ft across. E and the F chains are not clamped at 0 (H >= 0 absorbs a negative value: the same results).
// The last thread keeps lane 0's state (b, n_b, gmax, te, qe). On a dumped lane chunk (65,536 jobs, all fields
// identical): RTX 5070 Ti 31 -> 18 ms, L40S 21 -> 12 ms (kernel time). Register cap: 12 blocks of 64 per SM (80
// registers), best on both.
// ------------------------------------------------------------------------------------------------
#define MSW_WAVE_MINB 12
template <int SL, bool SAT, typename QF, typename TF>
__device__ MRes ksw_u8_wave_t(int qlen, int tlen, QF qbase, TF tbase, int xtra, int t, unsigned mask, uint32_t *b, int bcap) {
    const int o_del = c_par[0], e_del = c_par[1], o_ins = c_par[2], e_ins = c_par[3];
    int mn = 127, mx = 0;
    for (int a = 0; a < 25; ++a) { mn = c_mat[a] < mn ? c_mat[a] : mn; mx = c_mat[a] > mx ? c_mat[a] : mx; }
    const int shift = (uint8_t)(256 - (int)(uint8_t)(int8_t)mn);
    const int qmax = mx;
    const int sA = (uint8_t)(int8_t)(c_mat[0] + shift), sB = (uint8_t)(int8_t)(c_mat[1] + shift), sN = (uint8_t)(int8_t)(c_mat[4] + shift);
    const int sAs = sA - shift, sNs = sN - shift;
    int qv[SL], base[SL], baseN[SL];
#pragma unroll
    for (int j = 0; j < SL; ++j) {
        int k = j + t * SL;
        qv[j] = k >= qlen ? 255 : qbase(k);
        base[j] = qv[j] == 255 ? 0 : (qv[j] > 3 ? sNs : sB - shift);
        baseN[j] = qv[j] == 255 ? 0 : sNs;
    }
    const int minsc = (xtra & 0x40000) ? (xtra & 0xffff) : 0x10000;
    const int endsc = (xtra & 0x20000) ? (xtra & 0xffff) : 0x10000;
    const int oe_del = (uint8_t)(o_del + e_del), ed = (uint8_t)e_del, oe_ins = (uint8_t)(o_ins + e_ins), ei = (uint8_t)e_ins;
    int E[SL], Hp[SL];
#pragma unroll
    for (int j = 0; j < SL; ++j) { E[j] = 0; Hp[j] = 0; }
    // the last thread's state (as lane 0 / the reductions of ksw_u8_pass_t)
    int qe_save = 0, n_b = 0, te = -1, gmax = 0;
    uint32_t last_b = 0;
    bool over = false, done = false;
    // boundary from thread t-1: w = H(8) | F(8) << 8 | rowmax(8) << 16 | tc(3) << 24; key = H << 9 | (511 - pos)
    unsigned w_out = 0, key_out = 0;
    int diag_in = 0;   // H(i-1, my first column - 1)
    int gl = 0;        // the last thread's gmax as of a few rows back (a lower bound of gmax before this row)
    const int nsteps = tlen + 15;
    int tc_next = (t == 0 && tlen > 0) ? tbase(0) : 0;
    for (int s = 0; s < nsteps; ++s) {
        unsigned w_in = __shfl_up_sync(mask, w_out, 1, 16);
        unsigned key_in = __shfl_up_sync(mask, key_out, 1, 16);
        const int i = s - t;
        if (t == 0) {
            w_in = (unsigned)tc_next << 24; key_in = 0;
            if (i + 1 < tlen) tc_next = tbase(i + 1);
        }
        if (i >= 0 && i < tlen) {
            const int tc = (w_in >> 24) & 7;
            int d = diag_in;
            diag_in = w_in & 255;
            int ft = (w_in >> 8) & 255, fs = 0, rmax = (w_in >> 16) & 255;
            const int ft_in = ft;
#define WV_CELL(SC) { \
                int h = SAT ? max(min(d + (SC) + shift, 255) - shift, 0) : max(d + (SC), 0); \
                const int e = E[j]; \
                h = max(h, e); \
                const int hf = max(h, fs), hF = max(h, ft); \
                rmax = max(rmax, hf); \
                E[j] = max(e - ed, hf - oe_del); \
                fs = max(fs - ei, hf - oe_ins); \
                ft = max(ft - ei, hF - oe_ins); \
                d = Hp[j]; Hp[j] = hF; }
            if (tc < 4) {
#pragma unroll
                for (int j = 0; j < SL; ++j) WV_CELL(qv[j] == tc ? sAs : base[j])
            } else {
#pragma unroll
                for (int j = 0; j < SL; ++j) WV_CELL(baseN[j])
            }
#undef WV_CELL
            w_out = (unsigned)Hp[SL - 1] | (unsigned)max(ft, 0) << 8 | (unsigned)rmax << 16 | (unsigned)tc << 24;
            // the row's (largest final H, first position) is read only when the row's maximum beats gmax, and then
            // only a segment whose largest H beats it can hold it: the others pass the incoming key on. A segment's
            // largest final H is at most max(its largest first-sweep H, the F coming in) <= max(rmax, ft_in).
            const int mF = max(rmax, ft_in);
            unsigned key = key_in;
            if (mF > gl) {
#pragma unroll
                for (int j = 0; j < SL; ++j) key = max(key, (unsigned)Hp[j] << 9 | (unsigned)(511 - (j + t * SL)));
            }
            key_out = key;
            if (t == 15 && !done) {
                const int imax = rmax;
                if (imax >= minsc) {
                    if (n_b == 0 || (int)(last_b & 0xfffff) + 1 != i) {
                        last_b = (uint32_t)imax << 20 | (uint32_t)i;
                        if (n_b < bcap) b[n_b] = last_b; else over = true;
                        n_b++;
                    } else if ((int)(last_b >> 20) < imax) {
                        last_b = (uint32_t)imax << 20 | (uint32_t)i;
                        if (n_b <= bcap) b[n_b - 1] = last_b;
                    }
                }
                if (imax > gmax) {
                    gmax = imax; te = i;
                    qe_save = 511 - (int)(key & 511);
                    if (gmax + shift >= 255 || gmax >= endsc) done = true;
                }
            }
        }
        if ((s & 7) == 7) {
            const int g = __shfl_sync(mask, done ? -1 : gmax, 15, 16);
            if (g < 0) break;
            gl = g;
        }
    }
    MRes r; r.score = gmax + shift < 255 ? gmax : 255; r.te = te; r.qe = -1; r.score2 = -1; r.te2 = -1; r.over = false;
    if (t == 15 && r.score != 255) {
        r.qe = te >= 0 ? qe_save : 0;
        if (n_b > 0) {
            int wd = (r.score + qmax - 1) / qmax;
            int low = te - wd, high = te + wd;
            int nb = n_b < bcap ? n_b : bcap;
            for (int k = 0; k < nb; ++k) {
                int e = (int)(b[k] & 0xfffff);
                if ((e < low || e > high) && (int)(b[k] >> 20) > r.score2) { r.score2 = (int)(b[k] >> 20); r.te2 = e; }
            }
            r.over = over;
        }
    }
    // everyone gets the last thread's result
    r.score = __shfl_sync(mask, r.score, 15, 16); r.te = __shfl_sync(mask, r.te, 15, 16); r.qe = __shfl_sync(mask, r.qe, 15, 16);
    r.score2 = __shfl_sync(mask, r.score2, 15, 16); r.te2 = __shfl_sync(mask, r.te2, 15, 16); r.over = __shfl_sync(mask, (int)r.over, 15, 16);
    return r;
}

template <int SL, typename QF, typename TF>
__device__ MRes ksw_u8_wave_s(int qlen, int tlen, QF qbase, TF tbase, int xtra, int t, unsigned mask, uint32_t *b, int bcap) {
    int mn = 127, mx = 0;
    for (int a = 0; a < 25; ++a) { mn = c_mat[a] < mn ? c_mat[a] : mn; mx = c_mat[a] > mx ? c_mat[a] : mx; }
    const int shift = (uint8_t)(256 - (int)(uint8_t)(int8_t)mn);
    const bool sat = SL * 16 * mx + mx + shift > 255;
    return sat ? ksw_u8_wave_t<SL, true>(qlen, tlen, qbase, tbase, xtra, t, mask, b, bcap)
               : ksw_u8_wave_t<SL, false>(qlen, tlen, qbase, tbase, xtra, t, mask, b, bcap);
}

template <typename QF, typename TF>
__device__ MRes ksw_u8_wave(int qlen, int tlen, QF qbase, TF tbase, int xtra, int t, unsigned mask, uint32_t *b, int bcap) {
    switch ((qlen + 15) / 16) {
#define WS_CASE(n) case n: return ksw_u8_wave_s<n>(qlen, tlen, qbase, tbase, xtra, t, mask, b, bcap);
        WS_CASE(1) WS_CASE(2) WS_CASE(3) WS_CASE(4) WS_CASE(5) WS_CASE(6) WS_CASE(7) WS_CASE(8)
        WS_CASE(9) WS_CASE(10) WS_CASE(11) WS_CASE(12) WS_CASE(13) WS_CASE(14) WS_CASE(15) WS_CASE(16)
#undef WS_CASE
    }
    MRes r; r.score = 0; r.te = -1; r.qe = -1; r.score2 = -1; r.te2 = -1; r.over = true;
    return r;
}

__global__ void __launch_bounds__(MSW_BS, MSW_WAVE_MINB) k_matesw_wave(const MJob *jobs, int n, const uint8_t *reads, const int64_t *off, const uint8_t *pac, int64_t l_pac, MOut *out, uint32_t *bscratch) {
    int gt = blockIdx.x * blockDim.x + threadIdx.x;
    int job = gt >> 4, lane = threadIdx.x & 15;
    unsigned mask = 0xffffu << (threadIdx.x & 16);
    if (job >= n) return;
    const MJob jb = jobs[job];
    const uint8_t *rd = reads + off[jb.read];
    const int lq = jb.qlen;
    auto qf = [&](int k) -> int { int c = jb.rev ? rd[lq - 1 - k] : rd[k]; return jb.rev ? (c < 4 ? 3 - c : 4) : c; };
    auto tf = [&](int i) -> int { return ref_base(pac, l_pac, jb.rb + i); };
    uint32_t *b = bscratch + jb.boff;
    const int bcap = jb.tlen + 1;
    MRes r = ksw_u8_wave(lq, jb.tlen, qf, tf, jb.xtra, lane, mask, b, bcap);
    MOut o; o.score = r.score; o.te = r.te; o.qe = r.qe; o.score2 = r.score2; o.te2 = r.te2; o.tb = -1; o.qb = -1; o.ok = 1;
    if (r.over) o.ok = 0;
    const int xtra = jb.xtra;
    if (!((xtra & 0x80000) == 0 || ((xtra & 0x40000) && r.score < (xtra & 0xffff)))) {
        const int qe = r.qe, te = r.te;
        auto qf2 = [&](int k) -> int { return qf(qe - k); };
        auto tf2 = [&](int i) -> int { return i <= te ? tf(te - i) : tf(i); };
        MRes rr = ksw_u8_wave(qe + 1, jb.tlen, qf2, tf2, 0x20000 | r.score, lane, mask, b, bcap);
        if (r.score == rr.score) { o.tb = r.te - rr.te; o.qb = r.qe - rr.qe; }
    }
    if (lane == 15) out[job] = o;
}


static MJob *d_mjobs = 0; static size_t cap_mjobs = 0;
static uint8_t *d_reads2 = 0; static size_t cap_reads2 = 0;
static int64_t *d_off2 = 0;   static size_t cap_off2 = 0;

// The SAM stage's batch (its own buffer: the alignment stage uploads the next batch meanwhile)
extern "C" int aln_gpu_set_reads2(const uint8_t *reads, int64_t reads_bytes, const int64_t *off, int n_reads) {
    if (grow(&d_reads2, &cap_reads2, (size_t)reads_bytes + 1) || grow(&d_off2, &cap_off2, (size_t)n_reads + 1)) return -1;
    CK(cudaMemcpy(d_reads2, reads, reads_bytes, cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_off2, off, (n_reads + 1) * sizeof(int64_t), cudaMemcpyHostToDevice));
    return 0;
}
static MOut *d_mout = 0;  static size_t cap_mout = 0;
static uint32_t *d_bscr = 0; static size_t cap_bscr = 0;
static double t_msw = 0;
static thread_local double t_msw_call[3] = {0, 0, 0};

// Mate-rescue alignments of the reads set by aln_gpu_set_reads (u8 jobs only); ok = 0: the CPU redoes the job
// jobs[k].boff: offset of job k's per-row maxima in the scratch (sum of tlen + 1 before it); `scratch` its size
extern "C" int aln_gpu_matesw(const MJob *jobs, int n, int64_t scratch, MOut *out) {
    if (n == 0) return 0;
    SERIAL;
    if (grow(&d_mjobs, &cap_mjobs, (size_t)n) || grow(&d_mout, &cap_mout, (size_t)n) || grow(&d_bscr, &cap_bscr, (size_t)scratch)) return -1;
    cudaEvent_t e0 = tev()[0], e1 = tev()[1]; cudaEventRecord(e0);
    CK(cudaMemcpy(d_mjobs, jobs, (size_t)n * sizeof(MJob), cudaMemcpyHostToDevice));
    int bs = MSW_BS;  // 16 threads per job
    int64_t nb = ((int64_t)n * 16 + bs - 1) / bs;
    cudaEvent_t ek0 = tev()[2], ek1 = tev()[3]; cudaEventRecord(ek0);
    if (g_msw_wave) k_matesw_wave<<<(unsigned)nb, bs>>>(d_mjobs, n, d_rs[2], d_os[2], d_pac, g_l_pac, d_mout, d_bscr);
    else switch (g_msw_minb) {
    case 12: k_matesw<12><<<(unsigned)nb, bs>>>(d_mjobs, n, d_rs[2], d_os[2], d_pac, g_l_pac, d_mout, d_bscr); break;
    case 24: k_matesw<24><<<(unsigned)nb, bs>>>(d_mjobs, n, d_rs[2], d_os[2], d_pac, g_l_pac, d_mout, d_bscr); break;
    default: k_matesw<18><<<(unsigned)nb, bs>>>(d_mjobs, n, d_rs[2], d_os[2], d_pac, g_l_pac, d_mout, d_bscr);
    }
    CK(cudaGetLastError());
    cudaEventRecord(ek1);
    CK(cudaMemcpy(out, d_mout, (size_t)n * sizeof(MOut), cudaMemcpyDeviceToHost));
    cudaEventRecord(e1); cudaEventSynchronize(e1);
    float ms_; cudaEventElapsedTime(&ms_, e0, e1); t_msw += ms_ * 1e-3;
    add_kern(3, ek0, ek1);
    float a, b, c;
    cudaEventElapsedTime(&a, e0, ek0); cudaEventElapsedTime(&b, ek0, ek1); cudaEventElapsedTime(&c, ek1, e1);
    t_msw_call[0] += a * 1e-3; t_msw_call[1] += b * 1e-3; t_msw_call[2] += c * 1e-3;
    return 0;
}
// ALIGN_RS_TRACE: this thread's mate-SW device time since the last call (s): upload, kernel, download
extern "C" void aln_gpu_msw_take(double *out) {
    for (int i = 0; i < 3; ++i) { out[i] = t_msw_call[i]; t_msw_call[i] = 0; }
}

extern "C" void aln_gpu_report3() {
    fprintf(stderr, "[align-rs-gpu] mate SW %.2f s (incl. transfers)\n", t_msw);
    fprintf(stderr, "[align-rs-gpu] kernels only: extension %.2f s, seeding %.2f s, SA %.2f s, mate SW %.2f s%s\n", t_kern[0], t_kern[1], t_kern[2], t_kern[3],
            g_serial() ? " (ALIGN_RS_GPU_SERIAL: one at a time)" : " (concurrent kernels share the device)");
}

// The device's settings (environment first) and whether this build has kernels for it (-1 if not). Defaults from
// the device: extension blocks per SM so that the resident threads' DP rows (~170 bytes each in use) fit the
// L1 (L1 + shared per SM ~ the shared-memory maximum + 28 KB on sm_75 .. sm_120); the mate-SW register cap for
// three quarters of the SM's thread slots. Tuned on an RTX 5070 Ti (sm_120: 4 / 6 / 18); elsewhere a starting
// point (docs/notes/align-rs.md: GPU settings).
static int env_int(const char *k, int def) { const char *e = getenv(k); return e && *e ? atoi(e) : def; }
static int gpu_settings(int dev) {
    cudaDeviceProp pr;
    CK(cudaGetDeviceProperties(&pr, dev));
    cudaFuncAttributes fa;
    if (cudaFuncGetAttributes(&fa, k_extend) != cudaSuccess) {
        cudaGetLastError();
        fprintf(stderr, "[align-rs-gpu] %s (sm_%d%d): this build has no kernels for it (cuda/Makefile ARCH)\n", pr.name, pr.major, pr.minor);
        return -1;
    }
    g_n_sm = pr.multiProcessorCount;
    const int l1_kb = (int)(pr.sharedMemPerMultiprocessor >> 10) + 28;
    g_ext_bpsm = env_int("ALIGN_RS_EXT_BPSM", l1_kb / 21 > 2 ? l1_kb / 21 : 2);
    int want = pr.maxThreadsPerMultiProcessor * 3 / 4 / MSW_BS;  // 1536 -> 18, 2048 -> 24, 1024 -> 12
    int minb = want < 15 ? 12 : want < 21 ? 18 : 24;
    minb = env_int("ALIGN_RS_MSW_MINB", minb);
    if (minb != 12 && minb != 18 && minb != 24) { fprintf(stderr, "[align-rs-gpu] ALIGN_RS_MSW_MINB=%d: one of 12, 18, 24; using 18\n", minb); minb = 18; }
    g_msw_minb = minb;
    { const char *m = getenv("ALIGN_RS_MSW"); g_msw_wave = !(m && !strcmp(m, "striped")); }
    int g = env_int("ALIGN_RS_SMEM_G", 4);
    if (g != 1 && g != 2 && g != 4 && g != 8 && g != 16 && g != 32) { fprintf(stderr, "[align-rs-gpu] ALIGN_RS_SMEM_G=%d: one of 1, 2, 4, 8, 16, 32; using 4\n", g); g = 4; }
    g_smem_g = g;
    cudaFuncAttributes fm;
    int regs = 0;
    if ((g_msw_wave ? cudaFuncGetAttributes(&fm, k_matesw_wave) : minb == 12 ? cudaFuncGetAttributes(&fm, k_matesw<12>) : minb == 24 ? cudaFuncGetAttributes(&fm, k_matesw<24>) : cudaFuncGetAttributes(&fm, k_matesw<18>)) == cudaSuccess) regs = fm.numRegs;
    fprintf(stderr, "[align-rs-gpu] %s (sm_%d%d, %d SMs, ~%d KB L1/SM, %d threads/SM): seeding %d lanes per read, extension %d blocks/SM, mate SW %s %d registers (ALIGN_RS_SMEM_G, ALIGN_RS_EXT_BPSM, ALIGN_RS_MSW, ALIGN_RS_MSW_MINB)\n",
            pr.name, pr.major, pr.minor, g_n_sm, l1_kb, pr.maxThreadsPerMultiProcessor, g_smem_g, g_ext_bpsm, g_msw_wave ? "wavefront" : "striped", regs);
    return 0;
}

// ---- align-rs tune: settings, kernel times, device identity
// a setting < 0 is kept; invalid values are refused (-1)
extern "C" int aln_gpu_set_settings(int smem_g, int ext_bpsm, int msw_minb) {
    if (smem_g >= 0 && smem_g != 1 && smem_g != 2 && smem_g != 4 && smem_g != 8 && smem_g != 16 && smem_g != 32) return -1;
    if (msw_minb >= 0 && msw_minb != 12 && msw_minb != 18 && msw_minb != 24) return -1;
    if (smem_g >= 0) g_smem_g = smem_g;
    if (ext_bpsm >= 0) g_ext_bpsm = ext_bpsm;
    if (msw_minb >= 0) g_msw_minb = msw_minb;
    return 0;
}
// 1: mate SW runs k_matesw_wave, whose register cap is fixed (ALIGN_RS_MSW_MINB applies to k_matesw only)
extern "C" int aln_gpu_msw_wave() { return g_msw_wave; }
extern "C" void aln_gpu_get_settings(int *out) { out[0] = g_smem_g; out[1] = g_ext_bpsm; out[2] = g_msw_minb; }
// seconds of kernel time so far: extension, seeding, SA, mate SW (exact when one call runs at a time)
extern "C" void aln_gpu_kernel_times(double *out) {
    cudaDeviceSynchronize();
    std::lock_guard<std::mutex> g(g_tk_mu);
    for (int i = 0; i < 4; ++i) out[i] = t_kern[i];
}
extern "C" int aln_gpu_device(char *name, int cap, int *major, int *minor, int *driver) {
    int dev = 0;
    cudaDeviceProp pr;
    CK(cudaGetDevice(&dev));
    CK(cudaGetDeviceProperties(&pr, dev));
    snprintf(name, cap, "%s", pr.name);
    *major = pr.major; *minor = pr.minor;
    CK(cudaDriverGetVersion(driver));
    return 0;
}
