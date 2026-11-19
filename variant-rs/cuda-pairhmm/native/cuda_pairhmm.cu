// CUDA PairHMM for GATK HaplotypeCaller.
//
// Computes what Intel GKL's AVX PairHMM (com.intel.gkl.pairhmm.IntelPairHmm) computes for GATK's
// HaplotypeCaller / Mutect2 (all reads x all haplotypes of active regions), bit for bit: float first, double
// fallback when the float result underflows below 1e-28. variant-rs calls it through the C ABI at the end.
//
// Parallelization strategy
// ------------------------
//  * Inter-task: one CUDA warp per (read, haplotype) pair. A batch of R reads and
//    H haplotypes launches R*H warps.
//  * Intra-task: the DP matrix (rows = read bases, cols = haplotype bases) is cut
//    into horizontal stripes of 32 rows. Lane t of the warp owns row (32*s + t + 1).
//    The warp sweeps anti-diagonals: at step k, lane t computes column k - t + 1.
//    Dependencies:
//        M/X/Y(row-1, col-1)  -> lane t-1, two steps ago   (__shfl_up, delayed)
//        M/X  (row-1, col  )  -> lane t-1, previous step   (__shfl_up)
//        M/Y  (row,   col-1)  -> this lane, previous step  (registers)
//    Lane 31 writes its row into a shared-memory line buffer; lane 0 of the next
//    stripe reads it back as its "row-1". This is exactly what GKL's AVX code does
//    with 8/16 SIMD lanes and its shiftOut arrays - the GPU just uses 32 lanes and
//    runs thousands of pairs at once.
//
// Recurrence (GATK LoglessPairHMM; i = read row, j = haplotype col):
//   M[i][j] = prior(i,j) * (M[i-1][j-1]*mm[i] + (X[i-1][j-1] + Y[i-1][j-1])*gapm[i])
//   X[i][j] = M[i-1][j]*mx[i] + X[i-1][j]*xx[i]          (insertion)
//   Y[i][j] = M[i][j-1]*my[i] + Y[i][j-1]*yy[i]          (deletion)
//   init:  Y[0][j] = INITIAL_CONSTANT / haplen, everything else in row/col 0 is 0
//   result = log10( sum_j M[R][j] + X[R][j] ) - log10(INITIAL_CONSTANT)
//
// Probability tables replicate GKL Context.h (MIT License, Copyright (c) 2016-2024
// Intel Corporation) so results match the AVX implementation closely.

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <vector>
#include <chrono>
#include <algorithm>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <cstdlib>

#define CUDA_CHECK(call)                                                          \
  do {                                                                            \
    cudaError_t _e = (call);                                                      \
    if (_e != cudaSuccess) {                                                      \
      fprintf(stderr, "[cuda-pairhmm] CUDA error %s at %s:%d: %s\n", #call,      \
              __FILE__, __LINE__, cudaGetErrorString(_e));                        \
      return false;                                                               \
    }                                                                             \
  } while (0)

static constexpr int   WARP          = 32;
static constexpr int   WARPS_PER_BLK = 4;
static constexpr float MIN_ACCEPTED  = 1e-28f;
static constexpr int   NQ            = 128;                 // quals are masked with & 127
static constexpr int   MM_TABLE_SIZE = (NQ * (NQ + 1)) / 2; // triangular (max,min) table

// ---------------------------------------------------------------------------
// Host-side probability tables (bit-compatible with GKL's Context<NUMBER>)
// ---------------------------------------------------------------------------
#define MAX_QUAL 254
#define MAX_JACOBIAN_TOLERANCE 8.0
#define JACOBIAN_LOG_TABLE_STEP 0.0001
#define JACOBIAN_LOG_TABLE_INV_STEP (1.0 / JACOBIAN_LOG_TABLE_STEP)
#define JACOBIAN_LOG_TABLE_SIZE ((int)(MAX_JACOBIAN_TOLERANCE / JACOBIAN_LOG_TABLE_STEP) + 1)

template <class T> struct HostTables {
  T ph2pr[NQ];
  T mm[MM_TABLE_SIZE];
  T initialConstant, log10InitialConstant;
  std::vector<T> jacobian;

  int fastRound(T d) { return (d > (T)0.0) ? (int)(d + (T)0.5) : (int)(d - (T)0.5); }
  T approxLog10SumLog10(T small, T big) {
    if (small > big) std::swap(small, big);
    if (std::isinf(small) || std::isinf(big)) return big;
    T diff = big - small;
    if (diff >= (T)MAX_JACOBIAN_TOLERANCE) return big;
    int ind = fastRound((T)(diff * ((T)JACOBIAN_LOG_TABLE_INV_STEP)));
    return big + jacobian[ind];
  }
  void init(bool isFloat) {
    jacobian.resize(JACOBIAN_LOG_TABLE_SIZE);
    for (int k = 0; k < JACOBIAN_LOG_TABLE_SIZE; k++)
      jacobian[k] = (T)(log10(1.0 + pow(10.0, -((double)k) * JACOBIAN_LOG_TABLE_STEP)));
    const double INV_LN10 = 0.434294;
    std::vector<T> full(((MAX_QUAL + 1) * (MAX_QUAL + 2)) >> 1);
    for (int i = 0, offset = 0; i <= MAX_QUAL; offset += ++i)
      for (int j = 0; j <= i; j++) {
        double log10Sum = approxLog10SumLog10((T)-0.1 * (T)i, (T)-0.1 * (T)j);
        double mm10 = log1p(-std::min(1.0, pow(10, log10Sum))) * INV_LN10;
        full[offset + j] = (T)(pow(10, mm10));
      }
    std::copy(full.begin(), full.begin() + MM_TABLE_SIZE, mm);
    for (int x = 0; x < NQ; x++)
      ph2pr[x] = isFloat ? (T)powf(10.f, -((float)x) / 10.f) : (T)pow(10.0, -((double)x) / 10.0);
    if (isFloat) { initialConstant = (T)ldexpf(1.f, 120); log10InitialConstant = (T)log10f((float)initialConstant); }
    else         { initialConstant = (T)ldexp(1.0, 1020); log10InitialConstant = (T)log10((double)initialConstant); }
  }
};

// Device copies of the tables
__constant__ float  c_ph2pr_f[NQ];
__constant__ double c_ph2pr_d[NQ];
__constant__ float  c_init_f, c_log10init_f;
__constant__ double c_init_d, c_log10init_d;
__device__ float  g_mm_f[MM_TABLE_SIZE];
__device__ double g_mm_d[MM_TABLE_SIZE];

template <class T> struct DevTables;
template <> struct DevTables<float> {
  __device__ static float ph2pr(int q) { return c_ph2pr_f[q]; }
  __device__ static float mm(int idx) { return g_mm_f[idx]; }
  __device__ static float init() { return c_init_f; }
};
template <> struct DevTables<double> {
  __device__ static double ph2pr(int q) { return c_ph2pr_d[q]; }
  __device__ static double mm(int idx) { return g_mm_d[idx]; }
  __device__ static double init() { return c_init_d; }
};

// ---------------------------------------------------------------------------
// Batch layout on the device
// ---------------------------------------------------------------------------
struct Batch {
  // A launch may contain several GATK batches ("groups", one per active region).
  // Group g owns pairs [pairStart[g], pairStart[g+1]) laid out read-major, exactly
  // like the likelihood array GATK expects: pair = pairStart[g] + r*nHaps[g] + h.
  int nGroups;
  const int* pairStart; // [nGroups+1]
  const int* readStart; // [nGroups]   first read index of the group
  const int* hapStart;  // [nGroups]   first hap index of the group
  const int* nHaps;     // [nGroups]
  const int* readOff;   // [nReads]  offset into read arrays
  const int* readLen;   // [nReads]
  const int* hapOff;    // [nHaps]
  const int* hapLen;    // [nHaps]
  const uint8_t* rs;    // read bases, converted to 0..4 (A,C,T,G,N)
  const uint8_t* q;     // base quals
  const uint8_t* ins;   // insertion GOP
  const uint8_t* del;   // deletion GOP
  const uint8_t* gcp;   // gap continuation
  const uint8_t* hap;   // haplotype bases, converted to 0..4
};

// One warp computes one (read, hap) pair. Returns the (scaled) likelihood
// sum_j M[R][j] + X[R][j] in every lane.
// Map a global pair index to (global read index, global hap index).
__device__ __forceinline__ void locate_pair(const Batch& b, int pair, int& r, int& h) {
  int lo = 0, hi = b.nGroups - 1;  // find last g with pairStart[g] <= pair
  while (lo < hi) {
    const int mid = (lo + hi + 1) >> 1;
    if (b.pairStart[mid] <= pair) lo = mid; else hi = mid - 1;
  }
  const int local = pair - b.pairStart[lo];
  const int nh = b.nHaps[lo];
  r = b.readStart[lo] + local / nh;
  h = b.hapStart[lo] + local % nh;
}

// Flush-to-zero. GKL sets MXCSR.FTZ, which flushes denormal results of float AND double operations.
// nvcc -ftz=true covers float only, so every double result is flushed explicitly.
__device__ __forceinline__ float  FZ(float x)  { return x; }
__device__ __forceinline__ double FZ(double x) { return fabs(x) < 2.2250738585072014e-308 ? copysign(0.0, x) : x; }
// Explicit fused multiply-add (the file is compiled with -fmad=false, so nothing else is contracted).
__device__ __forceinline__ float  FMA(float a, float b, float c)    { return __fmaf_rn(a, b, c); }
__device__ __forceinline__ double FMA(double a, double b, double c) { return FZ(__fma_rn(a, b, c)); }
__device__ __forceinline__ float  MUL(float a, float b)  { return __fmul_rn(a, b); }
__device__ __forceinline__ double MUL(double a, double b) { return FZ(__dmul_rn(a, b)); }
__device__ __forceinline__ float  ADD(float a, float b)  { return __fadd_rn(a, b); }
__device__ __forceinline__ double ADD(double a, double b) { return FZ(__dadd_rn(a, b)); }
// Explicit IEEE round-to-nearest division and subtraction (GKL divides with divps/divpd).
__device__ __forceinline__ float  DIV(float a, float b)   { return __fdiv_rn(a, b); }
__device__ __forceinline__ double DIV(double a, double b) { return FZ(__ddiv_rn(a, b)); }
__device__ __forceinline__ float  SUB(float a, float b)   { return __fsub_rn(a, b); }
__device__ __forceinline__ double SUB(double a, double b) { return __dsub_rn(a, b); }

template <class T, bool DBG = false>
__device__ T pairhmm_warp(const Batch& b, int r, int h, T* bufM, T* bufX, T* bufY, int lane, T* dbg = nullptr) {
  const unsigned FULL = 0xffffffffu;
  const int R = b.readLen[r];
  const int H = b.hapLen[h];
  const uint8_t* hap = b.hap + b.hapOff[h];
  const int ro = b.readOff[r];
  const T initY = DIV(DevTables<T>::init(), (T)H);

  const int nStripes = (R + WARP - 1) / WARP;
  const int lastLane = (R - 1) % WARP;
  T sumM = 0, sumX = 0;   // GKL accumulates the last row's M and X separately

  for (int s = 0; s < nStripes; ++s) {
    const int row = s * WARP + lane + 1;  // 1-based DP row this lane owns
    const bool lastStripe = (s == nStripes - 1);

    // Per-row transition / emission probabilities (GKL initializeVectors)
    T pMM = 0, pGAPM = 0, pMX = 0, pXX = 0, pMY = 0, pYY = 0, dist1 = 0, dist3 = 0;
    int rbase = 255;
    if (row <= R) {
      const int p = ro + row - 1;
      const int qi = b.ins[p] & 127, qd = b.del[p] & 127, qc = b.gcp[p] & 127, qq = b.q[p] & 127;
      const int mx = max(qi, qd), mn = min(qi, qd);
      pMM   = DevTables<T>::mm(((mx * (mx + 1)) >> 1) + mn);
      pGAPM = SUB((T)1.0, DevTables<T>::ph2pr(qc));
      pMX   = DevTables<T>::ph2pr(qi);
      pXX   = DevTables<T>::ph2pr(qc);
      pMY   = DevTables<T>::ph2pr(qd);
      pYY   = DevTables<T>::ph2pr(qc);
      const T distm = DevTables<T>::ph2pr(qq);
      dist1 = SUB((T)1.0, distm);   // prior if bases match
      dist3 = DIV(distm, (T)3.0);   // prior if they mismatch
      rbase = b.rs[p];
    }

    // Own values at (row, col-1); start as column 0 = 0
    T M1 = 0, X1 = 0, Y1 = 0;
    // Values at (row-1, col-1) for the next cell. For lane 0 of stripe 0 that is
    // DP cell (0,0) whose Y is initY (free deletions before the read starts).
    T dM = 0, dX = 0, dY = (s == 0 && lane == 0) ? initY : (T)0;

    const int nSteps = H + WARP - 1;
    for (int k = 0; k < nSteps; ++k) {
      const int c = k - lane + 1;  // haplotype column (1-based) this lane works on

      // Neighbour from the row above: lane-1's result from the previous step = (row-1, c)
      T uM = __shfl_up_sync(FULL, M1, 1);
      T uX = __shfl_up_sync(FULL, X1, 1);
      T uY = __shfl_up_sync(FULL, Y1, 1);
      if (lane == 0) {
        if (s == 0) { uM = 0; uX = 0; uY = initY; }          // DP row 0
        else if (c <= H) { uM = bufM[c]; uX = bufX[c]; uY = bufY[c]; }  // last row of prev stripe
      }

      T nM = 0, nX = 0, nY = 0;
      if (c >= 1 && c <= H) {
        const int hb = hap[c - 1];
        const bool match = (rbase == hb) || rbase == 4 || hb == 4;  // 'N' matches anything
        const T prior = match ? dist1 : dist3;
        // Same operation order and FMA contraction as GKL's AVX-512 build (gcc -mfma):
        nM = MUL(FMA(dY, pGAPM, FMA(dX, pGAPM, MUL(dM, pMM))), prior);
        nX = FMA(uX, pXX, MUL(uM, pMX));
        nY = FMA(Y1, pYY, MUL(M1, pMY));
        if (lastStripe && lane == lastLane) { sumM = ADD(sumM, nM); sumX = ADD(sumX, nX); }
        if (DBG && row <= R) { const size_t o = 3 * ((size_t)row * (H + 1) + c); dbg[o] = nM; dbg[o + 1] = nX; dbg[o + 2] = nY; }
        if (lane == WARP - 1 && !lastStripe) { bufM[c] = nM; bufX[c] = nX; bufY[c] = nY; }
      }
      dM = uM; dX = uX; dY = uY;   // (row-1, c) becomes the diagonal for column c+1
      M1 = nM; X1 = nX; Y1 = nY;
      __syncwarp();
    }
    __syncwarp();
  }
  return __shfl_sync(FULL, ADD(sumM, sumX), lastLane);
}

// Pass 1: float for every pair. Pairs whose float result underflows are queued for pass 2.
__global__ void pairhmm_float_kernel(Batch b, int nPairs, double* out, unsigned char* path, int* redo, int* nRedo, int bufStride) {
  extern __shared__ unsigned char smem[];
  const int warpInBlk = threadIdx.x / WARP, lane = threadIdx.x % WARP;
  const int pair = blockIdx.x * (blockDim.x / WARP) + warpInBlk;  // warps per block: set per launch
  if (pair >= nPairs) return;
  float* buf = reinterpret_cast<float*>(smem) + warpInBlk * 3 * bufStride;
  int r, h; locate_pair(b, pair, r, h);
  if (b.readLen[r] == 0) { if (lane == 0) { out[pair] = 0.0; path[pair] = 2; } return; }
  const float res = pairhmm_warp<float>(b, r, h, buf, buf + bufStride, buf + 2 * bufStride, lane);
  if (lane == 0) {
    if (res < MIN_ACCEPTED) { redo[atomicAdd(nRedo, 1)] = pair; }
    else { out[pair] = (double)res; path[pair] = 0; }   // raw float sum; host takes log10f
  }
}

// Pass 2: double precision for the (rare) underflowing pairs, or for everything if requested.
// nPairsDev (if non-null) holds the device-side count of queued pairs, so the host does
// not need to synchronize between pass 1 and pass 2.
__global__ void pairhmm_double_kernel(Batch b, const int* pairs, int nPairs, const int* nPairsDev, double* out, unsigned char* path, int bufStride) {
  extern __shared__ unsigned char smem[];
  const int warpInBlk = threadIdx.x / WARP, lane = threadIdx.x % WARP;
  const int i = blockIdx.x * (blockDim.x / WARP) + warpInBlk;
  if (i >= (nPairsDev ? *nPairsDev : nPairs)) return;
  const int pair = pairs ? pairs[i] : i;
  double* buf = reinterpret_cast<double*>(smem) + warpInBlk * 3 * bufStride;
  int r, h; locate_pair(b, pair, r, h);
  if (b.readLen[r] == 0) { if (lane == 0) { out[pair] = 0.0; path[pair] = 2; } return; }
  const double res = pairhmm_warp<double>(b, r, h, buf, buf + bufStride, buf + 2 * bufStride, lane);
  if (lane == 0) { out[pair] = res; path[pair] = 1; }   // raw double sum; host takes log10
}

// ---------------------------------------------------------------------------
// The device's one-time setup (tables, kernel attributes), shared by the engines.
// ---------------------------------------------------------------------------
struct Device {
  bool ready = false;
  bool useDouble = false;
  float log10InitF = 0; double log10InitD = 0;
  uint8_t conv[256];
  int maxSmem = 48 * 1024;  // opt-in dynamic shared memory per block (set in init)

  bool init(bool dbl) {
    useDouble = dbl;
    int n = 0;
    if (cudaGetDeviceCount(&n) != cudaSuccess || n == 0) return false;
    // block (instead of spinning) while waiting for the GPU: the host threads have CPU work to do
    cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync);
    CUDA_CHECK(cudaSetDevice(0));
    CUDA_CHECK(cudaDeviceGetAttribute(&maxSmem, cudaDevAttrMaxSharedMemoryPerBlockOptin, 0));
    CUDA_CHECK(cudaFuncSetAttribute(pairhmm_float_kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, maxSmem));
    CUDA_CHECK(cudaFuncSetAttribute(pairhmm_double_kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, maxSmem));

    static HostTables<float> tf; static HostTables<double> td;
    tf.init(true); td.init(false);
    CUDA_CHECK(cudaMemcpyToSymbol(c_ph2pr_f, tf.ph2pr, sizeof(tf.ph2pr)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_ph2pr_d, td.ph2pr, sizeof(td.ph2pr)));
    CUDA_CHECK(cudaMemcpyToSymbol(g_mm_f, tf.mm, sizeof(tf.mm)));
    CUDA_CHECK(cudaMemcpyToSymbol(g_mm_d, td.mm, sizeof(td.mm)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_init_f, &tf.initialConstant, sizeof(float)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_log10init_f, &tf.log10InitialConstant, sizeof(float)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_init_d, &td.initialConstant, sizeof(double)));
    CUDA_CHECK(cudaMemcpyToSymbol(c_log10init_d, &td.log10InitialConstant, sizeof(double)));
    log10InitF = tf.log10InitialConstant; log10InitD = td.log10InitialConstant;

    memset(conv, 0, sizeof(conv));  // GKL ConvertChar: unknown -> 0 ('A')
    conv['A'] = 0; conv['C'] = 1; conv['T'] = 2; conv['G'] = 3; conv['N'] = 4;
    ready = true;
    return true;
  }

  // Longest haplotype one warp per block can handle (double pass: 3 rows of H + 1 doubles).
  int max_hap_len() const { return maxSmem / (3 * (int)sizeof(double)) - 1; }
};

static Device g_dev;

// ---------------------------------------------------------------------------
// Host engine: packs a batch into its pinned buffer, ships it, launches, returns. One caller at a time per engine;
// several engines (their own stream and buffers) let one caller pack while another's kernels run.
// ---------------------------------------------------------------------------
struct Engine {
  cudaStream_t stream{};
  cudaEvent_t evStart{}, evH2D{}, evKernel{}, evEnd{};

  // pinned host staging + device mirror (same layout)
  uint8_t* hIn = nullptr; uint8_t* dIn = nullptr; size_t inCap = 0;
  double* hOut = nullptr; double* dOut = nullptr; size_t outCap = 0;
  int* dRedo = nullptr; int* dNRedo = nullptr; int* hNRedo = nullptr; size_t redoCap = 0;
  unsigned char* dPath = nullptr; unsigned char* hPath = nullptr;

  // statistics (this engine's; one caller at a time)
  long long nsPack = 0, nsGpu = 0, nBatches = 0, nLaunches = 0, nPairs = 0, nCells = 0, nDouble = 0;
  double msH2D = 0, msKernel = 0, msD2H = 0;

  bool create() {
    CUDA_CHECK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    CUDA_CHECK(cudaEventCreate(&evStart)); CUDA_CHECK(cudaEventCreate(&evH2D));
    CUDA_CHECK(cudaEventCreate(&evKernel)); CUDA_CHECK(cudaEventCreate(&evEnd));
    CUDA_CHECK(cudaMallocHost(&hNRedo, sizeof(int)));
    CUDA_CHECK(cudaMalloc(&dNRedo, sizeof(int)));
    return true;
  }

  bool ensure(size_t inBytes, size_t nOut) {
    if (inBytes > inCap) {
      if (hIn) { cudaFreeHost(hIn); cudaFree(dIn); }
      inCap = std::max(inBytes, inCap * 2);
      CUDA_CHECK(cudaMallocHost(&hIn, inCap));
      CUDA_CHECK(cudaMalloc(&dIn, inCap));
    }
    if (nOut > outCap) {
      if (hOut) { cudaFreeHost(hOut); cudaFree(dOut); cudaFree(dRedo); cudaFreeHost(hPath); cudaFree(dPath); }
      outCap = std::max(nOut, outCap * 2);
      CUDA_CHECK(cudaMallocHost(&hOut, outCap * sizeof(double)));
      CUDA_CHECK(cudaMalloc(&dOut, outCap * sizeof(double)));
      CUDA_CHECK(cudaMalloc(&dRedo, outCap * sizeof(int)));
      CUDA_CHECK(cudaMallocHost(&hPath, outCap));
      CUDA_CHECK(cudaMalloc(&dPath, outCap));
    }
    return true;
  }

  static size_t align(size_t x) { return (x + 15) & ~size_t(15); }

  // One GATK batch (one active region): R reads x H haplotypes.
  struct HostGroup {
    std::vector<const int8_t*> rb, rq, ri, rd, rc, hb;
    std::vector<int> rl, hl;
    double* out;
  };

  // Computes all groups in a single GPU launch.
  bool compute(const std::vector<HostGroup>& groups) {
    auto t0 = std::chrono::steady_clock::now();
    const int G = (int)groups.size();
    long long totR = 0, totH = 0, cells = 0; int maxH = 0, nR = 0, nH = 0, nPairsB = 0;
    for (const auto& g : groups) {
      long long gr = 0, gh = 0;
      for (int x : g.rl) gr += x;
      for (int x : g.hl) { gh += x; maxH = std::max(maxH, x); }
      totR += gr; totH += gh; cells += gr * gh;
      nR += (int)g.rl.size(); nH += (int)g.hl.size(); nPairsB += (int)(g.rl.size() * g.hl.size());
    }
    if (nPairsB == 0) return true;

    // Layout: pairStart | readStart | hapStart | nHaps | readOff | readLen | hapOff | hapLen | rs | q | ins | del | gcp | hap
    size_t o_pairStart = 0;
    size_t o_readStart = align(o_pairStart + (G + 1) * sizeof(int));
    size_t o_hapStart  = align(o_readStart + G * sizeof(int));
    size_t o_nHaps     = align(o_hapStart + G * sizeof(int));
    size_t o_readOff   = align(o_nHaps + G * sizeof(int));
    size_t o_readLen = align(o_readOff + nR * sizeof(int));
    size_t o_hapOff  = align(o_readLen + nR * sizeof(int));
    size_t o_hapLen  = align(o_hapOff + nH * sizeof(int));
    size_t o_rs      = align(o_hapLen + nH * sizeof(int));
    size_t o_q = align(o_rs + totR), o_i = align(o_q + totR), o_d = align(o_i + totR), o_c = align(o_d + totR);
    size_t o_hap = align(o_c + totR);
    size_t total = align(o_hap + totH);
    if (!ensure(total, nPairsB)) return false;

    int* pairStart = (int*)(hIn + o_pairStart); int* readStart = (int*)(hIn + o_readStart);
    int* hapStart = (int*)(hIn + o_hapStart);   int* nHapsA = (int*)(hIn + o_nHaps);
    int* readOff = (int*)(hIn + o_readOff); int* readLen = (int*)(hIn + o_readLen);
    int* hapOff = (int*)(hIn + o_hapOff);   int* hapLen = (int*)(hIn + o_hapLen);
    long long roff = 0, hoff = 0; int ri_ = 0, hi_ = 0, ps = 0;
    for (int gi = 0; gi < G; ++gi) {
      const auto& g = groups[gi];
      pairStart[gi] = ps; readStart[gi] = ri_; hapStart[gi] = hi_; nHapsA[gi] = (int)g.hl.size();
      ps += (int)(g.rl.size() * g.hl.size());
      for (size_t r = 0; r < g.rl.size(); ++r, ++ri_) {
        const int L = g.rl[r];
        readOff[ri_] = (int)roff; readLen[ri_] = L;
        uint8_t* dst = hIn + o_rs + roff;
        for (int k = 0; k < L; ++k) dst[k] = g_dev.conv[(uint8_t)g.rb[r][k]];
        memcpy(hIn + o_q + roff, g.rq[r], L); memcpy(hIn + o_i + roff, g.ri[r], L);
        memcpy(hIn + o_d + roff, g.rd[r], L); memcpy(hIn + o_c + roff, g.rc[r], L);
        roff += L;
      }
      for (size_t h = 0; h < g.hl.size(); ++h, ++hi_) {
        const int L = g.hl[h];
        hapOff[hi_] = (int)hoff; hapLen[hi_] = L;
        uint8_t* dst = hIn + o_hap + hoff;
        for (int k = 0; k < L; ++k) dst[k] = g_dev.conv[(uint8_t)g.hb[h][k]];
        hoff += L;
      }
    }
    pairStart[G] = ps;
    auto t1 = std::chrono::steady_clock::now();

    Batch b;
    b.nGroups = G;
    b.pairStart = (const int*)(dIn + o_pairStart); b.readStart = (const int*)(dIn + o_readStart);
    b.hapStart = (const int*)(dIn + o_hapStart);   b.nHaps = (const int*)(dIn + o_nHaps);
    b.readOff = (const int*)(dIn + o_readOff); b.readLen = (const int*)(dIn + o_readLen);
    b.hapOff = (const int*)(dIn + o_hapOff);   b.hapLen = (const int*)(dIn + o_hapLen);
    b.rs = dIn + o_rs; b.q = dIn + o_q; b.ins = dIn + o_i; b.del = dIn + o_d; b.gcp = dIn + o_c; b.hap = dIn + o_hap;

    const int bufStride = maxH + 1;
    // Each warp needs 3 rows of (maxH + 1) values in shared memory. Long haplotypes get fewer warps per
    // block (same results: every warp computes one pair on its own); beyond one warp's limit
    // (max_hap_len()) the caller must not send the batch.
    int wpb = WARPS_PER_BLK;
    while (wpb > 1 && (size_t)wpb * 3 * bufStride * sizeof(double) > (size_t)g_dev.maxSmem) wpb /= 2;
    if ((size_t)wpb * 3 * bufStride * sizeof(double) > (size_t)g_dev.maxSmem) {
      fprintf(stderr, "[cuda-pairhmm] haplotype of %d bases exceeds the shared-memory limit (max %d)\n", maxH, g_dev.max_hap_len());
      return false;
    }
    const int blocks = (nPairsB + wpb - 1) / wpb;
    const size_t smemF = (size_t)wpb * 3 * bufStride * sizeof(float);
    const size_t smemD = (size_t)wpb * 3 * bufStride * sizeof(double);

    CUDA_CHECK(cudaEventRecord(evStart, stream));
    CUDA_CHECK(cudaMemcpyAsync(dIn, hIn, total, cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaEventRecord(evH2D, stream));
    if (g_dev.useDouble) {
      pairhmm_double_kernel<<<blocks, wpb * WARP, smemD, stream>>>(b, nullptr, nPairsB, nullptr, dOut, dPath, bufStride);
    } else {
      CUDA_CHECK(cudaMemsetAsync(dNRedo, 0, sizeof(int), stream));
      pairhmm_float_kernel<<<blocks, wpb * WARP, smemF, stream>>>(b, nPairsB, dOut, dPath, dRedo, dNRedo, bufStride);
      // Pass 2 is always launched; warps beyond the device-side redo count exit immediately.
      pairhmm_double_kernel<<<blocks, wpb * WARP, smemD, stream>>>(b, dRedo, nPairsB, dNRedo, dOut, dPath, bufStride);
      CUDA_CHECK(cudaMemcpyAsync(hNRedo, dNRedo, sizeof(int), cudaMemcpyDeviceToHost, stream));
    }
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaEventRecord(evKernel, stream));
    CUDA_CHECK(cudaMemcpyAsync(hOut, dOut, nPairsB * sizeof(double), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaMemcpyAsync(hPath, dPath, nPairsB, cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaEventRecord(evEnd, stream));
    CUDA_CHECK(cudaStreamSynchronize(stream));
    // log10 on the host with glibc, exactly as GKL: (double)(log10f(f) - LOG10_INIT_F) or log10(d) - LOG10_INIT_D
    for (int p = 0; p < nPairsB; ++p) {
      if (hPath[p] == 0) hOut[p] = (double)(log10f((float)hOut[p]) - g_dev.log10InitF);
      else if (hPath[p] == 1) hOut[p] = log10(hOut[p]) - g_dev.log10InitD;
    }
    {
      size_t o = 0;
      for (const auto& g : groups) {
        const size_t n = g.rl.size() * g.hl.size();
        memcpy(g.out, hOut + o, n * sizeof(double));
        o += n;
      }
    }
    auto t2 = std::chrono::steady_clock::now();

    float a = 0, k = 0, d = 0;
    cudaEventElapsedTime(&a, evStart, evH2D); cudaEventElapsedTime(&k, evH2D, evKernel); cudaEventElapsedTime(&d, evKernel, evEnd);
    msH2D += a; msKernel += k; msD2H += d;
    nsPack += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
    nsGpu += std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
    nBatches += G; nLaunches++; nPairs += nPairsB; nCells += cells; nDouble += g_dev.useDouble ? nPairsB : *hNRedo;
    return true;
  }

};

// The engines (VARIANT_RS_GPU_STREAMS, default 8; 1: one at a time, as before): a caller takes a free one for the
// whole call. Each pair is computed on its own, so which engine and stream run it changes nothing in the results.
// HaplotypeCaller, 30x whole genome, L40S, 64 threads, 8 contig workers: 1 engine 190 s, 2 168, 4 168, 8 160 (the
// callers had waited 2,543 s in all for the one engine); with 16 workers, 4 engines 124 s, 8 121.
static std::mutex g_pool_mu;
static std::condition_variable g_pool_cv;
static std::vector<Engine*> g_pool, g_free;
static std::atomic<long long> g_nsTotal{0}, g_nsWait{0};

static int pool_size() {
  const char* e = getenv("VARIANT_RS_GPU_STREAMS");
  int n = e ? atoi(e) : 8;
  return n < 1 ? 1 : n > 64 ? 64 : n;
}
static Engine* acquire() {
  std::unique_lock<std::mutex> lk(g_pool_mu);
  g_pool_cv.wait(lk, [] { return !g_free.empty(); });
  Engine* e = g_free.back();
  g_free.pop_back();
  return e;
}
static void release(Engine* e) {
  { std::lock_guard<std::mutex> lk(g_pool_mu); g_free.push_back(e); }
  g_pool_cv.notify_one();
}

static void report_all() {
  long long nsPack = 0, nsGpu = 0, nBatches = 0, nLaunches = 0, nPairs = 0, nCells = 0, nDouble = 0;
  double msH2D = 0, msKernel = 0, msD2H = 0;
  for (Engine* e : g_pool) {
    nsPack += e->nsPack; nsGpu += e->nsGpu; nBatches += e->nBatches; nLaunches += e->nLaunches; nPairs += e->nPairs;
    nCells += e->nCells; nDouble += e->nDouble; msH2D += e->msH2D; msKernel += e->msKernel; msD2H += e->msD2H;
  }
  if (nBatches == 0) return;
  fprintf(stderr,
          "[cuda-pairhmm] batches=%lld launches=%lld pairs=%lld cells=%.3e doubleFallback=%lld\n"
          "[cuda-pairhmm] %d engine(s), summed over them: pack=%.3fs  GPU wall=%.3fs (H2D=%.3fs kernel=%.3fs D2H=%.3fs)  total=%.3fs; waiting for an engine %.3fs\n"
          "[cuda-pairhmm] kernel throughput=%.1f GCUPS (summed kernel time), avg pairs/launch=%.0f\n",
          nBatches, nLaunches, nPairs, (double)nCells, nDouble,
          (int)g_pool.size(), nsPack * 1e-9, nsGpu * 1e-9, msH2D * 1e-3, msKernel * 1e-3, msD2H * 1e-3, g_nsTotal * 1e-9, g_nsWait * 1e-9,
          nCells / (msKernel * 1e-3) * 1e-9, (double)nPairs / nLaunches);
}


// ---------------------------------------------------------------------------
// Plain C ABI (used by the Rust HaplotypeCaller port). Several regions ("groups")
// are computed in ONE launch. Reads and haplotypes are passed as arrays of pointers,
// flattened across groups in group order; `out` receives, for each group in order,
// a row-major n_reads x n_haps matrix of log10 likelihoods.
// ---------------------------------------------------------------------------
extern "C" {

int cuda_pairhmm_init(int use_double) {
  if (g_dev.ready) return 1;
  if (!g_dev.init(use_double != 0)) return 0;
  for (int i = 0, n = pool_size(); i < n; ++i) {
    Engine* e = new Engine();
    if (!e->create()) { g_dev.ready = false; return 0; }
    g_pool.push_back(e);
    g_free.push_back(e);
  }
  return 1;
}

int cuda_pairhmm_compute(int n_groups, const int* group_n_reads, const int* group_n_haps,
                         const uint8_t* const* read_bases, const uint8_t* const* read_quals,
                         const uint8_t* const* read_ins, const uint8_t* const* read_del,
                         const uint8_t* const* read_gcp, const int* read_lens,
                         const uint8_t* const* hap_bases, const int* hap_lens, double* out) {
  if (!g_dev.ready) return 0;
  std::vector<Engine::HostGroup> groups(n_groups);
  size_t ri = 0, hi = 0, oo = 0;
  for (int g = 0; g < n_groups; ++g) {
    auto& G = groups[g];
    const int nr = group_n_reads[g], nh = group_n_haps[g];
    for (int r = 0; r < nr; ++r, ++ri) {
      G.rb.push_back((const int8_t*)read_bases[ri]); G.rq.push_back((const int8_t*)read_quals[ri]);
      G.ri.push_back((const int8_t*)read_ins[ri]);   G.rd.push_back((const int8_t*)read_del[ri]);
      G.rc.push_back((const int8_t*)read_gcp[ri]);   G.rl.push_back(read_lens[ri]);
    }
    for (int h = 0; h < nh; ++h, ++hi) { G.hb.push_back((const int8_t*)hap_bases[hi]); G.hl.push_back(hap_lens[hi]); }
    G.out = out + oo;
    oo += (size_t)nr * nh;
  }
  auto tw = std::chrono::steady_clock::now();
  Engine* e = acquire();
  auto t0 = std::chrono::steady_clock::now();
  const bool ok = e->compute(groups);
  release(e);
  g_nsWait += std::chrono::duration_cast<std::chrono::nanoseconds>(t0 - tw).count();
  g_nsTotal += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - t0).count();
  return ok ? 1 : 0;
}

void cuda_pairhmm_report(void) { report_all(); }

// Longest haplotype the GPU path accepts (0 if not initialized); longer ones need the CPU PairHMM.
int cuda_pairhmm_max_hap_len(void) { return g_dev.ready ? g_dev.max_hap_len() : 0; }

// Debug: float pass for ONE pair, writing every DP cell (M,X,Y interleaved, (R+1) x (H+1)) to `cells`
// and the raw float sum to `*raw`.  Reuses the production code path (pairhmm_warp<float, true>).
}  // extern "C"
__global__ void pairhmm_debug_kernel(Batch b, float* cells, float* raw) {
  extern __shared__ unsigned char smem[];
  const int bufStride = b.hapLen[0] + 1;
  float* buf = reinterpret_cast<float*>(smem);
  const float res = pairhmm_warp<float, true>(b, 0, 0, buf, buf + bufStride, buf + 2 * bufStride, threadIdx.x, cells);
  if (threadIdx.x == 0) *raw = res;
}
extern "C" {
int cuda_pairhmm_debug_cells(const uint8_t* rb, const uint8_t* rq, const uint8_t* ri, const uint8_t* rd, const uint8_t* rc,
                             int R, const uint8_t* hb, int H, float* cells, float* raw) {
  if (!g_dev.ready) return 0;
  std::vector<int> ints = {0, 1, 0, 0, 1, 0, R, 0, H};  // pairStart[2] readStart hapStart nHaps readOff readLen hapOff hapLen
  std::vector<uint8_t> bytes(5 * R + H);
  for (int k = 0; k < R; ++k) { bytes[k] = g_dev.conv[rb[k]]; bytes[R + k] = rq[k]; bytes[2 * R + k] = ri[k]; bytes[3 * R + k] = rd[k]; bytes[4 * R + k] = rc[k]; }
  for (int k = 0; k < H; ++k) bytes[5 * R + k] = g_dev.conv[hb[k]];
  int* dI; uint8_t* dB; float* dC; float* dR;
  const size_t nc = 3 * (size_t)(R + 1) * (H + 1);
  cudaMalloc(&dI, ints.size() * sizeof(int)); cudaMalloc(&dB, bytes.size()); cudaMalloc(&dC, nc * sizeof(float)); cudaMalloc(&dR, sizeof(float));
  cudaMemcpy(dI, ints.data(), ints.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(dB, bytes.data(), bytes.size(), cudaMemcpyHostToDevice);
  cudaMemset(dC, 0, nc * sizeof(float));
  Batch b; b.nGroups = 1; b.pairStart = dI; b.readStart = dI + 2; b.hapStart = dI + 3; b.nHaps = dI + 4;
  b.readOff = dI + 5; b.readLen = dI + 6; b.hapOff = dI + 7; b.hapLen = dI + 8;
  b.rs = dB; b.q = dB + R; b.ins = dB + 2 * R; b.del = dB + 3 * R; b.gcp = dB + 4 * R; b.hap = dB + 5 * R;
  pairhmm_debug_kernel<<<1, WARP, 3 * (H + 1) * sizeof(float)>>>(b, dC, dR);
  cudaMemcpy(cells, dC, nc * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(raw, dR, sizeof(float), cudaMemcpyDeviceToHost);
  cudaFree(dI); cudaFree(dB); cudaFree(dC); cudaFree(dR);
  return cudaGetLastError() == cudaSuccess ? 1 : 0;
}

}  // extern "C"
