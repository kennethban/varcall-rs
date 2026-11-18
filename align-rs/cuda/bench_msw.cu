// Mate-SW kernel benchmark: runs aln_gpu_matesw on a chunk dumped by align-rs (ALIGN_RS_DUMP_MSW=FILE) and checks every
// result against the dumped one (the kernel that produced it gives bwa-mem2's output).
//   nvcc ... -o bench_msw bench_msw.cu && ./bench_msw REF.fa DUMP [REPS]
#include "aln_ext.cu"
#include <vector>
#include <string>
#include <chrono>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

static std::vector<char> slurp(const char *p) {
    FILE *f = fopen(p, "rb"); if (!f) { perror(p); exit(1); }
    fseek(f, 0, SEEK_END); long n = ftell(f); fseek(f, 0, SEEK_SET);
    std::vector<char> v(n); if (fread(v.data(), 1, n, f) != (size_t)n) exit(1); fclose(f); return v;
}
int main(int argc, char **argv) {
    std::string ref = argv[1];
    int reps = argc > 3 ? atoi(argv[3]) : 3;
    int64_t l_pac = 0;
    { FILE *a = fopen((ref + ".ann").c_str(), "r"); if (fscanf(a, "%ld", &l_pac) != 1) exit(1); fclose(a); }
    int fd = open((ref + ".pac").c_str(), O_RDONLY); struct stat st; fstat(fd, &st);
    const uint8_t *pac = (const uint8_t *)mmap(0, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    int8_t mat[25];
    for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) mat[i * 5 + j] = (i < 4 && j < 4) ? (i == j ? 1 : -4) : -1;
    if (aln_gpu_init(pac, (st.st_size), l_pac, mat, 6, 1, 6, 1, 100)) { fprintf(stderr, "init failed\n"); return 1; }
    std::vector<char> d = slurp(argv[2]);
    const char *p = d.data();
    auto rd = [&](size_t n) { const char *q = p; p += n; return q; };
    uint64_t n_reads = *(const uint64_t *)rd(8), total = *(const uint64_t *)rd(8);
    const int64_t *off = (const int64_t *)rd(8 * (n_reads + 1));
    const uint8_t *reads = (const uint8_t *)rd(total);
    uint64_t n = *(const uint64_t *)rd(8), scratch = *(const uint64_t *)rd(8);
    const MJob *jobs = (const MJob *)rd(n * sizeof(MJob));
    const MOut *want = (const MOut *)rd(n * sizeof(MOut));
    if (aln_gpu_set_reads(2, 1, reads, total, off, (int)n_reads)) return 1;
    if (aln_gpu_reads_cached(2, 1) != 1) return 1;
    std::vector<MOut> out(n);
    double cells = 0; for (uint64_t i = 0; i < n; ++i) cells += (double)jobs[i].qlen * jobs[i].tlen;
    for (int r = 0; r < reps; ++r) {
        auto t0 = std::chrono::steady_clock::now();
        if (aln_gpu_matesw(jobs, (int)n, (int64_t)scratch, out.data())) { fprintf(stderr, "matesw failed\n"); return 1; }
        double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        printf("%lu jobs, %.2e cells: %.1f ms, %.0f GCUPS (incl. transfers)\n", n, cells, s * 1e3, cells / s * 1e-9);
    }
    uint64_t bad = 0;
    for (uint64_t i = 0; i < n; ++i) {
        const MOut &a = out[i], &b = want[i];
        if (a.ok != b.ok || (a.ok && (a.score != b.score || a.te != b.te || a.qe != b.qe || a.score2 != b.score2 || a.te2 != b.te2 || a.tb != b.tb || a.qb != b.qb))) {
            if (bad++ < 5) fprintf(stderr, "job %lu: got %d %d %d %d %d %d %d, want %d %d %d %d %d %d %d\n", i, a.score, a.te, a.qe, a.score2, a.te2, a.tb, a.qb, b.score, b.te, b.qe, b.score2, b.te2, b.tb, b.qb);
        }
    }
    printf("%s: %lu of %lu results differ\n", bad ? "DIFFERENT" : "identical", bad, n);
    aln_gpu_report3();
    return bad != 0;
}
