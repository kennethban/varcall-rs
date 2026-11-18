// align-rs <-> htslib: SAM line -> BAM record (sam_parse1, as samtools does), and the sorted BAM writer with its
// on-the-fly index (sam_idx_init, as `samtools sort --write-index` does).
//
// A record crosses into Rust as an opaque blob: bam1_core_t followed by bam1_t.data (l_data bytes).
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <stdio.h>

#define FAIL(what) do { fprintf(stderr, "[align-rs] htslib: %s failed (%s)\n", what, path); return NULL; } while (0)

typedef struct {
    samFile *fp;
    sam_hdr_t *h;
    bam1_t *b;
    char *idx;
} ah_writer;

// A header for parsing: the reference name hash is built here, so that sam_parse1 only reads it afterwards
// (thread-safe)
sam_hdr_t *ah_hdr(const char *text, size_t len) {
    sam_hdr_t *h = sam_hdr_parse(len, text);
    if (!h) return NULL;
    if (sam_hdr_nref(h) > 0 && sam_hdr_name2tid(h, sam_hdr_tid2name(h, 0)) != 0) return NULL;
    return h;
}

int ah_core_size(void) { return (int)sizeof(bam1_core_t); }

static __thread bam1_t *tb;
static __thread kstring_t ks;

// Parses one SAM line (no newline). Returns the blob size (written to `out` when it fits in `cap`), or -1 on a
// parse error. tid / pos / flag: for sorting.
int64_t ah_parse(sam_hdr_t *h, const char *line, size_t len, uint8_t *out, size_t cap, int32_t *tid, int64_t *pos, uint16_t *flag) {
    if (!tb) tb = bam_init1();
    ks.l = 0;
    if (kputsn(line, len, &ks) < 0) return -1;
    if (sam_parse1(&ks, h, tb) < 0) return -1;
    size_t n = sizeof(bam1_core_t) + (size_t)tb->l_data;
    *tid = tb->core.tid;
    *pos = tb->core.pos;
    *flag = tb->core.flag;
    if (n <= cap) {
        memcpy(out, &tb->core, sizeof(bam1_core_t));
        memcpy(out + sizeof(bam1_core_t), tb->data, tb->l_data);
    }
    return (int64_t)n;
}

// mode: "wb1" etc.; the header text gets @HD SO:coordinate as samtools sort sets it
ah_writer *ah_open(const char *path, const char *mode, const char *text, size_t len, int threads) {
    ah_writer *w = calloc(1, sizeof(ah_writer));
    w->h = sam_hdr_parse(len, text);
    if (!w->h) FAIL("sam_hdr_parse");
    if (sam_hdr_count_lines(w->h, "HD") == 0) {
        if (sam_hdr_add_line(w->h, "HD", "VN", SAM_FORMAT_VERSION, "SO", "coordinate", NULL) < 0) FAIL("sam_hdr_add_line");
    } else if (sam_hdr_change_HD(w->h, "SO", "coordinate") < 0) FAIL("sam_hdr_change_HD");
    w->fp = sam_open(path, mode);
    if (!w->fp) FAIL("sam_open");
    if (threads > 1 && hts_set_threads(w->fp, threads) < 0) FAIL("hts_set_threads");
    if (sam_hdr_write(w->fp, w->h) < 0) FAIL("sam_hdr_write");
    size_t pl = strlen(path);
    w->idx = malloc(pl + 5);
    memcpy(w->idx, path, pl);
    memcpy(w->idx + pl, ".bai", 5);
    if (sam_idx_init(w->fp, w->h, 0, w->idx) < 0) FAIL("sam_idx_init");
    w->b = bam_init1();
    return w;
}

int ah_write(ah_writer *w, const uint8_t *blob, size_t n) {
    bam1_t *b = w->b;
    size_t l = n - sizeof(bam1_core_t);
    if (l > b->m_data) {
        uint8_t *d = realloc(b->data, l);
        if (!d) return -1;
        b->data = d;
        b->m_data = l;
    }
    memcpy(&b->core, blob, sizeof(bam1_core_t));
    memcpy(b->data, blob + sizeof(bam1_core_t), l);
    b->l_data = (int)l;
    return sam_write1(w->fp, w->h, b) < 0 ? -1 : 0;
}

int ah_close(ah_writer *w) {
    int r = 0;
    if (sam_idx_save(w->fp) < 0) r = -1;
    if (sam_close(w->fp) < 0) r = -1;
    bam_destroy1(w->b);
    sam_hdr_destroy(w->h);
    free(w->idx);
    free(w);
    return r;
}
