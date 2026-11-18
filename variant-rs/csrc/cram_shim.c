// Minimal htslib wrapper for variant-rs CRAM input: open + index, multi-interval query, and each record
// serialized exactly as htslib's bam_write1 would put it in a BAM file (the bytes after block_size).
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/cram.h>

typedef struct {
    htsFile *fp;
    sam_hdr_t *hdr;
    hts_idx_t *idx;
    hts_itr_t *itr;
    bam1_t *b;
    uint8_t *buf;
    size_t cap;
} crh;

crh *crh_open(const char *path, const char *fasta, crh *share) {
    crh *h = calloc(1, sizeof(crh));
    if (!h) return NULL;
    h->fp = hts_open(path, "r");
    if (!h->fp || h->fp->format.format != cram) goto fail;
    if (hts_set_fai_filename(h->fp, fasta) != 0) goto fail;
    // No CRAM_OPT_SHARED_REF: a shared cache keeps only the last-used contig, so handles working on
    // different contigs evict each other and reload whole chromosomes per slice. Each handle keeps
    // its current contig (<= ~250 MB) instead.
    (void)share;
    h->hdr = sam_hdr_read(h->fp);
    if (!h->hdr) goto fail;
    h->idx = sam_index_load(h->fp, path);
    if (!h->idx) goto fail;
    h->b = bam_init1();
    return h;
fail:
    if (h->fp) hts_close(h->fp);
    if (h->hdr) sam_hdr_destroy(h->hdr);
    free(h);
    return NULL;
}

void crh_close(crh *h) {
    if (!h) return;
    if (h->itr) hts_itr_destroy(h->itr);
    if (h->b) bam_destroy1(h->b);
    if (h->idx) hts_idx_destroy(h->idx);
    if (h->hdr) sam_hdr_destroy(h->hdr);
    if (h->fp) hts_close(h->fp);
    free(h->buf);
    free(h);
}

const char *crh_header_text(crh *h, size_t *len) {
    *len = sam_hdr_length(h->hdr);
    return sam_hdr_str(h->hdr);
}
int crh_nref(crh *h) { return sam_hdr_nref(h->hdr); }
const char *crh_ref_name(crh *h, int tid) { return sam_hdr_tid2name(h->hdr, tid); }
int64_t crh_ref_len(crh *h, int tid) { return sam_hdr_tid2len(h->hdr, tid); }

// Records of `tid` overlapping any of the n intervals iv[2k]..iv[2k+1] (1-based, inclusive), each once,
// in file order. Returns 0, or -1 on error.
int crh_query(crh *h, int tid, const int64_t *iv, int n) {
    if (h->itr) { hts_itr_destroy(h->itr); h->itr = NULL; }
    if (n <= 0) return -1;
    hts_reglist_t *rl = calloc(1, sizeof(hts_reglist_t));
    hts_pair_pos_t *p = malloc(sizeof(hts_pair_pos_t) * n);
    if (!rl || !p) return -1;
    for (int k = 0; k < n; k++) { p[k].beg = iv[2 * k] - 1; p[k].end = iv[2 * k + 1]; }
    rl->reg = sam_hdr_tid2name(h->hdr, tid);
    rl->intervals = p;
    rl->tid = tid;
    rl->count = n;
    rl->min_beg = p[0].beg;
    rl->max_end = p[n - 1].end;
    h->itr = sam_itr_regions(h->idx, h->hdr, rl, 1); // takes ownership of rl
    return h->itr ? 0 : -1;
}

// The unplaced unmapped records (RNAME '*', after all placed ones), in file order. Returns 0, or -1 on error.
int crh_query_unplaced(crh *h) {
    if (h->itr) { hts_itr_destroy(h->itr); h->itr = NULL; }
    h->itr = sam_itr_queryi(h->idx, HTS_IDX_NOCOOR, 0, 0);
    return h->itr ? 0 : -1;
}

// Next record as BAM bytes (no block_size). Returns its length (*out valid until the next call),
// -1 at the end, < -1 on error.
int64_t crh_next(crh *h, const uint8_t **out) {
    int r = sam_itr_next(h->fp, h->itr, h->b);
    if (r < 0) return r < -1 ? r : -1;
    const bam1_t *b = h->b;
    const bam1_core_t *c = &b->core;
    if (c->n_cigar > 0xffff || c->pos > INT32_MAX || c->mpos > INT32_MAX) return -3;
    uint32_t lq = c->l_qname - c->l_extranul;
    size_t len = 32 + lq + (b->l_data - c->l_qname);
    if (len > h->cap) {
        h->cap = len * 2;
        h->buf = realloc(h->buf, h->cap);
        if (!h->buf) return -4;
    }
    uint8_t *x = h->buf;
    int32_t v[8];
    v[0] = c->tid;
    v[1] = (int32_t)c->pos;
    v[2] = (int32_t)((uint32_t)c->bin << 16 | (uint32_t)c->qual << 8 | lq);
    v[3] = (int32_t)((uint32_t)c->flag << 16 | c->n_cigar);
    v[4] = c->l_qseq;
    v[5] = c->mtid;
    v[6] = (int32_t)c->mpos;
    v[7] = (int32_t)c->isize;
    memcpy(x, v, 32); // little-endian host
    memcpy(x + 32, b->data, lq);
    memcpy(x + 32 + lq, b->data + c->l_qname, b->l_data - c->l_qname);
    *out = h->buf;
    return (int64_t)len;
}
