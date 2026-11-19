// align-rs <-> htslib: SAM line -> BAM record (sam_parse1, as samtools does), serialized as bam_write1 writes it;
// the sorted BAM's header (sam_hdr_write) and raw BGZF blocks. The blocks are made in Rust (bgzf_compress, in
// parallel) and the .bai with htslib's index builder (hts_idx_push), as `samtools sort --write-index` does.
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <stdio.h>
#include <limits.h>

#if !defined(__BYTE_ORDER__) || __BYTE_ORDER__ != __ORDER_LITTLE_ENDIAN__
#error "align-rs writes BAM records in host byte order: little-endian hosts only"
#endif

#define FAIL(what) do { fprintf(stderr, "[align-rs] htslib: %s failed (%s)\n", what, path); return NULL; } while (0)

typedef struct {
    samFile *fp;
    sam_hdr_t *h;
} ah_writer;

// A header for parsing: the reference name hash is built here, so that sam_parse1 only reads it afterwards
// (thread-safe)
sam_hdr_t *ah_hdr(const char *text, size_t len) {
    sam_hdr_t *h = sam_hdr_parse(len, text);
    if (!h) return NULL;
    if (sam_hdr_nref(h) > 0 && sam_hdr_name2tid(h, sam_hdr_tid2name(h, 0)) != 0) return NULL;
    return h;
}

static __thread bam1_t *tb;
static __thread kstring_t ks;

// The record as BAM stores it, byte for byte what bam_write1 (htslib 1.25, sam.c) writes: block_size, the 32-byte
// fixed part, the read name without its extra NULs, the rest of the data; more than 65535 CIGAR operations become
// <l_qseq>S<ref length>N with the real CIGAR in a CG:B,I tag at the end. Returns its size (written when it fits in
// `cap`), or -1 where bam_write1 refuses the record.
static int64_t ondisk(const bam1_t *b, uint8_t *out, size_t cap) {
    const bam1_core_t *c = &b->core;
    uint32_t block_len = b->l_data - c->l_extranul + 32;
    if (c->l_qname - c->l_extranul > BAM_MAX_QNAME_LEN + 1) return -1;
    if (c->n_cigar > 0xffff) block_len += 16;
    if (c->pos > INT_MAX || c->mpos > INT_MAX || c->isize < INT_MIN || c->isize > INT_MAX) return -1;
    size_t n = 4 + (size_t)block_len;
    if (n > cap) return (int64_t)n;
    uint32_t x[8];
    x[0] = c->tid;
    x[1] = c->pos;
    x[2] = (uint32_t)c->bin << 16 | c->qual << 8 | (c->l_qname - c->l_extranul);
    if (c->n_cigar > 0xffff) x[3] = (uint32_t)c->flag << 16 | 2;
    else x[3] = (uint32_t)c->flag << 16 | (c->n_cigar & 0xffff);
    x[4] = c->l_qseq;
    x[5] = c->mtid;
    x[6] = c->mpos;
    x[7] = c->isize;
    uint8_t *p = out;
    memcpy(p, &block_len, 4); p += 4;
    memcpy(p, x, 32); p += 32;
    memcpy(p, b->data, c->l_qname - c->l_extranul); p += c->l_qname - c->l_extranul;
    if (c->n_cigar <= 0xffff) {
        memcpy(p, b->data + c->l_qname, b->l_data - c->l_qname);
    } else {
        hts_pos_t cigreflen = bam_cigar2rlen(c->n_cigar, bam_get_cigar(b));
        if (cigreflen >= (1 << 28)) return -1;
        uint32_t cigar_st = (uint8_t *)bam_get_cigar(b) - b->data, cigar_en = cigar_st + c->n_cigar * 4;
        uint32_t cigar[2] = {(uint32_t)c->l_qseq << 4 | BAM_CSOFT_CLIP, (uint32_t)cigreflen << 4 | BAM_CREF_SKIP};
        memcpy(p, cigar, 8); p += 8;
        memcpy(p, b->data + cigar_en, b->l_data - cigar_en); p += b->l_data - cigar_en;
        memcpy(p, "CGBI", 4); p += 4;
        uint32_t nc = c->n_cigar;
        memcpy(p, &nc, 4); p += 4;
        memcpy(p, b->data + cigar_st, c->n_cigar * 4);
    }
    return (int64_t)n;
}

// Parses one SAM line (no newline). Returns the BAM record's size (written to `out` when it fits in `cap`), -1 on a
// parse error, -2 if BAM cannot store it. tid / pos / flag: for sorting; end: bam_endpos, for the index.
int64_t ah_parse(sam_hdr_t *h, const char *line, size_t len, uint8_t *out, size_t cap, int32_t *tid, int64_t *pos, uint16_t *flag, int64_t *end) {
    if (!tb) tb = bam_init1();
    ks.l = 0;
    if (kputsn(line, len, &ks) < 0) return -1;
    if (sam_parse1(&ks, h, tb) < 0) return -1;
    *tid = tb->core.tid;
    *pos = tb->core.pos;
    *flag = tb->core.flag;
    *end = bam_endpos(tb);
    int64_t n = ondisk(tb, out, cap);
    return n < 0 ? -2 : n;
}

// mode: "wb1" etc.; the header text gets @HD SO:coordinate as samtools sort sets it. The header is written and
// flushed (bam_hdr_write ends its block): the records' blocks follow with ah_raw.
ah_writer *ah_open(const char *path, const char *mode, const char *text, size_t len) {
    ah_writer *w = calloc(1, sizeof(ah_writer));
    w->h = sam_hdr_parse(len, text);
    if (!w->h) FAIL("sam_hdr_parse");
    if (sam_hdr_count_lines(w->h, "HD") == 0) {
        if (sam_hdr_add_line(w->h, "HD", "VN", SAM_FORMAT_VERSION, "SO", "coordinate", NULL) < 0) FAIL("sam_hdr_add_line");
    } else if (sam_hdr_change_HD(w->h, "SO", "coordinate") < 0) FAIL("sam_hdr_change_HD");
    w->fp = sam_open(path, mode);
    if (!w->fp) FAIL("sam_open");
    if (sam_hdr_write(w->fp, w->h) < 0) FAIL("sam_hdr_write");
    return w;
}

// the virtual offset after the header (the first record's)
uint64_t ah_tell(ah_writer *w) { return (uint64_t)bgzf_tell(w->fp->fp.bgzf); }

// whole BGZF blocks, after the header
int ah_raw(ah_writer *w, const uint8_t *data, size_t n) {
    return bgzf_raw_write(w->fp->fp.bgzf, data, n) == (ssize_t)n ? 0 : -1;
}

// writes the EOF block and closes
int ah_close(ah_writer *w) {
    int r = sam_close(w->fp) < 0 ? -1 : 0;
    sam_hdr_destroy(w->h);
    free(w);
    return r;
}
