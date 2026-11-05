/* magibu-meth: per-CpG methylation and entropy from a modified-base BAM (MM/ML tags).
   Written for Magibu with htslib; the MM/ML tags follow the public SAM tags specification (hts-specs).
   Replaces  modkit extract full --cpg --mapped-only  +  ModkitResorter.sh  (reshape, sort,
   ParseModkit.pl of PoreMeth2).

   Output, one row per CpG (0-based C of the CpG, strands combined), same columns as
   *_sorted.entropy.file.tsv:
     chr  pos  entropy  reads_with_3_CpG_pattern  beta  coverage
   Rules of ParseModkit.pl: per read call p5mC, p5hmC = (ML + 0.5) / 256 (0 for implicit calls);
   methylated if p5mC > 0.85, unmethylated if p5mC + p5hmC < 0.3, otherwise ambiguous.
   beta = methylated / (methylated + unmethylated). The entropy of a CpG uses the reads that
   carry it and the next two CpGs (8 possible patterns), normalised by 3. A CpG is written when
   both the reads with the full pattern and the confident calls are more than 2.
   Primary alignments only, chromosomes 1-22, X, Y (as ModkitResorter.sh).
   Without -r: whole genome (PoreMeth2); with -r: only those regions (Magibu, much faster).
   Entropy and beta are written with 4 decimals; --full-precision writes them as Perl does.
   --no-entropy: beta only (same rules), every CpG with at least --min-cov confident calls
   (default 3, as ParseModkit.pl), entropy and pattern reads written as NA.
   --sum-5hmc: 5mC + 5hmC counted as methylation (methylated if p5mC + p5hmC > 0.85).
   --split-5hmc: 5mC and 5hmC kept apart. Each call is 5mC if p5mC > 0.85, 5hmC if p5hmC > 0.85,
   unmodified C if p5mC + p5hmC < 0.3, otherwise ambiguous; three columns are added after the six
   above (which do not change):
     beta_5mC = n5mC / n  beta_5hmC = n5hmC / n  n = n5mC + n5hmC + nC
   beta_5hmC is NA when no read at the CpG carries 5hmC calls (BAM basecalled with a 5mC-only model).
   CpG sites: from the reference with --ref (as modkit --cpg --reference), otherwise from the
   read (C followed by G), dropping sites seen as CpG in less than --min-cpg-frac of the reads.

   usage: magibu-meth -o out.entropy.file.tsv [-r regions.bed] [--ref ref.fasta] [-t threads]
                      [--min-cpg-frac F] [--full-precision] [--no-entropy] [--min-cov N]
                      [--sum-5hmc] [--split-5hmc] in.bam
   (--entropy is accepted and ignored: it is the only mode) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include <math.h>
#include <ctype.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>

#define MAX_SPAN 200000
#define MAX_REG_PER_CHUNK 256

typedef struct { int tid; hts_pos_t beg, end; } region_t;            /* 0-based, half open */
typedef struct { int tid; hts_pos_t beg, end; int r0, r1; } chunk_t;  /* regions r0..r1-1 */

static region_t *regs; static int nregs;
static chunk_t *chunks; static int nchunks;
static const char *bam_path;
static double min_frac = 0.2;      /* min fraction of covering reads that show the CpG */
static int done_chunks = 0;
static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
static long n_reads_used = 0, n_reads_bad = 0;

static void die(const char *m, const char *a) { fprintf(stderr, "ERROR: %s %s\n", m, a ? a : ""); exit(1); }

/* ---- regions ---- */
static int cmp_reg(const void *a, const void *b) {
    const region_t *x = a, *y = b;
    if (x->tid != y->tid) return x->tid < y->tid ? -1 : 1;
    return x->beg < y->beg ? -1 : (x->beg > y->beg);
}

static int name2tid(sam_hdr_t *h, const char *name) {
    int tid = sam_hdr_name2tid(h, name);
    if (tid >= 0) return tid;
    char alt[300];
    if (strncmp(name, "chr", 3) == 0) snprintf(alt, sizeof alt, "%s", name + 3);
    else snprintf(alt, sizeof alt, "chr%s", name);
    tid = sam_hdr_name2tid(h, alt);
    if (tid < 0 && (strcmp(name, "chrM") == 0 || strcmp(name, "MT") == 0))
        tid = sam_hdr_name2tid(h, strcmp(name, "chrM") == 0 ? "MT" : "chrM");
    return tid;
}

static void read_regions(const char *path, sam_hdr_t *h) {
    FILE *f = fopen(path, "r");
    if (!f) die("cannot open regions file", path);
    char line[4096], chr[256]; long long b, e; int cap = 65536, skipped = 0;
    regs = malloc(cap * sizeof *regs);
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#' || strncmp(line, "track", 5) == 0 || strncmp(line, "browser", 7) == 0) continue;
        if (sscanf(line, "%255s %lld %lld", chr, &b, &e) != 3) continue;
        int tid = name2tid(h, chr);
        if (tid < 0) { skipped++; continue; }
        if (b < 0) b = 0;
        if (e <= b) continue;
        if (nregs == cap) { cap *= 2; regs = realloc(regs, cap * sizeof *regs); }
        regs[nregs].tid = tid; regs[nregs].beg = b; regs[nregs].end = e; nregs++;
    }
    fclose(f);
    if (skipped) fprintf(stderr, "regions on contigs absent from the BAM: %d (skipped)\n", skipped);
    if (!nregs) die("no usable regions in", path);
    qsort(regs, nregs, sizeof *regs, cmp_reg);
    int m = 0;                                       /* merge overlapping/adjacent */
    for (int i = 1; i < nregs; i++) {
        if (regs[i].tid == regs[m].tid && regs[i].beg <= regs[m].end) {
            if (regs[i].end > regs[m].end) regs[m].end = regs[i].end;
        } else regs[++m] = regs[i];
    }
    nregs = m + 1;
    int extra = 0;                                   /* split regions longer than a chunk */
    for (int i = 0; i < nregs; i++) extra += (int) ((regs[i].end - regs[i].beg - 1) / MAX_SPAN);
    if (extra) {
        region_t *nr = malloc((nregs + extra) * sizeof *nr); int k = 0;
        for (int i = 0; i < nregs; i++)
            for (hts_pos_t s = regs[i].beg; s < regs[i].end; s += MAX_SPAN) {
                nr[k].tid = regs[i].tid; nr[k].beg = s; nr[k].end = s + MAX_SPAN < regs[i].end ? s + MAX_SPAN : regs[i].end; k++;
            }
        free(regs); regs = nr; nregs = k;
    }
    chunks = malloc(nregs * sizeof *chunks);
    for (int i = 0; i < nregs; ) {
        chunk_t c = { regs[i].tid, regs[i].beg, regs[i].end, i, i + 1 };
        int j = i + 1;
        while (j < nregs && regs[j].tid == c.tid && regs[j].end - c.beg <= MAX_SPAN && j - i < MAX_REG_PER_CHUNK) {
            c.end = regs[j].end; j++;
        }
        c.r1 = j; chunks[nchunks++] = c; i = j;
    }
}

/* one base more on the right: a reverse strand read that starts on the G of the last CpG
   carries the call of that CpG */
static void reg_string(int tid, hts_pos_t beg, hts_pos_t end, sam_hdr_t *h, char *buf, size_t n) {
    snprintf(buf, n, "{%s}:%lld-%lld", sam_hdr_tid2name(h, tid), (long long) beg + 1, (long long) end + 1);
}

/* ---- calls ---- */
/* st: status for beta and entropy (1 methylated, 0 unmethylated, -1 ambiguous; --sum-5hmc rule).
   cls: split status (CLS_*), plus CLS_H_SEEN when the read carries 5hmC calls. */
enum { CLS_C = 0, CLS_M = 1, CLS_H = 2, CLS_AMB = 3, CLS_H_SEEN = 4 };
typedef struct { hts_pos_t pos; uint64_t uid; int8_t st, cls; } call_t;
typedef struct { hts_pos_t pos; int n, cap, cov, meth, nm, nh, nc, hseen; uint64_t *uid; int8_t *st; } wpos_t;
typedef struct { int tid; int c0, c1; } chrom_t;          /* chunks c0..c1-1 */

static const char *ref_path = NULL;
static chrom_t *echr; static int nechr;
static int next_chrom = 0;
static char **tmp_names;
static long ent_rows = 0;

static uint64_t name_uid(const char *s) {               /* FNV-1a 64 */
    uint64_t hsh = 1469598103934665603ULL;
    for (; *s; s++) { hsh ^= (unsigned char) *s; hsh *= 1099511628211ULL; }
    return hsh;
}

static int cmp_call(const void *a, const void *b) {
    const call_t *x = a, *y = b;
    if (x->pos != y->pos) return x->pos < y->pos ? -1 : 1;
    return x->uid < y->uid ? -1 : (x->uid > y->uid);
}

typedef struct { call_t *v; size_t n, cap; } callvec_t;
static inline void cv_push(callvec_t *c, hts_pos_t pos, uint64_t uid, int st, int cls) {
    if (c->n == c->cap) { c->cap = c->cap ? 2 * c->cap : 65536; c->v = realloc(c->v, c->cap * sizeof(call_t)); }
    c->v[c->n].pos = pos; c->v[c->n].uid = uid; c->v[c->n].st = (int8_t) st; c->v[c->n].cls = (int8_t) cls; c->n++;
}

static int ent_sum_h = 0;             /* --sum-5hmc: methylated if p5mC + p5hmC > 0.85 */
static int ent_split = 0;             /* --split-5hmc: add beta_5mC, beta_5hmC, n */
static inline int call_status(double ml, double hl) {
    if ((ent_sum_h ? ml + hl : ml) > 0.85) return 1;
    if (ml + hl < 0.3) return 0;
    return -1;
}
static inline int call_class(double ml, double hl, int has_h) {
    int c = ml > 0.85 ? CLS_M : hl > 0.85 ? CLS_H : ml + hl < 0.3 ? CLS_C : CLS_AMB;
    return c | (has_h ? CLS_H_SEEN : 0);
}

/* calls of one read inside [beg, end) (and in mask if given). ref = reference for [beg-1, end+1)
   upper case, or NULL (read-based CpG). depth: per position read depth (read-based mode only). */
static int entropy_read(bam1_t *b, hts_base_mod_state *st, hts_pos_t beg, hts_pos_t end, const uint8_t *mask,
                        const char *ref, int32_t *depth, callvec_t *cv) {
    if (bam_parse_basemod(b, st) < 0) return -1;
    int ntypes = 0, *types = bam_mods_recorded(st, &ntypes);
    int has_m = 0, has_h = 0, imp_m = 0;
    for (int k = 0; k < ntypes; k++) {
        int strand, implicit; char canon;
        if (bam_mods_query_type(st, types[k], &strand, &implicit, &canon) < 0) continue;
        if (canon != 'C' || strand != 0) continue;
        if (types[k] == 'm') { has_m = 1; imp_m = implicit; }
        if (types[k] == 'h') has_h = 1;
    }
    if (!has_m) return 0;
    int rev = bam_is_rev(b), qlen = b->core.l_qseq;
    const uint8_t *seq = bam_get_seq(b);
    uint64_t uid = name_uid(bam_get_qname(b));
    hts_base_mod mods[8]; int qp, k, ne = 0, ncap = 1024;
    int *eq = malloc(ncap * sizeof(int)), *em = malloc(ncap * sizeof(int)), *eh = malloc(ncap * sizeof(int));
    while ((k = bam_next_basemod(b, st, mods, 8, &qp)) > 0) {
        int qm = -1, qh = -1;
        for (int j = 0; j < k && j < 8; j++) {
            if (mods[j].canonical_base != 'C' || mods[j].strand != 0) continue;
            if (mods[j].modified_base == 'm') qm = mods[j].qual; else if (mods[j].modified_base == 'h') qh = mods[j].qual;
        }
        if (qm < 0 && qh < 0) continue;
        if (ne == ncap) { ncap *= 2; eq = realloc(eq, ncap * sizeof(int)); em = realloc(em, ncap * sizeof(int)); eh = realloc(eh, ncap * sizeof(int)); }
        eq[ne] = qp; em[ne] = qm; eh[ne] = qh; ne++;
    }
    const uint32_t *cig = bam_get_cigar(b);
    hts_pos_t r = b->core.pos; int q = 0, ei = 0;
    for (uint32_t c = 0; c < b->core.n_cigar; c++) {
        int op = bam_cigar_op(cig[c]), len = bam_cigar_oplen(cig[c]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            for (int i = 0; i < len; i++) {
                hts_pos_t rp = r + i; int qq = q + i;
                if (rp + 1 < beg || rp - 1 >= end) continue;
                if (depth && rp >= beg && rp < end) depth[rp - beg]++;
                int base = bam_seqi(seq, qq), site; hts_pos_t pos;
                if (!rev) {
                    if (base != 2) continue;
                    if (ref) site = ref[rp - beg + 1] == 'C' && ref[rp - beg + 2] == 'G';
                    else site = qq + 1 < qlen && bam_seqi(seq, qq + 1) == 4;
                    pos = rp;
                } else {
                    if (base != 4) continue;
                    if (ref) site = ref[rp - beg + 1] == 'G' && ref[rp - beg] == 'C';
                    else site = qq > 0 && bam_seqi(seq, qq - 1) == 2;
                    pos = rp - 1;
                }
                if (!site || pos < beg || pos >= end || (mask && !mask[pos - beg])) continue;
                while (ei < ne && eq[ei] < qq) ei++;
                int qm = -1, qh = -1;
                if (ei < ne && eq[ei] == qq) { qm = em[ei]; qh = eh[ei]; }
                double ml, hl;
                if (qm >= 0) ml = (qm + 0.5) / 256.0; else if (imp_m) ml = 0; else continue;
                if (qh >= 0) hl = (qh + 0.5) / 256.0; else hl = 0;
                cv_push(cv, pos, uid, call_status(ml, hl), call_class(ml, hl, has_h));
            }
        }
        if (bam_cigar_type(op) & 1) q += len;
        if (bam_cigar_type(op) & 2) r += len;
    }
    free(eq); free(em); free(eh);
    return 1;
}

/* ---- output ---- */
static int ent_decimals = 4;          /* -1 = full precision, as Perl prints (%.15g) */
static int ent_window = 1;            /* 0 = --no-entropy: beta only, no window of 3 CpGs */
static int ent_min_cov = 3;           /* confident calls needed to write a CpG (ParseModkit.pl: > 2) */
static void put_num(FILE *f, double x) {
    if (ent_decimals < 0) { if (x == 0) x = 0; fprintf(f, "%.15g", x); return; }
    char buf[64]; snprintf(buf, sizeof buf, "%.*f", ent_decimals, x);
    char *e = buf + strlen(buf) - 1;                     /* 0.5000 -> 0.5, 1.0000 -> 1 */
    if (strchr(buf, '.')) { while (*e == '0') *e-- = 0; if (*e == '.') *e = 0; }
    fputs(strcmp(buf, "-0") ? buf : "0", f);
}

/* counts of one CpG, one call per read (a read seen twice at the same CpG counts once) */
static void add_call(wpos_t *p, const call_t *c) {
    if (c->st >= 0) { p->cov++; p->meth += c->st; }
    int k = c->cls & 3;
    if (c->cls & CLS_H_SEEN) p->hseen = 1;
    if (k == CLS_M) p->nm++; else if (k == CLS_H) p->nh++; else if (k == CLS_C) p->nc++;
}

/* end of a row: beta, coverage and, with --split-5hmc, beta_5mC, beta_5hmC, n */
static void put_tail(FILE *f, const wpos_t *p) {
    put_num(f, (double) p->meth / p->cov); fprintf(f, "\t%d", p->cov);
    if (ent_split) {
        int n = p->nm + p->nh + p->nc;
        if (n == 0) fputs("\tNA\tNA\t0", f);
        else {
            fputc('\t', f); put_num(f, (double) p->nm / n); fputc('\t', f);
            if (p->hseen) put_num(f, (double) p->nh / n); else fputs("NA", f);
            fprintf(f, "\t%d", n);
        }
    }
    fputc('\n', f);
}

/* window of 3 consecutive CpGs: output for w[0] */
static int window_out(FILE *f, const char *chr, wpos_t *w) {
    long cnt[8] = {0}; long sum = 0;
    int i = 0, j = 0, k = 0;
    while (i < w[0].n && j < w[1].n && k < w[2].n) {
        uint64_t a = w[0].uid[i], b2 = w[1].uid[j], c = w[2].uid[k];
        if (a == b2 && b2 == c) {
            int s0 = w[0].st[i], s1 = w[1].st[j], s2 = w[2].st[k];
            if (s0 >= 0 && s1 >= 0 && s2 >= 0) { cnt[s0 * 4 + s1 * 2 + s2]++; sum++; }
            i++; j++; k++;
        } else {
            uint64_t mx = a > b2 ? a : b2; if (c > mx) mx = c;
            if (a < mx) i++;
            if (b2 < mx) j++;
            if (c < mx) k++;
        }
    }
    if (sum > 2 && w[0].cov >= ent_min_cov) {
        double e = 0;
        for (int p = 0; p < 8; p++) if (cnt[p]) { double fr = (double) cnt[p] / sum; e += fr * (log(fr) / log(2)); }
        e = -e * (1.0 / 3);
        fprintf(f, "%s\t%lld\t", chr, (long long) w[0].pos);
        put_num(f, e); fprintf(f, "\t%ld\t", sum);
        put_tail(f, &w[0]);
        return 1;
    }
    return 0;
}

static void *entropy_worker(void *arg) {
    (void) arg;
    samFile *in = sam_open(bam_path, "r");
    sam_hdr_t *h = sam_hdr_read(in);
    hts_idx_t *idx = sam_index_load(in, bam_path);
    faidx_t *fai = ref_path ? fai_load(ref_path) : NULL;
    bam1_t *b = bam_init1();
    hts_base_mod_state *st = hts_base_mod_state_alloc();
    callvec_t cv = {0};
    uint8_t *mask = malloc(MAX_SPAN + 2);
    int32_t *depth = fai ? NULL : malloc((MAX_SPAN + 2) * sizeof(int32_t));
    char **rs = malloc(MAX_REG_PER_CHUNK * sizeof(char *));
    for (int i = 0; i < MAX_REG_PER_CHUNK; i++) rs[i] = malloc(300);
    wpos_t w[4]; memset(w, 0, sizeof w);
    for (;;) {
        pthread_mutex_lock(&lock); int ci = next_chrom++; pthread_mutex_unlock(&lock);
        if (ci >= nechr) break;
        chrom_t *cr = &echr[ci];
        const char *chr = sam_hdr_tid2name(h, cr->tid);
        FILE *f = fopen(tmp_names[ci], "wb");
        if (!f) die("cannot write", tmp_names[ci]);
        int nw = 0; long rows = 0;
        for (int c = cr->c0; c < cr->c1; c++) {
            chunk_t *ch = &chunks[c];
            hts_pos_t span = ch->end - ch->beg;
            int nr = ch->r1 - ch->r0, use_mask = nr > 0;
            if (use_mask) {
                memset(mask, 0, span);
                for (int r = ch->r0; r < ch->r1; r++) {
                    memset(mask + (regs[r].beg - ch->beg), 1, regs[r].end - regs[r].beg);
                    reg_string(regs[r].tid, regs[r].beg, regs[r].end, h, rs[r - ch->r0], 300);
                }
            } else { nr = 1; reg_string(ch->tid, ch->beg, ch->end, h, rs[0], 300); }
            char *ref = NULL;
            if (fai) {
                hts_pos_t rl = 0, b0 = ch->beg > 0 ? ch->beg - 1 : 0;
                char *s = faidx_fetch_seq64(fai, chr, b0, ch->end, &rl);   /* [beg-1, end] */
                ref = calloc(span + 3, 1);
                if (s) {
                    int off = ch->beg > 0 ? 0 : 1;
                    for (hts_pos_t i = 0; i < rl && i + off < span + 2; i++) ref[i + off] = toupper((unsigned char) s[i]);
                    free(s);
                }
            }
            long nbad = 0, nok = 0;
            /* without reference: first pass finds the CpGs of the sample (sites seen as CpG by at
               least min_frac of the reads covering them) and builds a reference with only those */
            for (int pass = ref ? 1 : 0; pass < 2; pass++) {
                if (pass == 0) memset(depth, 0, span * sizeof(int32_t));
                cv.n = 0; nbad = nok = 0;
                hts_itr_t *it = sam_itr_regarray(idx, h, rs, nr);
                if (it) {
                    while (sam_itr_next(in, it, b) >= 0) {
                        if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
                        if (entropy_read(b, st, ch->beg, ch->end, use_mask ? mask : NULL, ref, pass == 0 ? depth : NULL, &cv) < 0) nbad++; else nok++;
                    }
                    hts_itr_destroy(it);
                }
                if (pass == 0) {
                    qsort(cv.v, cv.n, sizeof(call_t), cmp_call);
                    ref = calloc(span + 3, 1);
                    for (size_t i = 0; i < cv.n; ) {
                        size_t j = i; while (j < cv.n && cv.v[j].pos == cv.v[i].pos) j++;
                        hts_pos_t x = cv.v[i].pos - ch->beg;
                        if ((double) (j - i) >= min_frac * depth[x]) { ref[x + 1] = 'C'; ref[x + 2] = 'G'; }
                        i = j;
                    }
                }
            }
            free(ref);
            qsort(cv.v, cv.n, sizeof(call_t), cmp_call);
            for (size_t i = 0; i < cv.n; ) {
                size_t j = i; while (j < cv.n && cv.v[j].pos == cv.v[i].pos) j++;
                int n = (int) (j - i);
                if (!ent_window) {                               /* beta only: write the CpG now */
                    wpos_t p; memset(&p, 0, sizeof p);
                    for (size_t k = i; k < j; k++) {
                        if (k > i && cv.v[k].uid == cv.v[k - 1].uid) continue;
                        add_call(&p, &cv.v[k]);
                    }
                    if (p.cov >= ent_min_cov && p.cov > 0) {
                        fprintf(f, "%s\t%lld\tNA\tNA\t", chr, (long long) cv.v[i].pos);
                        put_tail(f, &p);
                        rows++;
                    }
                    i = j; continue;
                }
                if (nw == 3) {                                   /* new CpG: the oldest one is complete */
                    rows += window_out(f, chr, w);
                    wpos_t t = w[0]; w[0] = w[1]; w[1] = w[2]; w[2] = t; nw = 2;
                }
                wpos_t *p = &w[nw++];
                if (p->cap < n) { p->cap = n; p->uid = realloc(p->uid, n * sizeof(uint64_t)); p->st = realloc(p->st, n); }
                p->pos = cv.v[i].pos; p->n = 0; p->cov = p->meth = p->nm = p->nh = p->nc = p->hseen = 0;
                for (size_t k = i; k < j; k++) {
                    if (p->n && p->uid[p->n - 1] == cv.v[k].uid) continue;     /* same read twice */
                    p->uid[p->n] = cv.v[k].uid; p->st[p->n] = cv.v[k].st; p->n++;
                    add_call(p, &cv.v[k]);
                }
                i = j;
            }
            pthread_mutex_lock(&lock);
            n_reads_used += nok; n_reads_bad += nbad; done_chunks++;
            if (done_chunks % (nchunks / 100 + 1) == 0 || done_chunks == nchunks)
                fprintf(stderr, "PROGRESS %d/%d\n", done_chunks, nchunks);
            pthread_mutex_unlock(&lock);
        }
        pthread_mutex_lock(&lock); ent_rows += rows; pthread_mutex_unlock(&lock);
        fclose(f);                    /* the last CpGs of a chromosome have no complete window (as ParseModkit.pl) */
    }
    for (int i = 0; i < 4; i++) { free(w[i].uid); free(w[i].st); }
    for (int i = 0; i < MAX_REG_PER_CHUNK; i++) free(rs[i]);
    free(rs); free(mask); free(depth); free(cv.v);
    if (fai) fai_destroy(fai);
    hts_base_mod_state_free(st); bam_destroy1(b); hts_idx_destroy(idx); sam_hdr_destroy(h); sam_close(in);
    return NULL;
}

/* chromosomes 1-22, X, Y in this order (with or without "chr") */
static int canon_tids(sam_hdr_t *h, int *tids) {
    int n = 0; char nm[16];
    for (int i = 1; i <= 24; i++) {
        const char *base = i <= 22 ? NULL : (i == 23 ? "X" : "Y");
        if (base) snprintf(nm, sizeof nm, "chr%s", base); else snprintf(nm, sizeof nm, "chr%d", i);
        int t = sam_hdr_name2tid(h, nm);
        if (t < 0) t = sam_hdr_name2tid(h, nm + 3);
        if (t >= 0) tids[n++] = t;
    }
    return n;
}

static int run_entropy(const char *rfile, const char *ofile, int nth) {
    samFile *in = sam_open(bam_path, "r");
    if (!in) die("cannot open BAM", bam_path);
    sam_hdr_t *h = sam_hdr_read(in);
    if (!h) die("cannot read BAM header", bam_path);
    hts_idx_t *idx = sam_index_load(in, bam_path);
    if (!idx) die("BAM index (.bai) not found next to", bam_path);
    if (ref_path) { faidx_t *f = fai_load(ref_path); if (!f) die("cannot open reference FASTA (needs .fai)", ref_path); fai_destroy(f); }
    int tids[24], nt = canon_tids(h, tids);
    if (!nt) die("no chromosome 1-22, X, Y in", bam_path);
    if (rfile) read_regions(rfile, h);
    /* chunks grouped by chromosome, in canonical order */
    chunk_t *all = malloc(sizeof(chunk_t) * (rfile ? nchunks + 1 : 4000));
    int na = 0, cap = rfile ? nchunks + 1 : 4000;
    echr = malloc(nt * sizeof(chrom_t)); nechr = 0;
    for (int k = 0; k < nt; k++) {
        int t = tids[k], c0 = na;
        if (rfile) {
            for (int c = 0; c < nchunks; c++) if (chunks[c].tid == t) all[na++] = chunks[c];
        } else {
            hts_pos_t len = sam_hdr_tid2len(h, t);
            for (hts_pos_t s = 0; s < len; s += MAX_SPAN) {
                if (na == cap) { cap *= 2; all = realloc(all, cap * sizeof(chunk_t)); }
                chunk_t c = { t, s, s + MAX_SPAN < len ? s + MAX_SPAN : len, 0, 0 };
                all[na++] = c;
            }
        }
        if (na > c0) { echr[nechr].tid = t; echr[nechr].c0 = c0; echr[nechr].c1 = na; nechr++; }
    }
    if (rfile) free(chunks);
    chunks = all; nchunks = na;
    fprintf(stderr, "%s: %s, %d chromosomes, %d chunks, CpG sites from %s, minimum coverage %d, %s%s\n",
            ent_window ? "entropy mode" : "beta mode (no entropy)", rfile ? "model regions" : "whole genome",
            nechr, nchunks, ref_path ? "the reference" : "the reads", ent_min_cov,
            ent_sum_h ? "beta = 5mC + 5hmC" : "beta = 5mC", ent_split ? ", 5mC and 5hmC also split" : "");
    tmp_names = malloc(nechr * sizeof(char *));
    for (int i = 0; i < nechr; i++) {
        tmp_names[i] = malloc(strlen(ofile) + 32);
        sprintf(tmp_names[i], "%s.part%02d.tmp", ofile, i);
    }
    pthread_t *tt = malloc(nth * sizeof(pthread_t));
    for (int i = 0; i < nth; i++) pthread_create(&tt[i], NULL, entropy_worker, NULL);
    for (int i = 0; i < nth; i++) pthread_join(tt[i], NULL);
    FILE *out = fopen(ofile, "wb");
    if (!out) die("cannot write", ofile);
    char *buf = malloc(1 << 20);
    for (int i = 0; i < nechr; i++) {
        FILE *f = fopen(tmp_names[i], "rb"); size_t n;
        if (f) { while ((n = fread(buf, 1, 1 << 20, f)) > 0) fwrite(buf, 1, n, out); fclose(f); }
        remove(tmp_names[i]);
    }
    fclose(out); free(buf);
    fprintf(stderr, "done: %ld rows, %ld reads used, %ld reads skipped (methylation tags not readable)\n", ent_rows, n_reads_used, n_reads_bad);
    hts_idx_destroy(idx); sam_hdr_destroy(h); sam_close(in);
    return 0;
}

int main(int argc, char **argv) {
    const char *rfile = NULL, *ofile = NULL; int nth = 4;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-r") && i + 1 < argc) rfile = argv[++i];
        else if (!strcmp(argv[i], "-o") && i + 1 < argc) ofile = argv[++i];
        else if (!strcmp(argv[i], "-t") && i + 1 < argc) nth = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--entropy")) ;                    /* the only mode, kept for old command lines */
        else if (!strcmp(argv[i], "--full-precision")) ent_decimals = -1;
        else if (!strcmp(argv[i], "--no-entropy")) ent_window = 0;
        else if (!strcmp(argv[i], "--sum-5hmc")) ent_sum_h = 1;
        else if (!strcmp(argv[i], "--split-5hmc")) ent_split = 1;
        else if (!strcmp(argv[i], "--min-cov") && i + 1 < argc) ent_min_cov = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--ref") && i + 1 < argc) ref_path = argv[++i];
        else if (!strcmp(argv[i], "--min-cpg-frac") && i + 1 < argc) min_frac = atof(argv[++i]);
        else if (!strcmp(argv[i], "--version")) { printf("magibu-meth 0.3 (htslib %s)\n", hts_version()); return 0; }
        else if (argv[i][0] != '-' && !bam_path) bam_path = argv[i];
        else die("unknown argument", argv[i]);
    }
    if (!bam_path || !ofile) {
        fprintf(stderr, "usage: magibu-meth -o out.entropy.file.tsv [-r regions.bed] [--ref ref.fasta] [-t threads]\n"
                        "       [--min-cpg-frac F] [--full-precision] [--no-entropy] [--min-cov N]\n"
                        "       [--sum-5hmc] [--split-5hmc] in.bam\n"
                        "  --sum-5hmc    beta (column 5) = 5mC + 5hmC; default 5mC only\n"
                        "  --split-5hmc  also write beta_5mC, beta_5hmC and n (columns 7-9)\n");
        return 1;
    }
    if (nth < 1) nth = 1;
    hts_set_log_level(HTS_LOG_OFF);
    return run_entropy(rfile, ofile, nth);
}
