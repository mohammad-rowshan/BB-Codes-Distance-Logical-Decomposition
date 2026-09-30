/*
 * cluster_lb.c
 *
 * Anchored cluster search for bivariate bicycle (BB) codes, the lower-bound
 * engine behind Table "cluster certificates" of the manuscript
 * "Logical Operator Decomposition for Distance Analysis of Bivariate Bicycle
 * Codes" (Rowshan, Devitt).
 *
 * What it does, in one paragraph: we grow a support T one qubit at a time,
 * starting from an anchor qubit, and at every step we only branch on the
 * qubits of the lowest-index unsatisfied X-check. Lemma "syndrome
 * connectivity" (Sec. V) says every proper subset of a minimum-weight logical
 * has nonzero syndrome, so this never loses a minimum-weight logical.
 * Translations fix the anchor (left origin, plus right origin for supports
 * living only on the right block), see Theorem "exactness of the cluster search".
 * If the search at radius W ends without a hit, that is a proof that d_Z >= W+1,
 * and d_X = d_Z for every BB code (Proposition "X/Z symmetry").
 *
 * The idea of growing clusters in the Tanner graph goes back to
 * Dumer, Kovalev and Pryadko, IEEE Trans. Inf. Theory 63(7), 2017.
 * Everything here is exact integer / bit arithmetic, no floating point.
 *
 * Build:   gcc -O3 -march=native -o cluster_lb cluster_lb.c
 *
 * Usage:   ./cluster_lb [options] l m W "a(x,y)" "b(x,y)"
 *   polynomials are written like  "x^3+y+y^2"  or  "1+x^2+x^7"  or "x^2*y^3"
 *
 * Options:
 *   -e        enumerate: do not stop at the first hit, print every hit
 *             (used for the minimum-weight census, Table "census")
 *   -c        colon mode: a leaf only counts if pi(z) = v + (a) is nonzero,
 *             and zero-syndrome leaves with pi = 0 are pruned
 *             (Corollary "colon connectivity"); certifies d_C >= W+1
 *   -k        kernel mode: any nonzero element of K counts, stabilizers too.
 *             Gives w_K, the minimum nonzero weight of K, needed by
 *             Proposition "irreducible range"
 *   -n NODES  node budget; if exceeded we say so and certify nothing
 *
 * Output (one item per line, easy to parse from Python):
 *   # info lines
 *   LOGICAL <weight> <q_1> ... <q_w>     qubit q < N is left block, q >= N right
 *   RESULT found|none|enumerated|budget ...
 *
 * Qubit indexing matches the manuscript: monomial x^i y^j <-> index i*m + j,
 * left block first (0..N-1), right block second (N..2N-1).
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

/* ---------------- small helpers ---------------- */

static int L_, M_, N_, n_;          /* l, m, N = l*m, n = 2N */
static int W_;                      /* search radius */

static int idx(int i, int j) {      /* x^i y^j -> column index, exponents mod l, m */
    i %= L_; if (i < 0) i += L_;
    j %= M_; if (j < 0) j += M_;
    return i * M_ + j;
}

/* Parse "x^3+y+y^2" style polynomials into a list of (i,j) exponents.
 * Repeated monomials cancel, since we are over F_2. Returns number of terms. */
static int parse_poly(const char *s, int (*terms)[2], int maxterms) {
    int cnt = 0;
    const char *p = s;
    while (*p) {
        int ei = 0, ej = 0, seen = 0;
        /* one term: a product of 1, x^e, y^e with optional '*' */
        while (*p && *p != '+') {
            if (isspace((unsigned char)*p) || *p == '*') { p++; continue; }
            if (*p == '1' && !seen) { p++; seen = 1; continue; }
            if (*p == 'x' || *p == 'y') {
                char v = *p++; int e = 1;
                if (*p == '^') { p++; e = (int)strtol(p, (char **)&p, 10); }
                if (v == 'x') ei += e; else ej += e;
                seen = 1; continue;
            }
            fprintf(stderr, "cannot parse polynomial near '%s'\n", p); exit(2);
        }
        if (*p == '+') p++;
        if (!seen) continue;
        ei %= L_; ej %= M_;
        /* XOR the term into the list (x + x = 0 over F_2) */
        int k, dup = -1;
        for (k = 0; k < cnt; k++) if (terms[k][0] == ei && terms[k][1] == ej) dup = k;
        if (dup >= 0) { terms[dup][0] = terms[cnt-1][0]; terms[dup][1] = terms[cnt-1][1]; cnt--; }
        else {
            if (cnt >= maxterms) { fprintf(stderr, "too many terms\n"); exit(2); }
            terms[cnt][0] = ei; terms[cnt][1] = ej; cnt++;
        }
    }
    return cnt;
}

/* ---------------- GF(2) row space with highest-bit pivots ---------------- */
/* Used twice: for S = im[B;A] (n bits) and, in colon mode, for (a) (N bits). */

typedef struct {
    int nbits, words, rows;
    uint64_t *basis;   /* rows * words */
    int *piv;
} space_t;

static void space_init(space_t *S, int nbits, int maxrows) {
    S->nbits = nbits; S->words = (nbits + 63) / 64; S->rows = 0;
    S->basis = calloc((size_t)maxrows * S->words, sizeof(uint64_t));
    S->piv = calloc(maxrows, sizeof(int));
}

static void space_reduce(const space_t *S, uint64_t *v) {
    for (int r = 0; r < S->rows; r++) {
        int p = S->piv[r];
        if (v[p >> 6] >> (p & 63) & 1ULL) {
            const uint64_t *b = S->basis + (size_t)r * S->words;
            for (int w = 0; w < S->words; w++) v[w] ^= b[w];
        }
    }
}

static int is_zero(const uint64_t *v, int words) {
    for (int w = 0; w < words; w++) if (v[w]) return 0;
    return 1;
}

/* add v to the space (fully reduced echelon form, so reduce() is canonical) */
static void space_add(space_t *S, uint64_t *v) {
    space_reduce(S, v);
    int p = -1;
    for (int w = S->words - 1; w >= 0 && p < 0; w--)
        if (v[w]) p = w * 64 + 63 - __builtin_clzll(v[w]);
    if (p < 0) return;                      /* already in the span */
    for (int r = 0; r < S->rows; r++) {     /* clear the new pivot elsewhere */
        uint64_t *b = S->basis + (size_t)r * S->words;
        if (b[p >> 6] >> (p & 63) & 1ULL)
            for (int w = 0; w < S->words; w++) b[w] ^= v[w];
    }
    memcpy(S->basis + (size_t)S->rows * S->words, v, S->words * sizeof(uint64_t));
    S->piv[S->rows++] = p;
}

/* ---------------- the code ---------------- */

static int wa, wb, ta[32][2], tb[32][2];
static int *qchk, *qdeg;            /* checks touched by each qubit (<= max(wa,wb)) */
static int *cq, *cdeg;              /* qubits in each check (<= wa+wb) */
static int CMAX, QMAX;              /* max check weight, max qubit degree */
static space_t Sspace, Aideal;      /* stabilizers, and the ideal (a) for colon mode */

/* search state */
static uint64_t *syn;               /* syndrome bitset over the N X-checks */
static int swords, nunsat;
static char *inT;
static int *T, tsz;
static int rightonly;
static int mode_enum, mode_colon, mode_kernel;
static long long nodes, budget = -1, hits;
static int stop_all, budget_hit, best_found;

static void flip(int q) {
    inT[q] ^= 1;
    for (int k = 0; k < qdeg[q]; k++) {
        int c = qchk[q * QMAX + k];
        syn[c >> 6] ^= 1ULL << (c & 63);
        nunsat += (syn[c >> 6] >> (c & 63) & 1ULL) ? 1 : -1;
    }
}

static int T_in_S(void) {
    uint64_t v[Sspace.words];
    memset(v, 0, sizeof v);
    for (int t = 0; t < tsz; t++) v[T[t] >> 6] |= 1ULL << (T[t] & 63);
    space_reduce(&Sspace, v);
    return is_zero(v, Sspace.words);
}

/* pi(z) = v + (a): is the right-block part of T outside the ideal (a)? */
static int T_pi_nonzero(void) {
    uint64_t v[Aideal.words];
    memset(v, 0, sizeof v);
    for (int t = 0; t < tsz; t++) if (T[t] >= N_) {
        int q = T[t] - N_; v[q >> 6] |= 1ULL << (q & 63);
    }
    space_reduce(&Aideal, v);
    return !is_zero(v, Aideal.words);
}

static void report_hit(void) {
    hits++;
    if (!best_found || tsz < best_found) best_found = tsz;
    printf("LOGICAL %d", tsz);
    for (int t = 0; t < tsz; t++) printf(" %d", T[t]);
    printf("\n");
    if (!mode_enum) stop_all = 1;
}

static void grow(void) {
    if (stop_all) return;
    if (budget >= 0 && nodes >= budget) { budget_hit = 1; stop_all = 1; return; }
    nodes++;

    if (nunsat == 0) {
        /* T is in K. Decide whether it counts as a hit, otherwise prune.
         * Pruning is safe: a zero-syndrome proper subset of a minimum-weight
         * target cannot exist (Lemma "syndrome connectivity", or the colon
         * version, or minimality of w_K in kernel mode). */
        if (mode_kernel) report_hit();
        else if (mode_colon) { if (T_pi_nonzero()) report_hit(); }
        else if (!T_in_S()) report_hit();
        return;
    }
    /* parity prune: each extra qubit fixes at most QMAX unsatisfied checks */
    if (tsz + (nunsat + QMAX - 1) / QMAX > W_) return;

    /* canonical choice: the lowest-index unsatisfied check */
    int c = -1;
    for (int w = 0; w < swords; w++) if (syn[w]) { c = w * 64 + __builtin_ctzll(syn[w]); break; }

    for (int k = 0; k < cdeg[c]; k++) {
        int q = cq[c * CMAX + k];
        if (inT[q]) continue;
        if (rightonly && q < N_) continue;
        T[tsz++] = q; flip(q);
        grow();
        flip(q); tsz--;
        if (stop_all) return;
    }
}

int main(int argc, char **argv) {
    int ai = 1;
    while (ai < argc && argv[ai][0] == '-' && !isdigit((unsigned char)argv[ai][1])) {
        if (!strcmp(argv[ai], "-e")) mode_enum = 1;
        else if (!strcmp(argv[ai], "-c")) mode_colon = 1;
        else if (!strcmp(argv[ai], "-k")) mode_kernel = 1;
        else if (!strcmp(argv[ai], "-n") && ai + 1 < argc) budget = atoll(argv[++ai]);
        else { fprintf(stderr, "unknown option %s\n", argv[ai]); return 2; }
        ai++;
    }
    if (argc - ai != 5) {
        fprintf(stderr, "usage: %s [-e] [-c] [-k] [-n budget] l m W \"a\" \"b\"\n", argv[0]);
        return 2;
    }
    L_ = atoi(argv[ai]); M_ = atoi(argv[ai+1]); W_ = atoi(argv[ai+2]);
    N_ = L_ * M_; n_ = 2 * N_;
    wa = parse_poly(argv[ai+3], ta, 32);
    wb = parse_poly(argv[ai+4], tb, 32);
    if (wa == 0 || wb == 0) { fprintf(stderr, "a and b must be nonzero\n"); return 2; }

    CMAX = wa + wb; QMAX = wa > wb ? wa : wb;
    qchk = calloc((size_t)n_ * QMAX, sizeof(int)); qdeg = calloc(n_, sizeof(int));
    cq = calloc((size_t)N_ * CMAX, sizeof(int)); cdeg = calloc(N_, sizeof(int));

    /* H_X = [A | B] with column-vector convention (Remark "stabilizer
     * convention"): left qubit j sits in checks j + supp(a), right qubit j
     * in checks j + supp(b). */
    for (int j = 0; j < N_; j++) {
        int ji = j / M_, jj = j % M_;
        for (int k = 0; k < wa; k++) {
            int g = idx(ji + ta[k][0], jj + ta[k][1]);
            qchk[j * QMAX + qdeg[j]++] = g; cq[g * CMAX + cdeg[g]++] = j;
        }
        for (int k = 0; k < wb; k++) {
            int g = idx(ji + tb[k][0], jj + tb[k][1]);
            qchk[(N_ + j) * QMAX + qdeg[N_ + j]++] = g; cq[g * CMAX + cdeg[g]++] = N_ + j;
        }
    }

    /* S = {(b r, a r)}: generators r = x^i y^j (Sec. II-B) */
    space_init(&Sspace, n_, N_);
    {
        int words = Sspace.words;
        uint64_t v[words];
        for (int r = 0; r < N_; r++) {
            int ri = r / M_, rj = r % M_;
            memset(v, 0, sizeof v);
            for (int k = 0; k < wb; k++) { int q = idx(ri + tb[k][0], rj + tb[k][1]); v[q >> 6] ^= 1ULL << (q & 63); }
            for (int k = 0; k < wa; k++) { int q = N_ + idx(ri + ta[k][0], rj + ta[k][1]); v[q >> 6] ^= 1ULL << (q & 63); }
            space_add(&Sspace, v);
        }
    }
    /* for BB codes rank H_Z = rank H_X = dim S, so k = n - 2 dim S */
    int k = n_ - 2 * Sspace.rows;

    if (mode_colon) {           /* the ideal (a) = span{a x^r}, as N-bit words */
        space_init(&Aideal, N_, N_);
        uint64_t v[Aideal.words];
        for (int r = 0; r < N_; r++) {
            int ri = r / M_, rj = r % M_;
            memset(v, 0, sizeof v);
            for (int t = 0; t < wa; t++) { int q = idx(ri + ta[t][0], rj + ta[t][1]); v[q >> 6] ^= 1ULL << (q & 63); }
            space_add(&Aideal, v);
        }
    }

    printf("# BB code l=%d m=%d n=%d k=%d dimS=%d check_weight=%d radius=%d mode=%s%s\n",
           L_, M_, n_, k, Sspace.rows, CMAX, W_,
           mode_kernel ? "kernel" : (mode_colon ? "colon" : "logical"),
           mode_enum ? "+enumerate" : "");
    if (k == 0 && !mode_kernel) { printf("RESULT nologicals k=0\n"); return 0; }

    swords = (N_ + 63) / 64;
    syn = calloc(swords, sizeof(uint64_t));
    inT = calloc(n_, 1);
    T = calloc(n_ + 1, sizeof(int));

    clock_t t0 = clock();
    /* two anchors are enough (proof of Theorem "exactness of the cluster search") */
    int anchors[2] = {0, N_};
    for (int A = 0; A < 2 && !stop_all; A++) {
        rightonly = (A == 1);
        T[0] = anchors[A]; tsz = 1; flip(anchors[A]);
        grow();
        flip(anchors[A]); tsz = 0;
    }
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;

    if (budget_hit)
        printf("RESULT budget radius=%d nodes=%lld time=%.2f\n", W_, nodes, secs);
    else if (mode_enum)
        printf("RESULT enumerated radius=%d hits=%lld min_weight=%d nodes=%lld time=%.2f\n",
               W_, hits, best_found, nodes, secs);
    else if (hits)
        printf("RESULT found radius=%d weight=%d nodes=%lld time=%.2f\n", W_, best_found, nodes, secs);
    else
        printf("RESULT none radius=%d nodes=%lld time=%.2f  => %s >= %d\n", W_, nodes, secs,
               mode_kernel ? "w_K" : (mode_colon ? "d_C" : "d"), W_ + 1);
    return 0;
}
