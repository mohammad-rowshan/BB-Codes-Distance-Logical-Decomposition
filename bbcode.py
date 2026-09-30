#!/usr/bin/env python3
"""
bbcode.py -- decomposition and distance bounds for bivariate bicycle codes.

Companion code for
  M. Rowshan and S. Devitt, "Logical Operator Decomposition for Distance
  Analysis of Bivariate Bicycle Codes".

References in the comments give the number used in the manuscript together
with its LaTeX label in square brackets, e.g. "Lemma 2 [lem:conn]", so they
can be matched even if the numbering moves in a later version.

Everything is exact linear algebra over F_2. Vectors are Python ints used as
bitsets: bit (i*m + j) is the coefficient of x^i y^j, and a two-block vector
(u, v) is stored as u | (v << N), left block in the low bits (Sec. II-B).

Shell use, for any BB code:

    python3 bbcode.py -l 12 -m 6 -a "x^3+y+y^2" -b "y^3+x+x^2"
    python3 bbcode.py -l 12 -m 12 -a "x^3+y^2+y^7" -b "y^3+x+x^2" --seconds 60

It prints n, k, the exact-sequence profile of Theorem 1 [thm:exact] and a
pair of bounds d_L <= d <= d_U. Upper bounds are explicit logicals found by the
screens of Algorithm "search for low-weight logical witnesses" [alg:screen].
Lower bounds come from the anchored cluster search [alg:cluster], justified by
Lemma 2 [lem:conn] and Theorem 4 [thm:cluster], run at increasing radius until
the bounds meet (Sec. V-D) or the time budget runs out. When the compiled C
engine ./cluster_lb is present it is used; otherwise the pure-Python search
below runs the same algorithm with the same node counts.
"""

import argparse
import os
import random
import re
import shutil
import subprocess
import sys
import time
from itertools import combinations

popcount = int.bit_count          # Python >= 3.10


# ---------------------------------------------------------------------------
# GF(2) helpers.  An echelon basis is kept as a dict {pivot_bit: row}, fully
# reduced, so reduce() returns a canonical coset representative. The census
# (Table V [tab:census]) relies on this: two logicals are in the same class
# of K/S iff they reduce to the same vector modulo S.
# ---------------------------------------------------------------------------

class Span:
    """Row space over F_2 in fully reduced echelon form."""

    def __init__(self, rows=()):
        self.rows = {}                      # pivot -> row
        for r in rows:
            self.add(r)

    def reduce(self, v):
        for p, r in self.rows.items():
            if v >> p & 1:
                v ^= r
        return v

    def add(self, v):
        v = self.reduce(v)
        if not v:
            return False
        p = v.bit_length() - 1
        for q in list(self.rows):           # keep it fully reduced
            if self.rows[q] >> p & 1:
                self.rows[q] ^= v
        self.rows[p] = v
        return True

    def __contains__(self, v):
        return self.reduce(v) == 0

    def __len__(self):
        return len(self.rows)

    def basis(self):
        return list(self.rows.values())


def nullspace(images, nbits_in):
    """Kernel of the linear map sending basis vector e_i to images[i].

    Stack [image | identity] and eliminate on the image part; rows whose
    image part vanishes carry kernel vectors in the identity part. This is
    the "Gaussian elimination" of Sec. II-C, e.g. ann(a) = ker A.
    """
    rows = [(img << nbits_in) | (1 << i) for i, img in enumerate(images)]
    sp = Span(rows)
    mask = (1 << nbits_in) - 1
    return [r & mask for r in sp.basis() if r >> nbits_in == 0]


# ---------------------------------------------------------------------------
# The ring R = F_2[x,y]/(x^l - 1, y^m - 1)  (Sec. II-A)
# ---------------------------------------------------------------------------

def parse_poly(s, l, m):
    """'x^3+y+y^2' -> list of exponent pairs; repeated terms cancel mod 2."""
    terms = set()
    for raw in s.replace(" ", "").split("+"):
        if not raw:
            continue
        if raw == "1":
            i = j = 0
        else:
            ex = re.fullmatch(r"(x(\^\d+)?)?\*?(y(\^\d+)?)?", raw)
            if not ex or not raw.strip("*"):
                raise ValueError("cannot parse term %r" % raw)
            i = (int(ex.group(2)[1:]) if ex.group(2) else 1) if ex.group(1) else 0
            j = (int(ex.group(4)[1:]) if ex.group(4) else 1) if ex.group(3) else 0
        terms ^= {(i % l, j % m)}
    return sorted(terms)


class BBCode:
    """BB code of (a, b) on the l x m torus (Sec. II-B).

    H_X = [A | B] and H_Z = [B^T | A^T], with A the matrix of multiplication
    by a in the column-vector convention of Remark 1 [rem:conv]. Then
    K = {(u,v): a u + b v = 0} and S = {(b r, a r)}.
    """

    def __init__(self, l, m, a, b, name=None):
        self.l, self.m = l, m
        self.N = l * m
        self.n = 2 * self.N
        self.a_terms = parse_poly(a, l, m) if isinstance(a, str) else list(a)
        self.b_terms = parse_poly(b, l, m) if isinstance(b, str) else list(b)
        self.a = self.poly(self.a_terms)
        self.b = self.poly(self.b_terms)
        self.name = name or "BB(l=%d,m=%d)" % (l, m)
        self._build()

    # -- ring arithmetic ---------------------------------------------------
    def idx(self, i, j):
        return (i % self.l) * self.m + (j % self.m)

    def poly(self, terms):
        v = 0
        for i, j in terms:
            v ^= 1 << self.idx(i, j)
        return v

    def monomials(self, f):
        """exponent pairs of the monomials in f"""
        out, k = [], 0
        while f:
            if f & 1:
                out.append(divmod(k, self.m))
            f >>= 1
            k += 1
        return out

    def shift(self, f, di, dj):
        """multiply f by the monomial x^di y^dj, a torus translation"""
        r = 0
        for i, j in self.monomials(f):
            r |= 1 << self.idx(i + di, j + dj)
        return r

    def mul(self, f, g):
        """product in R: cyclic convolution in both indices (Sec. II-A)"""
        if popcount(f) > popcount(g):
            f, g = g, f
        r = 0
        for i, j in self.monomials(f):
            r ^= self.shift(g, i, j)
        return r

    def fmt(self, f):
        """pretty print an element of R, e.g. x^3y^2+y"""
        if not f:
            return "0"
        parts = []
        for i, j in self.monomials(f):
            s = ("x" if i else "") + ("^%d" % i if i > 1 else "") + \
                ("y" if j else "") + ("^%d" % j if j > 1 else "")
            parts.append(s or "1")
        return "+".join(parts)

    # -- Table I [tab:notation] and Table II [tab:glance] -------------------
    def _build(self):
        N = self.N
        mono = [1 << q for q in range(N)]
        self.Acols = [self.mul(self.a, e) for e in mono]     # columns of A
        self.Bcols = [self.mul(self.b, e) for e in mono]     # columns of B
        self.ideal_a = Span(self.Acols)                       # (a) = im A
        self.ideal_b = Span(self.Bcols)                       # (b) = im B
        # S = im [B; A], eq. (KSmat): generators (b x^r, a x^r)
        self.S = Span(self.Bcols[r] | (self.Acols[r] << N) for r in range(N))
        # K = ker [A | B] = ker H_X, eq. (KSmat)
        self.K = nullspace(self.Acols + self.Bcols, self.n)
        self.k = len(self.K) - len(self.S)
        # Tanner-graph data for the cluster search: qubit q -> X-checks it meets
        self.qmask = [0] * self.n
        for q in range(N):
            self.qmask[q] = self.Acols[q]
            self.qmask[N + q] = self.Bcols[q]
        self.check_qubits = [[] for _ in range(N)]
        for q in range(self.n):
            c = self.qmask[q]
            while c:
                g = (c & -c).bit_length() - 1
                self.check_qubits[g].append(q)
                c &= c - 1
        # c = max{wt(a), wt(b)}: most X-checks meeting one qubit (Sec. V-A)
        self.qdeg = max(len(self.a_terms), len(self.b_terms))

    # -- membership tests ---------------------------------------------------
    def split(self, z):
        return z & ((1 << self.N) - 1), z >> self.N

    def join(self, u, v):
        return u | (v << self.N)

    def in_K(self, z):
        u, v = self.split(z)
        return self.mul(self.a, u) ^ self.mul(self.b, v) == 0

    def is_logical(self, z):
        """z in K \\ S, i.e. a nontrivial Z-logical (eq. (dz))"""
        return z != 0 and self.in_K(z) and z not in self.S

    def pi(self, z):
        """pi[(u,v)] = v + (a) of Theorem 1 [thm:exact], as a reduced vector"""
        return self.ideal_a.reduce(self.split(z)[1])

    def component(self, z):
        """'A' (annihilator, pi = 0) or 'C' (colon, pi != 0), Sec. IV-A"""
        return "A" if self.pi(z) == 0 else "C"

    def class_key(self, z):
        """canonical representative of z + S, used to group the census"""
        return self.S.reduce(z)

    def translate(self, z, di, dj):
        """action of the translation group G, Proposition 1 [prop:equiv]"""
        u, v = self.split(z)
        return self.join(self.shift(u, di, dj), self.shift(v, di, dj))

    # -- Theorem 1 [thm:exact], Corollary 1 [cor:basis], Lemma 1 [lem:frob] --
    def decomposition(self):
        N = self.N
        ann_a = nullspace(self.Acols, N)                           # ker A
        b_ann_a = Span(self.mul(self.b, t) for t in ann_a)        # B(ker A)
        # colon ideal (a : b) = B^{-1}(im A), eq. (colon): kernel of v -> b v mod (a)
        colon = nullspace([self.ideal_a.reduce(c) for c in self.Bcols], N)
        # Lemma 1 [lem:frob], eq. (frobcolon): (a : b) = ann(b ann(a)).
        # b ann(a) is itself an ideal, so its vector-space basis generates it.
        ann_J = self._ann_of(b_ann_a.basis())
        colon_span = Span(colon)
        frob_ok = len(ann_J) == len(colon) and all(x in colon_span for x in ann_J)
        rA = len(ann_a) - len(b_ann_a)
        rC = len(colon) - len(self.ideal_a)
        return {
            "dim_ann_a": len(ann_a), "dim_b_ann_a": len(b_ann_a),
            "dim_colon": len(colon), "dim_ideal_a": len(self.ideal_a),
            "r_A": rA, "r_C": rC, "k": self.k,
            "sequence_ok": rA + rC == self.k,          # eq. (dimk)
            "frobenius_ok": frob_ok and rA == rC,      # eq. (balanced-dim)
            "ann_a": ann_a, "b_ann_a": b_ann_a, "colon": colon,
        }

    def _ann_of(self, gens):
        """ann(J) for the ideal J spanned, as a vector space, by gens"""
        N = self.N
        if not gens:
            return [1 << q for q in range(N)]
        images = []
        for q in range(N):
            e = 1 << q
            img = 0
            for s, g in enumerate(gens):
                img |= self.mul(e, g) << (s * N)
            images.append(img)
        return nullspace(images, N)


# ---------------------------------------------------------------------------
# Upper bounds: Algorithm "search for low-weight logical witnesses"
# [alg:screen], following the search order of Sec. IV-B.
# ---------------------------------------------------------------------------

def gray_min(basis, forbidden, max_dim=22, samples=200000, rng=None):
    """Least weight in span(basis) outside the subspace `forbidden`.

    With basis = ann(a) and forbidden = b ann(a) this is w_ann of
    eq. (wann). Exhaustive by Gray code when the dimension allows (this is how
    the w_ann column of Table VII [tab:profile] is obtained); otherwise random
    combinations, which only give an upper bound, and the flag says so.
    """
    d = len(basis)
    best, arg = None, None
    if d <= max_dim:
        cur = 0
        for g in range(1, 1 << d):
            cur ^= basis[(g & -g).bit_length() - 1]
            w = popcount(cur)
            if (best is None or w < best) and cur not in forbidden:
                best, arg = w, cur
        return best, arg, True
    rng = rng or random.Random(1)
    for _ in range(samples):
        cur = 0
        for bvec in basis:
            if rng.random() < 0.5:
                cur ^= bvec
        if cur and (best is None or popcount(cur) < best) and cur not in forbidden:
            best, arg = popcount(cur), cur
    return best, arg, False


def one_sided_witnesses(code, max_dim=22):
    """Steps 1 and 2 of Sec. IV-B: (t,0) with t in ann(a) \\ b ann(a), and
    (0,s) with s in ann(b) \\ a ann(b). The right-only ones are labelled by
    evaluating pi, since they may be colon classes (Example 4 [ex:72])."""
    N = code.N
    out = []
    for left in (True, False):
        g = code.b if left else code.a
        cols = code.Acols if left else code.Bcols
        ann = nullspace(cols, N)
        forb = Span(code.mul(g, t) for t in ann)
        w, t, exact = gray_min(ann, forb, max_dim)
        if w is not None:
            z = code.join(t, 0) if left else code.join(0, t)
            assert code.is_logical(z)
            out.append((w, z, "one-sided " + ("left" if left else "right"), exact))
    return out


def coset_min(u0, ann, exhaustive_dim=10):
    """Lightest element of u0 + span(ann), i.e. lambda_a of eq. (lambda).

    Exact by Gray code for small ann; otherwise a greedy descent, which is an
    upper bound on lambda_a and is all a witness needs.
    """
    if len(ann) <= exhaustive_dim:
        best = cur = u0
        for gi in range(1, 1 << len(ann)):
            cur ^= ann[(gi & -gi).bit_length() - 1]
            if popcount(cur) < popcount(best):
                best = cur
        return best
    best, improved = u0, True
    while improved:
        improved = False
        for t in ann:
            if popcount(best ^ t) < popcount(best):
                best ^= t
                improved = True
    return best


def colon_lift_witnesses(code, vmax=3, max_ann_dim=10):
    """Step 3 of Sec. IV-B: sparse two-block relations a u = b v.

    v runs over sparse words with one monomial anchored at 1, which is enough
    by translation equivariance (Proposition 1 [prop:equiv]). For each v with
    b v in (a) we solve a u = b v and minimise wt(u) over u0 + ann(a), the
    coset-leader weight of eq. (dcol). The same is repeated with a and b
    exchanged, as in the last loop of [alg:screen].
    """
    N = code.N
    found = []
    for swap in (False, True):
        g = code.a if swap else code.b
        fcols = code.Bcols if swap else code.Acols
        ann = nullspace(fcols, N)
        solver = Span((c << N) | (1 << i) for i, c in enumerate(fcols))
        maskN = (1 << N) - 1
        for wv in range(1, vmax + 1):
            for rest in combinations(range(1, N), wv - 1):
                v = 1
                for q in rest:
                    v |= 1 << q
                w = code.mul(g, v)
                red = solver.reduce(w << N)
                if red >> N:                  # g v not in (f): v not in the colon ideal
                    continue
                u0 = red & maskN              # f u0 = g v
                best_u = coset_min(u0, ann, max_ann_dim)
                z = code.join(v, best_u) if swap else code.join(best_u, v)
                if code.is_logical(z):
                    found.append((popcount(z), z, "colon lift" + (" (b-side)" if swap else ""), False))
    found.sort(key=lambda t: t[0])
    return found[:5]


# ---------------------------------------------------------------------------
# Lower bound, baseline: Algorithm "orbit-reduced support exclusion"
# [alg:supp] and Theorem 3 [thm:supp]. Used in Sec. V-D for [[18,4,4]] and
# [[72,12,6]]; its cost grows like binom(n, d0-1), so it stops there.
# ---------------------------------------------------------------------------

def support_exclusion(code, d0):
    """True if K(E) = S(E) for every anchored support |E| < d0, i.e. d_Z >= d0.

    Anchoring (Sec. V, paragraph after Theorem 3): one representative per
    translation orbit, supports through the left origin plus right-only
    supports through the right origin. K(E) is the kernel of the restricted
    check columns; K(E) = S(E) iff every basis vector of K(E) lies in S.
    Returns (ok, offending_support_or_None, supports_tested).
    """
    N, n = code.N, code.n
    tested = 0
    for anchor, pool in ((0, range(1, n)), (N, range(N + 1, n))):
        for size in range(1, d0):
            for rest in combinations(pool, size - 1):
                E = (anchor,) + rest
                tested += 1
                for kvec in nullspace([code.qmask[q] for q in E], len(E)):
                    z = 0
                    for i, q in enumerate(E):
                        if kvec >> i & 1:
                            z |= 1 << q
                    if z not in code.S:
                        return False, E, tested
    return True, None, tested


# ---------------------------------------------------------------------------
# Lower bound, main method: Algorithm "anchored cluster search"
# [alg:cluster], Sec. V-A.
# ---------------------------------------------------------------------------

class Budget(Exception):
    pass


def cluster_search_py(code, W, mode="logical", enumerate_all=False, node_budget=None):
    """Pure-Python version of [alg:cluster].

    mode = "logical": a leaf counts if not in S. No hit certifies d >= W+1
                      (Theorem 4 [thm:cluster]).
    mode = "colon"  : a leaf counts if pi != 0. No hit certifies d_C >= W+1
                      (Corollary 3 [cor:colonconn]).
    mode = "kernel" : any nonzero element of K counts; gives w_K for
                      Proposition 4 [prop:irred].
    enumerate_all   : do not stop at the first hit (census, Table V
                      [tab:census]; complete at W = d, and for W < 2 w_K by
                      Proposition 4 [prop:irred]).
    Returns (status, hits, nodes), status in {"found","none","enumerated","budget"}.
    Nodes are calls of GROW over both anchors, as counted in Table IV
    [tab:clustercert].
    """
    N, n = code.N, code.n
    qmask, cq, c = code.qmask, code.check_qubits, code.qdeg
    hits, state = [], {"nodes": 0}
    inT = [False] * n
    T = []

    def leaf_counts(z):
        if mode == "kernel":
            return True
        if mode == "colon":
            return code.pi(z) != 0
        return z not in code.S

    def grow(syn, z, rightonly):
        state["nodes"] += 1
        if node_budget is not None and state["nodes"] > node_budget:
            raise Budget
        if syn == 0:
            # zero syndrome: output or prune. Pruning is safe because a
            # minimum-weight target has no zero-syndrome proper subset
            # (Lemma 2 [lem:conn]; colon and kernel analogues as above).
            if leaf_counts(z):
                hits.append(z)
                if not enumerate_all:
                    return True
            return False
        u = popcount(syn)
        if len(T) + (u + c - 1) // c > W:          # parity prune |T| + ceil(u/c) > W
            return False
        g = (syn & -syn).bit_length() - 1          # lowest-index unsatisfied check
        for q in cq[g]:                            # at most wt(a)+wt(b)-1 branches
            if inT[q] or (rightonly and q < N):
                continue
            inT[q] = True
            T.append(q)
            done = grow(syn ^ qmask[q], z | (1 << q), rightonly)
            T.pop()
            inT[q] = False
            if done:
                return True
        return False

    sys.setrecursionlimit(max(1000, 4 * W + 100))
    try:
        # two anchors, as in the proof of Theorem 4 [thm:cluster]
        for anchor, rightonly in ((0, False), (N, True)):
            inT[anchor] = True
            T.append(anchor)
            done = grow(qmask[anchor], 1 << anchor, rightonly)
            T.pop()
            inT[anchor] = False
            if done:
                break
    except Budget:
        return "budget", hits, state["nodes"]
    if enumerate_all:
        return "enumerated", hits, state["nodes"]
    return ("found" if hits else "none"), hits, state["nodes"]


def find_c_engine():
    here = os.path.dirname(os.path.abspath(__file__))
    for cand in (os.path.join(here, "cluster_lb"), shutil.which("cluster_lb") or ""):
        if cand and os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    return None


def cluster_search_c(code, W, mode="logical", enumerate_all=False, node_budget=None, exe=None):
    """Same search through the C engine; parses its LOGICAL / RESULT lines."""
    exe = exe or find_c_engine()
    args = [exe]
    if enumerate_all:
        args.append("-e")
    if mode == "colon":
        args.append("-c")
    if mode == "kernel":
        args.append("-k")
    if node_budget is not None:
        args += ["-n", str(int(node_budget))]
    fmt = lambda ts: "+".join(("x^%d" % i if i else "") + ("*" if i and j else "") +
                              ("y^%d" % j if j else "") or "1" for i, j in ts)
    args += [str(code.l), str(code.m), str(W), fmt(code.a_terms), fmt(code.b_terms)]
    out = subprocess.run(args, capture_output=True, text=True, check=True).stdout
    hits, status, nodes = [], None, 0
    for line in out.splitlines():
        if line.startswith("LOGICAL"):
            qs = list(map(int, line.split()[2:]))
            hits.append(sum(1 << q for q in qs))
        elif line.startswith("RESULT"):
            status = line.split()[1]
            m_ = re.search(r"nodes=(\d+)", line)
            nodes = int(m_.group(1)) if m_ else 0
    return status, hits, nodes


def cluster_search(code, W, mode="logical", enumerate_all=False, node_budget=None, engine="auto"):
    exe = find_c_engine() if engine in ("auto", "c") else None
    if engine == "c" and not exe:
        raise RuntimeError("C engine requested but ./cluster_lb not found (see README)")
    if exe:
        return cluster_search_c(code, W, mode, enumerate_all, node_budget, exe)
    return cluster_search_py(code, W, mode, enumerate_all, node_budget)


# ---------------------------------------------------------------------------
# Screening on random pairs, Sec. VI-B and Table VIII [tab:sweep]
# ---------------------------------------------------------------------------

def sweep_sample(l, m, seed, samples=200):
    """The sampler of Sec. VI-B: independent uniform weight-three a and b
    (distinct monomials within each, overlaps between a and b allowed), drawn
    with Python's random.Random(seed); seed 3 for (6,6), seed 5 for (9,6)."""
    mons = [(i, j) for i in range(l) for j in range(m)]
    rng = random.Random(seed)
    return [(rng.sample(mons, 3), rng.sample(mons, 3)) for _ in range(samples)]


def sweep_screens(code):
    """The two diagnostics reported in Table VIII [tab:sweep].

    w_one: exact w_ann of eq. (wann) on the a-side (one-sided annihilator screen).
    lop  : a sparse colon shortcut v in (a:b) \\ (a) with wt(v) <= 2.
    """
    ann = nullspace(code.Acols, code.N)
    w_one = gray_min(ann, Span(code.mul(code.b, t) for t in ann), max_dim=24)[0]
    lop = False
    for w in (1, 2):
        for S in combinations(range(code.N), w):
            v, bv = 0, 0
            for q in S:
                v |= 1 << q
                bv ^= code.Bcols[q]
            if code.ideal_a.reduce(bv) == 0 and code.ideal_a.reduce(v) != 0:
                lop = True
                break
        if lop:
            break
    return w_one, lop


# ---------------------------------------------------------------------------
# Putting it together: the procedure of Sec. VI-C [sec:workflow]
# ---------------------------------------------------------------------------

def distance_bounds(code, seconds=60.0, engine="auto", verbose=True):
    """Return (d_L, d_U, witness, notes).

    d_U: lightest explicit logical among the screens of [alg:screen].
    d_L: the cluster search at W = 1, 2, ...; a radius with no hit certifies
         d >= W+1 (Theorem 4 [thm:cluster]). The first radius with a hit gives d
         exactly, since all smaller radii were already excluded (Sec. V-D). We
         stop when the bounds meet or the time budget is spent. d_X = d_Z by
         Proposition 3 [prop:sym], so this is the code distance.
    """
    log = print if verbose else (lambda *a, **k: None)
    if code.k == 0:
        return None, None, None, ["k = 0, no logical qubits"]
    notes = []
    cands = one_sided_witnesses(code) + colon_lift_witnesses(code, vmax=3 if code.N > 60 else 4)
    cands.sort(key=lambda t: t[0])
    dU, wit, how, _ = cands[0] if cands else (code.n, None, "none", False)
    log("  upper bound from screens: %d (%s)" % (dU, how))

    t0, dL = time.time(), 1
    W = 1
    while dL < dU:
        left = seconds - (time.time() - t0)
        if left <= 0:
            notes.append("time budget spent at radius %d" % W)
            break
        # rough node budget from the remaining time: the C engine manages a few
        # times 1e7 nodes/s, Python roughly 2e6/s. Only used to stop runaways.
        rate = 3e7 if (engine != "py" and find_c_engine()) else 1.5e6
        status, hits, nodes = cluster_search(code, W, node_budget=int(rate * left), engine=engine)
        log("  cluster search radius %2d: %-5s nodes=%d" % (W, status, nodes))
        if status == "budget":
            notes.append("node budget hit at radius %d" % W)
            break
        if status == "found":
            z = min(hits, key=popcount)
            dL = dU = popcount(z)
            wit, how = z, "cluster search"
            break
        dL = W + 1
        W += 1
    return dL, dU, wit, notes


def describe_witness(code, z):
    u, v = code.split(z)
    return "u = %s ; v = %s ; split (%d,%d) ; component %s" % (
        code.fmt(u), code.fmt(v), popcount(u), popcount(v), code.component(z))


def main():
    ap = argparse.ArgumentParser(description="Distance bounds for a bivariate bicycle code.")
    ap.add_argument("-l", type=int, required=True)
    ap.add_argument("-m", type=int, required=True)
    ap.add_argument("-a", required=True, help='e.g. "x^3+y+y^2"')
    ap.add_argument("-b", required=True, help='e.g. "y^3+x+x^2"')
    ap.add_argument("--seconds", type=float, default=60.0, help="time budget for the lower bound")
    ap.add_argument("--engine", choices=["auto", "c", "py"], default="auto")
    args = ap.parse_args()

    code = BBCode(args.l, args.m, args.a, args.b)
    print("BB code on %dx%d torus: a = %s, b = %s" % (code.l, code.m, code.fmt(code.a), code.fmt(code.b)))
    print("n = %d, k = %d" % (code.n, code.k))
    if code.k == 0:
        print("no logical qubits, nothing to bound")
        return
    dec = code.decomposition()
    print("exact sequence (Theorem 1): dim ann(a)=%d, dim b ann(a)=%d, r_A=%d; dim(a:b)=%d, dim(a)=%d, r_C=%d"
          % (dec["dim_ann_a"], dec["dim_b_ann_a"], dec["r_A"], dec["dim_colon"], dec["dim_ideal_a"], dec["r_C"]))
    print("checks: r_A + r_C = k (Corollary 1) -> %s ; (a:b) = ann(b ann(a)) and r_A = r_C (Lemma 1) -> %s"
          % (dec["sequence_ok"], dec["frobenius_ok"]))
    dL, dU, wit, notes = distance_bounds(code, args.seconds, args.engine)
    print()
    if dL == dU:
        print("distance established: d = %d  (d_X = d_Z by Proposition 3)" % dL)
    else:
        print("bounds: %d <= d <= %d" % (dL, dU))
    for nt in notes:
        print("  note:", nt)
    if wit is not None:
        print("witness:", describe_witness(code, wit))


if __name__ == "__main__":
    main()
