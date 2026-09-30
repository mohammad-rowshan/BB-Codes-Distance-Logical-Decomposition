#!/usr/bin/env python3
"""
reproduce.py -- regenerate the numbers reported in
"Logical Operator Decomposition for Distance Analysis of Bivariate Bicycle Codes".

Manuscript items are cited by number and LaTeX label, e.g. Table V [tab:census].

  [1] Table VII [tab:profile], structural columns: k, dim ann(a), dim b ann(a),
      r_A = r_C, the exact-sequence count of Corollary 1 [cor:basis] and the
      identity (a:b) = ann(b ann(a)) of Lemma 1 [lem:frob]; exhaustive w_ann,
      eq. (wann).
  [2] Table III [tab:wit] and Table IX [tab:codedata]: each witness lies in
      K \\ S with weight d; its split and component. Example 4 [ex:72] (both
      one-sided witnesses) and Example 5 [ex:18] ((a,1) in C, (1,b) in A,
      w_ann = 6 > d_A = 4).
  [3] Table IV [tab:clustercert]: the anchored cluster search [alg:cluster] at
      radius d-1 finds nothing, so d >= d_U (Theorem 4 [thm:cluster]); node counts.
      For [[18,4,4]] and [[72,12,6]] also the support-exclusion test of
      Theorem 3 [thm:supp] / [alg:supp], as stated in Sec. V-D.
  [4] Table V [tab:census]: all minimum-weight logicals grouped by class,
      labelled by component and by shape (Definition 1 [def:shape], tau = 3);
      d_A and d_C of Table VII where the census attains them.
  [5] Example 6 [ex:108]: w_K = 6, no logical of weight 11 (Proposition 4
      [prop:irred]), d_C = 10 (Corollary 3 [cor:colonconn]), d_A = 12.
  [6] Table VIII [tab:sweep]: the random-pair screening of Sec. VI-B.
  [7] optional (--milp): the integer-program cross-check of Sec. V-C,
      eq. (milpobj), for n <= 108. Needs scipy.

Usage:
    gcc -O3 -march=native -o cluster_lb cluster_lb.c     # recommended
    python3 reproduce.py                 # everything except the MILP
    python3 reproduce.py --skip-288      # a few seconds
    python3 reproduce.py --milp          # add the MILP cross-check

Without the C engine the pure-Python search is used (about 2e6 nodes/s); the
[[288,12,18]] certificate then takes a few minutes and its census is skipped
unless --slow is given.
"""

import argparse
import sys
import time
from collections import Counter, defaultdict

from bbcode import (BBCode, Span, nullspace, gray_min, cluster_search, find_c_engine,
                    popcount, support_exclusion, sweep_sample, sweep_screens)

# The six standard codes (Table IX [tab:codedata]) with the values printed in
# the manuscript. dA, dC are the exact component minima of Table VII.
CODES = [
    dict(name="[[18,4,4]]",    l=3,  m=3,  a="x+y+xy",      b="x^2+y^2+x^2y^2",
         k=4,  d=4,  dimann=2,  dimbann=0,  rA=2, w_ann=6,  dA=4,  dC=4,  nodes=19),
    dict(name="[[72,12,6]]",   l=6,  m=6,  a="x^3+y+y^2",   b="y^3+x+x^2",
         k=12, d=6,  dimann=12, dimbann=6,  rA=6, w_ann=6,  dA=6,  dC=6,  nodes=171),
    dict(name="[[90,8,10]]",   l=15, m=3,  a="x^9+y+y^2",   b="1+x^2+x^7",
         k=8,  d=10, dimann=6,  dimbann=2,  rA=4, w_ann=10, dA=10, dC=10, nodes=29301),
    dict(name="[[108,8,10]]",  l=9,  m=6,  a="x^3+y+y^2",   b="y^3+x+x^2",
         k=8,  d=10, dimann=6,  dimbann=2,  rA=4, w_ann=12, dA=12, dC=10, nodes=24863),
    dict(name="[[144,12,12]]", l=12, m=6,  a="x^3+y+y^2",   b="y^3+x+x^2",
         k=12, d=12, dimann=12, dimbann=6,  rA=6, w_ann=12, dA=12, dC=12, nodes=284281),
    dict(name="[[288,12,18]]", l=12, m=12, a="x^3+y^2+y^7", b="y^3+x+x^2",
         k=12, d=18, dimann=24, dimbann=18, rA=6, w_ann=18, dA=18, dC=18, nodes=455504708),
]

# Table V [tab:census]: classes attaining d, then per component the counts of
# one-sided / lopsided / balanced classes.
CENSUS = {
    "[[18,4,4]]":    (15,  (0, 3, 0),  (0, 12, 0)),
    "[[72,12,6]]":   (84,  (36, 0, 0), (36, 12, 0)),
    "[[90,8,10]]":   (72,  (9, 0, 0),  (9, 18, 36)),
    "[[108,8,10]]":  (9,   (0, 0, 0),  (0, 0, 9)),
    "[[144,12,12]]": (246, (36, 0, 12), (0, 72, 126)),
    "[[288,12,18]]": (84,  (36, 0, 0), (36, 0, 12)),
}

# Table VIII [tab:sweep]: (torus, seed) -> (k>0 count, k histogram,
# percent with w_ann <= 6, percent with w_ann <= 12, percent lopsided)
SWEEP = {
    (6, 6, 3): (42, {4: 28, 8: 13, 24: 1}, 45, 71, 57),
    (9, 6, 5): (39, {4: 36, 8: 3}, 3, 28, 74),
}

FAILS = []


def check(label, got, want):
    ok = got == want
    if not ok:
        FAILS.append(label)
    print("    %-46s %-18s %s" % (label, str(got), "ok" if ok else "MISMATCH (manuscript: %s)" % (want,)))


def witnesses(code, name):
    """Explicit logicals of Table III [tab:wit] / Table IX [tab:codedata]."""
    P = code.poly
    a, b = code.a, code.b
    if name == "[[18,4,4]]":                    # Example 5 [ex:18]
        return [("(a, 1)", code.join(a, 1)), ("(1, b)", code.join(1, b))]
    if name == "[[72,12,6]]":                   # Example 4 [ex:72], both witnesses
        return [("((1+y^2)a, 0)", code.join(code.mul(P([(0, 0), (0, 2)]), a), 0)),
                ("(0, (1+x^2)b)", code.join(0, code.mul(P([(0, 0), (2, 0)]), b)))]
    if name == "[[90,8,10]]":
        u = code.mul(P([(0, 0), (0, 1)]), P([(3 * j, 0) for j in range(5)]))
        return [("((1+y) sum x^{3j}, 0)", code.join(u, 0))]
    if name == "[[108,8,10]]":                  # eq. (w108)
        u = P([(1, 0), (1, 1), (4, 4), (5, 0), (7, 3), (8, 4)])
        v = P([(1, 4), (3, 1), (4, 2), (6, 5)])
        return [("eq. (w108)", code.join(u, v))]
    if name == "[[144,12,12]]":
        u = code.mul(code.mul(P([(0, 0), (6, 0)]), P([(0, 0), (0, 2)])), a)
        return [("((1+x^6)(1+y^2)a, 0)", code.join(u, 0))]
    if name == "[[288,12,18]]":                 # eq. (w288)
        u = P([(0, 1), (0, 3), (0, 4), (0, 6), (0, 8), (0, 9), (0, 10), (0, 11),
               (3, 2), (3, 4), (3, 6), (3, 8), (6, 2), (6, 6), (6, 7), (6, 11), (9, 0), (9, 4)])
        return [("eq. (w288)", code.join(u, 0))]
    return []


def census(code, d, engine):
    """All weight-d logicals -> classes -> (component, shape) tally.

    Enumeration mode of [alg:cluster] lists every minimum-weight logical through
    an anchor (paragraph after Theorem 4 [thm:cluster]); we close the list under
    the l*m translations and group by coset of S. Shape follows Definition 1
    [def:shape] with tau = max(wt a, wt b) = 3.
    """
    status, hits, nodes = cluster_search(code, d, enumerate_all=True, engine=engine)
    vecs = set()
    for z in hits:
        if popcount(z) != d:
            continue
        for i in range(code.l):
            for j in range(code.m):
                vecs.add(code.translate(z, i, j))
    classes = defaultdict(list)
    for z in vecs:
        classes[code.class_key(z)].append(z)
    tau = max(len(code.a_terms), len(code.b_terms))
    tally = {"A": [0, 0, 0], "C": [0, 0, 0]}
    for zs in classes.values():
        comp = code.component(zs[0])          # constant on a class (Theorem 1)
        sides = [min(popcount(code.split(z)[0]), popcount(code.split(z)[1])) for z in zs]
        if 0 in sides:
            shape = 0                          # one-sided
        elif min(sides) <= tau:
            shape = 1                          # tau-lopsided
        else:
            shape = 2                          # tau-balanced
        tally[comp][shape] += 1
    return len(vecs), len(classes), tuple(tally["A"]), tuple(tally["C"]), nodes


def milp_distance(code, time_limit=600):
    """eq. (milpobj): min 1^T z s.t. H_X z = 2p, L z - y = 2q, sum y >= 1.

    L has rows vanishing on S but independent on K, so ker(L|_K) = S as
    required in Sec. V-C. Floating point (HiGHS), so a cross-check only.
    """
    import numpy as np
    from scipy.optimize import milp, LinearConstraint, Bounds
    n, N = code.n, code.N
    bits = lambda v: [(v >> i) & 1 for i in range(n)]
    HX = [[(code.qmask[q] >> g) & 1 for q in range(n)] for g in range(N)]
    Sperp = nullspace([sum(((s >> q) & 1) << r for r, s in enumerate(code.S.basis()))
                       for q in range(n)], n)
    Kperp = Span(nullspace([sum(((kv >> q) & 1) << r for r, kv in enumerate(code.K))
                            for q in range(n)], n))
    L = [f for f in Sperp if Kperp.add(f)]
    k = len(L)
    nv = n + k + N + k
    rows, lo, hi = [], [], []
    for g in range(N):
        r = np.zeros(nv); r[:n] = HX[g]; r[n + k + g] = -2
        rows.append(r); lo.append(0); hi.append(0)
    for i, f in enumerate(L):
        r = np.zeros(nv); r[:n] = bits(f); r[n + i] = -1; r[n + k + N + i] = -2
        rows.append(r); lo.append(0); hi.append(0)
    r = np.zeros(nv); r[n:n + k] = 1
    rows.append(r); lo.append(1); hi.append(np.inf)
    cost = np.zeros(nv); cost[:n] = 1
    ub = np.concatenate([np.ones(n + k), np.full(N + k, n)])
    res = milp(cost, constraints=LinearConstraint(np.array(rows), lo, hi),
               integrality=np.ones(nv), bounds=Bounds(0, ub),
               options={"time_limit": time_limit, "mip_rel_gap": 0})
    return (round(res.fun) if res.status == 0 else None), res.message


def run_sweep():
    """[6] Table VIII [tab:sweep], Sec. VI-B."""
    print("\n=== Table VIII: screening on random pairs (Sec. VI-B) ===")
    pct = lambda x, n: int(round(100.0 * x / n))
    for (l, m, seed), want in SWEEP.items():
        rows = []
        for a, b in sweep_sample(l, m, seed):
            code = BBCode(l, m, a, b)
            if code.k > 0:
                rows.append((code.k,) + sweep_screens(code))
        n = len(rows)
        hist = dict(sorted(Counter(r[0] for r in rows).items()))
        le6 = sum(r[1] is not None and r[1] <= 6 for r in rows)
        le12 = sum(r[1] is not None and r[1] <= 12 for r in rows)
        lop = sum(r[2] for r in rows)
        print("  torus (%d,%d), seed %d, 200 samples" % (l, m, seed))
        check("pairs with k > 0", n, want[0])
        check("k histogram", hist, want[1])
        check("%% with w_ann <= 6   (%d of %d)" % (le6, n), pct(le6, n), want[2])
        check("%% with w_ann <= 12  (%d of %d)" % (le12, n), pct(le12, n), want[3])
        check("%% with sparse colon shortcut (%d of %d)" % (lop, n), pct(lop, n), want[4])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-288", action="store_true", help="leave out [[288,12,18]]")
    ap.add_argument("--engine", choices=["auto", "c", "py"], default="auto")
    ap.add_argument("--milp", action="store_true", help="MILP cross-check of Sec. V-C for n <= 108")
    ap.add_argument("--slow", action="store_true",
                    help="with the Python engine, also run the [[288]] census (20+ min)")
    ap.add_argument("--no-sweep", action="store_true", help="skip Table VIII")
    args = ap.parse_args()

    have_c = find_c_engine() is not None
    py_only = (not have_c) or args.engine == "py"
    print("cluster engine:", "pure Python" if py_only else "C (./cluster_lb)")
    codes = [c for c in CODES if not (args.skip_288 and c["name"] == "[[288,12,18]]")]
    if py_only and not args.skip_288:
        print("  (pure Python: the [[288]] certificate takes a few minutes; its census\n"
              "   is skipped unless --slow. Compile cluster_lb.c to make both quick.)")

    for spec in codes:
        t_code = time.time()
        code = BBCode(spec["l"], spec["m"], spec["a"], spec["b"], spec["name"])
        name, d = spec["name"], spec["d"]
        big = name == "[[288,12,18]]"
        print("\n=== %s   a = %s, b = %s on %dx%d ===" % (name, spec["a"], spec["b"], spec["l"], spec["m"]))

        # [1] Table VII, structural columns ---------------------------------
        print("  [1] Table VII: exact-sequence profile")
        dec = code.decomposition()
        check("k", code.k, spec["k"])
        check("dim ann(a)", dec["dim_ann_a"], spec["dimann"])
        check("dim b ann(a)", dec["dim_b_ann_a"], spec["dimbann"])
        check("(r_A, r_C), Lemma 1", (dec["r_A"], dec["r_C"]), (spec["rA"], spec["rA"]))
        check("eq. (dimk) and eq. (frobcolon) hold", dec["sequence_ok"] and dec["frobenius_ok"], True)
        # exhaustive over ann(a); for [[288]] that is 2^24 elements (Sec. VI-A)
        w, _, exact = gray_min(dec["ann_a"], dec["b_ann_a"], max_dim=24)
        check("w_ann, eq. (wann) (exhaustive=%s)" % exact, w, spec["w_ann"])

        # [2] Table III / IX witnesses ----------------------------------------
        print("  [2] Table III: explicit minimum-weight witnesses")
        for label, z in witnesses(code, name):
            u, v = code.split(z)
            check("%s in K\\S, weight" % label, (code.is_logical(z), popcount(z)), (True, d))
            print("        split (%d,%d), component %s" % (popcount(u), popcount(v), code.component(z)))
        if name == "[[18,4,4]]":
            # Example 5 [ex:18]: (1,b) is two-sided yet pi = 0
            check("Example 5: (1,b) lies in A", code.component(code.join(1, code.b)), "A")
        if name == "[[72,12,6]]":
            # Example 4 [ex:72]: the right-only witness is a colon class
            s = code.mul(code.poly([(0, 0), (2, 0)]), code.b)
            check("Example 4: (0,(1+x^2)b) lies in C", code.component(code.join(0, s)), "C")

        # [3] Table IV certificate, plus Sec. V-D support exclusion ----------
        print("  [3] Table IV: anchored cluster search at radius d-1")
        t0 = time.time()
        status, hits, nodes = cluster_search(code, d - 1, engine=args.engine)
        check("radius %d: no logical  => d >= %d" % (d - 1, d), status, "none")
        check("nodes", nodes, spec["nodes"])
        print("        (%.1f s on this machine)" % (time.time() - t0))
        if name in ("[[18,4,4]]", "[[72,12,6]]"):
            t0 = time.time()
            ok, E, tested = support_exclusion(code, d)
            check("Theorem 3, all |E| <= %d: K(E) = S(E)" % (d - 1), ok, True)
            print("        %d anchored supports, %.1f s" % (tested, time.time() - t0))

        # [4] Table V census -----------------------------------------------
        print("  [4] Table V: census of minimum-weight classes")
        if big and py_only and not args.slow:
            print("      skipped in pure Python (use --slow, or the C engine)")
        else:
            nvec, ncls, tA, tC, _ = census(code, d, args.engine)
            want = CENSUS[name]
            check("classes attaining d", ncls, want[0])
            check("A: one-sided/lopsided/balanced", tA, want[1])
            check("C: one-sided/lopsided/balanced", tC, want[2])
            print("        (%d minimum-weight vectors)" % nvec)
            # a component present in the census attains d (Table VII)
            if sum(tA):
                check("d_A = d (Table VII)", d, spec["dA"])
            if sum(tC):
                check("d_C = d (Table VII)", d, spec["dC"])

        # [5] Example 6 -------------------------------------------------------
        if name == "[[108,8,10]]":
            print("  [5] Example 6: component minima of [[108,8,10]]")
            wK = None
            for W in range(2, 8):              # kernel mode: minimum nonzero weight of K
                st, _, _ = cluster_search(code, W, mode="kernel", engine=args.engine)
                if st == "found":
                    wK = W
                    break
            check("w_K", wK, 6)
            # Proposition 4 [prop:irred]: radius 2 w_K - 1 = 11 lists every logical up to weight 11
            st, hits, _ = cluster_search(code, 2 * wK - 1, enumerate_all=True, engine=args.engine)
            weights = sorted({popcount(z) for z in hits})
            compA = [z for z in hits if code.component(z) == "A"]
            translates = {code.translate(z, i, j) for z in hits for i in range(code.l) for j in range(code.m)}
            check("weights of all logicals of weight <= 11", weights, [10])
            check("number of weight-10 logicals", len(translates), 54)
            check("annihilator logicals of weight <= 11", len(compA), 0)
            check("hence d_A = w_ann", 12 if not compA else None, spec["dA"])
            st, _, _ = cluster_search(code, 9, mode="colon", engine=args.engine)
            check("Corollary 3: colon search, radius 9", st, "none")

        # [7] MILP cross-check -------------------------------------------------
        if args.milp and code.n <= 108:
            try:
                val, _ = milp_distance(code)
                check("MILP optimum, eq. (milpobj) (cross-check)", val, d)
            except ImportError:
                print("      scipy not available, MILP skipped")
        print("  (%.1f s for this code)" % (time.time() - t_code))

    if not args.no_sweep:
        run_sweep()

    print("\n%s" % ("all checks passed" if not FAILS else "MISMATCHES: " + ", ".join(FAILS)))
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
