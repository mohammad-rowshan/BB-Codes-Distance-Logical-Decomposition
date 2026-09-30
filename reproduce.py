#!/usr/bin/env python3
"""
reproduce.py -- regenerate the main numbers of the manuscript
"Logical Operator Decomposition for Distance Analysis of Bivariate Bicycle Codes".

What is reproduced (section / table names follow the revised manuscript):

  [1] Table "profile": dim ann(a), dim b ann(a), r_A, r_C, k, w_ann, and the
      two structural checks (Corollary "basis": r_A + r_C = k; Lemma
      "Frobenius balance": (a:b) = ann(b ann(a)), r_A = r_C).
  [2] Table "witnesses" and Appendix B: every listed logical is in K \\ S with
      the stated weight, split and component; Example "[[18,4,4]]: component
      is not shape" ((a,1) in C, (1,b) in A, w_ann = 6 > d_A = 4).
  [3] Table "cluster certificates": the anchored cluster search at radius
      d-1 finds nothing, so d >= d_U; with the witnesses, d is exact.
  [4] Table "census": all minimum-weight logicals, grouped into classes and
      labelled by component (pi) and shape (Definition "representative shape").
  [5] Example "[[108,8,10]]": w_K = 6, no logical of weight 11, d_C = 10,
      d_A = 12 (Proposition "irreducible range", Corollary "colon connectivity").
  [6] optional: the parity-constrained MILP cross-check (needs scipy), eq. (milp).

Usage:
    gcc -O3 -march=native -o cluster_lb cluster_lb.c     # strongly recomended
    python3 reproduce.py                 # everything except the MILP
    python3 reproduce.py --skip-288      # quick run, a few seconds
    python3 reproduce.py --milp          # add the MILP cross-check (slow for n >= 90)

Without the C engine the pure-Python search is used. It does roughly 2e6
nodes/s, so everything is quick except [[288,12,18]]: its certificate needs
~4.6e8 nodes (a few minutes) and its census ~2.5e9 (skipped unless --slow).
"""

import argparse
import sys
import time
from collections import defaultdict

from bbcode import (BBCode, Span, nullspace, gray_min, cluster_search,
                    find_c_engine, popcount)

# The six standard codes (Table "code data", Appendix B) with the values the
# manuscript reports. d_A, d_C are the exact component minima of Table "profile".
CODES = [
    dict(name="[[18,4,4]]",    l=3,  m=3,  a="x+y+xy",      b="x^2+y^2+x^2y^2",
         k=4,  d=4,  dimann=2,  dimbann=0,  rA=2, w_ann=6,  dA=4,  dC=4),
    dict(name="[[72,12,6]]",   l=6,  m=6,  a="x^3+y+y^2",   b="y^3+x+x^2",
         k=12, d=6,  dimann=12, dimbann=6,  rA=6, w_ann=6,  dA=6,  dC=6),
    dict(name="[[90,8,10]]",   l=15, m=3,  a="x^9+y+y^2",   b="1+x^2+x^7",
         k=8,  d=10, dimann=6,  dimbann=2,  rA=4, w_ann=10, dA=10, dC=10),
    dict(name="[[108,8,10]]",  l=9,  m=6,  a="x^3+y+y^2",   b="y^3+x+x^2",
         k=8,  d=10, dimann=6,  dimbann=2,  rA=4, w_ann=12, dA=12, dC=10),
    dict(name="[[144,12,12]]", l=12, m=6,  a="x^3+y+y^2",   b="y^3+x+x^2",
         k=12, d=12, dimann=12, dimbann=6,  rA=6, w_ann=12, dA=12, dC=12),
    dict(name="[[288,12,18]]", l=12, m=12, a="x^3+y^2+y^7", b="y^3+x+x^2",
         k=12, d=18, dimann=24, dimbann=18, rA=6, w_ann=18, dA=18, dC=18),
]

# Table "census": number of minimum-weight classes and, per component,
# the counts of one-sided / lopsided / balanced classes (tau = 3).
CENSUS = {
    "[[18,4,4]]":    (15,  (0, 3, 0),  (0, 12, 0)),
    "[[72,12,6]]":   (84,  (36, 0, 0), (36, 12, 0)),
    "[[90,8,10]]":   (72,  (9, 0, 0),  (9, 18, 36)),
    "[[108,8,10]]":  (9,   (0, 0, 0),  (0, 0, 9)),
    "[[144,12,12]]": (246, (36, 0, 12), (0, 72, 126)),
    "[[288,12,18]]": (84,  (36, 0, 0), (36, 0, 12)),
}

FAILS = []


def check(label, got, want):
    ok = got == want
    if not ok:
        FAILS.append(label)
    print("    %-44s %-18s %s" % (label, str(got), "ok" if ok else "MISMATCH (manuscript: %s)" % (want,)))


def witnesses(code, name):
    """The explicit logicals of Table "witnesses" / Appendix B, as vectors."""
    P = code.poly
    a, b = code.a, code.b
    if name == "[[18,4,4]]":
        return [("(a, 1)", code.join(a, 1)), ("(1, b)", code.join(1, b))]
    if name in ("[[72,12,6]]",):
        return [("((1+y^2)a, 0)", code.join(code.mul(P([(0, 0), (0, 2)]), a), 0))]
    if name == "[[90,8,10]]":
        u = code.mul(P([(0, 0), (0, 1)]), P([(3 * j, 0) for j in range(5)]))
        return [("((1+y) sum x^{3j}, 0)", code.join(u, 0))]
    if name == "[[108,8,10]]":   # eq. (w108)
        u = P([(1, 0), (1, 1), (4, 4), (5, 0), (7, 3), (8, 4)])
        v = P([(1, 4), (3, 1), (4, 2), (6, 5)])
        return [("eq. (w108)", code.join(u, v))]
    if name == "[[144,12,12]]":
        u = code.mul(code.mul(P([(0, 0), (6, 0)]), P([(0, 0), (0, 2)])), a)
        return [("((1+x^6)(1+y^2)a, 0)", code.join(u, 0))]
    if name == "[[288,12,18]]":  # eq. (w288)
        u = P([(0, 1), (0, 3), (0, 4), (0, 6), (0, 8), (0, 9), (0, 10), (0, 11),
               (3, 2), (3, 4), (3, 6), (3, 8), (6, 2), (6, 6), (6, 7), (6, 11), (9, 0), (9, 4)])
        return [("eq. (w288)", code.join(u, 0))]
    return []


def census(code, d, engine):
    """All weight-d logicals -> classes -> (component, shape) tally.

    The cluster search in enumeration mode lists every minimum-weight logical
    through an anchor (the completeness part of Theorem "exactness"), then we
    close the list under the l*m translations and group by coset of S.
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
        comp = code.component(zs[0])          # the same for the whole class
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
    """eq. (milp): min 1^T z s.t. H_X z = 2p, L z - y = 2q, sum y >= 1."""
    import numpy as np
    from scipy.optimize import milp, LinearConstraint, Bounds
    n, N = code.n, code.N
    bits = lambda v: [(v >> i) & 1 for i in range(n)]
    # rows of H_X = [A | B]: check g touches left q iff bit g of A-column q
    HX = [[(code.qmask[q] >> g) & 1 for q in range(n)] for g in range(N)]
    # L: functionals that vanish on S but not on K, one per logical qubit.
    # (Any basis of X-logicals would do; we build it from S^perp mod K^perp.)
    Sperp = nullspace([sum(((s >> q) & 1) << r for r, s in enumerate(code.S.basis()))
                       for q in range(n)], n)
    Kperp = Span(nullspace([sum(((kv >> q) & 1) << r for r, kv in enumerate(code.K))
                            for q in range(n)], n))
    L = []
    for f in Sperp:
        if Kperp.add(f):
            L.append(f)
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-288", action="store_true", help="leave out the length-288 code")
    ap.add_argument("--engine", choices=["auto", "c", "py"], default="auto")
    ap.add_argument("--milp", action="store_true", help="MILP cross-check for n <= 108")
    ap.add_argument("--slow", action="store_true",
                    help="with the Python engine, also run the [[288]] census (~20+ min)")
    args = ap.parse_args()

    have_c = find_c_engine() is not None
    print("cluster engine:", "C (./cluster_lb)" if have_c and args.engine != "py" else "pure Python")
    codes = [c for c in CODES if not (args.skip_288 and c["name"] == "[[288,12,18]]")]
    if (not have_c or args.engine == "py") and not args.skip_288:
        print("  (pure Python: the [[288]] certificate takes a few minutes, its census is\n"
              "   skipped unless --slow; compile cluster_lb.c to make both quick)")

    for spec in codes:
        t_code = time.time()
        code = BBCode(spec["l"], spec["m"], spec["a"], spec["b"], spec["name"])
        name, d = spec["name"], spec["d"]
        big = name == "[[288,12,18]]"
        py_only = (not have_c) or args.engine == "py"
        print("\n=== %s   a = %s, b = %s on %dx%d ===" % (name, spec["a"], spec["b"], spec["l"], spec["m"]))

        # [1] structural part of Table "profile" ---------------------------
        print("  [1] exact-sequence profile")
        dec = code.decomposition()
        check("k", code.k, spec["k"])
        check("dim ann(a)", dec["dim_ann_a"], spec["dimann"])
        check("dim b ann(a)", dec["dim_b_ann_a"], spec["dimbann"])
        check("r_A = r_C (Lemma Frobenius)", (dec["r_A"], dec["r_C"]), (spec["rA"], spec["rA"]))
        check("r_A + r_C = k and (a:b) = ann(b ann(a))",
              dec["sequence_ok"] and dec["frobenius_ok"], True)
        # w_ann is an exhaustive minimum over ann(a) \ b ann(a); for the
        # 288 code that is 2^24 elements, a few seconds to a minute in Python
        w, t, exact = gray_min(dec["ann_a"], dec["b_ann_a"], max_dim=24)
        check("w_ann (exhaustive=%s)" % exact, w, spec["w_ann"])

        # [2] witnesses -------------------------------------------------------
        print("  [2] explicit logicals (upper bounds)")
        for label, z in witnesses(code, name):
            u, v = code.split(z)
            ok = code.is_logical(z)
            check("%s in K\\S, weight" % label, (ok, popcount(z)), (True, d))
            print("        split (%d,%d), component %s" % (popcount(u), popcount(v), code.component(z)))
        if name == "[[18,4,4]]":
            # Example "component is not shape": (1,b) has pi = 0 although it is
            # two-sided, and its one-sided representative weighs 6
            check("pi[(1,b)] = 0, i.e. (1,b) in A", code.component(code.join(1, code.b)), "A")

        # [3] lower bound certificate (Table "cluster certificates") --------
        print("  [3] anchored cluster search at radius d-1")
        t0 = time.time()
        status, hits, nodes = cluster_search(code, d - 1, engine=args.engine)
        check("radius %d: no logical  => d >= %d" % (d - 1, d), status, "none")
        print("        nodes = %d, %.1f s" % (nodes, time.time() - t0))

        # [4] census (Table "census") -------------------------------------
        print("  [4] census of minimum-weight classes")
        if big and py_only and not args.slow:
            print("      skipped in pure Python (use --slow, or the C engine)")
        else:
            nvec, ncls, tA, tC, nodes = census(code, d, args.engine)
            want = CENSUS[name]
            check("classes attaining d", ncls, want[0])
            check("A: one-sided/lopsided/balanced", tA, want[1])
            check("C: one-sided/lopsided/balanced", tC, want[2])
            print("        (%d minimum-weight vectors)" % nvec)
            # a component present in the census attains d; an absent one is > d
            if sum(tA):
                check("d_A = d (census)", d, spec["dA"])
            if sum(tC):
                check("d_C = d (census)", d, spec["dC"])

        # [5] the one code where the components differ ---------------------
        if name == "[[108,8,10]]":
            print("  [5] Example [[108,8,10]]: component minima")
            # w_K by the kernel-mode search (minimal nonzero elements of K
            # have no zero-syndrome proper subset, so the search is exact)
            wK = None
            for W in range(2, 8):
                st, _, _ = cluster_search(code, W, mode="kernel", engine=args.engine)
                if st == "found":
                    wK = W
                    break
            check("w_K (minimum nonzero weight of K)", wK, 6)
            # Proposition "irreducible range": radius 11 < 2 w_K lists every
            # logical of weight <= 11
            st, hits, _ = cluster_search(code, 2 * wK - 1, enumerate_all=True, engine=args.engine)
            weights = sorted({popcount(z) for z in hits})
            compA = [z for z in hits if code.component(z) == "A"]
            check("weights of all logicals of weight <= 11", weights, [10])
            check("annihilator logicals of weight <= 11", len(compA), 0)
            check("so d_A = w_ann", 12 if not compA else None, spec["dA"])
            st, _, _ = cluster_search(code, 9, mode="colon", engine=args.engine)
            check("colon search radius 9 (Corollary colon)", st, "none")
            print("        => d_C = 10 (witness eq. w108), d_A = 12")

        # [6] MILP cross-check ----------------------------------------------
        if args.milp and code.n <= 108:
            try:
                val, msg = milp_distance(code)
                check("MILP optimum (floating point, cross-check)", val, d)
            except ImportError:
                print("      scipy not available, MILP skipped")
        print("  (%.1f s for this code)" % (time.time() - t_code))

    print("\n%s" % ("all checks passed" if not FAILS else "MISMATCHES: " + ", ".join(FAILS)))
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
