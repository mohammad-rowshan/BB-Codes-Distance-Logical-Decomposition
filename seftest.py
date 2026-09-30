#!/usr/bin/env python3
"""
selftest.py -- the consistency checks behind Sec. V-C and Appendix A.

  1. For the six standard codes (Table IX [tab:codedata]): every basis vector
     of S lies in K, every translated stabilizer generator (b x^i y^j, a x^i y^j)
     reduces to zero modulo the stored basis of S, and dim K - dim S = k.
  2. Regression: on random small BB codes, the distance returned by the
     anchored cluster search [alg:cluster] equals the brute-force minimum over
     K \\ S. This is the "200 random BB codes of length at most 32" of Sec. V-C.
     Both engines are tested when the C binary exists.
  3. The C and Python engines visit exactly the same number of nodes, so the
     node counts of Table IV [tab:clustercert] do not depend on the engine.

    python3 selftest.py            # 200 random codes, a few seconds
    python3 selftest.py --count 1000 --seed 7
"""

import argparse
import random
import sys

from bbcode import BBCode, cluster_search, find_c_engine, popcount
from reproduce import CODES


def brute_force_distance(code):
    """min weight over K \\ S, eq. (dz), by Gray-code enumeration of K (small codes only)"""
    best, cur = None, 0
    K = code.K
    for g in range(1, 1 << len(K)):
        cur ^= K[(g & -g).bit_length() - 1]
        w = popcount(cur)
        if (best is None or w < best) and cur not in code.S:
            best = w
    return best


def cluster_distance(code, engine):
    """radius 1, 2, ...; the first radius with a hit is d (Theorem 4 [thm:cluster])"""
    for W in range(1, code.n + 1):
        status, hits, _ = cluster_search(code, W, engine=engine)
        if status == "found":
            return W
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--count", type=int, default=200, help="number of random codes")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--max-n", type=int, default=32, help="largest block length for brute force")
    args = ap.parse_args()
    fails = 0

    # 1. structural checks on the standard codes -----------------------------
    print("[1] stabilizer checks on the six standard codes")
    for spec in CODES:
        code = BBCode(spec["l"], spec["m"], spec["a"], spec["b"])
        gens_ok = True
        for r in range(code.N):
            z = code.join(code.Bcols[r], code.Acols[r])      # (b x^r, a x^r)
            if not code.in_K(z) or z not in code.S:
                gens_ok = False
        basis_ok = all(code.in_K(s) for s in code.S.basis())
        k_ok = code.k == spec["k"]
        ok = gens_ok and basis_ok and k_ok
        fails += not ok
        print("    %-15s S in K: %s, generators reduce to 0 mod S: %s, k = %d: %s"
              % (spec["name"], basis_ok, gens_ok, code.k, "ok" if ok else "FAIL"))

    # 2. + 3. regression against brute force, and engine agreement --------------
    engines = ["py"] + (["c"] if find_c_engine() else [])
    print("\n[2] cluster search vs brute force on %d random codes (engines: %s)"
          % (args.count, ", ".join(engines)))
    rng = random.Random(args.seed)
    tori = [(l, m) for l in range(2, 9) for m in range(2, 9) if 2 * l * m <= args.max_n]
    tested = node_mismatch = dist_mismatch = 0
    while tested < args.count:
        l, m = rng.choice(tori)
        mons = [(i, j) for i in range(l) for j in range(m)]
        wa, wb = rng.choice([2, 3]), rng.choice([2, 3])
        code = BBCode(l, m, rng.sample(mons, wa), rng.sample(mons, wb))
        if code.k == 0:
            continue
        d_true = brute_force_distance(code)
        for eng in engines:
            d_cl = cluster_distance(code, eng)
            if d_cl != d_true:
                dist_mismatch += 1
                print("    MISMATCH (%s engine): l=%d m=%d a=%s b=%s brute=%s cluster=%s"
                      % (eng, l, m, code.fmt(code.a), code.fmt(code.b), d_true, d_cl))
        if len(engines) == 2:
            W = max(1, d_true - 1)
            n_py = cluster_search(code, W, engine="py")[2]
            n_c = cluster_search(code, W, engine="c")[2]
            node_mismatch += n_py != n_c
        tested += 1
    print("    %d codes tested, %d distance mismatches" % (tested, dist_mismatch))
    fails += dist_mismatch
    if len(engines) == 2:
        print("[3] C vs Python node counts differ on %d of %d codes" % (node_mismatch, tested))
        fails += node_mismatch

    print("\n%s" % ("self-test passed" if fails == 0 else "self-test FAILED (%d problems)" % fails))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
