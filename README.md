# BB-Codes-Distance-Logical-Decomposition

Companion code for

> M. Rowshan and S. Devitt, *Logical Operator Decomposition for Distance
> Analysis of Bivariate Bicycle Codes*.

Repository: <https://github.com/mohammad-rowshan/BB-Codes-Distance-Logical-Decomposition>

The code does two things. It reproduces the numbers reported in the paper
(exact-sequence profiles, witnesses, distance certificates, the census of
minimum-weight logicals), and it gives upper and lower bounds on the distance
of **any** bivariate bicycle (BB) code you hand it, proving the distance
outright whenever the two bounds meet.

---

## Contents

1. [Files](#1-files)
2. [Requirements and installation](#2-requirements-and-installation)
3. [Quick start](#3-quick-start)
4. [Bounds for your own BB code: `bbcode.py`](#4-bounds-for-your-own-bb-code-bbcodepy)
5. [The search engine: `cluster_lb`](#5-the-search-engine-cluster_lb)
6. [Reproducing the paper: `reproduce.py`](#6-reproducing-the-paper-reproducepy)
7. [Self-test: `selftest.py`](#7-self-test-selftestpy)
8. [Using the Python library directly](#8-using-the-python-library-directly)
9. [Conventions](#9-conventions)
10. [How the bounds are certified](#10-how-the-bounds-are-certified)
11. [Performance and runtimes](#11-performance-and-runtimes)
12. [Limitations](#12-limitations)
13. [Troubleshooting](#13-troubleshooting)
14. [Citation](#14-citation)

---

## 1. Files

| file | role | paper |
|---|---|---|
| `bbcode.py` | Python library and command-line tool: builds the code, computes the exact sequence, finds witnesses, runs the cluster search, reports bounds | Sec. II to VI |
| `cluster_lb.c` | fast C implementation of the anchored cluster search, with enumeration, colon and kernel modes | Sec. V, Algorithm "anchored cluster search" |
| `reproduce.py` | regenerates the tables of the paper and compares every value with the published one | Tables "profile", "witnesses", "cluster certificates", "census", Example [[108,8,10]] |
| `selftest.py` | stabilizer sanity checks and a regression test against brute force | Appendix A |
| `README.md` | this file | |

## 2. Requirements and installation

* Python 3.10 or newer (the code uses `int.bit_count`). No third-party packages
  are needed for anything except the optional MILP cross-check.
* A C compiler (gcc or clang) for `cluster_lb.c`. Optional, but it is 15 to 20
  times faster than the Python fallback and makes the length-288 code quick.
* Optional: `numpy` and `scipy >= 1.9` for the MILP cross-check (`reproduce.py --milp`).

```sh
git clone https://github.com/mohammad-rowshan/BB-Codes-Distance-Logical-Decomposition.git
cd BB-Codes-Distance-Logical-Decomposition
gcc -O3 -march=native -o cluster_lb cluster_lb.c
python3 selftest.py          # should end with "self-test passed"
```

On macOS use `clang -O3 -o cluster_lb cluster_lb.c` if `-march=native` is
rejected. On Windows, WSL is the simplest route; with MinGW the same gcc line
works and produces `cluster_lb.exe` (rename it to `cluster_lb` or put it on
your `PATH`).

The Python scripts look for the binary next to `bbcode.py` first and then on
the `PATH`. If none is found they silently use the pure-Python engine.

## 3. Quick start

```sh
# distance of the Gross code, proved
python3 bbcode.py -l 12 -m 6 -a "x^3+y+y^2" -b "y^3+x+x^2"

# regenerate the paper's tables (about two minutes with the C engine)
python3 reproduce.py

# the same without the length-288 code (a few seconds)
python3 reproduce.py --skip-288
```

## 4. Bounds for your own BB code: `bbcode.py`

```
python3 bbcode.py -l L -m M -a "POLY" -b "POLY" [--seconds T] [--engine auto|c|py]
```

| option | meaning | default |
|---|---|---|
| `-l`, `-m` | torus sizes, R = F2[x,y]/(x^l - 1, y^m - 1) | required |
| `-a`, `-b` | the two defining polynomials | required |
| `--seconds` | time budget for the lower-bound search | 60 |
| `--engine` | `auto` uses C if available, `c` insists on it, `py` forces Python | `auto` |

### Polynomial syntax

Terms are joined by `+`. A term is `1`, a power of `x`, a power of `y`, or a
product of the two, with or without `*`:

```
x^3+y+y^2        1+x^2+x^7        x^2y^3+xy        x^2*y^3 + y^5
```

Exponents are reduced modulo `l` and `m`, and repeated terms cancel, since the
coefficients live in F2 (so `x+x` is `0`).

### What it does

1. Builds A, B, H_X = [A | B], H_Z = [B^T | A^T], K = ker H_X and S = im H_Z^T,
   and reports n and k.
2. Computes the exact sequence of Theorem 1: dim ann(a), dim b ann(a), r_A,
   dim (a:b), dim (a), r_C, and checks r_A + r_C = k (Corollary "basis") and
   (a:b) = ann(b ann(a)) with r_A = r_C (Lemma "Frobenius balance").
3. **Upper bound.** Runs the algebraic screens of Algorithm 1: least-weight
   one-sided logicals (t,0) with t in ann(a) \ b ann(a), and (0,s) on the
   other block, plus sparse colon lifts a u = b v with the coset-leader weight
   minimised over u0 + ann(a). The lightest logical found is d_U.
4. **Lower bound.** Runs the anchored cluster search at radius W = 1, 2, 3, ...
   A radius with no hit certifies d >= W+1. The first radius with a hit gives
   the distance exactly, since all smaller radii were already excluded. It
   stops when the bounds meet or when the time budget is spent.

### Reading the output

```
BB code on 12x6 torus: a = y+y^2+x^3, b = y^3+x+x^2
n = 144, k = 12
exact sequence (Theorem 1): dim ann(a)=12, dim b ann(a)=6, r_A=6; dim(a:b)=66, dim(a)=60, r_C=6
checks: r_A + r_C = k -> True ; Frobenius (a:b)=ann(b ann(a)) and r_A=r_C -> True
  upper bound from screens: 12 (one-sided left)
  cluster search radius  1: none  nodes=2
  ...
  cluster search radius 11: none  nodes=284281

distance proved: d = 12  (d_X = d_Z by Proposition 'X/Z symmetry')
witness: u = ... ; v = 0 ; split (12,0) ; component A
```

* `distance proved: d = ...` means an explicit logical meets a certified lower
  bound.
* `bounds only: d_L <= d <= d_U` means the budget ran out first. The lower
  bound is still a proof; the upper bound is the lightest logical found.
* The witness line gives the two blocks as polynomials, the weight split and
  the component under pi[(u,v)] = v + (a): `A` (annihilator, pi = 0) or `C`
  (colon, pi != 0).

### Example on a longer code

```sh
python3 bbcode.py -l 30 -m 6 -a "x^9+y+y^2" -b "y^3+x^25+x^26" --seconds 600
```

With 60 s this gives `18 <= d <= 30`. The lower bound grows with the budget;
the upper bound from the screens is weak on this code (a weight-24 logical is
known), see [Limitations](#12-limitations).

## 5. The search engine: `cluster_lb`

```
./cluster_lb [-e] [-c] [-k] [-n BUDGET] l m W "a" "b"
```

| option | meaning | paper |
|---|---|---|
| (none) | stop at the first logical of weight <= W, or certify d >= W+1 | Theorem "exactness of the cluster search" |
| `-e` | enumerate: print every logical of weight <= W reached through an anchor | census, Table "census" |
| `-c` | colon mode: a hit must have pi != 0; certifies d_C >= W+1 | Corollary "colon connectivity" |
| `-k` | kernel mode: any nonzero element of K is a hit; gives w_K | Proposition "irreducible range" |
| `-n BUDGET` | stop after BUDGET nodes and certify nothing | |

### Output format

```
# BB code l=9 m=6 n=108 k=8 dimS=50 check_weight=6 radius=10 mode=logical
LOGICAL 10 0 58 1 39 46 67 22 74 24 89
RESULT found radius=10 weight=10 nodes=3133 time=0.00
```

* `#` lines are informational.
* `LOGICAL w q1 ... qw` lists a support. Qubit `q < N` is the left-block
  monomial with index q, `q >= N` the right-block monomial with index q - N,
  where N = l*m and index i*m + j stands for x^i y^j.
* `RESULT none ... => d >= W+1` is a certificate.
* `RESULT found`, `RESULT enumerated`, `RESULT budget` are self-explanatory.
  A budget stop certifies nothing about radius W.

### Typical calls

```sh
./cluster_lb 12 12 17 "x^3+y^2+y^7" "y^3+x+x^2"   # [[288,12,18]]: d >= 18
./cluster_lb 12 12 18 "x^3+y^2+y^7" "y^3+x+x^2"   # finds a weight-18 logical
./cluster_lb -e 9 6 10 "x^3+y+y^2" "y^3+x+x^2"    # all weight-10 logicals of [[108,8,10]]
./cluster_lb -c 9 6 9 "x^3+y+y^2" "y^3+x+x^2"     # d_C >= 10
./cluster_lb -k 9 6 6 "x^3+y+y^2" "y^3+x+x^2"     # w_K = 6
```

Enumeration mode lists supports through the anchors only. To get all
minimum-weight logicals, close the list under the l*m translations, as
`reproduce.py` does.

Enumeration is exact at radius W = d. For W > d it is exact as long as
W < 2 w_K (Proposition "irreducible range"); beyond that, logicals whose
support splits into two zero-syndrome pieces can be missed.

## 6. Reproducing the paper: `reproduce.py`

```
python3 reproduce.py [--skip-288] [--engine auto|c|py] [--milp] [--slow]
```

| option | meaning |
|---|---|
| `--skip-288` | leave out [[288,12,18]] |
| `--engine` | as for `bbcode.py` |
| `--milp` | add the parity-constrained MILP cross-check of eq. (milp) for n <= 108 (needs scipy; floating point, so a cross-check only) |
| `--slow` | with the Python engine, also run the [[288]] census (about 20 minutes or more) |

For each of the six standard codes the script prints one line per quantity,
with `ok` or `MISMATCH (manuscript: ...)`, and exits with status 0 only if
everything matches:

| block | checks | paper |
|---|---|---|
| `[1]` | k, dim ann(a), dim b ann(a), r_A = r_C, the exact-sequence and Frobenius identities, exhaustive w_ann | Table "profile", Lemma "Frobenius balance" |
| `[2]` | every listed witness is in K \ S with weight d; split and component; for [[18,4,4]], (1,b) lies in A | Table "witnesses", Appendix B, Example "component is not shape" |
| `[3]` | cluster search at radius d-1 finds nothing; node count | Table "cluster certificates" |
| `[4]` | number of minimum-weight classes and their component/shape counts; d_A, d_C where attained | Table "census", Table "profile" |
| `[5]` | [[108,8,10]]: w_K = 6, no logical of weight 11, d_C = 10, d_A = 12 | Example [[108,8,10]], Proposition "irreducible range" |
| `[6]` | MILP optimum (only with `--milp`) | Section "Results of the lower-bound computations" |

Not reproduced: the random-pair screening statistics of Table "sweep" (the
exact sampler of the paper is not part of this package).

## 7. Self-test: `selftest.py`

```
python3 selftest.py [--count 200] [--seed 2026] [--max-n 32]
```

1. For the six standard codes: every basis vector of S lies in K, every
   translated stabilizer generator (b x^i y^j, a x^i y^j) reduces to zero
   modulo the stored basis of S, and dim K - dim S = k.
2. On `--count` random BB codes with n <= `--max-n` (random weight-2 and
   weight-3 polynomials on small tori), the cluster search distance equals the
   brute-force minimum over K \ S, for each available engine.
3. The C and Python engines visit exactly the same number of nodes on every
   test code, since they use the same canonical branching rule.

Run it after changing any of the code.

## 8. Using the Python library directly

```python
from bbcode import BBCode, cluster_search, distance_bounds, describe_witness

code = BBCode(12, 6, "x^3+y+y^2", "y^3+x+x^2")   # l, m, a, b
print(code.n, code.k)                            # 144 12

dec = code.decomposition()                       # Theorem 1 data
print(dec["r_A"], dec["r_C"], dec["frobenius_ok"])   # 6 6 True

status, hits, nodes = cluster_search(code, 11)   # radius 11
print(status, nodes)                             # none 284281  -> d >= 12

status, hits, nodes = cluster_search(code, 12)   # finds a weight-12 logical
z = hits[0]
print(describe_witness(code, z))                 # polynomials, split, component
print(code.is_logical(z), code.component(z))     # True, 'A' or 'C'

dL, dU, witness, notes = distance_bounds(code, seconds=30, verbose=False)
print(dL, dU)                                    # 12 12
```

Useful pieces:

| call | returns |
|---|---|
| `BBCode(l, m, a, b)` | the code; `a`, `b` as strings or lists of exponent pairs |
| `code.K`, `code.S` | basis of K (list of ints) and the echelon `Span` of S |
| `code.decomposition()` | dict with the dimensions of Theorem 1, `r_A`, `r_C`, the bases of ann(a) and (a:b), and the two identity checks |
| `code.split(z)`, `code.join(u, v)` | convert between a length-n vector and its two blocks |
| `code.is_logical(z)` | z in K \ S |
| `code.pi(z)`, `code.component(z)` | v + (a) reduced, and `'A'` or `'C'` |
| `code.class_key(z)` | canonical representative of z + S, to group logicals by class |
| `code.translate(z, i, j)` | multiply both blocks by x^i y^j |
| `code.mul(f, g)`, `code.fmt(f)` | product in R and pretty printing |
| `cluster_search(code, W, mode, enumerate_all, node_budget, engine)` | `(status, hits, nodes)`; `mode` is `"logical"`, `"colon"` or `"kernel"` |
| `one_sided_witnesses(code)`, `colon_lift_witnesses(code)` | lists of `(weight, z, description, exact_flag)` |
| `distance_bounds(code, seconds, engine, verbose)` | `(d_L, d_U, witness, notes)` |

Vectors are Python integers used as bitsets. Bit `i*m + j` is the coefficient
of x^i y^j, and a two-block vector (u, v) is `u | (v << N)`.

## 9. Conventions

* **Indexing.** x^i y^j is coordinate i*m + j; the left block occupies
  0..N-1 and the right block N..2N-1, with N = l*m.
* **Check matrices.** H_X = [A | B] and H_Z = [B^T | A^T], with A the matrix
  of multiplication by a in the column-vector convention (Remark "stabilizer
  convention"). Hence K = {(u,v) : a u + b v = 0} and S = {(b r, a r)}.
  Some papers use the row convention, which replaces a, b by their reciprocals;
  distances and all dimensions are unchanged.
* **Which distance.** Everything is computed for Z-type logicals. For BB codes
  d_X = d_Z (Proposition "X/Z symmetry"), so d_Z is the code distance.
* **Components.** A and C refer to the fixed projection pi[(u,v)] = v + (a).
  Swapping the roles of a and b gives the companion labels of Corollary
  "symmetric sequence".
* **Shapes.** In the census, a class is one-sided if some minimum-weight
  representative has an empty block, tau-lopsided if not and some minimum
  representative has a block of weight <= tau, and tau-balanced otherwise,
  with tau = max(wt a, wt b) (Definition "representative shape").

## 10. How the bounds are certified

* **Upper bounds** are explicit vectors z with H_X z = 0 and z not in S. Both
  properties are checked by exact elimination before a witness is reported.
* **Lower bounds** rest on Lemma "syndrome connectivity": every proper subset
  of a minimum-weight logical has nonzero syndrome. The search therefore grows
  a support from an anchor qubit and always branches on the qubits of the
  lowest-index unsatisfied X-check, which never loses a minimum-weight logical.
  Translation symmetry needs only two anchors: the left-block origin, plus the
  right-block origin for supports that live only on the right block.
* Two prunes keep the tree small. A zero-syndrome support that lies in S is
  dropped, by the lemma. A branch is also dropped when |T| + ceil(u/c) > W,
  where u is the number of unsatisfied checks and c the largest number of
  checks meeting one qubit.
* Everything is exact bit and integer arithmetic. The MILP of `--milp` uses a
  floating-point solver and is reported only as a cross-check.

The cluster idea is due to Dumer, Kovalev and Pryadko (IEEE Trans. Inf.
Theory 63(7), 2017). The anchoring, the parity prune, and the colon and kernel
variants are specific to this work.

## 11. Performance and runtimes

Node counts are machine independent; times below are for one core and will
vary.

| code | radius | nodes | C | Python |
|---|---|---|---|---|
| [[18,4,4]] | 3 | 19 | < 0.1 s | < 0.1 s |
| [[72,12,6]] | 5 | 171 | < 0.1 s | < 0.1 s |
| [[90,8,10]] | 9 | 29,301 | < 0.1 s | < 0.1 s |
| [[108,8,10]] | 9 | 24,863 | < 0.1 s | < 0.1 s |
| [[144,12,12]] | 11 | 284,281 | < 0.1 s | about 0.2 s |
| [[288,12,18]] | 17 | 455,504,708 | 7 to 13 s | about 4 min |

The worst-case tree size is 2 * 5^(W-1) for weight-six checks, independent of
n, but the parity prune keeps the practical count far smaller. Each extra unit
of radius costs roughly a factor of 2 to 6 on these codes, so distances in the
low twenties are feasible on a workstation.

Other runtimes:

* `reproduce.py` in full takes about two minutes with the C engine. Most of
  that time goes to the [[288]] census, about 2.5e9 nodes.
* `selftest.py` with the default 200 codes takes a few seconds.
* The exhaustive w_ann for [[288]] (2^24 elements) takes a few seconds in
  Python.

## 12. Limitations

* **Upper bounds on long codes.** The algebraic screens only look at one-sided
  logicals and sparse colon lifts. On long codes the lightest logical can be
  two-sided and dense on both blocks, and then d_U from the screens is loose
  until the cluster search itself reaches d. For better upper bounds on such
  codes, pair this package with QDistRnd or a BP-OSD sampler and use the
  cluster search for the lower bound.
* **Exponential cost in d.** The search is exponential in the target distance.
  Past d of about 22 to 24 you will want a larger budget, several cores (the
  two anchors and the first-level branches are independent jobs), or
  reduction by the full automorphism group.
* **Enumeration beyond d.** Complete only up to radius 2 w_K - 1, see
  Section 5.
* **d_A in general.** There is no annihilator analogue of the colon
  connectivity result. d_A is read off the census when it equals d, and
  otherwise obtained with Proposition "irreducible range" plus a one-sided
  witness, as for [[108,8,10]].
* **Size.** The C engine accepts any l, m and polynomials with up to 32 terms
  each. Memory is small; the Python engine keeps the recursion depth at W.

## 13. Troubleshooting

| symptom | fix |
|---|---|
| `cluster engine: pure Python` although you compiled | the binary must be named `cluster_lb` and sit next to `bbcode.py` or on the `PATH`, with execute permission |
| `cannot parse polynomial` | use only `x`, `y`, `^`, `*`, digits and `+`; no spaces inside exponents, no minus signs |
| `k = 0, no logical qubits` | this pair (a, b) encodes nothing on this torus; there is nothing to bound |
| `RESULT budget` or `bounds only` | raise `--seconds` (or `-n` for `cluster_lb`), or use the C engine |
| MILP check prints `scipy not available` | `pip install numpy scipy`, or run without `--milp` |
| `AttributeError: 'int' object has no attribute 'bit_count'` | Python older than 3.10; upgrade |

## 14. Citation

```bibtex
@article{rowshan_devitt_bbdistance,
  author  = {Rowshan, Mohammad and Devitt, Simon},
  title   = {Logical Operator Decomposition for Distance Analysis of Bivariate Bicycle Codes},
  journal = {submitted},
  year    = {2026}
}
```

Please also cite I. Dumer, A. A. Kovalev and L. P. Pryadko, "Distance
verification for classical and quantum LDPC codes," IEEE Trans. Inf. Theory
63(7), 2017, when using the cluster search.
