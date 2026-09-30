# bb-distance-decomposition: reproduction code

Scripts for *Logical Operator Decomposition for Distance Analysis of Bivariate
Bicycle Codes* (M. Rowshan, S. Devitt).

| file | what it is |
|---|---|
| `cluster_lb.c` | fast engine for the anchored cluster search (Sec. V, Algorithm "anchored cluster search") |
| `bbcode.py` | library + command line tool: exact sequence, witnesses, distance bounds for any BB code |
| `reproduce.py` | regenerates Tables "profile", "witnesses", "cluster certificates", "census" and Example [[108,8,10]] |

## Build and run

```sh
gcc -O3 -march=native -o cluster_lb cluster_lb.c   # optional but much faster
python3 reproduce.py              # all six codes, about 2 min with the C engine
python3 reproduce.py --skip-288   # a few seconds
python3 reproduce.py --milp       # adds the MILP cross-check (needs scipy)
```

Bounds for an arbitrary BB code:

```sh
python3 bbcode.py -l 12 -m 6 -a "x^3+y+y^2" -b "y^3+x+x^2"
python3 bbcode.py -l 30 -m 6 -a "x^9+y+y^2" -b "y^3+x^25+x^26" --seconds 600
```

The tool prints `distance proved: d = ...` when the explicit witness meets the
certified lower bound, and `bounds only: d_L <= d <= d_U` when the time budget
runs out first. The C engine can also be called directly:

```sh
./cluster_lb 12 12 17 "x^3+y^2+y^7" "y^3+x+x^2"      # certifies d >= 18
./cluster_lb -e 9 6 10 "x^3+y+y^2" "y^3+x+x^2"       # list all weight-10 logicals
./cluster_lb -c 9 6 9 "x^3+y+y^2" "y^3+x+x^2"        # colon component: d_C >= 10
./cluster_lb -k 9 6 6 "x^3+y+y^2" "y^3+x+x^2"        # w_K, min weight of K
```

## Conventions

Monomial x^i y^j is coordinate `i*m + j`; the left block comes first, so qubit
`q < N` is left and `q >= N` is right (N = l*m). H_X = [A | B] in the
column-vector convention of the manuscript (Remark "stabilizer convention"),
so S = {(b r, a r)} and pi[(u,v)] = v + (a). Only d_Z is computed; d_X = d_Z
for every BB code (Proposition "X/Z symmetry").

## Reference runtimes (one core)

The [[288,12,18]] certificate at radius 17 visits 455,504,708 nodes, about
7 to 13 s in C and about 4 min in pure Python. Node counts are machine
independent and are the numbers to compare against the manuscript.
