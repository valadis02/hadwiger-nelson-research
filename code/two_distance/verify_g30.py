"""Independent verification of the 30-vertex 6-chromatic {1,d}-graph G30 (exact arithmetic).

1. edges = exactly the pairs at distance 1 or d;   2. the stored colouring is a proper 6-colouring;
3. no proper 5-colouring (SAT refutation + DRUP proof checked by rup_checker);
4. vertex-critical: G30 - v has a proper 5-colouring for every v (each checked edge by edge)."""
import json, itertools, os, subprocess
from fractions import Fraction as Fr
from qf import *
from pysat.solvers import Solver

HERE = os.path.dirname(os.path.abspath(__file__))
g = json.load(open(os.path.join(HERE, "g30.json")))
P = [C(R(*[Fr(x, 8) for x in v["a"]]), R(*[Fr(x, 8) for x in v["b"]])) for v in g["vertices"]]
n = len(P); E = set(tuple(sorted(e)) for e in g["edges"])
found = {(a, b) for a, b in itertools.combinations(range(n), 2) if norm2(csub(P[a], P[b])) in FORBIDDEN}
assert found == E
assert all(norm2(csub(P[a], P[b])) != R(0) for a, b in itertools.combinations(range(n), 2))
units = sum(1 for a, b in E if norm2(csub(P[a], P[b])) == ONE)
print(f"G30: {n} vertices, {len(E)} edges ({units} of length 1, {len(E)-units} of length d) = all pairs at 1 or d  OK")
col = g["colouring6"]; assert all(col[a] != col[b] for a, b in E) and max(col) < 6
print("stored 6-colouring proper  OK")

def cnf(verts, k):
    vs = list(verts); ix = {v: i for i, v in enumerate(vs)}
    var = lambda v, j: ix[v] * k + j + 1
    F = [[var(v, j) for j in range(k)] for v in vs]
    F += [[-var(a, j), -var(b, j)] for a, b in E if a in ix and b in ix for j in range(k)]
    return F, var

F, _ = cnf(range(n), 5)
with open(os.path.join(HERE, "g30_k5.cnf"), "w") as f:
    f.write(f"p cnf {n*5} {len(F)}\n" + "".join(" ".join(map(str, c)) + " 0\n" for c in F))
s = Solver(name="glucose4", bootstrap_with=F, with_proof=True)
assert s.solve() is False
open(os.path.join(HERE, "g30_k5.drup"), "w").write("\n".join(s.get_proof()) + "\n")
subprocess.run(["gcc", "-O2", "-o", os.path.join(HERE, "rup_checker"), os.path.join(HERE, "rup_checker.c")], check=True)
out = subprocess.run([os.path.join(HERE, "rup_checker"), os.path.join(HERE, "g30_k5.cnf"), os.path.join(HERE, "g30_k5.drup")],
                     capture_output=True, text=True).stdout
assert "VERIFIED UNSAT" in out
print("no proper 5-colouring: DRUP proof checked  OK")
for v in range(n):
    rest = [w for w in range(n) if w != v]
    F, var = cnf(rest, 5)
    s = Solver(name="glucose4", bootstrap_with=F); assert s.solve()
    m = set(l for l in s.get_model() if l > 0)
    c5 = {w: next(j for j in range(5) if var(w, j) in m) for w in rest}
    assert all(c5[a] != c5[b] for a, b in E if a != v and b != v)
print("vertex-critical: 30/30 proper 5-colourings of G30 - v  OK")
