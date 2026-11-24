"""Independent verification of Theorem chi({1,(sqrt2+sqrt6)/2}) >= 7.

Checks, in exact arithmetic (qf.py, Fractions only):
  1. the 8-point gadget (K8 minus a perfect matching) of the averaging step;
  2. every point of cert7.json is distinct;
  3. the edge list is exactly the set of pairs at distance 1 or d (nothing missing, nothing extra);
  4. every triangle (P;Q,R) is rigid: |PQ| = |PR| in {sqrt2, 1+sqrt3} and |QR| in {1, d};
and writes the CNF whose unsatisfiability is the finite part of the proof (conf7.cnf).
Then run  ./prove_chi7.sh  to refute it with a SAT solver and check the DRUP proof."""
import json, itertools, sys, os
from qf import *

HERE = os.path.dirname(os.path.abspath(__file__))

def gadget_check():
    i = C(ZERO, ONE)
    X = [CZERO, CONE, i, cadd(CONE, ZETA), cadd(CONE, i), cpow(ZETA, 5), cpow(ZETA, 10), cadd(cpow(ZETA, 2), i)]
    matched = {(1, 2), (4, 0), (3, 5), (6, 7)}
    matched |= {(b, a) for a, b in matched}
    legs = {R(2), R(4, 2)}            # sqrt2^2, (1+sqrt3)^2
    for a, b in itertools.combinations(range(8), 2):
        n = norm2(csub(X[a], X[b]))
        if (a, b) in matched: assert n in legs, (a, b)
        else: assert n in FORBIDDEN, (a, b)
    print("gadget: 24 pairs at distance 1 or d, 4 matched pairs at sqrt2 / 1+sqrt3  OK")

def main():
    gadget_check()
    c = json.load(open(os.path.join(HERE, "cert7.json")))
    G = generators()
    P = []
    for idx in c["points"]:
        z = CZERO
        for j in idx: z = cadd(z, G[j])
        P.append(z)
    n = len(P)
    E = set(tuple(sorted(e)) for e in c["edges"])
    found = set(); dup = 0
    for a, b in itertools.combinations(range(n), 2):
        nn = norm2(csub(P[a], P[b]))
        if nn == R(0): dup += 1
        if nn in FORBIDDEN: found.add((a, b))
    assert dup == 0, "duplicate points"
    assert found == E, f"edge list mismatch: {len(found ^ E)} pairs differ"
    print(f"points: {n} distinct; edges: {len(E)} = exactly the pairs at distance 1 or d  OK")
    legs = {R(2), R(4, 2)}
    for p, q, r in c["triangles"]:
        l1 = norm2(csub(P[p], P[q])); l2 = norm2(csub(P[p], P[r])); base = norm2(csub(P[q], P[r]))
        assert l1 == l2 and l1 in legs and base in FORBIDDEN, (p, q, r)
    print(f"triangles: {len(c['triangles'])} rigid  OK")
    # CNF: x_{v,j} exactly one of 6 colours; edges proper; y_{uv,j} -> x_uj & x_vj; triangle clauses
    K = 6; var = lambda v, j: v * K + j + 1; nv = n * K; F = []
    for v in range(n):
        F.append([var(v, j) for j in range(K)])
        F += [[-var(v, a), -var(v, b)] for a in range(K) for b in range(a + 1, K)]
    for a, b in sorted(E):
        F += [[-var(a, j), -var(b, j)] for j in range(K)]
    eqv = {}
    def eq(p, q):
        nonlocal nv
        key = (min(p, q), max(p, q))
        if key not in eqv:
            ys = []
            for j in range(K):
                nv += 1; ys.append(nv); F.append([-nv, var(p, j)]); F.append([-nv, var(q, j)])
            eqv[key] = ys
        return eqv[key]
    for p, q, r in c["triangles"]: F.append(eq(p, q) + eq(p, r))
    a, b = c["edges"][0]
    F += [[var(a, 0)], [var(b, 1)]]      # colour symmetry: all constraints are invariant under permuting colours
    out = os.path.join(HERE, "conf7.cnf")
    with open(out, "w") as f:
        f.write(f"p cnf {nv} {len(F)}\n")
        for cl in F: f.write(" ".join(map(str, cl)) + " 0\n")
    print(f"CNF written: {out}  ({nv} variables, {len(F)} clauses)")

if __name__ == "__main__":
    main()
