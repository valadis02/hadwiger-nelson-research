"""
h8_impossible.py — step D3.b3(b): the n = 8 candidate is not realisable.

By R(3,4) = 9 the only combinatorial candidate on 8 vertices is the complement
of the Wagner (Moebius-Kantor) graph, and it is unique up to isomorphism.  So a
single impossibility proof closes h <= 7.

The elimination.  A triangle of the graph is three unit vectors with pairwise
inner product 1/2; their Gram is M = (I+J)/2, which is invertible, so they are a
basis of R^3.  Write t_v = (<v,a>, <v,b>, <v,c>) for each remaining vertex.
Then every inner product is a polynomial in these:

        <u,v> = t_u^T M^{-1} t_v ,      |v|^2 = t_v^T M^{-1} t_v = 1

Entries of t_v are 1/2 where v is adjacent to a triangle vertex and unknown
otherwise, so the whole realisability question becomes a small polynomial
system.  We solve it exactly with Groebner bases: an empty variety is a proof of
non-realisability, no numerics involved.
"""
import itertools, sys
import sympy as sp

# ---- the unique candidate: complement of the Wagner graph V8 --------------
V8 = set()
for i in range(8):
    V8.add(frozenset((i, (i+1) % 8)))
for i in range(4):
    V8.add(frozenset((i, i+4)))
G = [p for p in itertools.combinations(range(8), 2) if frozenset(p) not in V8]
Eset = set(map(frozenset, G))
adj = {i: {j for j in range(8) if j != i and frozenset((i, j)) in Eset}
       for i in range(8)}
print(f"candidate: {len(G)} edges, degrees {[len(adj[i]) for i in range(8)]}")

# sanity: alpha <= 2 and omega <= 3
assert all(any(frozenset(p) in Eset for p in itertools.combinations(t, 2))
           for t in itertools.combinations(range(8), 3)), "alpha>2"
assert not any(all(frozenset(p) in Eset for p in itertools.combinations(q, 2))
               for q in itertools.combinations(range(8), 4)), "K4 present"
print("  alpha <= 2 and omega <= 3 confirmed")

# ---- pick a triangle -----------------------------------------------------
tri = None
for t in itertools.combinations(range(8), 3):
    if all(frozenset(p) in Eset for p in itertools.combinations(t, 2)):
        tri = t; break
assert tri, "no triangle"
print(f"  triangle basis: {tri}")

M = sp.Matrix([[1, sp.Rational(1,2), sp.Rational(1,2)],
               [sp.Rational(1,2), 1, sp.Rational(1,2)],
               [sp.Rational(1,2), sp.Rational(1,2), 1]])
Minv = M.inv()

rest = [v for v in range(8) if v not in tri]
T = {}
unknowns = []
for v in rest:
    comps = []
    for k, a in enumerate(tri):
        if a in adj[v]:
            comps.append(sp.Rational(1, 2))
        else:
            s = sp.Symbol(f"x_{v}_{a}", real=True)
            unknowns.append(s); comps.append(s)
    T[v] = sp.Matrix(comps)
for k, a in enumerate(tri):
    e = [0, 0, 0]; e[k] = 1
    T[a] = M * sp.Matrix(e)          # inner products of a with the triangle

print(f"  unknowns: {len(unknowns)}  ({[str(u) for u in unknowns]})")

def ip(u, v):
    return sp.expand((T[u].T * Minv * T[v])[0, 0])

eqs = []
for v in rest:
    eqs.append(sp.expand(ip(v, v) - 1))              # unit length
for i, j in itertools.combinations(range(8), 2):
    if i in tri and j in tri:
        continue
    if frozenset((i, j)) in Eset:
        eqs.append(sp.expand(ip(i, j) - sp.Rational(1, 2)))   # edge
eqs = [sp.nsimplify(e) for e in eqs if e != 0]
print(f"  equations: {len(eqs)}")
sys.stdout.flush()

print("\ncomputing Groebner basis ...")
gb = sp.groebner(eqs, *unknowns, order='grevlex')
print(f"  basis size {len(gb.exprs)}")
if list(gb.exprs) == [sp.Integer(1)] or gb.exprs == (sp.Integer(1),):
    print("\n  GROEBNER BASIS = {1}  ->  the system has NO solution, even over C")
    print("  => the 8-vertex candidate is NOT realisable, hence h <= 7.")
else:
    print("  basis is not {1}; the complex variety is non-empty.")
    for g in gb.exprs[:8]:
        print("   ", g)
    print("  -> non-realisability must then be argued over the REALS")
    print("     (or via positive semidefiniteness), not by emptiness alone.")
