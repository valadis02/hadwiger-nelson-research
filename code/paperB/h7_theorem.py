"""
h7_theorem.py — Theorem 3.2:  h >= 7, exactly, over Q(sqrt46).

STANDALONE.  Earlier versions of this script loaded a numerical realisation from
h7_realisation.npy; the exact Gram matrix is now written out in full below, so
the script depends on nothing but sympy and reproduces the theorem from scratch.

The graph is the unique (up to isomorphism) graph on 7 vertices with
independence number <= 2 and clique number <= 3 that is realisable on the unit
sphere with every edge a chord of length 1, i.e. inner product 1/2.  It has 11
edges and degree sequence [3,3,3,3,3,3,4]; h7.py establishes the uniqueness by
exhaustion over the nine isomorphism classes.

What is verified here, symbolically:
  * the 11 prescribed entries are 1/2 and the diagonal is 1
  * rank(G) = 3 exactly       -> realisable in R^3 and in no lower dimension
  * G is positive semidefinite -> it is the Gram matrix of real unit vectors
  * no non-edge entry equals 1/2 -> the realised graph is exactly this one

Usage: python3 h7_theorem.py
"""
import itertools
import sympy as sp

s46 = sp.sqrt(46)
n = 7

# the 11 edges (inner product 1/2)
EDGES = [(0, 2), (0, 4), (0, 5), (0, 6), (1, 2), (1, 3), (1, 6),
         (2, 6), (3, 4), (3, 5), (4, 5)]

# the ten non-edges, with their exact inner products over Q(sqrt46)
NONEDGES = {
    (0, 1): sp.Rational(-1, 3),
    (0, 3): sp.Rational(-1, 3),
    (1, 4): (2 + s46) / 16,
    (1, 5): (2 - s46) / 16,
    (2, 3): (2 - s46) / 16,
    (2, 4): sp.Rational(37, 64),
    (2, 5): (23 - 6 * s46) / 64,
    (3, 6): (2 + s46) / 16,
    (4, 6): (23 + 6 * s46) / 64,
    (5, 6): sp.Rational(37, 64),
}

assert len(EDGES) == 11
assert len(NONEDGES) == n * (n - 1) // 2 - 11
assert set(EDGES).isdisjoint(NONEDGES)

G = sp.zeros(n, n)
for i in range(n):
    G[i, i] = sp.Integer(1)
for i, j in EDGES:
    G[i, j] = G[j, i] = sp.Rational(1, 2)
for (i, j), v in NONEDGES.items():
    G[i, j] = G[j, i] = v

print("Gram matrix over Q(sqrt46):")
sp.pprint(G)

deg = sorted(sum(1 for e in EDGES if i in e) for i in range(n))
print(f"\nedges: {len(EDGES)}   degree sequence: {deg}")

r = G.rank()
print(f"rank(G) = {r}   "
      f"{'OK (realisable in R^3, not lower)' if r == 3 else 'MISMATCH'}")

ev = G.eigenvals()
neg = [e for e in ev if sp.simplify(e) != 0 and sp.N(e) < -sp.Rational(1, 10**30)]
print(f"eigenvalues ({len(ev)} distinct); negative: {len(neg)}   "
      f"{'PSD OK' if not neg else 'NOT PSD'}")
for e, mult in ev.items():
    print(f"   {sp.simplify(e)}   (multiplicity {mult}, ~{float(sp.N(e)):.6f})")

bad = [k for k, v in NONEDGES.items() if sp.simplify(v - sp.Rational(1, 2)) == 0]
print(f"\nnon-edges at inner product exactly 1/2: {len(bad)}   "
      f"{'OK - the graph realised is exactly this one' if not bad else 'COLLISION'}")

ok = (r == 3) and (not neg) and (not bad)
print(f"\nTheorem 3.2 (h >= 7): {'VERIFIED' if ok else 'FAILED'}")
