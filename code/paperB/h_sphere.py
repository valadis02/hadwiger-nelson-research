"""
h_sphere.py  -  step D3.b2:  the value of h on the unit sphere S^2(1).

h = max number of points on the unit sphere, edges = chords of length 1
    (i.e. 60 degrees, inner product 1/2), with independence number alpha <= 2.

TWO STRUCTURAL FACTS, both proved, that make the search finite and tiny:

  (F1) K4-free.  Three unit vectors with pairwise inner product 1/2 have Gram
       matrix (I+J)/2, positive definite, hence linearly independent; a fourth
       cannot be added in R^3.  So omega(G) <= 3.                (report sec 14.4)

  (F2) alpha(G) <= 2  <=>  complement(G) is triangle-free.

  Together: G has no K4 and no independent 3-set, so by R(3,4) = 9 we get
       h <= 8,
  which recovers the recorded upper bound as a one-line Ramsey argument rather
  than an assumption.  Moreover the (3,4)-Ramsey graph on 8 vertices is unique
  (the complement of the Wagner / Moebius-Kantor graph V8), so at n = 8 there is
  exactly ONE combinatorial candidate to test for spherical realisability.

This script tests realisability of the candidates numerically with many random
restarts; a residual driven to ~1e-16 is a (numerical) realisation, a floor far
from zero over many restarts is evidence of non-realisability, to be upgraded to
a proof by exact means for whichever n turns out to be the answer.
"""
import itertools, numpy as np
from scipy.optimize import least_squares

def alpha_le2(n, edges):
    E = set(map(frozenset, edges))
    for t in itertools.combinations(range(n), 3):
        if all(frozenset(p) not in E for p in itertools.combinations(t, 2)):
            return False
    return True

def omega_le3(n, edges):
    E = set(map(frozenset, edges))
    for q in itertools.combinations(range(n), 4):
        if all(frozenset(p) in E for p in itertools.combinations(q, 2)):
            return False
    return True

def realisable(n, edges, tries=300, seed=0):
    """min residual over restarts; vectors unit, inner product 1/2 on edges."""
    rng = np.random.default_rng(seed)
    E = list(edges)
    NE = [p for p in itertools.combinations(range(n), 2) if frozenset(p) not in set(map(frozenset, E))]
    def resid(x):
        V = x.reshape(n, 3)
        r = [np.dot(V[i], V[i]) - 1 for i in range(n)]
        r += [np.dot(V[i], V[j]) - 0.5 for i, j in E]
        return np.array(r)
    best = None; bestV = None
    for t in range(tries):
        x0 = rng.normal(size=3 * n)
        s = least_squares(resid, x0, xtol=1e-15, ftol=1e-15, gtol=1e-15)
        c = float(np.max(np.abs(s.fun)))
        if best is None or c < best:
            best, bestV = c, s.x.reshape(n, 3)
    # check the non-edges really avoid inner product 1/2 (else it is a different graph)
    sep = None
    if best is not None and best < 1e-10 and NE:
        sep = min(abs(float(np.dot(bestV[i], bestV[j])) - 0.5) for i, j in NE)
    return best, sep

# ---- Wagner graph V8 = C8 plus the four diameters; G = its complement
V8 = set()
for i in range(8):
    V8.add(frozenset((i, (i + 1) % 8)))
for i in range(4):
    V8.add(frozenset((i, i + 4)))
G8 = [p for p in itertools.combinations(range(8), 2) if frozenset(p) not in V8]
print(f"n=8 candidate: complement of Wagner graph, {len(G8)} edges, "
      f"alpha<=2: {alpha_le2(8, G8)}, omega<=3: {omega_le3(8, G8)}")
b, s = realisable(8, G8, tries=200, seed=1)
print(f"   best residual over 200 restarts: {b:.3e}   -> "
      f"{'REALISABLE' if b < 1e-10 else 'no realisation found'}")
