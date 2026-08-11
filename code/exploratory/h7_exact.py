"""
h7_exact.py — step D3.b3, part (a): lift h = 7 from numerics to exact.

h7.py found that exactly one graph on 7 vertices (up to isomorphism) with
alpha <= 2 and omega <= 3 is realisable on the unit sphere with every edge a
chord of length 1 (inner product 1/2).  It has 11 edges and the numerical
solution had residual 2.2e-16 with non-edges separated from 1/2 by 0.049.

Numerics are not a theorem.  This script recovers the graph, solves for a
realisation, and inspects the Gram matrix: if the non-edge inner products come
out at recognisable algebraic values, the exact realisation follows directly.
"""
import itertools, sys
import numpy as np
from scipy.optimize import least_squares
from fractions import Fraction

n = 7
P = list(itertools.combinations(range(n), 2))
pidx = {frozenset(p): i for i, p in enumerate(P)}
def bit(i, j): return 1 << pidx[frozenset((i, j))]
TRIM = [bit(a, b) | bit(a, c) | bit(b, c) for a, b, c in itertools.combinations(range(n), 3)]
QUADM = [sum(bit(x, y) for x, y in itertools.combinations(q, 2))
         for q in itertools.combinations(range(n), 4)]

cands = []
for m in range(1 << len(P)):
    if any(m & t == 0 for t in TRIM):   # alpha <= 2
        continue
    if any(m & q == q for q in QUADM):  # omega <= 3
        continue
    cands.append(m)

perms = list(itertools.permutations(range(n)))
def canon(m):
    best = None
    for pm in perms:
        r = 0
        for i, (a, b) in enumerate(P):
            if m >> i & 1:
                r |= bit(pm[a], pm[b])
        if best is None or r < best:
            best = r
    return best

seen = {}
for m in cands:
    seen.setdefault(canon(m), m)
reps = list(seen.values())
print(f"{len(cands)} labelled, {len(reps)} up to isomorphism")


def realise(edges, tries=400, seed=1):
    rng = np.random.default_rng(seed)
    E = list(edges)
    def resid(x):
        V = x.reshape(n, 3)
        r = [float(np.dot(V[i], V[i]))-1 for i in range(n)]
        r += [float(np.dot(V[i], V[j]))-0.5 for i, j in E]
        return np.array(r)
    best = None; bv = None
    for _ in range(tries):
        s = least_squares(resid, rng.normal(size=3*n),
                          xtol=1e-15, ftol=1e-15, gtol=1e-15)
        c = float(np.max(np.abs(s.fun)))
        if best is None or c < best:
            best, bv = c, s.x.reshape(n, 3)
    return best, bv


good = []
for r in reps:
    E = [P[i] for i in range(len(P)) if r >> i & 1]
    b, V = realise(E)
    if b < 1e-10:
        good.append((r, E, V, b))
print(f"realisable up to isomorphism: {len(good)}")

for r, E, V, b in good:
    print(f"\ngraph with {len(E)} edges, residual {b:.2e}")
    G = V @ V.T
    NE = [p for p in P if frozenset(p) not in set(map(frozenset, E))]
    print("  non-edge inner products:")
    vals = sorted({round(float(G[i, j]), 9) for i, j in NE})
    for v in vals:
        # try to recognise the value
        guess = ""
        for den in range(1, 25):
            for num in range(-den*3, den*3+1):
                if abs(v - num/den) < 1e-9:
                    guess = f"= {num}/{den}"
                    break
            if guess: break
        if not guess:
            for expr, name in ((-1/3, "-1/3"), (1/np.sqrt(5), "1/sqrt5"),
                               (-1/np.sqrt(5), "-1/sqrt5"),
                               ((np.sqrt(5)-1)/4, "(sqrt5-1)/4"),
                               (-(np.sqrt(5)+1)/4, "-(sqrt5+1)/4")):
                if abs(v - expr) < 1e-9:
                    guess = f"= {name}"; break
        print(f"    {v:+.9f}  {guess}")
    print("  degree sequence:",
          sorted(sum(1 for e in E if i in e) for i in range(n)))
    np.save("h7_realisation.npy", V)
