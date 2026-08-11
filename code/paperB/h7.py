"""enumerate all graphs on n=7 with alpha<=2 and omega<=3, reduce by isomorphism,
   then test spherical realisability (unit vectors, inner product 1/2 on edges)."""
import itertools, numpy as np
from scipy.optimize import least_squares

n=7
P=list(itertools.combinations(range(n),2))
TRI=[ (i,j,k) for i,j,k in itertools.combinations(range(n),3)]
QUAD=list(itertools.combinations(range(n),4))
pidx={frozenset(p):i for i,p in enumerate(P)}
def bit(i,j): return 1<<pidx[frozenset((i,j))]
TRIM=[bit(a,b)|bit(a,c)|bit(b,c) for a,b,c in TRI]
QUADM=[sum(bit(x,y) for x,y in itertools.combinations(q,2)) for q in QUAD]

cands=[]
for m in range(1<<len(P)):
    ok=True
    for t in TRIM:                 # alpha<=2 : every triple has an edge
        if m & t == 0: ok=False; break
    if not ok: continue
    for q in QUADM:                # omega<=3 : no K4
        if m & q == q: ok=False; break
    if not ok: continue
    cands.append(m)
print(f"n=7: {len(cands)} labelled graphs with alpha<=2, omega<=3")

# isomorphism reduction by canonical form under all 7! permutations
import math
perms=list(itertools.permutations(range(n)))
def canon(m):
    best=None
    for pm in perms:
        r=0
        for i,(a,b) in enumerate(P):
            if m>>i & 1: r |= bit(pm[a],pm[b])
        if best is None or r<best: best=r
    return best
seen={}
for m in cands:
    c=canon(m)
    seen.setdefault(c,m)
reps=list(seen.values())
print(f"      {len(reps)} up to isomorphism")

def edges_of(m): return [P[i] for i in range(len(P)) if m>>i & 1]
def realisable(edges, tries=120, seed=0):
    rng=np.random.default_rng(seed); E=list(edges)
    def resid(x):
        V=x.reshape(n,3)
        r=[float(np.dot(V[i],V[i]))-1 for i in range(n)]
        r+=[float(np.dot(V[i],V[j]))-0.5 for i,j in E]
        return np.array(r)
    best=None;bv=None
    for t in range(tries):
        s=least_squares(resid,rng.normal(size=3*n),xtol=1e-15,ftol=1e-15,gtol=1e-15)
        c=float(np.max(np.abs(s.fun)))
        if best is None or c<best: best,bv=c,s.x.reshape(n,3)
    sep=None
    if best<1e-10:
        NE=[p for p in P if frozenset(p) not in set(map(frozenset,E))]
        if NE: sep=min(abs(float(np.dot(bv[i],bv[j]))-0.5) for i,j in NE)
    return best,sep

found=[]
for r in reps:
    E=edges_of(r)
    b,s=realisable(E)
    if b<1e-10:
        found.append((r,len(E),s))
        print(f"   REALISABLE: {len(E)} edges, residual {b:.2e}, non-edge separation {s}")
print(f"\nn=7 realisable candidates: {len(found)}")
