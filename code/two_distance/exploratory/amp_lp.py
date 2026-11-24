"""Idea 1: 'amplification' LP for chi({1,d}) in Z[zeta12].

Variables p_delta (averaged probability that a pair at distance delta is
monochromatic, w.r.t. an invariant mean).  Constraints:
  * p = 0 on forbidden distances (1, d)
  * gadgets X: (p_{|x-y|})_{pairs} in conv(partition patterns of X into <= k
    classes, each class containing no forbidden pair)
  * triangles: p_a + p_b - p_c <= 1
Objective: minimise M >= p_delta for all delta > 1/2.
M* >= 1/2  contradicts the strict lemma p_delta < 1/2  =>  chi > k.
"""
import itertools, math, sys, time
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix

# ---------- exact arithmetic in Z[zeta], zeta^4 = zeta^2 - 1 ----------
def mul(a, b):
    c = [0] * 7
    for i in range(4):
        for j in range(4):
            c[i + j] += a[i] * b[j]
    for k in range(6, 3, -1):        # zeta^k = zeta^{k-2} - zeta^{k-4}
        c[k - 2] += c[k]; c[k - 4] -= c[k]; c[k] = 0
    return tuple(c[:4])

ZETA = (0, 1, 0, 0)
UNITS = [(1, 0, 0, 0)]
for _ in range(11):
    UNITS.append(mul(UNITS[-1], ZETA))
assert mul(UNITS[-1], ZETA) == (1, 0, 0, 0)

def add(a, b): return tuple(x + y for x, y in zip(a, b))
def sub(a, b): return tuple(x - y for x, y in zip(a, b))

E1 = [np.exp(1j * math.pi / 6 * k) for k in range(4)]      # zeta -> e^{i pi/6}
E2 = [np.exp(5j * math.pi / 6 * k) for k in range(4)]      # zeta -> zeta^5 (sqrt3 -> -sqrt3)

def normkey(a):
    """|a|^2 = m + n sqrt3 exactly, returned as (m, n)."""
    z1 = sum(c * e for c, e in zip(a, E1)); z2 = sum(c * e for c, e in zip(a, E2))
    N1, N2 = abs(z1) ** 2, abs(z2) ** 2
    m = (N1 + N2) / 2; n = (N1 - N2) / (2 * math.sqrt(3))
    mi, ni = round(m), round(n)
    assert abs(m - mi) < 1e-7 and abs(n - ni) < 1e-7, (a, m, n)
    return (mi, ni)

def kval(key): return key[0] + key[1] * math.sqrt(3)

def ball(depth):
    pts = {(0, 0, 0, 0)}
    frontier = set(pts)
    for _ in range(depth):
        new = set()
        for p in frontier:
            for u in UNITS:
                q = add(p, u)
                if q not in pts:
                    new.add(q)
        pts |= new; frontier = new
    return sorted(pts)

# ---------- gadgets ----------
def find_gadgets(S, forb, size, maxdef, limit=200000):
    """Sets X (containing 0) of given size in which every vertex has at most
    `maxdef` non-forbidden partners inside X.  Returned modulo duplicates."""
    idx = {p: i for i, p in enumerate(S)}
    n = len(S)
    adj = [set() for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if normkey(sub(S[i], S[j])) in forb:
                adj[i].add(j); adj[j].add(i)
    zero = idx[(0, 0, 0, 0)]
    out = set()

    def rec(X, cands):
        if len(out) >= limit: return
        if len(X) == size:
            out.add(tuple(sorted(X))); return
        for c in sorted(cands):
            if c <= X[-1] and len(X) > 1:
                continue
            # deficiency check
            newX = X + [c]
            ok = True
            for v in newX:
                d = sum(1 for w in newX if w != v and w not in adj[v])
                if d > maxdef: ok = False; break
            if not ok: continue
            # candidates: vertices adjacent to all but <= maxdef of newX
            nc = {w for w in cands if w > c and
                  sum(1 for v in newX if w not in adj[v]) <= maxdef}
            if len(nc) + len(newX) >= size:
                rec(newX, nc)

    base = [zero]
    cands0 = {w for w in range(n) if w != zero}
    # first extra vertex must be forbidden-adjacent to 0 or a defect; pick all
    for c in sorted(cands0):
        X = [zero, c]
        if (c not in adj[zero]) and maxdef == 0: continue
        nc = {w for w in cands0 if w > c and
              sum(1 for v in X if w not in adj[v]) <= maxdef}
        rec(X, nc)
    return [[S[i] for i in g] for g in out], adj

def patterns(X, forb, k):
    """All set partitions of X into <= k classes with no forbidden pair inside."""
    n = len(X)
    bad = [[normkey(sub(X[i], X[j])) in forb if i != j else False for j in range(n)] for i in range(n)]
    res = []
    def rec(i, classes):
        if i == n:
            res.append([c[:] for c in classes]); return
        for c in classes:
            if not any(bad[i][j] for j in c):
                c.append(i); rec(i + 1, classes); c.pop()
        if len(classes) < k:
            classes.append([i]); rec(i + 1, classes); classes.pop()
    rec(0, [])
    return res

# ---------- LP ----------
def run(dkey, k, depth=2, gsizes=None, maxdef=1, tri_depth=None, verbose=True):
    t0 = time.time()
    forb = {(1, 0), dkey}
    S = ball(depth)
    keys = {}
    def var(key):
        if key not in keys: keys[key] = len(keys)
        return keys[key]
    rows = []      # (dict var->coef, sense, rhs)  sense: 'le' or 'eq'
    lam_count = [0]

    # gadgets
    gsizes = gsizes or [k + 1, k + 2]
    ngad = 0
    for sz in gsizes:
        gads, _ = find_gadgets(S, forb, sz, maxdef)
        for X in gads:
            pats = patterns(X, forb, k)
            if not pats:   # X itself not k-colourable -> chi > k outright
                return dict(status="gadget not k-colourable", X=X)
            # lambda variables for this gadget
            lam0 = lam_count[0]; lam_count[0] += len(pats)
            pairs = [(i, j) for i in range(len(X)) for j in range(i + 1, len(X))
                     if normkey(sub(X[i], X[j])) not in forb]
            # sum lambda = 1
            rows.append(({('l', lam0 + t): 1 for t in range(len(pats))}, 'eq', 1))
            for (i, j) in pairs:
                key = normkey(sub(X[i], X[j]))
                r = {('p', var(key)): -1}
                for t, P in enumerate(pats):
                    if any(i in c and j in c for c in P):
                        r[('l', lam0 + t)] = r.get(('l', lam0 + t), 0) + 1
                rows.append((r, 'eq', 0))
            ngad += 1

    # triangles with a vertex at 0
    T = ball(tri_depth if tri_depth is not None else depth)
    tris = set()
    nz = [p for p in T if any(p)]
    for a in range(len(nz)):
        ka = normkey(nz[a])
        for b in range(a + 1, len(nz)):
            kb = normkey(nz[b]); kc = normkey(sub(nz[a], nz[b]))
            tris.add(tuple(sorted((ka, kb, kc))))
    for (a, b, c) in tris:
        for (x, y, z) in ((a, b, c), (a, c, b), (b, c, a)):
            r = {}
            for kk, co in ((x, 1), (y, 1), (z, -1)):
                r[('p', var(kk))] = r.get(('p', var(kk)), 0) + co
            rows.append((r, 'le', 1))

    # forbidden distances
    for f in forb:
        if f in keys:
            rows.append(({('p', keys[f]): 1}, 'eq', 0))

    # M >= p_delta for delta > 1/2
    big = [kk for kk in keys if kval(kk) > 0.25 + 1e-12 and kk not in forb]
    for kk in big:
        rows.append(({('p', keys[kk]): 1, ('M', 0): -1}, 'le', 0))

    nP = len(keys); nL = lam_count[0]; nv = nP + nL + 1
    def col(v):
        return v[1] if v[0] == 'p' else (nP + v[1] if v[0] == 'l' else nP + nL)
    ub = [r for r in rows if r[1] == 'le']; eq = [r for r in rows if r[1] == 'eq']
    A_ub = lil_matrix((len(ub), nv)); b_ub = np.zeros(len(ub))
    for i, (r, _, rhs) in enumerate(ub):
        for v, co in r.items(): A_ub[i, col(v)] = co
        b_ub[i] = rhs
    A_eq = lil_matrix((len(eq), nv)); b_eq = np.zeros(len(eq))
    for i, (r, _, rhs) in enumerate(eq):
        for v, co in r.items(): A_eq[i, col(v)] = co
        b_eq[i] = rhs
    c = np.zeros(nv); c[-1] = 1
    res = linprog(c, A_ub=A_ub.tocsr(), b_ub=b_ub, A_eq=A_eq.tocsr(), b_eq=b_eq,
                  bounds=[(0, 1)] * nv, method="highs")
    out = dict(status=res.message, M=res.fun if res.success else None,
               nS=len(S), nkeys=nP, ngad=ngad, ntri=len(tris), nlam=nL,
               time=round(time.time() - t0, 1))
    if res.success:
        p = res.x[:nP]
        out["top"] = sorted(((round(p[keys[kk]], 4), round(math.sqrt(kval(kk)), 4)) for kk in big),
                            reverse=True)[:6]
    return out

if __name__ == "__main__":
    D = {"sqrt3": (3, 0), "d": (2, 1), "2": (4, 0), "sqrt2": (2, 0)}
    dname, k, depth = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    maxdef = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    tri = int(sys.argv[5]) if len(sys.argv) > 5 else None
    print(dname, "k=", k, "depth=", depth, "maxdef=", maxdef, "tri_depth=", tri)
    print(run(D[dname], k, depth, maxdef=maxdef, tri_depth=tri))
