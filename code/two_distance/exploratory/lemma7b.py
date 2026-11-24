"""Same test as lemma7.py, on Minkowski balls in Q(i, sqrt3, sqrt7) (float exploration).
Constraints for a 6-colouring under the hypothesis p_sqrt2 = p_{1+sqrt3} = 1/2:
  proper for distances 1, d;  for every P and Q, R with |PQ| = |PR| = delta in {sqrt2, 1+sqrt3}
  and |QR| in {1, d}:  c(P) = c(Q) or c(P) = c(R)."""
import sys, time, math, itertools, multiprocessing as mp
import numpy as np
from scipy.spatial import cKDTree
import ball37 as B

D2 = 2 + math.sqrt(3)
FORB = [1.0, D2]
DELT = [2.0, (1 + math.sqrt(3)) ** 2]

def classify(r2, targets, tol=1e-7):
    for t in targets:
        if abs(r2 - t) < tol: return True
    return False

def build(P, k, use_tri=True, deltas=DELT):
    xy = np.c_[P.real, P.imag]; T = cKDTree(xy); n = len(P)
    def pairs_at(r2):
        out = []
        for i, j in T.query_pairs(math.sqrt(r2) + 1e-6):
            if abs(((xy[i] - xy[j]) ** 2).sum() - r2) < 1e-7: out.append((i, j))
        return out
    edges = pairs_at(1.0) + pairs_at(D2)
    forbset = set(edges) | set((j, i) for i, j in edges)
    var = lambda v, c: v * k + c + 1
    nv = n * k; cls = []
    for v in range(n):
        cls.append([var(v, c) for c in range(k)])
        for a in range(k):
            for b in range(a + 1, k): cls.append([-var(v, a), -var(v, b)])
    for i, j in edges:
        for c in range(k): cls.append([-var(i, c), -var(j, c)])
    eqv = {}
    def eq(p, q):
        nonlocal nv
        key = (min(p, q), max(p, q))
        if key not in eqv:
            ys = []
            for c in range(k):
                nv += 1; ys.append(nv)
                cls.append([-nv, var(p, c)]); cls.append([-nv, var(q, c)])
            eqv[key] = ys
        return eqv[key]
    ntri = 0
    if use_tri:
        for dl in deltas:
            nb = {i: [] for i in range(n)}
            for i, j in pairs_at(dl):
                nb[i].append(j); nb[j].append(i)
            for p in range(n):
                for q, r in itertools.combinations(nb[p], 2):
                    if (q, r) in forbset:
                        cls.append(eq(p, q) + eq(p, r)); ntri += 1
    return edges, cls, ntri, var

def worker(cls, q):
    from pysat.solvers import Solver
    s = Solver(name="cadical195", bootstrap_with=cls)
    r = s.solve(); q.put((r, s.get_model() if r else None))

def run(P, k, use_tri, label, tl=240):
    t0 = time.time()
    edges, cls, ntri, var = build(P, k, use_tri)
    if edges:
        i, j = edges[0]; cls = cls + [[var(i, 0)], [var(j, 1)]]
    q = mp.Queue(); p = mp.Process(target=worker, args=(cls, q)); p.start(); p.join(tl)
    if p.is_alive():
        p.kill(); res = "timeout"
    else:
        r, m = q.get()
        if r:
            ms = set(l for l in m if l > 0)
            col = [next(c for c in range(k) if var(v, c) in ms) for v in range(len(P))]
            assert all(col[a] != col[b] for a, b in edges)
        res = "SAT" if r else "UNSAT"
    print(f"{label}: |F|={len(P)} forb-edges={len(edges)} triangles={ntri} -> {res} ({time.time()-t0:.0f}s)", flush=True)
    return res

if __name__ == "__main__":
    mode = sys.argv[1]
    U = B.units(1); Dv = [(1 + B.zeta) * u for u in U]
    if mode == "U1":      # sums of <= 2 unit vectors only
        P = B.ball(U, 2)
    elif mode == "UD1":   # depth 1 of units + d-vectors
        P = B.ball(U + Dv, 1)
    elif mode == "UD2":
        P = B.ball(U + Dv, 2)
    elif mode == "U3":
        P = B.ball(U, 3)
    run(P, 6, False, mode + " control(no tri)")
    run(P, 6, True, mode + " with tri")
