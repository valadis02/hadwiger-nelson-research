import itertools, math, pickle, time, random, sys
import numpy as np
from scipy.spatial import cKDTree
from pysat.solvers import Solver
import ball37 as B
import lemma7b as L

K = 6
U = B.units(1); Dv = [(1 + B.zeta) * u for u in U]
P = B.ball(U + Dv, 2); n = len(P)
xy = np.c_[P.real, P.imag]; T = cKDTree(xy)
def pairs_at(r2):
    return [(i, j) for i, j in T.query_pairs(math.sqrt(r2) + 1e-6) if abs(((xy[i] - xy[j]) ** 2).sum() - r2) < 1e-7]
E = pairs_at(1.0) + pairs_at(L.D2)
Eset = set(E) | set((j, i) for i, j in E)
TRI = []
for dl in L.DELT:
    nb = {i: [] for i in range(n)}
    for i, j in pairs_at(dl): nb[i].append(j); nb[j].append(i)
    for p in range(n):
        for q, r in itertools.combinations(nb[p], 2):
            if (q, r) in Eset: TRI.append((p, q, r))
log = lambda *a: print(*a, f"[{time.time()-t0:.0f}s]", flush=True)
t0 = time.time()
log("points", n, "edges", len(E), "triangles", len(TRI))

def encode(verts, tri_sel=False, vert_sel=False, tris=None):
    vs = set(verts); var = lambda v, c: v * K + c + 1
    nv = n * K; cls = []; tsel = {}; vsel = {}
    if vert_sel:
        for v in verts: nv += 1; vsel[v] = nv
    for v in verts:
        g = [-vsel[v]] if vert_sel else []
        cls.append(g + [var(v, c) for c in range(K)])
        for a in range(K):
            for b in range(a + 1, K): cls.append([-var(v, a), -var(v, b)])
    for i, j in E:
        if i in vs and j in vs:
            for c in range(K): cls.append([-var(i, c), -var(j, c)])
    eqv = {}
    def eq(p, q):
        nonlocal nv
        key = (min(p, q), max(p, q))
        if key not in eqv:
            ys = []
            for c in range(K):
                nv += 1; ys.append(nv); cls.append([-nv, var(p, c)]); cls.append([-nv, var(q, c)])
            eqv[key] = ys
        return eqv[key]
    for t in (tris if tris is not None else TRI):
        p, q, r = t
        if p in vs and q in vs and r in vs:
            g = []
            if tri_sel: nv += 1; tsel[t] = nv; g = [-nv]
            if vert_sel: g += [-vsel[p], -vsel[q], -vsel[r]]
            cls.append(g + eq(p, q) + eq(p, r))
    # symmetry breaking on one edge inside verts
    for i, j in E:
        if i in vs and j in vs:
            cls.append([var(i, 0)]); cls.append([var(j, 1)]); break
    return cls, tsel, vsel

stage = sys.argv[1]
if stage == "1":
    verts = list(range(n))
    cls, tsel, _ = encode(verts, tri_sel=True)
    s = Solver(name="cadical195", bootstrap_with=cls)
    r = s.solve(assumptions=list(tsel.values()))
    log("full with triangle selectors:", "SAT" if r else "UNSAT")
    inv = {l: t for t, l in tsel.items()}
    Tc = [inv[l] for l in s.get_core() if l in inv]
    log("core triangles:", len(Tc))
    pickle.dump(Tc, open("tricore.pkl", "wb"))
elif stage == "2":
    Tc = pickle.load(open("tricore.pkl", "rb"))
    for it in range(4):
        cls, tsel, _ = encode(list(range(n)), tri_sel=True, tris=Tc)
        s = Solver(name="cadical195", bootstrap_with=cls)
        r = s.solve(assumptions=list(tsel.values()))
        if r: log("SAT?!"); break
        inv = {l: t for t, l in tsel.items()}
        Tn = [inv[l] for l in s.get_core() if l in inv]
        log(f"iter {it}: triangles {len(Tc)} -> {len(Tn)}")
        if len(Tn) >= len(Tc): break
        Tc = Tn
    pickle.dump(Tc, open("tricore.pkl", "wb"))
    # now vertex selectors, triangles fixed to Tc
    verts = list(range(n))
    for it in range(6):
        cls, _, vsel = encode(verts, vert_sel=True, tris=Tc)
        s = Solver(name="cadical195", bootstrap_with=cls)
        r = s.solve(assumptions=[vsel[v] for v in verts])
        if r: log("SAT?!"); break
        inv = {l: v for v, l in vsel.items()}
        vn = sorted(inv[l] for l in s.get_core() if l in inv)
        log(f"vertex iter {it}: {len(verts)} -> {len(vn)}")
        if len(vn) >= len(verts): break
        verts = vn
    pickle.dump({"verts": verts, "tris": Tc}, open("core7_stage2.pkl", "wb"))
elif stage == "3":
    d = pickle.load(open("core7_stage2.pkl", "rb")); verts, Tc = d["verts"], d["tris"]
    def unsat(vs):
        cls, _, _ = encode(vs, tris=Tc)
        return not Solver(name="cadical195", bootstrap_with=cls).solve()
    random.seed(1)
    for v in sorted(verts, key=lambda v: random.random()):
        trial = [w for w in verts if w != v]
        if unsat(trial): verts = trial
    vs = set(verts)
    used_E = [(i, j) for i, j in E if i in vs and j in vs]
    used_T = [t for t in Tc if all(x in vs for x in t)]
    log("final vertices", len(verts), "edges", len(used_E), "triangles", len(used_T))
    pickle.dump({"P": P, "verts": verts, "edges": used_E, "tris": used_T}, open("core7.pkl", "wb"))
