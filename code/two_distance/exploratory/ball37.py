"""{1,d}-graphs, d = (sqrt2+sqrt6)/2, on Minkowski balls in Q(i, sqrt3, sqrt7).
Exploration in floating point (vertex hashing at 1e-9); results to be certified exactly later."""
import sys, time, threading, cmath, math, itertools
import numpy as np
from scipy.spatial import cKDTree
import networkx as nx
from pysat.solvers import Glucose4, Cadical153

D2 = 2 + math.sqrt(3)
zeta = cmath.exp(1j * math.pi / 6)
omega = (3 + 1j * math.sqrt(7)) / 4          # |omega| = 1, in Q(sqrt-7)
assert abs(abs(omega) - 1) < 1e-15 and abs(abs(1 + zeta) ** 2 - D2) < 1e-12

def units(level):
    rots = [zeta ** k for k in range(12)]
    ws = [omega ** a for a in range(-level, level + 1)]
    U = {}
    for r in rots:
        for w in ws:
            z = r * w
            U[(round(z.real, 9), round(z.imag, 9))] = z
    return list(U.values())

def ball(gens, depth):
    pts = {(0.0, 0.0): 0j}
    frontier = [0j]
    for _ in range(depth):
        new = []
        for p in frontier:
            for g in gens:
                q = p + g
                key = (round(q.real, 7), round(q.imag, 7))
                if key not in pts:
                    pts[key] = q; new.append(q)
        frontier = new
    return np.array(list(pts.values()))

def graph(P):
    xy = np.c_[P.real, P.imag]
    T = cKDTree(xy)
    E = set()
    for r2 in (1.0, D2):
        for i, j in T.query_pairs(math.sqrt(r2) + 1e-6):
            d = (xy[i, 0] - xy[j, 0]) ** 2 + (xy[i, 1] - xy[j, 1]) ** 2
            if abs(d - 1) < 1e-7 or abs(d - D2) < 1e-7:
                E.add((i, j))
    G = nx.Graph(); G.add_nodes_from(range(len(P))); G.add_edges_from(E)
    return G

def kcore_reduce(G, k):
    """vertices of degree < k can always be coloured last: drop them."""
    return nx.k_core(G, k)

def solve(G, k, tl):
    nodes = list(G.nodes()); ix = {v: i for i, v in enumerate(nodes)}
    var = lambda v, c: ix[v] * k + c + 1
    s = Glucose4()
    for v in nodes: s.add_clause([var(v, c) for c in range(k)])
    for u, v in G.edges():
        for c in range(k): s.add_clause([-var(u, c), -var(v, c)])
    cl = max(nx.find_cliques(G), key=len)
    for i, v in enumerate(cl[:k]): s.add_clause([var(v, i)])
    t = threading.Timer(tl, s.interrupt); t.start()
    r = s.solve_limited(expect_interrupt=True); t.cancel()
    if r is None:
        return "timeout", len(cl)
    if r:
        m = set(l for l in s.get_model() if l > 0)
        col = {v: next(c for c in range(k) if var(v, c) in m) for v in nodes}
        assert all(col[u] != col[v] for u, v in G.edges())
        return "SAT", len(cl)
    return "UNSAT", len(cl)

if __name__ == "__main__":
    level, depth, kk = int(sys.argv[1]), int(sys.argv[2]), [int(x) for x in sys.argv[3].split(",")]
    tl = float(sys.argv[4]) if len(sys.argv) > 4 else 100
    U = units(level)
    Dv = [(1 + zeta) * u for u in U]
    gens = U + Dv
    t0 = time.time()
    P = ball(gens, depth)
    G = graph(P)
    print(f"level={level} depth={depth} |U|={len(U)} gens={len(gens)} |V|={len(P)} |E|={G.number_of_edges()} ({time.time()-t0:.1f}s)", flush=True)
    for k in kk:
        H = kcore_reduce(G, k)
        r, om = solve(H, k, tl)
        print(f"  k={k}: core |V|={H.number_of_nodes()} |E|={H.number_of_edges()} omega={om} -> {r} ({time.time()-t0:.1f}s)", flush=True)
